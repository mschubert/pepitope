#include <Rcpp.h>

#include <zlib.h>

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <ostream>
#include <string>
#include <unordered_map>
#include <vector>

namespace {

struct Segment {
    int start;
    int width;
};

std::vector<Segment> as_segments(Rcpp::IntegerVector start, Rcpp::IntegerVector width) {
    std::vector<Segment> out;
    out.reserve(start.size());
    for (R_xlen_t i = 0; i < start.size(); ++i)
        out.push_back({start[i] - 1, width[i]});
    return out;
}

class FastqReader {
    const std::string path;
    gzFile fp;
    std::vector<char> buffer;
    std::size_t pos = 0, end = 0;

public:
    explicit FastqReader(const std::string &path) : path(path), fp(nullptr), buffer(1 << 20) {
        fp = gzopen(path.c_str(), "r");
        if (fp == nullptr)
            Rcpp::stop("Could not open FASTQ file '%s'", path);
        gzbuffer(fp, 1 << 20);
    }

    ~FastqReader() { gzclose(fp); }
    FastqReader(const FastqReader &) = delete;
    FastqReader &operator=(const FastqReader &) = delete;

    // Passing no output skips a line without copying its contents.
    bool line(std::string *out = nullptr) {
        if (out)
            out->clear();
        bool seen = false;
        for (;;) {
            if (pos == end) {
                const int size = gzread(fp, buffer.data(), buffer.size());
                // A truncated stream can return data and set an error in the same read.
                int status;
                const char *message = gzerror(fp, &status);
                if (status != Z_OK && status != Z_STREAM_END)
                    Rcpp::stop("Could not read FASTQ file '%s': %s", path, message);
                if (size <= 0)
                    break;
                pos = 0;
                end = size;
            }
            const char *begin = buffer.data() + pos;
            const char *newline = static_cast<const char *>(std::memchr(begin, '\n', end-pos));
            const std::size_t length = newline ? static_cast<std::size_t>(newline-begin) : end-pos;
            if (out)
                out->append(begin, length);
            seen = true;
            pos += length;
            if (newline) {
                ++pos;
                break;
            }
        }
        if (out)
            while (!out->empty() && out->back() == '\r')
                out->pop_back();
        return seen;
    }
};

bool append_segments(std::string &out, const std::vector<Segment> &segments, const std::string &seq) {
    for (const Segment &segment : segments) {
        if (seq.size() < static_cast<std::size_t>(segment.start + segment.width))
            return false;
        const auto previous = out.size();
        out.append(seq, segment.start, segment.width);
        std::transform(out.begin() + previous, out.end(), out.begin() + previous,
            [](unsigned char c) {
                return static_cast<char>(c >= 'a' && c <= 'z' ? c - ('a' - 'A') : c);
            });
    }
    return true;
}

std::unordered_map<std::string, int> make_index(const std::vector<std::string> &values) {
    std::unordered_map<std::string, int> index;
    index.reserve(values.size() * 2);
    for (std::size_t i = 0; i < values.size(); ++i)
        index.emplace(values[i], static_cast<int>(i));
    return index;
}

std::string format_reads(std::uint64_t value, Rcpp::Function formatter) {
    Rcpp::CharacterVector formatted = formatter(
        static_cast<double>(value),
        Rcpp::Named("big.mark") = ",",
        Rcpp::Named("scientific") = false,
        Rcpp::Named("trim") = true
    );
    return Rcpp::as<std::string>(formatted[0]);
}

} // namespace

// [[Rcpp::export]]
Rcpp::List count_fastq_barcodes_cpp(Rcpp::CharacterVector fq,
                                    Rcpp::CharacterVector sample_barcodes,
                                    Rcpp::CharacterVector construct_barcodes,
                                    Rcpp::IntegerVector sample_start,
                                    Rcpp::IntegerVector sample_width,
                                    Rcpp::IntegerVector construct_start,
                                    Rcpp::IntegerVector construct_width,
                                    bool verbose,
                                    bool progress) {
    const std::vector<std::string> samples = Rcpp::as<std::vector<std::string>>(sample_barcodes);
    const std::vector<std::string> constructs = Rcpp::as<std::vector<std::string>>(construct_barcodes);

    const std::vector<Segment> sample_segments = as_segments(sample_start, sample_width);
    const std::vector<Segment> construct_segments = as_segments(construct_start, construct_width);
    const auto sample_index = make_index(samples);
    const auto construct_index = make_index(constructs);

    Rcpp::NumericMatrix counts(constructs.size(), samples.size());
    Rcpp::NumericVector unassigned_reads(constructs.size());
    Rcpp::NumericVector total_reads(samples.size());
    double unmapped_reads = 0.0;
    std::uint64_t reads_seen = 0;
    Rcpp::Function format("format");

    const std::string path = Rcpp::as<std::string>(fq[0]);
    if (verbose)
        Rcpp::Rcout << "Counting barcodes in " << path << "\n";

    std::string sample_key;
    std::string construct_key;
    std::size_t sample_size = 0, construct_size = 0;
    for (const auto &segment : sample_segments)
        sample_size += segment.width;
    for (const auto &segment : construct_segments)
        construct_size += segment.width;
    sample_key.reserve(sample_size);
    construct_key.reserve(construct_size);

    FastqReader reader(path);
    std::string seq;
    while (reader.line()) {
        if (!reader.line(&seq) || !reader.line() || !reader.line())
            Rcpp::stop("Malformed FASTQ file '%s'", path);
        ++reads_seen;
        if (reads_seen % 1000000 == 0) {
            Rcpp::checkUserInterrupt();
            if (progress)
                Rcpp::Rcout << "\rProcessed " << format_reads(reads_seen, format) << " reads" << std::flush;
        }

        sample_key.clear();
        construct_key.clear();
        if (!append_segments(sample_key, sample_segments, seq) ||
                !append_segments(construct_key, construct_segments, seq))
            continue;

        const auto sample_hit = sample_index.find(sample_key);
        const auto construct_hit = construct_index.find(construct_key);
        if (sample_hit == sample_index.end()) {
            if (construct_hit != construct_index.end())
                ++unassigned_reads[construct_hit->second];
            else
                ++unmapped_reads;
            continue;
        }

        const int sample_i = sample_hit->second;
        ++total_reads[sample_i];
        if (construct_hit == construct_index.end())
            continue;

        ++counts(construct_hit->second, sample_i);
    }

    if (verbose)
        Rcpp::Rcout << "\rProcessed " << format_reads(reads_seen, format) << " reads\n";

    return Rcpp::List::create(
        Rcpp::_["counts"] = counts,
        Rcpp::_["unassigned_reads"] = unassigned_reads,
        Rcpp::_["total_reads"] = total_reads,
        Rcpp::_["unmapped_reads"] = unmapped_reads
    );
}
