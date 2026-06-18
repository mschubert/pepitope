#include <Rcpp.h>

#include <zlib.h>

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <ostream>
#include <string>
#include <unordered_map>
#include <vector>

namespace {

struct Segment {
    int start;
    int width;
};

std::string uppercase_subseq(const char *seq, int start, int width) {
    std::string out(seq + start, seq + start + width);
    std::transform(out.begin(), out.end(), out.begin(), [](unsigned char c) {
        return static_cast<char>(std::toupper(c));
    });
    return out;
}

std::vector<Segment> as_segments(Rcpp::IntegerVector start, Rcpp::IntegerVector width) {
    std::vector<Segment> out;
    out.reserve(start.size());
    for (R_xlen_t i = 0; i < start.size(); ++i)
        out.push_back({start[i] - 1, width[i]});
    return out;
}

bool read_gz_line(gzFile fp, std::string &out) {
    out.clear();
    char buffer[65536];
    while (gzgets(fp, buffer, sizeof(buffer)) != Z_NULL) {
        out += buffer;
        if (!out.empty() && out.back() == '\n')
            break;
    }
    if (out.empty())
        return false;
    while (!out.empty() && (out.back() == '\n' || out.back() == '\r'))
        out.pop_back();
    return true;
}

bool append_segments(std::string &out, const std::vector<Segment> &segments, const std::string &seq) {
    for (const Segment &segment : segments) {
        if (seq.size() < static_cast<std::size_t>(segment.start + segment.width))
            return false;
        out += uppercase_subseq(seq.data(), segment.start, segment.width);
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
    Rcpp::NumericVector total_reads(samples.size());
    Rcpp::NumericVector mapped_reads(samples.size());
    double unmatched_reads = 0.0;
    double too_short_reads = 0.0;
    std::uint64_t reads_seen = 0;
    Rcpp::Function format("format");

    const std::string path = Rcpp::as<std::string>(fq[0]);
    if (verbose)
        Rcpp::Rcout << "Counting barcodes in " << path << "\n";

    gzFile fp = gzopen(path.c_str(), "r");
    if (fp == nullptr)
        Rcpp::stop("Could not open FASTQ file '%s'", path);

    std::string header;
    std::string seq;
    std::string plus;
    std::string qual;
    while (read_gz_line(fp, header)) {
        if (!read_gz_line(fp, seq) || !read_gz_line(fp, plus) || !read_gz_line(fp, qual)) {
            gzclose(fp);
            Rcpp::stop("Malformed FASTQ file '%s'", path);
        }
        ++reads_seen;
        if (reads_seen % 1000000 == 0) {
            Rcpp::checkUserInterrupt();
            if (progress)
                Rcpp::Rcout << "\rProcessed " << format_reads(reads_seen, format) << " reads" << std::flush;
        }

        std::string sample_key;
        std::string construct_key;
        if (!append_segments(sample_key, sample_segments, seq) ||
                !append_segments(construct_key, construct_segments, seq)) {
            ++too_short_reads;
            continue;
        }

        const auto sample_hit = sample_index.find(sample_key);
        if (sample_hit == sample_index.end()) {
            ++unmatched_reads;
            continue;
        }

        const int sample_i = sample_hit->second;
        ++total_reads[sample_i];
        const auto construct_hit = construct_index.find(construct_key);
        if (construct_hit == construct_index.end())
            continue;

        ++mapped_reads[sample_i];
        ++counts(construct_hit->second, sample_i);
    }

    gzclose(fp);

    if (verbose)
        Rcpp::Rcout << "\rProcessed " << format_reads(reads_seen, format) << " reads\n";

    return Rcpp::List::create(
        Rcpp::_["counts"] = counts,
        Rcpp::_["total_reads"] = total_reads,
        Rcpp::_["mapped_reads"] = mapped_reads,
        Rcpp::_["unmatched_reads"] = unmatched_reads,
        Rcpp::_["too_short_reads"] = too_short_reads,
        Rcpp::_["reads_seen"] = static_cast<double>(reads_seen)
    );
}
