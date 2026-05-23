#include <Rcpp.h>

#include <zlib.h>
#include <htslib/kseq.h>

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

KSEQ_INIT(gzFile, gzread)

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

bool append_segments(std::string &out, const std::vector<Segment> &segments, const kseq_t *seq) {
    for (const Segment &segment : segments) {
        if (seq->seq.l < static_cast<std::size_t>(segment.start + segment.width))
            return false;
        out += uppercase_subseq(seq->seq.s, segment.start, segment.width);
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

} // namespace

// [[Rcpp::export]]
Rcpp::List count_fastq_barcodes_cpp(Rcpp::CharacterVector fq,
                                    Rcpp::CharacterVector sample_barcodes,
                                    Rcpp::CharacterVector construct_barcodes,
                                    Rcpp::IntegerVector sample_start,
                                    Rcpp::IntegerVector sample_width,
                                    Rcpp::IntegerVector construct_start,
                                    Rcpp::IntegerVector construct_width,
                                    bool verbose) {
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

    const std::string path = Rcpp::as<std::string>(fq[0]);
    if (verbose)
        Rcpp::Rcout << "Counting barcodes in " << path << "\n";

    gzFile fp = gzopen(path.c_str(), "r");
    if (fp == nullptr)
        Rcpp::stop("Could not open FASTQ file '%s'", path);
    kseq_t *seq = kseq_init(fp);

    while (kseq_read(seq) >= 0) {
        ++reads_seen;
        if (reads_seen % 1000000 == 0) {
            Rcpp::checkUserInterrupt();
            if (verbose)
                Rcpp::Rcout << "Processed " << reads_seen << " reads\n";
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

    kseq_destroy(seq);
    gzclose(fp);

    return Rcpp::List::create(
        Rcpp::_["counts"] = counts,
        Rcpp::_["total_reads"] = total_reads,
        Rcpp::_["mapped_reads"] = mapped_reads,
        Rcpp::_["unmatched_reads"] = unmatched_reads,
        Rcpp::_["too_short_reads"] = too_short_reads
    );
}
