#include <Rcpp.h>

#include <algorithm>
#include <array>
#include <limits>
#include <string>
#include <vector>

// For motif-aware codon optimization by dynamic programming, see:
// Condon & Thachuk (2012), Efficient codon optimization with motif engineering.
// https://www.cs.ubc.ca/labs/algorithms/Projects/codon-optimizer/downloads/jda2012-codon-opt.pdf
// This implementation uses motif-prefix states and minimizes nucleotide substitutions;
// the paper optimizes the codon adaptation index (CAI).

namespace {

const std::string bases = "ACGT";

int codon_index(const std::string &codon) {
    int index = 0;
    for (char base : codon) {
        const std::size_t digit = bases.find(base);
        if (digit == std::string::npos)
            Rcpp::stop("Sequences containing cut sites must use only A, C, G, and T");
        index = index * 4 + static_cast<int>(digit);
    }
    return index;
}

bool ends_with(const std::string &text, const std::string &suffix) {
    return text.size() >= suffix.size() &&
        text.compare(text.size() - suffix.size(), suffix.size(), suffix) == 0;
}

struct GeneticCode {
    std::array<std::string, 64> codons;
    std::array<std::vector<int>, 64> synonyms;

    explicit GeneticCode(Rcpp::CharacterVector code) {
        const Rcpp::CharacterVector names = code.names();
        std::array<std::string, 64> amino_acids;
        std::vector<int> order;
        for (R_xlen_t i = 0; i < code.size(); ++i) {
            const std::string codon = Rcpp::as<std::string>(names[i]);
            const int index = codon_index(codon);
            codons[index] = codon;
            amino_acids[index] = Rcpp::as<std::string>(code[i]);
            order.push_back(index);
        }
        for (int original = 0; original < 64; ++original)
            for (int candidate : order)
                if (amino_acids[original] == amino_acids[candidate])
                    synonyms[original].push_back(candidate);
    }
};

// A state remembers the longest suffix that is a proper prefix of a forbidden motif.
// Once a motif is completed, the transition is rejected, including within a codon.
class MotifMatcher {
    std::vector<std::string> motifs;
    std::vector<std::string> prefixes = {""};

public:
    std::vector<std::array<int, 64>> next_codon;

    MotifMatcher(Rcpp::CharacterVector patterns, const GeneticCode &code) {
        for (R_xlen_t i = 0; i < patterns.size(); ++i) {
            const std::string motif = Rcpp::as<std::string>(patterns[i]);
            motifs.push_back(motif);
            for (std::size_t length = 1; length < motif.size(); ++length) {
                const std::string prefix = motif.substr(0, length);
                if (std::find(prefixes.begin(), prefixes.end(), prefix) == prefixes.end())
                    prefixes.push_back(prefix);
            }
        }
        next_codon.resize(prefixes.size());
        for (std::size_t state = 0; state < prefixes.size(); ++state)
            for (int codon = 0; codon < 64; ++codon)
                next_codon[state][codon] = advance(state, code.codons[codon]);
    }

    int advance(int state, const std::string &sequence) const {
        for (char base : sequence) {
            const std::string suffix = prefixes[state] + base;
            for (const std::string &motif : motifs)
                if (ends_with(suffix, motif))
                    return -1;
            state = 0;
            for (std::size_t i = 1; i < prefixes.size(); ++i)
                if (prefixes[i].size() > prefixes[state].size() &&
                    ends_with(suffix, prefixes[i]))
                    state = i;
        }
        return state;
    }
};

struct Trace {
    int previous = -1;
    int codon = -1;
};

std::string optimize(const std::string &sequence, const GeneticCode &code,
                     const MotifMatcher &matcher, R_xlen_t sequence_index) {
    const std::size_t n_codons = sequence.size() / 3;
    const std::size_t n_states = matcher.next_codon.size();
    const int unreachable = std::numeric_limits<int>::max();
    std::vector<int> costs(n_states, unreachable);
    std::vector<Trace> trace(n_codons * n_states);
    costs[0] = 0;

    for (std::size_t position = 0; position < n_codons; ++position) {
        if (position % 1024 == 0)
            Rcpp::checkUserInterrupt();
        const std::string original = sequence.substr(3 * position, 3);
        const std::vector<int> &synonyms = code.synonyms[codon_index(original)];
        std::array<int, 64> substitutions = {};
        for (int codon : synonyms)
            for (int base = 0; base < 3; ++base)
                substitutions[codon] += original[base] != code.codons[codon][base];

        std::vector<int> next_costs(n_states, unreachable);
        for (std::size_t state = 0; state < n_states; ++state) {
            if (costs[state] == unreachable)
                continue;
            for (int codon : synonyms) {
                const int next = matcher.next_codon[state][codon];
                if (next < 0)
                    continue;
                const int cost = costs[state] + substitutions[codon];
                if (cost < next_costs[next]) {
                    next_costs[next] = cost;
                    trace[position * n_states + next] = {static_cast<int>(state), codon};
                }
            }
        }
        costs.swap(next_costs);
    }

    // Preserve incomplete terminal codons, checking motifs across their boundary too.
    const std::string tail = sequence.substr(3 * n_codons);
    int best_state = -1;
    int best_cost = unreachable;
    for (std::size_t state = 0; state < n_states; ++state)
        if (costs[state] < best_cost && matcher.advance(state, tail) >= 0) {
            best_state = state;
            best_cost = costs[state];
        }
    if (best_state < 0)
        Rcpp::stop("No synonymous sequence avoids the cut sites (sequence %d)",
                   static_cast<int>(sequence_index + 1));

    std::string result = sequence;
    for (std::size_t position = n_codons; position > 0; --position) {
        const Trace &step = trace[(position - 1) * n_states + best_state];
        result.replace(3 * (position - 1), 3, code.codons[step.codon]);
        best_state = step.previous;
    }
    return result;
}

} // namespace

// [[Rcpp::export(rng = false)]]
Rcpp::CharacterVector remove_cutsite_nuc_cpp(Rcpp::CharacterVector nuc,
                                            Rcpp::CharacterVector sites,
                                            Rcpp::CharacterVector genetic_code) {
    const GeneticCode code(genetic_code);
    const MotifMatcher matcher(sites, code);
    Rcpp::CharacterVector result(nuc.size());
    for (R_xlen_t i = 0; i < nuc.size(); ++i) {
        Rcpp::checkUserInterrupt();
        result[i] = optimize(Rcpp::as<std::string>(nuc[i]), code, matcher, i);
    }
    return result;
}
