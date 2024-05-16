#include <utility>
#include <algorithm>

#include <vector>
#include <mutex>

#include <immintrin.h>

#include <parlay/sequence.h>

template <typename index_t = size_t>
struct memoization_table {
    using prefix_t = std::pair<index_t, index_t>;

    parlay::sequence<parlay::sequence<prefix_t>> diagonals;
    std::vector<std::mutex> locks;

    memoization_table(index_t n) {
        diagonals = parlay::sequence<parlay::sequence<prefix_t>>(n);
        locks = std::vector<std::mutex>(n);
    }

    index_t manual_lcp(char *str, index_t len, index_t a, index_t b, index_t n) {
        if (a == b) {
            return len - a;
        }

        index_t lcp = 0;
        for (index_t i = 0; i < n; i += 32) {
            __m256i suff_a = _mm256_loadu_si256((__m256i*)&str[a + i]);
            __m256i suff_b = _mm256_loadu_si256((__m256i*)&str[b + i]);

            __m256i cmp = _mm256_cmpeq_epi8(suff_a, suff_b);
            uint32_t mask = _mm256_movemask_epi8(cmp);

            if (mask != 0xFFFFFFFF) {
                return std::min(n, lcp + __builtin_ctz(~mask));
            }
            lcp += 32;
        }
        return std::min(n, lcp);
    }

    index_t get_lcp(char *str, size_t len, index_t a, index_t b) {
        index_t early_lcp = manual_lcp(str, len, a, b, 32);
        if (early_lcp < 32) return early_lcp;

        if (a == b) {
            return len - a;
        }

        // TODO: start with one step of naive string comparison
        // TODO: (possible) save the results of this naive string comparisons and start from that point after consulting mem table
        if (a > b) {
            std::swap(a, b);
        }

        index_t diagonal_index = b - a;
        if (diagonals[diagonal_index].size() > 1) {
            int extra_cmps = 31 - __builtin_clz(diagonals[diagonal_index].size());
            int extra_lcp = manual_lcp(str, len, a + 32, b + 32, 32 * extra_cmps);
            if (extra_lcp < 32 * extra_cmps) return early_lcp + extra_lcp;
            else early_lcp += extra_lcp;
        }

        locks[diagonal_index].lock();

        // Find the first prefix in the diagonal that comes after (a, b)
        auto i = std::find_if(diagonals[diagonal_index].begin(), diagonals[diagonal_index].end(), [&] (prefix_t& prefix) -> bool {
                return prefix.first > a;
            });
        index_t index_in_diagonal = i - diagonals[diagonal_index].begin();

        // If there is a previous LCP that covers our queried LCP
        if (index_in_diagonal > 0 && a - diagonals[diagonal_index][index_in_diagonal - 1].first <= diagonals[diagonal_index][index_in_diagonal - 1].second) {
            locks[diagonal_index].unlock();
            return diagonals[diagonal_index][index_in_diagonal - 1].second - (a - diagonals[diagonal_index][index_in_diagonal - 1].first);
        }

        // TODO: (possible) when inserting a new LCP, replace old LCPs of length <= 64 * log(diagonal.size())
        index_t limit = (index_in_diagonal == diagonals[diagonal_index].size() ? len - a : diagonals[diagonal_index][index_in_diagonal].first - a);
        index_t next_lcp_len = (index_in_diagonal == diagonals[diagonal_index].size() ? 0 : diagonals[diagonal_index][index_in_diagonal].second);
        index_t lcp = early_lcp;
        if (lcp < limit) {
            lcp += manual_lcp(str, len, a + early_lcp, b + early_lcp, limit - early_lcp);
        }

        if (lcp >= limit) {
            diagonals[diagonal_index][index_in_diagonal].first = a;
            diagonals[diagonal_index][index_in_diagonal].second = limit + next_lcp_len;

            locks[diagonal_index].unlock();
            return diagonals[diagonal_index][index_in_diagonal].second;
        }
        
        if (lcp > 32) {
            diagonals[diagonal_index].resize(diagonals[diagonal_index].size() + 1);
            std::move_backward(diagonals[diagonal_index].begin() + index_in_diagonal, diagonals[diagonal_index].end() - 1, diagonals[diagonal_index].end());
            diagonals[diagonal_index][index_in_diagonal].first = a;
            diagonals[diagonal_index][index_in_diagonal].second = lcp;
        }

        locks[diagonal_index].unlock();
        return lcp;
    }
};
