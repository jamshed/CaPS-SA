#include <utility>
#include <algorithm>

#include <vector>
#include <mutex>

#include <parlay/sequence.h>

#define BITS_PER_CHAR 2

template <typename index_t = size_t>
struct memoization_table {
    using prefix_t = std::pair<index_t, index_t>;

    parlay::sequence<parlay::sequence<prefix_t>> diagonals;
    std::vector<std::mutex> locks;

    memoization_table(index_t n) {
        diagonals = parlay::sequence<parlay::sequence<prefix_t>>(n);
        locks = std::vector<std::mutex>(n);
    }

    index_t get_lcp(char *str, size_t len, index_t a, index_t b) {
        if (a == b) {
            return len - a;
        }

        // TODO: start with one step of naive string comparison
        // TODO: (possible) save the results of this naive string comparisons and start from that point after consulting mem table
        if (a > b) {
            std::swap(a, b);
        }

        index_t diagonal_index = b - a;
        locks[diagonal_index].lock();

        // Find the first prefix in the diagonal that comes after (a, b)
        auto i = find_if(diagonals[diagonal_index], [&] (prefix_t& prefix) -> bool {
                return prefix.first > a;
            });
        index_t index_in_diagonal = i - diagonals[diagonal_index].begin();

        // If there is a previous LCP that covers our queried LCP
        if (index_in_diagonal > 0 && a - diagonals[diagonal_index][index_in_diagonal - 1].first <= diagonals[diagonal_index][index_in_diagonal - 1].second) {
            locks[diagonal_index].unlock();
            return diagonals[diagonal_index][index_in_diagonal - 1].second - (a - diagonals[diagonal_index][index_in_diagonal - 1]);
        }

        // TODO: (possible) when inserting a new LCP, replace old LCPs of length <= 64 * log(diagonal.size())
        index_t limit = (index_in_diagonal == diagonals[diagonal_index].size() ? len - a : diagonals[diagonal_index][index_in_diagonal].first - a);
        index_t next_lcp_len = (index_in_diagonal == diagonals[diagonal_index].size() ? 0 : diagonals[diagonal_index][index_in_diagonal].second);
        for (int j = 0; j < limit; j++) {
            if (str[a + j] != str[b + j]) {
                // Reached the end of the current LCP
                if (j > 64) {
                    diagonals[diagonal_index].push_back((index_t)0);
                    std::memmove(&diagonals[diagonal_index][index_in_diag], &diagonals[diagonal_index][index_in_diag + 1], (diagonals[diagonal_index].size() - index_in_diag + 1) * sizeof(prefix_t));
                    diagonals[diagonal_index][index_in_diag].first = a;
                    diagonals[diagonal_index][index_in_diag].second = j;
                }
                
                locks[diagonal_index].unlock();
                return j;
            }
        }

        // Reached the start of the next LCP
        diagonals[diagonal_index][index_in_diag].first = a;
        diagonals[diagonal_index][index_in_diag].second = limit + next_lcp_len;

        locks[diagonal_index].unlock();
        return diagonals[diagonal_index][index_in_diag].second;
    }
};
