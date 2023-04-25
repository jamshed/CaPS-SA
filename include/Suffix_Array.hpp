
#ifndef THEMIS_SUFFIX_ARRAY_HPP
#define THEMIS_SUFFIX_ARRAY_HPP



#include <cstdint>
#include <cstddef>
#include <atomic>
#include <cstdlib>
#include <chrono>
#include <immintrin.h>
#include <iostream>
// =============================================================================

namespace themis
{

// The Suffix Array (SA) and the Longest Common Prefix (LCP) array constructor
// class for some given sequence.
// TODO: templatize.
class Suffix_Array
{
private:

    typedef std::size_t idx_t;

    const char* const T_;   // The input text.
    const idx_t n_; // Length of the input text.
    idx_t* const SA_;   // The suffix array.
    idx_t* const LCP_;  // The LCP array.
    idx_t* SA_w;    // Working space for the SA construction.
    idx_t* LCP_w;   // Working space for the LCP construction.
    const std::size_t p_;   // Count of subproblems used in construction.
    idx_t* pivot_;  // Pivots for the global suffix array.
    const idx_t pivot_per_part_;    // Number of pivots to sample per subarray.
    idx_t* part_size_scan_; // Inclusive scan (prefix sum) of the sizes of the pivoted final partitions containing appropriate sorted sub-subarrays.
    idx_t* part_ruler_; // "Ruler" for the partitions—contains the indices of each sub-subarray in each partition.
    std::atomic_uint64_t solved_;   // Progress tracker—number of subproblems solved in some step.

    static constexpr std::size_t default_subproblem_count = 1024;   // Default subproblem-count to use in construction.

    // Fields for profiling time.
    typedef std::chrono::high_resolution_clock::time_point time_point_t;
    constexpr static auto now = std::chrono::high_resolution_clock::now;
    constexpr static auto duration = [](const std::chrono::nanoseconds& d) { return std::chrono::duration_cast<std::chrono::duration<double>>(d).count(); };


    // Returns the LCP length of `x` and `y`, where `min_len` is the length of
    // the shorter of `x` and `y`.
    static idx_t lcp(const char* x, const char* y, idx_t min_len);

    // Returns the LCP length of `x` and `y`, where `min_len` is the length of
    // the shorter of `x` and `y`. Optimized with some poor man's vectorization.
    static idx_t lcp_opt(const char* x, const char* y, idx_t min_len);

    // Returns the LCP length of `x` and `y`, where `min_len` is the length of
    // the shorter of `x` and `y`. Optimized with some poor man's vectorization.
    static idx_t lcp_opt_avx(const char* x, const char* y, idx_t min_len);


    // Merges the sorted collections of suffixes, `X` and `Y`, with lengths
    // `len_x` and `len_y` and LCP arrays `LCP_x` and `LCP_y` respectively, into
    // `Z`. Also constructs `Z`'s LCP array in `LCP_z`.
    void merge(const idx_t* X, idx_t len_x, const idx_t* Y, idx_t len_y, const idx_t* LCP_x, const idx_t* LCP_y, idx_t* Z, idx_t* LCP_z) const;

    // Merge-sorts the suffix collection `X` of length `n` into `Y`. Also
    // constructs the LCP array of `X` in `LCP`, using `W` as working space.
    // A necessary precondition is that `Y` must be equal to `X`.  `X` may
    // not remain the same after the sort.
    void merge_sort(idx_t* X, idx_t* Y, idx_t n, idx_t* LCP, idx_t* W) const;

    // Initializes internal data structures for the construction algorithm.
    void initialize();

    // Sorts uniform-sized subarrays independently.
    void sort_subarrays();

    // Samples `m` pivots from the sorted suffix collection `X` of size `n`
    // into `P`.
    static void sample_pivots(const idx_t* X, idx_t n, idx_t m, idx_t* P);

    // Selects pivots for parallel merging of the sorted subarrays.
    void select_pivots();

    // Locates the positions (upper-bounds) of the selected pivots in the sorted
    // subarrays and flattens them in `P`. Besides these pivots, two flanking
    // pivots, `0` and `|X_i|`, for each subarray `X_i` are also placed.
    void locate_pivots(idx_t* P) const;

    // Returns the first index `idx` into the sorted suffix collection `X` of
    // length `n` such that `X[idx]` is strictly greater than the query pattern
    // `P` of length `P_len`.
    std::size_t upper_bound(const idx_t* X, idx_t n, const char* P, std::size_t P_len) const;

    // Collates the sub-subarrays delineated by the pivot locations in each
    // sorted subarray, present in `P`, into appropriate partitions.
    void partition_sub_subarrays(const idx_t* P);

    // Merges the sorted sub-subarrays laid flat together in each partition.
    void merge_sub_subarrays();

    // Computes the LCPs at the partition boundaries, specifically at the
    // starting index of each partition in their flat collection.
    void compute_partition_boundary_lcp();

    // Merge-sorts the collection `X` that contains `n` sorted arrays of
    // suffixes laid flat together, into `Y`. `S` contains the delimiter indices
    // of the `n` arrays in `X`. The LCP array of sorted `X` is constructed in
    // `LCP_y`; `LCP_x` contains the LCP arrays of the `n` arrays of `X`.
    // A necessary precondition is that `Y` must be equal to `X`, and `LCP_y` to
    // `LCP_x`. `X` and `LCP_x` may not remain the same after the sort.
    void sort_partition(idx_t* X, idx_t* Y, idx_t n, const idx_t* S, idx_t* LCP_x, idx_t* LCP_y);

    // Cleans up after the construction algorithm.
    void clean_up();

    // Returns pointer to a memory-allocation for `size` elements of type `T_`.
    template <typename T_>
    static T_* allocate(std::size_t size) { return static_cast<T_*>(std::malloc(size * sizeof(T_))); }

    // Returns true iff `X` is a valid (partial) suffix array with size `n`.
    bool is_sorted(const idx_t* X, idx_t n) const;

public:

    // Constructs a suffix array object for the input text `T` of size
    // `n`. Optionally, the number of subproblems to decompose the original
    // construction problem into can be provided with `subproblem_count`.
    Suffix_Array(const char* T, std::size_t n, std::size_t subproblem_count = 0);

    // Copy constructs the suffix array object from `other`.
    Suffix_Array(const Suffix_Array& other);

    ~Suffix_Array();

    const Suffix_Array& operator=(const Suffix_Array& rhs) = delete;

    // Returns the text.
    const char* T() const { return T_; }

    // Returns the length of the text.
    std::size_t n() const { return n_; }

    // Returns the suffix array.
    const idx_t* SA() const { return SA_; }

    // Returns the LCP array.
    const idx_t* LCP() const { return LCP_; }

    // Constructs the suffix array and the LCP array.
    void construct();
};


inline Suffix_Array::idx_t Suffix_Array::lcp(const char* const x, const char* const y, const idx_t min_len)
{
    idx_t l = 0;
    while(l < min_len && x[l] == y[l])
        l++;

    return l;
}

inline Suffix_Array::idx_t Suffix_Array::lcp_opt_avx(const char* str1, const char* str2, const idx_t len_in) {
  int64_t i = 0;
  int64_t len = static_cast<int64_t>(len_in);
  if (len >= 32) {
    for (; i <= len - 32; i += 32) {
      __m256i v1 = _mm256_loadu_si256((__m256i*)(str1 + i));
      __m256i v2 = _mm256_loadu_si256((__m256i*)(str2 + i));
      __m256i cmp = _mm256_cmpeq_epi8(v1, v2);
      int mask = _mm256_movemask_epi8(cmp);
      if (mask != 0xFFFFFFFF) {
        int j = __builtin_ctz(~mask) + i;
        return static_cast<idx_t>(j);
      }
    }
  }
  for (; i < len; i++) {
    if (str1[i] != str2[i]) {
      break;
    }
  }
  return static_cast<idx_t>(i);
}
/*
    const size_t len = min_len;

    const __m256i* p1 = reinterpret_cast<const __m256i*>(str1);
    const __m256i* p2 = reinterpret_cast<const __m256i*>(str2);
    const __m256i* end = reinterpret_cast<const __m256i*>(str1 + (len & ~0x1F)); // round down to nearest multiple of 32 bytes

    for (; p1 < end; p1++, p2++) {
        __m256i v1 = _mm256_loadu_si256(p1);
        __m256i v2 = _mm256_loadu_si256(p2);
        __m256i cmp = _mm256_cmpeq_epi8(v1, v2);
        int mask = _mm256_movemask_epi8(cmp);
        if (mask != 0xFFFFFFFF) {
            // found a mismatch, determine the position of the first differing byte
            uint32_t pos = __builtin_ctz(~mask);
            return (reinterpret_cast<const char*>(p1) - str1) + pos;
        }
    }

    // handle the remaining bytes (less than 32)
    for (; reinterpret_cast<const char*>(p1) < str1 + len; p1++, p2++) {
        if (*reinterpret_cast<const char*>(p1) != *reinterpret_cast<const char*>(p2)) {
            return reinterpret_cast<const char*>(p1) - str1;
        }
    }

    // strings are identical up to the length of the shorter string
    return len;
}
*/

inline Suffix_Array::idx_t Suffix_Array::lcp_opt(const char* const x, const char* const y, const idx_t min_len)
{
    auto const X = reinterpret_cast<const uint64_t*>(x);
    auto const Y = reinterpret_cast<const uint64_t*>(y);
    const auto word_count = (min_len >> 3);

    idx_t i = 0;
    while(i < word_count && X[i] == Y[i])
        i++;

    return (i << 3) + lcp(x + (i << 3), y + (i << 3), min_len - (i << 3));
}

}



#endif
