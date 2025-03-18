
#ifndef THEMIS_SUFFIX_ARRAY_HPP
#define THEMIS_SUFFIX_ARRAY_HPP



#include "x86/avx2.h"

#include <cstdint>
#include <cstddef>
#include <cstdlib>
#include <fstream>
#include <chrono>

#include <immintrin.h>

// =============================================================================

namespace CaPS_SA
{

// The Suffix Array (SA) and the Longest Common Prefix (LCP) array constructor
// class for some given sequence.
template <typename T_idx_>
class Suffix_Array
{
private:

    typedef T_idx_ idx_t;   // Integer-type for indexing the input text.

    const char* const T_;   // The input text.
    const idx_t n_; // Length of the input text.
    idx_t* const SA_;   // The suffix array.
    idx_t* const LCP_;  // The LCP array.
    idx_t* SA_w;    // Working space for the SA construction.
    idx_t* LCP_w;   // Working space for the LCP construction.
    const idx_t p_; // Count of subproblems used in construction.
    const idx_t max_context;    // Maximum prefix-context length for comparing suffixes.
    idx_t* pivot_;  // Pivots for the global suffix array.
    const idx_t pivot_per_part_;    // Number of pivots to sample per subarray.
    idx_t* part_size_scan_; // Inclusive scan (prefix sum) of the sizes of the pivoted final partitions containing appropriate sorted sub-subarrays.
    idx_t* part_ruler_; // "Ruler" for the partitions—contains the indices of each sub-subarray in each partition.

    static constexpr idx_t default_subproblem_count = 8192; // Default subproblem-count to use in construction.
    static constexpr idx_t nested_par_grain_size = (1lu << 13); // Granularity for nested parallelism to kick in.

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

    // Returns the LCP length of `x` and `y`, where `min_len` is the length of
    // the shorter of `x` and `y`. Optimized with some poor man's vectorization.
    // NOTE: hand unrolled version of `lcp_opt_avx`.
    static idx_t lcp_opt_avx_unrolled(const char* x, const char* y, idx_t min_len);

    // Returns the LCP length of `x` and `y`, where `min_len` is the length of
    // the shorter of `x` and `y`. `N x 32` bytes of prefix comparisons are
    // loop-unrolled.
    template <std::size_t N = 8>
    static idx_t LCP(const char* x, const char* y, idx_t min_len);

    // Returns the LCP length of the `32 x N`-bytes prefix of `x` and `y`.
    template <std::size_t N> static idx_t LCP_unrolled(const char* x, const char* y);

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

    // Populates the suffix array with some permutation of `[0, len)`.
    void permute();

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
    idx_t upper_bound(const idx_t* X, idx_t n, const char* P, idx_t P_len) const;

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
    static T_* allocate(idx_t size) { return static_cast<T_*>(std::malloc(size * sizeof(T_))); }

    // Deallocates the pointer `ptr`, allocated with `allocate`.
    template <typename T_>
    static void deallocate(T_* const ptr) { std::free(ptr); }

    // Returns true iff `X` is a valid (partial) suffix array with size `n`, and
    // `L` is its LCP-array.
    bool is_sorted(const idx_t* X, const idx_t* L, idx_t n) const;

public:

    // Constructs a suffix array object for the input text `T` of size
    // `n`. Optionally, the number of subproblems to decompose the original
    // construction problem into can be provided with `subproblem_count`, and
    // the maximum prefix-context length for the suffixes can be bounded by
    // `max_context`.
    Suffix_Array(const char* T, idx_t n, idx_t subproblem_count = 0, idx_t max_context = 0);

    Suffix_Array(const Suffix_Array&) = delete;
    Suffix_Array& operator=(const Suffix_Array&) = delete;
    Suffix_Array(Suffix_Array&&) = delete;
    Suffix_Array& operator=(Suffix_Array&&) = delete;

    ~Suffix_Array();

    // Returns the text.
    const char* T() const { return T_; }

    // Returns the length of the text.
    idx_t n() const { return n_; }

    // Returns the suffix array.
    const idx_t* SA() const { return SA_; }

    // Returns the LCP array.
    const idx_t* LCP() const { return LCP_; }

    // Constructs the suffix array and the LCP array.
    void construct();

    // Dumps the suffix array and the LCP array into the stream `output`.
    void dump(std::ofstream& output);
};


template <typename T_idx_>
inline T_idx_ Suffix_Array<T_idx_>::lcp(const char* const x, const char* const y, const idx_t min_len)
{
    idx_t l = 0;
    while(l < min_len && x[l] == y[l])
        l++;

    return l;
}


template <typename T_idx_>
template <std::size_t N>
inline T_idx_ Suffix_Array<T_idx_>::LCP(const char* const x, const char* const y, const idx_t min_len)
{
    idx_t lcp = 0;

    if constexpr(N == 1)
    {
        for(; lcp < min_len; ++lcp)
            if(x[lcp] != y[lcp])
                break;

        return lcp;
    }
    else
    {
        while((min_len - lcp) >= N * 32)
        {
            const auto l = LCP_unrolled<N>(x + lcp, y + lcp);
            lcp += l;
            if(l < N * 32)
                return lcp;
        }

        return lcp + LCP<N - 1>(x + lcp, y + lcp, min_len - lcp);
    }
}


template <typename T_idx_>
template <std::size_t N>
inline T_idx_ Suffix_Array<T_idx_>::LCP_unrolled(const char* const x, const char* const y)
{
    if constexpr(N == 0)
        return 0;
    else
    {
        const auto v1 = simde_mm256_loadu_si256(reinterpret_cast<const simde__m256i*>(x));
        const auto v2 = simde_mm256_loadu_si256(reinterpret_cast<const simde__m256i*>(y));
        const auto cmp = simde_mm256_cmpeq_epi8(v1, v2);
        const auto mask = static_cast<uint32_t>(simde_mm256_movemask_epi8(cmp));
        if(mask != 0xFFFFFFFF)
            return __builtin_ctz(~mask);

        return 32 + LCP_unrolled<N - 1>(x + 32, y + 32);
    }
}

}



#endif
