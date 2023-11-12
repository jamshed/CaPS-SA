
#ifndef THEMIS_SUFFIX_ARRAY_HPP
#define THEMIS_SUFFIX_ARRAY_HPP



#include "Bit_Packed_Text.hpp"
#include "utilities.hpp"

#include <cstdint>
#include <cstddef>
#include <atomic>
#include <chrono>
#include <immintrin.h>
#include <iostream>
#include <vector>


// =============================================================================

namespace CaPS_SA
{

struct PrefixLookupTab {
  static constexpr uint16_t context_len = 8;
  std::vector<uint64_t> entries;
  std::vector<int64_t> breakpoints;
  bool empty = true;
  uint16_t last = 0;
  int64_t max_breakpoint = 0;

  PrefixLookupTab() {
    entries.resize(65536, std::numeric_limits<uint64_t>::max());
    breakpoints.push_back(0);
  }

  inline bool insert(uint16_t u, uint64_t offset) {
    if ((u > last) or empty) {
      entries[u] = breakpoints.size();
      breakpoints.push_back(offset);
      last = u;
      empty = false;
      return true;
    }
    return false;
  }
  
  void finish(uint64_t n) { 
    breakpoints.push_back(n); 
    max_breakpoint = n;
  }

  void fill() {
    uint64_t prev = 0;
    for (size_t i = 0; i < entries.size(); ++i) {
      // fill in sentinel values with the previous value
      if (entries[i] == std::numeric_limits<uint64_t>::max()) {
        entries[i] = prev;
      } else {
        prev = entries[i];
      }
    }
  }

  std::pair<uint64_t, uint64_t> get(uint16_t u) const {
    auto i = entries[u];
    return std::make_pair(breakpoints[i], breakpoints[i+1]);
  }

  inline std::pair<int64_t, int64_t> get_expanded(uint16_t u) const {
    auto i = entries[u];
    auto start = breakpoints[i]; --start; 
    auto stop = breakpoints[i+1];
    stop = std::min(stop + 1, max_breakpoint);//(stop <  breakpoints.back()) ? stop + 1 : breakpoints.back();
    return std::make_pair(start, stop);
  }

  size_t size() const { return entries.size(); }
};

// The Suffix Array (SA) and the Longest Common Prefix (LCP) array constructor
// class for some given sequence.
template <typename T_idx_>
class Suffix_Array
{
private:

    typedef T_idx_ idx_t;   // Integer-type for indexing the input text.

    const char* const T_;   // The input text.
    const idx_t n_; // Length of the input text.
    Bit_Packed_Text B;  // Bit-packed representation of the text.
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
    std::atomic_uint64_t solved_;   // Progress tracker—number of subproblems solved in some step.

    static constexpr idx_t default_subproblem_count = 8192; // Default subproblem-count to use in construction.
    static constexpr idx_t nested_par_grain_size = (1lu << 13); // Granularity for nested parallelism to kick in.

    // Fields for profiling time.
    typedef std::chrono::high_resolution_clock::time_point time_point_t;
    constexpr static auto now = std::chrono::high_resolution_clock::now;
    constexpr static auto duration = [](const std::chrono::nanoseconds& d) { return std::chrono::duration_cast<std::chrono::duration<double>>(d).count(); };


    // Returns the LCP length of `x` and `y`, where `min_len` is the length of
    // the shorter of `x` and `y`.
    static idx_t lcp(const char* x, const char* y, idx_t min_len);

    // Returns the LCP length of `x` and `y`, where `ctx` is the length of their
    // context.
    idx_t lcp(idx_t x, idx_t y, idx_t ctx) const;

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

    // Merges the sorted collections of suffixes, `X` and `Y`, with lengths
    // `len_x` and `len_y` and LCP arrays `LCP_x` and `LCP_y` respectively, into
    // `Z`. Also constructs `Z`'s LCP array in `LCP_z`.
    void merge(const idx_t* X, idx_t len_x, const idx_t* Y, idx_t len_y, const idx_t* LCP_x, const idx_t* LCP_y, idx_t* Z, idx_t* LCP_z) const;

    // Merge-sorts the suffix collection `X` of length `n` into `Y`. Also
    // constructs the LCP array of `Y` in `LCP`, using `W` as working space.
    // A necessary precondition is that `Y` must be equal to `X`.  `X` may
    // not remain the same after the sort.
    void merge_sort(idx_t* X, idx_t* Y, idx_t n, idx_t* LCP, idx_t* W) const;

    // Insertion-sorts the suffix collection `X` of length `n` into `Y`. Also
    // constructs the LCP-array of `Y` in `LCP`.
    void insertion_sort(idx_t* X, idx_t* Y, idx_t n, idx_t* LCP) const;

    // Returns `true` iff the suffix `x` is lesser than the suffix `y`.
    bool suf_less(idx_t x, idx_t y) const;

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
    idx_t upper_bound(const idx_t* X, idx_t n, const char* P, idx_t P_len) const;
    idx_t upper_bound_with_lookup(const idx_t* X, idx_t n, const char* P, idx_t P_len, 
                                  const PrefixLookupTab& lookup) const;


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


public:

    // Constructs a suffix array object for the input text `T` of size
    // `n`. Optionally, the number of subproblems to decompose the original
    // construction problem into can be provided with `subproblem_count`, and
    // the maximum prefix-context length for the suffixes can be bounded by
    // `max_context`.
    Suffix_Array(const char* T, idx_t n, idx_t subproblem_count = 0, idx_t max_context = 0);

    // Copy constructs the suffix array object from `other`.
    Suffix_Array(const Suffix_Array& other);

    ~Suffix_Array();

    const Suffix_Array& operator=(const Suffix_Array& rhs) = delete;

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

    // Returns true iff `X` is a valid (partial) suffix array with size `n`.
    // TODO: add LCP-array check.
    bool is_sorted(const idx_t* X, idx_t n) const;
};


template <typename T_idx_>
inline bool Suffix_Array<T_idx_>::suf_less(const idx_t x, const idx_t y) const
{
    const idx_t max_n = n_ - std::max(x, y);  // Length of the shorter suffix.
    const idx_t context = std::min(max_context, max_n); // Prefix-context length for the suffixes.
    const idx_t lcp = lcp_opt_avx_unrolled(T_ + x, T_ + y, context);

    if(lcp == max_n)
        return x > y;

    return T_[x + lcp] < T_[y + lcp];
}


template <typename T_idx_>
inline T_idx_ Suffix_Array<T_idx_>::lcp(const char* const x, const char* const y, const idx_t min_len)
{
    idx_t l = 0;
    while(l < min_len && x[l] == y[l])
        l++;

    return l;
}


template <typename T_idx_>
inline T_idx_ Suffix_Array<T_idx_>::lcp(const idx_t x, const idx_t y, const idx_t ctx) const
{
    const auto v_x = *reinterpret_cast<const uint64_t*>(T_ + x);
    const auto v_y = *reinterpret_cast<const uint64_t*>(T_ + y);
    if(v_x != v_y)
        return __builtin_ctzll(v_x ^ v_y) >> 3;

    return B.LCP(x, y, ctx);
}


#define LCPCMP(N, IDX_T) \
      __m256i v1 ## N = _mm256_loadu_si256((__m256i*)(str1 + i + N));\
      __m256i v2 ## N = _mm256_loadu_si256((__m256i*)(str2 + i + N));\
      __m256i cmp ## N = _mm256_cmpeq_epi8(v1##N, v2##N);\
      int mask ## N = _mm256_movemask_epi8(cmp##N);\
      if (mask ## N != static_cast<int>(0xFFFFFFFF)) {\
        int j = __builtin_ctz(~mask ## N) + i + N;\
        return static_cast<IDX_T>(j);\
      } 

#define M_REPEAT_1(X, T) X(0, T)
#define M_REPEAT_2(X, T) X(0, T) X(32, T)
#define M_REPEAT_3(X, T) M_REPEAT_2(X, T) X(64, T)
#define M_REPEAT_4(X, T) M_REPEAT_3(X, T) X(96, T)
#define M_REPEAT_5(X, T) M_REPEAT_4(X, T) X(128, T)
#define M_REPEAT_6(X, T) M_REPEAT_5(X, T) X(160, T)
#define M_REPEAT_7(X, T) M_REPEAT_6(X, T) X(192, T)
#define M_REPEAT_8(X, T) M_REPEAT_7(X, T) X(224, T)

template <typename T_idx_>
inline T_idx_ Suffix_Array<T_idx_>::lcp_opt_avx_unrolled(const char* str1, const char* str2, const idx_t len_in) {
  int64_t i = 0;
  int64_t len = static_cast<int64_t>(len_in);

  if ((len - i)>= 160) {
    for (; i <= len - 160; i += 160) {
      M_REPEAT_5(LCPCMP, idx_t);
    } 
  }

  if ((len - i)>= 128) {
    for (; i <= len - 128; i += 128) {
      M_REPEAT_4(LCPCMP, idx_t);
    } 
  }

  if ((len - i)>= 96) {
    for (; i <= len - 96; i += 96) {
      M_REPEAT_3(LCPCMP, idx_t);
    }
  }

  if ((len - i)>= 64) {
    for (; i <= len - 64; i += 64) {
      M_REPEAT_2(LCPCMP, idx_t);
    }
  }

  if ((len - i) >= 32) {
    for (; i <= len - 32; i += 32) {
      __m256i v1 = _mm256_loadu_si256((__m256i*)(str1 + i));
      __m256i v2 = _mm256_loadu_si256((__m256i*)(str2 + i));
      __m256i cmp = _mm256_cmpeq_epi8(v1, v2);
      int mask = _mm256_movemask_epi8(cmp);
      if (mask != static_cast<int>(0xFFFFFFFF)) {
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

template <typename T_idx_>
inline T_idx_ Suffix_Array<T_idx_>::lcp_opt_avx(const char* str1, const char* str2, const idx_t len_in) {
  int64_t i = 0;
  int64_t len = static_cast<int64_t>(len_in);
  if (len >= 32) {
    for (; i <= len - 32; i += 32) {
      __m256i v1 = _mm256_loadu_si256((__m256i*)(str1 + i));
      __m256i v2 = _mm256_loadu_si256((__m256i*)(str2 + i));
      __m256i cmp = _mm256_cmpeq_epi8(v1, v2);
      int mask = _mm256_movemask_epi8(cmp);
      if (mask != static_cast<int>(0xFFFFFFFF)) {
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


template <typename T_idx_>
inline T_idx_ Suffix_Array<T_idx_>::lcp_opt(const char* const x, const char* const y, const idx_t min_len)
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
