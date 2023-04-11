
#ifndef THEMIS_SUFFIX_ARRAY_HPP
#define THEMIS_SUFFIX_ARRAY_HPP



#include <cstddef>
#include <cstdlib>
#include <chrono>


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

    const char* const str_; // The input string.
    const idx_t n_; // Length of the input string.
    idx_t* const SA_;   // The suffix array.
    idx_t* const LCP_;  // The LCP array.
    idx_t* SA_w;    // Working space for the SA construction.
    idx_t* LCP_w;   // Working space for the LCP construction.
    const std::size_t p_;   // Count of threads used in construction.
    idx_t* P;   // Pivots for the global suffix array.
    const idx_t pivot_per_part; // Number of pivots to sample per subarray.

    // Fields for profiling time.
    typedef std::chrono::high_resolution_clock::time_point time_point_t;
    constexpr static auto now = std::chrono::high_resolution_clock::now;
    constexpr static auto duration = [](const std::chrono::nanoseconds& d) { return std::chrono::duration_cast<std::chrono::duration<double>>(d).count(); };


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

    // Cleans up after the construction algorithm.
    void clean_up();

    // Returns pointer to a memory-allocation for `size` elements of type `T_`.
    template <typename T_>
    static T_* allocate(std::size_t size) { return static_cast<T_*>(std::malloc(size * sizeof(T_))); }

public:

    // Constructs a suffix array object for the input string `str` of size
    // `n`.
    Suffix_Array(const char* str, std::size_t n);

    // Copy constructs the suffix array object from `other`.
    Suffix_Array(const Suffix_Array& other);

    ~Suffix_Array();

    const Suffix_Array& operator=(const Suffix_Array& rhs) = delete;

    // Returns the length of the string.
    std::size_t n() const { return n_; }

    // Returns the suffix array.
    const idx_t* SA() const { return SA_; }

    // Returns the LCP array.
    const idx_t* LCP() const { return LCP_; }

    // Constructs the suffix array and the LCP array.
    void construct();
};

}



#endif
