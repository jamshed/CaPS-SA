
#ifndef THEMIS_SUFFIX_ARRAY_HPP
#define THEMIS_SUFFIX_ARRAY_HPP



#include <cstddef>


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
    const idx_t n;  // Length of the input string.
    idx_t* const SA;  // The suffix array.
    idx_t* const LCP;   // The LCP array.


    // Merges the sorted collections of suffixes, `X` and `Y`, with lengths
    // `len_x` and `len_y` and LCP arrays `LCP_x` and `LCP_y` respectively, into
    // `Z`. Also constructs `Z`'s LCP array in `LCP_z`.
    void merge(const idx_t* X, idx_t len_x, const idx_t* Y, idx_t len_y, const idx_t* LCP_x, const idx_t* LCP_y, idx_t* Z, idx_t* LCP_z) const;

    // Merge-sorts the suffix collection `X` of length `n` into `Y`. Also
    // constructs the LCP array of `X` in `LCP`, using `W` as working space.
    // A necessary precondition is that `Y` must be equal to `X`.  `X` may
    // not remain the same after the sort.
    void merge_sort(idx_t* X, idx_t* Y, idx_t n, idx_t* LCP, idx_t* W) const;


public:

    // Constructs a suffix array object for the input string `str` of size
    // `n`.
    Suffix_Array(const char* str, std::size_t n);

    // Copy constructs the suffix array object from `other`.
    Suffix_Array(const Suffix_Array& other);

    ~Suffix_Array();

    const Suffix_Array& operator=(const Suffix_Array& rhs) = delete;

    // Constructs the suffix array and the LCP array.
    void construct();
};

}



#endif
