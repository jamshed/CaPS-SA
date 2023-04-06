
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
    const idx_t len_;   // Length of the input string.
    idx_t* const suf_arr_;  // The suffix array.
    idx_t* const lcp_arr_;  // The LCP array.


public:

    // Constructs a suffix array object for the input string `str` of size
    // `len`.
    Suffix_Array(const char* str, std::size_t len);

    ~Suffix_Array();
};

}



#endif
