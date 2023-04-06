
#include "Suffix_Array.hpp"

#include <cstdlib>
#include <algorithm>

namespace themis
{

Suffix_Array::Suffix_Array(const char* const str, const std::size_t len):
    str_(str),
    len_(len),
    suf_arr_(static_cast<idx_t*>(std::malloc(len * sizeof(idx_t)))),
    lcp_arr_(static_cast<idx_t*>(std::malloc(len * sizeof(idx_t))))
{}


Suffix_Array::~Suffix_Array()
{
    std::free(suf_arr_);
    std::free(lcp_arr_);
}

}
