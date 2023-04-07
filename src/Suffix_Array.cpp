
#include "Suffix_Array.hpp"

#include <cstdlib>
#include <algorithm>

namespace themis
{

Suffix_Array::Suffix_Array(const char* const str, const std::size_t n):
    str_(str),
    n(n),
    SA(static_cast<idx_t*>(std::malloc(n * sizeof(idx_t)))),
    LCP(static_cast<idx_t*>(std::malloc(n * sizeof(idx_t))))
{}


Suffix_Array::~Suffix_Array()
{
    std::free(SA);
    std::free(LCP);
}


void Suffix_Array::merge(const idx_t* X, idx_t len_x, const idx_t* Y, idx_t len_y, const idx_t* LCP_x, const idx_t* LCP_y, idx_t* Z, idx_t* LCP_z) const
{
    idx_t m = 0;    // LCP of the last compared pair.
    idx_t l_x;  // LCP(X_i, X_{i - 1}).
    idx_t i = 0;    // Index into `X`.
    idx_t j = 0;    // Index into `Y`.
    idx_t k = 0;    // Index into `Z`.

    while(i < len_x && j < len_y)
    {
        l_x = LCP_x[i];

        if(l_x > m)
            Z[k] = X[i],
            LCP_z[k] = l_x,
            m = m;
        else if(l_x < m)
            Z[k] = Y[j],
            LCP_z[k] = m,
            m = l_x;
        else    // Compute LCP of X_i and Y_j through linear scan.
        {
            idx_t n = m;    // LCP(X_i, Y_j).
            const idx_t max_n = n - std::max(X[i], Y[j]);   // Length of the shorter suffix.
            while(n < max_n && str_[X[i] + n] == str_[Y[j] + n])
                n++;

            // Whether the shorter suffix is a prefix of the longer one.
            Z[k] = (n == max_n ?    std::max(X[i], Y[j]) :
                                    (str_[X[i] + n] < str_[Y[j] + n] ? X[i] : Y[j]));
            LCP_z[k] = (Z[k] == X[i] ? l_x : m);
            m = n;
        }


        if(Z[k] == X[i])
            i++;
        else    // Swap X and Y (and associated data structures) when Y_j gets pulled into Z.
        {
            j++;
            std::swap(X, Y),
            std::swap(len_x, len_y),
            std::swap(LCP_x, LCP_y),
            std::swap(i, j);
        }

        k++;
    }


    for(; i < len_x; ++i, ++k)  // Copy rest of the data from X to Z.
        Z[k] = X[i], LCP_z[k] = LCP_x[i];

    if(j < len_y)   // Copy rest of the data from Y to Z.
    {
        Z[k] = Y[j], LCP_z[k] = m;
        for(j++, k++; j < len_y; ++j, ++k)
            Z[k] = Y[j], LCP_z[k] = LCP_y[j];
    }
}


Suffix_Array::~Suffix_Array()
{
    std::free(suf_arr_);
    std::free(lcp_arr_);
}

}
