
#include "Suffix_Array.hpp"
#include "parlay/parallel.h"

#include <cstdlib>
#include <cstring>
#include <iostream>
#include <cstdlib>
#include <algorithm>
#include <cassert>

namespace themis
{

Suffix_Array::Suffix_Array(const char* const str, const std::size_t n):
    str_(str),
    n_(n),
    SA_(static_cast<idx_t*>(std::malloc(n * sizeof(idx_t)))),
    LCP_(static_cast<idx_t*>(std::malloc(n * sizeof(idx_t)))),
    SA_w(nullptr),
    LCP_w(nullptr),
    p_(std::getenv("PARLAY_NUM_THREADS") == nullptr ? 0 : std::atoi(std::getenv("PARLAY_NUM_THREADS")))
{
    if(p_ == 0)
    {
        std::cerr << "The environment variable `PARLAY_NUM_THREADS` needs to be set. Aborting.\n";
        std::exit(EXIT_FAILURE);
    }
}


Suffix_Array::Suffix_Array(const Suffix_Array& other): Suffix_Array(other.str_, other.n_)
{
    std::memcpy(SA_, other.SA_, n_ * sizeof(idx_t));
    std::memcpy(LCP_, other.LCP_, n_ * sizeof(idx_t));
}


Suffix_Array::~Suffix_Array()
{
    std::free(SA_);
    std::free(LCP_);
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
            const idx_t max_n = n_ - std::max(X[i], Y[j]);  // Length of the shorter suffix.
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


void Suffix_Array::merge_sort(idx_t* const X, idx_t* const Y, const idx_t n, idx_t* const LCP, idx_t* const W) const
{
    assert(std::memcmp(X, Y, n * sizeof(idx_t)) == 0);

    if(n == 1)
        LCP[0] = 0;
    else
    {
        const idx_t m = n / 2;
        merge_sort(Y, X, m, W, LCP);
        merge_sort(Y + m, X + m, n - m, W + m, LCP + m);
        merge(X, m, X + m, n - m, W, W + m, Y, LCP);
    }
}


void Suffix_Array::initialize()
{
    SA_w = static_cast<idx_t*>(std::malloc(n_ * sizeof(idx_t)));    // Working space for the SA construction.
    LCP_w = static_cast<idx_t*>(std::malloc(n_ * sizeof(idx_t)));   // Working space for the LCP construction.

    const auto idx_init = [SA_ = SA_, SA_w = SA_w](const std::size_t i){ SA_[i] = SA_w[i] = i; };
    parlay::parallel_for(0, n_, idx_init);
}


void Suffix_Array::clean_up()
{
    std::free(SA_w);
    std::free(LCP_w);
}


void Suffix_Array::construct()
{
    const auto t_start = now();

    initialize();
    merge_sort(SA_w, SA_, n_, LCP_, LCP_w);

    clean_up();

    const auto t_end = now();
    std::cerr << "Constructed the suffix array. Time taken: " << duration(t_end - t_start) << " seconds.\n";
}

}
