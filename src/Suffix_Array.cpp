
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
    SA_(allocate<idx_t>(n_)),
    LCP_(allocate<idx_t>(n_)),
    SA_w(nullptr),
    LCP_w(nullptr),
    p_(std::getenv("PARLAY_NUM_THREADS") == nullptr ? 0 : std::atoi(std::getenv("PARLAY_NUM_THREADS"))),
    P(nullptr),
    pivot_per_part(p_ - 1)
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
            const idx_t max_n = n_ - std::max(X[i], Y[j]);  // Length of the shorter suffix.
            const idx_t n = m + lcp(str_ + (X[i] + m), str_ + (Y[j] + m), max_n - m);   // LCP(X_i, Y_j)

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
    SA_w = allocate<idx_t>(n_); // Working space for the SA construction.
    LCP_w = allocate<idx_t>(n_);    // Working space for the LCP construction.

    const auto sample_count = p_ * pivot_per_part;
    P = allocate<idx_t>(sample_count);

    const auto idx_init = [SA_ = SA_, SA_w = SA_w](const std::size_t i){ SA_[i] = SA_w[i] = i; };
    parlay::parallel_for(0, n_, idx_init);
}


void Suffix_Array::sort_subarrays()
{
    const auto subarr_size = n_ / p_;   // Size of each subarray to be sorted independently.
    const auto sort_subarr =
        [&](const std::size_t i)
        {
            merge_sort( SA_w + i * subarr_size, SA_ + i * subarr_size,
                        subarr_size + (i < p_ - 1 ? 0 : n_ % p_),
                        LCP_ + i * subarr_size, LCP_w + i * subarr_size);
        };

    parlay::parallel_for(0, p_, sort_subarr, 1);
}


void Suffix_Array::sample_pivots(const idx_t* const X, const idx_t n, const idx_t m, idx_t* const P)
{
    const auto gap = n / (m + 1);   // Distance-gap between pivots.
    for(idx_t i = 0; i < m; ++i)
        P[i] = X[(i + 1) * gap - 1];
}


void Suffix_Array::select_pivots()
{
    const auto sample_count = p_ * pivot_per_part;  // Total number of samples to select pivots from.
    idx_t* const P_w = allocate<idx_t>(sample_count);   // Working space to sample pivots.
    const auto subarr_size = n_ / p_;   // Size of each sorted subarray.

    for(idx_t i = 0; i < p_; ++i)
        sample_pivots(  SA_ + i * subarr_size, subarr_size + (i < p_ - 1 ? 0 : n_ % p_),
                        pivot_per_part, P + i * pivot_per_part);

    auto const temp_1 = allocate<idx_t>(sample_count), temp_2 = allocate<idx_t>(sample_count);

    std::memcpy(P_w, P, sample_count * sizeof(idx_t));
    merge_sort(P, P_w, sample_count, temp_1, temp_2);

    sample_pivots(P_w, sample_count, p_ - 1, P);

    std::free(P_w), std::free(temp_1), std::free(temp_2);
}


std::size_t Suffix_Array::upper_bound(const idx_t* const X, const idx_t n, const char* const q, const std::size_t q_len) const
{
    // Invariant: SA[l] < s < SA[r].

    int64_t l = -1, r = n;  // (Exclusive-) Range of the iterations in the binary search.
    idx_t c;    // Midpoint in each iteration.
    idx_t soln = n; // Solution of the search.
    idx_t lcp_l = 0, lcp_r = 0; // LCP(s, SA[l]) and LCP(s, SA[r]).

    while(r - l > 1)    // Candidate matches exist.
    {
        c = (l + r) / 2;
        const char* const suf = str_ + X[c];    // The suffix at the middle.
        const auto suf_len = n_ - X[c]; // Length of the suffix.

        idx_t lcp_c = std::min(lcp_l, lcp_r);   // LCP(X[c], q).
        const auto max_lcp = std::min(suf_len, q_len);    // Maximum possible LCP, i.e. length of the shorter string.
        lcp_c += lcp(suf + lcp_c, q + lcp_c, max_lcp - lcp_c);  // Skip an informed number of character comparisons.

        if(lcp_c == max_lcp)    // One is a prefix of the other.
        {
            if(lcp_c == q_len)  // q is a prefix of the suffix.
            {
                if(q_len == suf_len)  // The query is the suffix itself, i.e. q = X[c]
                    return c + 1;
                else    // q < X[c]
                    r = c, lcp_r = lcp_c, soln = c;
            }
            else    // The suffix is a prefix of the query, so X[c] < q; technically impossible if the text terminates with $.
                l = c, lcp_l = lcp_c;
        }
        else    // Neither is a prefix of the other.
            if(suf[lcp_c + 1] < q[lcp_c + 1])   // X[c] < q
                l = c, lcp_l = lcp_c;
            else    // q < X[c]
                r = c, lcp_r = lcp_c, soln = c;
    }


    return soln;
}


void Suffix_Array::clean_up()
{
    std::free(SA_w);
    std::free(LCP_w);

    std::free(P);
}


void Suffix_Array::construct()
{
    const auto t_start = now();

    initialize();

    // merge_sort(SA_w, SA_, n_, LCP_, LCP_w);  // Monolithic construction.
    sort_subarrays();

    select_pivots();

    clean_up();

    const auto t_end = now();
    std::cerr << "Constructed the suffix array. Time taken: " << duration(t_end - t_start) << " seconds.\n";
}

}
