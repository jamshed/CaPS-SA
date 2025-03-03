
#include "Suffix_Array.hpp"
#include "parlay/parallel.h"

#include <cstring>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <cassert>

namespace CaPS_SA
{

template <typename T_idx_>
Suffix_Array<T_idx_>::Suffix_Array(const char* const T, const idx_t n, const idx_t subproblem_count, const idx_t max_context):
    T_(T),
    n_(n),
    SA_(allocate<idx_t>(n_)),
    LCP_(allocate<idx_t>(n_)),
    SA_w(nullptr),
    LCP_w(nullptr),
    p_(std::min(subproblem_count > 0 ? subproblem_count : default_subproblem_count, n / 16)),   // TODO: fix subproblem-count for small `n`.
    max_context(max_context ? max_context : n_),
    pivot_(nullptr),
    pivot_per_part_(std::min(static_cast<idx_t>(std::ceil(32.0 * std::log(n_))), n_ / p_ - 1)), // (c \ln n) or (|subarray| - 1)
    part_size_scan_(nullptr),
    part_ruler_(nullptr)
{
    assert(n_ >= 16);   // TODO: fix subproblem-count for small `n`.

    if(p_ > n_)
    {
        std::cerr << "Incompatible subproblem-count. Aborting.\n";
        std::exit(EXIT_FAILURE);
    }
}


template <typename T_idx_>
Suffix_Array<T_idx_>::~Suffix_Array()
{
    deallocate(SA_), deallocate(LCP_);
}


template <typename T_idx_>
void Suffix_Array<T_idx_>::merge(const idx_t* X, idx_t len_x, const idx_t* Y, idx_t len_y, const idx_t* LCP_x, const idx_t* LCP_y, idx_t* Z, idx_t* LCP_z) const
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
            const idx_t context = std::min(max_context, max_n); // Prefix-context length for the suffixes.
            const idx_t n = m + lcp_opt_avx_unrolled(T_ + (X[i] + m), T_ + (Y[j] + m), context - m); // LCP(X_i, Y_j)

            // Whether the shorter suffix is a prefix of the longer one.
            Z[k] = (n == max_n ?    std::max(X[i], Y[j]) :
                                    (T_[X[i] + n] < T_[Y[j] + n] ? X[i] : Y[j]));
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


template <typename T_idx_>
void Suffix_Array<T_idx_>::merge_sort(idx_t* const X, idx_t* const Y, const idx_t n, idx_t* const LCP, idx_t* const W) const
{
    assert(std::memcmp(X, Y, n * sizeof(idx_t)) == 0);

    if(n == 1)
        LCP[0] = 0;
    else
    {
        const idx_t m = n / 2;
        const auto f = [&](){ merge_sort(Y, X, m, W, LCP); };
        const auto g = [&](){ merge_sort(Y + m, X + m, n - m, W + m, LCP + m); };

        m < nested_par_grain_size ?
            (f(), g()) : parlay::par_do(f, g);
        merge(X, m, X + m, n - m, W, W + m, Y, LCP);
    }
}


template <typename T_idx_>
void Suffix_Array<T_idx_>::initialize()
{
    const auto t_s = now();

    SA_w = allocate<idx_t>(n_); // Working space for the SA construction.
    LCP_w = allocate<idx_t>(n_);    // Working space for the LCP construction.

    const auto sample_count = p_ * pivot_per_part_;
    pivot_ = allocate<idx_t>(sample_count);

    const auto t_e = now();
    std::cerr << "Initialized required data structures. Time taken: " << duration(t_e - t_s) << " seconds.\n";
}


template <typename T_idx_>
void Suffix_Array<T_idx_>::permute()
{
    const auto t_s = now();

    const auto populate = [&](const idx_t i){ SA_[i] = SA_w[i] = i; };
    parlay::parallel_for(0, n_, populate);

    const auto t_e = now();
    std::cerr << "Populated the suffix array with a permutation. Time taken: " << duration(t_e - t_s) << " seconds.\n";
}


template <typename T_idx_>
void Suffix_Array<T_idx_>::sort_subarrays()
{
    const auto t_s = now();

    const auto subarr_size = n_ / p_;   // Size of each subarray to be sorted independently.
    std::atomic_uint64_t solved_ = 0;   // Progress tracker—number of subproblems solved in some step.
    const auto sort_subarr =
        [&](const idx_t i)
        {
            merge_sort( SA_w + i * subarr_size, SA_ + i * subarr_size,
                        subarr_size + (i < p_ - 1 ? 0 : n_ % p_),
                        LCP_ + i * subarr_size, LCP_w + i * subarr_size);

            if(++solved_ % 8 == 0)
                std::cerr << "\rSorted " << solved_ << " subarrays.";
        };

    parlay::parallel_for(0, p_, sort_subarr, 1);
    std::cerr << "\n";

    const auto t_e = now();
    std::cerr << "Sorted the subarrays independently. Time taken: " << duration(t_e - t_s) << " seconds.\n";
}


template <typename T_idx_>
void Suffix_Array<T_idx_>::sample_pivots(const idx_t* const X, const idx_t n, const idx_t m, idx_t* const P)
{
    assert(m <= n);
    const auto gap = n / m; // Distance-gap between pivots.
    for(idx_t i = 0; i < m; ++i)
        P[i] = X[(i + 1) * gap - 1];
}


template <typename T_idx_>
void Suffix_Array<T_idx_>::select_pivots()
{
    const auto t_s = now();

    const auto sample_count = p_ * pivot_per_part_; // Total number of samples to select pivots from.
    idx_t* const pivot_w = allocate<idx_t>(sample_count);   // Working space to sample pivots.
    const auto subarr_size = n_ / p_;   // Size of each sorted subarray.

    assert(pivot_per_part_ < subarr_size);
    for(idx_t i = 0; i < p_; ++i)
        sample_pivots(  SA_ + i * subarr_size, subarr_size + (i < p_ - 1 ? 0 : n_ % p_),
                        pivot_per_part_, pivot_ + i * pivot_per_part_);

    auto const temp_1 = allocate<idx_t>(sample_count), temp_2 = allocate<idx_t>(sample_count);

    std::memcpy(pivot_w, pivot_, sample_count * sizeof(idx_t));
    merge_sort(pivot_, pivot_w, sample_count, temp_1, temp_2);

    sample_pivots(pivot_w, sample_count, p_ - 1, pivot_);

    deallocate(pivot_w), deallocate(temp_1), deallocate(temp_2);

    const auto t_e = now();
    std::cerr << "Selected the global pivots. Time taken: " << duration(t_e - t_s) << " seconds.\n";
}


template <typename T_idx_>
void Suffix_Array<T_idx_>::locate_pivots(idx_t* const P) const
{
    const auto t_s = now();

    const auto subarr_size = n_ / p_;   // Size of each independent sorted subarray.

    const auto locate =
        [&](const idx_t i)
        {
            const auto X_i = SA_ + i * subarr_size; // The i'th subarray.
            const auto P_i = P + i * (p_ + 1);  // Pivot locations in `X_i` are to be placed in `P_i`.

            const auto tot_subarr_size = subarr_size + (i < p_ - 1 ? 0 : n_ % p_);
            P_i[0] = 0, P_i[p_] = tot_subarr_size; // The two flanking pivot indices.

            for(idx_t j = 0; j < p_ - 1; ++j) // TODO: try parallelizing this loop too; observe performance diff.
                P_i[j + 1] = upper_bound(X_i, P_i[p_], T_ + pivot_[j], n_ - pivot_[j]);
        };

    parlay::parallel_for(0, p_, locate, 1);

    const auto t_e = now();
    std::cerr << "Located the pivots in each sorted subarray. Time taken: " << duration(t_e - t_s) << " seconds.\n";
}


template <typename T_idx_>
T_idx_ Suffix_Array<T_idx_>::upper_bound(const idx_t* const X, const idx_t n, const char* const P, const idx_t P_len) const
{
    // Invariant: SA[l] < P < SA[r].

    int64_t l = -1, r = n;  // (Exclusive-) Range of the iterations in the binary search.
    idx_t c;    // Midpoint in each iteration.
    idx_t soln = n; // Solution of the search.
    idx_t lcp_l = 0, lcp_r = 0; // LCP(P, SA[l]) and LCP(P, SA[r]).
	constexpr idx_t cutoff = 65536; // TODO: better tune and document.

    while(r - l > 1)    // Candidate matches exist.
    {
        c = (l + r) / 2;
        auto const suf = T_ + X[c]; // The suffix at the middle.
        const auto suf_len = n_ - X[c]; // Length of the suffix.

        idx_t lcp_c = std::min(lcp_l, lcp_r);   // LCP(X[c], P).
        lcp_c = std::min(lcp_c, cutoff);
        auto max_lcp = std::min(std::min(suf_len, P_len), max_context); // Maximum possible LCP, i.e. length of the shorter string.
        max_lcp = std::min(max_lcp, cutoff);
        lcp_c += lcp_opt_avx_unrolled(suf + lcp_c, P + lcp_c, max_lcp - lcp_c);  // Skip an informed number of character comparisons.

        if(lcp_c == max_lcp)    // One is a prefix of the other, or they align at least up-to the context- or the cutoff-length.
        {
            if(lcp_c == P_len)  // P is a prefix of the suffix.
            {
                if(P_len == suf_len)    // The query is the suffix itself, i.e. P = X[c]
                    return c + 1;
                else    // P < X[c]
                    r = c, lcp_r = lcp_c, soln = c;
            }
            else    // The suffix is a prefix of the query, so X[c] < P; technically impossible if the text terminates with $.
                    // Or, their relevant prefixes align; moving to the right, as searching for the upper-bound.
                l = c, lcp_l = lcp_c;
        }
        else    // They mismatch within their relevant prefixes.
            if(suf[lcp_c] < P[lcp_c])   // X[c] < P
                l = c, lcp_l = lcp_c;
            else    // P < X[c]
                r = c, lcp_r = lcp_c, soln = c;
    }


    return soln;
}


template <typename T_idx_>
void Suffix_Array<T_idx_>::partition_sub_subarrays(const idx_t* const P)
{
    const auto t_s = now();

    part_size_scan_ = allocate<idx_t>(p_ + 1);

    const auto collect_size =   // Collects the size of the `j`'th partition.
        [&](const idx_t j)
        {
            part_size_scan_[j] = 0;
            for(idx_t i = 0; i < p_; ++i)   // For subarray `i`.
            {
                const auto P_i = P + i * (p_ + 1);  // Pivot collection of subarray `i`.
                part_size_scan_[j] += (P_i[j + 1] - P_i[j]);    // Collect its `j`'th sub-subarray's size.
            }
        };

    parlay::parallel_for(0, p_, collect_size, 1);   // Collect the individual size of each partition.


    // Compute inclusive-scan (prefix sum) of the partition sizes.
    idx_t curr_sum = 0;
    for(idx_t j = 0; j < p_; ++j) // For partition `j`.
    {
        const auto part_size = part_size_scan_[j];

        part_size_scan_[j] = curr_sum;
        curr_sum += part_size;
    }

    part_size_scan_[p_] = curr_sum;
    assert(part_size_scan_[p_] == n_);


    // Collate the sorted sub-subarrays to appropriate partitions.
    part_ruler_ = allocate<idx_t>(p_ * (p_ + 1));
    const idx_t subarr_size = n_ / p_;
    const auto collate =    // Collates the `j`'th sub-subarray from each sorted subarray to partition `j`.
        [&](const idx_t j)
        {
            auto const Y_j = SA_w + part_size_scan_[j]; // Memory-base for partition `j`.
            auto const LCP_Y_j = LCP_w + part_size_scan_[j];    // Memory-base for LCPs of partition `j`.
            auto const sub_subarr_idx = part_ruler_ + j * (p_ + 1); // Index of the sorted sub-subarrays in `Y_j`.
            idx_t curr_idx = 0; // Current index into `Y_j`.

            for(idx_t i = 0; i < p_; ++i)   // Subarray `i`.
            {
                const auto X_i = SA_ + i * subarr_size; // `i`'th sorted subarray.
                const auto LCP_X_i = LCP_ + i * subarr_size;    // LCP array of `X_i`.
                const auto P_i = P + i * (p_ + 1);  // Pivot collection of subarray `i`.

                const auto sub_subarr_size = P_i[j + 1] - P_i[j];   // Size of the `j`'th sub-subarray of subarray `i`.
                sub_subarr_idx[i] = curr_idx;
                std::memcpy(Y_j + sub_subarr_idx[i], X_i + P_i[j], sub_subarr_size * sizeof(idx_t));
                std::memcpy(LCP_Y_j + sub_subarr_idx[i], LCP_X_i + P_i[j], sub_subarr_size * sizeof(idx_t));
                LCP_Y_j[sub_subarr_idx[i]] = 0;
                curr_idx += sub_subarr_size;
            }

            sub_subarr_idx[p_] = curr_idx;
            assert(curr_idx == part_size_scan_[j + 1] - part_size_scan_[j]);
        };

    parlay::parallel_for(0, p_, collate, 1);

    const auto t_e = now();
    std::cerr << "Collated the sorted sub-subarrays into partitions. Time taken: " << duration(t_e - t_s) << " seconds.\n";
}


template <typename T_idx_>
void Suffix_Array<T_idx_>::merge_sub_subarrays()
{
    const auto t_s = now();

    const auto mem_init =
        [&](const idx_t j)
        {
            const auto part_size = part_size_scan_[j + 1] - part_size_scan_[j];
            std::memcpy(SA_ + part_size_scan_[j], SA_w + part_size_scan_[j], part_size * sizeof(idx_t));
            std::memcpy(LCP_ + part_size_scan_[j], LCP_w + part_size_scan_[j], part_size * sizeof(idx_t));
        };

    parlay::parallel_for(0, p_, mem_init, 1);   // Fulfill `sort_partition`'s precondition.


    std::atomic_uint64_t solved_ = 0;   // Progress tracker—number of subproblems solved in some step.
    const auto sort_part =
        [&](const idx_t j)
        {
            const auto part_idx = part_size_scan_[j];   // Index of the partition in the partitions' flat collection.
            auto const X_j = SA_w + part_idx;   // Memory-base for partition `j`.
            auto const Y_j = SA_ + part_idx;    // Location to sort partition `j`.
            auto const LCP_X_j = LCP_w + part_idx;  // Memory-base for the LCP-arrays of partition `j`.
            auto const LCP_Y_j = LCP_ + part_idx;   // LCP array of `Y_j`.
            auto const sub_subarr_idx = part_ruler_ + j * (p_ + 1); // Indices of the sorted subarrays in `X_i`.

            sort_partition(X_j, Y_j, p_, sub_subarr_idx, LCP_X_j, LCP_Y_j);

            if(++solved_ % 8 == 0)
                std::cerr << "\rMerged " << solved_ << " partitions.";
        };

    parlay::parallel_for(0, p_, sort_part, 1);  // Merge the sorted subarrays in each partitions.
    std::cerr << "\n";

    const auto t_e = now();
    std::cerr << "Merged the sorted subarrays in each partition. Time taken: " << duration(t_e - t_s) << " seconds.\n";
}


template <typename T_idx_>
void Suffix_Array<T_idx_>::sort_partition(idx_t* const X, idx_t* const Y, const idx_t n, const idx_t* const S, idx_t* const LCP_x, idx_t* const LCP_y)
{
    if(n == 1)
        return;

    const auto m = n / 2;
    const auto flat_count_l = S[m] - S[0];
    const auto flat_count_r = S[n] - S[m];

    const auto f = [&](){ sort_partition(Y, X, m, S, LCP_y, LCP_x); };
    const auto g = [&](){ sort_partition(Y + flat_count_l, X + flat_count_l, n - m, S + m, LCP_y + flat_count_l, LCP_x + flat_count_l); };

    (flat_count_l < nested_par_grain_size || flat_count_r < nested_par_grain_size) ?
        (f(), g()) : parlay::par_do(f, g);
    merge(X, flat_count_l, X + flat_count_l, flat_count_r, LCP_x, LCP_x + flat_count_l, Y, LCP_y);
}


template <typename T_idx_>
void Suffix_Array<T_idx_>::compute_partition_boundary_lcp()
{
    const auto t_s = now();

    const auto compute_boundary_lcp =
        [&](const idx_t j)
        {
          const auto part_idx = part_size_scan_[j];
          LCP_[part_idx] = lcp_opt_avx_unrolled(T_ + SA_[part_idx - 1], T_ + SA_[part_idx], n_ - std::max(SA_[part_idx - 1], SA_[part_idx]));
        };

    parlay::parallel_for(1, p_, compute_boundary_lcp, 1);

    const auto t_e = now();
    std::cerr << "Computed the LCPs at the partition boundaries. Time taken: " << duration(t_e - t_s) << " seconds.\n";
}


template <typename T_idx_>
void Suffix_Array<T_idx_>::clean_up()
{
    const auto t_s = now();

    deallocate(SA_w), deallocate(LCP_w);

    deallocate(pivot_);
    deallocate(part_size_scan_);
    deallocate(part_ruler_);

    const auto t_e = now();
    std::cerr << "Released the temporary data structures. Time taken: " << duration(t_e - t_s) << " seconds.\n";
}


template <typename T_idx_>
void Suffix_Array<T_idx_>::construct()
{
    const auto t_start = now();

    initialize();

    permute();

    // merge_sort(SA_w, SA_, n_, LCP_, LCP_w);  // Monolithic construction.
    sort_subarrays();

    select_pivots();

    idx_t* const P = allocate<idx_t>(p_ * (p_ + 1));  // Collection of pivot locations in the subarrays.
    locate_pivots(P);
    partition_sub_subarrays(P);
    deallocate(P);

    merge_sub_subarrays();

    compute_partition_boundary_lcp();

    clean_up();

    const auto t_end = now();
    std::cerr << "Constructed the suffix array. Time taken: " << duration(t_end - t_start) << " seconds.\n";
}


template <typename T_idx_>
void Suffix_Array<T_idx_>::dump(std::ofstream& output)
{
    const auto t_start = now();

    const std::size_t n = n_;
    output.write(reinterpret_cast<const char*>(&n), sizeof(std::size_t));
    output.write(reinterpret_cast<const char*>(SA_), n_ * sizeof(idx_t));
    output.write(reinterpret_cast<const char*>(LCP_), n_ * sizeof(idx_t));

    const auto t_end = now();
    std::cerr << "Dumped the suffix array. Time taken: " << duration(t_end - t_start) << " seconds.\n";
}


template <typename T_idx_>
bool Suffix_Array<T_idx_>::is_sorted(const idx_t* const X, const idx_t n) const
{
    for(idx_t i = 1; i < n; ++i)
    {
        const auto x = T_ + X[i - 1], y = T_ + X[i];
        const auto l = std::min(n_ - X[i - 1], n_ - X[i]);

        for(idx_t i = 0; i < l; ++i)
            if(x[i] < y[i])
                break;
            else if(x[i] > y[i])
                return false;
    }

    return true;
}

}



// Template instantiations for the required instances.
template class CaPS_SA::Suffix_Array<uint32_t>;
template class CaPS_SA::Suffix_Array<uint64_t>;
