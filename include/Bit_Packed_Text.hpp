
#ifndef BIT_PACKED_TEXT_HPP
#define BIT_PACKED_TEXT_HPP



#include "utilities.hpp"

#include <cstdint>
#include <cstddef>
#include <cstring>
#include <immintrin.h>
#include <cassert>


// =============================================================================

namespace CaPS_SA
{

// A class to to represent genomic texts in a 2-bit-packed manner.
class Bit_Packed_Text
{
private:

    const char* const T;    // The input text.
    const std::size_t n;    // Length of the input text.

    const std::size_t pack_sz;  // Size of the bit-pack.
    uint8_t* const B;   // The bit-packed representation.


    // Returns the 124-nucleobase block from the `i` 'th character, in 256-bits
    // little-endian. No guarantees are provided on the highest byte.
    __m256i load(std::size_t idx) const;


public:

    // Constructs an object for 2-bit-packing a genomic text `T` of length `n`.
    Bit_Packed_Text(const char* T = nullptr, std::size_t n = 0);

    // Constructs the bit-packed representation.
    void construct();

    // Returns the 28-nucleobase block from the `i`'th character, in 64-bits
    // little-endian. No guarantees are provided on the highest byte. 
    uint64_t load28(std::size_t i) const;

    uint64_t loadSmall(const std::size_t i, uint32_t ctx) const;

    // Returns the LCP length of the suffixes `x` and `y`, where `ctx` is the
    // context length.
    std::size_t LCP(std::size_t x, std::size_t y, std::size_t ctx) const;

    // Prints the bit-pack, with left-to-right corresponding to high-to-low
    // indices and bits.
    void print() const;
};


inline __m256i Bit_Packed_Text::load(const std::size_t i) const
{
    assert(i + 124 <= n);

    const auto base_idx = i / 4;    // Base word's index.
    const auto blk = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(B + base_idx));

    const auto base_trail = (i & 3);    // Number of unwanted bases trailing in the base word.
    if(!base_trail)
        return blk;


    const auto to_clear_r = _mm256_set1_epi64x(base_trail * 2);
    const auto trail_cleared = _mm256_srlv_epi64(blk, to_clear_r);
    const auto r_shifted_words = _mm256_permute4x64_epi64(blk, 0b11'11'10'01);
    const auto to_clear_l = _mm256_set1_epi64x((32 - base_trail) * 2);
    const auto lost_bits = _mm256_sllv_epi64(r_shifted_words, to_clear_l);
    const auto restored = _mm256_or_si256(trail_cleared, lost_bits);

    return restored;
}

inline uint64_t Bit_Packed_Text::loadSmall(const std::size_t i, uint32_t ctx) const {
    const auto w_idx = i / 4;
    const auto base_trail = (i & 3);

    uint64_t w;
    std::memcpy(&w, B + w_idx, 8);

    return (w >> (base_trail * 2)) & ((1lu << ((32 - ctx) * 2)) -1);
}

inline uint64_t Bit_Packed_Text::load28(const std::size_t i) const
{
    assert(i + 28 <= n);

    const auto w_idx = i / 4;
    const auto base_trail = (i & 3);

    uint64_t w;
    std::memcpy(&w, B + w_idx, 8);

    return w >> (base_trail *2);
}


inline std::size_t Bit_Packed_Text::LCP(const std::size_t x, const std::size_t y, const std::size_t ctx) const
{
    assert(ctx <= n);

    constexpr std::size_t blk_sz = 124;
    constexpr __mmask32 clear_MSB_mask = ~(__mmask32(1) << 31); // To clear out the unguaranteed top byte from `load`.

    std::size_t i = x, j = y;
    std::size_t lcp = 0;
    for(std::size_t compared = 0; compared + blk_sz <= ctx; compared += blk_sz)
    {
        const auto X = load(i);
        const auto Y = load(j);

        const auto neq_mask = (~_mm256_movemask_epi8(_mm256_cmpeq_epi8(X, Y))) & clear_MSB_mask;
                            //  _mm256_cmpneq_epi8_mask(X, Y) & clear_MSB_mask; // AVX512
        if(neq_mask)
        {
            const auto bytes_eq = __builtin_ctz(neq_mask);
            alignas(32) uint8_t X_bytes[32];
            _mm256_store_si256(reinterpret_cast<__m256i*>(X_bytes), X);
            alignas(32) uint8_t Y_bytes[32];
            _mm256_store_si256(reinterpret_cast<__m256i*>(Y_bytes), Y);
            assert(X_bytes[bytes_eq] != Y_bytes[bytes_eq]);
            const auto bits_eq  = __builtin_ctz(X_bytes[bytes_eq] ^ Y_bytes[bytes_eq]);

            lcp += (bytes_eq << 2) + (bits_eq >> 1);
            return lcp;
        }

        i += blk_sz, j += blk_sz, lcp += blk_sz;
    }

    while(lcp < ctx && T[x + lcp] == T[y + lcp])
        lcp++;

    return lcp;
}

}



#endif
