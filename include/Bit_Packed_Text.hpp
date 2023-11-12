
#ifndef BIT_PACKED_TEXT_HPP
#define BIT_PACKED_TEXT_HPP



#include "utilities.hpp"

#include <cstdint>
#include <cstddef>
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

}



#endif
