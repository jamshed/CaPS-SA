
#ifndef BIT_PACKED_TEXT_HPP
#define BIT_PACKED_TEXT_HPP



#include "utilities.hpp"

#include <cstdint>
#include <cstddef>


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


public:

    // Constructs an object for 2-bit-packing a genomic text `T` of length `n`.
    Bit_Packed_Text(const char* T = nullptr, std::size_t n = 0);

    // Constructs the bit-packed representation.
    void construct();

    // Prints the bit-pack, with left-to-right corresponding to high-to-low
    // indices and bits.
    void print();
};

}



#endif
