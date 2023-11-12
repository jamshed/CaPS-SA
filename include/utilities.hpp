
#ifndef UTILITIES_HPP
#define UTILITIES_HPP



#include <cstdint>
#include <cstddef>
#include <climits>
#include <immintrin.h>
#include <cstdlib>
#include <iostream>


// =============================================================================

namespace CaPS_SA
{
    // Returns pointer to a memory-allocation for `size` elements of type `T_`.
    template <typename T_>
    inline T_* allocate(std::size_t size) { return static_cast<T_*>(std::malloc(size * sizeof(T_))); }


namespace Debug
{
    // Prints the binary representation of `x`, with left-to-right corresponding
    // to high-to-low bits.
    template <typename T_>
    inline void print_bin(T_ x)
    {
        for(int i = sizeof(T_) * CHAR_BIT - 1; i >= 0; --i)
            std::cerr << ((x & (uint64_t(1) << i)) ? 1 : 0);
    }


    // Prints the hex representation of `x`.
    template <typename T_>
    inline void print_hex(T_ x) { std::cerr << std::hex << x; }


    // Prints the hex representation of the 256-bit register `val`.
    inline void print256_hex(__m256i val)
    {
        alignas(32) int64_t W[4];
        _mm256_store_si256(reinterpret_cast<__m256i*>(W), val);

        for(int i = 3; i >= 0; --i)
        {
            print_hex(W[i]);
            std::cerr << " ";
        }

        std::cerr << "\n";
    }

    // Prints the binary representation of the 256-bit register `val`.
    inline void print256_bin(__m256i val)
    {
        alignas(32) int64_t W[4];
        _mm256_store_si256(reinterpret_cast<__m256i*>(W), val);

        for(int i = 3; i >= 0; --i)
        {
            print_bin(W[i]);
            std::cerr << " ";
        }

        std::cerr << "\n";
    }
}
};



#endif