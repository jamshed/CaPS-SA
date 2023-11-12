
#include "utilities.hpp"


namespace CaPS_SA
{

namespace Debug
{

void print256_hex(__m256i val)
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


void print256_bin(__m256i val)
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

}
