
#include "Bit_Packed_Text.hpp"
#include "parlay/parallel.h"


namespace CaPS_SA
{

Bit_Packed_Text::Bit_Packed_Text(const char* T, std::size_t n):
      T(T)
    , n(n)
    , pack_sz((n + 3) / 4)
    , B(T != nullptr ? allocate<uint8_t>(pack_sz) : nullptr)
{}


void Bit_Packed_Text::construct()
{
    const auto base_code =
        [](char ch)
        {
            ch &= ~32;  // To upper-case.
            assert(ch == 'A' || ch == 'C' || ch == 'G' || ch == 'T');

            return ((ch >> 2) ^ (ch >> 1)) & 0b11;
        };

    const auto pack =
        [&](const std::size_t i)
        {
            assert(i < pack_sz);
            assert(4 * i + 3 < n);

            B[i]  = (base_code(T[i << 2])) |
                    (base_code(T[(i << 2) | 0b01]) << 2) |
                    (base_code(T[(i << 2) | 0b10]) << 4) |
                    (base_code(T[(i << 2) | 0b11]) << 6);
        };


    parlay::parallel_for(0, n / 4, pack);
    if(n & 3)
    {
        B[n / 4] = 0;
        for(std::size_t i = 0; i < (n & 3); ++i)
            B[n / 4] |= (base_code(T[(n / 4) * 4 + i]) << (i * 2));
    }
}

}
