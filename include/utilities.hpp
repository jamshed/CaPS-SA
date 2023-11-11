
#ifndef UTILITIES_HPP
#define UTILITIES_HPP



#include <cstddef>
#include <cstdlib>


namespace CaPS_SA
{
    // Returns pointer to a memory-allocation for `size` elements of type `T_`.
    template <typename T_>
    T_* allocate(std::size_t size) { return static_cast<T_*>(std::malloc(size * sizeof(T_))); }
};



#endif
