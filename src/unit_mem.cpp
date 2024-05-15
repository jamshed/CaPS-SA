#include <iostream>

#include <cstring>

#include <openssl/rand.h>

#include "memoization_table.h"

#define BITMASK(X) ((1ull << X) - 1)
#define DIV_ROUND_UP(X, Y) ((X - 1) / Y + 1)

using index_t = uint32_t;

index_t lcp(char *str, index_t a, index_t b) {
    index_t a_old = a;
    while (str[a++] == str[b++]);
    return a - a_old - 1;
}

int main() {
    const index_t n = 20;

    char string[n + 1];
    string[n] = '\0';

    const char nucl_map[4] = {'A', 'C', 'G', 'T'};
    RAND_bytes((unsigned char*)&string[0], DIV_ROUND_UP(n, 4));
    for (int i = n - 1; i >= 0; i--) {
        string[i] = nucl_map[string[i >> 2] & ((0b11 << (i & 0b11)) >> (i & 0b11))];
    }

    memoization_table<index_t> mem_table(n);
    
    

    std::cout << "Done" << std::endl;
    return 0;
}
