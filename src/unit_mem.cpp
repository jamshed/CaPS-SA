#include <iostream>

#include <cstring>

#include <openssl/rand.h>

#include "memoization_table.h"

#define BITMASK(X) ((1ull << X) - 1)
#define DIV_ROUND_UP(X, Y) ((X - 1) / Y + 1)
#define ROUND_UP(X, Y) (DIV_ROUND_UP(X, Y) * Y)

using index_t = uint32_t;

index_t lcp(char *str, index_t len, index_t a, index_t b) {
    if (a == b) return len - a;
    index_t a_old = a;
    while (str[a++] == str[b++]);
    return a - a_old - 1;
}

int main() {
    const index_t n = 1000;

    char string[ROUND_UP(n + 1, 32)];
    string[n] = '\0';

    const char nucl_map[4] = {'A', 'C', 'G', 'T'};
    RAND_bytes((unsigned char*)&string[0], DIV_ROUND_UP(n, 4));
    for (int i = n - 1; i >= 0; i--) {
        string[i] = 'A';//nucl_map[(string[i >> 2] & (0b11 << (i & 0b11))) >> (i & 0b11)];
    }

    memoization_table<index_t> mem_table(n);
    
    const int num_tests = 10000;
    srand(time(NULL));
    for (int i = 0; i < num_tests; i++) {
        int a = rand() % n, b = rand() % n;
        int exp_lcp = lcp(string, n, a, b);
        int obt_lcp = mem_table.get_lcp(string, n, a, b);
        if (exp_lcp != obt_lcp) {
            std::cout << "lcp(" << a << ", " << b << ") is wrong" << std::endl;
            mem_table.get_lcp(string, n, a, b);
        }
    }

    std::cout << "Done" << std::endl;
    return 0;
}
