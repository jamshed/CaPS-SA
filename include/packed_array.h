#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <bit>

#include <iostream>
#include <fstream>

#include <immintrin.h>

#define BITS_PER_CHAR 2

#define BITMASK(X) ((1ull << X) - 1)
#define DIV_ROUND_UP(X, Y) ((X - 1) / Y) + 1
#define BYTE_ROUND_UP(BITS) ((BITS + 7) / 8)
#define LOG2_ROUND_UP(X) (X ? (63 - __builtin_clzll(X - 1)) + 1 : 0)

static inline uint64_t extract(unsigned char *bitvector, int i, int n) {
	return (*((uint64_t*)&bitvector[i / 8]) >> (i % 8)) & BITMASK(n);
}

static inline int popcnt(uint64_t x) {
	asm("popcnt %[x], %[x]"
			: [x] "+r" (x)
			:
			: "cc");
	return x;
}

static inline int64_t bsr(uint64_t x) {
	if (!x) return -1;
	asm("bsr %[x], %[x]"
			: [x] "+r" (x)
			:
			: "cc");
	return x;
}

class packed_array {
public:
	int arr_size;
	int int_size;
	int num_bytes;
	unsigned char *data;

	packed_array() {
		arr_size = 0;
		int_size = 0;
		data = NULL;
	}
	packed_array(size_t _arr_size, int _int_size) : arr_size(_arr_size), int_size(_int_size) {
		num_bytes = BYTE_ROUND_UP(arr_size * int_size);
		data = new unsigned char[num_bytes];
		std::memset(data, (char)0, num_bytes);
	}

	std::string see_bits(size_t k) const {
		std::string s = "";
		for (int i = 0; i < 8; i++) {
			for (int j = 0; j < 8; j++) {
				s += (data[k + i] & (1 << j) ? "1" : "0");
			}
		}
		return s;
	}

	inline uint64_t get(size_t i) const {
		return extract(data, i * int_size, int_size);
	}

	inline void set(const size_t i, const uint64_t x) {
		uint64_t prior_bits = i * int_size;
		uint64_t bit_offset = prior_bits % 8;
		uint64_t byte_offset = prior_bits / 8;
		uint64_t *value = (uint64_t*)&data[byte_offset];
		*value = *value & ~(BITMASK(int_size) << bit_offset) | ((x & BITMASK(int_size)) << bit_offset);
	}

    size_t lcp(index_t a, index_t b) {
        size_t lcp = 0;

        __m256i suff_a = _mm256_loadu_si256((__m256i*)&data[a]);

        __m256i eq = _mm256_cmpeq_epi8(suff_a, suff_b);
        int mask = _mm256_movemask_epi8(eq);

    }
};

