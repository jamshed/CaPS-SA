#include <iostream>
#include <fstream>
#include <cstring>
#include <algorithm>

#include "openssl/rand.h"

//#include "parlay/parallel.h"

#define BITMASK(X) ((1ull << X) - 1)
#define DIV_ROUND_UP(X, Y) ((X - 1) / Y + 1)

int main(int argc, char **argv) {
	if (argc < 3) {
		std::cout << "Usage: " << argv[0] << " [output file] [gen length]" << std::endl;
		exit(0);
	}

	std::ofstream writer(argv[1]);
	if (!writer.is_open()) {
		std::cout << "Unable to open file " << argv[1] << std::endl;
		exit(0);
	}

	int gen_length = std::atoi(argv[2]);
	int num_bytes = DIV_ROUND_UP(gen_length, 4);

	unsigned char *rand = new unsigned char[num_bytes];
	RAND_bytes(rand, num_bytes);

	const int buffer_size = 256;
	const char nucl_map[4] = {'A', 'C', 'G', 'T'};
	char buffer[buffer_size + 1];

	for (int i = 0; i < gen_length;) {
		int j = i;
		for (int k = std::min(i + buffer_size, gen_length); i < k; i++) {
			int a = 0b11 << (i & BITMASK(2));
			int b = (rand[i / 4] & (0b11 << (i & BITMASK(2)))) >> (i & BITMASK(2));
			char c = nucl_map[rand[i / 4] & (0b11 << (i & BITMASK(2)))];
			buffer[i - j] = nucl_map[(rand[i / 4] & (0b11 << (i & BITMASK(2)))) >> (i & BITMASK(2))];
		}
		buffer[i - j] = '\0';
		writer << buffer;
	}

	writer.close();
}
