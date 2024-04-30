#include <iostream>
#include <fstream>
#include <cstring>
#include <algorithm>

#include "openssl/rand.h"

//#include "parlay/parallel.h"

#define BITMASK(X) ((1ull << X) - 1)
#define DIV_ROUND_UP(X, Y) ((X - 1) / Y + 1)

int main(int argc, char **argv) {
    if (argc < 5) {
        std::cout << "Usage: " << argv[0] << " [output file] [gen length] [num words] [word length]" << std::endl;
        exit(0);
    }

    std::ofstream writer(argv[1]);
    if (!writer.is_open()) {
        std::cout << "Unable to open file " << argv[1] << std::endl;
        exit(0);
    }

    int gen_length = std::atoi(argv[2]);
    int num_bytes = DIV_ROUND_UP(gen_length, 4);

    int num_words = std::atoi(argv[3]);
    int word_len = std::atoi(argv[4]);

    unsigned char *rand = new unsigned char[num_words * (word_len + 1)];
    RAND_bytes(rand, num_words * word_len);

    const char nucl_map[4] = {'A', 'C', 'G', 'T'};
    
    char *words = new char[num_words * (word_len + 1)];
    for (int i = 0; i < num_words * (word_len + 1); i++) {
        words[i] = nucl_map[(rand[i / 4] & (0b11 << (i & BITMASK(2)))) >> (i & BITMASK(2))];
    }
    for (int i = 0; i < num_words; i++) {
        words[(i + 1) * (word_len + 1) - 1] = '\0';
    }

    int num_chunks = (gen_length - 1) / word_len + 1;
    unsigned int *word_choices = new unsigned int[num_chunks];
    RAND_bytes((unsigned char*)word_choices, num_chunks * sizeof(int));

    for (int i = 0; i < num_chunks; i++) {
        int word_ptr = (word_choices[i] % num_words) * (word_len + 1);
        int remaining_chars = gen_length - (i * word_len);
        if (remaining_chars < word_len) {
            words[word_ptr + remaining_chars] = '\0';
        }
        writer << &words[word_ptr];
    }

    writer.close();
}
