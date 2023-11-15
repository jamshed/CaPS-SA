#include "CLI11.hpp"
#include "Suffix_Array.hpp"
#include "parlay/parallel.h"

#include <algorithm>
#include <string>
#include <cstdlib>
#include <limits>
#include <fstream>
#include <iostream>
#include <filesystem>


void read_input(const std::string& ip_path, std::string& text)
{
    std::error_code ec;
    const auto file_size = std::filesystem::file_size(ip_path, ec);

    if(ec)
    {
        std::cerr << ip_path << " : " << ec.message() << "\n";
        std::exit(EXIT_FAILURE);
    }

    // pad by 8 here so we can always read
    // suffixes into `uint64_t`s without worrying
    // about going out of bounds.
    text.resize(file_size + 8);
    std::ifstream input(ip_path);
    input.read(text.data(), file_size);
    input.close();
}


template <typename T_idx_>
void pretty_print(const CaPS_SA::Suffix_Array<T_idx_>& suf_arr, std::ofstream& output)
{
    const std::size_t n = suf_arr.n();
    for(std::size_t i = 0; i < n; ++i)
        output << suf_arr.SA()[i] << " \n"[i == n - 1];
    for(std::size_t i = 0; i < n; ++i)
        output << suf_arr.LCP()[i] << " \n"[i == n - 1];
}


int main(int argc, char* argv[]) {
    CLI::App app{"Build the suffix array using the divsufsort algorithm"};
    std::string ip_path;
    std::string op_path;

    app.add_option("-i,--input", ip_path, "input file on which to build sa")->required();
    app.add_option("-o,--output", op_path, "ouput file where SA should be written");
    CLI11_PARSE(app, argc, argv);

    std::string text;
    read_input(ip_path, text);
    constexpr char lookup[4] = {'A', 'C', 'T', 'G'};
    // the actual size of the text is 8 less than 
    // text.length() because of the padding.
    std::size_t n = text.length() - 8;
    parlay::blocked_for(0, n, 65536, 
      [&, n](size_t i, size_t start, size_t end) {
        (void)i;
        for (size_t j = start; j < std::min(end, n); ++j) {
          char c = text[j];
          text[j] = lookup[((std::toupper(c) & 0x6) >> 1)];
        };
    });
  std::ofstream output(op_path);
  if (output.is_open()) {
    output.write(reinterpret_cast<const char*>(text.data()), n);
    if (output.good()) {
      std::cerr << "Dumped the converted string to " << op_path << "\n";
    } else {
      std::cerr << "Error writing converted string to " << op_path << "\n";
    }
  } else {
    std::cerr << "Error opening output file " << op_path << "\n";
  }
  return 0;
} 
