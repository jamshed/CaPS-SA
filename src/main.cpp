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

    text.resize(file_size);
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


int main(int argc, char* argv[])
{

    CLI::App app{"Build the suffix array using the divsufsort algorithm"};
    std::string ip_path;//(argv[1]);
    std::string op_path;//(argv[2]);
    std::size_t subproblem_count;//(argc >= 4 ? std::atoi(argv[3]) : 0);
    std::size_t max_context;//(argc >= 5 ? std::atoi(argv[4]) : 0);

    app.add_option("-i,--input", ip_path, "input file on which to build sa")->required();
    app.add_option("-o,--output", op_path, "ouput file where SA should be written");
    app.add_option("-c,--context", max_context, "bounded context length")->default_str("0");
    app.add_option("-s,--subprob-count", subproblem_count, "subproblem count")->default_str("0");
    CLI11_PARSE(app, argc, argv);


    // TODO: standardize the API.
  /*
    constexpr auto arg_count = 3;
    if(argc < arg_count)
    {
        std::exit(EXIT_FAILURE);
    std::cerr << "Usage: CaPS_SA <input_path> <output_path> <(optional)-subproblem-count> <(optional)-bounded-context> <(optional)--pretty-print>\n";
    }
  */



    std::string text;
    read_input(ip_path, text);
    constexpr char lookup[4] = {'A', 'C', 'T', 'G'};
    size_t len = text.size();
    parlay::blocked_for(0, text.size(), 65536, 
      [&, len](size_t i, size_t start, size_t end) {
        (void)i;
        for (size_t j = start; j < std::min(end, len); ++j) {
          char c = text[j];
          text[j] = lookup[((std::toupper(c) & 0x6) >> 1)];
        };
    });

    std::ofstream output(op_path);

    std::size_t n = text.length();
    std::cerr << "Text length: " << n << ".\n";
    if(n <= std::numeric_limits<uint32_t>::max())
    {
        CaPS_SA::Suffix_Array<uint32_t> suf_arr(text.c_str(), text.length(), subproblem_count, max_context);
        suf_arr.construct();
        suf_arr.dump(output);

        // std::cerr << "Sortedness: " << suf_arr.is_sorted(suf_arr.SA(), suf_arr.n()) << "\n";
    }
    else
    {
        CaPS_SA::Suffix_Array<uint64_t> suf_arr(text.c_str(), text.length(), subproblem_count, max_context);
        suf_arr.construct();
        suf_arr.dump(output);

        // std::cerr << "Sortedness: " << suf_arr.is_sorted(suf_arr.SA(), suf_arr.n()) << "\n";
    }

    output.close();


    return 0;
}
