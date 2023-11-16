
#include "Suffix_Array.hpp"
#include "parlay/parallel.h"
#include "CLI11/CLI11.hpp"

#include <algorithm>
#include <string>
#include <cstdlib>
#include <limits>
#include <fstream>
#include <iostream>
#include <filesystem>


// Reads the input string at file-path `ip_path` into `text`, pads it with 8
// empty characters, and returns the size of the read string.
std::size_t read_and_pad_input(const std::string& ip_path, char*& text)
{
    const auto t_i = CaPS_SA::now();

    std::error_code ec;
    const auto file_size = std::filesystem::file_size(ip_path, ec);

    if(ec)
    {
        std::cerr << ip_path << " : " << ec.message() << "\n";
        std::exit(EXIT_FAILURE);
    }

    text = CaPS_SA::allocate<char>(file_size + 8);
    std::ifstream input(ip_path);
    input.read(text, file_size);
    input.close();

    const auto t_e = CaPS_SA::now();
    std::cerr << "Time to read input: " << CaPS_SA::duration(t_e - t_i) << "s.\n";

    return file_size;
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
    std::size_t subproblem_count(0);//(argc >= 4 ? std::atoi(argv[3]) : 0);
    std::size_t max_context(0);//(argc >= 5 ? std::atoi(argv[4]) : 0);

    app.add_option("-i,--input", ip_path, "input file on which to build sa")->required();
    app.add_option("-o,--output", op_path, "ouput file where SA should be written");
    // app.add_option("-c,--context", max_context, "bounded context length")->default_str("0");
    // app.add_option("-s,--subprob-count", subproblem_count, "subproblem count")->default_str("0");
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



    char* text = nullptr;
    const std::size_t n = read_and_pad_input(ip_path, text);

    std::ofstream output(op_path);

    std::cerr << "Text length: " << n << ".\n";
    if(n <= std::numeric_limits<uint32_t>::max())
    {
        CaPS_SA::Suffix_Array<uint32_t> suf_arr(text, n, subproblem_count, max_context);
        suf_arr.construct();
        suf_arr.dump(output);

        assert(suf_arr.is_sorted(suf_arr.SA(), suf_arr.n()));
        // std::cerr << "Sortedness: " << suf_arr.is_sorted(suf_arr.SA(), suf_arr.n()) << "\n";
    }
    else
    {
        CaPS_SA::Suffix_Array<uint64_t> suf_arr(text, n, subproblem_count, max_context);
        suf_arr.construct();
        suf_arr.dump(output);

        assert(suf_arr.is_sorted(suf_arr.SA(), suf_arr.n()));
        // std::cerr << "Sortedness: " << suf_arr.is_sorted(suf_arr.SA(), suf_arr.n()) << "\n";
    }

    output.close();
    std::free(text);


    return 0;
}
