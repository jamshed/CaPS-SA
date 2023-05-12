
#include "Suffix_Array.hpp"

#include <string>
#include <cstdlib>
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


void pretty_print(const CaPS_SA::Suffix_Array& suf_arr, std::ofstream& output)
{
    const std::size_t n = suf_arr.n();
    for(std::size_t i = 0; i < n; ++i)
        output << suf_arr.SA()[i] << " \n"[i == n - 1];
    for(std::size_t i = 0; i < n; ++i)
        output << suf_arr.LCP()[i] << " \n"[i == n - 1];
}


int main(int argc, char* argv[])
{
    // TODO: standardize the API.
    constexpr auto arg_count = 3;
    if(argc < arg_count)
    {
        std::cerr << "Usage: CaPS_SA <input_path> <output_path> <(optional)-subproblem-count> <(optional)-bounded-context> <(optional)--pretty-print>\n";
        std::exit(EXIT_FAILURE);
    }


    const std::string ip_path(argv[1]);
    const std::string op_path(argv[2]);
    const std::size_t subproblem_count(argc >= 4 ? std::atoi(argv[3]) : 0);
    const std::size_t max_context(argc >= 5 ? std::atoi(argv[4]) : 0);
    const bool print(argc >= 6 ? (std::string(argv[5]) == "--pretty-print") : false);

    std::string text;
    read_input(ip_path, text);

    CaPS_SA::Suffix_Array suf_arr(text.c_str(), text.length(), subproblem_count, max_context);
    suf_arr.construct();

    // TODO: time the o/p.
    std::ofstream output(op_path);
    print ? pretty_print(suf_arr, output) : suf_arr.dump(output);
    output.close();


    return 0;
}
