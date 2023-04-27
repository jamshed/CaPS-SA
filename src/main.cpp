
#include "Suffix_Array.hpp"

#include <string>
#include <cstdlib>
#include <fstream>
#include <iostream>


void pretty_print(const themis::Suffix_Array& suf_arr, std::ofstream& output)
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
        std::cerr << "Usage: themis <input_path> <output_path> <(optional)-subproblem-count> <(optional)--pretty-print>\n";
        std::exit(EXIT_FAILURE);
    }


    const std::string ip_path(argv[1]);
    const std::string op_path(argv[2]);
    const std::size_t subproblem_count(argc == 4 ? std::atoi(argv[3]) : 0);
    const bool print(argc == 5 ? (std::string(argv[4]) == "--pretty-print") : false);
    std::string str;

    std::getline(std::ifstream(ip_path), str);

    themis::Suffix_Array suf_arr(str.c_str(), str.length(), subproblem_count);
    suf_arr.construct();

    std::ofstream output(op_path);
    print ? pretty_print(suf_arr, output) : suf_arr.dump(output);
    output.close();


    return 0;
}
