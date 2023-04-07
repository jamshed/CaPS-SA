
#include "Suffix_Array.hpp"

#include <string>
#include <cstdlib>
#include <fstream>
#include <iostream>

int main(int argc, char* argv[])
{
    // TODO: standardize the API.
    constexpr auto arg_count = 3;
    if(argc != arg_count)
    {
        std::cerr << "Usage: themis <input_path> <output_path>\n";
        std::exit(EXIT_FAILURE);
    }


    const std::string ip_path(argv[1]);
    const std::string op_path(argv[2]);
    std::string str;

    std::ifstream(ip_path) >> str;

    themis::Suffix_Array suf_arr(str.c_str(), str.length());
    suf_arr.construct();

    std::ofstream output(op_path);

    const std::size_t n = suf_arr.n();
    for(std::size_t i = 0; i < n; ++i)
        output << suf_arr.SA()[i] << " \n"[i == n - 1];
    for(std::size_t i = 0; i < n; ++i)
        output << suf_arr.LCP()[i] << " \n"[i == n - 1];

    output.close();


    return 0;
}
