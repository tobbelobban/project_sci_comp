#include <iostream>
#include <omp.h>
#include "csr.hpp"
#include "sellcs.hpp"

int main(int argc, char const *argv[]) {
    if(argc < 2) {
        std::cout << "Error. Need filename with graph." << std::endl;
    }

    std::string file_path(argv[1]);

    csr_graph csr_g;
    read_csr_graph_from_file(file_path, csr_g);
    print_csr_graph(csr_g);
    delete_csr(csr_g);

    std::cout << std::endl;

    sellcs sellcs_g;
    sellcs_g.C = 2;
    read_sellcs_graph_from_file(file_path, sellcs_g);
    print_sellcs_graph(sellcs_g);
    delete_sellcs(sellcs_g);
    
    std::cout << std::endl;
    


    return 0;

}
