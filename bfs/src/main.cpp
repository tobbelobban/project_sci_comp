#include <iostream>
#include <omp.h>
#include "csr.hpp"
#include "sellcs.hpp"

int main(int argc, char const *argv[]) {
    if(argc < 2) {
        std::cout << "Error. Need filename with graph." << std::endl;
    }
    int32_t root = 2;
    std::string file_path(argv[1]);
    double time;

    csr_graph csr_g;
    time = omp_get_wtime();
    read_csr_graph_from_file(file_path, csr_g);
    time = omp_get_wtime() - time;
    std::cout << "CSR read time: " << time << " s" << std::endl;
    //print_csr_graph(csr_g);
    time = omp_get_wtime();
    auto csr_res = csr_bfs(csr_g, root);
    time = omp_get_wtime() - time;
    std::cout << "CSR solve time: " << time << " s" << std::endl;
    //print_vector(csr_res);
    delete_csr(csr_g);

    std::cout << std::endl;

    sellcs sellcs_g;
    sellcs_g.C = 4;
    sellcs_g.sigma = 1024;
    time = omp_get_wtime();
    read_sellcs_graph_from_file(file_path, sellcs_g);
    time = omp_get_wtime() - time;
    std::cout << "SELL-C-sigma read time: " << time << " s" << std::endl;
    std::cout << "SELL-C-sigma beta = " << sellcs_g.beta << std::endl;
    //print_sellcs_graph(sellcs_g);
    time = omp_get_wtime();
    auto sellcs_res = sellcs_bfs(sellcs_g, get_permutated_vid(root, sellcs_g));
    time = omp_get_wtime() - time;
    std::cout << "SELL-" << sellcs_g.C << '-' << sellcs_g.sigma << " solve time: " << time << " s" << std::endl;
    permutate_solution(sellcs_res, sellcs_g);
    //print_vector(sellcs_res);
    std::cout << "Same solution? " << (same(sellcs_res, csr_res) ? "yes" : "no") << std::endl;
    delete_sellcs(sellcs_g);

    return 0;

}
