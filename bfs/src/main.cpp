#include <iostream>
#include <getopt.h>
#include <cmath>
#include <omp.h>
#include "csr.hpp"
#include "sellcs.hpp"

int main(int argc, char *const *argv) {

    int32_t root = 0, chunk_size=8, sigma=1,SCALE=-1,EDGEFACTOR=-1;
    std::string file_path = "none";
    if(argc < 2) {
        std::cout << "Error. Need filename with graph." << std::endl;
        exit(0);
    }

    // parse arguments
    int opt;
    while ((opt = getopt(argc,argv,":C:s:S:E:r:f:")) != -1)
    {
        switch (opt)
        {
            case 'C':
                chunk_size = std::atoi(optarg);
                break;
            case 's':
                sigma = std::atoi(optarg);            
                break;
            case 'S':
                SCALE = std::atoi(optarg);
                break;
            case 'E':
                EDGEFACTOR = std::atoi(optarg);
                break;
            case 'f':
                file_path = optarg;
                break;
            case 'r':
                root = std::atoi(optarg);
                break;
            default:
                std::cout << "Unrecognized flag: " << opt << "\nTerminating!" << std::endl;
                exit(0);
                break;
        }
    }

    // check arguments
    if(SCALE == -1) {
        std::cout << "ERROR! SCALE must be supplied with: -S X, where X is 2-logarithm of number of vertices." << std::endl;
        exit(0);
    }
    if(EDGEFACTOR == -1) {
        std::cout << "ERROR! EDGEFACTOR must be supplied with: -E X, where X is #edges / #vertices." << std::endl;
        exit(0);
    }
    if(root < 0 || root >= pow(2,SCALE)) {
        std::cout << "ERROR! root = " << root << " is not valid root. Root must be 0 <= root < 2^SCALE. Supply root with: -r X" << std::endl;
        exit(0);
    }
    if(chunk_size != 4 && chunk_size != 8) {
        std::cout << "ERROR! chunk size must be =4 or =8. Modify chunk size with: -C X" << std::endl;
        exit(0);
    }
    if(sigma < 0 || sigma > SCALE) {
        std::cout << "ERROR! sigma must be 0 <= sigma <= SCALE. Modify sigma with: -s X" << std::endl;
        exit(0);
    }
    if(file_path == "none") {
        std::cout << "ERROR! file path containing edges must be specified! Use: -f path/to/file/" << std::endl; 
        exit(0);
    }

    // csr_graph csr_g;
    // // time = omp_get_wtime();
    // read_csr_graph_from_file(file_path, csr_g);
    // // time = omp_get_wtime() - time;
    // // std::cout << "CSR read time: " << time << " s" << std::endl;
    // // //print_csr_graph(csr_g);
    // // time = omp_get_wtime();
    // auto csr_res = csr_bfs(csr_g, root);
    // // time = omp_get_wtime() - time;
    // // std::cout << "CSR solve time: " << time << " s" << std::endl;
    // //print_vector(csr_res);
    // delete_csr(csr_g);

    // std::cout << std::endl;


    sellcs sellcs_g;
    sellcs_g.C = chunk_size;
    sellcs_g.sigma = pow(2,sigma);

    read_sellcs_graph_from_file(file_path, sellcs_g, SCALE, EDGEFACTOR);
    //print_sellcs_graph(sellcs_g);
    auto res = sellcs_bfs(sellcs_g, root);
    // permutate_solution(res, sellcs_g);
    // for(auto i = 0; i < sellcs_g.nverts; ++i) {
    //     if (res[i] != pow(2,SCALE)) {
    //         std::cout << i << " " << res[i]<< std::endl;
    //     }
    // }
    //print_vector(res);
    delete_sellcs(sellcs_g);

    return 0;

}
