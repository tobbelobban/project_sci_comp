#include <iostream>
#include <getopt.h>
#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <string>

#include "make_graph.h"
#include "graph.hpp"

int main(int argc, char* const* argv) {
    uint64_t scale = 10, edge_factor = 16, s1 = 2, s2 = 32;

    std::string f_path;
    int option, tmp;
    bool in_fn = false;
    while((option = getopt(argc, argv, "s:e:f:")) != -1) {
        switch (option)
        {
        case 's':
            tmp = atoi(optarg);
            scale = tmp > 0 ? tmp : scale;
            break;
        case 'e':
            tmp = atoi(optarg);
            edge_factor = tmp > 0 ? tmp : edge_factor;
            break;
        case 'f':
            f_path = optarg;
            in_fn = true;
            break;
        default:
            std::cout << "Unrecognized argument: " << option << std::endl;
            break;
        }
    }
    if(!in_fn)
        f_path = "../bin/" + std::to_string(scale) + '_' + std::to_string(edge_factor) + ".bin"; 
    std::cout << "GENERATING GRAPH TO FILE: " << f_path << std::endl;
    generate_graph_to_file(f_path, scale, edge_factor,s1,s2);
    std::cout << "--GENERATION DONE!--" << std::endl;
    return 0;
}
