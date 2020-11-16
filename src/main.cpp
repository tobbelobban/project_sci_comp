#include<iostream>

#include "make_graph.h"
#include "graph.hpp"
#include "csr.hpp"


int main(int argc, char const *argv[]) {
    if(argc < 2) {
        std::cout << "Error. Need filename with graph." << std::endl;
    }

    graph g;
    read_graph_from_file(argv[1], g);
    sort_graph_edges(g);
    write_graph_to_file(argv[1], g);
    //print_graph_edges(g);
    delete_graph(g);

    csr_graph csr_g;
    read_csr_graph_from_file(argv[1], csr_g);
    print_csr_graph(csr_g);
    return 0;

}
