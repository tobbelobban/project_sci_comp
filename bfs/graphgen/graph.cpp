#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <errno.h>
#include <algorithm>
#include <iostream>
#include <cmath>

#include "make_graph.h"
#include "graph.hpp"

void read_graph_from_file(const std::string& f_path, graph & g) {        
    // attempt to open file
    FILE *file_ptr;
    file_ptr = fopen(f_path.c_str(), "rb");
    if(file_ptr == NULL) {
        std::cout << "ERROR! Failed to open file: " << f_path << std::endl;
        exit(0);
    }

    //read SCALE and EDGEFACTOR from file path
    // filename must be of format "path/to/file/SCALE_EDFEFACTOR_blablabla"
    errno = 0;
    auto bottom_dir = f_path.find_last_of('/');
    std::string dir = f_path.substr(0,bottom_dir);
    const char* file = f_path.substr(bottom_dir, f_path.length()).c_str();
    char * end;
    const int64_t nverts = (int) pow(2,strtol(file, &end, 10));
    end++;
    const int64_t nedges = strtol(end, &end, 10) * nverts;
    if(errno == ERANGE) {
        printf("Range error occured while extracting scale and edge factor\n");
        exit(0);
    }
    packed_edge * edges = new packed_edge[nedges];
    int64_t edge_buffer[2];
    for(int i = 0; i < nedges; i++) {
        fread(edge_buffer, sizeof(packed_edge), 1, file_ptr);
        packed_edge pe = {edge_buffer[0], edge_buffer[1]};
        edges[i] = pe;
    }
    fclose(file_ptr);
    g.nedges = nedges;
    g.nverts = nverts;
    g.edges = edges;
}

bool compare_edges(const packed_edge& e1, const packed_edge& e2) {
    int64_t min_e1 = std::min(e1.v0,e1.v1);
    int64_t min_e2 = std::min(e2.v0,e2.v1);
    if(min_e1 == min_e2) {
        return std::max(e1.v0,e1.v1) < std::max(e2.v0,e2.v1);
    }
    return min_e1 < min_e2;
}

void sort_graph_edges(graph& g) {
    std::sort(g.edges, g.edges+g.nedges, compare_edges);
}

void delete_graph(graph& g) {
    free(g.edges);
}

void print_graph_edges(const graph& g) {
    for(int i = 0; i < g.nedges; i++) {
        std::cout << "Edge " << i << ": " << g.edges[i].v0 << " -- " << g.edges[i].v1 << std::endl;
    }
}

void write_graph_to_file(const std::string& f_path, const graph& g) {
    // attempt to open file
    FILE *file_ptr;
    file_ptr = fopen(f_path.c_str(), "wb");
    if(file_ptr == NULL) {
        std::cout << "ERROR! Failed to open file: " << f_path << std::endl;
        exit(0);
    }
    // convert 64-bit nodes to 32-bit
    int32_t buffer[2];
    for(int i = 0; i < g.nedges; ++i) {
        buffer[0] = g.edges[i].v0;
        buffer[1] = g.edges[i].v1;
        size_t f = fwrite(buffer, sizeof(int32_t), 2, file_ptr);
        if(f != 2) {
            std::cout << "Error writing graph... exiting." << std::endl;
            fclose(file_ptr);
            exit(0);
        }
    }
        
    fclose(file_ptr);
}

void generate_graph_to_file(const std::string& f_path, const uint64_t scale, const uint64_t edge_factor, uint64_t s1, uint64_t s2) {
    graph g;
    int64_t nverts = (uint32_t)pow(2,scale);
    int64_t nedges = edge_factor * nverts;
    make_graph(scale, nedges, s1, s2, &g.nedges, &g.edges);
    sort_graph_edges(g);
    write_graph_to_file(f_path, g);
    delete_graph(g);
}