#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <errno.h>
#include <algorithm>
#include <iostream>
#include <cmath>

#include "make_graph.h"
#include "graph.hpp"

void read_graph_from_file(const char * filename, graph & g) {        
    // attempt to open file
    FILE *file_ptr;
    file_ptr = fopen(filename, "rb");
    if(file_ptr == NULL) {
        printf("FAILED TO OPEN FILE: %s\n", filename);
        exit(0);
    }

    //read SCALE and EDGEFACTOR
    // filename must be of format "SCALE_EDGEFACTOR_*"
    errno = 0;
    char * end;
    const int64_t nverts = (int) pow(2,strtol(filename, &end, 10));
    end++;
    const int16_t nedges = strtol(end, &end, 10) * nverts;
    if(errno == ERANGE) {
        printf("Range error occured while extracting scale and edge factor\n");
        exit(0);
    }
    packed_edge * edges = new packed_edge[nedges];
    int64_t edge_buffer[2];
    for(int i = 0; i < nedges; i++) {
        fread(edge_buffer, 16, 1, file_ptr);
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

void write_graph_to_file(const char* filename, const graph& g) {
    // attempt to open file
    FILE *file_ptr;
    file_ptr = fopen(filename, "wb");
    if(file_ptr == NULL) {
        printf("FAILED TO OPEN FILE: %s\n", filename);
        exit(0);
    }
    fwrite(g.edges, sizeof(packed_edge)*g.nedges, 1, file_ptr);
    fclose(file_ptr);
}