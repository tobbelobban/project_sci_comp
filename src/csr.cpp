#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <errno.h>
#include <iostream>
#include <cmath>
#include <queue>
#include "csr.hpp"

void read_csr_graph_from_file(const char* filename, csr_graph& csr_g) {
    // assuming that the edges are sorted in file

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
    std::vector<std::queue<int64_t>> edges(nverts,std::queue<int64_t>());
    int64_t edge_buffer[2];
    int64_t prev_min = 0, prev_max = -1, min, max, count = 0;
    for(int i = 0; i < nedges; i++) {
        fread(edge_buffer, sizeof(int64_t)*2, 1, file_ptr);
        min = edge_buffer[0];
        max = edge_buffer[1];
        if(min == max) continue; // skip self-loops
        if(min > max) {
            min = max;
            max = edge_buffer[0];
        } 
        if(min == prev_min && max == prev_max) continue; // skip duplicates
        edges[min].push(max);
        edges[max].push(min);
        count += 2;
        prev_min = min;
        prev_max = max;
        // std::cout << "Edge " << i << ": " << min << " -- " << max << std::endl;
    }
    fclose(file_ptr);
    csr_g.cols = new int64_t[count];
    csr_g.rows = new int64_t[nverts+1];
    csr_g.nverts = nverts;
    csr_g.nedges = count;
    count = 0;
    for(int v = 0; v < nverts; ++v) {
        csr_g.rows[v] = count;
        while(!edges[v].empty()) {
            csr_g.cols[count++] = edges[v].front();
            edges[v].pop();
        }
    }
    csr_g.rows[nverts] = count;
}

void print_csr_graph(const csr_graph& csr_g) {
    for(int64_t from = 0; from < csr_g.nverts; ++from) {
        for(int64_t to = csr_g.rows[from]; to < csr_g.rows[from+1]; ++to) {
            std::cout << csr_g.cols[to] << " ";
        }
        std::cout << std::endl;
    }
}
