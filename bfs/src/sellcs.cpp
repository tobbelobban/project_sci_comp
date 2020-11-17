#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <errno.h>
#include <iostream>
#include <cmath>
#include <queue>
#include <cstring>

#include <immintrin.h>
#include "sellcs.hpp"

void read_sellcs_graph_from_file(const std::string& file_path, sellcs& sellcs_g) {
    // assuming that the edges are sorted in row-major order in file, 32-bit vert ids
    auto bottom_dir = file_path.find_last_of('/')+1;
    std::string dir = file_path.substr(0,bottom_dir);
    const char* file = file_path.substr(bottom_dir, file_path.length()).c_str();
    // attempt to open file
    FILE *file_ptr;
    file_ptr = fopen(file_path.c_str(), "rb");
    if(file_ptr == NULL) {
        std::cout << "ERROR! Failed to open file: " << file_path << std::endl;
        exit(0);
    }

    //read SCALE and EDGEFACTOR
    // file_path must be of format "path/to/file/SCALE_EDGEFACTOR_blablabla"
    errno = 0;
    char * end;
    const uint32_t nverts = (uint32_t) pow(2,strtol(file, &end, 10));
    end++;
    const uint32_t nedges = strtol(end, &end, 10) * nverts;
    if(errno == ERANGE) {
        printf("Range error occured while extracting scale and edge factor\n");
        exit(0);
    }
    std::vector<std::queue<uint32_t>> edges(nverts,std::queue<uint32_t>());
    const uint32_t num_chunks = nverts/sellcs_g.C; 
    sellcs_g.nverts = nverts;
    sellcs_g.n_chunks = num_chunks;
    sellcs_g.cl = new uint32_t[num_chunks];
    std::memset(sellcs_g.cl, 0, num_chunks*sizeof(uint32_t));
    uint32_t edge_buffer[2];
    uint32_t row_len = 0;
    uint32_t prev_min = 0, prev_max = 0, min, max;
    for(uint32_t i = 0; i < nedges; ++i) {
        fread(edge_buffer, sizeof(uint32_t), 2, file_ptr);
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
        if(min != prev_min) {
            uint32_t c = prev_min / sellcs_g.C;
            sellcs_g.cl[c] = std::max(sellcs_g.cl[c], row_len);
            row_len = edges[min].size();
        } else {
            ++row_len;
        }
        prev_min = min;
        prev_max = max;
    }
    sellcs_g.cl[num_chunks-1] = std::max(row_len, sellcs_g.cl[num_chunks-1]);
    fclose(file_ptr);

    // for(int i = 0; i < nverts; ++i) {
    //     while(!edges[i].empty()) {
    //         std::cout << edges[i].front() << " ";
    //         edges[i].pop();
    //     }
    //     std::cout << std::endl;
    // }

    uint32_t size = 0;
    sellcs_g.cs = new uint32_t[num_chunks];
    for(uint32_t c = 0; c < num_chunks; ++c) {
        sellcs_g.cs[c] = size;
        size += sellcs_g.cl[c]*sellcs_g.C;
    }
    sellcs_g.cols = new uint32_t[size];
    uint32_t vid;
    for(uint32_t c = 0; c < num_chunks; ++c) {
        for(uint32_t l = 0; l < sellcs_g.cl[c]; ++l) {
            for(uint32_t r = 0; r < sellcs_g.C; ++r) {
                vid = c*sellcs_g.C+r;
                if(edges[vid].empty()) {
                    sellcs_g.cols[sellcs_g.cs[c]+l*sellcs_g.C+r] = 0;    
                } else {
                    sellcs_g.cols[sellcs_g.cs[c]+l*sellcs_g.C+r] = edges[vid].front();
                    edges[vid].pop();
                }
            } 
        }
    }
}

void print_sellcs_graph(const sellcs& sellcs_g) {
    for(uint32_t c = 0; c < sellcs_g.n_chunks; ++c) {
        for(uint32_t r = 0; r < sellcs_g.C; ++r) {
            for(uint32_t l = 0; l < sellcs_g.cl[c]; ++l) {        
                std::cout << sellcs_g.cols[sellcs_g.cs[c]+l*sellcs_g.C+r] << " ";
            } 
            std::cout << std::endl;
        }
    }
}

void delete_sellcs(sellcs& sellcs_g) {
    delete sellcs_g.cols;
    delete sellcs_g.cs;
    delete sellcs_g.cl;
}