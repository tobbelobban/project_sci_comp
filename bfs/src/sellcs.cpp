#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <errno.h>
#include <iostream>
#include <cmath>
#include <queue>
#include <cstring>

#include <immintrin.h>
#include "csr.hpp"
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
    const int32_t nverts = (int32_t) pow(2,strtol(file, &end, 10));
    end++;
    const int32_t nedges = strtol(end, &end, 10) * nverts;
    if(errno == ERANGE) {
        printf("Range error occured while extracting scale and edge factor\n");
        exit(0);
    }
    std::vector<std::queue<int32_t>> edges(nverts,std::queue<int32_t>());
    const int32_t num_chunks = nverts/sellcs_g.C; 
    sellcs_g.nverts = nverts;
    sellcs_g.n_chunks = num_chunks;
    sellcs_g.cl = new int32_t[num_chunks];
    std::memset(sellcs_g.cl, 0, num_chunks*sizeof(int32_t));
    int32_t edge_buffer[2];
    int32_t prev_min = 0, prev_max = 0, min, max;
    for(int32_t i = 0; i < nedges; ++i) {
        fread(edge_buffer, sizeof(int32_t), 2, file_ptr);
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
        sellcs_g.cl[min/sellcs_g.C] = std::max(sellcs_g.cl[min/sellcs_g.C],(int32_t)edges[min].size());
        sellcs_g.cl[max/sellcs_g.C] = std::max(sellcs_g.cl[max/sellcs_g.C],(int32_t)edges[max].size());
        prev_min = min;
        prev_max = max;
    }
    fclose(file_ptr);
    
    int32_t size = 0;
    sellcs_g.cs = new int32_t[num_chunks];
    for(int32_t c = 0; c < num_chunks; ++c) {
        sellcs_g.cs[c] = size;
        size += sellcs_g.cl[c]*sellcs_g.C;
    }
    sellcs_g.cols = new int32_t[size];
    int32_t vid;

    for(int32_t c = 0; c < num_chunks; ++c) {
        for(int32_t l = 0; l < sellcs_g.cl[c]; ++l) {
            for(int32_t r = 0; r < sellcs_g.C; ++r) {
                vid = c*sellcs_g.C+r;
                if(edges[vid].empty()) {
                    sellcs_g.cols[sellcs_g.cs[c]+l*sellcs_g.C+r] = -1;    
                } else {
                    sellcs_g.cols[sellcs_g.cs[c]+l*sellcs_g.C+r] = edges[vid].front();
                    edges[vid].pop();
                }
            } 
        }
    }
}

void tropical_sellcs_mv_mult(std::vector<int32_t>& y, const sellcs& g, const std::vector<int32_t>& x) {
    // for(int i = 0; i < g.n_chunks; ++i) {
    //     for(int j = 0; j < g.cl[i]; ++j) {
    //         y[i*g.C+0] = std::min(y[i*g.C+0], );
    //         y[i*g.C+1] = std::min(y[i*g.C+1],);
    //         y[i*g.C+2] = std::min(y[i*g.C+2],);
    //         y[i*g.C+3] = std::min(y[i*g.C+3],);
    //     }
    // }
}

std::vector<int32_t> sellcs_bfs(const sellcs& g, const int32_t r) {
    std::vector<int32_t> dists(g.nverts, g.nverts);
    if(r >= g.nverts) return dists;
    dists[r] = 0;
    std::vector<int32_t> prev_dists(g.nverts,0);
    while(!same(dists,prev_dists)) {
        prev_dists = dists;
        tropical_sellcs_mv_mult(dists, g, prev_dists);
    }
    return dists;
}

void print_sellcs_graph(const sellcs& sellcs_g) {
    for(int32_t c = 0; c < sellcs_g.n_chunks; ++c) {
        for(int32_t r = 0; r < sellcs_g.C; ++r) {
            for(int32_t l = 0; l < sellcs_g.cl[c]; ++l) {        
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