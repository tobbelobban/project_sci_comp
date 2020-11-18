#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <errno.h>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <queue>
#include <cstring>
#include <omp.h>
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
    std::vector<std::queue<int32_t>> edges(nverts, std::queue<int32_t>());
    std::vector<vertex> vertex_degs(nverts, vertex());
    for(int32_t vid = 0; vid < nverts; ++vid) vertex_degs[vid].vid = vid;
    const int32_t num_chunks = nverts/sellcs_g.C; 
    sellcs_g.nverts = nverts;
    sellcs_g.n_chunks = num_chunks;
    sellcs_g.cl = new int32_t[num_chunks];
    std::memset(sellcs_g.cl, 0, num_chunks*sizeof(int32_t));
    int32_t edge_buffer[2];
    int32_t prev_min = 0, prev_max = 0, min, max;
    int64_t e_count = 0;
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
        e_count += 2;
        prev_min = min;
        prev_max = max;
    }
    fclose(file_ptr);
    for(int32_t vid = 0; vid < nverts; ++vid) vertex_degs[vid].degree = edges[vid].size();
    for(int32_t i = 0; i < nverts/sellcs_g.sigma; ++i) {
        int32_t offs_front = i*sellcs_g.sigma;
        int32_t offs_back =  (i+1)*sellcs_g.sigma;
        std::partial_sort(vertex_degs.begin()+offs_front,vertex_degs.begin()+offs_back, vertex_degs.begin()+offs_back,[](const vertex& v1, const vertex& v2){return v1.degree > v2.degree;});
    }
    //std::sort(vertex_degs.begin(), vertex_degs.end(), [](const vertex& v1, const vertex& v2){return v1.degree > v2.degree;});
    sellcs_g.permuts = new int32_t[nverts];
    for(int32_t vid = 0; vid < nverts; ++vid) sellcs_g.permuts[vertex_degs[vid].vid] = vid;
    int32_t size = 0;
    sellcs_g.cs = new int32_t[num_chunks];
    for(int32_t c = 0; c < num_chunks; ++c) {
        sellcs_g.cl[c] = vertex_degs[c*sellcs_g.C].degree;
        sellcs_g.cs[c] = size;
        size += sellcs_g.cl[c]*sellcs_g.C;
    }
    sellcs_g.cols = new int32_t[size];
    sellcs_g.beta = (double)e_count / (double)size;
    int32_t vid;
    for(int32_t c = 0; c < num_chunks; ++c) {
        for(int32_t l = 0; l < sellcs_g.cl[c]; ++l) {
            for(int32_t r = 0; r < sellcs_g.C; ++r) {
                vid = vertex_degs[c*sellcs_g.C+r].vid;
                if(edges[vid].empty()) {
                    sellcs_g.cols[sellcs_g.cs[c]+l*sellcs_g.C+r] = -1;    
                } else {
                    sellcs_g.cols[sellcs_g.cs[c]+l*sellcs_g.C+r] = sellcs_g.permuts[edges[vid].front()];
                    edges[vid].pop();
                }
            } 
        }
    }
}

void tropical_sellcs_mv_mult(std::vector<int32_t>& y, const sellcs& g, const std::vector<int32_t>& x) {
    __m128i ones = _mm_set_epi32(1,1,1,1), m_ones = _mm_set_epi32(-1,-1,-1,-1), infs = _mm_set_epi32(g.nverts, g.nverts, g.nverts, g.nverts);
    __m128i tmps, col, vals, rhs;
    int32_t c_offs;
    for(int32_t i = 0; i < g.n_chunks; ++i) {
        tmps = _mm_load_si128((__m128i*)&x[i*g.C]); // load chunk from frontier (x)
        c_offs = g.cs[i];
        for(int32_t j = 0; j < g.cl[i]; ++j) {
            col = _mm_load_si128((__m128i*)&g.cols[c_offs]);
            vals = _mm_cmpeq_epi32(m_ones,col);
            vals = _mm_blendv_epi8(ones,infs,vals);
            rhs = _mm_set_epi32(x[g.cols[c_offs+3]], x[g.cols[c_offs+2]], x[g.cols[c_offs+1]], x[g.cols[c_offs+0]]);
            tmps = _mm_min_epi32(_mm_add_epi32(rhs,vals),tmps);
            c_offs += g.C;
        }
        _mm_store_si128((__m128i*)&y[i*g.C], tmps);
    }
}

std::vector<int32_t> sellcs_bfs(const sellcs& g, const int32_t r) {
    std::vector<int32_t> dists(g.nverts, g.nverts);
    if(r >= g.nverts) return dists;
    dists[r] = 0;
    double total_t = 0, t;
    std::vector<int32_t> prev_dists(g.nverts,0);
    int32_t its = 0;
    while(!same(dists,prev_dists)) {
        prev_dists = dists;
        t = omp_get_wtime();
        tropical_sellcs_mv_mult(dists, g, prev_dists);
        total_t += omp_get_wtime() - t;
        ++its;
    }
    std::cout << "avg time SELL-C-S = " << total_t/its << std::endl;
    std::cout << "iterations = " << its << std::endl;
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

void permutate_solution(std::vector<int32_t>& solv, const sellcs& g) {
    std::vector<int32_t> solv_cp(solv);
    for(int32_t oldp_v = 0; oldp_v < g.nverts; ++oldp_v) {
        solv[oldp_v] = solv_cp[g.permuts[oldp_v]];
    }
}

int32_t get_permutated_vid(const int32_t vid, sellcs& g) {
    if(vid < 0 || vid > g.nverts) return -1;
    return g.permuts[vid];
} 

void delete_sellcs(sellcs& sellcs_g) {
    delete sellcs_g.cols;
    delete sellcs_g.cs;
    delete sellcs_g.cl;
    delete sellcs_g.permuts;
}