#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <errno.h>
#include <iostream>
#include <cmath>
#include <queue>

#include <immintrin.h>
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
    const int64_t nverts = (int64_t) pow(2,strtol(filename, &end, 10));
    end++;
    const int64_t nedges = strtol(end, &end, 10) * nverts;
    if(errno == ERANGE) {
        printf("Range error occured while extracting scale and edge factor\n");
        exit(0);
    }
    std::vector<std::queue<int64_t>> edges(nverts,std::queue<int64_t>());
    int64_t edge_buffer[2];
    int64_t prev_min = 0, prev_max = -1, min, max, count = 0;
    for(int i = 0; i < nedges; ++i) {
        fread(edge_buffer, sizeof(int64_t), 2, file_ptr);
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
    for(int64_t v = 0; v < nverts; ++v) {
        csr_g.rows[v] = count;
        while(!edges[v].empty()) {
            csr_g.cols[count++] = edges[v].front();
            edges[v].pop();
        }
    }
    csr_g.rows[nverts] = count;
}

void delete_csr(csr_graph& csr_g) {
    delete csr_g.cols;
    delete csr_g.rows;
}

void print_csr_graph(const csr_graph& csr_g) {
    for(int64_t from = 0; from < csr_g.nverts; ++from) {
        for(int64_t to = csr_g.rows[from]; to < csr_g.rows[from+1]; ++to) {
            std::cout << csr_g.cols[to] << " ";
        }
        std::cout << std::endl;
    }
}

bool same(const std::vector<int64_t>& d1, const std::vector<int64_t>& d2) {
    const int64_t n = d1.size();
    for(int64_t i = 0; i < n; ++i) {
        if(d1[i] != d2[i]) return false;
    }
    return true;
}

// void tropical_iin_csr_mv_mult(std::vector<int64_t>& y, const csr_graph& csr_g, const std::vector<int64_t>& x) {
//     for(int64_t i = 0; i < csr_g.nverts; ++i) {
//         int64_t j;
//         __m256i tmps = _mm256_set_epi32(csr_g.nverts,csr_g.nverts,csr_g.nverts,csr_g.nverts,0,0,0,0);
//         __m256i xs;
//         for(j = csr_g.rows[i]; j < csr_g.rows[i+1]-3; j+=4) {
//             xs = _mm256_set_epi32(1+x[csr_g.cols[j+0]], 1+x[csr_g.cols[j+1]], 1+x[csr_g.cols[j+2]], 1+x[csr_g.cols[j+3]],0,0,0,0);
//             tmps = _mm256_min_epu32(xs,tmps);
//         }
        
//         y[i] = std::min( y[i], (int64_t)std::min(
//                     std::min(
//                         _mm256_extract_epi32(tmps,4),
//                         _mm256_extract_epi32(tmps,5)
//                     ),
//                     std::min(
//                         _mm256_extract_epi32(tmps,7),
//                         _mm256_extract_epi32(tmps,6)
//                     )));
//         for(; j < csr_g.rows[i+1]; ++j) {  // clean up last part of loop
//             y[i] = std::min(y[i], 1+x[csr_g.cols[j]]);
//         }
//     }
// }

void tropical_csr_mv_mult(std::vector<int64_t>& y, const csr_graph& csr_g, const std::vector<int64_t>& x) {
    int64_t tmp0, tmp1, tmp2, tmp3;
    for(int64_t i = 0; i < csr_g.nverts; ++i) {
        tmp0 = tmp1 = tmp2 = tmp3 = csr_g.nverts;
        int64_t j;
        for(j = csr_g.rows[i]; j < csr_g.rows[i+1]-3; j+=4) {
            tmp0 = std::min(tmp0, 1+x[csr_g.cols[j+0]]);
            tmp1 = std::min(tmp1, 1+x[csr_g.cols[j+1]]);
            tmp2 = std::min(tmp2, 1+x[csr_g.cols[j+2]]);
            tmp3 = std::min(tmp3, 1+x[csr_g.cols[j+3]]);
            
        }
        y[i] = std::min(y[i],std::min(std::min(tmp0,tmp1), std::min(tmp2,tmp3))); // reduce
        for(; j < csr_g.rows[i+1]; ++j) {  // clean up last part of loop
            y[i] = std::min(y[i], 1+x[csr_g.cols[j]]);
        }
    }
}

// std::vector<int64_t> csr_bfs_iin(const csr_graph& csr_g, const int64_t r) {
//     std::vector<int64_t> dists(csr_g.nverts, csr_g.nverts);
//     if(r < 0 || r >= csr_g.nverts) return dists;
//     dists[r] = 0;
//     std::vector<int64_t> prev_dists(csr_g.nverts,0);
//     while(!same(dists,prev_dists)) {
//         prev_dists = dists;
//         tropical_iin_csr_mv_mult(dists, csr_g, prev_dists);
//     }
//     return dists;
// }

std::vector<int64_t> csr_bfs(const csr_graph& csr_g, const int64_t r) {
    std::vector<int64_t> dists(csr_g.nverts, csr_g.nverts);
    if(r < 0 || r >= csr_g.nverts) return dists;
    dists[r] = 0;
    std::vector<int64_t> prev_dists(csr_g.nverts,0);
    while(!same(dists,prev_dists)) {
        prev_dists = dists;
        tropical_csr_mv_mult(dists, csr_g, prev_dists);
    }
    return dists;
}

void print_vector(const std::vector<int64_t>& v) {
    for(uint i = 0; i < v.size(); ++i)
        std::cout << v[i];
    std::cout << std::endl;
}