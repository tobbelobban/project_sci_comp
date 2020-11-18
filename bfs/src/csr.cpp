#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <errno.h>
#include <iostream>
#include <cmath>
#include <queue>
#include <omp.h>
#include <immintrin.h>
#include "csr.hpp"

void read_csr_graph_from_file(const std::string& file_path, csr_graph& csr_g) {
    // assuming that the edges are sorted in file, 32-bit vert ids
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
    int32_t edge_buffer[2];
    int32_t prev_min = 0, prev_max = 0, min, max, count = 0;
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
        count += 2;
        prev_min = min;
        prev_max = max;
    }
    fclose(file_ptr);
    csr_g.cols = new int32_t[count];
    csr_g.rows = new int32_t[nverts+1];
    csr_g.nverts = nverts;
    csr_g.nedges = count;
    count = 0;
    for(int32_t v = 0; v < nverts; ++v) {
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
    for(int32_t from = 0; from < csr_g.nverts; ++from) {
        for(int32_t to = csr_g.rows[from]; to < csr_g.rows[from+1]; ++to) {
            std::cout << csr_g.cols[to] << " ";
        }
        std::cout << std::endl;
    }
}

bool same(const std::vector<int32_t>& d1, const std::vector<int32_t>& d2) {
    const int32_t n = d1.size();
    for(int32_t i = 0; i < n; ++i) {
        if(d1[i] != d2[i]) return false;
    }
    return true;
}

void tropical_iin_csr_mv_mult(std::vector<int32_t>& y, const csr_graph& csr_g, const std::vector<int32_t>& x) {
    for(int32_t i = 0; i < csr_g.nverts; ++i) {
        int64_t j;
        __m256i tmps = _mm256_set_epi32(csr_g.nverts,csr_g.nverts,csr_g.nverts,csr_g.nverts,csr_g.nverts,csr_g.nverts,csr_g.nverts,csr_g.nverts);
        __m256i xs;
        for(j = csr_g.rows[i]; j < (int64_t)csr_g.rows[i+1]-7; j+=8) {
            xs = _mm256_set_epi32(1+x[csr_g.cols[j+0]], 1+x[csr_g.cols[j+1]], 1+x[csr_g.cols[j+2]], 1+x[csr_g.cols[j+3]],1+x[csr_g.cols[j+4]],1+x[csr_g.cols[j+5]],1+x[csr_g.cols[j+6]],1+x[csr_g.cols[j+7]]);
            tmps = _mm256_min_epu32(xs,tmps);
        }
        
        y[i] = std::min( y[i], (int32_t)std::min(
                    (int32_t)std::min(
                    std::min(
                        _mm256_extract_epi32(tmps,0),
                        _mm256_extract_epi32(tmps,1)
                    ),
                    std::min(
                        _mm256_extract_epi32(tmps,2),
                        _mm256_extract_epi32(tmps,3)
                    )),
                    (int32_t)std::min(
                    std::min(
                        _mm256_extract_epi32(tmps,4),
                        _mm256_extract_epi32(tmps,5)
                    ),
                    std::min(
                        _mm256_extract_epi32(tmps,6),
                        _mm256_extract_epi32(tmps,7)
                    ))));
        for(; j < csr_g.rows[i+1]; ++j) {  // clean up last part of loop
            y[i] = std::min(y[i], 1+x[csr_g.cols[j]]);
        }
    }
}

void tropical_csr_mv_mult(std::vector<int32_t>& y, const csr_graph& csr_g, const std::vector<int32_t>& x) {
    for(int32_t i = 0; i < csr_g.nverts; ++i) {
        for(int64_t j = csr_g.rows[i]; j < (int64_t)csr_g.rows[i+1]; ++j) {
            y[i] = std::min(y[i],1+x[csr_g.cols[j]]);
        }
    }
}

std::vector<int32_t> csr_bfs_iin(const csr_graph& csr_g, const int32_t r) {
    std::vector<int32_t> dists(csr_g.nverts, csr_g.nverts);
    if(r >= csr_g.nverts) return dists;
    dists[r] = 0;
    std::vector<int32_t> prev_dists(csr_g.nverts,0);
    while(!same(dists,prev_dists)) {
        prev_dists = dists;
        tropical_iin_csr_mv_mult(dists, csr_g, prev_dists);
    }
    return dists;
}

std::vector<int32_t> csr_bfs(const csr_graph& csr_g, const int32_t r) {
    std::vector<int32_t> dists(csr_g.nverts, csr_g.nverts);
    if(r < 0 || r >= csr_g.nverts) return dists;
    dists[r] = 0;
    std::vector<int32_t> prev_dists(csr_g.nverts,0);
    double total_t = 0, t;
    int32_t its = 0;
    while(!same(dists,prev_dists)) {
        prev_dists = dists;
        t = omp_get_wtime();
        tropical_csr_mv_mult(dists, csr_g, prev_dists);
        total_t += omp_get_wtime() - t;
        ++its;
    }
    std::cout << "avg time CSR = " << total_t/its << std::endl;
    std::cout << "iterations = " << its << std::endl;
    return dists;
}

void print_vector(const std::vector<int32_t>& v) {
    for(uint32_t i = 0; i < v.size(); ++i)
        std::cout << v[i] << " ";
    std::cout << std::endl;
}