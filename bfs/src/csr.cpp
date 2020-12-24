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

void read_csr_graph_from_file(const std::string& file_path, csr_graph& csr_g, int32_t SCALE, int32_t EDGEFACTOR, double* const timer) {

    // assuming that the edges are sorted in file, 32-bit vert ids
    auto bottom_dir = file_path.find_last_of('/')+1;
    std::string dir = file_path.substr(0,bottom_dir);

    // set number of vertices and edges
    const int32_t nverts = (int32_t) pow(2,SCALE);
    const int32_t nedges = EDGEFACTOR * nverts;

    // queue for maintaining edges 
    //std::vector<std::queue<int32_t>> edges(nverts,std::queue<int32_t>());
    
    // attempt to open file
    FILE *file_ptr;
    file_ptr = fopen(file_path.c_str(), "rb");
    if(file_ptr == NULL) {
        std::cout << "ERROR! Failed to open file: " << file_path << std::endl;
        exit(0);
    }

    // buffer to read into
    int32_t* edge_buffer = new int32_t[nedges*2];
    int32_t* degrees = new int32_t[nverts];
    memset((void*)degrees, 0, sizeof(int32_t)*nverts);

    // read edges to buffer
    timer[0] = cpuSecond();
    auto read_res = fread(edge_buffer, 2*sizeof(int32_t), nedges, file_ptr);
    timer[0] = cpuSecond() - timer[0];
    fclose(file_ptr);
    if(read_res != (size_t)nedges) 
    {
        std::cout << "Error reading edges from file... exiting!" << std::endl;
        exit(0);
    }

    // process edges
    timer[1] = cpuSecond();
    int32_t prev_min = 0, prev_max = 0, min, max, e_count = 0;
    for(int32_t e = 0; e < nedges; ++e) {
        min = edge_buffer[e*2];
        max = edge_buffer[e*2+1];
        
        if(min > max)
        {
            min = max;
            max = edge_buffer[e*2];
        }

        // skip self-loops and duplicate edges
        if(min == max || (min == prev_min && max == prev_max)) 
        {
            edge_buffer[e*2] = -1; // if this edge is not valid, -1 will tell us so later    
            continue; 
        }

        ++degrees[min];
        ++degrees[max];
        e_count += 2;
        prev_min = min;
        prev_max = max;
    }
    timer[1] = cpuSecond() - timer[1];
    
    // set edges in CSR
    timer[2] = cpuSecond();
    
    csr_g.cols = new int32_t[e_count];
    csr_g.rows = new int32_t[nverts+1];
    csr_g.nverts = nverts;
    csr_g.nedges = e_count;
    
    int32_t* offsets = new int32_t[nverts];
    int32_t offset = 0;
    for(int32_t v = 0; v < nverts; ++v)
    {
        offsets[v] = offset;
        csr_g.rows[v] = offset;
        offset += degrees[v];
    }
    csr_g.rows[nverts] = offset;

    int32_t v1,v2;
    for(int32_t e = 0; e < nedges; ++e) {
        v1 = edge_buffer[e*2];
        if(v1 == -1) continue;
        v2 = edge_buffer[2*e+1];
        csr_g.cols[offsets[v1]++] = v2;
        csr_g.cols[offsets[v2]++] = v1;
    }
    timer[2] = cpuSecond() - timer[2];
    
    delete[] edge_buffer;
    delete[] offsets;
}

void delete_csr(csr_graph& csr_g) {
    delete[] csr_g.cols;
    delete[] csr_g.rows;
}

void print_csr_graph(const csr_graph& csr_g) {
    for(int32_t from = 0; from < csr_g.nverts; ++from) {
        for(int32_t to = csr_g.rows[from]; to < csr_g.rows[from+1]; ++to) {
            std::cout << csr_g.cols[to] << "\t";
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


void tropical_csr_mv_mult(std::vector<int32_t>& y, const csr_graph& csr_g, const std::vector<int32_t>& x) {
    #pragma omp parallel for
    for(int32_t i = 0; i < csr_g.nverts; ++i) {
        for(int64_t j = csr_g.rows[i]; j < (int64_t)csr_g.rows[i+1]; ++j) {
            y[i] = std::min(y[i],1+x[csr_g.cols[j]]);
        }
    }
}

int32_t csr_bfs(const csr_graph& csr_g, const int32_t r, std::vector<int32_t>& res) {
    // res must be initialzed s.t. res[r] = 0, res[v != r] = nverts
    if(r < 0 || r >= csr_g.nverts) return -1;
    std::vector<int32_t> prev_dists(csr_g.nverts,0);
    int32_t its = 0;
    while(!same(res,prev_dists)) {
        prev_dists = res;
        tropical_csr_mv_mult(res, csr_g, prev_dists);
        ++its;
    }
    return its;
}

void print_vector(const std::vector<int32_t>& v) {
    for(uint32_t i = 0; i < v.size(); ++i)
        std::cout << v[i] << " ";
    std::cout << std::endl;
}

int isisolated(int32_t root, const csr_graph& g){
    return g.rows[root] == g.rows[root+1];
}