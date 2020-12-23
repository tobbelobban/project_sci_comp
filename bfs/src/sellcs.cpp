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
#include <sys/time.h>

void local_sort(std::vector<vertex> &v, const sellcs &g) {
    std::vector<vertex> tmp(v);
    #pragma omp parallel for 
    for(int32_t i = 0; i < g.nverts/g.sigma; ++i) 
    {
        int32_t start = i*g.sigma;
        int32_t end =  (i+1)*g.sigma;
        sorter(tmp, start, end, v);
    }
}

void merger(std::vector<vertex>& v_tmp, int32_t start, int32_t mid, int32_t end, std::vector<vertex>& v_sorted) {
    int32_t i = start, j = mid;
    for(int32_t k = start; k < end; ++k) {
        if(i < mid && (j >= end || v_tmp[i].degree > v_tmp[j].degree)) {
            v_sorted[k] = v_tmp[i++];
        } else {
            v_sorted[k] = v_tmp[j++];
        }
    }
}

void sorter_serial(std::vector<vertex>& tmp, int32_t start, int32_t end, std::vector<vertex>& v) {
    if(end - start <= 1) {
        return;
    }
    int32_t mid = (start+end)/2;
    sorter(v, start, mid, tmp);
    sorter(v, mid, end, tmp);
    merger(tmp, start, mid, end, v);
}

void sorter(std::vector<vertex>& tmp, int32_t start, int32_t end, std::vector<vertex>& v) {
    if(end - start <= 1024) {
        sorter_serial(tmp, start, end, v);
        return;
    }
    int32_t mid = (start+end)/2;

    #pragma omp parallel
    {
        #pragma omp single nowait 
        {
            #pragma omp task shared(tmp, v, start, mid)
            sorter(v, start, mid, tmp);

            #pragma omp task shared(tmp, v, mid, end)
            sorter(v, mid, end, tmp);
            
            #pragma omp taskwait
            merger(tmp, start, mid, end, v);
        }
    }
}

void read_sellcs_graph_from_file(const std::string& file_path, sellcs& sellcs_g, int32_t SCALE, int32_t EDGEFACTOR, double* const timer) {
    
    // assuming that the edges are stored in row-major order in file, 32-bit vertex ids
    auto bottom_dir = file_path.find_last_of('/')+1;
    std::string dir = file_path.substr(0,bottom_dir);
    
    // set number of vertices and edges
    const int32_t nverts = (int32_t) pow(2,SCALE);
    const int32_t nedges = EDGEFACTOR * nverts;
    
    // vector for counting degree of each vertex
    std::vector<vertex> vertex_degs(nverts, vertex());
    for(int32_t vid = 0; vid < nverts; ++vid) 
    {
        vertex_degs[vid].vid = vid;
        vertex_degs[vid].degree = 0;
    }
    
    const int32_t num_chunks = nverts/sellcs_g.C; 

    // initialze sell-C-sigma members
    sellcs_g.nverts = nverts;
    sellcs_g.n_chunks = num_chunks;
    sellcs_g.cl = new int32_t[num_chunks];    
    sellcs_g.permuts = new int32_t[nverts];
    sellcs_g.cs = new int32_t[num_chunks];
    
    // attempt to open file
    FILE *file_ptr;
    file_ptr = fopen(file_path.c_str(), "rb");
    if(file_ptr == NULL) {
        std::cout << "ERROR! Failed to open file: " << file_path << std::endl;
        exit(0);
    }

    // prepare for reading edges
    int32_t* edge_buffer = new int32_t[nedges*2];

    // read edges to buffer
    timer[0]  = cpuSecond();
    if(fread(edge_buffer, 2*sizeof(int32_t), nedges, file_ptr) != (size_t)nedges) 
    {
        std::cout << "Error reading edges from file... exiting!" << std::endl;
        fclose(file_ptr);
        exit(0);
    }
    fclose(file_ptr);
    timer[0] = cpuSecond() - timer[0];

    // process edges
    timer[1] = cpuSecond();
    int32_t prev_min = 0, prev_max = 0, min, max, e_count = 0;
    for(int32_t e = 0; e < nedges; ++e) 
    {
        // get the vertices of the current edge
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
            edge_buffer[e*2] = -1; // mark the edge as invalid
            continue; 
        }
         
        vertex_degs[min].degree += 1;
        vertex_degs[max].degree += 1;
        e_count += 2;
        prev_min = min;
        prev_max = max;
    }
    timer[1] = cpuSecond() - timer[1];
    
    // sort vertices based on degree (descending order)
    timer[2] = cpuSecond(); 
    local_sort(vertex_degs, sellcs_g);
    timer[2] = cpuSecond() - timer[2];
    
    // compute the chunk lengths and size
    int32_t size=0, offset=0, v=0, currlen;
    int32_t* edge_offsets = new int32_t[nverts];

    for(int32_t c = 0; c < num_chunks; ++c)
    {
        int32_t maxlen = -1;
        for(int32_t r = 0; r < sellcs_g.C; ++r)
        {
            // store the permutated vertex mapping: permuts(original order) --> permutated order
            sellcs_g.permuts[vertex_degs[v].vid] = v;            
            currlen = vertex_degs[v].degree;
            edge_offsets[v] = offset;
            offset += currlen;
            if(currlen > maxlen) maxlen = currlen;
            ++v;
        }
        
        // set chunk length and size
        sellcs_g.cl[c] = maxlen;
        sellcs_g.cs[c] = size;
        size += maxlen*sellcs_g.C;
    }
    
    // store all valid edges in row-major order for permutated order
    int32_t* edges = new int32_t[e_count];
    int32_t v1,v2;
    for(int32_t e = 0; e < nedges; ++e)
    {
        v1 = edge_buffer[e*2];
        if(v1 == -1) continue;
        v1 = sellcs_g.permuts[v1];
        v2 = sellcs_g.permuts[edge_buffer[e*2+1]];
        edges[edge_offsets[v1]++] = v2;
        edges[edge_offsets[v2]++] = v1;
    }
    delete[] edge_buffer;
    
    // initialze SELL-C-Sigma col-array
    sellcs_g.cols = new int32_t[size];
    sellcs_g.beta = (double)e_count / (double)size;
    
    // store edges
    int32_t remaining;
    offset = 0;
    timer[3] = cpuSecond();
    for(int32_t c = 0; c < num_chunks; ++c) 
    {
        for(int32_t l = 0; l < sellcs_g.cl[c]; ++l) 
        {
            for(int32_t r = 0; r < sellcs_g.C; ++r) 
            {
                v = c*sellcs_g.C+r;  
                remaining = vertex_degs[v].degree;
                // if current vertex has no more edges, insert -1
                if(remaining == 0) {
                    sellcs_g.cols[offset] = -1;    
                } else {
                    sellcs_g.cols[offset] = edges[edge_offsets[v]-remaining];
                    --vertex_degs[c*sellcs_g.C+r].degree;
                }
                ++offset;
            } 
        }
    }
    timer[3] = cpuSecond() - timer[3];
    delete[] edges;
    delete[] edge_offsets;
    
}

void tropical_sellcs_mv_mult_w8(std::vector<int32_t>& y, const sellcs& g, const std::vector<int32_t>& x) {
    __m256i ones = _mm256_set1_epi32(1), m_ones = _mm256_set1_epi32(-1), infs = _mm256_set1_epi32(g.nverts);
    #pragma omp parallel for
    for(int32_t i = 0; i < g.n_chunks; ++i) {
        __m256i tmps, col, vals, rhs;
        int32_t c_offs;
        tmps = _mm256_loadu_si256((__m256i*)&x[i*g.C]); // load chunk from frontier (x)
        c_offs = g.cs[i];
        for(int32_t j = 0; j < g.cl[i]; ++j) {
            col = _mm256_loadu_si256((__m256i*)&g.cols[c_offs]);
            vals = _mm256_cmpeq_epi32(m_ones,col);
            vals = _mm256_blendv_epi8(ones,infs,vals);
            rhs = _mm256_set_epi32(x[g.cols[c_offs+7]],x[g.cols[c_offs+6]],x[g.cols[c_offs+5]],x[g.cols[c_offs+4]],x[g.cols[c_offs+3]], x[g.cols[c_offs+2]], x[g.cols[c_offs+1]], x[g.cols[c_offs+0]]);
            tmps = _mm256_min_epi32(_mm256_add_epi32(rhs,vals),tmps);
            c_offs += g.C;
        }
        _mm256_storeu_si256((__m256i*)&y[i*g.C], tmps);
    }
}

void tropical_sellcs_mv_mult_w4(std::vector<int32_t>& y, const sellcs& g, const std::vector<int32_t>& x) {
    __m128i ones = _mm_set1_epi32(1), m_ones = _mm_set1_epi32(-1), infs = _mm_set1_epi32(g.nverts);
    #pragma omp parallel for
    for(int32_t i = 0; i < g.n_chunks; ++i) {
        __m128i tmps, col, vals, rhs;
        int32_t c_offs;
        tmps = _mm_loadu_si128((__m128i*)&x[i*g.C]); // load chunk from frontier (x)
        c_offs = g.cs[i];
        for(int32_t j = 0; j < g.cl[i]; ++j) {
            col = _mm_loadu_si128((__m128i*)&g.cols[c_offs]);
            vals = _mm_cmpeq_epi32(m_ones,col);
            vals = _mm_blendv_epi8(ones,infs,vals);
            rhs = _mm_set_epi32(x[g.cols[c_offs+3]], x[g.cols[c_offs+2]], x[g.cols[c_offs+1]], x[g.cols[c_offs+0]]);
            tmps = _mm_min_epi32(_mm_add_epi32(rhs,vals),tmps);
            c_offs += g.C;
        }
        _mm_storeu_si128((__m128i*)&y[i*g.C], tmps);
    }
}

int32_t sellcs_bfs(const sellcs& g, const int32_t r, std::vector<int32_t>& res) {
    if(r >= g.nverts || r < 0) return -1;
    std::vector<int32_t> prev_dists(g.nverts,0);
    int32_t its = 0;
    if(g.C == 8) {
        while(!same(res,prev_dists)) {
            prev_dists = res;
            tropical_sellcs_mv_mult_w8(res, g, prev_dists);
            ++its;
        }
    } else {
        while(!same(res,prev_dists)) {
            prev_dists = res;
            tropical_sellcs_mv_mult_w4(res, g, prev_dists);
            ++its;
        }
    }
    return its;
}

void print_sellcs_graph(const sellcs& sellcs_g) {
    for(int32_t c = 0; c < sellcs_g.n_chunks; ++c) {
        for(int32_t r = 0; r < sellcs_g.C; ++r) {
            for(int32_t l = 0; l < sellcs_g.cl[c]; ++l) {        
                std::cout << sellcs_g.cols[sellcs_g.cs[c]+l*sellcs_g.C+r] << "\t";
            } 
            std::cout << std::endl;
        }
    }
}

void permutate_solution(std::vector<int32_t>& solution, const sellcs& g) {
    std::vector<int32_t> solv_cp(solution);
    for(int32_t oldp_v = 0; oldp_v < g.nverts; ++oldp_v) {
        solution[oldp_v] = solv_cp[g.permuts[oldp_v]];
    }
}

int32_t get_permutated_vid(const int32_t vid, const sellcs& g) {
    if(vid < 0 || vid > g.nverts) return -1;
    return g.permuts[vid];
} 

void delete_sellcs(sellcs& sellcs_g) {
    delete[] sellcs_g.cols;
    delete[] sellcs_g.cs;
    delete[] sellcs_g.cl;
    delete[] sellcs_g.permuts;
}