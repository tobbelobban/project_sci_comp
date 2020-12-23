#ifndef SELL_C_SIGMA_
#define SELL_C_SIGMA_

#include <cstdint>
#include <queue>
#include <string>
#include <vector>

typedef struct sellcs {
    int32_t C;
    int32_t nverts;
    int32_t n_chunks;
    int32_t* cs;
    int32_t* cl;
    int32_t* cols;
    int32_t* permuts;
    double beta;
    int32_t sigma;

} sellcs;

typedef struct vertex {
    int32_t degree;
    int32_t vid;
} vertex;

// mergesort functions
void sorter_serial(std::vector<vertex>& tmp, int32_t start, int32_t end, std::vector<vertex>& v);
void sorter(std::vector<vertex>& tmp, int32_t start, int32_t end, std::vector<vertex>& v);
void merger(std::vector<vertex>& v_tmp, int32_t start, int32_t mid, int32_t end, std::vector<vertex>& v_sorted);
void local_sort(std::vector<vertex>&, const sellcs&);

// bfs 
int32_t get_permutated_vid(const int32_t, const sellcs&);
void tropical_sellcs_mv_mult_w8(std::vector<int32_t>&, const sellcs&, const std::vector<int32_t>&);
void tropical_sellcs_mv_mult_w4(std::vector<int32_t>&, const sellcs&, const std::vector<int32_t>&);
int32_t sellcs_bfs(const sellcs&, const int32_t, std::vector<int32_t>&);
void permutate_solution(std::vector<int32_t>&, const sellcs&);

// misc
    double cpuSecond();
    void read_sellcs_graph_from_file(const std::string&, sellcs&, int, int, double* const);
    void print_sellcs_graph(const sellcs&);
    void delete_sellcs(sellcs&);

#endif