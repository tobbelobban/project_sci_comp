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

typedef struct vq {
    std::queue<int32_t> e;
    int32_t vid;
} vq;

void read_sellcs_graph_from_file(const std::string&, sellcs&);
void print_sellcs_graph(const sellcs&);
void tropical_sellcs_mv_mult(std::vector<int32_t>&, const sellcs&, const std::vector<int32_t>&);
std::vector<int32_t> sellcs_bfs(const sellcs&, const int32_t);
int32_t get_permutated_vid(const int32_t, sellcs&);
void permutate_solution(std::vector<int32_t>&, const sellcs&);
void delete_sellcs(sellcs&);

#endif