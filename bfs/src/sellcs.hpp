#ifndef SELL_C_SIGMA_
#define SELL_C_SIGMA_

#include <cstdint>
#include <string>
#include <vector>

typedef struct sellcs {
    int32_t C;
    int32_t nverts;
    int32_t n_chunks;
    int32_t* cs;
    int32_t* cl;
    int32_t* cols;

} sellcs;

void read_sellcs_graph_from_file(const std::string&, sellcs&);
void print_sellcs_graph(const sellcs&);
void tropical_sellcs_mv_mult(std::vector<int32_t>&, const sellcs&, const std::vector<int32_t>&);
void delete_sellcs(sellcs&);

#endif