#ifndef SELL_C_SIGMA_
#define SELL_C_SIGMA_

#include <cstdint>
#include <string>

typedef struct sellcs {
    uint32_t C;
    uint32_t nverts;
    uint32_t n_chunks;
    uint32_t* cs;
    uint32_t* cl;
    uint32_t* cols;

} sellcs;

void read_sellcs_graph_from_file(const std::string&, sellcs&);
void print_sellcs_graph(const sellcs&);
void delete_sellcs(sellcs&);

#endif