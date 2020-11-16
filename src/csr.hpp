#ifndef _CSR_GRAPH
#define _CSR_GRAPH
#include <cstdint>

typedef struct csr_graph {
    int64_t* cols;
    int64_t* rows;
    int64_t nverts;
    int64_t nedges;
} csr_graph;

void read_csr_graph_from_file(const char*, csr_graph&);
void print_csr_graph(const csr_graph& csr_g);
#endif