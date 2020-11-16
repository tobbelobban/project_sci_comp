#ifndef _CSR_GRAPH
#define _CSR_GRAPH
#include <cstdint>
#include <vector>

typedef struct csr_graph {
    int64_t* cols;
    int64_t* rows;
    int64_t nverts;
    int64_t nedges;
} csr_graph;

void delete_csr(csr_graph&);
void read_csr_graph_from_file(const char*, csr_graph&);
void print_csr_graph(const csr_graph& csr_g);
void tropical_csr_mv_mult(std::vector<int64_t>&, const csr_graph&, const std::vector<int64_t>&);
void tropical_iin_csr_mv_mult(std::vector<int64_t>&, const csr_graph&, const std::vector<int64_t>&);
std::vector<int64_t> csr_bfs(const csr_graph&, const int64_t);
std::vector<int64_t> csr_bfs_iin(const csr_graph&, const int64_t);
bool same(const std::vector<int64_t>&, const std::vector<int64_t>&);
void print_vector(const std::vector<int64_t>& v);

#endif