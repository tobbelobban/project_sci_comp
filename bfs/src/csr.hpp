#ifndef _CSR_GRAPH
#define _CSR_GRAPH
#include <cstdint>
#include <vector>

typedef struct csr_graph {
    uint32_t* cols;
    uint32_t* rows;
    uint32_t nverts;
    uint32_t nedges;
} csr_graph;

void delete_csr(csr_graph&);
void read_csr_graph_from_file(const std::string&, csr_graph&);
void print_csr_graph(const csr_graph& csr_g);
void tropical_csr_mv_mult(std::vector<uint32_t>&, const csr_graph&, const std::vector<uint32_t>&);
void tropical_iin_csr_mv_mult(std::vector<uint32_t>&, const csr_graph&, const std::vector<uint32_t>&);
std::vector<uint32_t> csr_bfs(const csr_graph&, const uint32_t);
std::vector<uint32_t> csr_bfs_iin(const csr_graph&, const uint32_t);
bool same(const std::vector<uint32_t>&, const std::vector<uint32_t>&);
void print_vector(const std::vector<uint32_t>& v);

#endif