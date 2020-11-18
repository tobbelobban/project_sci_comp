#ifndef _CSR_GRAPH
#define _CSR_GRAPH
#include <cstdint>
#include <vector>

typedef struct csr_graph {
    int32_t* cols;
    int32_t* rows;
    int32_t nverts;
    int32_t nedges;
} csr_graph;

void delete_csr(csr_graph&);
void read_csr_graph_from_file(const std::string&, csr_graph&);
void print_csr_graph(const csr_graph& csr_g);
void tropical_csr_mv_mult(std::vector<int32_t>&, const csr_graph&, const std::vector<int32_t>&);
void tropical_iin_csr_mv_mult(std::vector<int32_t>&, const csr_graph&, const std::vector<int32_t>&);
std::vector<int32_t> csr_bfs(const csr_graph&, const int32_t);
std::vector<int32_t> csr_bfs_iin(const csr_graph&, const int32_t);
bool same(const std::vector<int32_t>&, const std::vector<int32_t>&);
void print_vector(const std::vector<int32_t>& v);

#endif