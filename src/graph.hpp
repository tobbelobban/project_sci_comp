#ifndef _READ_GRAPH
#define _READ_GRAPH

/*
    GRAPH READER/WRITER
    EDGE SORTER
    EDGE PRINTER
    PACKED_EDGE: two int64_t 
*/

typedef struct graph {
    int64_t nverts;
    int64_t nedges;
    packed_edge* edges;
} graph;

void read_graph_from_file(const char*, graph&);
bool compare_edges(const packed_edge&, const packed_edge&);
void sort_graph_edges(graph&);
void delete_graph(graph&);
void print_graph_edges(const graph&);
void write_graph_to_file(const char*, const graph&);

#endif