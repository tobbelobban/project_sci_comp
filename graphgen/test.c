#include <stdio.h>
#include <math.h>
#include <omp.h>

#include "make_graph.h"

int main(int argc, char const *argv[]) {
    double start, end;

    int log_num_verts = 20;
    int64_t desired_nedges = 16 * pow(2, log_num_verts);
    uint64_t seed1 = 91;
    uint64_t seed2 = 142;
    int64_t nedges;
    struct packed_edge* result;
    
    start = omp_get_wtime();
    make_graph(log_num_verts, desired_nedges, seed1, seed2, &nedges, &result);
    end = omp_get_wtime();

    double time_taken = end - start;

    int id;
    #pragma omp parallel private(id) 
    {
        id = omp_get_thread_num();
        printf("My id: %i\n", id);    
    }

    printf("%li edge%s generated in %fs (%f Medges/s)\n", nedges, (nedges == 1 ? "" : "s"), time_taken, 1. * nedges / time_taken * 1.e-6);

    
    free(result);
    return 0;
}
