#include <iostream>
#include <getopt.h>
#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <string>

#include "make_graph.h"
#include "../src/csr.hpp"
#include "graph.hpp"
#include <sys/time.h>

double cpuSecond() {
    struct timeval tp;
    gettimeofday(&tp,NULL);
    return ((double)tp.tv_sec + (double)tp.tv_usec*1.e-6);
}

int main(int argc, char* const* argv) {
    uint64_t scale = 10, edge_factor = 16, s1 = 2, s2 = 32;
    int32_t num_roots = -1;
    std::string f_path;
    int option, tmp;
    bool in_fn = false;
    while((option = getopt(argc, argv, "s:e:f:t:y:r:")) != -1) {
        switch (option)
        {
        case 's':
            tmp = atoi(optarg);
            scale = tmp > 0 ? tmp : scale;
            break;
        case 'e':
            tmp = atoi(optarg);
            edge_factor = tmp > 0 ? tmp : edge_factor;
            break;
        case 'f':
            f_path = optarg;
            in_fn = true;
            break;
        case 't':
            tmp = atoi(optarg);
            s1 = tmp > 0 ? tmp : s1;
            break;
        case 'y':
            tmp = atoi(optarg);
            s2 = tmp > 0 ? tmp : s1;
            break;
        case 'r':
            tmp = atoi(optarg);
            num_roots = tmp > 0 ? tmp : num_roots;
            break;
        default:
            std::cout << "Unrecognized argument: " << option << std::endl;
            break;
        }
    }
    if(!in_fn)
        f_path = "../bin/" + std::to_string(scale) + '_' + std::to_string(edge_factor) + ".bin"; 
    std::cout << "GENERATING GRAPH TO FILE: " << f_path << std::endl;
    generate_graph_to_file(f_path, scale, edge_factor,s1,s2);
    std::cout << "--GENERATION DONE!--" << std::endl;
    if(num_roots > 0)
    {
        std::cout << "GENERATING ROOTS" << std::endl;
        std::flush(std::cout);
        int32_t* bfs_roots = (int32_t*)malloc(num_roots * sizeof(int32_t));
        csr_graph csr;
        double CSR_timer[5];
        read_csr_graph_from_file(f_path, csr, scale, edge_factor, CSR_timer);

        // this root generator is based off of generator used by the supplied Graph500 reference bfs
        // I take no credit for this solution
        int64_t counter = 0;
		int bfs_root_idx;
		for (bfs_root_idx = 0; bfs_root_idx < num_roots; ++bfs_root_idx) {
			int32_t root;
			while (1) {
				double d[2];
				make_random_numbers(2, s1, s2, counter, d);
				root = (int32_t)((d[0] + d[1]) * csr.nverts) % csr.nverts;
				counter += 2;
				if (counter > 2 * csr.nverts) break;
				int is_duplicate = 0;
				int i;
				for (i = 0; i < bfs_root_idx; ++i) {
					if (root == bfs_roots[i]) {
						is_duplicate = 1;
						break;
					}
				}
				if (is_duplicate) continue; /* Everyone takes the same path here */
				int root_bad = isisolated(root,csr);
				if (!root_bad) break;
			}
			bfs_roots[bfs_root_idx] = root;
		}
        delete_csr(csr);
		num_roots = bfs_root_idx;
        
        char buffer[40];
        if(sprintf(buffer, "../bin/s%li_ef%li_roots%i.bin", scale, edge_factor, num_roots) > 0)
        {
            FILE* file_ptr = fopen(buffer, "wb");
            if(fwrite(bfs_roots, sizeof(int32_t),num_roots, file_ptr) != (size_t)num_roots)
            {
                printf("Error writing roots to file %s\n", buffer);
            }
            fclose(file_ptr);
            free(bfs_roots);
        }   
        std::cout << "--GENERATION DONE!--" << std::endl;
        std::cout << "Roots in file: " << buffer << std::endl;
    }
    return 0;
}
