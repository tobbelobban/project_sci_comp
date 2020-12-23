#include <iostream>
#include <getopt.h>
#include <cmath>
#include <omp.h>
#include "csr.hpp"
#include "sellcs.hpp"
#include <sys/time.h>


void calc_statistics(double* stats, const double* times, int size) {
    
	double temp=0.0;
	int i;
	// compute the mean
	
	for (i = 0; i < size; ++i) temp += times[i];
    stats[0] = temp;
	temp /= size;
	stats[1] = temp;
	double mean = temp;
    temp = 0;

    // compte the standard deviation
	for (i = 0; i < size; ++i) temp += (times[i] - mean) * (times[i] - mean);
	temp /= size - 1;
	stats[2] = sqrt(temp);
}

int is_same_solution(const std::vector<int32_t>& res1, const std::vector<int32_t>& res2, const int size) {
    for(int32_t i = 0; i < size; ++i)
    {
        if(res1[i] != res2[i]) return 0;
    }
    return 1;
}

double cpuSecond() {
    struct timeval tp;
    gettimeofday(&tp,NULL);
    return ((double)tp.tv_sec + (double)tp.tv_usec*1.e-6);
}

int main(int argc, char *const *argv) {

    int32_t root=0, chunk_size=8, sigma=1,SCALE=-1,EDGEFACTOR=-1, mode=1;
    int32_t* roots;
    std::string file_path = "none", roots_file = "none";
    if(argc < 2) {
        std::cout << "Error. Need filename with graph." << std::endl;
        exit(0);
    }

    // parse arguments
    int opt;
    while ((opt = getopt(argc,argv,":C:s:S:E:r:f:m:R:")) != -1)
    {
        switch (opt)
        {
            case 'C':
                chunk_size = std::atoi(optarg);
                break;
            case 's':
                sigma = std::atoi(optarg);            
                break;
            case 'S':
                SCALE = std::atoi(optarg);
                break;
            case 'E':
                EDGEFACTOR = std::atoi(optarg);
                break;
            case 'f':
                file_path = optarg;
                break;
            case 'r':
                root = std::atoi(optarg);
                break;
            case 'm':
                mode = std::atoi(optarg);
                break;
            case 'R':
                roots_file = optarg;
                break;
            default:
                std::cout << "Unrecognized flag: " << opt << "\nTerminating!" << std::endl;
                exit(0);
                break;
        }
    }

    // check arguments
    if(SCALE == -1) {
        std::cout << "ERROR! SCALE must be supplied with: -S X, where X is 2-logarithm of number of vertices." << std::endl;
        exit(0);
    }
    if(EDGEFACTOR == -1) {
        std::cout << "ERROR! EDGEFACTOR must be supplied with: -E X, where X is #edges / #vertices." << std::endl;
        exit(0);
    }
    if(root < 0 || root >= pow(2,SCALE)) {
        std::cout << "ERROR! r = " << root << " is not valid root. Root must be 0 <= root < 2^SCALE. Supply root with: -r X" << std::endl;
        exit(0);
    }
    if(chunk_size != 4 && chunk_size != 8) {
        std::cout << "ERROR! chunk size must be =4 or =8. Modify chunk size with: -C X" << std::endl;
        exit(0);
    }
    if(sigma < 0 || sigma > SCALE) {
        std::cout << "ERROR! sigma must be 0 <= sigma <= SCALE. Modify sigma with: -s X" << std::endl;
        exit(0);
    }
    if(file_path == "none") {
        std::cout << "ERROR! file path containing edges must be specified! Use: -f path/to/file/" << std::endl; 
        exit(0);
    }
    if(mode < 1 || mode > 3) {
        std::cout << "ERROR! Mode must be one of the following:\n" << std::endl;
        std::cout << "mode = 1: CSR only.\n" << std::endl;
        std::cout << "mode = 2: SELL-C-Sigma only.\n" << std::endl;
        std::cout << "mode = 3: CSR and SELL-C-Sigma.\n" << std::endl;
        std::cout << "supply mode with flag m: -m <mode> \n" << std::endl;
        exit(0);
    }
    if(roots_file != "none") {
        FILE* roots_ptr = fopen(roots_file.c_str(), "rb");
        if(roots_ptr == NULL)
        {
            std::cout << "ERROR! Failed to open file: " << roots_file << std::endl;
            exit(0);
        }
        roots = new int32_t[root];
        if(fread(roots, sizeof(int32_t), root, roots_ptr) != (size_t)root) 
        {
            std::cout << "Error reading roots from file: " << roots_file << " ... exiting!" << std::endl;
            exit(0);
        }
    } else {
        roots = new int32_t[1];
        roots[0] = root;
        root = 1;
    }
    
    // timer indices have following defintions:
        // timer[0] = time to read graph to buffer
        // timer[1] = time to process edges in buffer
        // timer[2] = time to sort edges (does not exist for CSR)
        // timer[3] = time to create SELL-C-SIGMA datastructure

    double* CSR_create_timer, *SELLCS_create_timer, *CSR_compute_timer, *SELLCS_compute_timer, *SELLCS_permutate_timer;
    csr_graph csr_g;
    sellcs sellcs_g;
    int32_t nverts = pow(2,SCALE);

    if(mode & 1) 
    {    
        // timing for reading graph
        std::cout << "-- CSR --\n" << std::endl;
        std::cout << "Creating CSR from file, please wait...";
        std::flush(std::cout);

        CSR_create_timer = new double[3];
        CSR_compute_timer = new double[root];
        read_csr_graph_from_file(file_path, csr_g, SCALE, EDGEFACTOR, CSR_create_timer);
        //print_csr_graph(csr_g);

        std::cout << " done!" << std::endl;
        std::cout <<  "Buffer read time: \t\t" << CSR_create_timer[0] << " s"  << std::endl;
        std::cout << "Process buffer time: \t\t" << CSR_create_timer[1] << " s"  << std::endl;
        std::cout << "Create CSR time: \t\t" << CSR_create_timer[2] << " s"  << std::endl;
        std::cout << std::endl;
    }

    // BFS with SELL-C-SIGMA datastructure
    if(mode & 2) 
    {
        sellcs_g.C = chunk_size;
        sellcs_g.sigma = pow(2,sigma);

        // timing for reading graph
        std::cout << "-- SELL-C-SIGMA --\n" << std::endl;
        std::cout << "Creating SELL-C-SIGMA from file, please wait...";
        std::flush(std::cout);

        SELLCS_create_timer = new double[5];
        SELLCS_compute_timer = new double[root];
        SELLCS_permutate_timer = new double[root];
        read_sellcs_graph_from_file(file_path, sellcs_g, SCALE, EDGEFACTOR, SELLCS_create_timer);
        //print_sellcs_graph(sellcs_g);

        std::cout << " done!" << std::endl;
        std::cout <<  "Buffer read time: \t\t" << SELLCS_create_timer[0] << " s"  << std::endl;
        std::cout << "Process buffer time: \t\t" << SELLCS_create_timer[1] << " s"  << std::endl;
        std::cout << "Sort vertices time: \t\t" << SELLCS_create_timer[2] << " s"  << std::endl;
        std::cout << "Create SELL-C-SIGMA time: \t" << SELLCS_create_timer[3] << " s"  << std::endl;
        std::cout << "Beta = \t\t\t\t" << sellcs_g.beta << std::endl;
        std::cout << std::endl;
    }

    std::cout << "Computing BFS...";
    std::flush(std::cout);

    std::vector<int32_t> res_CSR, res_SELLCS;
    int32_t CSR_its=0, SELLCS_its=0, tmp, same_count=0;

    for(int32_t i = 0; i < root; ++i)
    {
        int32_t curr_root = roots[i];
        if(curr_root < 0 || curr_root > nverts) continue;

        if(mode & 1)
        {
            // warm up
            res_CSR = std::vector<int32_t>(csr_g.nverts,csr_g.nverts);
            res_CSR[curr_root] = 0;
            tmp = csr_bfs(csr_g, curr_root, res_CSR);

            // now compute with timing
            res_CSR = std::vector<int32_t>(csr_g.nverts,csr_g.nverts);
            res_CSR[curr_root] = 0;
            CSR_compute_timer[i] = cpuSecond();
            tmp = csr_bfs(csr_g, curr_root, res_CSR);
            CSR_compute_timer[i] = cpuSecond() - CSR_compute_timer[i];
            if(tmp == -1)
            {
                CSR_compute_timer[i] = 0;
                std::cout << "Failed to compute BFS with CSR, root = " << curr_root << ", i = " << i << std::endl;
            } 
            else
            {
                CSR_its += tmp;
            }
        }
        if(mode & 2)
        {            
            // warm up
            int32_t permutated_root = get_permutated_vid(curr_root, sellcs_g);
            res_SELLCS = std::vector<int32_t>(nverts,nverts);
            res_SELLCS[permutated_root] = 0;
            tmp = sellcs_bfs(sellcs_g, permutated_root, res_SELLCS);

            // compute with timing
            SELLCS_compute_timer[i] = cpuSecond();
            tmp = sellcs_bfs(sellcs_g, permutated_root, res_SELLCS);
            SELLCS_compute_timer[i] = cpuSecond() - SELLCS_compute_timer[i];
            
            SELLCS_permutate_timer[i] = cpuSecond();
            permutate_solution(res_SELLCS, sellcs_g);
            SELLCS_permutate_timer[i] = cpuSecond() - SELLCS_permutate_timer[i];
            if(tmp == -1)
            {
                SELLCS_compute_timer[i] = 0;
                std::cout << "Failed to compute BFS with SELLCS, root = " << curr_root << ", i = " << i << std::endl;
            } 
            else
            {
                SELLCS_its += tmp;
            }
        }
        if((mode & 1) && (mode & 2))
        {
            same_count += is_same_solution(res_CSR, res_SELLCS, nverts);    
        }
    }
    std::cout << " DONE!\n" << std::endl;

    // print results
    double stats[3];
    if(mode & 1)
    {
        calc_statistics(stats, CSR_compute_timer, root);
        delete_csr(csr_g);

        std::cout << "-- CSR RESULTS --" << std::endl;
        std::cout << "Total number matrix-vector mult.\t= " << CSR_its << std::endl;
        std::cout << "Total BFS-CSR time\t\t\t= " << stats[0] << " s" << std::endl;
        std::cout << "Mean time per matrix-vector mult.\t= " << stats[0]/(double)CSR_its << " s" << std::endl;
        std::cout << "Mean BFS-CSR time\t\t\t= " << stats[1] << " s" << std::endl;
        std::cout << "Std. deviation BFS-CSR\t\t\t= " << stats[2] << " s" << std::endl;
        std::cout << std::endl;
    }
    if(mode & 2)
    {
        calc_statistics(stats, SELLCS_compute_timer, root);
        delete_sellcs(sellcs_g);

        std::cout << "-- SELL-C-SIGMA RESULTS --" << std::endl;
        std::cout << "Total number matrix-vector mult. \t= " << SELLCS_its << std::endl;
        std::cout << "Total BFS-SELLCS time \t\t\t= " << stats[0] << " s" << std::endl;
        std::cout << "Mean time per matrix-vector mult. \t= " << stats[0]/(double)SELLCS_its << " s" << std::endl;
        std::cout << "Mean BFS-SELLCS time \t\t\t= " << stats[1] << " s" << std::endl;
        std::cout << "Std. deviation BFS-SELLCS \t\t= " << stats[2] << " s" << std::endl;
        std::cout << std::endl;
    }
    if(mode & 2 && mode & 1)
    {
        std::cout << "-- COMPARISON RESULTS --" << std::endl;
        std::cout << "BFS compute for " << root << (root > 1 ? " roots." : " root.") << std::endl;
        std::cout << "Number of same results between BFS-CSR and BFS-SELLCS:\t" << same_count << std::endl;
        std::cout << std::endl; 
    }

    delete[] roots;
    return 0;
}
