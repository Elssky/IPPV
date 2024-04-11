#include <iostream>
#include "Args.h"
#include "Graph.h"
#include <omp.h>


int main(int argc, char* argv[]) {
    EdgeList* el;
    setbuf(stdout, NULL);
    clock_t begin = clock();
    Args* args = new Args();
    args->parse_args(argc, argv);
    FILE* d_file = fopen(args->address, "r");
    omp_set_num_threads(4);

    printf("Reading edgelist from file %s\n", args->address);

    Graph g = Graph(d_file, args->NT, args->topk, args->h, args->p);

    clock_t io_end = clock();

    g.findLhCDS();
    g.output(args->ds_address);

    clock_t end = clock();
    double io_secs = double(io_end - begin) / CLOCKS_PER_SEC;
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    printf("io time: %.4f, total time: %.4f\n", io_secs, elapsed_secs);
    return 0;
}