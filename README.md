# "An Efficient and Exact Algorithm for Locally ℎ-Clique Densest Subgraph Discovery"

## Usage
`mkdir output` 

 Examples:  

`./LhCDScvx --args -g ./dataset/CA-CondMat.txt -h 3 -t 10 -k 20 -p 0` 

The input data format is: Fisrt row is `nodes_num, edges_num`, then each following row represents an edge. Node is should start from 0, e.g.
```csv
23133 93439
0 1
0 2
0 3
0 4
0 5
0 6
0 7
0 8
0 9
```

The following are the regular subgraph models we support.

``` cpp
if (p == 0)
    {
        h_cliques = get_kcliques(adj, v2cliques, k, selected);
    }
    else if (p == 1)
    {
        h_cliques = get_kloops(adj, v2cliques, k, selected);
    }
    else if (p == 2)
    {
        h_cliques = get_doubletriangle(adj, v2cliques, k, selected);
    }
    else if (p == 3)
    {
        h_cliques = get_tritriangle(adj, v2cliques, k, selected);
    }
    else if (p == 4)
    {
        h_cliques = get_kstar(adj, v2cliques, k - 1, selected);
    }
    else if (p == 5)
    {
        h_cliques = get_ctriangle(adj, v2cliques, k - 1, selected);
    }
    for (int i = 0; i < h_cliques.size(); i++) {
        h_cliques_no.push_back(i);
    }
```
