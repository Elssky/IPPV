# IPPV (Iterative Propose-Prune-and-Verify)
[SIGMOD2025] An Efficient and Exact Algorithm for Locally h-Clique Densest Subgraph Discovery

Detecting locally non-overlapping, near-clique densest subgraphs is essential for community search in social networks. Our work introduces an efficient and exact algorithm to tackle this challenge by identifying the top-k non-overlapping, locally h-clique densest subgraphs (LhCDS).

> 检测局部不重叠、接近派系最密集的子图对于社交网络中的社区搜索至关重要。我们的工作引入了一种高效且精确的算法，通过识别 top-k 非重叠、局部 h 团最稠密子图 (LhCDS) 来应对这一挑战。

Key highlights of our approach:

- Iterative Propose-Prune-and-Verify (IPPV) Pipeline: 
  - [x] Propose: Generate candidates using convex programming without missing any dense substructures.
  - [x] Prune: Handle overlapping cliques with an advanced graph decomposition technique.
  - [x] Verify: Use basic and fast verification methods to confirm LhCDS candidates, leveraging a flow network to enhance speed.
- Generalization: Our approach extends to broader densest subgraph detection problems while ensuring exactness and low computational complexity.
  - [x] Extensive experiments on real datasets demonstrate the effectiveness and efficiency of our IPPV algorithm.

[ArXiv](https://www.arxiv.org/abs/2408.14022)
![image](https://github.com/user-attachments/assets/84075e89-39bb-4122-9172-feacd6308a82)
Figure 1 shows the relationships between a subset of characters in “Harry Potter”. The top-1 and top-2 L3CDSes are the blue and green subgraphs, respectively. The top-1 L3CDS is a family named Weasley, and the top-2 L3CDS is an organization named Death Eaters, which indicate the potential of LℎCDS discovery for mining diverse dense communities.
## Code
`mkdir output` 

 Examples:  

`./LhCDScvx --args -g ./dataset/CA-CondMat.txt -h 3 -t 10 -k 20 -p 0` 

`h` means h-clique, `t` is iteration, `k` means top-k LhCDSm `p` is pattern.

The input data format is: 
Fisrt row is `nodes_num, edges_num`, then each following row represents an edge. Node is should start from 0, e.g.
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
## Datasets

## Experiments
### 1. Efficiency under Parameter Variations
![image](https://github.com/user-attachments/assets/bd66566b-08e4-4b8e-aa09-3c8021076222)
### 2. Efficiency v.s. Existing Algorithms
![image](https://github.com/user-attachments/assets/953d3a2e-921c-4539-91be-eb52b396dcc6)
### 3. Memory Overheads
![image](https://github.com/user-attachments/assets/ca39b3cf-f031-4e0f-962f-4232653fc9b2)



## Citation

```
@article{xu2024efficient,
  author = {Xiaojia Xu and Haoyu Liu and Xiaowei Lv and Yongcai Wang and Deying Li},
  title = {An Efficient and Exact Algorithm for Locally h-Clique Densest Subgraph Discovery},
  journal = {Proc. ACM Manag. Data},
  volume = {2},
  number = {N6},
  article = {225},
  year = {2024},
  month = {December},
  pages = {26},
  doi = {10.1145/3698800}
}
```
