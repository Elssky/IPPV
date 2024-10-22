# "An Efficient and Exact Algorithm for Locally ℎ-Clique Densest Subgraph Discovery [SIGMOD 2025]"
Detecting locally, non-overlapping, near-clique densest subgraphs is a crucial problem for community search in social networks. As a vertex may be involved in multiple overlapped local cliques, detecting locally densest sub-structures considering h-clique density, i.e., locally h-clique densest subgraph (LhCDS) attracts great interests. This paper investigates the LhCDS detection problem and proposes an efficient and exact algorithm to list the top-k non-overlapping, locally h-clique dense, and compact subgraphs. We in particular jointly consider h-clique compact number and LhCDS and design a new "Iterative Propose-Prune-and-Verify" pipeline (IPPV) for top-k LhCDS detection. (1) In the proposal part, we derive initial bounds for h-clique compact numbers; prove the validity, and extend a convex programming method to tighten the bounds for proposing LhCDS candidates without missing any. (2) Then a tentative graph decomposition method is proposed to solve the challenging case where a clique spans multiple subgraphs in graph decomposition. (3) To deal with the verification difficulty, both a basic and a fast verification method are proposed, where the fast method constructs a smaller-scale flow network to improve efficiency while preserving the verification correctness. The verified LhCDSes are returned, while the candidates that remained unsure reenter the IPPV pipeline. (4) We further extend the proposed methods to locally more general pattern densest subgraph detection problems. We prove the exactness and low complexity of the proposed algorithm. Extensive experiments on real datasets show the effectiveness and high efficiency of IPPV.
## Introduction
![image](https://github.com/user-attachments/assets/84075e89-39bb-4122-9172-feacd6308a82)
Figure 1 shows the relationships between a subset of characters in “Harry Potter”. The top-1 and top-2 L3CDSes are the blue and green subgraphs, respectively. The top-1 L3CDS is a family named Weasley, and the top-2 L3CDS is an organization named Death Eaters, which indicate the potential of LℎCDS discovery for mining diverse dense communities.
## Code
`mkdir output` 

 Examples:  

`./LhCDScvx --args -g ./dataset/CA-CondMat.txt -h 3 -t 10 -k 20 -p 0` 

`h` means h-clique, `t` is iteration, `k` means top-k LhCDSm `p` is pattern.

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
## Experiments
### Efficiency under Parameter Variations
![image](https://github.com/user-attachments/assets/bd66566b-08e4-4b8e-aa09-3c8021076222)
### Efficiency v.s. Existing Algorithms
![image](https://github.com/user-attachments/assets/953d3a2e-921c-4539-91be-eb52b396dcc6)
### Memory Overheads
![image](https://github.com/user-attachments/assets/ca39b3cf-f031-4e0f-962f-4232653fc9b2)



## Citation

`
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

`
