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
<p align="center">
  <img src="https://github.com/user-attachments/assets/857f3a91-5bfa-446c-add3-d228c8eda6fd" width="80%">
<!--   <img src="https://github.com/user-attachments/assets/47d16988-84e6-4c7e-bffa-2872c702d3c1" width="45%"> -->
</p>


Left figure shows the relationships between a subset of characters in “Harry Potter”. The top-1 and top-2 L3CDSes are the blue and green subgraphs, respectively. The top-1 L3CDS is a family named Weasley, and the top-2 L3CDS is an organization named Death Eaters, which indicate the potential of LℎCDS discovery for mining diverse dense communities. 

##  Data Preparation
The input data format is: 
Fisrt row is `nodes_num edges_num`, then each following row `src_node dst_node` represents an edge. Nodes should start from 0, e.g.
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

### Table: Datasets used in our experiments

| **Name**             | **Abbr.** | **\|V\|** | **\|E\|** | 3-clique nums | 5-cliques num|
|----------------------|-----------|-----------|-----------|------------------|------------------|
| soc-hamsterster       | HA        | 2,426     | 16,630    | 53,251           | 298,013          |
| CA-GrQc              | GQ        | 5,242     | 14,484    | 48,260           | 2,215,500        |
| fb-pages-politician   | PP        | 5,908     | 41,706    | 174,632          | 2,002,250        |
| fb-pages-company      | PC        | 14,113    | 52,126    | 56,005           | 207,829          |
| web-webbase-2001      | WB        | 16,062    | 25,593    | 21,115           | 382,674          |
| CA-CondMat            | CM        | 23,133    | 93,439    | 173,361          | 511,088          |
| soc-epinions          | EP        | 26,588    | 100,120   | 159,700          | 521,106          |
| Email-Enron           | EN        | 36,692    | 183,831   | 727,044          | 5,809,356        |
| loc-gowalla           | GW        | 196,591   | 950,327   | 2,273,138        | 14,570,875       |
| DBLP                 | DB        | 317,080   | 1,049,866 | 2,224,385        | 262,663,639      |
| Amazon               | AM        | 334,863   | 925,872   | 667,129          | 61,551           |
| soc-youtube           | YT        | 495,957   | 1,936,748 | 2,443,886        | 5,306,643        |
| soc-lastfm            | LF        | 1,191,805 | 4,519,330 | 3,946,207        | 10,404,656       |
| soc-flixster          | FX        | 2,523,386 | 7,918,801 | 7,897,122        | 96,315,278       |
| soc-wiki-talk         | WT        | 2,394,385 | 4,659,565 | 9,203,519        | 382,777,822      |

All the datasets can be downloaded from [network repository](https://networkrepository.com/)

## Installation
```
mkdir build
cd build/
cmake ..
make
```
## Pipeline
<img width="978" alt="image" src="https://github.com/user-attachments/assets/2e03e968-be10-4977-87c3-1e78774ba752">

 Examples:  

`./LhCDScvx --args -g ./dataset/CA-CondMat.txt -h 3 -t 10 -k 20 -p 0` 

`h` means h-clique, `t` is iteration, `k` means top-k LhCDSm `p` is pattern.


The following are the regular subgraph models we support.
<img width="800" alt="image" src="https://github.com/user-attachments/assets/e6567574-18e6-4b1e-b941-40988a3dc615">


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
### 1. Efficiency under Parameter Variations
![image](https://github.com/user-attachments/assets/bd66566b-08e4-4b8e-aa09-3c8021076222)
Execute the following shell file to obtain the results of this experiment (note: the specific time may vary due to differences in hardware devices)
```shell
k_list=(5 10 15 20)
h_list=(3 4 5)
v_list=(1 0)
t=20
p=0
data_list=("fb-pages-company" "soc-hamsterster" "soc-epinions" "Email-Enron" "loc-gowalla" "CA-CondMat" "CA-GrQc" "Amazon")

for dataset in ${data_list[@]}
do
    if [ -z "$dataset" ]; then
    echo "Error: The variable 'dataset' is not set."
    exit 1
    fi

    # Create the directory if it does not exist
    mkdir -p "./output/$dataset"

    for k in ${k_list[@]}
    do
        for h in ${h_list[@]}
        do
            for t in ${t_list[@]}
            do
                for s in ${s_list[@]}
                do
                    for v in ${v_list[@]}
                    do
                    (/usr/bin/time -v nohup ./LDScvx -g ./dataset/${dataset}.txt -h $h -k $k -t $t -p $p -v $v&) >&  ./output/$dataset/${dataset}_LhCDScvx_k=${k}_h=${h}_t=${t}_s=${s}_p=${p}_v=${v}.txt
                    done
                done
            done
        done
    done
done

```

### 2. Efficiency v.s. Existing Algorithms
![image](https://github.com/user-attachments/assets/953d3a2e-921c-4539-91be-eb52b396dcc6)
The code of LDSflow comes from the author of Locally Densest Subgraph Discovery.
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
