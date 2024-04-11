#ifndef LhCDSCVX_GRAPH_H
#define LhCDSCVX_GRAPH_H

#include <vector>
#include <stack>
#include <algorithm>
#include <queue>
#include <unordered_map>
#include "KClistCore.h"

#include "FlowNetwork.h"

const unsigned MAX_EDGES = 100000000; // Maximum number of edges for memory allocation; will increase if needed

inline unsigned max3(unsigned a, unsigned b, unsigned c) {
    a = (a > b) ? a : b;
    return (a > c) ? a : c;
}

using namespace std;

typedef struct {
    unsigned s;
    unsigned t;
} Edge;


typedef struct {
    unsigned n; // Number of nodes
    unsigned e; // Number of edges
    Edge* edges; // List of edges
    unsigned* rank; // Ranking of the nodes according to degeneracy ordering
} EdgeList;


EdgeList* ReadEdgeList(char* file_name);

class Graph {
public:
    Graph(FILE* pFile, int NT, int topk, int h, int p);
    void findLhCDS();

    void output(char* ds_address);

private:
    int n;
    unsigned long m;
    int h; // h-clique
    int p; // pattern 
    int hc_num; // number of h-cliques
    int NT;
    int CT;
    int topk;
    int nsg;
    int max_d;
    int num_verify;
    bool check_first = false;
    bool fn_baseline = false;
    //Edge* edges; // List of edges

    vector<vector<int>> adj;
    vector<bool> selected;
    vector<bool> active;
    vector<int> nodes;
    vector<int> h_cliques_no;
    vector<int> slt_nodes;
    vector<int> slt_edges;
    vector<int> slt_h_cliques;
    stack<vector<int>> stk_nodes;
    stack<vector<int>> stk_edges;
    stack<vector<int>> stk_h_cliques;
    stack<int> stk_CT;
    vector<vector<int>> lhcdses;
    vector<vector<int>> lhcds_edge;
    vector<double> lhcds_rho;
    vector<double>  lhcds_den;
    vector<int> lhcds_num;
    vector<int> deg;
    vector<int> core;
    vector<int> pos;
    vector<int> sg;
    vector<pair<int, int>> edges;
    vector<vector<int>> h_cliques;
    vector<vector<int>> v2cliques;
    vector<double> alpha;
    vector<std::unordered_map<int, double>> alpha_c;
    vector<double> r;
    vector<double> rho_u;
    vector<double> rho_gu;
    vector<double> rho_l;
    vector<double> val;
    vector<int> nag;
    vector<int> fa;
    //vector<pair<int, int>> cmpt;
    vector<std::tuple<int, int, int>> cmpt;
    vector<int> veri_vtx;
    //void clique_enum();
    void frank_wolfe();
    void pava();
    void check_sg();
    void pruning();
    void compute_core();
    void prune_by_core();
    bool verify_by_maxflow(vector<int>& nodes, double g);
    bool verify_LhCDS(vector<int>& nodes, double g);
    bool verify_LhCDS(vector<int>temp_nodes, vector<int>& nodes, double g);
    bool verify_LhCDS_small(vector<int>& nodes, double g, vector<int>& tmp_hcliques);
    int find_fa(int x);
    void connected_components();
    // compute h-clique-density
    void compute_h_clique_core();
    void frank_wolfe_h_clique();
    void frank_wolfe_h_clique_new();
    void compute_h_cliques();
    void enumerateKCliques(int k);
    void backtrack(int currentNode, int k, vector<int>& currentClique);
    bool isValidAddition(const vector<int>& currentClique, int newNode);
    void pava_h_clique();
    void check_sg_h_clique();
    bool Get_Maxflow();
    void LODA_hclique();
    void new_LODA();
    void LODA_sort();
};





double ComputeAverageClusteringCoefficient(Graph g, vector<vector<int>> groups);

#endif //LhCDSCVX_GRAPH_H