#include "Graph.h"
#include <ctime>
#include <iostream>
#include <unordered_map>
#include "MaxFlow.h"
#include <omp.h>

using namespace std;

#define DELTA 1e-5

#define EPS 1e-9

bool IsEqual(double a, double b)
{
    return (fabs(a - b) < EPS);
}

class SmallGraph {
public:
    std::unordered_map<int, vector<int>> adj; 

    void addEdge(int v, int w) {
        adj[v].push_back(w);
        adj[w].push_back(v);
    }

    double clusteringCoefficient(int v) {
        if (adj[v].size() < 2) return 0.0; 

        double possibleEdges = 0.0, existingEdges = 0.0;
        for (int i = 0; i < adj[v].size(); ++i) {
            for (int j = i + 1; j < adj[v].size(); ++j) {
                possibleEdges++;
                if (find(adj[adj[v][i]].begin(), adj[adj[v][i]].end(), adj[v][j]) != adj[adj[v][i]].end()) {
                    existingEdges++;
                }
            }
        }
        return (2.0 * existingEdges) / (adj[v].size() * (adj[v].size() - 1));
    }

    double averageClusteringCoefficient() {
        double sum = 0.0;
        for (auto& p : adj) {
            sum += clusteringCoefficient(p.first);
        }
        return sum / adj.size();
    }
};


Graph::Graph(FILE* file, int NT, int topk, int h, int p) {
    printf("Graph construction\n");
    fscanf(file, "%d%lu", &n, &m);

    adj.resize(n);
    r.resize(n);
    deg.resize(n);
    core.resize(n);
    selected.resize(n, true);
    active.resize(n, true);
    lhcds_num.resize(n, -1);
    num_verify = 0;
    veri_vtx.resize(n, 0);
    rho_gu.resize(n);
    rho_u.resize(n);
    rho_l.resize(n);
    val.resize(n);
    nag.resize(n);
    sg.resize(n, -1);
    fa.resize(n);
    slt_nodes.clear();
    slt_edges.clear();
    slt_h_cliques.clear();

    //    printf("ok1\n");
    for (int i = 0; i < m; i++) {
        int u, v;
        fscanf(file, "%d%d", &u, &v);
        edges.emplace_back(u, v);
        adj[u].push_back(v);
        adj[v].push_back(u);
    }
    alpha.resize(m);
    alpha_c.resize(hc_num);
    //    printf("ok2\n");

    this->NT = NT;
    this->topk = topk;
    this->h = h;
    this->p = p;
}

EdgeList* ReadEdgeList(char* file_name) {
    unsigned e1 = MAX_EDGES;
    EdgeList* el = (EdgeList*)malloc(sizeof(EdgeList));
    FILE* file = fopen(file_name, "r");
    el->n = 0;
    el->e = 0;
    el->edges = (Edge*)malloc(e1 * sizeof(Edge));
    while (fscanf(file, "%u %u", &(el->edges[el->e].s), &(el->edges[el->e].t)) == 2) { // Read one edge
        el->n = max3(el->n, el->edges[el->e].s, el->edges[el->e].t);
        ++el->e;
        if (el->e == e1) {
            e1 += MAX_EDGES;
            el->edges = (Edge*)realloc(el->edges, e1 * sizeof(Edge));
        }
    }
    fclose(file);
    ++el->n; // So the nodes are numbered from 0 to el->n - 1
    el->edges = (Edge*)realloc(el->edges, el->e * sizeof(Edge));
    return el;
}


void Graph::enumerateKCliques(int k) {
    h_cliques.clear();
    v2cliques.clear();
    for (int i = 0; i < n; i++) v2cliques.push_back(vector<int>());

    clock_t begin = clock();
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

    clock_t end = clock();

    double io_secs = double(end - begin) / CLOCKS_PER_SEC;
    cout << "kclist time is:" << io_secs << endl;
} 

void Graph::backtrack(int currentNode, int k, vector<int>& currentClique) {
    if (currentClique.size() == k) {
        h_cliques.push_back(currentClique);
        for (int node : currentClique) {
            //cout << node << " ";
            v2cliques[node].push_back(h_cliques.size() - 1);
        }
        //cout << endl;
        return;
    }

    for (int i = currentNode; i < n; ++i) {
        if (isValidAddition(currentClique, i)) {
            currentClique.push_back(i);
            backtrack(i + 1, k, currentClique);
            currentClique.pop_back();
        }
    }
}

bool Graph::isValidAddition(const vector<int>& currentClique, int newNode) {
    for (int node : currentClique) {
        bool found = false;
        for (int neighbor : adj[node]) {
            if (neighbor == newNode) {
                found = true;
                break;
            }
        }
        if (!found) {
            return false;
        }
    }
    return true;
}

void Graph::frank_wolfe() {
    printf("#node %lld #edge %lld\n", slt_nodes.size(), slt_edges.size());
    if (CT == 0) {
        for (auto& edge : slt_edges) {
            alpha[edge] = 0.5;
        }
    }

    for (auto& node : slt_nodes) {
        r[node] = 0;
    }
    for (auto& e : slt_edges) {
        r[edges[e].first] += alpha[e];
        r[edges[e].second] += 1 - alpha[e];
    }


    for (int t = CT + 1; t <= CT + NT; t++) {
        double gamma_t = 2.0 / (t + 2);
        for (auto& e : slt_edges) {
            //            alpha[e] *= 1 - gamma_t;
            if (r[edges[e].first] < r[edges[e].second]) {
                r[edges[e].first] += gamma_t * (1 - alpha[e]);
                r[edges[e].second] -= gamma_t * (1 - alpha[e]);
                alpha[e] = alpha[e] * (1 - gamma_t) + gamma_t;
            }
            else if (r[edges[e].first] > r[edges[e].second]) {
                r[edges[e].first] -= gamma_t * alpha[e];
                r[edges[e].second] += gamma_t * alpha[e];
                alpha[e] = alpha[e] * (1 - gamma_t);
            }
        }
    }
}

void Graph::frank_wolfe_h_clique() {
    printf("#node %lu #edge %lu\n", slt_nodes.size(), slt_edges.size());
    if (CT == 0) {
        alpha_c.resize(h_cliques.size());
        for (auto& hc : slt_h_cliques) {
            for (auto& hc_node : h_cliques[hc]) {
                alpha_c[hc][hc_node] = 1 / double(h);
            }
        }
    }

    for (auto& node : slt_nodes) {
        r[node] = 0;
    }

    for (auto& hc : slt_h_cliques) {
        for (auto& hc_node : h_cliques[hc]) {
            r[hc_node] += alpha_c[hc][hc_node];
        }
    }

    for (int t = CT + 1; t <= CT + NT; t++) {
        double gamma_t = 2.0 / double(t + 2);
        for (auto& hc : slt_h_cliques) {
            for (auto& hc_node : h_cliques[hc]) {
                alpha_c[hc][hc_node] *= 1 - gamma_t;
            }
        }
        for (auto& node : slt_nodes) {
            r[node] *= 1 - gamma_t;
        }
        for (int i = 0; i < slt_h_cliques.size(); i++) {
            int hc = slt_h_cliques[i];
            int node_index = h_cliques[hc][0];
            for (auto& hc_node : h_cliques[hc]) {
                if (r[hc_node] < r[node_index]) {
                    node_index = hc_node;
                }
            }
            r[node_index] += gamma_t;
            alpha_c[hc][node_index] += gamma_t;
        }
    }
}

void Graph::frank_wolfe_h_clique_new() {
    printf("#node %lu #edge %lu\n", slt_nodes.size(), slt_edges.size());
    if (CT == 0) {
        alpha_c.resize(h_cliques.size());
        for (auto& hc : slt_h_cliques) {
            for (auto& hc_node : h_cliques[hc]) {
                alpha_c[hc][hc_node] = 1 / double(h);
            }
        }
    }

    for (auto& node : slt_nodes) {
        r[node] = 0;
    }

    for (auto& hc : slt_h_cliques) {
        for (auto& hc_node : h_cliques[hc]) {
            r[hc_node] += alpha_c[hc][hc_node];
        }
    }

    for (int t = CT + 1; t <= CT + NT; t++) {
        double gamma_t = 2.0 / double(t + 2);
        for (auto& hc : slt_h_cliques) {
            for (auto& hc_node : h_cliques[hc]) {
                alpha_c[hc][hc_node] *= 1 - gamma_t;
            }
        }
        for (auto& node : slt_nodes) {
            r[node] *= 1 - gamma_t;
        }
        for (int i = 0; i < slt_h_cliques.size(); i++) {
            int hc = slt_h_cliques[i];
            int node_index = h_cliques[hc][0];
            int sub_node_index = h_cliques[hc][0];
            for (auto& hc_node : h_cliques[hc]) {
                if (r[hc_node] < r[node_index]) {
                    sub_node_index = node_index;
                    node_index = hc_node;
                }
                else if (r[hc_node] < r[sub_node_index])
                {
                    sub_node_index = hc_node;
                }
            }
            r[node_index] = r[node_index] + (h - 1.0)/(1.0 * h) * gamma_t;
            alpha_c[hc][node_index] = alpha_c[hc][node_index] + (h - 1.0) / (1.0 * h) *  gamma_t;

            r[sub_node_index] = r[sub_node_index] + 1.0 / (1.0 * h) * gamma_t;
            alpha_c[hc][sub_node_index] = alpha_c[hc][sub_node_index] + 1.0 / (1.0 * h) * gamma_t;
        }
    }
}

void Graph::LODA_hclique() {
    printf("#node %lu #edge %lu\n", slt_nodes.size(), slt_edges.size());
    if (CT==0) {
        alpha_c.resize(h_cliques.size());
        for (auto& hc : slt_h_cliques) {
            for (auto& hc_node : h_cliques[hc]) {
                alpha_c[hc][hc_node] = 1 / double(h);
            }
        }
    }

    for (auto& node : slt_nodes) {
        r[node] = 0;
    }

    for (auto& hc : slt_h_cliques) {
        for (auto& hc_node : h_cliques[hc]) {
            r[hc_node] += alpha_c[hc][hc_node];
        }
    }
    for (int t = 1; t <= NT; t++) {
        for (auto& hc : slt_h_cliques) {
            double gamma_avg = 0.0;
            for (auto& hc_node : h_cliques[hc]) {
                //if (IsEqual(alpha_c[hc][hc_node], 0))continue;
                gamma_avg += r[hc_node];
            }
            gamma_avg /= 1.0 * h;
            double sumD = 0.0;
            for (auto& hc_node : h_cliques[hc]) {
                double dif = 0.0;
                //if (IsEqual(alpha_c[hc][hc_node], 0))continue;
                if (r[hc_node] > gamma_avg)
                {
                    dif = r[hc_node] - gamma_avg;
                }
                else continue;
                if (dif > alpha_c[hc][hc_node])
                {
                    dif = alpha_c[hc][hc_node];
                }
                alpha_c[hc][hc_node] -= dif;
                r[hc_node] -= dif;
                sumD +=dif;
            }
            for (auto& hc_node : h_cliques[hc])
            {
                double dif = 0.0;
                if (r[hc_node] < gamma_avg)
                {
                    dif = gamma_avg - r[hc_node];
                    if (dif > 1 - alpha_c[hc][hc_node])
                    {
                        dif = 1 - alpha_c[hc][hc_node];
                    }
                    if (dif > sumD)
                    {
                        dif = sumD;
                    }
                    alpha_c[hc][hc_node] += dif;
                    r[hc_node] += dif;
                    sumD -= dif;
                    if (sumD == 0)break;
                }
            }
        }
    }
}

void Graph::new_LODA()
{
    printf("#node %lu #edge %lu\n", slt_nodes.size(), slt_edges.size());
    if (CT == 0) {
        alpha_c.resize(h_cliques.size());
        for (auto& hc : slt_h_cliques) {
            for (auto& hc_node : h_cliques[hc]) {
                alpha_c[hc][hc_node] = 1 / double(h);
            }
        }
    }

    for (auto& node : slt_nodes) {
        r[node] = 0;
    }

    for (auto& hc : slt_h_cliques) {
        for (auto& hc_node : h_cliques[hc]) {
            r[hc_node] += alpha_c[hc][hc_node];
        }
    }
    for (int t = 1; t <= NT; t++) {
        for (auto& hc : slt_h_cliques) {
            double max_gamma = 0.0, max_pos = -1;
            double min_gamma = 1e9, min_pos = -1;
            for (auto& hc_node : h_cliques[hc]) {
                if (r[hc_node] > max_gamma)
                {
                    max_gamma = r[hc_node];
                    max_pos = hc_node;
                }
                if (r[hc_node] < min_gamma)
                {
                    min_gamma = r[hc_node];
                    min_pos = hc_node;
                }
            }
            if (max_gamma < min_gamma + DELTA) continue;
            double d = min((max_gamma - min_gamma) / 2.0, alpha_c[hc][max_pos]);
            r[max_pos] -= d;
            alpha_c[hc][max_pos] -= d;
            r[min_pos] += d;
            alpha_c[hc][min_pos] += d;
        }
    }
}

void Graph::pava() {
   
    sort(slt_nodes.begin(), slt_nodes.end(), [this](int a, int b)->bool {
        return r[a] > r[b];
        });


    for (int i = 0; i < slt_nodes.size(); i++) {
        deg[slt_nodes[i]] = i;
    }


    vector<int> ne(slt_nodes.size(), 0);
    for (auto e : slt_edges) {
        int u = deg[edges[e].first];
        int v = deg[edges[e].second];
        ne[(u > v) ? u : v]++;
    }


    nag[0] = 1;
    val[0] = ne[0];
    nsg = 0;


    for (int i = 1; i < slt_nodes.size(); i++) {
        nsg += 1;            
        val[nsg] = ne[i];    
        nag[nsg] = 1;        

        while ((nsg > 0) && (val[nsg] > val[nsg - 1] - 1e-5)) {
          
            val[nsg - 1] = (nag[nsg] * val[nsg] + nag[nsg - 1] * val[nsg - 1]) / (nag[nsg] + nag[nsg - 1]);
            nag[nsg - 1] += nag[nsg]; 
            nsg--;                    
        }
    }
    nsg++; 


    printf("nsg %d\n", nsg);

    int cur = 0;
    for (int i = 0; i < nsg; i++) {
        for (int j = cur; j < cur + nag[i]; j++) {
            sg[slt_nodes[j]] = i;
        }
        cur += nag[i];
    }

}

int max(vector<int>* h_clique) {
    int max_node = 0;
    for (auto& hc_node : *h_clique) {
        if (hc_node > max_node) max_node = hc_node;
    }
    return max_node;
}


int max(vector<int>* sg, vector<int> *h_clique) {
    int max_sg = 0;
    for (auto& hc_node : *h_clique) {
        if (sg->at(hc_node) > max_sg) max_sg = sg->at(hc_node);
    }
    return max_sg;
}

int min(vector<int>* h_clique) {
    int min_node = INT_MAX;
    for (auto& hc_node : *h_clique) {
        if (hc_node < min_node) min_node = hc_node;
    }
    return min_node;
}


int min(vector<int>* sg, vector<int>* h_clique) {
    int min_sg = INT_MAX;
    for (auto& hc_node : *h_clique) {
        if (sg->at(hc_node) < min_sg) min_sg = sg->at(hc_node);
    }
    return min_sg;
}

void Graph::pava_h_clique() {

    sort(slt_nodes.begin(), slt_nodes.end(), [this](int a, int b)->bool {
        return r[a] > r[b];
        });


    for (int i = 0; i < slt_nodes.size(); i++) {
        deg[slt_nodes[i]] = i;
    }

  
    vector<int> nc(slt_nodes.size(), 0);
    for (auto& hc : slt_h_cliques) {
        int max_node = 0;
        for (auto& hc_node : h_cliques[hc]) {
            if (deg[hc_node] > max_node) max_node = deg[hc_node];
        }
        nc[max_node]++;
    }



    nag[0] = 1;
    val[0] = nc[0];
    nsg = 0;

    for (int i = 1; i < slt_nodes.size(); i++) {
        nsg += 1;           
        val[nsg] = nc[i];    
        nag[nsg] = 1;        

        while ((nsg > 0) && (val[nsg] > val[nsg - 1] - 1e-5)) {

            val[nsg - 1] = (nag[nsg] * val[nsg] + nag[nsg - 1] * val[nsg - 1]) / (nag[nsg] + nag[nsg - 1]);
            nag[nsg - 1] += nag[nsg]; 
            nsg--;                    
        }
    }
    nsg++;  

    printf("nsg %d\n", nsg);

    int cur = 0;
    for (int i = 0; i < nsg; i++) {
        for (int j = cur; j < cur + nag[i]; j++) {
            sg[slt_nodes[j]] = i;
        }
        cur += nag[i];
    }
}


void Graph::check_sg_h_clique() {
  
    if (nsg <= 1) return;

    for (auto hc : slt_h_cliques) {
        sort(h_cliques[hc].begin(), h_cliques[hc].end(), [this](int x, int y)->bool {
            return sg[x] < sg[y];
            });
    }
    sort(slt_h_cliques.begin(), slt_h_cliques.end(), [this](int a, int b)->bool {
        for (int i = 0; i < h; i++) {
            if (sg[h_cliques[a][i]] != sg[h_cliques[b][i]]) {
                return sg[h_cliques[a][i]] < sg[h_cliques[b][i]];
            }
            
        }
        return false;
    });
   
    vector<bool> valid(nsg);
    vector<int> bin(nsg + 1, 0);

    for (auto& hc : slt_h_cliques) {  
        int a = sg[h_cliques[hc][0]];
        ++bin[a];
    }

    int s = 0;
    for (int i = 0; i <= nsg; i++) {
        int tmp = s;
        s += bin[i];
        bin[i] = tmp;
    }

    vector<double> max_r(nsg);

    int cur = 0;
    for (int i = 0; i < nsg; i++) {
        max_r[i] = r[slt_nodes[cur]];

        for (int j = cur + 1; j < cur + nag[i]; j++) {
            max_r[i] = max(max_r[i], r[slt_nodes[j]]);
        }
        cur += nag[i];
    }

    double min_r = r[slt_nodes[0]]; 
    cur = 0;
    vector<int> cpt(bin);
    clock_t start = clock();

    for (int i = 0; i < nsg - 1; i++) { 
        for (int j = cur; j < cur + nag[i]; j++) { 
            min_r = min(min_r, r[slt_nodes[j]]);
        }
        double min_t = min_r; 

        vector<double> max_tmp(max_r.begin() + i + 1, max_r.end()); 
        double max_t = *max_element(std::begin(max_tmp), std::end(max_tmp)); 

        queue<std::tuple<int, int, std::unordered_map<int, double>>> q; 
        valid[i] = true; 

        for (int j = 0; j <= i; j++) { 

            while (cpt[j] < bin[j + 1]) { 
                int hc = slt_h_cliques[cpt[j]]; 
                int max_sg = 0;
                for (auto& hc_node : h_cliques[hc]) {
                    if (sg[hc_node] > max_sg) max_sg = sg[hc_node];
                }


                if (max_sg > i) { 
                    q.push(make_tuple(j, cpt[j], alpha_c[hc]));

                    double s = 0;
                    int num = 0;
                    for (auto& hc_node : h_cliques[hc]) {
                        if (sg[hc_node] != max_sg) {
                            s += alpha_c[hc][hc_node];
                            r[hc_node] -= alpha_c[hc][hc_node];
                            min_t = min(min_t, r[hc_node]);
                            alpha_c[hc][hc_node] = 0;
                            num++;
                        }
                    }
                    for (auto& hc_node : h_cliques[hc]) {
                        if (sg[hc_node] == max_sg) {
                            r[hc_node] += s / num;
                            max_r[sg[hc_node]] = max(max_r[sg[hc_node]], r[hc_node]);
                            max_t = max(max_t, r[hc_node]);
                            alpha_c[hc][hc_node] += s / num;
                        }
                    }
                    if (min_t <= max_t) {
                        valid[i] = false; 
                        break;
                    }
                }

                ++cpt[j]; 
            }

            if (!valid[i]) break;
        }

        if (valid[i]) {
            min_r = min(min_t, min_r); 
        }
        else {
            for (int j = i + 1; j < max_r.size(); j++) {
                max_r[j] = max_tmp[j - i - 1];
            }

            while (!q.empty()) {
                auto tp = q.front();
                q.pop();
                cpt[get<0>(tp)] = min(cpt[get<0>(tp)], get<1>(tp)); 

                int hc = slt_h_cliques[get<1>(tp)];

                int max_sg = 0;
                for (auto& hc_node : h_cliques[hc]) {
                    if (sg[hc_node] > max_sg) max_sg = sg[hc_node];
                }

                double s = 0;
                int num = 0;
                for (auto& hc_node : h_cliques[hc]) {
                    if (sg[hc_node] != max_sg) {
                        s += get<2>(tp)[hc_node];
                        r[hc_node] += get<2>(tp)[hc_node];
                        alpha_c[hc][hc_node] = get<2>(tp)[hc_node];
                        num++;
                    }
                }
                for (auto& hc_node : h_cliques[hc]) {
                    if (sg[hc_node] == max_sg) {
                        r[hc_node] -= s / num;
                        alpha_c[hc][hc_node] -= s / num;
                    }
                }
            }
        }

        clock_t cur_t = clock(); 
    }


    check_first = (nsg == 1) || valid[0];

    if (nsg > 1) {
        vector<int> n_nag;
        for (int i = 0; i < nsg; i++) {

            if (i == 0 || valid[i - 1]) {

                n_nag.push_back(nag[i]);
            }
            else {
             
                n_nag[n_nag.size() - 1] += nag[i];
            }

        }
        nag = n_nag; 

        cur = 0; 
        nsg = nag.size(); 
        double minr = m;


        for (int i = 0; i < nsg; i++) {
            double tmp_r = m; 
            for (int j = cur; j < cur + nag[i]; j++) {
               
                sg[slt_nodes[j]] = i; 
                tmp_r = min(tmp_r, r[slt_nodes[j]]); 
                minr = min(minr, r[slt_nodes[j]]); 
            }
            cur += nag[i]; 
        }
    }


    printf("updated nsg %d\n", nsg);
}

void Graph::check_sg() {

    if (nsg <= 1) return;

    sort(slt_edges.begin(), slt_edges.end(), [this](int a, int b)->bool {
        int a1 = min(sg[edges[a].first], sg[edges[a].second]);
        int a2 = max(sg[edges[a].first], sg[edges[a].second]);
        int b1 = min(sg[edges[b].first], sg[edges[b].second]);
        int b2 = max(sg[edges[b].first], sg[edges[b].second]);
        if (a1 != b1)
            return a1 < b1;
        else
            return a2 < b2;
        });


    vector<bool> valid(nsg);


    vector<int> bin(nsg + 1, 0);


    for (int i = 0; i < slt_edges.size(); i++) {

        int a = min(sg[edges[slt_edges[i]].first], sg[edges[slt_edges[i]].second]);
        ++bin[a];
    }

    int s = 0;
    for (int i = 0; i <= nsg; i++) {
        int tmp = s;
        s += bin[i];
        bin[i] = tmp;
    }


    vector<double> max_r(nsg);

    int cur = 0;
    for (int i = 0; i < nsg; i++) {

        max_r[i] = r[slt_nodes[cur]];


        for (int j = cur + 1; j < cur + nag[i]; j++) {
            max_r[i] = max(max_r[i], r[slt_nodes[j]]);
        }
        cur += nag[i];
    }


    double min_r = r[slt_nodes[0]]; 
    cur = 0; 
    vector<int> cpt(bin); 
    clock_t start = clock(); 

    for (int i = 0; i < nsg - 1; i++) { 
        for (int j = cur; j < cur + nag[i]; j++) { 
        }
        double min_t = min_r; 

        vector<double> max_tmp(max_r.begin() + i + 1, max_r.end()); 
        double max_t = *max_element(std::begin(max_tmp), std::end(max_tmp)); 

        queue<std::tuple<int, int, double>> q; 
        valid[i] = true; 

        for (int j = 0; j <= i; j++) { 
            while (cpt[j] < bin[j + 1]) { 
                int e = slt_edges[cpt[j]]; 

                if (max(sg[edges[e].first], sg[edges[e].second]) > i) { 
                    q.push(make_tuple(j, cpt[j], alpha[e]));

                    int u = edges[e].first;
                    int v = edges[e].second;

                    if (sg[u] > sg[v]) {
                        r[u] += 1 - alpha[e];
                        r[v] -= 1 - alpha[e];
                        max_r[sg[u]] = max(max_r[sg[u]], r[u]);
                        max_t = max(max_t, r[u]);
                        min_t = min(min_t, r[v]);
                        alpha[e] = 1;
                    }
                    else {
                        r[u] -= alpha[e];
                        r[v] += alpha[e];
                        min_t = min(min_t, r[u]);
                        max_r[sg[v]] = max(max_r[sg[v]], r[v]);
                        max_t = max(max_t, r[v]);
                        alpha[e] = 0;
                    }

                    if (min_t <= max_t) {
                        valid[i] = false; 
                        break;
                    }
                }

                ++cpt[j]; 
            }

            if (!valid[i]) break; 
        }

        if (valid[i]) {
            min_r = min(min_t, min_r);
        }
        else {
            
            for (int j = i + 1; j < max_r.size(); j++) {
                max_r[j] = max_tmp[j - i - 1];
            }

            while (!q.empty()) {
                auto tp = q.front();
                q.pop();
                cpt[get<0>(tp)] = min(cpt[get<0>(tp)], get<1>(tp));
                int e = slt_edges[get<1>(tp)];
                int u = edges[e].first, v = edges[e].second;
                r[u] += get<2>(tp) - alpha[e];
                r[v] -= get<2>(tp) - alpha[e];
                alpha[e] = get<2>(tp);
            }
        }

        clock_t cur_t = clock();
    }


    
    check_first = (nsg == 1) || valid[0];


    if (nsg > 1) {
        vector<int> n_nag;
        for (int i = 0; i < nsg; i++) {
          
            if (i == 0 || valid[i - 1]) {

                n_nag.push_back(nag[i]);
            }
            else {
               
                n_nag[n_nag.size() - 1] += nag[i];
            }

        }
        nag = n_nag; 

        cur = 0; 
        nsg = nag.size(); 
        double minr = m; 

    
        for (int i = 0; i < nsg; i++) {
            double tmp_r = m; 

      
            for (int j = cur; j < cur + nag[i]; j++) {
                sg[slt_nodes[j]] = i; 
                tmp_r = min(tmp_r, r[slt_nodes[j]]); 
                minr = min(minr, r[slt_nodes[j]]); 
            }
            cur += nag[i]; 
        }
    }


    printf("updated nsg %d\n", nsg);
}




void Graph::pruning() {  
    vector<double> min_r(nsg);
    vector<double> max_r(nsg);

    int cur = 0;
    for (int i = 0; i < nsg; i++) {
        min_r[i] = r[slt_nodes[cur]];
        max_r[i] = r[slt_nodes[cur]];
        for (int j = cur + 1; j < cur + nag[i]; j++) {
            min_r[i] = min(min_r[i], r[slt_nodes[j]]);
            max_r[i] = max(max_r[i], r[slt_nodes[j]]);
        }
        for (int j = cur; j < cur + nag[i]; j++) {
            /* int threadId = omp_get_thread_num();
             printf("Thread %d is processing j=%d\n", threadId, j);*/
            selected[slt_nodes[j]] = true;
            deg[slt_nodes[j]] = 0;
            rho_u[slt_nodes[j]] = min(rho_u[slt_nodes[j]], max_r[i]);
            rho_l[slt_nodes[j]] = max(rho_l[slt_nodes[j]], min_r[i]);
            if (rho_u[slt_nodes[j]] + DELTA < rho_l[slt_nodes[j]]) {
                selected[slt_nodes[j]] = false;
            }
            if (active[slt_nodes[j]]) {
                rho_gu[slt_nodes[j]] = min(rho_gu[slt_nodes[j]], max_r[i]);
            }
        }


        cur += nag[i];
    }


    for (auto e : slt_edges) {
        int u = edges[e].first;
        int v = edges[e].second;
        if (sg[u] > sg[v]) {
            selected[u] = false;
        }
        if (sg[u] < sg[v]) {
            selected[v] = false;
        }
    }

    int u;
    for (int i = 0; i < slt_nodes.size(); i++) {
        u = slt_nodes[i];
        if (selected[u]) {
            for (auto v : adj[u]) {
                if (rho_l[v] > rho_u[u] + DELTA) {
                    selected[u] = false;
                }
            }
        }
    }

    for (auto& hc : slt_h_cliques) {
        int flag = 1;
        for (auto& hc_node : h_cliques[hc]) {
            if (!selected[hc_node]) {
                flag = 0;
                break;
            }
        }
        if (flag) {
            for (auto& hc_node : h_cliques[hc]) {
                ++deg[hc_node];
            }
        }
    }

  

    queue<int> q;
    for (auto u : slt_nodes) {
        if (selected[u] && deg[u] < rho_l[u] - DELTA) {
            selected[u] = false;
            q.push(u);
        }
    }
    std::unordered_map<int, int> slt_h_cliques_map;
    for (int i = 0; i < slt_h_cliques.size(); i++) {
        slt_h_cliques_map.insert({ slt_h_cliques[i], 1 });
    }
    while (!q.empty()) {
        int u = q.front(); q.pop();
        for (int i = 0; i < v2cliques[u].size(); i++) {
            if (slt_h_cliques_map.find(v2cliques[u][i]) == slt_h_cliques_map.end()) continue;
        
            vector<int> hc = h_cliques[v2cliques[u][i]];
            for (int j = 0; j < hc.size(); j++) {
                int v = hc[j];
                if (selected[v])
                {
                    if (--deg[v] < rho_l[v] - DELTA)
                    {
                        selected[v] = false;
                        q.push(v);
                    }
                }
            }
        }
    }

    vector<int> tmp_nodes;
    vector<bool> inactive_sg(nsg, false);
    for (auto u : slt_nodes) {
        if (selected[u])
            tmp_nodes.push_back(u);
        else
            inactive_sg[sg[u]] = true;
    }
    slt_nodes = tmp_nodes;

    for (auto u : slt_nodes) {
        if (inactive_sg[sg[u]]) {
            active[u] = false;
        }
    }

    vector<int> tmp_edges;
    for (auto e : slt_edges) {
        if (selected[edges[e].first] && selected[edges[e].second]) {
            tmp_edges.push_back(e);
        }
    }


    slt_edges = tmp_edges;
    sort(slt_edges.begin(), slt_edges.end(), [this](int a, int b)->bool {
        return sg[edges[a].first] < sg[edges[b].first];
        });


    vector<int> tmp_h_cliques;
    for (auto hc : slt_h_cliques) {
        int flag = 1;
        for (auto hc_node : h_cliques[hc]) {
            if (!selected[hc_node]) {
                flag = 0;
                break;
            }
        }
        if (flag == 1) tmp_h_cliques.push_back(hc);
    }

    slt_h_cliques = tmp_h_cliques;

    sort(slt_h_cliques.begin(), slt_h_cliques.end(), [this](int a, int b)->bool {
        for (int i = 0; i < h; i++) {
            if (sg[h_cliques[a][i]] != sg[h_cliques[b][i]]) {
                return sg[h_cliques[a][i]] < sg[h_cliques[b][i]];
            }
            
        }
        return false;
        });

    while (slt_nodes.size() > 0 && sg[slt_nodes.back()] > 0) {
        int n_nodes = slt_nodes.size() - 1;
        while (n_nodes > 0 && sg[slt_nodes[n_nodes]] == sg[slt_nodes.back()]) {
            --n_nodes;
        }
        n_nodes++;
        vector<int> t_nodes(slt_nodes.begin() + n_nodes, slt_nodes.end());
        stk_nodes.push(t_nodes);
        slt_nodes.resize(n_nodes);


        if (sg[t_nodes[0]] > sg[edges[slt_edges.back()].first]) {
            stk_nodes.pop();
            continue;
        }
        int n_edges = slt_edges.size() - 1;
        while (n_edges > 0 && sg[edges[slt_edges[n_edges]].first] == sg[edges[slt_edges.back()].first]) {
            --n_edges;
        }
        n_edges++;
        vector<int> t_edges(slt_edges.begin() + n_edges, slt_edges.end());
        stk_edges.push(t_edges);
        slt_edges.resize(n_edges);

        if ( sg[t_nodes[0]] > sg[h_cliques[slt_h_cliques.back()][0]]) {
            stk_nodes.pop();
            stk_edges.pop();
            continue;
        }

        int n_h_cliques = slt_h_cliques.size() - 1;

        cout << min(&sg, &h_cliques[slt_h_cliques.back()]) << endl;
        while (n_h_cliques > 0 && min(&sg, &h_cliques[slt_h_cliques[n_h_cliques]]) == min(&sg, &h_cliques[slt_h_cliques.back()])) {
            --n_h_cliques;
        }
        n_h_cliques++;
        vector<int> t_h_cliques(slt_h_cliques.begin() + n_h_cliques, slt_h_cliques.end());
        //cout << "cliques_num:" << t_h_cliques.size() << endl;
        stk_h_cliques.push(t_h_cliques);
        slt_h_cliques.resize(n_h_cliques);


        printf("sg %d %d %lu %lu\n", sg[t_nodes[0]], sg[edges[t_edges[0]].first], t_nodes.size(), t_edges.size());
        stk_CT.push(CT);
    }
}

bool Graph::Get_Maxflow() {
    Network network;
    vector<int> id_in_network(n, -1);
    vector<int> idx_in_slt_nodes(slt_nodes.size(), 0);
    unsigned long long m = 0;
    vector<Network::Vertex> R;
    Network::Vertex s = network.AddVertex(), t = network.AddVertex();
    for (unsigned i = 0; i < slt_nodes.size(); ++i) {
        id_in_network[slt_nodes[i]] = i;
        R.push_back(network.AddVertex());
    }


    for (auto hc : slt_h_cliques) {
        ++m;
        Network::Vertex v = network.AddVertex();
        for (auto hc_node : h_cliques[hc]) {
            network.AddEdge(v, R[id_in_network[hc_node]], slt_nodes.size());
        }
        network.AddEdge(s, v, slt_nodes.size());
    }

    for (unsigned i = 0; i < slt_nodes.size(); ++i)
        network.AddEdge(R[i], t, m);
    return network.MaxFlow(s, t) >= m * slt_nodes.size();
}

void Graph::findLhCDS() {

    //#pragma omp threadprivate(h_cliques)
    for (int i = 0; i < n; i++) {
        if (p == 0)
        {
            if (adj[i].size() < h-1)
            {
                selected[i] = false;
                continue;
            }
        }
        slt_nodes.push_back(i);
        /*selected.push_back(true);*/
    }
    for (int i = 0; i < m; i++) {
        if (p == 0)
        {
            if (adj[edges[i].first].size() < h - 1 || adj[edges[i].second].size() < h-1)
            {
                continue;
            }
        }
        slt_edges.push_back(i);
    }
    nodes = slt_nodes;
    
    double fw_time = 0;
    double sg_time = 0;
    double pr_time = 0;
    double mf_time = 0;
    int i = 0;
    clock_t start = clock();
    enumerateKCliques(h);
    printf("#slt node %lu #slt edge %lu #h_cliques %lu\n", slt_nodes.size(), slt_edges.size(), h_cliques.size());
    CT = 0;
    clock_t t_cc1 = clock();
    selected.clear();
    compute_h_clique_core();
    prune_by_core();
    clock_t t_cc2 = clock();
    printf("compute core and prune time: %.4f\n", double(t_cc2 - t_cc1) / CLOCKS_PER_SEC);


    while (!slt_nodes.empty() || !stk_nodes.empty()) {
        //std::cout << "Iteration: " << i++ << "  " << "slt_nodes.size(): " << slt_nodes.size() << "  " << "stk_nodes.size() : " << stk_nodes.size() << std::endl;
        clock_t start = clock();
        if (slt_nodes.empty()) {
            slt_nodes = stk_nodes.top(); stk_nodes.pop();
            slt_edges = stk_edges.top(); stk_edges.pop();
            slt_h_cliques = stk_h_cliques.top(); stk_h_cliques.pop();
            CT = stk_CT.top(); stk_CT.pop();
        }
        //frank_wolfe();
        frank_wolfe_h_clique(); //frank_wolfe Algorithm. 
        //frank_wolfe_h_clique_new();
        //LODA_hclique();
        clock_t t_fw = clock();
        fw_time += double(t_fw - start) / CLOCKS_PER_SEC;
        printf("fw time: %.4f\n", double(t_fw - start) / CLOCKS_PER_SEC);
        CT += NT;
        pava_h_clique();  //Pool Adjacent Violators Algorithm. 
        clock_t t_pv = clock();
        printf("pava time: %.4f\n", double(t_pv - t_fw) / CLOCKS_PER_SEC);
        //check_sg();
        check_sg_h_clique();
        clock_t t_check_sg = clock();
        sg_time += double(t_check_sg - t_fw) / CLOCKS_PER_SEC;
        printf("check sg time: %.4f\n", double(t_check_sg - t_pv) / CLOCKS_PER_SEC);
        pruning();  //pruning Algorithm. 
        clock_t t_prune = clock();
        pr_time += double(t_prune - t_check_sg) / CLOCKS_PER_SEC;
        printf("pruning time: %.4f\n", double(t_prune - t_check_sg) / CLOCKS_PER_SEC);
        sort(slt_nodes.begin(), slt_nodes.end());

        if (!slt_nodes.empty()) {     // && check_first     
            //double g = (double) slt_h_cliques.size() / slt_nodes.size(); 
            vector<pair<int, int>> tmp_edges;
            for (auto e : slt_edges) {
                tmp_edges.push_back(edges[e]);
            }
            //TODO check by edge number? and connected component

           /* FlowNetwork fn = FlowNetwork(tmp_edges, g);
            double max_flow = fn.get_maxflow(0, fn.n-1);
            std::cout << "abs(max_flow - slt_edges.size()): " << abs(max_flow - slt_h_cliques.size()) << std::endl;*/
            // IsDensestMF
            bool flag = Get_Maxflow();
            //if (abs(max_flow - slt_h_cliques.size()) <= 1e-3) {
            if (flag) {
                connected_components();
                int cur_u = 0;
                int cur_e = 0;
                int cur_hc = 0;
                for (auto pr : cmpt) {
                    vector<int> tmp_nodes(slt_nodes.begin() + cur_u, slt_nodes.begin() + std::get<0>(pr));
                    vector<int> tmp_cliques(slt_h_cliques.begin() + cur_hc, slt_h_cliques.begin() + std::get<2>(pr));
                    printf("lhcdses candidate: #nodes %lu #edges %d #h-cliques %d\n", tmp_nodes.size(), std::get<1>(pr) - cur_e, std::get<2>(pr) - cur_hc);
                    double g = (double)(std::get<2>(pr) - cur_hc) / tmp_nodes.size();

                    if (lhcdses.size() == 0 || verify_LhCDS_small(tmp_nodes, g, tmp_cliques)) {
                        lhcdses.push_back(tmp_nodes);
                        lhcds_rho.push_back(g);
                        vector<int> tmp_edges(slt_edges.begin() + cur_e, slt_edges.begin() + std::get<1>(pr));
                        lhcds_edge.push_back(tmp_edges);
                        lhcds_den.push_back((double)(std::get<1>(pr) - cur_e) / tmp_nodes.size());
                        if (lhcdses.size() >= topk)
                            break;
                    }

                    cur_u = std::get<0>(pr);
                    cur_e = std::get<1>(pr);
                    cur_hc = std::get<2>(pr);
                }
                //lhcdses.push_back(slt_nodes);
                if (lhcdses.size() >= topk)
                    break;
                slt_nodes.clear();
                slt_edges.clear();
            }
        }
        clock_t t_verify_LhCDS = clock();
        mf_time += double(t_verify_LhCDS - t_prune) / CLOCKS_PER_SEC;
        printf("verifyLhCDS time: %.4f\n", double(t_verify_LhCDS - t_prune) / CLOCKS_PER_SEC);
    }

    for (int i = 0; i < lhcdses.size(); i++) {
        sort(lhcdses[i].begin(), lhcdses[i].end());
        for (int j = 0; j < lhcdses[i].size(); j++) {
            cout << lhcdses[i][j] << " ";
        }
        cout << "nodes_num:" << lhcdses[i].size() << " ";
        cout << "edges_num:" << lhcds_edge[i].size() << " ";
        cout << "cli-den:" << lhcds_rho[i] << " ";
        cout << "density: " << lhcds_den[i] << endl;
    }

    SmallGraph sg;
  
    for (int i = 0; i < lhcdses.size(); i++) {
        for (int j = 0; j < lhcds_edge[i].size(); j++) {
            sg.addEdge(edges[lhcds_edge[i][j]].first, edges[lhcds_edge[i][j]].second);
        } 
    }
    double avgClusteringCoeff = sg.averageClusteringCoefficient();
    cout << "Average Clustering Coefficient: " << avgClusteringCoeff << endl;


    printf("fw %.4f sec, sg %.4f sec, pr %.4f sec, mf %.4f sec\n", fw_time, sg_time, pr_time, mf_time);

}

void Graph::compute_h_clique_core() {
    max_d = 0;
    for (int i = 0; i < n; i++) {
        deg[i] = v2cliques[i].size();
        max_d = max(deg[i], max_d);
    }
    printf("max_d %d\n", max_d);
    vector<int> bin(max_d + 1, 0);
    pos.resize(n);
    vector<int> vert(n);
    for (int i = 0; i < n; i++) {
        ++bin[deg[i]];
    }
    int start = 0;
    for (int d = 0; d <= max_d; d++) {
        int num = bin[d];
        bin[d] = start;
        start += num;
    }
    for (int i = 0; i < n; i++) {
        pos[i] = bin[deg[i]]++;
        vert[pos[i]] = i;
    }
    for (int d = max_d; d > 0; d--) {
        bin[d] = bin[d - 1];
    }
    bin[0] = 0;
    bool* vis_clique = new bool[h_cliques.size()];
    memset(vis_clique, false, sizeof(bool) * h_cliques.size());
    for (int i = 0; i < n; i++) {
        int u = vert[i];
        for (int k = 0; k < v2cliques[u].size(); k++) {
            if (vis_clique[v2cliques[u][k]])continue;
            vis_clique[v2cliques[u][k]] = true;
            vector<int> hc = h_cliques[v2cliques[u][k]];
            for (int j = 0; j < hc.size(); j++) {
                int v = hc[j];
                if (deg[v] > deg[u])
                {
                    int dv = deg[v], pv = pos[v];
                    int pw = bin[dv], w = vert[pw];
                    if (v != w) {
                        pos[v] = pw; vert[pv] = w;
                        pos[w] = pv; vert[pw] = v;
                    }
                    ++bin[dv]; --deg[v];
                }
            }
        }
    }
}



void Graph::prune_by_core() {

    for (int u = 0; u < n; u++) {
        rho_l[u] = deg[u] / double(h);
        rho_u[u] = deg[u];
        rho_gu[u] = deg[u];
    }
    vector<double> tmp_rho_l(n);
    queue<int> q;
    for (int u = 0; u < n; u++) {
        tmp_rho_l[u] = rho_l[u];
        for (auto v : adj[u]) {
            tmp_rho_l[u] = max(tmp_rho_l[u], rho_l[v]);
        }
        if (tmp_rho_l[u] > rho_u[u] + DELTA) {
            selected[u] = false;
            q.push(u);
        }
    }


    vector<int> vtx;
    while (!q.empty()) {
        int u = q.front(); q.pop();
        vtx.push_back(u);
        for (int i = 0; i < v2cliques[u].size(); i++) {
            vector<int> hc = h_cliques[v2cliques[u][i]];
            for (int j = 0; j < hc.size(); j++) {
                int v = hc[j];
                if (selected[v] && deg[u] > deg[v]) {
                    --deg[v];
                    if (tmp_rho_l[v] > deg[v] + DELTA) {
                        selected[v] = false;
                        q.push(v);
                    }
                }
            }
        }
    }

    slt_nodes.clear();
    slt_edges.clear();

    vector<int> selected_cliques(h_cliques.size(), 1);
    for (int i = 0; i < n; i++) {
        if (selected[i]) slt_nodes.push_back(i);
        else if (!selected[i]) {
            for (int j = 0; j < v2cliques[i].size(); j++) {
                selected_cliques[v2cliques[i][j]] = 0;
            }
        }
    }
    for (int i = 0; i < m; i++) {
        if (selected[edges[i].first] && selected[edges[i].second]) {
            slt_edges.push_back(i);
        }
    }

    for (int i = 0; i < h_cliques.size(); i++) {
        if (selected_cliques[i]) {
            slt_h_cliques.push_back(i);
        }
    }
    printf("#slt node %lu #slt edge %lu #slt_h_cliques %lu\n", slt_nodes.size(), slt_edges.size(), slt_h_cliques.size());

    sort(vtx.begin(), vtx.end(), [this](int a, int b)->bool {
        return rho_gu[a] > rho_gu[b];
        });
    for (auto u : vtx) {
        if (active[u]) {
            active[u] = false;
            double threshold = rho_gu[u] + DELTA;
            q.push(u);
            while (!q.empty()) {
                int v = q.front(); q.pop();
                for (auto w : adj[v]) {
                    if (active[w] && rho_l[w] <= threshold) {
                        active[w] = false;
                        q.push(w);
                    }
                }
            }
        }
    }

    int cnt = 0;
    for (int i = 0; i < n; i++) {
        if (active[i]) ++cnt;
    }
    printf("#active nodes %d\n", cnt);
}


bool Graph::verify_by_maxflow(vector<int>& nodes, double g){

}

bool Graph::verify_LhCDS(vector<int>temp_nodes, vector<int>& nodes, double g) {
    for (auto u : temp_nodes) {
        for (auto v : adj[u]) {
            if (rho_l[v] > g) {
                printf("validate failed\n");
                return false;
            }
        }
    }

    Network network;
    vector<int> id_in_network(n, -1);
    vector<int> idx_in_slt_nodes(nodes.size(), 0);
    unsigned long long m = 0;
    vector<Network::Vertex> R;
    //vector<int> S;
    double g_temp = g - 1.0 / n / n;

    //Add a source node s and a sink node
    Network::Vertex s = network.AddVertex(), t = network.AddVertex();
    //Add vertices in G
    for (unsigned i = 0; i < nodes.size(); ++i) {
        id_in_network[nodes[i]] = i;
        R.push_back(network.AddVertex());
    }

    for (unsigned i = 0; i < nodes.size(); ++i) {
        network.AddEdge(s, R[i], v2cliques[nodes[i]].size());
        network.AddEdge(R[i], t, g_temp * h);
    }

    for (auto hc : h_cliques) {
        Network::Vertex v = network.AddVertex();
        for (auto hc_node : hc) {
            network.AddEdge(R[id_in_network[hc_node]], v, 1);
            network.AddEdge(v, R[id_in_network[hc_node]], h - 1);
        }
    }

    network.MaxFlow(s, t);
    //S.clear();
    std::unordered_map<int, int> S;
    std::vector<bool> visited(nodes.size() + h_cliques_no.size() + 2, false);
    network.dfs(0, visited);
    for (int i = 0; i < visited.size(); ++i) {
        if (visited[i]) {
            S.insert({ i,1 });
        }
    }
   
    std::unordered_map<int, int> mp;

    for (int i = 0; i < temp_nodes.size(); i++)
    {
        mp.insert({ temp_nodes[i],1 });
    }
    for (int i = 0; i < temp_nodes.size(); i++)
    {
        int u = temp_nodes[i];
        for (auto w : adj[u])
        {
            if (u == 1144 && w == 2693) {
                cout << "here" << endl;
            }
            if (mp.find(w) == mp.end() && id_in_network[w] != -1 && S.find(R[id_in_network[w]]) != S.end()) //w ²»ÊôÓÚnodes£¬µ«ÊÇÊôÓÚS
            {

                //cout << u << " "<< w;
                return false;
            }
        }
    }


    return true;

}

bool Graph::verify_LhCDS_small(vector<int>& nodes, double g, vector<int>& tmp_hcliques) {

    for (auto u : nodes) {
        for (auto v : adj[u]) {
            if (rho_l[v] > g) {
                printf("validate failed\n");
                return false;
            }
        }
    }
    for (auto u : nodes) {
        lhcds_num[u] = lhcdses.size();
    }
    bool flag = true;

    ++num_verify;
    queue<int> q;
    vector<int> tmp_nodes = nodes;
    std::unordered_map<int, int> vis_tmp_nodes;
    for (auto tn : nodes)
    {
        vis_tmp_nodes.insert({ tn,1 });
    }
    //vector<int> tmp_cliques;
    std::unordered_map<int, int> vis_cliques;


    std::unordered_map<int, vector<int>>tmp_cliques;
    std::unordered_map<int, int>in_cliques;
    vector<vector<int>> tmp_cliques_in_node;
    vector<int> clq;
    std::unordered_map<int, int>vis_clique;

    for (auto u : nodes) {
        if (veri_vtx[u] != num_verify) {
            q.push(u);
            veri_vtx[u] = num_verify;
        }
        while (!q.empty()) {
            int v = q.front(); q.pop();
            for (auto hc : v2cliques[v]) {
             
                clq = h_cliques[hc];
                bool Vaild_cq = true;

                if (vis_clique.find(hc) != vis_clique.end()) Vaild_cq = false;
                if (Vaild_cq)
                {
                    for (int i = 0; i < h; i++) {
                        int w = clq[i];
                        if (rho_gu[w] + DELTA < g) {
                            Vaild_cq = false;
                            break;
                        }
                    }
                }
                if (vis_clique.find(hc) == vis_clique.end())
                {
                    vis_clique.insert({ hc,1 });
                }
                if (Vaild_cq == true) {
                    vector<int> in_node = {v};
                    int cnt = 1;
                    for (int i = 0; i < h; i++) {
                        int w = clq[i];
                        if (w == v) continue;
                        if (veri_vtx[w] != num_verify && lhcds_num[w] != -1 && lhcds_num[w] < lhcds_num[u]) {  
                            flag = false;
                        }

                        if (rho_l[w] <= g) {
                            if (veri_vtx[w] != num_verify) {
                                veri_vtx[w] = num_verify;
                                q.push(w);
                            }
                                
                            if (vis_tmp_nodes.find(w) == vis_tmp_nodes.end()) {
                                //in_node.push_back(w);
                                vis_tmp_nodes.insert({w, 1 });
                                tmp_nodes.push_back(w);   
                            }
                            in_node.push_back(w);
                            cnt++;
                        }
                        else {
                            /*flag = false;*/
                        }
                      
                    }

                    if (cnt == h && in_cliques.find(hc) == in_cliques.end()) {
                        in_cliques.insert({ hc,1 });
                    }
                    else if (cnt != h && tmp_cliques.find(hc) == tmp_cliques.end()) {
                        tmp_cliques.insert({ hc, in_node });
                        flag = false;
                    }

                
                }
            }
            for (auto w : adj[v]) {
                if (rho_gu[w] - DELTA >= g) {
                    if (veri_vtx[w] != num_verify) {
                        if (lhcds_num[w] != -1 && lhcds_num[w] < lhcds_num[u]) {
                            flag = false;
                        }
                        if (rho_l[w] + DELTA <= g) {
                            veri_vtx[w] = num_verify;
                            q.push(w);
                            if (vis_tmp_nodes.find(w) == vis_tmp_nodes.end()) {
                                vis_tmp_nodes.insert({ w, 1 });
                                tmp_nodes.push_back(w);
                            }
                        }
                        else {
                            flag = false;
                        }
                    }
                }
            }
        }
    }

    if (flag) return flag;

    Network network;
    vector<int> id_in_network(n, -1);
    vector<int> idx_in_slt_nodes(nodes.size(), 0);
    unsigned long long m = 0;
    vector<Network::Vertex> R;

    double g_temp = g - 1.0 / tmp_nodes.size() / tmp_nodes.size();

    //Add a source node s and a sink node 
    Network::Vertex s = network.AddVertex(), t = network.AddVertex();

    for (unsigned i = 0; i < tmp_nodes.size(); ++i) {
        id_in_network[tmp_nodes[i]] = i;
        R.push_back(network.AddVertex());
    }
    std::unordered_map<int, double> deg_new;

    for (auto _hc : in_cliques) {
        Network::Vertex v = network.AddVertex();
        for (auto hc_node : h_cliques[_hc.first]) {
            if (deg_new.find(hc_node) == deg_new.end())
            {
                deg_new.insert({ hc_node,1.0 });
            }
            else deg_new[hc_node] += 1.0;
            network.AddEdge(R[id_in_network[hc_node]], v, 1);
            network.AddEdge(v, R[id_in_network[hc_node]], h - 1);
        }
    }
    std::unordered_map<int, vector<int>>::iterator iter;
    for (iter = tmp_cliques.begin(); iter != tmp_cliques.end(); iter++) {
        Network::Vertex v = network.AddVertex();
        vector<int> in_nodes = iter->second;
        for (auto hc_node : in_nodes) {
            if (deg_new.find(hc_node) == deg_new.end())
            {
                deg_new.insert({ hc_node,1.0 + (h - in_nodes.size()) * 1.0 / in_nodes.size() });
            }
            else deg_new[hc_node] += 1.0 + (h - in_nodes.size()) * 1.0 / in_nodes.size();
            network.AddEdge(R[id_in_network[hc_node]], v, 1.0 + (h - in_nodes.size()) * 1.0 / in_nodes.size());
            network.AddEdge(v, R[id_in_network[hc_node]], h - 1);
        }
    }
    for (unsigned i = 0; i < tmp_nodes.size(); ++i) {
        network.AddEdge(s, R[i], deg_new[tmp_nodes[i]]);  //
        network.AddEdge(R[i], t, g_temp * h);
    }

    network.MaxFlow(s, t);
    std::unordered_map<int, int> S;
    std::vector<bool> visited(tmp_nodes.size() + in_cliques.size() + tmp_cliques.size() + 2, false);
    network.dfs(0, visited);
    for (int i = 0; i < visited.size(); ++i) {
        if (visited[i]) {
            S.insert({ i,1 });
        }
    }
    std::unordered_map<int, int>mp;

    for (int i = 0; i < nodes.size(); i++)
    {
        mp.insert({ nodes[i],1 });
    }
    for (int i = 0; i < nodes.size(); i++)
    {
        int u = nodes[i];
        for (auto w : adj[u])
        {
            if (u == 1144 && w == 2693) {
                cout << "here" << endl;
            }
            if (mp.find(w) == mp.end() && id_in_network[w] != -1 && S.find(R[id_in_network[w]]) != S.end())
            {
                
                cout << u << " " << w;
                return false;
            }
        }
    }

    return true;

}

bool Graph::verify_LhCDS(vector<int>& nodes, double g) {   //
    for (auto u : nodes) {
        for (auto v : adj[u]) {
            if (rho_l[v] > g) {
                printf("validate failed\n");
                return false;
            }
        }
    }
    for (auto u : nodes) {
        lhcds_num[u] = lhcdses.size();
    }
    bool flag = true;


    ++num_verify;
    queue<int> q;
    vector<pair<int, int>> tmp_edges;
    for (auto u : nodes) {
        if (veri_vtx[u] != num_verify) { 
            q.push(u);
            veri_vtx[u] = num_verify;
        }
        while (!q.empty()) {
            int v = q.front(); q.pop();
            for (auto w : adj[v]) {
                if (fn_baseline) {
                    if (rho_gu[w] >= g) {
                        if (veri_vtx[w] != num_verify) {
                            if (lhcds_num[w] != -1 && lhcds_num[w] < lhcds_num[u]) {
                                flag = false;
                            }
                            veri_vtx[w] = num_verify;
                            q.push(w);
                        }
                        if (v < w)
                            tmp_edges.emplace_back(v, w);
                    }
                }
                else {
                    if (rho_gu[w] >= g) {
                        if (veri_vtx[w] != num_verify) { // w not belong to T
                            if (lhcds_num[w] != -1 && lhcds_num[w] < lhcds_num[u]) {  
                                flag = false;
                            }
                            //                        veri_vtx[w] = num_verify;
                            //                        q.push(w);
                            if (rho_l[w] <= g) {
                                veri_vtx[w] = num_verify;
                                q.push(w);
                            }
                            else {
                                flag = false;
                                tmp_edges.emplace_back(v, v);
                            }
                        }
                        if (v < w && rho_l[w] <= g) 
                            tmp_edges.emplace_back(v, w);
                    }
                }
            }
        }
    }
    printf("size of tmp edges %lu\n", tmp_edges.size());
    if (flag) return flag;

    flag = true;
    FlowNetwork fn = FlowNetwork(tmp_edges, g - 1.0 / n / n, true);
    vector<int> tmp_nodes;
    fn.get_mincut(0, fn.n - 1, tmp_nodes); //tmp_nodes = S
    ++num_verify;
    printf("size of tmp nodes %lu\n", tmp_nodes.size());
    for (auto u : tmp_nodes) {
        //        printf("%d ", u);
        veri_vtx[u] = num_verify;
    }
    printf("\n");
    for (auto u : nodes) {
        for (auto v : adj[u]) {
            if (lhcds_num[v] != lhcds_num[u] && veri_vtx[v] == num_verify) {
                flag = false;
                break;
            }
        }
        if (!flag) break;
    }

    if (flag) return flag;
    printf("validate failed\n");

    for (auto u : nodes) {
        lhcds_num[u] = -1;
    }
    return false;
}

void Graph::output(char* ds_address) {
    printf("num of LhCDS %lu\n", lhcdses.size());
    FILE* dsFile = fopen(ds_address, "w");
    for (int i = 0; i < lhcdses.size(); i++) {

        fprintf(dsFile, "%lu %.4f\n", lhcdses[i].size(), lhcds_rho[i]);
        for (auto& u : lhcdses[i]) {
            fprintf(dsFile, "%d ", u);
        }
        fprintf(dsFile, "\n");
    }

    fclose(dsFile);
}

int Graph::find_fa(int x) {
    if (x != fa[x])
        fa[x] = find_fa(fa[x]);
    return fa[x];
}

void Graph::connected_components() {
    for (auto u : slt_nodes) {
        fa[u] = u;
    }

    for (auto e : slt_edges) {
        //        if (find_fa(edges[e].first) != find_fa(edges[e].second))
        int x = find_fa(edges[e].first);
        int y = find_fa(edges[e].second);
        //        printf("%d %d %d %d\n", edges[e].first, edges[e].second, x, y);
        fa[x] = y;
    }
    for (auto u : slt_nodes) {
        find_fa(u);
    }

    sort(slt_nodes.begin(), slt_nodes.end(), [this](int a, int b)->bool {
        return fa[a] > fa[b];
        });
    sort(slt_edges.begin(), slt_edges.end(), [this](int a, int b)->bool {
        return fa[edges[a].first] > fa[edges[b].first];
        });
    sort(slt_h_cliques.begin(), slt_h_cliques.end(), [this](int a, int b)->bool {
        return fa[h_cliques[a][0]] > fa[h_cliques[b][0]];
        });

    //    printf("fa ");
    //    for (int i = 0; i < n; i++) {
    //        printf("%d ", find_fa(i));
    //    }
    //    printf("\n");

    cmpt.clear();
    int x = 1, y = 1, z = 1;
    while (x < slt_nodes.size() && y < slt_edges.size()) {
        while (x < slt_nodes.size() && fa[slt_nodes[x]] == fa[slt_nodes[x - 1]]) x++;
        while (y < slt_edges.size() && fa[edges[slt_edges[y]].first] == fa[edges[slt_edges[y - 1]].first]) y++;
        while (z < slt_h_cliques.size() && fa[h_cliques[slt_h_cliques[z]][0]] == fa[h_cliques[slt_h_cliques[z - 1]][0]]) z++;
        cmpt.emplace_back(x, y, z);
        ++x;
        ++y;
        ++z;
    }
    //    for (auto pr : cmpt) {
    //        printf("%d %d\n", pr.first, pr.second);
    //    }
}
