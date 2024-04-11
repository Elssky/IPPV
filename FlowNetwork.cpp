#include "FlowNetwork.h"

FlowNetwork::FlowNetwork(vector<pair<int, int>> edges, double g, bool LhCDSvalidate) {
    n = 0;
    ori_id.push_back(0);
    for (auto & edge : edges) {
        if (mapping.find(edge.first) == mapping.end()) {
            ++n;
            ori_id.push_back(edge.first);
            mapping.insert(make_pair(edge.first, n));
        }
        if (mapping.find(edge.second) == mapping.end()) {
            ++n;
            ori_id.push_back(edge.second);
            mapping.insert(make_pair(edge.second, n));
        }
    }
    n += 2;
    adj.resize(n);
    vector<int> loops(n, 0);
    for (auto &edge : edges) {
        int u = mapping.find(edge.second)->second;
        int v = mapping.find(edge.first)->second;
        if (u != v) {
            add_edge(u, v, 1);
            add_edge(v, u, 1);
        } else {
            ++loops[u];
        }
    }

    if (LhCDSvalidate) {
        for (int u = 1; u < n - 1; u++) {
            add_edge(u, n - 1, edges.size() + 2 * g - adj[u].size() / 2.0 - 2.0 * loops[u]); 
        }
        for (int u = 1; u < n - 1; u++) {
            add_edge(0, u, edges.size()); 
        }
    } else {
        for (int u = 1; u < n - 1; u++) {
            add_edge(0, u, adj[u].size() / 4.0 + loops[u]); // s->C
        }
        for (int u = 1; u < n - 1; u++) {
            add_edge(u, n - 1, g); // u->t
        }
    }
}


FlowNetwork::FlowNetwork(vector<pair<int, int>> edges, vector<pair<int, int>>clique_edges, vector<pair<int, int>>dtr, vector<pair<int, double>>cap, double g) {
    n = 0;
    ori_id.push_back(0);
    for (auto& edge : edges) {
        if (mapping.find(edge.first) == mapping.end()) {
            ++n;
            ori_id.push_back(edge.first);
            mapping.insert(make_pair(edge.first, n));
        }
        if (mapping.find(edge.second) == mapping.end()) {
            ++n;
            ori_id.push_back(edge.second);
            mapping.insert(make_pair(edge.second, n));
        }
    }

    for (auto& edge : clique_edges) {
        if (mapping.find(edge.first) == mapping.end()) {
            ++n;
            ori_id.push_back(edge.first);
            mapping.insert(make_pair(edge.first, n));
        }
        if (mapping.find(edge.second) == mapping.end()) {
            ++n;
            ori_id.push_back(edge.second);
            mapping.insert(make_pair(edge.second, n));
        }
    }
    for (auto& edge : dtr) {
        if (mapping.find(edge.first) == mapping.end()) {
            ++n;
            ori_id.push_back(edge.first);
            mapping.insert(make_pair(edge.first, n));
        }
    }

    n += 2;
    adj.resize(n);
    for (auto& edge : edges) {
        int u = mapping.find(edge.second)->second;
        int v = mapping.find(edge.first)->second;
        add_edge(v, u, 1.0);
    }
    for (auto& edge : clique_edges) {
        int u = mapping.find(edge.second)->second;
        int v = mapping.find(edge.first)->second;
        add_edge(v, u, 2.0);
    }

    for (auto& edge : dtr) {
        int dtrv = edge.second;
        int v = mapping.find(edge.first)->second;
        add_edge(0, v, (double)dtrv);
    }

    for (auto& edge : cap) {
        double checks = edge.second;
        int v = mapping.find(edge.first)->second;
        add_edge(v, n-1, checks);
    }
}

void FlowNetwork::add_edge(int from, int to, double cap) {
    adj[from].push_back(EdgeFN(from, to, cap, 0, adj[to].size()));
    if (from == to) {
        adj[from].back().index++;
    }
    adj[to].push_back(EdgeFN(to, from, 0, 0, adj[from].size() - 1));
}

void FlowNetwork::enqueue (int v) {
    if (!active[v] && excess[v] > 0 && dist[v] < n) {
        active[v] = true;
        B[dist[v]].push_back(v);
        b = std::max(b, dist[v]);
    }
}

void FlowNetwork::push (EdgeFN &e) {
    double amt = std::min(excess[e.from], e.cap - e.flow);
    if (dist[e.from] == dist[e.to] + 1 && amt > 0) {
        e.flow += amt;
        adj[e.to][e.index].flow -= amt;
        excess[e.to] += amt;
        excess[e.from] -= amt;
        enqueue(e.to);
    }
}

void FlowNetwork::gap (int k) {
    for (int v = 0; v < n; v++) {
        if (dist[v] >= k) {
            count[dist[v]]--;
            dist[v] = std::max(dist[v], n);
            count[dist[v]]++;
            enqueue(v);
        }
    }
}

void FlowNetwork::relabel(int v) {
    count[dist[v]]--;
    dist[v] = n;
    for (auto e : adj[v]) if (e.cap - e.flow > 0) {
            dist[v] = std::min(dist[v], dist[e.to] + 1);
        }
    count[dist[v]]++;
    enqueue(v);
}

void FlowNetwork::discharge(int v) {
    for (auto &e : adj[v]) {
        if (excess[v] > 0) {
            push(e);
        } else {
            break;
        }
    }

    if (excess[v] > 0) {
        if (count[dist[v]] == 1) {
            gap(dist[v]);
        } else {
            relabel(v);
        }
    }
}

double FlowNetwork::get_maxflow(int s, int t, bool need_initial) {
    if (need_initial) {
        dist = std::vector<int>(n, 0);
        excess = std::vector<double>(n, 0);
        count = std::vector<int>(n + 1, 0);
        active = std::vector<bool>(n, false);
        B = std::vector<std::vector<int>>(n);
        b = 0;

        for (auto &e: adj[s]) {
            excess[s] += e.cap;
        }

        count[0] = n;
        enqueue(s);
        active[t] = true;
    }

    while (b >= 0) {
        if (!B[b].empty()) {
            int v = B[b].back();
            B[b].pop_back();
            active[v] = false;
            discharge(v);
        } else {
            b--;
        }
    }

    return excess[t];
}

double
FlowNetwork::get_mincut(int s, int t, std::vector<int> &S, bool need_initial) {
    auto ret = get_maxflow(s, t, need_initial);
    S.clear();

    for (int v = 0; v < n; v++) {
//        printf("dist [%d] %d\n", v, dist[v]);
        if (dist[v] >= n) {
            if (v == s)
                continue;
            if (v > 0 && v < n) {
                S.push_back(ori_id[v]);
            }
        }
    }

    return ret;
}