#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/push_relabel_max_flow.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/stoer_wagner_min_cut.hpp>
#include "MaxFlow.h"
using namespace std;
using namespace boost;

Network::Network() : g(), rev(get(edge_reverse, g)) {}

Network::Vertex Network::AddVertex() {
	return add_vertex(g);
}

void Network::AddEdge(Vertex& v1, Vertex& v2, const long capacity) {
	Traits::edge_descriptor e1 = add_edge(v1, v2, g).first;
	Traits::edge_descriptor e2 = add_edge(v2, v1, g).first;
	put(edge_capacity, g, e1, capacity);
	rev[e1] = e2;
	rev[e2] = e1;
}

unsigned long long Network::MaxFlow(Vertex& s, Vertex& t) {
	return push_relabel_max_flow(g, s, t);
}

void Network::dfs(int v, std::vector<bool>& visited) {
	visited[v] = true;
	Graph::out_edge_iterator ei, e_end;
	for (boost::tie(ei, e_end) = out_edges(v, g); ei != e_end; ++ei) {
		//cout << target(*ei, g);
		if (!visited[target(*ei, g)] && get(boost::edge_residual_capacity, g, *ei) > 0) {
			dfs(target(*ei, g), visited);
		}
	}
}