#ifndef _MAX_FLOW_HPP_
#define _MAX_FLOW_HPP_
#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_utility.hpp>

using namespace std;
using namespace boost;

class Network {
private:
    typedef adjacency_list_traits < vecS, vecS, directedS > Traits;

public:
    typedef Traits::vertex_descriptor Vertex;
    Network();
    Vertex AddVertex();
    void AddEdge(Vertex& v1, Vertex& v2, const long capacity);
    void dfs(int v, std::vector<bool>& visited);
    unsigned long long MaxFlow(Vertex& s, Vertex& t);

private:
    typedef adjacency_list < vecS, vecS, directedS,
        property < vertex_name_t, std::string,
        property < vertex_index_t, long,
        property < vertex_color_t, boost::default_color_type,
        property < vertex_distance_t, long,
        property < vertex_predecessor_t, Traits::edge_descriptor > > > > >,

        property < edge_capacity_t, long long,
        property < edge_residual_capacity_t, long long,
        property < edge_reverse_t, Traits::edge_descriptor > > > > Graph;
    Graph g;
    property_map < Graph, edge_reverse_t >::type rev;
};

#endif