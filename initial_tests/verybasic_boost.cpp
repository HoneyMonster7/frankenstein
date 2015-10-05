

#include <boost/graph/adjacency_list.hpp>

int main (){

using namespace boost;
adjacency_list<> g;
// adds four vertices to the graph
adjacency_list<>::vertex_descriptor v1 = add_vertex(g);
adjacency_list<>::vertex_descriptor v2 = add_vertex(g);
adjacency_list<>::vertex_descriptor v3 = add_vertex(g);
adjacency_list<>::vertex_descriptor v4 = add_vertex(g); 

bfs_visitor<null_visitor> vis;

}
