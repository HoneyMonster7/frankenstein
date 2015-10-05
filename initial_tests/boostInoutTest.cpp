

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/visitors.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <cstdlib>
#include<string>


int main (int argc, char* argv[]){
using namespace boost;
const int matsize=100;

std::cout<<"Prob(out of 10000): ";
int prob;
//std::cin>>prob;


adjacency_list<vecS,vecS,undirectedS> g;

//bfs_visitor<my_discover_visitor> vis;

srand( time(NULL) );
int edgeN=0;
int nrTries=1e9;


/*for (int i=0; i<100; i++){
	adjacency_list<>::vertex_descriptor v = add_vertex(g)
}*/

for (int i=0; i<nrTries; i++){
        int x=rand()%100;
        int y=rand()%100;
        if (x!=y){ add_edge(x,y,g);
                   remove_edge(x,y,g);}

}




std::cout<<"Edges added, trial ended."<<std::endl;


}
