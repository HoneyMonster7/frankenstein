

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
std::cin>>prob;


adjacency_list<vecS,vecS,undirectedS> g;

//bfs_visitor<my_discover_visitor> vis;

srand( time(NULL) );
int edgeN=0;
for (int i=0; i<matsize; i++){
	for (int j=i+1; j<matsize; j++){
	if (std::rand()%10000 <prob){
		add_edge(i,j,g);
		edgeN++;
		/*std::cout<<"Added ("<<i<<","<<j<<")"<<std::endl;*/ }
	}

}

std::cout<<edgeN<<" edges added of the possible "<<(matsize*(matsize-1))/2<<"."<<std::endl;


std::array<int, matsize> distances = {{ 0 }}; 
std::array<int, matsize> predecessors = {{ 0 }}; 
predecessors[4]=4;
/*breadth_first_search(g, 4,
 visitor(
 make_bfs_visitor(
 record_distances(distances.data(),
 on_tree_edge())))); */


breadth_first_search(g,4,
 visitor(
 make_bfs_visitor(
 record_predecessors(predecessors.data(),
 on_tree_edge())))); 


int p=91; //where we want to reach
int prev=91;
while (p !=4){

std::cout<<p<<",";
p=predecessors[p];
if (p==prev){ std::cout<<std::endl<<"No route found."<<std::endl; break;}
else {prev=p;}
}
std::cout<<p<<std::endl;

/*for (auto d : distances ) {
std::cout<<d<<","; }
std::cout<<std::endl;*/



}
