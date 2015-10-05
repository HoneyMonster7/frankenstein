

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/visitors.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <cstdlib>
#include<string>
#include<time.h>

int main (int argc, char* argv[]){
using namespace boost;
const int matsize=100;

std::cout<<"Prob(out of 10000): ";
int prob;
std::cin>>prob;

int numberofTries=100;
int success=0;
int edgetot=0;
bool printRoute=false;

for (int cnt=0; cnt<numberofTries; cnt++){


adjacency_list<vecS,vecS,undirectedS> g;

//bfs_visitor<my_discover_visitor> vis;

struct timespec ts;
clock_gettime(CLOCK_MONOTONIC, &ts);

srand((time_t)ts.tv_nsec );
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
bool reachedGoal=true;
while (p !=4){


if(printRoute) {std::cout<<p<<",";}
p=predecessors[p];
if (p==prev){
		if (printRoute){std::cout<<std::endl<<"No route found."<<std::endl;} 
		reachedGoal=false;
		break;}
else {prev=p;}
}
if(printRoute){std::cout<<p<<std::endl;}
if(reachedGoal){success++; edgetot+=edgeN;}


}

std::cout<<"Out of "<<numberofTries<<" trials "<<success<<" were successful, and on average these had "<<edgetot/success<<" edges."<<std::endl;

/*for (auto d : distances ) {
std::cout<<d<<","; }
std::cout<<std::endl;*/



}
