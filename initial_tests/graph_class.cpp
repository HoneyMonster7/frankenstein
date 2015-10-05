
#include<set>
#include<iostream>
#include<vector>
using namespace std;



class graph {
   	static const int msize = 150;
	bool adjmat [msize][msize] ;
	
public: 
   int edgeN=0;

   const int matsize=msize;
   void addEdge (int,int);
   void printStruct ();
   void connectedTo(int);
   void depthS(int,set<int> *visited );
   
} testgraph;


void graph::addEdge (int x, int y) {

 if (( x == y) || (x>msize) || (y>msize)) { std::cout<<"Can't connect those."<<endl;}
 else {
	adjmat[x][y]=true;
	adjmat[y][x]=true; 
//	std::cout<<"Edge ("<<x<<","<<y<<") added."<<endl; 
	edgeN++;}

}

void graph::remEdge (int x, int y) {
	 if (( x == y) || (x>msize) || (y>msize)) { std::cout<<"Can't connect those."<<endl;}
	 else {adjmat[x][y]=false;
	       adjmat[y][x]=false;
	       edgeN--;
}


void graph::printStruct() {
cout<<"Printing the graph with "<<edgeN<<" edges."<<endl;
	for (auto& i  : adjmat) {
		for (auto& edge : i) {
			if (edge) {std::cout<<"1";}
			else {std::cout<<"0";}
	}
		std::cout<<endl;
}
}

void graph::connectedTo(int origV){
	set<int> visited;
        visited.insert(origV);
	cout<<"Original inserted."<<endl;
	testgraph.depthS(origV,&visited);
 	cout<<"The edges connected to "<<origV<<" are:";
	std::set<int>::iterator it;
	for (it=visited.begin(); it !=visited.end(); it++){
		int curr=*it;
		cout<<curr<<",";
}

}

void graph::depthS( int calledon, set<int> *visited){
	cout<<"DS called on:"<<calledon<<"... ";
	int vcount=0;
//	vector<int> vertexVis;
	for (int i=0; i<msize; i++){
		if ((adjmat[calledon][i]) && ((*visited).find(i)==(*visited).end())){
			(*visited).insert(i);
		/*	vertexVis[vcount]=i; vcount++;*/
			testgraph.depthS(i,visited);
			
			
	cout<<"DS on "<<calledon<<" finished."<<endl;
			
		}
	}
}

int main () {

cout<<"Starting..."<<endl;
srand( time( NULL ) );


for (int i=0; i<testgraph.matsize; i++){

	for (int j=0; j<i; j++){
	/*	cout<<"("<<i<<","<<j<<")  "<<rand() % 2<<endl; */
        if ((rand() % 100) < 1) {testgraph.addEdge(i,j);}
	}
}

testgraph.printStruct();

testgraph.connectedTo(3);

std::cout<<"Hopefully all went well."<<endl;

}
