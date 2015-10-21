

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/visitors.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <cstdlib>
#include <boost/graph/bipartite.hpp>

#include<string>
#include<iostream>
#include<fstream>
#include<utility>
#include<algorithm>

struct environment{

	double temperature=293; //temperature in kelvins
	double atpCont=1;
	double adpCont=1;
	double ampCont=1;
	double nadplusCont=1;
	double nadhCont=1;
	double co2Cont=1;
	double piCont=1;
	double ppiCont=1;
	double h2oCont=1;
	double nh3Cont=1;
	double glutCont=1;
	double oxo2Cont=1;
};



struct substrate{
	//properties of each substrate
	int index;
	double freeOfCreation;
	std::string molecule;
	std::string name;
	int charge;
};



	struct SubReac; 	

	 typedef boost::adjacency_list<boost::listS,boost::vecS,boost::directedS,SubReac,bool> ReactionNetwork;
	 typedef boost::graph_traits<ReactionNetwork>::vertex_descriptor Vertex;


class reaction{
	//properties of each reaction, should include an int (possible shorter) for each internal metabolite to aid the recalculation of freeEchange for different environments
//	std::string type;
	std::vector<int> substrates;
	std::vector<int> products;
	int [9] internalMets;
	double freeEChange;
	//std::string humanReadable;	
	double currentFreeEChange=0;
	public:

	static void readReactions (std::string fileName, std::vector<reaction> &reacPointer, ReactionNetwork &lofasz,  std::vector<Vertex> &vertexList, const std::vector<Vertex> &compoundVList);

	static void readCompounds (std::string fileName, ReactionNetwork &lofasz, std::vector<Vertex> &compoundlist);

	reaction(std::vector<int> tmpsubstrates, std::vector<int> tmproducts, int [] internalMets ); 
	reaction();
	void printReaction ();
	double freeEchange();
	int getsubstrateIndex();
	int getproductIndex();
	int getnrATP();
	int getnrNADH();
	std::string gettype();
	void recalcEchange(environment env);
//	~reaction();
};

struct SubReac{
substrate sub;
reaction reac;
};


	 //template <typename Graph >

int main (int argc, char* argv[]){
using namespace boost;


std::cout<<"Tests begin."<<std::endl;

std::vector<reaction> reacVector;



	 std::vector<Vertex> reacVList,compoundVList;
	 ReactionNetwork lofasz;


         //reaction::readCompounds("compounds_list__4C_v3_2_2_ext_100.dat",lofasz,compoundVList);
         reaction::readCompounds("newsortedcompounds.txt",lofasz,compoundVList);

	 std::cout<<"length of the vector is: "<<reacVector.size()<<std::endl;
	 reaction::readReactions("reactions__4C_v3_2_2_ext_100.dat", reacVector,lofasz,reacVList,compoundVList);
	 std::cout<<"length of the vector is: "<<reacVector.size()<<std::endl;
	 reacVector[13].printReaction();
	 reacVector[299].printReaction();
	
	 

	 graph_traits<ReactionNetwork>::vertex_iterator vi, vi_end;
	 int count=0;
	 for (boost::tie(vi, vi_end)=vertices(lofasz); vi!=vi_end; ++vi){
	 count++;
	 }

	 std::cout<<"The graph has "<<count<<" vertices."<<std::endl;


	 typedef graph_traits <ReactionNetwork> traits;
	   typename traits::vertex_iterator vertex_iter, vertex_end;
		 
	   bool bipartiteee = boost::is_bipartite(lofasz);

	 if (bipartiteee) {std::cout<<"The graph is bipartite"<<std::endl;}
	 else {std::cout<<"The graph is not bipartite."<<std::endl;}

	 reaction newreaction(" ",0,0,0,0,0," ");

	 reaction evennewer;
	 evennewer.printReaction();
	 evennewer=newreaction;

	 newreaction.printReaction();
	 evennewer.printReaction();
	 std::cin.ignore();



std::cout<<"Tests completed."<<std::endl;
}

void reaction::readCompounds (std::string fileName, ReactionNetwork  &graph,std::vector<Vertex> &vertexList){
	  
	  
	  
	  
	          //lofasz[vertexList[vertexList.size()-1]].reac=newreaction;
		   
		  std::ifstream inFile(fileName);
		  std::string line,nameChem,namenorm;
		  int index, charge; 
		  double formFreeE;
		  while (std::getline(inFile,line)){
		  	std::stringstream iss(line);
		  	iss>>index;
		  	iss>>formFreeE;
		  	iss>>nameChem;
		  	iss>>namenorm;
		  	iss>>charge;
		 vertexList.emplace_back(boost::add_vertex(graph));
			}
}
void reaction::readReactions (std::string fileName, std::vector<reaction> &reacPointer,ReactionNetwork  &graph,std::vector<Vertex> &vertexList,const std::vector<Vertex> &compoundVList){




	//lofasz[vertexList[vertexList.size()-1]].reac=newreaction;

	std::ifstream inFile(fileName);
	std::string line,tmpsubs,tmpprods;

	std::vector<int> tmpsubstrates, tmproducts;
       	int tmpnrATP=0, tmpnrPPi=0, tmpnrPi=0 ,  tmpnrADP=0,tmpnrNAD_red=0,tmpnrNAD_ox=0, tmpnrCO2=0,tmpnrH2O=0;
	int [9] tmpinternalMets= {};

	double tmpfreeE;

	typedef boost::graph_traits<ReactionNetwork>::edge_descriptor Edge;

	while (std::getline(inFile,line)){
		std::stringstream iss(line);
		iss>>tmpfreeE;
		std::getline(iss,tmpsubs, '>');
		std::getline(iss,tmpprods);
		
		std::istringstream subss(tmpsubs);
		std::istringstream prodss(tmpprods);	

		for (int tmp; subss>>tmp;){

			tmpsubstrates.push_back(tmp);
		}

		for (int tmp; prodss>>tmp;){
			tmproducts.push_back(tmp);
		}

		//finding if any of the internal metabolites appear on any side of the reaction
		for (int interMIndex=-9; i<0; i++){
			if(std::find(tmpsubstrates.begin(), tmpsubstrates.end(),i)){
				tmpinternalMets[i+9]--;
			}
	
			if(std::find(tmproducts.begin(), tmproducts.end(),i)){
				tmpinternalMets[i+9]++;
			}
		}


	 	reacPointer.emplace_back(tmpsubstrates,tmproducts,tmpinternalMets);	


		vertexList.emplace_back(boost::add_vertex(graph));
		graph[vertexList[vertexList.size()-1]].reac=reaction(tmpsubstrates,tmproducts,tmpinternalMets);


		Edge e1;
		//not adding edges yet
		//e1=(boost::add_edge(compoundVList[tmpsubI],vertexList[vertexList.size()-1],graph)).first;




//		tmpSR.reac=new reaction(tmpType,tmpsubI,tmpProdI,tmpnrATP,tmpnrNADH,tmpfreeE,tmpHumRead);


		//Possibly also create the bipartite graph at this stage
		//graph structure: substrate number n is vertex number 2n
		//		   reaction number n is vertex number 2n+1
		//vertex properties are a pair <substrate,reaction> with the
		//non-important one being blank

		//std::cout<<"Type: "<<tmpType<<std::endl;
		//std::cout<<"Numbers: "<<tmpsubI<<" "<<tmpProdI<<" "<<tmpnrATP<<" "<<tmpnrNADH<<" "<<tmpfreeE<<std::endl;
		//std::cout<<"Humanreadable: "<<tmpHumRead<<std::endl;

	}

}

reaction::reaction(std::vector<int> tmpsubstrates, std::vector<int> tmproducts, int [9] tmpInternalMets )  {

	substrates=tmpsubstrates;
	products=tmproducts;

	internalMets=tmpinternalMets;



}
void reaction::printReaction(){
	std::cout<<"From: ";
	for(auto i =substrates.begin(); i!=substrates.end(); ++i) std::cout<<*i<<' ';
	std::cout<<std::endl;
	std::cout<<"To: ";
	for(auto i=products.begin(); i!=products.end(); ++i) std::cout<<*i<<' ';
	std::cout<<std::endl;
	std::cout<<"Free energy change: "<<freeEChange<<std::endl;

}

double reaction::freeEchange(){ return freeEChange;}
int reaction::getsubstrateIndex(){return substrateIndex;}
int reaction::getproductIndex(){return productIndex;}
int reaction::getnrATP(){return nrATP;}
int reaction::getnrNADH(){return nrNADH;}
std::string reaction::gettype(){return type;}

void recalcEchange(environment env){
	//recalculate the freeEchang using the formula G=G0+RT*ln(([C]^c*[D]^d)/([A]^a*[B]^b))
	//need the reactions database changed for this

}
reaction::reaction(){
	substrates=new std::vector<int>;
	products=new std::vector<int>;

	int [9] tmpint={};
	internalMets=tmpint;
}

/*reaction::~reaction(){
	delete &type;
	delete &substrateIndex;
	delete &productIndex;
	delete &nrATP;
	delete &nrNADH;
	delete &freeEChange;
	delete &humanReadable;



}*/

