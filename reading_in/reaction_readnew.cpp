

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

	//substrates from newsortedcompounds.txt in this order
	//in the internalMets array they are kept in this order

	double ppiCont=1;
	double piCont=1;
	double atpCont=1;
	double adpCont=1;
	double ampCont=1;
	double nadredcont=1;
	double nadoxcont=1;
	double co2cont=1;
	double h2ocont=1;


	double temperature=293; //temperature in kelvins
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
	int  internalMets [9];
	double freeEChange;
	//std::string humanReadable;	
	double currentFreeEChange=0;
	public:

	static void readReactions (std::string fileName, std::vector<reaction> &reacPointer, ReactionNetwork &lofasz,  std::vector<Vertex> &vertexList, const std::vector<Vertex> &compoundVList);

	static void readCompounds (std::string fileName, ReactionNetwork &lofasz, std::vector<Vertex> &compoundlist);

	reaction(double tmpfreechange, std::vector<int> tmpsubstrates, std::vector<int> tmproducts,int (&internalMets)[9] ); 
	reaction();
	void printReaction ();
	double freeEchange();
	std::vector<int> getsubstrates();
	std::vector<int> getproducts();
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


         //reaction::readCompounds(DATA_PATH "compounds_list__4C_v3_2_2_ext_100.dat",lofasz,compoundVList);
         reaction::readCompounds(DATA_PATH "newsortedcompounds.txt",lofasz,compoundVList);

	 std::cout<<"length of the vector is: "<<reacVector.size()<<std::endl;
	 //reaction::readReactions(DATA_PATH "reactions__4C_v3_2_2_ext_100.dat", reacVector,lofasz,reacVList,compoundVList);
	 reaction::readReactions(DATA_PATH "newreactions.txt", reacVector,lofasz,reacVList,compoundVList);
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

	 reaction newreaction;

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

       	//int tmpnrATP=0, tmpnrPPi=0, tmpnrPi=0 ,  tmpnrADP=0,tmpnrNAD_red=0,tmpnrNAD_ox=0, tmpnrCO2=0,tmpnrH2O=0;
	int tmpinternalMets [9] = {};

	double tmpfreeE;

	typedef boost::graph_traits<ReactionNetwork>::edge_descriptor Edge;

	while (std::getline(inFile,line)){
		
		std::vector<int> tmpsubstrates, tmproducts;
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
		for (int i=-9; i<0; i++){
			if(std::find(tmpsubstrates.begin(), tmpsubstrates.end(),i)!=tmpsubstrates.end()){
				tmpinternalMets[i+9]--;
			}
	
			if(std::find(tmproducts.begin(), tmproducts.end(),i)!=tmproducts.end()){
				tmpinternalMets[i+9]++;
			}
		}


	 	reacPointer.emplace_back(tmpfreeE,tmpsubstrates,tmproducts,tmpinternalMets);	


		vertexList.emplace_back(boost::add_vertex(graph));
		graph[vertexList[vertexList.size()-1]].reac=reaction(tmpfreeE,tmpsubstrates,tmproducts,tmpinternalMets);


		for (int i : tmpsubstrates){
		Edge e1;
		e1=(boost::add_edge(compoundVList[i+9],vertexList[vertexList.size()-1],graph)).first;

		}


		for (int i: tmproducts){

			Edge e1;
			e1=(boost::add_edge(vertexList[vertexList.size()-1],compoundVList[i+9],graph)).first;

				}

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

reaction::reaction(double tmpfreechange, std::vector<int> tmpsubstrates, std::vector<int> tmproducts, int (&tmpInternalMets) [9] )  {

	freeEChange=tmpfreechange;
	substrates=tmpsubstrates;
	products=tmproducts;

	std::copy(std::begin(tmpInternalMets),std::end(tmpInternalMets),std::begin(internalMets));



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
std::vector<int> reaction::getsubstrates() {return substrates;}
std::vector<int> reaction::getproducts() {return products ;}


void reaction::recalcEchange(environment env){
	//recalculate the freeEchang using the formula G=G0+RT*ln(([C]^c*[D]^d)/([A]^a*[B]^b))
	//need the reactions database changed for this

	double insideLog=(std::pow(env.ppiCont,internalMets[0])*std::pow(env.piCont,internalMets[1])*std::pow(env.atpCont,internalMets[2])*std::pow(env.adpCont,internalMets[3])*std::pow(env.ampCont,internalMets[4])*std::pow(env.nadredcont,internalMets[5])*std::pow(env.nadoxcont,internalMets[6])*std::pow(env.co2cont,internalMets[7])*std::pow(env.h2ocont,internalMets[8]));

	currentFreeEChange=freeEChange+8.3144598*env.temperature*std::log(insideLog);
}
reaction::reaction(){
	std::vector<int> substrates;
	std::vector<int> products;

	int tmpint [9] = {};
	std::copy(std::begin(tmpint),std::end(tmpint),std::begin(internalMets));
	
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

