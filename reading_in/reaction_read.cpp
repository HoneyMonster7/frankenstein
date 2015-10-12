

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/visitors.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <cstdlib>
#include<string>
#include<iostream>
#include<fstream>
#include<utility>


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

class reaction{
	//properties of each reaction, should include an int (possible shorter) for each internal metabolite to aid the recalculation of freeEchange for different environments
	std::string type;
	int substrateIndex;
        int productIndex;
	int nrATP;
	int nrNADH;
	double freeEChange;
	std::string humanReadable;	
	double currentFreeEChange=0;
	public:

	static void readReactions (std::string fileName, std::vector<reaction> &reacPointer);
	reaction(std::string tmpType,int tmpsubI, int tmpProdI, int tmpnrATP, int tmpnrNADH, double tmpfreeE, std::string tmpHumRead);
	reaction();
	void printReaction ();
	double freeEchange();
	int getsubstrateIndex();
	int getproductIndex();
	int getnrATP();
	int getnrNADH();
	std::string gettype();
	void recalcEchange(environment env);
};

struct SubReac{
substrate sub;
reaction reac;
};


int main (int argc, char* argv[]){
using namespace boost;


std::cout<<"Tests begin."<<std::endl;

std::vector<reaction> reacVector;


	 std::cout<<"length of the vector is: "<<reacVector.size()<<std::endl;
	 reaction::readReactions("reactions__4C_v3_2_2_ext_100.dat", reacVector);
	 std::cout<<"length of the vector is: "<<reacVector.size()<<std::endl;
	 reacVector[13].printReaction();
	 reacVector[299].printReaction();
	

	 reaction newreaction(" ",0,0,0,0,0," ");

	 reaction evennewer;

	 evennewer=newreaction;

	 newreaction.printReaction();
	 evennewer.printReaction();

std::cout<<"Tests completed."<<std::endl;
}



void reaction::readReactions (std::string fileName, std::vector<reaction> &reacPointer){

	//defining the boost graph (to be biparite)
	typedef boost::adjacency_list<boost::listS,boost::vecS,boost::directedS,SubReac,bool> ReactionNetwork;
	reaction emptyReac("Empty reaction",0,0,0,0,0,"Empty Humanreadable");
	substrate emptySubstrate;
	emptySubstrate.molecule="EmptyMolecule";
	emptySubstrate.index=0;
	emptySubstrate.charge=0;
	emptySubstrate.name="EmptyMolecule name";
	emptySubstrate.freeOfCreation=0;

	SubReac emptySubReac;
	emptySubReac.reac=emptyReac;
	emptySubReac.sub=emptySubstrate;




	std::ifstream inFile(fileName);
	std::string tmpType, tmpHumRead,line,tmpstring;
	int tmpsubI,tmpProdI,tmpnrATP,tmpnrNADH;
	double tmpfreeE;
	while (std::getline(inFile,line)){
		std::stringstream iss(line);
		iss>>tmpType;
		iss>>tmpsubI;	
	        iss>>tmpProdI;	
	        iss>>tmpnrATP;	
	        iss>>tmpnrNADH;	
	        iss>>tmpfreeE;	
		std::getline(iss, tmpHumRead);

	 	reacPointer.emplace_back(tmpType,tmpsubI,tmpProdI,tmpnrATP,tmpnrNADH,tmpfreeE,tmpHumRead);	

		SubReac tmpSR;
		tmpSR.sub=emptySubstrate;
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

reaction::reaction(std::string tmpType,int tmpsubI, int tmpProdI, int tmpnrATP, int tmpnrNADH, double tmpfreeE, std::string tmpHumRead) {
	type=tmpType;
	substrateIndex=tmpsubI;
	productIndex=tmpProdI;
	nrATP=tmpnrATP;
	nrNADH=tmpnrNADH;
	freeEChange=tmpfreeE;
	humanReadable=tmpHumRead;
}

void reaction::printReaction(){
	std::cout<<"Type: "<<type<<std::endl;
	std::cout<<"Numbers: "<<substrateIndex<<" "<<productIndex<<" "<<nrATP<<" "<<nrNADH<<" "<<freeEChange<<std::endl;
	
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
	type="Empty";
	substrateIndex=0;
	productIndex=0;
	nrATP=0;
	nrNADH=0;
	freeEChange=0;
	humanReadable="EmptyHuman";

}
