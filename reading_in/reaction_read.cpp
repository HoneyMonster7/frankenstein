

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/visitors.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <cstdlib>
#include<string>
#include<iostream>
#include<fstream>


class reaction{
	std::string type;
	int substrateIndex;
        int productIndex;
	int nrATP;
	int nrNADH;
	double freeEChange;
	std::string humanReadable;	
	public:

	static void readReactions (std::string fileName, std::vector<reaction> &reacPointer);
	reaction(std::string tmpType,int tmpsubI, int tmpProdI, int tmpnrATP, int tmpnrNADH, double tmpfreeE, std::string tmpHumRead);

};


int main (int argc, char* argv[]){
using namespace boost;


std::cout<<"Tests begin."<<std::endl;

std::vector<reaction> reacVector;


	 std::cout<<"length of the vector is: "<<reacVector.size()<<std::endl;
	 reaction::readReactions("reactions__4C_v3_2_2_ext_100.dat", reacVector);
	 std::cout<<"length of the vector is: "<<reacVector.size()<<std::endl;
std::cout<<"Tests completed."<<std::endl;
}



void reaction::readReactions (std::string fileName, std::vector<reaction> &reacPointer){
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

