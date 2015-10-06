

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/visitors.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <cstdlib>
#include<string>
#include<iostream>
#include<fstream>


struct reaction{
	std::string type;
	int substrateIndex;
        int productIndex;
	int nrATP;
	int nrNADH;
	double freeEChange;
	std::string humanReadable;	
};

void readReactions (std::string fileName); 

int main (int argc, char* argv[]){
using namespace boost;


std::cout<<"Tests begin."<<std::endl;
readReactions("reactions__4C_v3_2_2_ext_100.dat");
std::cout<<"Tests completed."<<std::endl;
}

void readReactions (std::string fileName){
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



		std::cout<<"Type: "<<tmpType<<std::endl;
		std::cout<<"Numbers: "<<tmpsubI<<" "<<tmpProdI<<" "<<tmpnrATP<<" "<<tmpnrNADH<<" "<<tmpfreeE<<std::endl;
		std::cout<<"Humanreadable: "<<tmpHumRead<<std::endl;

	}

}
