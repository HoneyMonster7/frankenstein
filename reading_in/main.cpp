#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/visitors.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <cstdlib>
#include <boost/random.hpp>
#include <boost/graph/bipartite.hpp>

#include <string>
#include <iostream>
#include <fstream>
#include <utility>
#include <algorithm>
#include <ctime>
#include <cstdlib>

#include "reaction.h"
#include "cell.h"

using namespace boost;

int main(int argc, char* argv[])
{
	//defined in cell.h
	RandomGeneratorType generator(1);

	//for automatically creating an output file named from current time
	time_t now=time(0); 
	tm *ltm = localtime(&now);
	std::ostringstream forDateTime;
	forDateTime<<1900+ltm->tm_year<<"."<<1+ltm->tm_mon<<"."<<ltm->tm_mday<<"_"<<1+ltm->tm_hour<<":"<<1+ltm->tm_min<<":"<<1+ltm->tm_sec;
	std::string dateForFileName=forDateTime.str();

	std::string dirCommand="mkdir "+dateForFileName;
	const int dir_err = system(dirCommand.c_str());
	if (-1 == dir_err)
	{
		std::cout<<"Error creating directory!"<<std::endl;
		exit(1);
	}

	std::ofstream improvementlog;
	improvementlog.open(dateForFileName+"/"+dateForFileName+".fitt");
	//set the number of internalmetabolites here:
	
	int nrOfInternalMetabolites=13;
	reaction::nrOfInternalMetabolites=nrOfInternalMetabolites;

	std::cout<<"Tests begin."<<std::endl;

	std::vector<reaction> reacVector;
	std::vector<substrate> substrateVector;


	try
	{
		//reaction::readCompounds(DATA_PATH "compounds_list__4C_v3_2_2_ext_100.dat",lofasz,compoundVList);
		reaction::readCompounds(DATA_PATH "fullnewcompounds.dat", substrateVector);

		std::cout<<"length of the reaction vector is: "<<reacVector.size()<<std::endl;
		//reaction::readReactions(DATA_PATH "reactions__4C_v3_2_2_ext_100.dat", reacVector,lofasz,reacVList,compoundVList);
		reaction::readReactions(DATA_PATH "fullnewreactions.dat", reacVector, substrateVector);
		std::cout<<"length of the reaction vector is: "<<reacVector.size()<<std::endl;
		std::cout<<"length of the substrate vector is: "<<substrateVector.size()<<std::endl;
	}
	catch(std::runtime_error& e)
	{
		std::cout << e.what() << std::endl;
		return 1;
	}

	if(reacVector.size() < 300)
	{
		std::cout << "not enough data read." << std::endl;
		return 1;
	}

	reacVector[13].printReaction();
	reacVector[209].printReaction();

	std::cout<<"Building the neighbour list..."<<std::endl;
	substrate::buildNeighbourList(reacVector,substrateVector);
	std::cout<<"Neighbour list built."<<std::endl;

	
	reacVector[129].printNeighbours();
	reacVector[5108].printNeighbours();

	environment currentEnvironment;
	currentEnvironment.atpCont=1e-1;
	currentEnvironment.adpCont=1e-2;
	currentEnvironment.ampCont=1e-4;
	currentEnvironment.nadoxcont=1e-1;
	currentEnvironment.nadredcont=1e-4;
	currentEnvironment.piCont=1e-3;
	currentEnvironment.ppiCont=1e-3;
	currentEnvironment.co2cont=1e-5;
	currentEnvironment.nh3aqCont=1e-5;
	currentEnvironment.glutCont=1e-2;
	currentEnvironment.oxo2Cont=1e-3;


	for (reaction& current : reacVector){

		current.recalcEchange(currentEnvironment);
	}


	std::vector<int> subset= {417,884,1070,448,816,629};


	cell::nrOfInternalMetabolites=nrOfInternalMetabolites;
	cell::reactionVector=reacVector;
	cell::substrateVector=substrateVector;
	//this is the k value for the fitness function
	cell::smallKforFitness=1e-2;
	//don't add nrofinternalmetabolites here
	cell::sourceSubstrate=104;
	cell::sinkSubstrate=54;


	cell trialcell(subset);

	std::vector<cell> cellVector(100,trialcell);
	//int compsize=compoundVList.size();
	double previousFittness=trialcell.getPerformance();
	std::vector<int> previousNetwork=trialcell.getReacs();


	std::string fileName="initial";
		trialcell.printXGMML(fileName);

		std::cout<<"Adding&deleting tests."<<std::endl;
		for (int k=0; k<2000; k++){
			//for testing
			//std::cout<<"Current reactions:";
			//std::vector<int> currentreacs=trialcell.getReacs();
			//for(auto i:currentreacs){std::cout<<i<<" ";}
			//std::cout<<k<<", "<<std::endl;
			if (k%1000 ==0){
				std::cout<<k<<" iterations done..."<<std::endl;
				std::cout<<"Performance is; "<<trialcell.getPerformance()<<std::endl;
			}
		trialcell.mutate(generator);
		//trialcell.printHumanReadable(compoundVList);

		//trialcell.calcThroughput();
		}
		std::cout<<"Final fitness is: "<<trialcell.getPerformance()<<std::endl;


		//std::vector<int> needMore, needLess, currentReactions;
		//int targetCompound = 54; //target is pyruvate - 54

		//needLess.emplace_back((int)-6);
		//needMore.emplace_back((int)-7);
		//currentReactions.emplace_back(11790);
		//reacVector[11790].printReaction();

		// this is the pathfinding algorithm, not used now
		//cell::findThePaths(needMore, needLess, currentReactions, targetCompound, reacVector, substrateVector, dateForFileName);


		for (int k=0; k<8000000; k++){
			cell::mutatePopulation(cellVector,generator);
			if (k%2500==0){
				
			std::cout<<k<<": ";
			cell currentBest=cell::printNFittest(cellVector,10);

			if (previousFittness+0.5<currentBest.getPerformance()){
				cell previousBest=cell(previousNetwork);
				std::ostringstream fname;
				fname<<dateForFileName<<"/"<<dateForFileName<<"_before_step_"<<k;
				previousBest.printXGMML(fname.str());
				fname.str("");
				fname.clear();
				fname<<dateForFileName<<"/"<<dateForFileName<<"_after_step_"<<k;
				currentBest.printXGMML(fname.str());
			}
			previousFittness=currentBest.getPerformance();
			previousNetwork=currentBest.getReacs();
			improvementlog<<k<<" "<<previousFittness<<std::endl;
			
			//cell::printPopulationFittnesses(cellVector);

			}

		}

		cell::printPopulationFittnesses(cellVector);
		int N=10;
		std::vector<cell> bestCells=cell::getBestNCells(cellVector,N);
		for (int i=0; i<N; i++){
			std::ostringstream forFileName;
			forFileName<<dateForFileName<<"/NR"<<i+1<<"cell";
			bestCells[i].printXGMML(forFileName.str());
		}
		

		fileName="final";
		trialcell.printXGMML(fileName);

		improvementlog.close();
  std::cout<<"Tests completed."<<std::endl;
}
