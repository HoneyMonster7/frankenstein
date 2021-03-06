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
#include <queue>
//for the command line flags
#include <unistd.h>
#include "reaction.h"
#include "cell.h"
using namespace boost;
int main(int argc, char **argv)
{
	int seedForGenerator=1;
	int c;
	std::string jobName;
	double probOfPointMut=1.0;
	double probOfHorizGene=0;
	double smallK=1e-3;
	double probOfSinkMutation=0;
	while ((c = getopt(argc,argv,"hs:j:p:g:k:a:")) != -1)
		switch(c)
		{
			case 's':
				seedForGenerator=std::stoi(optarg);
				break;
			case 'p':
				probOfPointMut=std::stod(optarg);
				break;
			case 'g':
				probOfHorizGene=std::stod(optarg);
				break;
			case 'k':
				smallK=std::stod(optarg);
				break;
			case 'a':
				probOfSinkMutation=std::stod(optarg);
				break;
			case 'h':
				std::cout<<"Accepted options:"<<std::endl;
				std::cout<<"\t -s [intSeed] seed for the MersenneTwister, default is 1, max is 2147483647"<<std::endl;
				std::cout<<"\t -p [probOfPointMut] probablity of point mutations (double) should be beteen 0 and 1. Default:1"<<std::endl;
				std::cout<<"\t -g [probOfHorizGene] probablity of horzontal gene transfer (double) should be beteen 0 and 1, preferably lower than probOfPointMut. Default:0"<<std::endl;
				std::cout<<"\t -k [smallKforFitness] cost of a reaction within the reaction network. Default:0.001"<<std::endl;
				std::cout<<"\t -a [probOfSinkMutation] probability of a sink mutation (addition or deletion with equal probability). Default:0"<<std::endl;
				std::cout<<"\t -j [jobName] name for the folder to output the results into. If skipped a timestamped directory will be used for this."<<std::endl;
				exit(0);
				break;
			case 'j':
				jobName=optarg;
				break;
			case '?':
				if (optopt =='s')
					std::cout<<"Option -s requires an integer argument"<<std::endl;
				else if (optopt =='j')
					std::cout<<"Option -j requires the desired folder name (jobName)"<<std::endl;
				else if (optopt =='p')
					std::cout<<"Option -p requires the desired probality of point mutations (double)"<<std::endl;
				else if (optopt =='g')
					std::cout<<"Option -g requires the desired probality of horizontal gene transfer"<<std::endl;
				else if (optopt =='k')
					std::cout<<"Option -k requires the desired cost for reactions"<<std::endl;
				else if (optopt =='a')
					std::cout<<"Option -a requires the desired probability"<<std::endl;
				else
					std::cout<<"Unknown option, try -h for allowed options."<<std::endl;
				return 1;
			default:
				exit(1);
		}
	RandomGeneratorType generator(seedForGenerator);
	//for automatically creating an output file named from current time
	time_t now=time(0); 
	tm *ltm = localtime(&now);
	std::ostringstream forDateTime;
	forDateTime<<1900+ltm->tm_year<<"."<<1+ltm->tm_mon<<"."<<ltm->tm_mday<<"_"<<1+ltm->tm_hour<<":"<<1+ltm->tm_min<<":"<<1+ltm->tm_sec;
	std::string dateForFileName=forDateTime.str();
	std::string actualFilename;
	if (!jobName.empty()){
		actualFilename=jobName;
	}
	else {
		actualFilename=dateForFileName;
	}
	std::string dirCommand="mkdir "+actualFilename;
	const int dir_err = system(dirCommand.c_str());
	if (-1 == dir_err)
	{
		std::cout<<"Error creating directory!"<<std::endl;
		exit(1);
	}
	std::cout<<"The chosen probablity for point mutations is: "<<probOfPointMut<<" and for horizontal gene transfer it is: "<<probOfHorizGene<<std::endl;
	std::ofstream improvementlog;
	improvementlog.open(actualFilename+"/"+actualFilename+".fitt");
	//set the number of internalmetabolites here:
	int nrOfInternalMetabolites=13;
	reaction::nrOfInternalMetabolites=nrOfInternalMetabolites;
	std::cout<<"Tests begin."<<std::endl;
	std::vector<reaction> reacVector;
	std::vector<substrate> substrateVector;
	try
	{
		reaction::readCompounds(DATA_PATH "fullnewcompounds.dat", substrateVector);
		std::cout<<"length of the reaction vector is: "<<reacVector.size()<<std::endl;
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
	std::cout<<"Building the neighbour list..."<<std::endl;
	substrate::buildNeighbourList(reacVector,substrateVector);
	std::cout<<"Neighbour list built."<<std::endl;
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
	//setting static member variables 
	cell::nrOfInternalMetabolites=nrOfInternalMetabolites;
	cell::reactionVector=reacVector;
	cell::substrateVector=substrateVector;
	//this is the k value for the fitness function
	cell::smallKforFitness=smallK;
	//setting the probabilities for mutations
	cell::probabilityOfMutation=probOfPointMut;
	cell::probabilityOfHorizontalGenetransfer=probOfHorizGene;
	cell::sourceSubstrate.push_back(104);
	cell::sinkSubstrate.push_back(54);
	//creating the original cell here
	cell trialcell(subset);
	//Restructuring the population mutations 
	int numberOfCells=100;
	std::vector<int> populationIndex(numberOfCells,0);
	std::vector<int> howManyOfEach(numberOfCells,0);
	//cell type vector
	std::vector<cell> cellVector(numberOfCells);
	cellVector[0]=trialcell;
	howManyOfEach[0]=100;

	std::string fileName="initial";
	trialcell.printXGMML(fileName);
	std::cout<<"Initial performance is; "<<trialcell.getPerformance()<<std::endl;
	// this is the pathfinding algorithm, not used now
	//cell::findThePaths(needMore, needLess, currentReactions, targetCompound, reacVector, substrateVector, actualFilename);
		double probabilityOfAnyMutation=cell::probabilityOfMutation+cell::probabilityOfHorizontalGenetransfer;
		int numberOfMutationsToSimulate=60000000;
		double dcheckpointLenght=(numberOfMutationsToSimulate/probabilityOfAnyMutation)/(10*numberOfCells);
		int NRofCheckpoints=10;
		int checkPointLength=(int)dcheckpointLenght;
		const int generationsPerWriteout=1000;
		double previousAvgFittness=cellVector[0].getPerformance();
		std::cout<<"To simulate "<<numberOfMutationsToSimulate<<" mutations, the checkPointLength has been set to "<<checkPointLength<<std::endl;
		//defining the queues here for the intermittent steps between logging
		double maxFittQueue [generationsPerWriteout];
		double avgFittQueue [generationsPerWriteout];
		double entropyQueue [generationsPerWriteout];
		int bestNetSizeQueue [generationsPerWriteout];
		double avgNetSizeQueue [generationsPerWriteout];
		int bestUsedReacsQueue [generationsPerWriteout];
		double avgUsedReacsQueue [generationsPerWriteout];
		int arrayPos=0;
		std::time_t startTime=std::time(nullptr);
		//outer loop is there in order to save networks every 10% of the simulation
		for (int outerLoop=0; outerLoop<NRofCheckpoints; outerLoop++){
			for (int k=0; k<checkPointLength; k++){
				for (int iter=0; iter<numberOfCells; ++iter){
					//running a generation here - as many mutations as many cells there are in the population
					cell::mutatePopulation(populationIndex,howManyOfEach,cellVector,generator);
				}
				cell::printProgressFile(populationIndex,cellVector,howManyOfEach,k,outerLoop,generationsPerWriteout,checkPointLength,improvementlog,previousAvgFittness,maxFittQueue,avgFittQueue,entropyQueue,bestNetSizeQueue,avgNetSizeQueue,bestUsedReacsQueue,avgUsedReacsQueue);
			}
			std::time_t currentTime=std::time(nullptr);
			double timeElapsed=std::difftime(currentTime,startTime);
			std::cout<<"Time for Checkpointing! "<<timeElapsed<<" seconds have passed since the start, that's "<<(outerLoop+1)*checkPointLength/timeElapsed<<" generations per seconds."<<std::endl;
			cell currentBest=cell::printNFittest(populationIndex,cellVector,numberOfCells);
			std::vector<cell> bestCells=cell::getBestNCells(populationIndex,cellVector,cellVector.size());
			if (outerLoop != NRofCheckpoints){
				//the checkpoints will be saved into separate folders
				//creating a directory for saving the checkpoints into
				std::ostringstream forCheckpointFolder;
				forCheckpointFolder<<actualFilename<<"/CP"<<outerLoop+1;
				std::string dirCommand="mkdir "+forCheckpointFolder.str();
				const int dir_err = system(dirCommand.c_str());
				for (int i=0; i<bestCells.size(); i++){
					std::ostringstream forFileName;
					forFileName<<forCheckpointFolder.str()<<"/"<<actualFilename<<"CP"<<outerLoop+1<<"NR"<<i+1<<"cell";
					bestCells[i].printXGMML(forFileName.str());
				}
			}
		}
		improvementlog.close();
  std::cout<<"Tests completed."<<std::endl;
}
