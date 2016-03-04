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
#include <time.h>
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

	while ((c = getopt(argc,argv,"hs:j:")) != -1)
		switch(c)
		{
			case 's':
				//std::cout<<"Input is: "<<optarg<<std::endl;
				seedForGenerator=std::stoi(optarg);
				//std::cout<<"Parsed, it is: "<<seedForGenerator<<std::endl;
				break;
			case 'h':
				std::cout<<"Accepted options:"<<std::endl;
				std::cout<<"\t -s [intSeed] seed for the MersenneTwister, default is 1, max is 2147483647"<<std::endl;
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
				else
					std::cout<<"Unknown option, try -h for allowed options."<<std::endl;
				return 1;
			default:
				exit(1);
		}



		

	//defined in cell.h
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
	cell::smallKforFitness=1e-3;
	//setting the probabilities for mutations
	cell::probabilityOfMutation=0.01;
	cell::probabilityOfHorizontalGenetransfer=0;
	//don't add nrofinternalmetabolites here
	cell::sourceSubstrate.push_back(104);
	//cell::sourceSubstrate.push_back(267);
	cell::sinkSubstrate.push_back(54);


	cell trialcell(subset);

	//std::vector<cell> cellVector(100,trialcell);

	//Restructuring the population mutations 
	int numberOfCells=100;
	std::vector<int> populationIndex(numberOfCells,0);
	std::vector<int> howManyOfEach(numberOfCells,0);

	std::vector<cell> cellVector(numberOfCells);
	cellVector[0]=trialcell;
	howManyOfEach[0]=100;
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
		//cell::findThePaths(needMore, needLess, currentReactions, targetCompound, reacVector, substrateVector, actualFilename);


		int NRofCheckpoints=10;
		int checkPointLength=6000000;
		const int generationsPerWriteout=10000;
		//43200 seconds is 12 hours
		double timeBetweenTwoWriteouts=43200;
		double previousAvgFittness=cellVector[0].getPerformance();


		//defining the queues here
		double maxFittQueue [generationsPerWriteout];
		double avgFittQueue [generationsPerWriteout];
		double entropyQueue [generationsPerWriteout];
		int bestNetSizeQueue [generationsPerWriteout];
		double avgNetSizeQueue [generationsPerWriteout];
		int bestUsedReacsQueue [generationsPerWriteout];
		double avgUsedReacsQueue [generationsPerWriteout];
		int arrayPos=0;


		int outerLoop=0;


		std::time_t previousWriteoutTime=std::time(nullptr);
		//outer loop is there in order to save networks every 10% of the simulation
		//for (int outerLoop=0; outerLoop<NRofCheckpoints; outerLoop++){
		while (true){

			for (int k=0; k<checkPointLength; k++){

				for (int iter=0; iter<numberOfCells; ++iter){

					//running a generation here - as many mutations as many cells there are in the population
					cell::mutatePopulation(populationIndex,howManyOfEach,cellVector,generator);
				}
				//if (k%generationsPerWriteout==0){
				//	
				//std::cout<<k+outerLoop*checkPointLength<<": ";
				//cell currentBest=cell::printNFittest(populationIndex,cellVector,10);

				//if (previousFittness+0.5<currentBest.getPerformance()){
				//	cell previousBest=cell(previousNetwork);
				//	std::ostringstream fname;
				//	fname<<actualFilename<<"/"<<actualFilename<<"_before_step_"<<k;
				//	previousBest.printXGMML(fname.str());
				//	fname.str("");
				//	fname.clear();
				//	fname<<actualFilename<<"/"<<actualFilename<<"_after_step_"<<k;
				//	currentBest.printXGMML(fname.str());
				//}
				//previousFittness=currentBest.getPerformance();
				//previousNetwork=currentBest.getReacs();


				//double enthropy=0;

				//for(auto element:howManyOfEach){

				//	if(element!=0){
				//		enthropy+=element*std::log(element);
				//	}
				//}

				//improvementlog<<k+outerLoop*checkPointLength<<" "<<previousFittness<<" "<<-1*enthropy<<std::endl;
				cell::printProgressFile(populationIndex,cellVector,howManyOfEach,k,outerLoop,generationsPerWriteout,checkPointLength,improvementlog,previousAvgFittness,maxFittQueue,avgFittQueue,entropyQueue,bestNetSizeQueue,avgNetSizeQueue,bestUsedReacsQueue,avgUsedReacsQueue);
				
				//cell::printPopulationFittnesses(cellVector);

				//}

			}

			//cell::printPopulationFittnesses(populationIndex,cellVector);

			cell currentBest=cell::printNFittest(populationIndex,cellVector,numberOfCells);
			//int N=10;
			std::vector<cell> bestCells=cell::getBestNCells(populationIndex,cellVector,cellVector.size());


			std::time_t currentTime=std::time(nullptr);
			double timeElapsed=std::difftime(currentTime,previousWriteoutTime);
			if (timeElapsed>timeBetweenTwoWriteouts){

				std::cout<<timeElapsed<<" seconds have passed ("<<timeElapsed/(double)3600<<" hours), time for Checkpointing!"<<std::endl;
				//the checkpoints will be saved into separate folders
				
				//creating a directory for saving the checkpoints into
				std::ostringstream forCheckpointFolder;
				forCheckpointFolder<<actualFilename<<"/CP"<<outerLoop+1;
				std::string dirCommand="mkdir "+forCheckpointFolder.str();
				const int dir_err = system(dirCommand.c_str());

				//the script starting reaction will make sure the folder is empty
				//
				//const int dir_err = system(dirCommand.c_str());
				//if (-1 == dir_err)
				//{
				//	std::cout<<"Error creating directory!"<<std::endl;
				//	exit(1);
				//}
				for (int i=0; i<bestCells.size(); i++){
					std::ostringstream forFileName;
					forFileName<<forCheckpointFolder.str()<<"/"<<actualFilename<<"CP"<<outerLoop+1<<"NR"<<i+1<<"cell";
					bestCells[i].printXGMML(forFileName.str());
				}
				previousWriteoutTime=currentTime;
			}
			else {
				//final network doesn't need a checkpoint folder - Yes it does.
				//for (int i=0; i<N; i++){
				//	std::ostringstream forFileName;
				//	forFileName<<actualFilename<<"/"<<actualFilename<<"CP"<<outerLoop+1<<"NR"<<i+1<<"cell";
				//	bestCells[i].printXGMML(forFileName.str());
				//}
			}
			++outerLoop;
		}
		

		fileName="final";
		//trialcell.printXGMML(fileName);

		improvementlog.close();
  std::cout<<"Tests completed."<<std::endl;
}
