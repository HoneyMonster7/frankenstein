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

#include "reaction.h"
//#include "cell.h"

using namespace boost;

int main(int argc, char* argv[])
{
	//defined in cell.h
	RandomGeneratorType generator(1);

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

	substrate::buildNeighbourList(reacVector,substrateVector);

	
	reacVector[129].printNeighbours();
	reacVector[5109].printNeighbours();

//	auto verts = vertices(lofasz);
//	size_t count = std::distance(verts.first, verts.second);
//
//	std::cout<<"The graph has "<<count<<" vertices."<<std::endl;
//
//
//	typedef graph_traits <ReactionNetwork> traits;
//	typename traits::vertex_iterator vertex_iter, vertex_end;
//
//	bool bipartiteee = boost::is_bipartite(lofasz);
//
//	if (bipartiteee) {std::cout<<"The graph is bipartite"<<std::endl;}
//	else {std::cout<<"The graph is not bipartite."<<std::endl;}
//
//
//
//	environment currentEnvironment;
//	currentEnvironment.atpCont=1e-1;
//	currentEnvironment.adpCont=1e-2;
//	currentEnvironment.ampCont=1e-4;
//	currentEnvironment.nadoxcont=1e-2;
//	currentEnvironment.nadredcont=1e-2;
//	currentEnvironment.piCont=1e-3;
//	currentEnvironment.ppiCont=1e-3;
//	currentEnvironment.co2cont=1e-5;
//	currentEnvironment.nh3aqCont=1e-5;
//	currentEnvironment.glutCont=1e-2;
//	currentEnvironment.oxo2Cont=1e-3;
//
//
//	for (Vertex current : reacVList){
//
//		lofasz[current].reac.recalcEchange(currentEnvironment);
//	}
//
//
//	std::vector<int> subset= {130,285,5107,5109};
//	std::vector<Vertex> testReacList = cell::subsetVertices(subset,reacVList);
//
//
//	cell::allTheReactions=lofasz;
//	cell::reactionVertexList=reacVList;
//	cell::substrateVertexList=compoundVList;
//	cell::internalMetaboliteVList=internals;
//	cell::nrOfInternalMetabolites=nrOfInternalMetabolites;
//	//this is the k value for the fitness function
//	cell::smallKforFitness=1e-2;
//	cell trialcell(subset);
//
//	//int compsize=compoundVList.size();
//
//		trialcell.printCytoscape(internals);
//
//		std::cout<<"Adding&deleting tests."<<std::endl;
//		for (int k=0; k<20; k++){
//			//for testing
//			//std::cout<<"Current reactions:";
//			//std::vector<int> currentreacs=trialcell.getReacs();
//			//for(auto i:currentreacs){std::cout<<i<<" ";}
//			//std::cout<<k<<", "<<std::endl;
//		trialcell.mutate(generator,internals);
//		//trialcell.printHumanReadable(compoundVList);
//
//		trialcell.calcThroughput();
//		}
//
//		trialcell.printCytoscape(internals);
//	//cell::calcThroughput(compsize,lofasz,testReacList);
//
	std::cout<<"Tests completed."<<std::endl;
}
