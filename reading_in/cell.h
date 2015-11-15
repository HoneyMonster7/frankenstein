#pragma once


#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/visitors.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/generator_iterator.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/uniform_int.hpp>
#include <cstdlib>
#include <boost/graph/bipartite.hpp>



#include <string>
#include <iostream>
#include <fstream>
#include <utility>
#include <algorithm>
#include <array>
#include <math.h>

#include "reaction.h"

typedef boost::mt19937 RandomGeneratorType;
typedef boost::uniform_int<> UniIntDistType;
typedef boost::variate_generator<RandomGeneratorType&, UniIntDistType> Gen_Type;
typedef boost::uniform_real<> UniRealDistType;
typedef boost::variate_generator<RandomGeneratorType&, UniRealDistType> RealGenType;

class cell{

	std::vector<int> availableReactions;
	double performance;
	private: double firstPerformance;
	public: static ReactionNetwork allTheReactions;

			static std::vector<Vertex> reactionVertexList;
			static std::vector<Vertex> substrateVertexList;
			static std::vector<Vertex> internalMetaboliteVList;
			static int nrOfInternalMetabolites;
			static double smallKforFitness;

	//public: 
	cell(std::vector<int>& availableReactons);
	cell();
	
	inline std::vector<int> getReacs() {return availableReactions;}
	inline double getPerformance() {return performance;}

	void mutate(RandomGeneratorType& generator , std::vector<Vertex>& internals);

	void printReacs();

	 std::vector<int> canBeAdded(std::vector<Vertex>& internals);


	double calcThroughput();
	static std::vector<Vertex> subsetVertices( std::vector<int> vertexIDs, std::vector<Vertex> reacList);
	static std::string niceSubstrateName(Vertex currentVertex);


	void printHumanReadable(std::vector<Vertex>& substrateList);
	void printCytoscape(std::vector<Vertex> internals);

	private: int randomIntInRange(RandomGeneratorType& generator, int maxNumber);
			 double randomRealInRange(RandomGeneratorType& generator, double maxNumber);

};
