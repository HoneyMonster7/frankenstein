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
			 std::vector<double> fluxThroughReacs;
	public: 

			static int sourceSubstrate,sinkSubstrate;
			static std::vector<reaction> reactionVector;
			static std::vector<substrate> substrateVector;
			static int nrOfInternalMetabolites;
			static double smallKforFitness;

	//public: 
	cell(std::vector<int>& availableReactons);
	cell();
	
	inline std::vector<int> getReacs() {return availableReactions;}
	inline double getPerformance() {return performance;}

	inline std::vector<double> getFluxes() {return fluxThroughReacs;}
	void mutate(RandomGeneratorType& generator);
	cell  mutateAndReturn(RandomGeneratorType& generator);

	void printReacs();

	 std::vector<int> canBeAdded();


	void setFluxes(std::vector<double>& fluxVector);
	double calcThroughput();
	//static std::vector<Vertex> subsetVertices( std::vector<int> vertexIDs, std::vector<Vertex> reacList);

	static void mutatePopulation(std::vector<cell>& population, RandomGeneratorType& generator);
	static std::vector<double> getPopulationFittness(std::vector<cell>& population);
	static void printPopulationFittnesses(std::vector<cell>& population);

	void printHumanReadable();
	void printXGMML(std::string fileName);

	private: static int randomIntInRange(RandomGeneratorType& generator, int maxNumber);
			 static double randomRealInRange(RandomGeneratorType& generator, double maxNumber);

};
