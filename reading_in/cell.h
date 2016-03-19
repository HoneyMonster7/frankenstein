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
#include <queue>
#include <unordered_set>
#include "reaction.h"
typedef boost::mt19937 RandomGeneratorType;
typedef boost::uniform_int<> UniIntDistType;
typedef boost::variate_generator<RandomGeneratorType&, UniIntDistType> Gen_Type;
typedef boost::uniform_real<> UniRealDistType;
typedef boost::variate_generator<RandomGeneratorType&, UniRealDistType> RealGenType;
class cell{
	std::vector<int> availableReactions;
	double performance;
	std::vector<int> additionalSinks;
	private: double firstPerformance;
			 std::vector<double> fluxThroughReacs;
	public: 
			 static std::vector<int> sourceSubstrate,sinkSubstrate;
			 static std::vector<reaction> reactionVector;
			 static std::vector<substrate> substrateVector;
			 static int nrOfInternalMetabolites;
			 static double smallKforFitness;
			 static double probabilityOfMutation;
			 static double probabilityOfHorizontalGenetransfer;
			 static double probabilityOfSinkMutation;
			 cell(std::vector<int>& availableReactons);
			 //other constructor that works similarily but won't calculate fittness automatically
			 cell(std::vector<int>& availableReactons, int notUsed);
			 cell();
			 bool operator<(const cell& other) const;
			 inline std::vector<int> getReacs() {return availableReactions;}
			 inline double getPerformance() const {return performance;}
			 inline std::vector<double> getFluxes() {return fluxThroughReacs;}
			 inline std::vector<int> getAddSinks() {return additionalSinks;}
			 cell  mutateAndReturn(RandomGeneratorType& generator);
			 void printReacs();
			 std::vector<int> canBeAdded();
			 void setFluxes(std::vector<double>& fluxVector);
			 double calcThroughput();
			 static void mutatePopulation(std::vector<int>& population,std::vector<int>& howManyOfEach, std::vector<cell>& cellVector, RandomGeneratorType& generator);
			 static std::vector<double> getPopulationFittness(std::vector<int>& population, std::vector<cell>& cellVector);
			 static void printPopulationFittnesses(std::vector<int>& population, std::vector<cell>& cellVector);
			 static cell printNFittest(std::vector<int>& population, std::vector<cell>& cellVector, int N);
			 static std::vector<cell> getBestNCells(std::vector<int>& population, std::vector<cell>& cellVector, int N);
			 static void findThePaths(std::vector<int> needMore, std::vector<int> needLess, std::vector<int> currentReactions, int TargetCompound, std::vector<reaction>& ReactionVector, std::vector<substrate>& SubstrateVector, std::string fnameString);
			 void printHumanReadable();
			 void printXGMML(std::string fileName);
			 static void printProgressFile(std::vector<int>& population, std::vector<cell>& cellVector, std::vector<int>& howManyOfEach, int k, int outerloop,const int generationsPerWriteout, int checkPointLength, std::ofstream& fileToWrite,double& previousFittness, double maxFittQueue [], double avgFittQueue [], double entropyQueue [], int bestNetSizeQueue [], double avgNetSizeQueue [], int bestUsedReacsQueue [], double avgUsedReacsQueue []);
			 void deleteThisSink(int thisSinkIsToBeGone);
			 void addThisSink(int thisSinkIsNowYours);
			 void theseAreYourAddSinks(std::vector<int> theseSinksAreNowYours);
			 //Only use in mutations. take care with this, remember to calc throughput
			 void setReacs(std::vector<int> theseAreYourNewReacs);
	private: static int randomIntInRange(RandomGeneratorType& generator, int maxNumber);
			 static double randomRealInRange(RandomGeneratorType& generator, double maxNumber);
};
