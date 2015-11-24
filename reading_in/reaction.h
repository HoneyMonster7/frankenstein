#pragma once


// header file for the reaction class


#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/visitors.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <cstdlib>
#include <boost/graph/bipartite.hpp>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/uniform_int.hpp>
//for the glpk library
#include "stdio.h"
#include "stdlib.h"
#include "glpk.h"

#include <string>
#include <iostream>
#include <fstream>
#include <utility>
#include <algorithm>
#include <array>

typedef boost::mt19937 RandomGeneratorType;
typedef boost::uniform_int<> UniIntDistType;
typedef boost::variate_generator<RandomGeneratorType&, UniIntDistType> Gen_Type;
typedef boost::uniform_real<> UniRealDistType;
typedef boost::variate_generator<RandomGeneratorType&, UniRealDistType> RealGenType;


struct environment
{

	//substrates from newsortedcompounds.txt in this order
	//in the internalMets array they are kept in this order

	double nh2acceptorCont=1;
	double nh2donorCont=1;
	double conh22Cont=1;
	double nh3aqCont=1;

	double ppiCont=1;
	double piCont=1;
	double atpCont=1;
	double adpCont=1;
	double ampCont=1;
	double nadredcont=1;
	double nadoxcont=1;
	double co2cont=1;
	double h2ocont=1;


	double temperature=293; //temperature in kelvins
	double glutCont=1;
	double oxo2Cont=1;
};

class reaction;

class substrate
{
	//properties of each substrate
	int index;
	double freeOfCreation;
	std::string molecule;
	std::string name;
	int charge;
	std::vector<int> involvedInReacs;

	public:

	substrate(int index, double freeOfCreation,std::string molecule, std::string name, int charge);
	substrate();
	inline std::vector<int> getinvolved() { return involvedInReacs;}
	void addInvolved(int addThisReac);

	static void buildNeighbourList(std::vector<reaction>& reacList, std::vector<substrate>& substrateVector);
	void printInvolved();
	std::string niceSubstrateName();
};




typedef std::array<int, 13> InternalMetsT;


class reaction
{
	//properties of each reaction, should include an int (possible shorter) for each internal metabolite to aid   the recalculation of freeEchange for different environments
	//      std::string type;

	int listNR;
	std::vector<int> substrates;
	std::vector<int> products;
	InternalMetsT internalMets = {};
	double freeEChange = 0.0;
	//std::string humanReadable;
	double currentFreeEChange=0;
	std::vector<int> neighbourOf;
	public: static int nrOfInternalMetabolites;

public:

	reaction(int listNR,double tmpfreechange, const std::vector<int>& tmpsubstrates, const std::vector<int>& tmproducts, const InternalMetsT& tmpInternalMets);
	reaction();

	static void readReactions(std::string fileName, std::vector<reaction>& reacPointer, std::vector<substrate>& substrateVector);

	static void readCompounds(std::string fileName,  std::vector<substrate>& substrateVector);

	void printReaction();
	double freeEchange();
	void recalcEchange(const environment& env);
	void printNeighbours();

	// inlining the return of big vectors may help the compiler
	// to get rid of some needless copying
	inline std::vector<int> getsubstrates() { return substrates; }
	inline std::vector<int> getproducts() { return products; }
	inline std::vector<int> getNeighbours() { return neighbourOf; }
	inline double getCurrentFreeEChange() {return currentFreeEChange;}
	int getListNr();
	void setNeighbours( std::vector<int>& neighbourList);
};

