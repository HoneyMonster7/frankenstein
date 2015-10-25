#pragma once


// header file for the reaction class


#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/visitors.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <cstdlib>
#include <boost/graph/bipartite.hpp>

#include <string>
#include <iostream>
#include <fstream>
#include <utility>
#include <algorithm>
#include <array>

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

struct substrate
{
	//properties of each substrate
	int index;
	double freeOfCreation;
	std::string molecule;
	std::string name;
	int charge;
};


struct SubReac;


typedef boost::adjacency_list<boost::listS, boost::vecS, boost::directedS, SubReac, bool> ReactionNetwork;
typedef boost::graph_traits<ReactionNetwork>::vertex_descriptor Vertex;

typedef std::array<int, 13> InternalMetsT;


class reaction
{
	//properties of each reaction, should include an int (possible shorter) for each internal metabolite to aid   the recalculation of freeEchange for different environments
	//      std::string type;

	std::vector<int> substrates;
	std::vector<int> products;
	InternalMetsT internalMets = {};
	double freeEChange = 0.0;
	//std::string humanReadable;
	double currentFreeEChange=0;

public:

	reaction(double tmpfreechange, const std::vector<int>& tmpsubstrates, const std::vector<int>& tmproducts, const InternalMetsT& tmpInternalMets);
	reaction();

	static void readReactions(std::string fileName, std::vector<reaction>& reacPointer, ReactionNetwork& lofasz, std::vector<Vertex>& vertexList, const std::vector<Vertex>& compoundVList);

	static void readCompounds(std::string fileName, ReactionNetwork& lofasz, std::vector<Vertex>& compoundlist);

	void printReaction();
	double freeEchange();
	void recalcEchange(const environment& env);

	// inlining the return of big vectors may help the compiler
	// to get rid of some needless copying
	inline std::vector<int> getsubstrates() { return substrates; }
	inline std::vector<int> getproducts() { return products; }
};

struct SubReac
{
	substrate sub;
	reaction reac;
};
