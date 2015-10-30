#pragma once


#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/visitors.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random.hpp>
#include <boost/random/uniform_real.hpp>
#include <cstdlib>
#include <boost/graph/bipartite.hpp>



#include <string>
#include <iostream>
#include <fstream>
#include <utility>
#include <algorithm>
#include <array>

#include "reaction.h"

typedef boost::mt19937 RandomGeneratorType;

class cell{


	std::vector<int> availableReactions;


	public: 
	cell(std::vector<int>& availableReactons);
	cell();
	
	inline std::vector<int> getReacs() {return availableReactions;}

	void mutate(double probToAdd, double probToDel, ReactionNetwork& allReacs, RandomGeneratorType& generator );

	void printReacs();

	 std::vector<int> canBeAdded(ReactionNetwork& allReacs, std::vector<Vertex>& Vertexlist, std::vector<Vertex>& internals);
};
