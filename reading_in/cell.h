#pragma once


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

#include "reaction.h"


class cell{


	std::vector<int> availableReactions;


	public: 
	cell(std::vector<int>& availableReactons);
	cell();
	
	inline std::vector<int> getReacs() {return availableReactions;}

	void mutate(double probToAdd, double probToDel, ReactionNetwork& allReacs);

	void printReacs();

	 std::vector<int> canBeAdded(ReactionNetwork& allReacs, std::vector<Vertex>& Vertexlist, std::vector<Vertex>& internals);
};
