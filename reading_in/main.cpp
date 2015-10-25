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

#include "reaction.h"

using namespace boost;

int main(int argc, char* argv[])
{
	std::cout<<"Tests begin."<<std::endl;

	std::vector<reaction> reacVector;

	std::vector<Vertex> reacVList, compoundVList;
	ReactionNetwork lofasz;

	try
	{
		//reaction::readCompounds(DATA_PATH "compounds_list__4C_v3_2_2_ext_100.dat",lofasz,compoundVList);
		reaction::readCompounds(DATA_PATH "fullnewcompounds.dat", lofasz, compoundVList);

		std::cout<<"length of the vector is: "<<reacVector.size()<<std::endl;
		//reaction::readReactions(DATA_PATH "reactions__4C_v3_2_2_ext_100.dat", reacVector,lofasz,reacVList,compoundVList);
		reaction::readReactions(DATA_PATH "fullnewreactions.dat", reacVector, lofasz, reacVList, compoundVList);
		std::cout<<"length of the vector is: "<<reacVector.size()<<std::endl;
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
	reacVector[299].printReaction();


	auto verts = vertices(lofasz);
	size_t count = std::distance(verts.first, verts.second);

	std::cout<<"The graph has "<<count<<" vertices."<<std::endl;


	typedef graph_traits <ReactionNetwork> traits;
	typename traits::vertex_iterator vertex_iter, vertex_end;

	bool bipartiteee = boost::is_bipartite(lofasz);

	if (bipartiteee) {std::cout<<"The graph is bipartite"<<std::endl;}
	else {std::cout<<"The graph is not bipartite."<<std::endl;}

	reaction newreaction;

	reaction evennewer;
	evennewer.printReaction();
	evennewer=newreaction;

	newreaction.printReaction();
	evennewer.printReaction();


	std::cout<<"Tests completed."<<std::endl;
}
