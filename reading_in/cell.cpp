#include "cell.h"




cell::cell(std::vector<int>& tmpAvailReacs)

		 :availableReactions(tmpAvailReacs)

{}

cell::cell() {}


void cell::printReacs() {

	std::cout<<"Current reactions: ";
	for (auto i : availableReactions){

		std::cout<<i<<", ";
	}
	std::cout<<" END"<<std::endl;
}

std::vector<int> cell::canBeAdded(ReactionNetwork& allReacs, std::vector<Vertex>& Vertexlist, std::vector<Vertex>& internals){



	//set to keep the substrate vertices that are currently used in the network
	std::set<Vertex> substrateSet;

	for (int i=0; i<availableReactions.size();i++){
		std::cout<<"Reaction number:"<<i<<std::endl;

	typename boost::graph_traits<ReactionNetwork>::adjacency_iterator vi, vi_end;


	//for testing
	allReacs[Vertexlist[availableReactions[i]-1]].reac.printReaction();
	std::cout<<"Compounds:";

	//iterating through the vertices connected to the current reaction (the compounds taking part)
	for (boost::tie(vi,vi_end)=boost::adjacent_vertices(Vertexlist[availableReactions[i]-1],allReacs); vi!=vi_end; ++vi){

		//adding the substrate vertex to the list (only unique elements are kept)
		substrateSet.insert(vi.dereference());

	}

	//testin
	std::cout<<std::endl;
	}

	//getting rid of the internal metabolites in there
	for(auto metab:internals){ substrateSet.erase(metab);}


	//the set that the possible new (and existing) reactions will be stored in
	std::set<Vertex> setOfNewReactions;

	while(!substrateSet.empty()){
		int substrateid;
		substrateid=allReacs[*substrateSet.begin()].sub.index;
		std::cout<<", "<<substrateid ;

		//looping through the substrates currently in play
		typename boost::graph_traits<ReactionNetwork>::adjacency_iterator inner, inner_end;


		for(boost::tie(inner,inner_end)=boost::adjacent_vertices(*substrateSet.begin(),allReacs); inner!=inner_end; ++inner){

			setOfNewReactions.insert(inner.dereference());
		}

		substrateSet.erase(substrateSet.begin());
	}



	std::vector<int> tobereturned;
	while(!setOfNewReactions.empty()){

		//testing
		//allReacs[*setOfNewReactions.begin()].reac.printReaction();
		int reacnumber=allReacs[*setOfNewReactions.begin()].reac.getListNr();
		tobereturned.emplace_back(reacnumber);
		setOfNewReactions.erase(setOfNewReactions.begin());
	}
		return tobereturned;

}

void cell::mutate( double probToAdd, double probToDel, ReactionNetwork& allReacs){

	boost::random::mt19937 gen{static_cast<std::uint32_t>(1)};

	std::cout<<"Random numbers:";
	for (int i=1; i<50; i++){
		std::cout<<gen()<<", ";
	}
	std::cout<<std::endl;
	double doWeAdd;
}

