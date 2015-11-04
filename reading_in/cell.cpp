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
	//	std::cout<<"Reaction number:"<<i<<std::endl;

	typename boost::graph_traits<ReactionNetwork>::adjacency_iterator vi, vi_end;


	//for testing
	//allReacs[Vertexlist[availableReactions[i]-1]].reac.printReaction();
	//std::cout<<"Compounds:";

	//iterating through the vertices connected to the current reaction (the compounds taking part)
	for (boost::tie(vi,vi_end)=boost::adjacent_vertices(Vertexlist[availableReactions[i]-1],allReacs); vi!=vi_end; ++vi){

		//adding the substrate vertex to the list (only unique elements are kept)
		substrateSet.insert(vi.dereference());

	}

	//testin
	//std::cout<<std::endl;
	}

	//getting rid of the internal metabolites in there
	for(auto metab:internals){ substrateSet.erase(metab);}


	//the set that the possible new (and existing) reactions will be stored in
	std::set<Vertex> setOfNewReactions;

	while(!substrateSet.empty()){
		int substrateid;
		substrateid=allReacs[*substrateSet.begin()].sub.index;
		//std::cout<<", "<<substrateid ;

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
	//uncomment for printing out the line numbers of possible reactions to be added 
	//for (auto whatever:tobereturned){
	//	std::cout<<whatever+1<<std::endl;
	//}
		return tobereturned;

}

void cell::mutate( double probToAdd, double probToDel, ReactionNetwork& allReacs, RandomGeneratorType& generator, std::vector<Vertex>& Vertexlist, std::vector<Vertex>& internals ){


//	std::cout<<"Random numbers:";
//	for (int i=1; i<50; i++){
//		std::cout<<gen()%availableReactions.size()<<", ";
//	}


	double doWeAdd=randomRealInRange(generator, 1);
	double doWeDelete=randomRealInRange(generator,1);

	if(doWeAdd<=probToAdd){
		std::vector<int> whatCanWeAdd = canBeAdded(allReacs, Vertexlist, internals);
		int whichOneToAdd=randomIntInRange(generator,whatCanWeAdd.size()-1);
		std::cout<<"We add reaction nr: "<<whatCanWeAdd[whichOneToAdd]<<"from a possible "<<whatCanWeAdd.size()<<"reactions"<<std::endl;
		availableReactions.push_back(whatCanWeAdd[whichOneToAdd]);
	}

	if(doWeDelete<=probToDel){

		int whichOneToDel=randomIntInRange(generator,availableReactions.size()-1);
		std::cout<<"We delete nr: "<<availableReactions[whichOneToDel]<<std::endl;
		availableReactions.erase(availableReactions.begin()+whichOneToDel);
	}
	std::cout<<"Currently in network: "<<availableReactions.size()<<" reactions."<<std::endl;
}

 int cell::randomIntInRange(RandomGeneratorType& generator, int maxNumber){

	Gen_Type intgenerator(generator, UniIntDistType(0,maxNumber));
	//next line probs not needed
	boost::generator_iterator<Gen_Type> randomInt(&intgenerator);
	int currentRandomNumber=intgenerator();
	//was used for testing
	//std::cout<<currentRandomNumber<<", ";

	return currentRandomNumber;
}

double cell::randomRealInRange(RandomGeneratorType& generator, double maxNumber){
		RealGenType realgenerator(generator, UniRealDistType(0,maxNumber));

		double currentRandomNumber=realgenerator();

		return currentRandomNumber;
		}




void cell::calcThroughput(const int NrCompounds,ReactionNetwork& graph, std::vector<Vertex> reacList){

	glp_prob *lp;
	lp=glp_create_prob();
	glp_set_prob_name(lp,"network_throughput");
	glp_set_obj_dir(lp,GLP_MAX);
	
	glp_add_rows(lp,NrCompounds);

	for (int i=1; i<=NrCompounds; i++){
		//specifying that the rows must sum to zero (flux vector in the nullspace of S matrix)
		glp_set_row_bnds(lp,i,GLP_FX,0.0,0.0);
	}
	//extra column for the imaginary reaction getting rid of the final compound (objective function)
	int listSize=reacList.size();
	glp_add_cols(lp,listSize+4);

	//for(int i=1; i<=listSize+4; i++){ glp_set_col_bnds(lp,i,GLP_LO,0.0,0.0);}

	std::vector<int> ia,ja;
	std::vector<double> ar;

	for (int i=1;i<=listSize;i++){


		//preparing the sparse matrix's values 
		reaction tmpreac=graph[reacList[i-1]].reac;


		//if(tmpreac.currentFreeEChange<0){
		glp_set_col_bnds(lp,i,GLP_DB,0.0,2.0);
		//}
		//else{glp_set_col_bnds(lp,i,GLP_UP,0.0,0.0);}


		std::vector<int> tmpsubs=tmpreac.getsubstrates();
		std::vector<int> tmpprods=tmpreac.getproducts();

		//ia.push_back(1); ja.push_back(1); ar.push_back(0.0);
		for (int j: tmpsubs){
			//using i+14 as the column numbering starts from 1, and there are 13 internal metabolites
			//with negative substrate indices
			ia.push_back(j+14);
			ja.push_back(i);
			ar.push_back(-1.0);
		}

		for (int j: tmpprods){
			//using i+14 as the column numbering starts from 1, and there are 13 internal metabolites
			//with negative substrate indices
			ia.push_back(j+14);
			ja.push_back(i);
			ar.push_back(1.0);
		}

	}

	//bounding the imaginary reactions, only certain amount of material can be taken at once
	glp_set_col_bnds(lp,listSize+1,GLP_DB,0.0,10.0);
	glp_set_col_bnds(lp,listSize+2,GLP_DB,0.0,10.0);
	glp_set_col_bnds(lp,listSize+3,GLP_DB,0.0,10.0);

	glp_set_col_bnds(lp,listSize+4,GLP_DB,-10.0,10.0);
	//add imaginary reaction here:
	ia.push_back(908+14);	ja.push_back(listSize+1); ar.push_back(1.0);
	ia.push_back(-1+14);	ja.push_back(listSize+2); ar.push_back(1.0);
	ia.push_back(-2+14);	ja.push_back(listSize+3); ar.push_back(-1.0);
	ia.push_back(14);	ja.push_back(listSize+4); ar.push_back(-1.0);

	//ia.push_back(43+14);	ja.push_back(listSize+2); ar.push_back(1.0);
	//ia.push_back(88+14);	ja.push_back(listSize+3); ar.push_back(1.0);
	//
	//target is to maximize the imaginary reactions throughput
	glp_set_obj_coef(lp,listSize+4,1.0);


	//creating the arrays now
	int length=ia.size();
	int iarray[length+1],jarray[length+1];
	double ararray[length+1];
	std::copy(ia.begin(),ia.end(),iarray+1);
	std::copy(ja.begin(),ja.end(),jarray+1);
	std::copy(ar.begin(),ar.end(),ararray+1);


	for (int i=0; i<=length; i++){


		std::cout<<"Matrix element: "<<iarray[i]<<", "<<jarray[i]<<", "<<ararray[i]<<std::endl;
	}

	std::cout<<"The length of the vectors are: "<<length<<", "<<ja.size()<<", "<<ar.size()<<std::endl;

	std::cout<<"Last elements are: "<<iarray[length]<<", "<<jarray[length]<<", "<<ararray[length]<<std::endl;
	glp_load_matrix(lp,length,iarray,jarray,ararray);

	glp_simplex(lp,NULL);

	double goodness=glp_get_obj_val(lp);

	std::cout<<"Objective function value: "<<goodness<<std::endl;
	std::cout<<"Coefficients:";
	for (int i=0; i<(listSize+4); i++){
		double colvalue=glp_get_col_prim(lp,i+1);
		std::cout<<colvalue<<", ";
	}
	std::cout<<std::endl;
}

 std::vector<Vertex> cell::subsetVertices( std::vector<int> vertexIDs, std::vector<Vertex> reacList){
	 //sort not required? 
	 std::sort(vertexIDs.begin(),vertexIDs.end());

	 std::vector<Vertex> tobeReturned;

	 for(int i :vertexIDs){
		 //uses i-1 in order to correspond to line numbers in the input reaction list
		 tobeReturned.push_back(reacList[i-1]);
	 }
	 return tobeReturned;
 }

void cell::printHumanReadable(ReactionNetwork& graph, std::vector<Vertex>& Vertexlist, std::vector<Vertex>& substratlist){

	for (int i: availableReactions){

		reaction currentReac=graph[Vertexlist[i]].reac;

	}

}
