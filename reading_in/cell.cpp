#include "cell.h"


// static class variables
ReactionNetwork cell::allTheReactions;

std::vector<Vertex> cell::reactionVertexList;
std::vector<Vertex> cell::substrateVertexList;
std::vector<Vertex> cell::internalMetaboliteVList;
int cell::nrOfInternalMetabolites;
double cell::smallKforFitness;

cell::cell(std::vector<int>& tmpAvailReacs)

		 :availableReactions(tmpAvailReacs)

{
	double initialPerformance=calcThroughput();
	performance=initialPerformance;
	firstPerformance=initialPerformance;
}


cell::cell() {}


void cell::printReacs() {

	std::cout<<"Current reactions: ";
	for (auto i : availableReactions){

		std::cout<<i<<", ";
	}
	std::cout<<" END"<<std::endl;
}

void cell::printCytoscape(std::vector<Vertex> internals){

	ReactionNetwork allReacs=allTheReactions;
	std::vector<Vertex> Vertexlist=reactionVertexList;
	std::vector<Vertex> substrateList=substrateVertexList;



	for (int i: availableReactions){

		reaction currentReac=allTheReactions[reactionVertexList[i-1]].reac;
		std::vector<int> currentsubs=currentReac.getsubstrates();
		std::vector<int> currentproducts=currentReac.getproducts();
		int reacNR=currentReac.getListNr();


		for (int sub:currentsubs){
			std::string substrateName=cell::niceSubstrateName(substrateVertexList[sub+nrOfInternalMetabolites]);

			std::cout<<substrateName<<" cr "<<reacNR<<std::endl;
		}

		for (int prod:currentproducts){
			std::string productName=cell::niceSubstrateName(substrateVertexList[prod+nrOfInternalMetabolites]);

			std::cout<<reacNR<<" rc "<<productName<<std::endl;

		}
	}





	}






std::vector<int> cell::canBeAdded(std::vector<Vertex>& internals){


	ReactionNetwork allReacs =allTheReactions;
	std::vector<Vertex> Vertexlist=reactionVertexList;

	//set to keep the substrate vertices that are currently used in the network
	std::set<Vertex> substrateSet;
	int nrOfAvailableReactions=availableReactions.size();

	for (int i=0; i<nrOfAvailableReactions;i++){
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

void cell::mutate( RandomGeneratorType& generator, std::vector<Vertex>& internals){


	ReactionNetwork allReacs=allTheReactions;
	std::vector<Vertex> Vertexlist=reactionVertexList;
	//int compoundSize= substrateVertexList.size();
//	std::cout<<"Random numbers:";
//	for (int i=1; i<50; i++){
//		std::cout<<gen()%availableReactions.size()<<", ";
//	}

	int desiredNetworkSize=10;
	int currentReactionNumber = availableReactions.size();


	int deleteThisOne,addThisOne;

	std::vector<int> trialNewCell = availableReactions;

	//double addProb= 0.1*exp(desiredNetworkSize-currentReactionNumber);
	//double delProb=0.1*exp(currentReactionNumber-desiredNetworkSize);

	double addProb=0.5;
	double delProb=0.5;

	double doWeAdd=randomRealInRange(generator, 1);
	double doWeDelete=randomRealInRange(generator,1);
	double doWeAccept=randomRealInRange(generator,1);


	//calculating current throughput
	bool areWeAdding= doWeAdd<=addProb;
	if(areWeAdding){
		std::vector<int> whatCanWeAdd = canBeAdded(internals);
		int whichOneToAdd=randomIntInRange(generator,whatCanWeAdd.size()-1);
		std::cout<<"We add reaction nr: "<<whatCanWeAdd[whichOneToAdd]<<"from a possible "<<whatCanWeAdd.size()<<"reactions"<<std::endl;
		//+1 required as the file from which we read starts with line 1, but vectors number from 0 
		addThisOne=whatCanWeAdd[whichOneToAdd]+1;
		trialNewCell.push_back(addThisOne);
	}

	else{

		int whichOneToDel=randomIntInRange(generator,trialNewCell.size()-1);
		std::cout<<"We delete nr: "<<availableReactions[whichOneToDel]<<std::endl;
		deleteThisOne=whichOneToDel;
		trialNewCell.erase(trialNewCell.begin()+deleteThisOne);
	}

	
	cell tryIfWorks(trialNewCell);

	double proposedThroughput=tryIfWorks.getPerformance();

	//if the network is not severly paralyzed we implement the changes
	//if (performance*0.5<proposedThroughput){
	//	if(areWeAdding){availableReactions.push_back(addThisOne);}
	//	else{availableReactions.erase(availableReactions.begin()+deleteThisOne);}
	//	performance=proposedThroughput;
	//}
	//else{ std::cout<<"Changes too destructive, not implemented."<<std::endl;}
	
	if (doWeAccept<exp(-1.0*firstPerformance*(performance-proposedThroughput)))
	{
		if(areWeAdding){availableReactions.push_back(addThisOne);}
		else{availableReactions.erase(availableReactions.begin()+deleteThisOne);}
		performance=proposedThroughput;
	}
	else
	{
		std::cout<<"Changes not implemented."<<std::endl;
	}

	std::cout<<"Currently in network: "<<availableReactions.size()<<" reactions, with a fittness of "<<performance<<"."<<std::endl;
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




double cell::calcThroughput(){


	ReactionNetwork graph=allTheReactions;
	std::vector<Vertex> allreacList=reactionVertexList;
	std::vector<Vertex> reacList=subsetVertices(availableReactions,allreacList);

	//need to find a way to only include compounds that are used in current reactions
	//
	
	//plan: create a vector<int> with the same length as substrateVertices that lists which row corresponds to which substrate
	
	//first nrOfInternalMetabolites is the internal metabolites (need to use variable instead of 13 later)
	std::vector<int> substrateIndex(substrateVertexList.size());

	for(int i=0; i<nrOfInternalMetabolites; i++){
		substrateIndex[i]=i+1;
	}

	//counter to keep track of which row of the GLPK problem the next substrate will belong to
	int nextRowNumber=14;
	std::set<Vertex> substrateSet;

	for (int i=0; i<availableReactions.size();i++){

		typename boost::graph_traits<ReactionNetwork>::adjacency_iterator vi, vi_end;

		for(boost::tie(vi,vi_end)=boost::adjacent_vertices(reactionVertexList[availableReactions[i]-1],allTheReactions); vi!=vi_end; ++vi){

			substrateSet.insert(vi.dereference());
		}
	}
	
	//erasing internal metabolites from the set, as we already have those at the beginning of the list
	for(auto metab:internalMetaboliteVList){substrateSet.erase(metab);}
		
	while(!substrateSet.empty()){
		int substrateid=allTheReactions[*substrateSet.begin()].sub.index;

		substrateIndex[substrateid+nrOfInternalMetabolites]=nextRowNumber;
		nextRowNumber++;
		substrateSet.erase(substrateSet.begin());
	}










	glp_prob *lp;
	lp=glp_create_prob();
	glp_set_prob_name(lp,"network_throughput");
	glp_set_obj_dir(lp,GLP_MAX);

	//silencing GLPK output
	glp_term_out(GLP_OFF);

	glp_add_rows(lp,nextRowNumber-1);

	for (int i=1; i<nextRowNumber; i++){
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


		double freeChange=tmpreac.getCurrentFreeEChange();
		//std::cout<<"FreeEChange: "<<freeChange<<std::endl;
		if(freeChange<0){
		glp_set_col_bnds(lp,i,GLP_DB,0.0,1.0);
		}
		else{glp_set_col_bnds(lp,i,GLP_UP,-1.0,0.0);}


		std::vector<int> tmpsubs=tmpreac.getsubstrates();
		std::vector<int> tmpprods=tmpreac.getproducts();

		//ia.push_back(1); ja.push_back(1); ar.push_back(0.0);
		for (int j: tmpsubs){
			//using i+14 as the column numbering starts from 1, and there are 13 internal metabolites
			//with negative substrate indices
			int rownumber=substrateIndex[j+nrOfInternalMetabolites];
			//ia.push_back(j+14);  previous method
			ia.push_back(rownumber);
			ja.push_back(i);
			ar.push_back(-1.0);
		}

		for (int j: tmpprods){
			//using i+14 as the column numbering starts from 1, and there are 13 internal metabolites
			//with negative substrate indices
			int rownumber=substrateIndex[j+nrOfInternalMetabolites];
			//ia.push_back(j+14);  previous method
			ia.push_back(rownumber);
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
	ia.push_back(substrateIndex[908+nrOfInternalMetabolites]);	ja.push_back(listSize+1); ar.push_back(1.0);
	ia.push_back(substrateIndex[-1+nrOfInternalMetabolites]);	ja.push_back(listSize+2); ar.push_back(1.0);
	ia.push_back(substrateIndex[-2+nrOfInternalMetabolites]);	ja.push_back(listSize+3); ar.push_back(-1.0);
	ia.push_back(substrateIndex[nrOfInternalMetabolites]);	ja.push_back(listSize+4); ar.push_back(-1.0);

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


	//for (int i=0; i<=length; i++){


	//	std::cout<<"Matrix element: "<<iarray[i]<<", "<<jarray[i]<<", "<<ararray[i]<<std::endl;
	//}

	//std::cout<<"The length of the vectors are: "<<length<<", "<<ja.size()<<", "<<ar.size()<<std::endl;

	//std::cout<<"Last elements are: "<<iarray[length]<<", "<<jarray[length]<<", "<<ararray[length]<<std::endl;
	
	
	glp_load_matrix(lp,length,iarray,jarray,ararray);

	glp_simplex(lp,NULL);

	double goodness=glp_get_obj_val(lp);

	//std::cout<<"Objective function value: "<<goodness<<std::endl;
	//std::cout<<"Coefficients:";
	//for (int i=0; i<(listSize+4); i++){
	//	double colvalue=glp_get_col_prim(lp,i+1);
	//	std::cout<<colvalue<<", ";
	//}
	//std::cout<<std::endl;

	return goodness-smallKforFitness*availableReactions.size();
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

std::string cell::niceSubstrateName(Vertex currentVertex){

	std::string emptyName=("---");
	std::string tobereturned;
	substrate currentSub=allTheReactions[currentVertex].sub;
	if (currentSub.name.compare(emptyName) ==0)
	{tobereturned=currentSub.molecule;}
	else
	{tobereturned=currentSub.name;}
	return tobereturned;

}

void cell::printHumanReadable(std::vector<Vertex>& substrateList){

	ReactionNetwork graph=allTheReactions;
	std::vector<Vertex> reacList=reactionVertexList;
	std::string emptyName = ("---");

	for (int i: availableReactions){

		int id=graph[reacList[i]].reac.getListNr();

		std::vector<int> products = graph[reacList[i-1]].reac.getproducts();
		std::vector<int> substrates = graph[reacList[i-1]].reac.getsubstrates();

		std::cout<<id<<": ";
		for (int j:substrates){
			//figuring out whether we have a meaningful name
			std::string toPrint,name,molecule;
			molecule=graph[substrateList[j+nrOfInternalMetabolites]].sub.molecule;
			name=graph[substrateList[j+nrOfInternalMetabolites]].sub.name;

			if (name.compare(emptyName) == 0){
				toPrint=molecule;
			}
			else {
				toPrint=name;
			}
			std::cout<<toPrint<<" + ";
			//std::cout<<graph[substrateList[j+13]].sub.index<<" + ";
		}


		std::cout<<" > ";

		for (int j:products){
			//figuring out whether we have a meaningful name
			std::string toPrint,name,molecule;
			molecule=graph[substrateList[j+nrOfInternalMetabolites]].sub.molecule;
			name=graph[substrateList[j+nrOfInternalMetabolites]].sub.name;


			if (name.compare(emptyName) == 0){
				toPrint=molecule;
			}
			else {
				toPrint=name;
			}
			std::cout<<toPrint<<" + ";
			//std::cout<<graph[substrateList[j+13]].sub.index<<" + ";

		}
		std::cout<<std::endl;
	
	}
}
