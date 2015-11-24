#include "cell.h"
#include <tuple>


// static class variables

std::vector<reaction> cell::reactionVector;
std::vector<substrate> cell::substrateVector;
int cell::nrOfInternalMetabolites;
int cell::sourceSubstrate;
int cell::sinkSubstrate;
double cell::smallKforFitness;

cell::cell(std::vector<int>& tmpAvailReacs)

		 :availableReactions(tmpAvailReacs)

{
	double initialPerformance=calcThroughput();
	std::vector<double> fluxThroughReacs(availableReactions.size());
	performance=initialPerformance;
	firstPerformance=initialPerformance/22;
}


cell::cell() {}


void cell::printReacs() {

	std::cout<<"Current reactions: ";
	for (auto i : availableReactions){

		std::cout<<i<<", ";
	}
	std::cout<<" END"<<std::endl;
}

void cell::printXGMML(std::string filename){


	std::ofstream outfile,typesfile,edgeFile,xgmmlFile;
	//for the node_types
	std::set<int> reacNumbers;
	std::set<int> compoundIDs;
	std::set<int> internalMetIDs;
	std::vector<std::tuple<int,double,int>> edgeVector;

	std::string fileToOpen=filename + ".xgmml";
	xgmmlFile.open(fileToOpen);


	//typesfile<<"Name Type"<<std::endl;
	xgmmlFile<<"<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>"<<std::endl;
	xgmmlFile<<"<graph label=\"justAGraph\" xmlns:cy=\"http://www.cytoscape.org\" xmlns=\"http://www.cs.rpi.edu/XGMML\" directed=\"1\">"<<std::endl;


	for (int j=0; j<availableReactions.size();j++){

		int i=availableReactions[j];
		double fluxOfCurrentReaction=fluxThroughReacs[j];
		//in order to get rid of the 1e-300 kind of fluxes
		if (std::abs(fluxOfCurrentReaction)<1e-6){fluxOfCurrentReaction=0;}
		reaction currentReac= reactionVector[i-1];
		std::vector<int> currentsubs=currentReac.getsubstrates();
		std::vector<int> currentproducts=currentReac.getproducts();
		int reacNR=currentReac.getListNr();

		//adding the current reaction to a set, will be used to generate the node-type file
		reacNumbers.insert(reacNR);

		for (int sub:currentsubs){
				std::string substrateName=substrateVector[sub+nrOfInternalMetabolites].niceSubstrateName();
					compoundIDs.insert(sub);
					//-1 here because that will be the ID of the reactions
					edgeVector.emplace_back(sub+nrOfInternalMetabolites,fluxOfCurrentReaction,-1*reacNR);

		}

		for (int prod:currentproducts){
				std::string productName=substrateVector[prod+nrOfInternalMetabolites].niceSubstrateName();

					compoundIDs.insert(prod);
					//-1 here because that will be the ID of the reactions
					edgeVector.emplace_back(-1*reacNR,fluxOfCurrentReaction,prod+nrOfInternalMetabolites);

		}
	}

	for (int i=0; i<nrOfInternalMetabolites; i++){
		internalMetIDs.insert(i-nrOfInternalMetabolites);
	}


	std::string sourceName=substrateVector[sourceSubstrate+nrOfInternalMetabolites].niceSubstrateName();
	std::string sinkName=substrateVector[sinkSubstrate+nrOfInternalMetabolites].niceSubstrateName();

	compoundIDs.erase(sourceSubstrate);
	compoundIDs.erase(sinkSubstrate);
	xgmmlFile<<"<node label=\""<<sourceName<<"\" id=\""<<sourceSubstrate+nrOfInternalMetabolites<<"\">"<<std::endl;
	xgmmlFile<<"\t <att name=\"Type\" type=\"string\" value=\"Source\"/>"<<std::endl;
	xgmmlFile<<"</node>"<<std::endl;
	xgmmlFile<<"<node label=\""<<sinkName<<"\" id=\""<<sinkSubstrate+nrOfInternalMetabolites<<"\">"<<std::endl;
	xgmmlFile<<"\t <att name=\"Type\" type=\"string\" value=\"Sink\"/>"<<std::endl;
	xgmmlFile<<"</node>"<<std::endl;


	//looping through all the internalMet names writing them into the type file, removing them from the substrate list
	while(!internalMetIDs.empty()){
		compoundIDs.erase(*internalMetIDs.begin());
		std::string currentName=substrateVector[*internalMetIDs.begin()+nrOfInternalMetabolites].niceSubstrateName();

		xgmmlFile<<"<node label=\""<<currentName<<"\" id=\""<<*internalMetIDs.begin()+nrOfInternalMetabolites<<"\">"<<std::endl;
		xgmmlFile<<"\t <att name=\"Type\" type=\"string\" value=\"InternalMet\"/>"<<std::endl;
		xgmmlFile<<"</node>"<<std::endl;
		internalMetIDs.erase(internalMetIDs.begin());
	}
	//doing the same with reacnubmers
	while(!reacNumbers.empty()){
		int currentReacNumber=*reacNumbers.begin();
		xgmmlFile<<"<node label=\""<<currentReacNumber<<"\" id=\""<<-1*currentReacNumber<<"\">"<<std::endl;
		xgmmlFile<<"\t <att name=\"Type\" type=\"string\" value=\"Reaction\"/>"<<std::endl;
		xgmmlFile<<"</node>"<<std::endl;
		reacNumbers.erase(reacNumbers.begin());
	}
	//now with the normal substrates NEED TO DIFFERENTIATE BETWEEN SOURCE AND SINK LATER
	

	while(!compoundIDs.empty()){
		std::string currentName=substrateVector[*compoundIDs.begin()+nrOfInternalMetabolites].niceSubstrateName();
		xgmmlFile<<"<node label=\""<<currentName<<"\" id=\""<<*compoundIDs.begin()+nrOfInternalMetabolites<<"\">"<<std::endl;
		xgmmlFile<<"\t <att name=\"Type\" type=\"string\" value=\"Compound\"/>"<<std::endl;
		xgmmlFile<<"</node>"<<std::endl;
		compoundIDs.erase(compoundIDs.begin());
	}

	for (std::tuple<int,double,int> currentEdge:edgeVector){
		int source, sink;
		double flux;
		std::tie (source,flux,sink) = currentEdge;

		xgmmlFile<<"<edge label=\""<<source<<" - "<<sink<<"\" source=\""<<source<<"\" target=\""<<sink<<"\">"<<std::endl;
		xgmmlFile<<"\t <att name=\"flux\" type=\"real\" value=\""<<flux<<"\"/>"<<std::endl;
		xgmmlFile<<"</edge>"<<std::endl;

	}

	xgmmlFile<<"</graph>"<<std::endl;
	
	xgmmlFile.close();
	}






std::vector<int> cell::canBeAdded(){



	int nrOfAvailableReactions=availableReactions.size();
	
	std::set<int> canAdd;

	for (int i=0; i<nrOfAvailableReactions;i++){
	//	std::cout<<"Reaction number:"<<i<<std::endl;

		std::vector<int> neighboursOfCurrentReac=reactionVector[availableReactions[i]-1].getNeighbours();
		//reactionVector[availableReactions[i]-1].printReaction();
		//reactionVector[availableReactions[i]-1].printNeighbours();

		for (int j:neighboursOfCurrentReac){
			canAdd.insert(j);
		}
	}


	for (int i:availableReactions){

		//getting rid of reactions that are already in the network
	
		int alreadyInReacNr=reactionVector[i-1].getListNr();
		canAdd.erase(alreadyInReacNr);
	}

	std::vector<int> tobereturned(canAdd.begin(),canAdd.end());
	
	//uncomment for printing out the line numbers of possible reactions to be added 
	//std::cout<<"We can add: ";
	//for (auto whatever:tobereturned){
	//	std::cout<<whatever<<", ";
	//}
	//std::cout<<std::endl;
		return tobereturned;

}

void cell::mutate( RandomGeneratorType& generator ){


	//int compoundSize= substrateVertexList.size();
//	std::cout<<"Random numbers:";
//	for (int i=1; i<50; i++){
//		std::cout<<gen()%availableReactions.size()<<", ";
//	}

	//int currentReactionNumber = availableReactions.size();


	int deleteThisOne,addThisOne;

	std::vector<int> trialNewCell = availableReactions;

	//double addProb= 0.1*exp(desiredNetworkSize-currentReactionNumber);
	//double delProb=0.1*exp(currentReactionNumber-desiredNetworkSize);

	double addProb=0.5;

	double doWeAdd=randomRealInRange(generator, 1);
	//double doWeDelete=randomRealInRange(generator,1);
	double doWeAccept=randomRealInRange(generator,1);


	//calculating current throughput
	bool areWeAdding= doWeAdd<=addProb;
	if(areWeAdding){
		std::vector<int> whatCanWeAdd = canBeAdded();
		int whichOneToAdd=randomIntInRange(generator,whatCanWeAdd.size()-1);
		//std::cout<<"We add reaction nr: "<<whatCanWeAdd[whichOneToAdd]<<"from a possible "<<whatCanWeAdd.size()<<"reactions"<<std::endl;
		//+1 required as the file from which we read starts with line 1, but vectors number from 0 
		//it is required, as the canbeAdded returns positions in the reactionVector and the
		//availableReactions needs positions in the original file (vector position +1)
		addThisOne=whatCanWeAdd[whichOneToAdd];
		trialNewCell.push_back(addThisOne+1);
	}

	else{

		int whichOneToDel=randomIntInRange(generator,trialNewCell.size()-1);
		//std::cout<<"We delete nr: "<<availableReactions[whichOneToDel]<<std::endl;
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
	
	if ((doWeAccept<exp(-1.0*(performance-proposedThroughput)/firstPerformance) && (availableReactions.size()>2)))
	{
		if(areWeAdding){availableReactions.push_back(addThisOne+1);}
		else{availableReactions.erase(availableReactions.begin()+deleteThisOne);}
		performance=proposedThroughput;
		std::vector<double> tmpFluxes=tryIfWorks.getFluxes();
		setFluxes(tmpFluxes);
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



	//need to find a way to only include compounds that are used in current reactions
	//
	
	//plan: create a vector<int> with the same length as substrateVector that lists which row corresponds to which substrate
	
	//first nrOfInternalMetabolites is the internal metabolites (need to use variable instead of 13 later)
	std::vector<int> substrateIndex(substrateVector.size());

	for(int i=0; i<nrOfInternalMetabolites; i++){
		substrateIndex[i]=i+1;
	}

	//counter to keep track of which row of the GLPK problem the next substrate will belong to
	int nextRowNumber=14;
	std::set<int> substrateSet;

	for (int i=0; i<availableReactions.size();i++){


		//CHECK IF THERE IS A +1 REQUIRED
		reaction currentReac=reactionVector[availableReactions[i]-1];

		//std::cout<<"Reaction "<<currentReac.getListNr()<<" has ";
		for (int j:currentReac.getsubstrates()){
			substrateSet.insert(j+nrOfInternalMetabolites);
		//	std::cout<<substrateVector[j+nrOfInternalMetabolites].niceSubstrateName()<<", ";
		}
		for (int j:currentReac.getproducts()){
			substrateSet.insert(j+nrOfInternalMetabolites);
		//	std::cout<<substrateVector[j+nrOfInternalMetabolites].niceSubstrateName()<<", ";
		}

		//std::cout<<std::endl;
	}
	
	
	//erasing internal metabolites from the set, as we already have those at the beginning of the list
	for(int metab=0; metab<nrOfInternalMetabolites; metab++){substrateSet.erase(metab);}

	//in order to always have the source and sink nodes
	substrateSet.insert(sinkSubstrate+nrOfInternalMetabolites);
	substrateSet.insert(sourceSubstrate+nrOfInternalMetabolites);
		
	while(!substrateSet.empty()){

		substrateIndex[*substrateSet.begin()]=nextRowNumber;
		//std::cout<<nextRowNumber<<" is "<<substrateVector[*substrateSet.begin()].niceSubstrateName()<<std::endl;
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
	int listSize=availableReactions.size();
	glp_add_cols(lp,listSize+5);

	//for(int i=1; i<=listSize+4; i++){ glp_set_col_bnds(lp,i,GLP_LO,0.0,0.0);}

	std::vector<int> ia,ja;
	std::vector<double> ar;

	for (int i=1;i<=listSize;i++){


		//preparing the sparse matrix's values 
		//final -1 because the availableReactions stores the id form the file
		//starting with 1, as opposed to the vector starting at 0
		reaction tmpreac=reactionVector[availableReactions[i-1]-1];


		double freeChange=tmpreac.getCurrentFreeEChange();
		//std::cout<<"FreeEChange of : "<<tmpreac.getListNr()<<" is "<<freeChange<<std::endl;
		//tmpreac.printReaction();
		if(freeChange<0){
		glp_set_col_bnds(lp,i,GLP_DB,0.0,1.0);
		}
		else{glp_set_col_bnds(lp,i,GLP_DB,-1.0,0.0);}


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
	glp_set_col_bnds(lp,listSize+5,GLP_DB,0.0,10.0);

	glp_set_col_bnds(lp,listSize+4,GLP_DB,-10.0,10.0);
	//add imaginary reaction here:
	ia.push_back(substrateIndex[sourceSubstrate+nrOfInternalMetabolites]);	ja.push_back(listSize+1); ar.push_back(1.0);
	ia.push_back(substrateIndex[-1+nrOfInternalMetabolites]);	ja.push_back(listSize+2); ar.push_back(1.0);
	ia.push_back(substrateIndex[-2+nrOfInternalMetabolites]);	ja.push_back(listSize+3); ar.push_back(-1.0);
	ia.push_back(substrateIndex[-8+nrOfInternalMetabolites]);	ja.push_back(listSize+5); ar.push_back(-1.0);
	ia.push_back(substrateIndex[sinkSubstrate+nrOfInternalMetabolites]);	ja.push_back(listSize+4); ar.push_back(-1.0);

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


	//for (int i=1; i<=length+1; i++){
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

	std::vector<double> tmpFluxes(availableReactions.size());

	//std::cout<<"Fluxes are: ";
	for (int i=0;i<availableReactions.size();i++){

		tmpFluxes[i]=glp_get_col_prim(lp,i+1);
	//	std::cout<<tmpFluxes[i]<<", ";
	}
	//std::cout<<std::endl;

	setFluxes(tmpFluxes);
	//std::cout<<"goodness is "<<goodness<<std::endl;;
	//std::cout<<"Fittness is "<<goodness-smallKforFitness*availableReactions.size()<<std::endl;;
	glp_delete_prob(lp);

	return goodness-smallKforFitness*availableReactions.size();
}


void cell::printHumanReadable(){


	for (int i: availableReactions){

		int id=reactionVector[i].getListNr();

		std::vector<int> products = reactionVector[i].getproducts();
		std::vector<int> substrates = reactionVector[i].getsubstrates();

		std::cout<<id<<": ";
		for (int j:substrates){
			//figuring out whether we have a meaningful name
			std::string toPrint=substrateVector[j].niceSubstrateName();

			std::cout<<toPrint<<" + ";
			//std::cout<<graph[substrateList[j+13]].sub.index<<" + ";
		}


		std::cout<<" > ";

		for (int j:products){
			//figuring out whether we have a meaningful name

			std::string toPrint=substrateVector[j].niceSubstrateName();

			std::cout<<toPrint<<" + ";
			//std::cout<<graph[substrateList[j+13]].sub.index<<" + ";

		}
		std::cout<<std::endl;
	
	}
}
void cell::setFluxes(std::vector<double>& fluxVector){

	fluxThroughReacs=fluxVector;
}
