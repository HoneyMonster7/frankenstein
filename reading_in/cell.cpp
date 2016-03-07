#include "cell.h"
#include <tuple>


// static class variables

std::vector<reaction> cell::reactionVector;
std::vector<substrate> cell::substrateVector;
int cell::nrOfInternalMetabolites;
std::vector<int> cell::sourceSubstrate;
std::vector<int> cell::sinkSubstrate;
double cell::smallKforFitness;
double cell::probabilityOfMutation;
double cell::probabilityOfSinkMutation;
double cell::probabilityOfHorizontalGenetransfer;

cell::cell(std::vector<int>& tmpAvailReacs)

		 :availableReactions(tmpAvailReacs)

{
	double initialPerformance=calcThroughput();
	std::vector<double> fluxThroughReacs(availableReactions.size());
	performance=initialPerformance;
	firstPerformance=initialPerformance/22;
	std::vector<int> additionalSinks;
}

//this constructor just creates the vector availablereacs and initializes the other variables
cell::cell(std::vector<int>& tmpAvailReacs, int notUsed)

		 :availableReactions(tmpAvailReacs)

{
	double initialPerformance=0;
	std::vector<double> fluxThroughReacs(availableReactions.size());
	performance=initialPerformance;
	firstPerformance=initialPerformance/22;
	std::vector<int> additionalSinks;
}

cell::cell() {}

void cell::addThisSink(int thisSinkIsNowYours) {

	additionalSinks.emplace_back(thisSinkIsNowYours);
	performance=calcThroughput();
}

void cell::theseAreYourAddSinks(std::vector<int> theseSinksAreNowYours){

	additionalSinks=theseSinksAreNowYours;
	performance=calcThroughput();
}
void cell::deleteThisSink(int thisSinkIsToBeGone){

	//thisSinkIsToBeGone needs to be the index of the sink to remove
	additionalSinks.erase(additionalSinks.begin()+thisSinkIsToBeGone);
}

void cell::printReacs() {

	std::cout<<"Current reactions: ";
	for (auto i : availableReactions){

		std::cout<<i<<", ";
	}
	std::cout<<" END"<<std::endl;
}

bool cell::operator<(const cell& other) const {

	//orders such that the largest performance wil be [0]
	return getPerformance()>other.getPerformance();
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
	xgmmlFile<<"<graph label=\""<<filename<<"\""<<std::endl<<" xmlns:cy=\"http://www.cytoscape.org\" "<<std::endl<<"xmlns=\"http://www.cs.rpi.edu/XGMML\" "<<std::endl<<"directed=\"1\">"<<std::endl;

	xgmmlFile<<"\t <att name=\"Fittness\" type=\"real\" value=\""<<getPerformance()<<"\"/>"<<std::endl;
//<<"fittnesvalue="<<getPerformance()<<std::endl
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


	for (auto sourceSubstrate_one:sourceSubstrate){
		std::string sourceName=substrateVector[sourceSubstrate_one+nrOfInternalMetabolites].niceSubstrateName();

		compoundIDs.erase(sourceSubstrate_one);
		xgmmlFile<<"<node label=\""<<sourceName<<"\" id=\""<<sourceSubstrate_one+nrOfInternalMetabolites<<"\">"<<std::endl;
		xgmmlFile<<"\t <att name=\"Type\" type=\"string\" value=\"Source\"/>"<<std::endl;
		xgmmlFile<<"</node>"<<std::endl;
	}

	for (auto sinkSubstrate_one:sinkSubstrate){
		std::string sinkName=substrateVector[sinkSubstrate_one+nrOfInternalMetabolites].niceSubstrateName();

		compoundIDs.erase(sinkSubstrate_one);
		xgmmlFile<<"<node label=\""<<sinkName<<"\" id=\""<<sinkSubstrate_one+nrOfInternalMetabolites<<"\">"<<std::endl;
		xgmmlFile<<"\t <att name=\"Type\" type=\"string\" value=\"Sink\"/>"<<std::endl;
		xgmmlFile<<"</node>"<<std::endl;
	}
	for (auto sinkSubstrate_one:additionalSinks){
		std::string sinkName=substrateVector[sinkSubstrate_one+nrOfInternalMetabolites].niceSubstrateName();

		compoundIDs.erase(sinkSubstrate_one);
		xgmmlFile<<"<node label=\""<<sinkName<<"\" id=\""<<sinkSubstrate_one+nrOfInternalMetabolites<<"\">"<<std::endl;
		xgmmlFile<<"\t <att name=\"Type\" type=\"string\" value=\"Sink\"/>"<<std::endl;
		xgmmlFile<<"</node>"<<std::endl;
	}

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
	//	std::cout<<"Changes not implemented."<<std::endl;
	}

	//std::cout<<"Currently in network: "<<availableReactions.size()<<" reactions, with a fittness of "<<performance<<"."<<std::endl;
}


cell cell::mutateAndReturn( RandomGeneratorType& generator ){
	//this function is for the population mutation. in this case there is no acceptance
	//probability, every mutation is accepted, as survival of the fittest is dealt with 
	//in chosing the cells to reproduce


	int deleteThisOne,addThisOne;

	std::vector<int> trialNewCell = availableReactions;

	double addProb=0.5;

	double doWeAdd=randomRealInRange(generator, 1);


	//calculating current throughput
	bool areWeAdding= doWeAdd<=addProb;
	if(areWeAdding){
		std::vector<int> whatCanWeAdd = canBeAdded();
		int whichOneToAdd=randomIntInRange(generator,whatCanWeAdd.size()-1);

		addThisOne=whatCanWeAdd[whichOneToAdd];
		trialNewCell.push_back(addThisOne+1);
	}

	else{

		int whichOneToDel=randomIntInRange(generator,trialNewCell.size()-1);
		//std::cout<<"We delete nr: "<<availableReactions[whichOneToDel]<<std::endl;
		deleteThisOne=whichOneToDel;
		trialNewCell.erase(trialNewCell.begin()+deleteThisOne);
	}

	//not calculating flux just here
	cell tryIfWorks(trialNewCell,0);

	return tryIfWorks;

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


void cell::mutatePopulation(std::vector<int>& population,std::vector<int>& howManyOfEach, std::vector<cell>& cellVector,  RandomGeneratorType& generator){

	bool gotOneToMutate=false;
	double maxPossibleFittness=10;
	int whichOneToMutate;
	bool areWeMutating=false;
	bool areWeHorizontalTransferring=false;
	bool wasThereAChange=false;

	while(!gotOneToMutate){

		//implementing the moran process selection
		whichOneToMutate=cell::randomIntInRange(generator,population.size()-1);
		double compareFittnesWithThis=cell::randomRealInRange(generator,maxPossibleFittness);
		double fittnessOfCurrentCell=cellVector[population[whichOneToMutate]].getPerformance();
		//std::cout<<"Trying if the lucky one is nr "<<whichOneToMutate<<std::endl;
		if(fittnessOfCurrentCell>compareFittnesWithThis){gotOneToMutate=true;}

	}

	int whichCellDies=cell::randomIntInRange(generator,population.size()-1);
	//now figure out if we are mutating
	
	double luckyToMutate=randomRealInRange(generator,1);
	double luckyToHorizontalGene=randomRealInRange(generator,1);

	if (luckyToMutate<cell::probabilityOfMutation) {
		areWeMutating=true;	
	}

	if (luckyToHorizontalGene<cell::probabilityOfHorizontalGenetransfer) {
		areWeHorizontalTransferring=true;
	}


	//now mutate the selected cell, and put it in the place of the dying cell
	
	std::vector<int> additionalSinks=cellVector[population[whichOneToMutate]].getAddSinks();
	std::vector<int> originalReacs=cellVector[population[whichOneToMutate]].getReacs();

	cell tmpCellWithoutPerfCalc(originalReacs,0);

	if (areWeHorizontalTransferring){

		//chosing the cell that transfers a reaction to our current cell
		int whichCellDonatesGenes=randomIntInRange(generator,population.size()-1);

		std::vector<int> GenesOfDonor=cellVector[population[whichCellDonatesGenes]].getReacs();

		int whichReacToDonate=randomIntInRange(generator,GenesOfDonor.size()-1);

		std::vector<int> GenesOfCurrent=tmpCellWithoutPerfCalc.getReacs();

		//converting to unordered set to avoid double reactions
		std::unordered_set<int> setOfGenesOfCurrent(GenesOfCurrent.begin(), GenesOfCurrent.end());

		setOfGenesOfCurrent.insert(GenesOfDonor[whichReacToDonate]);

		//if the donated reaction was already present in the current cell don't bother with creating a new cell
		if (setOfGenesOfCurrent.size() != GenesOfCurrent.size()){
			//std::cout<<"There was a successful horziontal gene transfer."<<std::endl;
			//converting back to vector for initialization of the new cell
			std::vector<int> GenesWithTransferred(setOfGenesOfCurrent.begin(), setOfGenesOfCurrent.end());

			tmpCellWithoutPerfCalc.setReacs(GenesWithTransferred);
			wasThereAChange=true;
		}

		//std::vector<int> genesOfDonor=cellVector[population[
	}

	if (areWeMutating) {
		//creating the mutatnt cell
		cell offspringOfChosenCell=tmpCellWithoutPerfCalc.mutateAndReturn(generator);

		tmpCellWithoutPerfCalc.setReacs(offspringOfChosenCell.getReacs());
		wasThereAChange=true;


	}
	
	if(wasThereAChange){
		//since there was a change now we'll have to recalculate the performance
		//this is done when we re-add the additional sinks (if any)
		tmpCellWithoutPerfCalc.theseAreYourAddSinks(additionalSinks);

		//now we need to find a place to store this mutatnt cell in cellVector
		//
		//decrease the number of those cells that die
		--howManyOfEach[population[whichCellDies]];
		//if the for loop doesn't find any unused places in cellVector, the statement after will exit with an error (this should never happen in normal circumstances)
		int whichIsUnused=population.size();
		for (int i=0; i<howManyOfEach.size();i++){

			if (howManyOfEach[i]==0) {
				whichIsUnused=i;
				break;
			}
		}
		//we have an unused element in cellVector now
		cellVector[whichIsUnused]=tmpCellWithoutPerfCalc;

		//create the population element that points to the newborn cell
		population[whichCellDies]=whichIsUnused;
		//increase the number corresponding to this newborn cell
		++howManyOfEach[whichIsUnused];
	}
	else{
		//if there's only reproduction without mutation just change the pointer of the dying cell to the reproducing one and adjust the numbers in howManyOfEach
		--howManyOfEach[population[whichCellDies]];
		population[whichCellDies]=population[whichOneToMutate];
		++howManyOfEach[population[whichOneToMutate]];

	}

	
}


std::vector<double> cell::getPopulationFittness(std::vector<int>& population, std::vector<cell>& cellVector){

	std::vector<double> toReturn(population.size());

	for (int k=0; k<population.size();k++){

		toReturn[k]=cellVector[population[k]].getPerformance();
	}


	return toReturn;
}

void cell::printPopulationFittnesses(std::vector<int>& population, std::vector<cell>& cellVector){

	std::vector<double> fittnesses=getPopulationFittness(population, cellVector);
	for (auto i:fittnesses){
		std::cout<<i<<", ";
	}
	std::cout<<std::endl;
}

cell cell::printNFittest(std::vector<int>& population, std::vector<cell>& cellVector,int N){


	
	std::sort(population.begin(),population.end(),[&,cellVector,population](int first, int second) {return cellVector[first].getPerformance() > cellVector[second].getPerformance();});


			std::cout<<"The best "<<N<<" are:";
			for (int i=0;i<N;i++){
			std::cout<<cellVector[population[i]].getPerformance()<<", ";
			}
			std::cout<<std::endl;
			return cellVector[population[0]];

}

std::vector<cell> cell::getBestNCells(std::vector<int>& population, std::vector<cell>& cellVector, int N){

	std::vector<cell> toReturn(N);
	//std::vector<double> allFittnesses=getPopulationFittness(population);
	//std::vector<double> sortedFittness(allFittnesses.size());
	//std::partial_sort_copy(allFittnesses.begin(),allFittnesses.end(),sortedFittness.begin(),sortedFittness.end());

	//for (int i=0; i<N; i++){

	//	double currentFittness=*(sortedFittness.end()-i-1);
	//	int position=find(allFittnesses.begin(),allFittnesses.end(),currentFittness)-allFittnesses.begin();
	//	if(position<allFittnesses.size()){
	//		toReturn[i]=population[position];
	//	}
	//	else {
	//		std::cout<<"There was a problem in finding the best "<<N<<" elements of the population."<<std::endl;
	//	}

	//0}
	//doesn't matter if population is sorted even when the Moran process is running, random selection 
	//doesn't care for ordered vector
	std::sort(population.begin(),population.end(),[&,cellVector,population](int first, int second) {return cellVector[first].getPerformance() > cellVector[second].getPerformance();});

	for (int i=0; i<N; i++){

		toReturn[i]=cellVector[population[i]];
		//std::cout<<"Added a cell with fittness: "<<population[i].getPerformance()<<std::endl;
	}
	return toReturn;
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
	for(int sinkSubstrate_one:sinkSubstrate){
		substrateSet.insert(sinkSubstrate_one+nrOfInternalMetabolites);
	}
	for(int sinkSubstrate_one:additionalSinks){
		substrateSet.insert(sinkSubstrate_one+nrOfInternalMetabolites);
	}
	for(int sourceSubstrate_one:sourceSubstrate){
		substrateSet.insert(sourceSubstrate_one+nrOfInternalMetabolites);
	}
		
	while(!substrateSet.empty()){

		substrateIndex[*substrateSet.begin()]=nextRowNumber;
		//std::cout<<nextRowNumber<<" is "<<substrateVector[*substrateSet.begin()].niceSubstrateName()<<std::endl;
		nextRowNumber++;
		substrateSet.erase(substrateSet.begin());
	}


	int nrOfNormalSinks=sinkSubstrate.size();
	int nrOfNormalSource=sourceSubstrate.size();
	int nrOfAdditionalSinks=additionalSinks.size();

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
	glp_add_cols(lp,listSize+4+nrOfNormalSinks+nrOfNormalSource+nrOfAdditionalSinks);

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
		if(freeChange<-10){
		glp_set_col_bnds(lp,i,GLP_DB,0.0,1.0);
		//std::cout<<"going normal"<<std::endl;
		}
		else if(freeChange>10){
			glp_set_col_bnds(lp,i,GLP_DB,-1.0,0.0); 
			//std::cout<<"going backwards"<<std::endl;
		}
		else {
			glp_set_col_bnds(lp,i,GLP_DB,-0.5,0.5);
			//std::cout<<"Going both ways"<<std::endl;
		}


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
	//water
	glp_set_col_bnds(lp,listSize+1,GLP_DB,0.0,40.0);
	//co2
	glp_set_col_bnds(lp,listSize+2,GLP_DB,0.0,40.0);
	//adp->atp auxilliary
	glp_set_col_bnds(lp,listSize+3,GLP_DB,-10.0,40.0);
	//nad_red->nad_ox
	glp_set_col_bnds(lp,listSize+4,GLP_DB,-40.0,40.0);
	//adding/removing phosphate
	glp_set_col_bnds(lp,listSize+4,GLP_DB,-40.0,40.0);


	//add imaginary reaction here:
	//adding water
	ia.push_back(substrateIndex[-1+nrOfInternalMetabolites]);	ja.push_back(listSize+1); ar.push_back(1.0);
	//adding or removing co2
	ia.push_back(substrateIndex[-2+nrOfInternalMetabolites]);	ja.push_back(listSize+2); ar.push_back(-1.0);
	//adding ADP
	ia.push_back(substrateIndex[-6+nrOfInternalMetabolites]);	ja.push_back(listSize+3); ar.push_back(1.0);
	//removing ATP
	ia.push_back(substrateIndex[-7+nrOfInternalMetabolites]);	ja.push_back(listSize+3); ar.push_back(-1.0);
	//adding NadOx
	ia.push_back(substrateIndex[-3+nrOfInternalMetabolites]);	ja.push_back(listSize+4); ar.push_back(1.0);
	//removing Nad_red
	ia.push_back(substrateIndex[-4+nrOfInternalMetabolites]);	ja.push_back(listSize+4); ar.push_back(-1.0);
	//adding phosphate
	ia.push_back(substrateIndex[-8+nrOfInternalMetabolites]);	ja.push_back(listSize+5); ar.push_back(1.0);
	//removing the sink substrate
	for (int sinkIndex=0; sinkIndex<nrOfNormalSinks; ++sinkIndex){
		ia.push_back(substrateIndex[sinkSubstrate[sinkIndex]+nrOfInternalMetabolites]);	ja.push_back(listSize+5+sinkIndex+1); ar.push_back(-1.0);
		//bounding the sinks
		glp_set_col_bnds(lp,listSize+5+sinkIndex+1,GLP_DB,0.0,10.0);
	}
	//the additional sinks (if any)
	for (int sinkIndex=0; sinkIndex<nrOfAdditionalSinks; ++sinkIndex){
		ia.push_back(substrateIndex[additionalSinks[sinkIndex]+nrOfInternalMetabolites]);	ja.push_back(listSize+5+sinkIndex+1+nrOfNormalSinks); ar.push_back(-1.0);
		//bounding the sinks
		glp_set_col_bnds(lp,listSize+5+nrOfNormalSinks+sinkIndex+1,GLP_DB,0.0,10.0);
	}
	//adding source substrates
	for (int sourceIndex=0; sourceIndex<nrOfNormalSource; ++sourceIndex){
		ia.push_back(substrateIndex[sourceSubstrate[sourceIndex]+nrOfInternalMetabolites]);	ja.push_back(listSize+5+nrOfNormalSinks+nrOfAdditionalSinks+sourceIndex+1); ar.push_back(1.0);
		//bounding the sources
		glp_set_col_bnds(lp,listSize+5+nrOfNormalSinks+nrOfAdditionalSinks+sourceIndex+1,GLP_DB,0.0,10.0);
	}

	//ia.push_back(43+14);	ja.push_back(listSize+2); ar.push_back(1.0);
	//ia.push_back(88+14);	ja.push_back(listSize+3); ar.push_back(1.0);
	//
	//target is to maximize the imaginary reactions throughput of the ADP->ATP reaction
	glp_set_obj_coef(lp,listSize+3,1.0);

	//target is to maximize the imaginary reactions throughput of the ADP->ATP reaction
	//glp_set_obj_coef(lp,listSize+3,0.5);
	//glp_set_obj_coef(lp,listSize+5,1);

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

	//WARNING! this is ALL the fluxes, also the ones for the auxilliary reactions
	//care must be taken to extract the fluxes you actually want
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

void cell::findThePaths(std::vector<int> needMore, std::vector<int> needLess, std::vector<int> currentReactions, int TargetCompound, std::vector<reaction>& ReactionVector, std::vector<substrate>& SubstrateVector, std::string fnameString){

	int writeoutcounter=1;
	std::vector<int> doesntHaveToBalance = {-3, -4, -1, -2, 94};

	doesntHaveToBalance.emplace_back(TargetCompound);

	std::set<int> canWeAdd;
	//first figure out what can be added
	for (int i=0; i<currentReactions.size();i++)
	{
		std::vector<int> neighboursOfCurrentReac=ReactionVector[currentReactions[i]].getNeighbours();
		for (int j:neighboursOfCurrentReac){
			canWeAdd.insert(j);
		}
	}

	for (int j:currentReactions){
		canWeAdd.erase(j);
	}

	std::vector<int> PossibleReactionsToAdd(canWeAdd.begin(),canWeAdd.end());

	//std::cout<<"We can add:";
	//for (auto whatever:PossibleReactionsToAdd){
	//	std::cout<<whatever<<", ";
	//}
	//std::cout<<std::endl;

	for (int j:canWeAdd){

		std::vector<int> substrates, products;
		if (reactionVector[j].getCurrentFreeEChange()>0){
			//std::cout<<"For "<<j<<" freeechange is greater than zero."<<std::endl;
			substrates=ReactionVector[j].getproducts();
			products=ReactionVector[j].getsubstrates();
		}
		else{
			substrates=ReactionVector[j].getsubstrates();
			products=ReactionVector[j].getproducts();
		}

		bool doesItSolveANeed=false;
		bool doesItSolveALess=false;

		for (int currProduct:products){
			doesItSolveANeed= doesItSolveANeed || (std::find(needMore.begin(), needMore.end(),currProduct)!=needMore.end());
			if (std::find(needMore.begin(), needMore.end(), currProduct) != needMore.end()){
				doesItSolveANeed=true;
			//	std::cout<<j<<" solves a more."<<std::endl;
			}	
		}
		for (int currSubs:substrates){
			//doesItSolveALess=doesItSolveALess || (std::find(needLess.begin(), needLess.end(), currSubs) != needLess.end());
		
			if (std::find(needLess.begin(), needLess.end(), currSubs) != needLess.end()){
				doesItSolveALess=true;
				//std::cout<<j<<" solves a less."<<std::endl;
			}
		}

		//now if adding the reaction solves a need or a too much problem, we add it

		//std::cout<<"solveless: "<<doesItSolveALess<<" solvemore: "<<doesItSolveANeed<<std::endl;
		if (doesItSolveALess || doesItSolveANeed){

			std::vector<int> currentReactionsInLoop=currentReactions;
			std::vector<int> needMoreInLoop=needMore;
			std::vector<int> needLessInLoop=needLess;

			//std::cout<<"Needmore: ";
			//for (int m:needMoreInLoop){
			//	std::cout<<m<<", ";
			//}
			//std::cout<<std::endl;
			//
			//std::cout<<"Needless: ";
			//for (int m:needLessInLoop){
			//	std::cout<<m<<", ";
			//}
			//std::cout<<std::endl;



			currentReactionsInLoop.emplace_back(j);
			// and we note the problems it creates, and solves
	
			for (int noBalance:doesntHaveToBalance){

				substrates.erase(std::remove(substrates.begin(), substrates.end(), noBalance), substrates.end());
				products.erase(std::remove(products.begin(), products.end(), noBalance), products.end());
			}	

			for (int k:substrates){

				if (std::find(needLessInLoop.begin(), needLessInLoop.end(), k)!=needLessInLoop.end()){
					needLessInLoop.erase(std::remove(needLessInLoop.begin(), needLessInLoop.end(), k), needLessInLoop.end());
					//std::cout<<"We no longer need less "<<k<<" thanks to "<<j<<std::endl;
				}
				else{
				needMoreInLoop.emplace_back(k);
					//std::cout<<"We now need more "<<k<<" thanks to "<<j<<std::endl;
				}
			}
			for (int k:products){
				if (std::find(needMoreInLoop.begin(), needMoreInLoop.end(), k) != needMoreInLoop.end()){
					needMoreInLoop.erase(std::remove(needMoreInLoop.begin(), needMoreInLoop.end(), k), needMoreInLoop.end());
					//std::cout<<"We no longer need more "<<k<<" thanks to "<<j<<std::endl;
				}
				else{needLessInLoop.emplace_back(k);
				//std::cout<<"We now need less "<<k<<" thanks to "<<j<<std::endl;
				}
			

				}

			//now writing out the current state, and calling it again

			int thingsInNeed=needMoreInLoop.size()+needLessInLoop.size();

			if(thingsInNeed<2){
				std::cout<<std::endl<<"Current state of the system:"<<std::endl;

				std::cout<<"Reactions: ";
				for (int k:currentReactionsInLoop){std::cout<<k<<", ";}
				std::cout<<std::endl<<"Needmore: ";
				for (int k:needMoreInLoop){std::cout<<k<<", ";}
				std::cout<<std::endl<<"Needless: ";
				for (int k:needLessInLoop){std::cout<<k<<", ";}
				std::cout<<std::endl;
				if (thingsInNeed<1){

					std::ostringstream forFileName;
					forFileName<<fnameString<<"/NR"<<writeoutcounter<<"cell";


					cell tmpcell(currentReactionsInLoop);
					tmpcell.printXGMML(forFileName.str());

					writeoutcounter++;
				}
			}

			if (currentReactionsInLoop.size()<7){
				cell::findThePaths(needMoreInLoop,needLessInLoop, currentReactionsInLoop, TargetCompound, ReactionVector, SubstrateVector, fnameString);
			}
		}
	}


}



void cell::printProgressFile(std::vector<int>& population, std::vector<cell>& cellVector, std::vector<int>& howManyOfEach,int k, int outerloop,const int generationsPerWriteout, int checkPointLength, std::ofstream& fileToWrite,double& previousAvgFittness, double maxFittQueue [], double avgFittQueue [], double entropyQueue [], int bestNetSizeQueue [], double avgNetSizeQueue [], int bestUsedReacsQueue [], double avgUsedReacsQueue []){

	double maxFittness=0;
	double totFittness=0;
	int bestNetworkSize;
	int totalNetworkSize=0;
	double enthropy=0;
	int bestUsedReacs=0;
	int totUsedReacs=0;

	//now runnign through the population to find the desired quantities
	//only running through the vector containing how many of each cells there are in the population
	//in a population containing many similar cells this can be much faster than running through each of the cells (in the population vector)
	for (int i = 0; i < population.size(); ++i) {
		int element=howManyOfEach[i];	
		if(element!=0){
			enthropy+=element*std::log(element);

			double currFittness=cellVector[i].getPerformance();
			int currReactionSize=cellVector[i].getReacs().size();

			int fluxcounter=0;
			std::vector<double> currentFluxes=cellVector[i].getFluxes();
			//the resizing is to make sure that we only count the nonzero real reactions not the auxilliary ones
			currentFluxes.resize(currReactionSize);
			for (double inspectedFlux:currentFluxes)
			{
				if (std::abs(inspectedFlux)>1e-5) {
					++fluxcounter;
				}
			}

			totUsedReacs+=fluxcounter*element;

			totFittness+=currFittness*element;
			totalNetworkSize+=currReactionSize*element;

			if (maxFittness<currFittness) {
				
				maxFittness=currFittness;
				bestNetworkSize=currReactionSize;
				bestUsedReacs=fluxcounter;
			}

		}
	}
	int remainderOfK=k%generationsPerWriteout;
	//now writing the required stuff in the file
	if (remainderOfK==0){

		double avgFittness=totFittness/(double)cellVector.size();
		//was used for testing
		//fileToWrite<<"Current avgfittness is "<<avgFittness<<" previousAvgFittness is "<<previousAvgFittness<<" difference is "<<avgFittness-previousAvgFittness<<std::endl;
		//if there was a large enough improvement since the last writeout, write out the steps inbetween
		if (avgFittness>0.3+previousAvgFittness){

			for (int i=1; i<generationsPerWriteout; ++i){

				long long int number=k-generationsPerWriteout+i+outerloop*checkPointLength;
				fileToWrite<<number<<" "<<maxFittQueue[i]<<" "<<entropyQueue[i]<<" "<<avgFittQueue[i]<<" "<<bestNetSizeQueue[i]<<" "<<avgNetSizeQueue[i]<<" "<<bestUsedReacsQueue[i]<<" "<<avgUsedReacsQueue[i]<<std::endl;

			}
		}
	
		long long int generationNr=k+outerloop*checkPointLength;
		fileToWrite<<generationNr<<" "<<maxFittness<<" "<<-1*enthropy<<" "<<avgFittness<<" "<<bestNetworkSize<<" "<<totalNetworkSize/(double)cellVector.size()<<" "<<bestUsedReacs<<" "<<totUsedReacs/(double)cellVector.size()<<std::endl;
		
		//setting the previous avg network fittness to the current value
		previousAvgFittness=avgFittness;
	}

	//Popoulating the next positions of the queue
	maxFittQueue[remainderOfK]=maxFittness;
	avgFittQueue[remainderOfK]=totFittness/cellVector.size();
	entropyQueue[remainderOfK]=-1*enthropy;
	bestNetSizeQueue[remainderOfK]=bestNetworkSize;
	avgNetSizeQueue[remainderOfK]=totalNetworkSize/(double)cellVector.size();
	bestUsedReacsQueue[remainderOfK]=bestUsedReacs;
	avgUsedReacsQueue[remainderOfK]=totUsedReacs/(double)cellVector.size();
}

void cell::setReacs(std::vector<int> theseAreYourNewReacs){

	//Note: no performance recalc was done here. We do that after all the mutations are done when re-adding the additional sinks
	availableReactions=theseAreYourNewReacs;
}
	


