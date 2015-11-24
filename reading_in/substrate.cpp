#include "reaction.h"

substrate::substrate(int tmpindex, double tmpfreeOfCreation, std::string tmpmolecule, std::string tmpname, int tmpcharge)
   : index(tmpindex)
   , freeOfCreation(tmpfreeOfCreation)
   , molecule(tmpmolecule)
   , name (tmpname)
   , charge(tmpcharge)

{std::vector<int> involvedInReacs;}

substrate::substrate() {}


void substrate::addInvolved(int addThisReac){
	int tmpadd=addThisReac;
	//involvedInReacs.push_back(tmpadd);
	involvedInReacs.emplace_back(tmpadd);
}

void substrate::printInvolved(){
	std::cout<<niceSubstrateName()<<" is involved in: ";
	for (int i: involvedInReacs){
		std::cout<<i<<", ";
	}

	std::cout<<std::endl;
}

void substrate::buildNeighbourList(std::vector<reaction>& reacList, std::vector<substrate>& substrateVector){


	for (reaction& currentReac : reacList){

		std::set<int> neighbourSet;

		std::vector<int> currSubs=currentReac.getsubstrates();
		std::vector<int> currProds=currentReac.getproducts();

		//std::cout<<currentReac.getListNr()<<" has neighbours: ";


		for (int i:currSubs)
		{
			std::vector<int> involvedIn=substrateVector[i+reaction::nrOfInternalMetabolites].getinvolved();
			for (int j:involvedIn){
				neighbourSet.insert(j);
				//std::cout<<j<<", ";
			}
		}
		for (int i:currProds)
		{
			std::vector<int> involvedIn=substrateVector[i+reaction::nrOfInternalMetabolites].getinvolved();
			for (int j:involvedIn){
				neighbourSet.insert(j);
				//std::cout<<j<<", ";
			}
		}

		//std::cout<<std::endl;

		std::vector<int> toBeNeighbours(neighbourSet.begin(), neighbourSet.end());
		currentReac.setNeighbours(toBeNeighbours);
		//int currnr=currentReac.getListNr();
		//std::cout<<currnr<<" has "<<toBeNeighbours.size()<<" neighbours, set has "<<neighbourSet.size()<<std::endl;
	}

}

std::string substrate::niceSubstrateName(){

		std::string emptyName=("---");
		std::string tobereturned;
		if (name.compare(emptyName) ==0)
			{tobereturned=molecule;}
		else
			{tobereturned=name;}
		return tobereturned;

}
