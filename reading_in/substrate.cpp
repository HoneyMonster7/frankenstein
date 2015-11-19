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
	involvedInReacs.emplace_back(addThisReac);
}

void substrate::buildNeighbourList(std::vector<reaction>& reacList, std::vector<substrate>& substrateVector){


	for (reaction& currentReac : reacList){

		std::set<int> neighbourSet;

		std::vector<int> currSubs=currentReac.getsubstrates();
		std::vector<int> currProds=currentReac.getproducts();

		for (int i:currSubs)
		{
			std::vector<int> involvedIn=substrateVector[i+reaction::nrOfInternalMetabolites].getinvolved();
			for (int j:involvedIn){neighbourSet.insert(j);}
		}
		for (int i:currProds)
		{
			std::vector<int> involvedIn=substrateVector[i+reaction::nrOfInternalMetabolites].getinvolved();
			for (int j:involvedIn){neighbourSet.insert(j);}
		}

		std::vector<int> toBeNeighbours(neighbourSet.begin(), neighbourSet.end());
		currentReac.setNeighbours(toBeNeighbours);
		int currnr=currentReac.getListNr();
		std::cout<<currnr<<" has "<<toBeNeighbours.size()<<" neighbours"<<std::endl;
	}

}
