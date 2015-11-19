#include "reaction.h"

int reaction::nrOfInternalMetabolites;


reaction::reaction(int reacNr,double tmpfreechange, const std::vector<int>& tmpsubstrates, const std::vector<int>& tmproducts, const InternalMetsT& tmpInternalMets)
	: listNR(reacNr)
	, substrates(tmpsubstrates)
	, products(tmproducts)
	, internalMets(tmpInternalMets)
	, freeEChange(tmpfreechange)
	  
{}

reaction::reaction() {}

void reaction::printReaction()
{
	std::cout<<"From: ";
	for(auto i =substrates.begin(); i!=substrates.end(); ++i) std::cout<<*i<<' ';
	std::cout<<std::endl;
	std::cout<<"To: ";
	for(auto i=products.begin(); i!=products.end(); ++i) std::cout<<*i<<' ';
	std::cout<<std::endl;
	std::cout<<"Free energy change: "<<freeEChange<<std::endl;
}

//possibly move this function to the substrate class
void reaction::readCompounds(std::string fileName,  std::vector<substrate>& substrateVector)
{
	std::ifstream inFile(fileName);
	std::string line, nameChem, namenorm;

	if(!inFile)
	{
		throw std::runtime_error("Can't open input file " + fileName);
	}

	while(std::getline(inFile, line))
	{

		int tmpindex, tmpcharge;
		double tmpfreeEofcreation;
		std::string tmpname, tmpmolecule;
		std::stringstream iss(line);
		iss>>tmpindex;
		iss>>tmpfreeEofcreation;
		iss>>tmpmolecule;
		iss>>tmpname;
		iss>>tmpcharge;

		//creating substrate and putting it into the vector
		substrate tmpsubstrate(tmpindex,tmpfreeEofcreation,tmpmolecule,tmpname, tmpcharge);

		substrateVector.emplace_back(tmpsubstrate);

	}
}

int reaction::getListNr(){ return listNR;}

void reaction::readReactions(std::string fileName, std::vector<reaction>& reacPointer, std::vector<substrate>& substrateVector)
{
	std::ifstream inFile(fileName);
	std::string line,tmpsubs,tmpprods;


	double tmpfreeE;


	if(!inFile)
	{
		throw std::runtime_error("Can't open input file " + fileName);
	}

	int counter=0;
	while(std::getline(inFile, line))
	{
	InternalMetsT tmpinternalMets = {};
		std::vector<int> tmpsubstrates, tmproducts;
		std::stringstream iss(line);
		iss>>tmpfreeE;
		std::getline(iss,tmpsubs, '>');
		std::getline(iss,tmpprods);

		std::istringstream subss(tmpsubs);
		std::istringstream prodss(tmpprods);

		for(int tmp; subss>>tmp;)
		{
			tmpsubstrates.push_back(tmp);
		}

		for(int tmp; prodss>>tmp;)
		{
			tmproducts.push_back(tmp);
		}

		//finding if any of the internal metabolites appear on any side of the reaction
		for(int i=-nrOfInternalMetabolites; i<0; i++)
		{
			if(std::find(tmpsubstrates.begin(), tmpsubstrates.end(), i) != tmpsubstrates.end())
			{
				tmpinternalMets[i+nrOfInternalMetabolites]--;
				//std::cout<<"Internal metabolite nr "<<i<<" found going in."<<std::endl;
			}

			if(std::find(tmproducts.begin(), tmproducts.end(), i) != tmproducts.end())
			{
				tmpinternalMets[i+nrOfInternalMetabolites]++;
			}
		}


		reacPointer.emplace_back(counter,tmpfreeE, tmpsubstrates, tmproducts, tmpinternalMets);




		for(int i : tmpsubstrates)
		{
		if(i>=0){
			substrateVector[i+nrOfInternalMetabolites].addInvolved(i);	
		}
		}


		for(int i: tmproducts)
		{
		if(i>=0){
			substrateVector[i+nrOfInternalMetabolites].addInvolved(i);	
		}
		}

		counter++;
	}
}

double reaction::freeEchange() { return freeEChange; }


void reaction::recalcEchange(const environment& env)
{
	//recalculate the freeEchang using the formula G=G0+RT*ln(([C]^c*[D]^d)/([A]^a*[B]^b))
	//need the reactions database changed for this

	double insideLog=(std::pow(env.nh2acceptorCont,internalMets[0])*std::pow(env.nh2donorCont,internalMets[1])*std::pow(env.conh22Cont,internalMets[2])*std::pow(env.nh3aqCont,internalMets[3])*std::pow(env.ppiCont,internalMets[4])*std::pow(env.piCont,internalMets[5])*std::pow(env.atpCont,internalMets[6])*std::pow(env.adpCont,internalMets[7])*std::pow(env.ampCont,internalMets[8])*std::pow(env.nadredcont,internalMets[9])*std::pow(env.nadoxcont,internalMets[10])*std::pow(env.co2cont,internalMets[11])*std::pow(env.h2ocont,internalMets[12]));

	currentFreeEChange=freeEChange+8.3144598e-3*env.temperature*std::log(insideLog);
	//std::cout<<"Inside the log: "<<insideLog<<"freechange now: "<<currentFreeEChange<<std::endl;
}

void reaction::setNeighbours(std::vector<int>& neighbourList){
neighbourOf=neighbourList;
}

void reaction::printNeighbours(){
	for (const int & i:neighbourOf){
		std::cout<<i<<", ";
	}
	std::cout<<std::endl;
}

