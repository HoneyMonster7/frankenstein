#include "reaction.h"


void reaction::readCompounds(std::string fileName, ReactionNetwork& graph, std::vector<Vertex>& vertexList)
{
	std::ifstream inFile(fileName);
	std::string line, nameChem, namenorm;
	int index, charge;
	double formFreeE;

	while (std::getline(inFile, line))
	{
		std::stringstream iss(line);
		iss>>index;
		iss>>formFreeE;
		iss>>nameChem;
		iss>>namenorm;
		iss>>charge;
		vertexList.emplace_back(boost::add_vertex(graph));
	}
}

void reaction::readReactions(std::string fileName, std::vector<reaction>& reacPointer, ReactionNetwork& graph, std::vector<Vertex>& vertexList, const std::vector<Vertex>& compoundVList)
{
	std::ifstream inFile(fileName);
	std::string line,tmpsubs,tmpprods;

	int tmpinternalMets[9] = {};

	double tmpfreeE;

	typedef boost::graph_traits<ReactionNetwork>::edge_descriptor Edge;

	while(std::getline(inFile, line))
	{
		std::vector<int> tmpsubstrates, tmproducts;
		std::stringstream iss(line);
		iss>>tmpfreeE;
		std::getline(iss,tmpsubs, '>');
		std::getline(iss,tmpprods);

		std::istringstream subss(tmpsubs);
		std::istringstream prodss(tmpprods);

		for (int tmp; subss>>tmp;)
		{
			tmpsubstrates.push_back(tmp);
		}

		for (int tmp; prodss>>tmp;)
		{
			tmproducts.push_back(tmp);
		}

		//finding if any of the internal metabolites appear on any side of the reaction
		for (int i=-9; i<0; i++)
		{
			if(std::find(tmpsubstrates.begin(), tmpsubstrates.end(), i) != tmpsubstrates.end())
			{
				tmpinternalMets[i+9]--;
			}

			if(std::find(tmproducts.begin(), tmproducts.end(), i) != tmproducts.end())
			{
				tmpinternalMets[i+9]++;
			}
		}


		reacPointer.emplace_back(tmpfreeE, tmpsubstrates, tmproducts, tmpinternalMets);


		vertexList.emplace_back(boost::add_vertex(graph));
		graph[vertexList[vertexList.size()-1]].reac=reaction(tmpfreeE,tmpsubstrates,tmproducts,tmpinternalMets);


		for (int i : tmpsubstrates)
		{
			Edge e1;
			e1=(boost::add_edge(compoundVList[i+9],vertexList[vertexList.size()-1],graph)).first;
		}


		for (int i: tmproducts)
		{
			Edge e1;
			e1=(boost::add_edge(vertexList[vertexList.size()-1],compoundVList[i+9],graph)).first;
		}

		//		tmpSR.reac=new reaction(tmpType,tmpsubI,tmpProdI,tmpnrATP,tmpnrNADH,tmpfreeE,tmpHumRead);


		//Possibly also create the bipartite graph at this stage
		//graph structure: substrate number n is vertex number 2n
		//		   reaction number n is vertex number 2n+1
		//vertex properties are a pair <substrate,reaction> with the
		//non-important one being blank

		//std::cout<<"Type: "<<tmpType<<std::endl;
		//std::cout<<"Numbers: "<<tmpsubI<<" "<<tmpProdI<<" "<<tmpnrATP<<" "<<tmpnrNADH<<" "<<tmpfreeE<<std::endl;
		//std::cout<<"Humanreadable: "<<tmpHumRead<<std::endl;

	}
}

reaction::reaction(double tmpfreechange, const std::vector<int>& tmpsubstrates, const std::vector<int>& tmproducts, int (&tmpInternalMets) [9] )  {

	freeEChange=tmpfreechange;
	substrates=tmpsubstrates;
	products=tmproducts;

	std::copy(std::begin(tmpInternalMets), std::end(tmpInternalMets), std::begin(internalMets));
}

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

double reaction::freeEchange() { return freeEChange; }


void reaction::recalcEchange(const environment& env)
{
	//recalculate the freeEchang using the formula G=G0+RT*ln(([C]^c*[D]^d)/([A]^a*[B]^b))
	//need the reactions database changed for this

	double insideLog=(std::pow(env.ppiCont,internalMets[0])*std::pow(env.piCont,internalMets[1])*std::pow(env.atpCont,internalMets[2])*std::pow(env.adpCont,internalMets[3])*std::pow(env.ampCont,internalMets[4])*std::pow(env.nadredcont,internalMets[5])*std::pow(env.nadoxcont,internalMets[6])*std::pow(env.co2cont,internalMets[7])*std::pow(env.h2ocont,internalMets[8]));

	currentFreeEChange=freeEChange+8.3144598*env.temperature*std::log(insideLog);
}

reaction::reaction()
{
	int tmpint [9] = {};
	std::copy(std::begin(tmpint), std::end(tmpint), std::begin(internalMets));
}
