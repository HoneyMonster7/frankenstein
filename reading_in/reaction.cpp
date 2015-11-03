#include "reaction.h"


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


void reaction::readCompounds(std::string fileName, ReactionNetwork& graph, std::vector<Vertex>& vertexList, std::vector<Vertex>& internals)
{
	std::ifstream inFile(fileName);
	std::string line, nameChem, namenorm;

	if(!inFile)
	{
		throw std::runtime_error("Can't open input file " + fileName);
	}

	while(std::getline(inFile, line))
	{
		substrate tmpsubstrate;

		std::stringstream iss(line);
		iss>>tmpsubstrate.index;
		iss>>tmpsubstrate.freeOfCreation;
		iss>>nameChem;
		iss>>tmpsubstrate.name;
		iss>>tmpsubstrate.charge;
		vertexList.emplace_back(boost::add_vertex(graph));

		graph[vertexList[vertexList.size()-1]].sub=tmpsubstrate;

		if(tmpsubstrate.index<0){
			internals.emplace_back(vertexList[vertexList.size()-1]);
		}

	}
}

int reaction::getListNr(){ return listNR;}

void reaction::readReactions(std::string fileName, std::vector<reaction>& reacPointer, ReactionNetwork& graph, std::vector<Vertex>& vertexList, const std::vector<Vertex>& compoundVList)
{
	std::ifstream inFile(fileName);
	std::string line,tmpsubs,tmpprods;

	InternalMetsT tmpinternalMets = {};

	double tmpfreeE;

	typedef boost::graph_traits<ReactionNetwork>::edge_descriptor Edge;

	if(!inFile)
	{
		throw std::runtime_error("Can't open input file " + fileName);
	}

	int counter=0;
	while(std::getline(inFile, line))
	{
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
		for(int i=-13; i<0; i++)
		{
			if(std::find(tmpsubstrates.begin(), tmpsubstrates.end(), i) != tmpsubstrates.end())
			{
				tmpinternalMets[i+13]--;
			}

			if(std::find(tmproducts.begin(), tmproducts.end(), i) != tmproducts.end())
			{
				tmpinternalMets[i+13]++;
			}
		}


		reacPointer.emplace_back(counter,tmpfreeE, tmpsubstrates, tmproducts, tmpinternalMets);


		vertexList.emplace_back(boost::add_vertex(graph));
		graph[vertexList[vertexList.size()-1]].reac=reaction(counter,tmpfreeE,tmpsubstrates,tmproducts,tmpinternalMets);


		for(int i : tmpsubstrates)
		{
			Edge e1;
			e1=(boost::add_edge(compoundVList[i+13],vertexList[vertexList.size()-1],graph)).first;
		}


		for(int i: tmproducts)
		{
			Edge e1;
			e1=(boost::add_edge(vertexList[vertexList.size()-1],compoundVList[i+13],graph)).first;
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

	currentFreeEChange=freeEChange+8.3144598*env.temperature*std::log(insideLog);
}


void reaction::calcThroughput(const int NrCompounds,ReactionNetwork& graph, std::vector<Vertex> reacList){

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

 std::vector<Vertex> reaction::subsetVertices( std::vector<int> vertexIDs, std::vector<Vertex> reacList){
	 //sort not required? 
	 std::sort(vertexIDs.begin(),vertexIDs.end());

	 std::vector<Vertex> tobeReturned;

	 for(int i :vertexIDs){
		 //uses i-1 in order to correspond to line numbers in the input reaction list
		 tobeReturned.push_back(reacList[i-1]);
	 }
	 return tobeReturned;
 }
