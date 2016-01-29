//to find the matching reactions between two metabolic networks
//the shell script extracts only the reaction numbers sorted in .jnk files
//this will read them in and compare them

#include <iostream>
#include <algorithm>
#include <vector>
#include <fstream>
#include <set>
#include<unistd.h>

int main(int argc, char **argv){

	std::string filenames;
	int nrReadNow=1;
	bool versusFirst=false;

	std::string listOfFiles="columns.list";

	int c;

	while ((c=getopt(argc,argv,"hl:b")) !=-1)
		switch(c)
		{
			case 'l':
				listOfFiles=optarg;
				break;
			case 'b':
				versusFirst=true;
				break;
			case 'h':
				std::cout<<"Accepted options:"<<std::endl;
				std::cout<<"\t -l [fileName] file containing the filenames that have the integers for the reactions in the given networks"<<std::endl;
				std::cout<<"\t -b compares fittnesses only against the first cell read in (you might want to order inputs accordingly)"<<std::endl;
				exit(0);
				break;
			case '?':
				if (optopt=='l')
					std::cout<<"Option -l requires a valid filename"<<std::endl;
				else
					std::cout<<"Unknown option. Try -h for allowed options."<<std::endl;
				return 1;
			default:
				exit(1);
		}




		


	std::ifstream infile(listOfFiles);

	std::vector<std::set<int>> reacMatrix;
	std::vector<std::string> names;

	if (infile)
	{
		while (infile >> filenames){

			std::ifstream nowReading(filenames);
			std::vector<int> tmpVector;
			std::set<int> tmpSet;
			int tmpInt;

			names.emplace_back(filenames);

			while (nowReading>>tmpInt){

				tmpVector.emplace_back(tmpInt);
				tmpSet.insert(tmpInt);
			}
			reacMatrix.emplace_back(tmpSet);

		}
	}	
	else {

		std::cout<<"File "<<listOfFiles<<" not found. Exiting now."<<std::endl;
		exit(1);
	}

	//for (auto writingOut:reacMatrix){


	//	std::cout<<"NR "<<nrReadNow<<std::endl;
	//	for (auto j:writingOut){

	//		std::cout<<j<<", ";
	//	}
	//	std::cout<<std::endl;
	//	nrReadNow++;
	//}
	
	if (versusFirst){

			std::set<int> firstSet=reacMatrix[0];
			for (int j=1; j<reacMatrix.size();j++)
			{

				std::set<int> secondSet=reacMatrix[j];

				std::vector<int> similarities;

				std::set_intersection(firstSet.begin(),firstSet.end(),secondSet.begin(),secondSet.end(),back_inserter(similarities));
				
				//std::cout<<"Between "<<1<<" and "<<j+1<<" there are "<<similarities.size()<<" similarities"<<std::endl;
				double coeff=similarities.size()/(double)std::max(firstSet.size(),secondSet.size());
				std::cout<<coeff<<std::endl;
				//std::cout<<coeff<<" "<<names[0]<<" "<<names[j]<<std::endl;
			}
	}
	else{
		for (int i=0;i<reacMatrix.size();i++){

			std::set<int> firstSet=reacMatrix[i];

			for (int j=0; j<reacMatrix.size();j++)
			{

				std::set<int> secondSet=reacMatrix[j];

				std::vector<int> similarities;

				std::set_intersection(firstSet.begin(),firstSet.end(),secondSet.begin(),secondSet.end(),back_inserter(similarities));
				
				//std::cout<<"Between "<<i+1<<" and "<<j+1<<" there are "<<similarities.size()<<" similarities"<<std::endl;
				double coeff=similarities.size()/(double)std::max(firstSet.size(),secondSet.size());
				std::cout<<coeff<<" ";
				//std::cout<<coeff<<" "<<names[i]<<" "<<names[j]<<std::endl;
			}
			std::cout<<std::endl;
		}
	}


}
