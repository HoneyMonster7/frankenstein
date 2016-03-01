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
	bool usedProvided=false;
	bool allProvided=false;

	std::string listOfAllFiles="columns.list";
	std::string listOfUsedFiles="columns.list";

	int c;

	while ((c=getopt(argc,argv,"hl:u:b")) !=-1)
		switch(c)
		{
			case 'l':
				listOfAllFiles=optarg;
				allProvided=true;
				break;
			case 'u':
				listOfUsedFiles=optarg;
				usedProvided=true;
				break;
			case 'b':
				versusFirst=true;
				break;
			case 'h':
				std::cout<<"Accepted options:"<<std::endl;
				std::cout<<"For the calculation of the similarity matrix both -l and -u must be specified."<<std::endl;
				std::cout<<"For the similarity array -l is needed."<<std::endl;
				std::cout<<"\t -l [fileName] file containing the filenames that have the integers for all the reactions in the given networks"<<std::endl;
				std::cout<<"\t -u [fileName] file containing the filenames that have the integers for the used reactions in the given networks"<<std::endl;
				std::cout<<"\t -b compares fittnesses only against the first cell read in (you might want to order inputs accordingly) This only uses -l."<<std::endl;
				exit(0);
				break;
			case '?':
				if (optopt=='l'){
					std::cout<<"Option -l requires a valid filename"<<std::endl;
					allProvided=false;
				}
				else
					if (optopt=='u'){
						std::cout<<"Option -u requires a valid filename"<<std::endl;
						usedProvided=false;
					}
					else
						std::cout<<"Unknown option. Try -h for allowed options."<<std::endl;
				return 1;
			default:
				exit(1);
		}




	if ((allProvided && usedProvided) || (versusFirst && allProvided)){
		//everything good
	}
		else
		{std::cout<<"There was an error in providing the right number of arguments. Try -h for help. Exiting now. "<<std::endl;
			exit(1);
		}

	//reading in the all reactions
	std::ifstream infile(listOfAllFiles);

	std::vector<std::set<int>> reacMatrixAll;
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
			reacMatrixAll.emplace_back(tmpSet);

		}
	}	
	else {

		std::cout<<"File "<<listOfAllFiles<<" not found. Exiting now."<<std::endl;
		exit(1);
	}

	//reading in the used reactions
	std::ifstream inUsedFile(listOfUsedFiles);

	std::vector<std::set<int>> reacMatrixUsed;

	if (inUsedFile)
	{
		while (inUsedFile >> filenames){

			std::ifstream nowReading(filenames);
			std::vector<int> tmpVector;
			std::set<int> tmpSet;
			int tmpInt;

			names.emplace_back(filenames);

			while (nowReading>>tmpInt){

				tmpVector.emplace_back(tmpInt);
				tmpSet.insert(tmpInt);
			}
			reacMatrixUsed.emplace_back(tmpSet);

		}
	}	
	else {

		std::cout<<"File "<<listOfUsedFiles<<" not found. Exiting now."<<std::endl;
		exit(1);
	}


	//for (auto writingOut:reacMatrixAll){


	//	std::cout<<"NR "<<nrReadNow<<std::endl;
	//	for (auto j:writingOut){

	//		std::cout<<j<<", ";
	//	}
	//	std::cout<<std::endl;
	//	nrReadNow++;
	//}
	
	if (versusFirst){

			std::set<int> firstSet=reacMatrixAll[0];
			for (int j=1; j<reacMatrixAll.size();j++)
			{

				std::set<int> secondSet=reacMatrixAll[j];

				std::vector<int> similarities;

				std::set_intersection(firstSet.begin(),firstSet.end(),secondSet.begin(),secondSet.end(),back_inserter(similarities));
				
				//std::cout<<"Between "<<1<<" and "<<j+1<<" there are "<<similarities.size()<<" similarities"<<std::endl;
				double coeff=similarities.size()/(double)std::max(firstSet.size(),secondSet.size());
				std::cout<<coeff<<std::endl;
				//std::cout<<coeff<<" "<<names[0]<<" "<<names[j]<<std::endl;
			}
	}
	else{
		for (int i=0;i<reacMatrixAll.size();i++){

			std::set<int> firstAllSet=reacMatrixAll[i];
			std::set<int> firstUsedSet=reacMatrixUsed[i];

			for (int j=0; j<reacMatrixAll.size();j++)
			{

				std::set<int> secondAllSet=reacMatrixAll[j];
				std::set<int> secondUsedSet=reacMatrixUsed[j];

				std::vector<int> similarities;
				double normalizer;

				// above the diagonal in the similarity matrix we'll have the all-all comparison, below the used-used ones
				if (i>j){
					std::set_intersection(firstAllSet.begin(),firstAllSet.end(),secondAllSet.begin(),secondAllSet.end(),back_inserter(similarities));
					normalizer=std::max(firstAllSet.size(),secondAllSet.size());
				}
				else
				{
					std::set_intersection(firstUsedSet.begin(),firstUsedSet.end(),secondUsedSet.begin(),secondUsedSet.end(),back_inserter(similarities));
					normalizer=std::max(firstUsedSet.size(),secondUsedSet.size());
				}
				
				//std::cout<<"Between "<<i+1<<" and "<<j+1<<" there are "<<similarities.size()<<" similarities"<<std::endl;
				double coeff=similarities.size()/normalizer;
				std::cout<<coeff<<" ";
				//std::cout<<coeff<<" "<<names[i]<<" "<<names[j]<<std::endl;
			}
			std::cout<<std::endl;
		}
	}


}
