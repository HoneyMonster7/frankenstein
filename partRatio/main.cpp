//to find the inverse participation ratio of the fittnesses 
//
//reads in a file with doubles

#include <iostream>
#include <algorithm>
#include <vector>
#include <fstream>
#include <set>
#include<unistd.h>
#include <math.h>

int main(int argc, char **argv){

	std::string filenames;
	int nrReadNow=1;
	int numberOfBins=10;

	std::string listOfFiles="columns.list";

	int c;

	while ((c=getopt(argc,argv,"hl:n:")) !=-1)
		switch(c)
		{
			case 'l':
				listOfFiles=optarg;
				break;
			case 'n':
				numberOfBins=std::stoi(optarg);
				break;
			case 'h':
				std::cout<<"Accepted options:"<<std::endl;
				std::cout<<"\t -l [fileName] file containing the fittness values"<<std::endl;
				std::cout<<"\t -n [NumberOfBins] the number of bins to use"<<std::endl;
				exit(0);
				break;
			case '?':
				if (optopt=='l')
					std::cout<<"Option -l requires a valid filename"<<std::endl;
				else if (optopt=='n')
					std::cout<<"Option -n requires a valid integer"<<std::endl;
				else
					std::cout<<"Unknown option. Try -h for allowed options."<<std::endl;
				return 1;
			default:
				exit(1);
		}




		


	std::ifstream infile(listOfFiles);

	std::vector<double> values;

	if (infile)
	{
			double tmpDouble;
		while (infile >> tmpDouble){

			values.emplace_back(tmpDouble);
		}
	}	
	else {

		std::cout<<"File "<<listOfFiles<<" not found. Exiting now."<<std::endl;
		exit(1);
	}

	std::sort(values.begin(),values.end());

	//	for (auto j:values){

	//		std::cout<<j<<", ";
	//	}
	//
	//	std::cout<<std::endl;

	double minValue=values[0];
	//need to use a bit more than the actual max, due to the errors made while doing maths 
	double maxValue=values.back()+1e-5;

	std::vector<int> numberInBins(numberOfBins);

	int currentBin=0;
	double binLength=(maxValue-minValue)/(double)numberOfBins;

	for (int i=0; i<values.size(); i++){

		if (values[i] < minValue+binLength*(currentBin+1)){
			numberInBins[currentBin]++;
		}
		else {
			currentBin++;
			numberInBins[currentBin]++;
		}
		
	}
	//std::cout<<"Bin size:"<<binLength<<std::endl;
	//std::cout<<"Numbers in bins:"<<std::endl;

	int sum=0;
	for (auto i:numberInBins){
		//std::cout<<i<<", ";
		sum+=i;
	}
	//std::cout<<std::endl;
	double dsum=(double)sum;

	std::vector<double> doubleBinNr(numberInBins.begin(),numberInBins.end());
	double sumOfSquares=0;

	for (int i=0; i<doubleBinNr.size(); i++){
		doubleBinNr[i]=doubleBinNr[i]/dsum;
		sumOfSquares+=pow(doubleBinNr[i],2);
	}
	
	double ratio=(double)1/sumOfSquares;
	//std::cout<<"The inverse participation ratio is: "<<ratio<<std::endl;
	std::cout<<ratio<<std::endl;

}
