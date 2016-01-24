#!/bin/bash


#first asking which job 

echo "Which job do you want analysed?"

read jobtoan

#jobtoan="newtars"

if [ ! -d $jobtoan ]; then
	echo "$jobtoan folder doesn't exist. "
	exit 1
fi

if [ ! $( ls $jobtoan | grep tar.gz) ]; then
	echo "Can't find jobs in $jobtoan. Are you sure it's the right folder?"
	exit 1
fi

#ls $jobtoan | grep tar.gz
for fname in $( ls $jobtoan | grep tar.gz); do


	jobnr=$(echo $fname | cut -d. -f2)

	tar -zxvf $jobtoan/$fname --wildcards -C $jobtoan --strip=1 job$jobnr/"job*CP10*xgmml"


done

if [ $? != 0 ]; then
	echo "There was a problem while extracting job results. Check if the folder contains jobs. Exiting."
	exit 1
fi

cp simChecker/similarityCalc.sh $jobtoan

cp simChecker/simMatrix $jobtoan

cp simChecker/plotter.gnup $jobtoan
cd $jobtoan

./similarityCalc.sh

#./simMatrix

#rm similarityCalc.sh
#rm simMatrix
#rm plotter.gnup
