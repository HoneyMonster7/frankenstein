#!/bin/bash


#first asking which job 

echo "Which job do you want analysed?"

#read jobtoan

jobtoan="newtars"

if [ ! -d $jobtoan ]; then
	echo "$jobtoan folder doesn't exist. "
	exit 1
fi

#ls $jobtoan | grep tar.gz
for fname in $( ls $jobtoan | grep tar.gz); do


	jobnr=$(echo $fname | cut -d. -f2)

	tar -zxvf $jobtoan/$fname --wildcards -C $jobtoan --strip=1 job$jobnr/"job*CP10*xgmml"


done

cp simChecker/similarityCalc.sh $jobtoan

cp simChecker/simMatrix $jobtoan

cd $jobtoan

./similarityCalc.sh

#./simMatrix

#rm similarityCalc.sh
#rm simMatrix
