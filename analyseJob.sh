#!/bin/bash

jobwasgiven=0
onlybest=0
onlyused=0
numberneeded=0

while getopts ":bj:hun:" opt; do


	case $opt in
		u)
			onlyused=1
			;;
		n)
			echo "-n was triggered, with $OPTARG"
			numberneeded=$OPTARG
			;;
		b)
			echo "-b was triggered"
			onlybest=1
			;;
		j)
			echo "-j was triggered"
			jobtoan=$OPTARG
			jobwasgiven=1
			;;
		h)
			echo "Script to analyse jobs that were run on the cplab."
			echo "Needs a job folder with tar.gz files collected from the computing nodes. This can either be provided with the -j [jobName] switch, or if omitted the script will ask for it."
			#echo "By default it uses the best 10 cells, if the -b flag is used only the best will be usedfrom each node"
			echo -e "\t -b use only the best network from each node. Default is the best 10."
			echo -e "\t -u calculate the similarity index using only the reactions with nonzero flux."
			echo -e "\t -j [jobName] jobs are provided in jobName  Default will ask user for folder."
			echo -e "\t -n [number] number of jobs to analyse. Will use the first [number] jobs. 0 for all jobs. Defaults to 0."
			exit 0
			;;
		\?)
			echo "Invalid option: -$OPTARG Try -h for help."
			exit 1
			;;
		:)
			echo "Option -$OPTARG requires an argument."
			exit 1
			;;
	esac
	

done


#first asking which job 

if [ "$jobwasgiven" == 0 ]; then

echo "Which job do you want analysed?"

read jobtoan
fi

#jobtoan="newtars"


if [ ! -d $jobtoan ]; then
	echo "$jobtoan folder doesn't exist. "
	exit 1
fi

#echo "jobtoan is $jobtoan within there we have $(ls $jobtoan | grep tar.gz)"
needToUntar=$(ls $jobtoan | grep tar.gz)

#if [ -z  $(ls $jobtoan | grep tar.gz) ]; then
if [ -z  "$needToUntar" ]; then
	echo "Can't find jobs in $jobtoan. Are you sure it's the right folder?"
	exit 1
fi

#ls $jobtoan | grep tar.gz
if [ "$numberneeded" > 0 ]; then
	#loopthroughthis=$(ls $jobtoan | grep tar.gz | head -n $numberneeded)
	needToUntar=$(echo "$needToUntar" | head -n $numberneeded)
#else
#	loopthroughthis=$(ls $jobtoan | grep tar.gz)
fi
#
#echo "we need to loop through $loopthroughthis"
#
#echo "needToUntar is $needToUntar"
for fname in $needToUntar; do
#for fname in $( ls $jobtoan | grep tar.gz); do


	jobnr=$(echo $fname | cut -d. -f2)

	if [ "$onlybest" == 0 ]; then

		tar -zxvf $jobtoan/$fname --wildcards -C $jobtoan --strip=1 job$jobnr/"job*CP10*xgmml"
	else
		tar -zxvf $jobtoan/$fname --wildcards -C $jobtoan --strip=1 job$jobnr/"job*CP10NR1cell.xgmml"
	fi


done

if [ $? != 0 ]; then
	echo "There was a problem while extracting job results. Check if the folder contains jobs. Exiting."
	exit 1
fi

cp simChecker/similarityCalc.sh $jobtoan

cp simChecker/simMatrix $jobtoan

if [ "$onlybest" == 0 ]; then
	cp simChecker/plotter.gnup $jobtoan
else
	cat simChecker/plotter.gnup | sed 's/i=5:words(XTICS):10/i=1:words(XTICS)/g' >$jobtoan/plotter.gnup

fi	
cd $jobtoan

if [ "$onlyused" == 1 ]; then

	./similarityCalc.sh -u 
else
	./similarityCalc.sh
fi

#./simMatrix

rm *.xgmml
#rm similarityCalc.sh
#rm simMatrix
#rm plotter.gnup
