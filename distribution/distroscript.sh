#!/bin/bash

# NEED: TAKE A NOTE OF THE RANDOM SEED, POTENTIALLY MAKE IT REALLY RANDOM


# SHOULD BE RUN IN THE FOLDER CONTAINING reading_in AND data WITH THE NECESSARY FILES

thishostsname=`hostname | cut -d. -f1`
thisdirectory=`pwd `
compNRLength=3

nrofmachines=10

firstToTry=88

echo "Find a job name:"
read jobname

if [ -d "$jobname" ]; then

	echo "Job name already used. Remove or find an other name. Exiting now. "
	exit
fi

if [ ! -e reading_in/build/reaction ]; then
	echo "Can't find the executable. Are you running in frankenstein?"
	exit
fi

if [ ! -d data ]; then
	echo "Can't find the data folder. Are you running in frankenstein?"
	exit
fi

mkdir $jobname

#first get the kerberos credentials right



#now get to work

for i in `seq 1 "$nrofmachines"`
do

	#this is to make new random numbers (not always the same)
	modrnd=$(($i + 50))
	#change the seed in main.cpp will need to change that file actually
	#sed -i  "s/RandomGeneratorType\ generator(1);/RandomGeneratorType\ generator($i);/" reading_in/main.cpp

	##compile reaction 

	#cd reading_in
	#cd build

	#make

	#cd ../../

	# CHANGE THE RANDOM SEED BACK SO NEXT CYCLE CAN USE UNMODIFIED MAIN.CPP
	#sed -i "s/RandomGeneratorType\ generator($i);/RandomGeneratorType\ generator(1);/" reading_in/main.cpp


	#echo "executable built, main.cpp changed back, now tarring"

	# tar the exectuable and the data files needed

	echo "The jobname is $jobname" > whichjob.txt
	tar acf backup.tar.gz data/fullnew*.dat reading_in/build/reaction whichjob.txt

	echo "tarring completed, now looking for a node that is up"


	#send the script that runs the system out


	# FIRST WE FIND THE NEXT AVAILABLE HOST

	#hostName=`printf "cplab%0*d\n" $compNRLength $firstToTry`
	hostname=`printf "cplab%0*d\n" $compNRLength $firstToTry`
	echo "now trying $hostname"

	isHostUp=`ssh -o BatchMode=yes -o ConnectTimeout=5 "$hostname" echo 1 2>&1`
	
	if  [["$isHostUp" -eq "1" ]]; then
		isReactionAlreadyRunning=`ssh "$hostname" ps ax -u s1134965 | grep reaction | grep -v grep`
		if [ -z "isReactionAlreadyRunning" ]; then
			isHostUp=0
			echo "Reaction is already running on $hostname"
		fi
	fi
	echo "tried $hostname"
	while [[ "$isHostUp" != "1" ]]
	do
	
		firstToTry=$((firstToTry+1))
		hostname=`printf "cplab%0*d\n" $compNRLength $firstToTry`
		echo " no luck there, trying $hostname"
		isHostUp=`ssh -o BatchMode=yes -o ConnectTimeout=5 "$hostname" echo 1 2>&1`
		if [[ "$isHostUp" -eq 1 ]]; then
			isReactionAlreadyRunning=`ssh "$hostname" ps ax -u s1134965 | grep reaction | grep -v grep`
			if [ -z "isReactionAlreadyRunning" ]; then
				isHostUp=0
				echo "Reaction is already running on $hostname"
			fi
		fi
		#if [ "$isHostUp" -eq 1 ].
		#then
		#	echo "$hostname is reachable"
		#
		#fi
	done
	echo " found $hostname personalizing the script for it"

	scriptname="oncurrentNode$jobname$i"

	echo "hostname is $hostname, thishostsname is $thishostsname, thisfolder is $thisdirectory"
	#echo "sed \"s/NODENR/$hostname/\" distribution/onNode.sh | sed  \"s/MOTHERHOST/$thishostsname/\"| sed  \"s#FOLDERTOCOLLECT#$thisdirectory#\" > oncurrentNode.sh"
	sed "s/JOBNR/$i/g" distribution/onNode.sh | sed  "s/MOTHERHOST/$thishostsname/g"| sed  "s#FOLDERTOCOLLECT#$thisdirectory/$jobname#g"| sed "s/RANDOMSEED/$modrnd/" > $scriptname


	#sed -i "s/MOTHERHOST/$thishostsname/g/" oncurrentNode.sh
	#sed -i "s/FOLDERTOCOLLECT/$thisdirectory/g/" oncurrentNode.sh

	chmod +x oncurrentNode.sh

	echo "$hostname is running job $i with random seed $modrnd" >> $jobname/joblist.txt
	echo "sending script to $hostname"
	rsync -aPhq {$scriptname,backup.tar.gz} "$hostname":/scratch/s1134965/frankenstein
	echo "running script at $hostname"
	ssh  -t "$hostname" 'bash -l -c /scratch/s1134965/frankenstein/$scripname >/scratch/s1134965/frankenstein/out 2>/scratch/s1134965/frankenstein/err </dev/null &' &

	mv $scriptname distribution/tmp/

	firstToTry=$((firstToTry+1))

	
done

