#!/bin/bash

jobwasgiven=0
onlybest=0
onlyused=0
numberneeded=0
fittnessgraph=0
onlylastcp=0
averaging=1

junkForLastCPStays=0

optionsforchecker=""

while getopts ":bj:hun:flsa" opt; do

	#things to do: 

	#get rid of -f switch, do it always
	#do -u alongside the normal always


	case $opt in
		u)
			onlyused=1
			if [ -z "$optionsforchecker" ]; then
				optionsforchecker="$optionsforchecker-u"
			else
				optionsforchecker="$optionsforchecker -u"
			fi
			;;
		f)
			fittnessgraph=1;
			if [ -z "$optionsforchecker" ]; then
				echo "optionsforchecker was empty"
				optionsforchecker="$optionsforchecker-f"
			else
				echo "optionsforchecker was notempty"
				optionsforchecker="$optionsforchecker -f"
			fi

		;;
		n)
			echo "-n was triggered, with $OPTARG"
			numberneeded=$OPTARG
			;;
		l)
			echo "-l was triggered, with $OPTARG"
			onlylastcp=1
			;;
		s)
			echo "-s was triggered. Junk stays for CP10."
			junkForLastCPStays=1
			;;
		a)
			echo "-a was triggered. No averaging of fitness"
			averaging=0
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
			echo -e "\t -f outputs the fittness graphs too. at the moment it doesn't work together with the -u switch."
			echo -e "\t -l extract and analyse the last checkpoint only"
			echo -e "\t -s don't remove .jnk files for the last CP (replotting the SimMatrix)"
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


echo "$optionsforchecker"
#read valami

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

./moveUnfinished.sh $jobtoan

#echo "jobtoan is $jobtoan within there we have $(ls $jobtoan | grep tar.gz)"
needToUntar=$(ls $jobtoan | grep tar.gz| sort -nt. -k2,2)

if [ -z  "$needToUntar" ]; then
	echo "Can't find jobs in $jobtoan. Are you sure it's the right folder?"
	exit 1
fi

#echo "needToUntar is $needToUntar"
if [ "$numberneeded" -gt 0 ]; then
	#loopthroughthis=$(ls $jobtoan | grep tar.gz | head -n $numberneeded)
	needToUntar=$(echo "$needToUntar" | head -n $numberneeded)
fi

jobnames=$(echo "$needToUntar" | cut -d. -f2)

echo "needToUntar is $needToUntar"

#if we need to do all the checkpoints do all of them
if [ ! "$onlylastcp" == 1 ]; then
	for checkpoint in `seq 1 9`; do

		for fname in $needToUntar; do

			jobnr=$(echo $fname | cut -d. -f2)
			#echo $jobnr

			#if we need to get all the cells get all of them, otherwise get the best ones only
			if [ "$onlybest" == 0 ]; then

				tar -zxf $jobtoan/$fname --wildcards -C $jobtoan --strip=1 *job$jobnr/CP$checkpoint/"*job*CP$checkpoint*xgmml"
			else
				tar -zxf $jobtoan/$fname --wildcards -C $jobtoan --strip=1 *job$jobnr/CP$checkpoint/"*job*CP"$checkpoint"NR1cell.xgmml"
			fi

		done
		#now that we have extracted the whole checkpoint into the folder CP$checkpoint we can parse the xgmml files
		
		#USE PATH FOR THESE SCRIPTS then no copying is necessary
		cp simChecker/similarityCalc.sh $jobtoan/CP$checkpoint
	

		cp simChecker/simMatrix $jobtoan/CP$checkpoint

		cp partRatio/partratio $jobtoan/CP$checkpoint

		cd $jobtoan/CP$checkpoint

		#remove but silence the errors
		rm $checkpoint.IPR.all $checkpoint.IPR.used 2>/dev/null

		#now we should have the jnk files, the fittness list, and the list of the jnk files, so let's generate the similarity indexes

		if [ "$onlyused" == 1 ]; then
			./similarityCalc.sh -u -b
		else
			./similarityCalc.sh -b
		fi


		#for job in $jobnames; do
		for job in `ls *.xgmml | grep NR11| sed 's/.*job/job/g' | cut -db -f2 | cut -dC -f1`; do

			#echo "Jobname is $job"
			grep "job$job"CP fittness.list | sort -n -k 3 | awk '{print $1}' >listofjob.jnk
			grep "job$job"CP fittness.list | sort -n -k 3 | awk '{print $2}' >listoffit.jnk
			#cat listofjob.jnk
			./simMatrix -l listofjob.jnk -b > "$job".all.array

			sed -i 's/all/used/g' listofjob.jnk

			./simMatrix -l listofjob.jnk -b > "$job".used.array

			partall=$(./partratio -l "$job".all.array)
			partused=$(./partratio -l "$job".used.array)
			partfit=$(./partratio -l listoffit.jnk)
			#echo "for all: $partall for used: $partused for fittness: $partfit"

			echo -n "$partall " >> $checkpoint.IPR.all
			echo -n "$partused " >> $checkpoint.IPR.used
			#need to look into why the fittness one doesn't work
			echo -n "$partfit " >> $checkpoint.IPR.fitt


			
		done
			#read vlaami
		rm *.jnk
		rm *.xgmml
		#rm *.array

		cd ../..

		#read valami

	done

fi
#read valami

#testing if the last checkpoint is in a folder in the tar file (format changed at the end of january)
testthistar=$(echo "$needToUntar" | head -n 1)
testresult=$(tar -zvtf $jobtoan/$testthistar | grep CP10\/$)
testresult="forcedtobetrue"
prefix=""
antiprefix=""
if [ -z "$testresult" ]; then
	echo "No CP10 folder found"
	prefix=""
	antiprefix="CP10/"
else
	echo "CP10 folder found"
	prefix=/CP10
	antiprefix=""
fi

#read valami
mkdir -p $jobtoan/CP10/
mkdir -p $jobtoan/progress/

#doing the last checkpoint separately, as that is not in a folder
for fname in $needToUntar; do
#for fname in $( ls $jobtoan | grep tar.gz); do


	jobnr=$(echo $fname | cut -d. -f2)
	#echo $jobnr

	if [ "$onlybest" == 0 ]; then

		tar -zxf $jobtoan/$fname --wildcards -C $jobtoan/$antiprefix --strip=1 *job$jobnr/"*job*CP10*xgmml"
	else
		tar -zxf $jobtoan/$fname --wildcards -C $jobtoan/$antiprefix --strip=1 *job$jobnr/"*job*CP10NR1cell.xgmml"
	fi

		tar -zxf $jobtoan/$fname --wildcards -C $jobtoan --strip=1 *job$jobnr/"*job*.fitt"

done


#if [ $? != 0 ]; then
#	echo "There was a problem while extracting job results. Check if the folder contains jobs. Exiting."
#	exit 1
#fi

cp simChecker/similarityCalc.sh $jobtoan/CP10

cp simChecker/simMatrix $jobtoan/CP10

cp simChecker/plotter.gnup $jobtoan/CP10

cp partRatio/partratio $jobtoan/CP10

if [ "$onlybest" == 1 ]; then
	sed -i  's/i=5:words(XTICS):10/i=1:words(XTICS)/g' $jobtoan/CP10/plotter.gnup
fi	

if [ "$onlyused" == 1 ]; then
	sed -i 's/WHICH/used/g' $jobtoan/CP10/plotter.gnup
	echo "should be used"
else
	sed -i 's/WHICH/all/g' $jobtoan/CP10/plotter.gnup
	echo "should be all"
fi


cd $jobtoan/CP10

echo "options are: $optionsforchecker"

./similarityCalc.sh -b 
./similarityCalc.sh -p "${optionsforchecker[@]}"

#this is specifically for CP10, output of main simulation should be changed to be able to handle this without being a special case
checkpoint=10

#rm *.jnk
#for job in $jobnames; do
		for job in `ls *.xgmml | grep NR11| sed 's/.*job/job/g' | cut -db -f2 | cut -dC -f1`; do

	#echo "Jobname is $job"
	grep "job$job"CP fittness.list | sort -n -k 3 | awk '{print $1}' >listofjob.jnk
	grep "job$job"CP fittness.list | sort -n -k 3 | awk '{print $2}' >listoffit.jnk
	#cat listofjob.jnk
	./simMatrix -l listofjob.jnk -b > "$job".all.array

	sed -i 's/all/used/g' listofjob.jnk

	./simMatrix -l listofjob.jnk -b > "$job".used.array

	partall=$(./partratio -l "$job".all.array)
	partused=$(./partratio -l "$job".used.array)
	partfit=$(./partratio -l listoffit.jnk)
	#echo "for all: $partall for used: $partused for fittness: $partfit"

	echo -n "$partall " >> $checkpoint.IPR.all
	echo -n "$partused " >> $checkpoint.IPR.used
	#need to look into why the fittness one doesn't work
	echo -n "$partfit " >> $checkpoint.IPR.fitt

done

for jobfile in `ls ../*.fitt `; do
	job="$(echo $jobfile| sed 's/.*job//g' | cut -d. -f1)"
	#now calculating the running average of the best fittnesses

	#sed -ni '/Current/!p' ../job$job.fitt

	awk '{print $1}' "$jobfile" > "../progress/$job.numbersonly"
	awk '{print $2}' "$jobfile" > "../progress/$job.fittonly"
	awk '{print $3}' "$jobfile" > "../progress/$job.enthonly"
	awk '{print $4}' "$jobfile" > "../progress/$job.avgfitonly"
	awk '{print $5}' "$jobfile" > "../progress/$job.netsizeonly"
	awk '{print $6}' "$jobfile" > "../progress/$job.avgnetsizeonly"
	awk '{print $7}' "$jobfile" > "../progress/$job.bestusedonly"
	awk '{print $8}' "$jobfile" > "../progress/$job.avgusedonly"

	#echo "job$job" > "../$job.fittavg"
	#echo "job$job" > "../$job.enthavg"
	#echo "job$job" 

	if [[ "$averaging" != 1 ]]; then
		cat "../progress/$job.fittonly" >> "../progress/$job.fittavg"
		cat "../progress/$job.enthonly" >> "../progress/$job.enthavg"
		cat "../progress/$job.avgfitonly" >> "../progress/$job.avgfitavg"
		cat "../progress/$job.netsizeonly" >> "../progress/$job.netsizeavg"
		cat "../progress/$job.avgnetsizeonly" >> "../progress/$job.avgnetsizeavg"
		cat "../progress/$job.bestusedonly" >> "../progress/$job.bestusedavg"
		cat "../progress/$job.avgusedonly" >> "../progress/$job.avgusedavg"
	else

		../../movAvg/movAvg 100 "../progress/$job.fittonly" >> "../progress/$job.fittavg"
		../../movAvg/movAvg 100 "../progress/$job.enthonly" >> "../progress/$job.enthavg"
		../../movAvg/movAvg 100 "../progress/$job.avgfitonly" >> "../progress/$job.avgfitavg"
		../../movAvg/movAvg 100 "../progress/$job.netsizeonly" >> "../progress/$job.netsizeavg"
		../../movAvg/movAvg 100 "../progress/$job.avgnetsizeonly" >> "../progress/$job.avgnetsizeavg"
		#../../movAvg/movAvg 100 "../progress/$job.bestusedonly" >> "../progress/$job.bestusedavg"
		cat "../progress/$job.bestusedonly" >> "../progress/$job.bestusedavg"
		../../movAvg/movAvg 100 "../progress/$job.avgusedonly" >> "../progress/$job.avgusedavg"
	fi
	
	paste ../progress/$job.numbersonly ../progress/$job.fittavg ../progress/$job.enthavg ../progress/$job.avgfitavg ../progress/$job.netsizeavg ../progress/$job.avgnetsizeavg ../progress/$job.bestusedavg ../progress/$job.avgusedonly| grep ^[0-9] > ../progress/$job.progress






	
done
	
# leaving the jnk files for the last checkpoint only
if [[ "$junkForLastCPStays" != 1 ]]; then
	rm *.jnk

	rm *.xgmml
	rm ../progress/*.fittonly
	rm ../progress/*.enthonly
	rm ../progress/*.avgfitonly
	rm ../progress/*.netsizeonly
	rm ../progress/*.avgnetsizeonly
	rm ../progress/*.bestusedonly
	rm ../progress/*.avgusedonly


	rm ../progress/*.numbersonly
	rm ../progress/*.fittavg
	rm ../progress/*.enthavg
	rm ../progress/*.avgfitavg
	rm ../progress/*.netsizeavg
	rm ../progress/*.avgnetsizeavg
	rm ../progress/*.bestusedavg
	rm ../progress/*.avgusedavg
fi
#paste ../*.fittavg > ../fittavgs.fitt
#paste ../*.enthavg > ../enthavgs.enth

#now that every checkpoint's IPR's have been calculated, let's gather them into a common file

cd ..

rm completeIPR.all

for job in $jobnames; do
	#creating the header for the file
	echo -n "$job " >> completeIPR.all
done

echo " " >> completeIPR.all
cp completeIPR.all completeIPR.used
cp completeIPR.all completeIPR.fittness

for cp in `seq 1 10`; do
	cat CP$cp/$cp.IPR.all >>completeIPR.all
	echo " " >> completeIPR.all
	cat CP$cp/$cp.IPR.used >>completeIPR.used
	echo " " >> completeIPR.used

	#cat CP$cp/$cp.IPR.fitt 
	cat CP$cp/$cp.IPR.fitt >>completeIPR.fittness
	echo " " >> completeIPR.fittness

done

cp ../simChecker/lineplotter.gnup progress
cp ../simChecker/IPRplotter.gnup .

cd progress

gnuplot -persist lineplotter.gnup

cd ..
gnuplot --persist IPRplotter.gnup

#if [ "$onlyused" == 1 ]; then
#
#	./similarityCalc.sh -u 
#else
#	./similarityCalc.sh
#fi

#./simMatrix

#rm *.xgmml 2> /dev/null
#rm *.fitt
#rm similarityCalc.sh
#rm simMatrix
#rm plotter.gnup
