#!/bin/bash

jobwasgiven=0
onlybest=0
onlyused=0
numberneeded=0
fittnessgraph=0
onlylastcp=0

optionsforchecker=""

while getopts ":bj:hun:fl" opt; do


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

#echo "jobtoan is $jobtoan within there we have $(ls $jobtoan | grep tar.gz)"
needToUntar=$(ls $jobtoan | grep tar.gz)

#if [ -z  $(ls $jobtoan | grep tar.gz) ]; then
if [ -z  "$needToUntar" ]; then
	echo "Can't find jobs in $jobtoan. Are you sure it's the right folder?"
	exit 1
fi

#echo "needToUntar is $needToUntar"
#ls $jobtoan | grep tar.gz
if [ "$numberneeded" -gt 0 ]; then
	#loopthroughthis=$(ls $jobtoan | grep tar.gz | head -n $numberneeded)
	needToUntar=$(echo "$needToUntar" | head -n $numberneeded)
#else
#	loopthroughthis=$(ls $jobtoan | grep tar.gz)
fi

jobnames=$(echo "$needToUntar" | cut -d. -f2)

#echo "jobnames are: $jobnames"
#
#echo "we need to loop through $loopthroughthis"
#
echo "needToUntar is $needToUntar"

if [ ! "$onlylastcp" == 1 ]; then
	for checkpoint in `seq 1 9`; do

		for fname in $needToUntar; do
		#for fname in $( ls $jobtoan | grep tar.gz); do


			jobnr=$(echo $fname | cut -d. -f2)
			echo $jobnr

			if [ "$onlybest" == 0 ]; then

				tar -zxvf $jobtoan/$fname --wildcards -C $jobtoan --strip=1 job$jobnr/CP$checkpoint/"job*CP$checkpoint*xgmml"
			else
				tar -zxvf $jobtoan/$fname --wildcards -C $jobtoan --strip=1 job$jobnr/CP$checkpoint/"job*CP"$checkpoint"NR1cell.xgmml"
			fi

		done
		#now that we have extracted the whole checkpoint into the folder CP$checkpoint we can parse the xgmml files
		cp simChecker/similarityCalc.sh $jobtoan/CP$checkpoint

		cp simChecker/simMatrix $jobtoan/CP$checkpoint

		cp partRatio/partratio $jobtoan/CP$checkpoint

		cd $jobtoan/CP$checkpoint

		rm $checkpoint.IPR.all $checkpoint.IPR.used

		#now we should have the jnk files, the fittness list, and the list of the jnk files, so let's generate the similarity indexes

		if [ "$onlyused" == 1 ]; then
			./similarityCalc.sh -u -b
		else
			./similarityCalc.sh -b
		fi


		for job in $jobnames; do

			echo "Jobname is $job"
			grep "job$job" fittness.list | sort -n -k 3 | awk '{print $1}' >listofjob.jnk
			grep "job$job" fittness.list | sort -n -k 3 | awk '{print $2}' >listoffit.jnk
			#cat listofjob.jnk
			./simMatrix -l listofjob.jnk -b > "$job".all.array

			sed -i 's/all/used/g' listofjob.jnk

			./simMatrix -l listofjob.jnk -b > "$job".used.array

			partall=$(./partratio -l "$job".all.array)
			partused=$(./partratio -l "$job".used.array)
			partfit=$(./partratio -l listoffit.jnk)
			echo "for all: $partall for used: $partused for fittness: $partfit"

			echo -n "$partall " >> $checkpoint.IPR.all
			echo -n "$partused " >> $checkpoint.IPR.used
			#need to look into why the fittness one doesn't work
			echo -n "$partfit " >> $checkpoint.IPR.fitt


			
		done
			#read vlaami
		rm *.jnk
		rm *.xgmml
		rm *.array

		cd ../..

		#read valami

	done

fi
#read valami
mkdir -p $jobtoan/CP10/

#doing the last checkpoint separately, as that is not in a folder
for fname in $needToUntar; do
#for fname in $( ls $jobtoan | grep tar.gz); do


	jobnr=$(echo $fname | cut -d. -f2)
	echo $jobnr

	if [ "$onlybest" == 0 ]; then

		tar -zxvf $jobtoan/$fname --wildcards -C $jobtoan/CP10/ --strip=1 job$jobnr/"job*CP10*xgmml"
	else
		tar -zxvf $jobtoan/$fname --wildcards -C $jobtoan/CP10/ --strip=1 job$jobnr/"job*CP10NR1cell.xgmml"
	fi

	if [ "$fittnessgraph" == 1 ]; then

		tar -zxvf $jobtoan/$fname --wildcards -C $jobtoan --strip=1 job$jobnr/"job*.fitt"
	fi

done


if [ $? != 0 ]; then
	echo "There was a problem while extracting job results. Check if the folder contains jobs. Exiting."
	exit 1
fi

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

for job in $jobnames; do

	echo "Jobname is $job"
	grep "job$job" fittness.list | sort -n -k 3 | awk '{print $1}' >listofjob.jnk
	grep "job$job" fittness.list | sort -n -k 3 | awk '{print $2}' >listoffit.jnk
	#cat listofjob.jnk
	./simMatrix -l listofjob.jnk -b > "$job".all.array

	sed -i 's/all/used/g' listofjob.jnk

	./simMatrix -l listofjob.jnk -b > "$job".used.array

	partall=$(./partratio -l "$job".all.array)
	partused=$(./partratio -l "$job".used.array)
	partfit=$(./partratio -l listoffit.jnk)
	echo "for all: $partall for used: $partused for fittness: $partfit"

	echo -n "$partall " >> $checkpoint.IPR.all
	echo -n "$partused " >> $checkpoint.IPR.used
	#need to look into why the fittness one doesn't work
	echo -n "$partfit " >> $checkpoint.IPR.fitt


	
done
	

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

	cat CP$cp/$cp.IPR.fitt 
	cat CP$cp/$cp.IPR.fitt >>completeIPR.fittness
	echo " " >> completeIPR.fittness

done



#if [ "$onlyused" == 1 ]; then
#
#	./similarityCalc.sh -u 
#else
#	./similarityCalc.sh
#fi

#./simMatrix

rm *.xgmml
rm *.fitt
#rm similarityCalc.sh
#rm simMatrix
#rm plotter.gnup
