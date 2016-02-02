#!/bin/bash

onlyUsed=0
fitplotsneeded=0
plotsneeded=0
simVsBest=0

rm columns.list
rm columns.all.list
rm columns.used.list
rm rows.list
rm fittness.list

while getopts ":ufhpb" opt; do


case $opt in
	u)
		echo "-u was triggered"
		onlyUsed=1
		;;
	b)
		echo "-b was triggered"
		simVsBest=1
		;;
	f)
		echo "-f was triggered"
		fitplotsneeded=1
		;;
	p)
		echo "-p fwas triggered"
		plotsneeded=1
		;;
	h)
		echo "Allowed options:"
		echo -e "\t -u output only the used relations"
		echo -e "\t -f also produces the plot for the fittnesses"
		echo -e "\t -p creates plots of the similarityMatrix"
		echo -e "\t -b calculates the fittness indices against the best performing network only. Array of indices produced instead of matrix."
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


counter=1;
for i in `ls *.xgmml`
do

	alljnkname=${i/xgmml/all.jnk}
	usedjnkname=${i/xgmml/used.jnk}
	jobname=${i/.xgmml/}
	numberofcelltmp=${i/cell.xgmml/}
	numberofcell=$(echo "$numberofcelltmp" | sed 's/job[0-9]*CP[0-9]*NR//g')
	#if [ "$onlyUsed" == 0 ]; then
		grep Reaction -B1 $i | grep label | cut -d\" -f2 | sort -n > $alljnkname
	#else
		grep -A 1 "source=\"\-[0-9]*" $i | grep -v "\-\-" | paste - - | grep -v "value=\"0\"" | cut -d\" -f4 | tr -d - | sort -n | uniq > $usedjnkname
	#fi

	fittness=$(grep Fittness $i | cut -d \" -f6)

	echo "$alljnkname $fittness $numberofcell" >> fittness.list

	echo $alljnkname >> columns.all.list
	echo $usedjnkname >> columns.used.list
	echo -n "$alljnkname " | sed 's/CP10//g' | sed 's/cell.all.jnk//g' >> rows.list


	counter=$(($counter +1))

done

cat fittness.list | sort -nr -k2 | awk '{print $1}' > columns.all.list.ordered
cat fittness.list | sort -nr -k2 | awk '{print $1}'| sed 's/all/used/g'  > columns.used.list.ordered
#these two below don't work for some reason...
sed 's/cell.all.jnk//g' fittness.list | awk '{print $1}' | tr '\r\n' ' ' > newrows.all.list
grep Fittness *.xgmml | cut -d \" -f6 | sort -nr > fittnessbulk.list

#this is now done by simMatrix much faster (doesn't have to read from files every iteration)

##for i in `ls *.jnk`
##do
##	ilength=`wc -l <$i`
##
##	
##
##	for j in `ls *.jnk`
##	do
##
##		jlength=`wc -l <$j`
##
##		#echo "ilength is $ilength, jlength is $jlength"
##
##
##		#cat $i $j|sort -n | uniq -c| grep 
##		matcoeff=`cat $i $j|sort -n | uniq -c | grep "^[[:space:]]*2" | wc -l`
##		
##		largerone=$((ilength > jlength ? ilength : jlength))
##		matdecim=`echo "scale=4; $matcoeff/$largerone" | bc`
##		
##		if [ "$matdecim" != "1.0000" ]
##		then
##		echo -n "$matdecim "
##		else 
##			echo -n "0 "
##		fi
##	done
##	echo 
##done

if [ ! -e simMatrix ]; then
	echo "Can't find the simMatrix executable in this folder. columns.list has been generated, but you have to run simMatrix by hand."
	exit
else
	if [ "$simVsBest" == 1 ]; then
		#only calculating the similarity index versus the best performing network	
		./simMatrix -b -l columns.all.list.ordered > similarityArray.all.dat
		./simMatrix -b -l columns.used.list.ordered > similarityArray.used.dat
	else
		#calculating the whole similarity matrix
		if [ "$onlyUsed" == 1 ]; then
			./simMatrix -l columns.used.list > similarityMatrix.dat
		else
			./simMatrix -l columns.all.list > similarityMatrix.dat
		fi
	fi
fi

echo "Removing .jnk files now."
#rm *.jnk

if [ "$plotsneeded" == 1 ]; then
	gnuplot -persist plotter.gnup
	if [ "$fitplotsneeded" == 1 ]; then

		gnuplot -persist lineplotter.plot
	fi
fi