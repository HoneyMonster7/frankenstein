#!/bin/bash

onlyUsed=0

rm columns.list
rm rows.list

while getopts ":uh" opt; do


case $opt in
	u)
		#echo "-u was triggered"
		onlyUsed=1
		;;
	h)
		echo "Allowed options:"
		echo -e "\t -u output only the used relations"
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

	jnkname=${i/xgmml/jnk}
	if [ "$onlyUsed" == 0 ]; then
		grep Reaction -B1 $i | grep label | cut -d\" -f2 | sort -n > $jnkname
	else
		grep -A 1 "source=\"\-[0-9]*" $i | grep -v "\-\-" | paste - - | grep -v "value=\"0\"" | cut -d\" -f4 | tr -d - | sort -n | uniq > $jnkname
	fi



	echo $jnkname >> columns.list
	echo -n "$jnkname " | sed 's/CP10//g' | sed 's/cell.jnk//g' >> rows.list


	counter=$(($counter +1))

done

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
./simMatrix > similarityMatrix.dat
fi

echo "Removing .jnk files now."
rm *.jnk

gnuplot plotter.gnup
