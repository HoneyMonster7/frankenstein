#!/bin/bash


rm columns.list

for i in `ls *.xgmml`
do

	jnkname=${i/xgmml/jnk}
	grep Reaction -B1 $i | grep label | cut -d\" -f2 | sort -n > $jnkname

	echo $jnkname >> columns.list
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
