#!/bin/bash



for i in `ls *.xgmml`
do

	jnkname=${i/xgmml/jnk}
	grep Reaction -B1 $i | grep label | cut -d\" -f2 | sort -n > $jnkname

done


for i in `ls *.jnk`
do
	ilength=`wc -l <$i`

	echo $i> columns.list
	for j in `ls *.jnk`
	do

		jlength=`wc -l <$j`

		#echo "ilength is $ilength, jlength is $jlength"


		#cat $i $j|sort -n | uniq -c| grep 
		matcoeff=`cat $i $j|sort -n | uniq -c | grep "^[[:space:]]*2" | wc -l`
		
		largerone=$((ilength > jlength ? ilength : jlength))
		matdecim=`echo "scale=4; $matcoeff/$largerone" | bc`
		
		echo -n "$matdecim "
	done
	echo 
done


rm *.jnk
