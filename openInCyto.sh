#!/bin/bash

if [[ $# -eq 0 ]] ; then
    echo 'Pass the jobname as argument'
    exit 0
fi

sed "s/longerLog/$1/g" cytoscript.sh > tmpcytoscript.sh 

for i in `seq 1 10` ; do

tar -zxf $1/response.$i.tar.gz --wildcards -C $1 --strip=2 "*job$i/CP10/*job"$i"CP10NR1cell.xgmml"

propername=$(ls $1/*job$iCP*.xgmml| grep job"$i"CP)

					echo $propername
sed -i "s;file=.*job$i""CP.*xgmml;file=$propername;" tmpcytoscript.sh

done

~/Cytoscape_v3.3.0/cytoscape.sh -S tmpcytoscript.sh
