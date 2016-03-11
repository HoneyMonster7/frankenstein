#!/bin/bash

if [[ $# -eq 0 ]] ; then
    echo 'Pass the jobname as argument'
    exit 0
fi

for i in `ls $1/*.tar.gz`; do


	mkdir -p $1/unfinished

	istherecp10=$(tar -ztf $i | grep CP10/$)

	if [ -z $istherecp10 ]; then
		echo "$i is unfinished"

		mv $i $1/unfinished
	fi


done
