#!/bin/bash

seedtostartwith=RANDOMSEED
source /etc/profile
#export KRB5CCNAME=/home/s1134965/krbrealm

#set >y
cd /scratch/s1134965/frankenstein

longname="SIMNjobJOBNR"

tar xvf backup.tar.gz

if [[ -d "$longname" ]]; then
	rm -rf $longname
fi


nohup reading_in/build/reaction -s $seedtostartwith -j $longname

echo "The random seed used to initiate the prng was $seedtostartwith" > $longname/originalseed.txt

latestfolder=$(ls -dt */| head -1)

tar acf response.JOBNR.tar.gz --exclude reaction  $longname


#scp response.from.NODENR.tar.gz MOTHERHOST:FOLDERTOCOLLECT


rsync -aPhq response.JOBNR.tar.gz MOTHERHOST:FOLDERTOCOLLECT

mv $longname previousJob
