#!/bin/bash

seedtostartwith=RANDOMSEED
source /etc/profile
#export KRB5CCNAME=/home/s1134965/krbrealm

#set >y
cd /scratch/s1134965/frankenstein

longname="jobJOBNR"

tar xvf backup.tar.gz

nohup reading_in/build/reaction -s $seedtostartwith -j $longname

latestfolder=$(ls -dt */| head -1)

tar acf response.JOBNR.tar.gz --exclude reaction  $longname

klist -f

#scp response.from.NODENR.tar.gz MOTHERHOST:FOLDERTOCOLLECT


rsync -aPhq response.JOBNR.tar.gz MOTHERHOST:FOLDERTOCOLLECT
