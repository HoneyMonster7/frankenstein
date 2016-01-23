#!/bin/bash

seedtostartwith=RANDOMSEED
source /etc/profile
#export KRB5CCNAME=/home/s1134965/krbrealm

#set >y
cd /scratch/s1134965/frankenstein


tar xvf backup.tar.gz

nohup reading_in/build/reaction -s $seedtostartwith

latestfolder=$(ls -dt */| head -1)

tar acf response.from.job.JOBNR.tar.gz --exclude reaction  $latestfolder

klist -f

#scp response.from.NODENR.tar.gz MOTHERHOST:FOLDERTOCOLLECT


rsync -aPhq response.from.NODENR.tar.gz MOTHERHOST:FOLDERTOCOLLECT
