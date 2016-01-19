#!/bin/bash

cd /scratch/s1134965/frankenstein


tar xvf backup.tar.gz

nohup reading_in/build/reaction

latestfolder=$(ls -dt */| head -1)

tar acf response.from.NODENR.tar.gz --exclude reaction  $latestfolder


scp response.from.NODENR.tar.gz MOTHERHOST:FOLDERTOCOLLECT
