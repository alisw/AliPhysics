#!/bin/bash

folder="/lustre/alice/jthaeder/data/compressionSGE/data_2011-08-06"
version=20a


#  -----------------------------------------------------------
#  -- make results
#  -----------------------------------------------------------

if [ ! -d results ] ; then 
    mkdir -p results
fi

for useFriends in 0 1 ; do
    aliroot -b -l -q readClusters.C'("'${folder}'_Pythia","Pythia","'${version}'",30000, 30020,'${useFriends}')'
done

exit

#  -----------------------------------------------------------
#  -- draw results
#  -----------------------------------------------------------

for scale in 1 0 ; do 
    for pad in -1 0 1 2; do     
	aliroot -b -l -q drawClusters.C'("results_friends_Pythia_'${version}'.root","Pythia","'${version}'",'${pad}','${scale}')'

    done
done
