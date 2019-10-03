#!/bin/bash

###########################
## Author : Beomkyu Kim  ##
## email  : kimb@cern.ch ##
###########################

##update your user location
basepath=/alice/cern.ch/user/k/kimb/

if [ -z "$1" ]
then
    echo "Wrong usage"
    echo "usage : ./run.sh [full] or [terminate]"
    exit 0;
fi


##taskname : task name
taskname=CD

##periods : remove some periods if you don't want to run all
periods="LHC10b LHC10c LHC10d LHC10e LHC15f"

if [ $1 = "full" ]
then 
    for i in ${periods}
    do 
        root -l -b -q runCD.C\(\"${taskname}\",\"${i}\"\)
    done
fi

if [ $1 = "terminate" ]
then
    alien-token-init  
    for i in ${periods}
    do
        perl ${ALICE_PHYSICS}/PWGUD/DIFFRACTIVE/macros/alien_cp.pl ${basepath}/${taskname}${i} root_archive.zip AnalysisResults.root
        cd alice
        Files=`find $PWD -name AnalysisResults.root`
        n=0
        nFiles=$(( $(echo $Files | perl -ne "print s/ //g;") + 1 )) 
        echo "Total files $nFiles"
        PartFiles=" "
        t=0
        FullFiles=" "
        for j in $Files
        do
            PartFiles=" $PartFiles $j "
            n=$(( $n + 1 ))
            t=$(( $t + 1 ))
            if [ $n -eq 25 ] || [ $t -eq $nFiles ]
            then
            n=0
            echo $ParFiles
            hadd AnalysisResults$t.root $PartFiles
            PartFiles=" "
            FullFiles=" $FullFiles AnalysisResults$t.root "
            fi
       done
       hadd AnalysisResults${i}.root $FullFiles
       mv AnalysisResults${i}.root ../
       cd ..
       mv alice tobedeleted${i}
    done
fi


