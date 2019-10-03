#! /bin/bash

# run it:
# ./copyfiles 2016 LHC16l muon_calo_pass1 752_20170220-2207 runlist.txt

# year=2016
# period=LHC16l
# pass=muon_calo_pass1
# train=752_20170220-2207
# list=runlist.txt

year=$1
period=$2
pass=$3
train=$4
list=$5

for i in `cat ${list}`
do
    
    runnumber=`basename --suffix=, $i`

    if [ ! -d ${runnumber} ]
    then
	mkdir ${runnumber}
    fi

    alien_cp -m -t 100 alien:///alice/data/${year}/${period}/000${runnumber}/${pass}/PWGPP/PP_EMCAL_Calibration/${train}/AnalysisResults.root ./${runnumber}/

done
    
