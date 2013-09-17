#!/bin/sh
#small script to download files from alien

#
#needed info: e.g. /alice/data/2013/LHC13f/000197258/pass2/event_stat.root
#
#dataPeriodFolder=/alice/data/2013/LHC13c/
#productionFileName=/ESDs/pass2/event_stat.root
#destinationFolder=`pwd`/LHC13c/
dataPeriodFolder=$1
productionFileName=$2
destinationFolder=$3


#
# prepare directories
#
mkdir $destinationFolder

#
# copy the files
#

counter=0
for i in `alien_ls $dataPeriodFolder`; do
    counter=`expr $counter + 1`;
    runNumber=`expr $i |  sed 's/\/000/\//g'`
    echo Copying file number $counter which corresponds to $runNumber to $destinationFolder$runNumber\_event_stat.root
    alien_cp  alien://$dataPeriodFolder$i$productionFileName file://$destinationFolder$runNumber\_event_stat.root
done
