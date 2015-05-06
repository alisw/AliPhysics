#!/bin/bash

#
# Aim: (Filtered tree --> raw event list) 
# File structures:
#     input: root file  --> root file with "highPt, Laser, V0s, Cosmics" trees:
#                           tree -> Scan() <esdFilePath> <eventNumberinFile> <triggerName> <GID> <timeStamp> <blabla>
#     output:ASCII file -->  <rawFilePath> <eventNumberinFile> <triggerName> <GID> <timeStamp> <blabla>
#                       E.g.: 
#                      /hera/alice/alien/raw/alice/data/2010/LHC10h/000137161/raw/10000137161001.10.root 240 V0s 52242539733 1289261418 
#
# how to run:  
#            source $ALICE_PHYSICS/../src/PWGPP/rawmerge/makeEventList.sh
#            makeEventListAllRuns FilteredESDs.list 
#            makeEventListPerRun  FilteredESDs.root
#

DATE=`date +%Y%m%d_%H_%M`
qsubCommand="qsub -cwd -V -l h_rt=24:0:0,h_rss=4G -b y -r y -o out.log -e err.log"
codeDir=$ALICE_PHYSICS/../src/PWGPP/rawmerge/
source $codeDir/makeEventList.config

###############################################################################
makeEventListAllRuns()
{
 
  # input --> <path>/FilteredESDs.list  
  filteredTreeList=$1  
  outbase=`pwd`
  rawFilteringDir=$outbase/FilteredRawEventLists_$DATE; DirCheckCreate $rawFilteringDir; cd $rawFilteringDir
  
  # loop over files and run makeEventListPerRun for each
  fileIndex=0
  for dataFile in $(cat $filteredTreeList); do
  
    fileIndex=$(($fileIndex+1))
    resultDirName=$rawFilteringDir/run_$fileIndex; mkdir $resultDirName; cd $resultDirName
    echo $fileIndex $dataFile
    cp $codeDir/*.C $codeDir/*.sh  $codeDir/*.config $codeDir/*.jdl .
    
    # submit jobs
    if [ $localRunning == 1 ]; then
       echo "local"
       ./makeEventList.sh  makeEventListPerRun $dataFile 
    elif [ $localRunning == 0 ]; then
       echo "on farm"
       eval $qsubCommand ./makeEventList.sh  makeEventListPerRun $dataFile
    fi
    
  done

  cd $outbase

}

###########################################################################################################
makeEventListPerRun()
{ 
  chmod +x makeEventList.config; source makeEventList.config

  file=$1
  ## dump computer info 
  (source $ALICE_PHYSICS/../src/PWGPP/scripts/utilities.sh; hostInfo) 1>hostInfo.log

  echo "*********************************** Input Information *****************************"  
  echo ptMinHighPt =$ptMinHighPt
  echo ptMinV0s    =$ptMinV0s
  echo inputFile   =$file  
  echo PWD         =$PWD
  echo aliroot -b -q "${scriptPath}/makeEventList.C+(\"${file}\",${ptMinHighPt},${ptMinV0s})"
  echo "***********************************************************************************" 
  echo "------------------------- raw event list is being produced ------------------------" 

  # run the macro and filter the output of TTree::Scan 
  aliroot -b -q "$ALICE_ROOT/../src/PWGPP/rawmerge/makeEventList.C(\"${file}\",${ptMinHighPt},${ptMinV0s})"

}
###############################################################################
DirCheckCreate()
{

dirName=$1
if [ -d "$dirName" ]; then
echo " !!! Attention !!! " $dirName "  exist already delete and recreate  "
rm -rf $dirName 
fi
mkdir $dirName

}
###############################################################################
main()
{
  eval "$@"
}
main "$@"
