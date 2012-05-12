#! /bin/bash

# Script to merge objects coming out of the calibration train and export of the OCDB:
# Arguments for local use:
#    1 - directory on which to look for the files to be merged 
#    2 - runNumber
#    3 - OCDB output path
#    [4 - input default OCDB]

#ALIEN setting
# $1 = directory where to perform the find 
# $2 = runNumber
# $3 = OCDB path

# if fourth argument given, its the default OCDB, otherwise use the default raw://
defaultOCDB=$4
[[ -z $4 ]] && defaultOCDB="raw://"

# init
path=$1
run=$2
ocdb=$3
isLocal=0
[[ -f $path ]] && isLocal=1

echo "***********************" 2>&1 | tee -a merge.log
echo mergeMakeOCDB.byComponent.sh started 2>&1 | tee -a merge.log
echo path = $path 2>&1 | tee -a merge.log
echo run  = $run 2>&1 | tee -a merge.log
echo ocdb = $ocdb 2>&1 | tee -a merge.log
echo defaultOCDB=$defaultOCDB | tee -a merge.log
echo isLocal = $isLocal 2>&1 | tee -a merge.log
echo "***********************" 2>&1 | tee -a merge.log

# setup components
components="TOF MeanVertex T0 SDD TRD TPCCalib  TPCAlign TPCCluster"

alienFileList="alien.list"
localFileList="local.list"
partialLocalFileListPrefix=${localFileList}_
partialAlienFileListPrefix=${alienFileList}_
runningMergeByComponentLockFile="runningMergeByComponent.lock"

mergeByComponent()
{
  # process by component
  # first argument is the file list to process
  # second argument can be "doCleanup" to delete merged files

  #lock
  touch $runningMergeByComponentLockFile

  previousResults=PREVIOUS_ITERATION_CalibObjects.root
  if [[ -f CalibObjects.root ]]; then
    mv CalibObjects.root $previousResults
    echo "$previousResults" >> $1
  fi

  echo "####DEBUG" | tee -a merge.log
  echo "####processed list $1" | tee -a merge.log
  cat $1

  for det in $components; do
    # merge
    echo "***********************" 2>&1 | tee -a merge.log
    echo merging $det data 2>&1 | tee -a merge.log
    echo "***********************" 2>&1 | tee -a merge.log
    echo aliroot -b -q "$ALICE_ROOT/PWGPP/CalibMacros/CPass0/mergeByComponent.C(\"$det\", \"$1\")" 2>&1 | tee -a merge.log
    aliroot -b -q "$ALICE_ROOT/PWGPP/CalibMacros/CPass0/mergeByComponent.C(\"$det\", \"$1\")" 2>&1 | tee -a merge.log
    mv syswatch.log syswatch_merge_$det.log
    mv CalibObjects.root CalibObjects_$det.root
  done
  
  rm -f $previousResults
  
  # global merge
  echo "***********************" 2>&1 | tee -a merge.log
  echo merging ALL data 2>&1 | tee -a merge.log
  echo "***********************" 2>&1 | tee -a merge.log
  partialCalibObjectsList="objects.list.${1}"
  ls -1 CalibObjects_*.root > $partialCalibObjectsList
  aliroot -b -q "$ALICE_ROOT/PWGPP/CalibMacros/CPass0/mergeByComponent.C(\"ALL\", \"$partialCalibObjectsList\")" 2>&1 | tee -a merge.log
  mv syswatch.log syswatch_ALL.log
  
  #cleanup
  if [[ "$2" == "doCleanup" ]]; then
    while read filename; do
      echo rm -f $filename | tee -a merge.log
      rm -f $filename
    done < $1
  fi
  rm -f CalibObjects_*.root

  #unlock
  rm -f $runningMergeByComponentLockFile
}

waitIfLocked()
{
  while [ 1 ]; do
    [[ ! -f $1 ]] && break
    sleep 1
  done
}

if [ $isLocal -eq 0 ]; then
    #with alien files copy them first to local
    echo "***********************" 2>&1 | tee -a merge.log
    echo copying files for run $run 2>&1 | tee -a merge.log
    echo from $path 2>&1 | tee -a merge.log
    echo "***********************" 2>&1 | tee -a merge.log
    aliroot -b -q "$ALICE_ROOT/PWGPP/CalibMacros/CPass0/mergeByComponent.C(\"MAKEALIENLIST\",\"$alienFileList\", \"$path\", \"AliESDfriends_v1.root\")" 2>&1 | tee -a merge.log
    split --numeric-suffixes --suffix-length=4 --lines=100 ${alienFileList} ${partialAlienFileListPrefix}
    counter=0
    rm -f $runningMergeByComponentLockFile
    for partialAlienFileList in ${partialAlienFileListPrefix}*
    do
      #copy the alien files to local
      partialAlienFileListPostfix=${partialAlienFileList#$partialAlienFileListPrefix}
      partialLocalFileList=${partialLocalFileListPrefix}${partialAlienFileListPostfix}
      aliroot -b -q "$ALICE_ROOT/PWGPP/CalibMacros/CPass0/mergeByComponent.C(\"COPY\",\"$partialAlienFileList\",\"noPath\",\"noPattern\",10,\"$partialLocalFileList\")" 2>&1 | tee -a merge.log
      [[ $counter -ne 0 ]] && sed -i '1d' syswatch.log && $(( $counter++ ))
      cat syswatch.log >> syswatch_copy.log
      rm syswatch.log

      #merge in parallel, use a simple lockfile
      waitIfLocked $runningMergeByComponentLockFile
      mergeByComponent $partialLocalFileList "doCleanup" &
    done
else
  #locally just use the list directly
  mergeByComponent $path
fi;

#cleanup
cat ${partialAlienFileListPrefix}* >> ${localFileList}
rm -f ${partialAlienFileListPrefix}*
rm -f ${partialLocalFileListPrefix}*

# make OCDB
echo "***********************" 2>&1 | tee -a ocdb.log
echo making $det OCDB 2>&1 | tee -a ocdb.log
echo "***********************" 2>&1 | tee -a ocdb.log
aliroot -b -q "$ALICE_ROOT/PWGPP/CalibMacros/CPass0/makeOCDB.C($run, \"$ocdb\", \"$defaultOCDB\")" 2>&1 | tee -a ocdb.log
mv syswatch.log syswatch_makeOCDB.log

# summary
echo "***********************" 2>&1 | tee -a ocdb.log
echo SUMMARY 2>&1 | tee -a ocdb.log
echo "***********************" 2>&1 | tee -a ocdb.log
ls -altr *CalibObjects.root *done 2>&1 | tee -a ocdb.log
