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

makeAbsolutePathsInList()
{
  #make sure the files in the list have absolute paths
  rm -f $2
  while read file; do
    readlink -f $file >> $2
  done < $1
}

mergeByComponent()
{
  # process by component
  # first argument is the file list to process
  # second argument can be "doCleanup" to delete merged files

  #lock
  touch $runningMergeByComponentLockFile

  # run inside a dedicated running directory
  # whic means copy the file list to process and prefic each line with ../
  # since the file names have no absolute paths!
  runningDirectory="${runningMergeByComponentLockFile}.${1}.dir"
  parentDirectory=$PWD
  fileList="$1.local"
  mkdir -p $runningDirectory
  makeAbsolutePathsInList $1 $runningDirectory/$fileList
  cd $runningDirectory

  #copy the macro from $parentDirectory 
  [[ -f $parentDirectory/mergeByComponent.C ]] && cp $parentDirectory/mergeByComponent.C $PWD

  previousResults=PREVIOUS_ITERATION_CalibObjects.root
  if [[ -f $parentDirectory/CalibObjects.root ]]; then
    mv -f $parentDirectory/CalibObjects.root $previousResults
    echo "$previousResults" >> $fileList
  fi

  for det in $components; do
    # merge
    echo "***********************" 
    echo merging $det data
    echo "***********************"
    echo aliroot -b -q "mergeByComponent.C(\"$det\", \"$fileList\")"
    aliroot -b -q "mergeByComponent.C(\"$det\", \"$fileList\")"
    mv syswatch.log syswatch_merge_$det.log
    mv CalibObjects.root CalibObjects_$det.root
  done
  
  # global merge
  echo "***********************"
  echo merging ALL data
  echo "***********************"
  partialCalibObjectsList="objects.list.${1}"
  ls -1 CalibObjects_*.root > $partialCalibObjectsList
  echo aliroot -b -q "mergeByComponent.C(\"ALL\", \"$partialCalibObjectsList\")"
  aliroot -b -q "mergeByComponent.C(\"ALL\", \"$partialCalibObjectsList\")"
  mv syswatch.log syswatch_ALL.log
  rm $partialCalibObjectsList
  
  #cleanup of the merged files
  if [[ "$2" == "doCleanup" ]]; then
    while read filename; do
      echo rm -f $filename 
      rm -f $filename
    done < $fileList
  fi
  rm -f $fileList
  rm -f CalibObjects_*.root

  #move stuff back to the parent dir and clean up
  #merge the syswatch logs
  for x in syswatch*log; do
    if [[ -f $parentDirectory/$x ]]
    then 
      echo "sed '1d' $x  >> $parentDirectory/$x"
      sed '1d' $x >> $parentDirectory/$x
      rm -f $x
    else 
      echo mv -f $x $parentDirectory/$x
      mv -f $x $parentDirectory/$x
    fi
  done
  #merge the other logs
  for x in *.log; do
    if [[ -f $parentDirectory/$x ]]
    then
      echo "cat $x >> $parentDirectory/$x"
      cat $x >> $parentDirectory/$x
    else
      echo mv -f $x $parentDirectory/$x
      mv -f $x $parentDirectory/$x
    fi
  done
  #move the CalibObjects
  echo mv -f CalibObjects.root $parentDirectory
  mv -f CalibObjects.root $parentDirectory

  #final cleanup
  cd $parentDirectory
  rm -rf $runningDirectory

  #unlock
  rm -f $runningMergeByComponentLockFile
  echo "***mergeByComponent() DONE"
}

waitIfLocked()
{
  while [ 1 ]; do
    [[ ! -f $1 ]] && break
    sleep 1
  done
}

copyScripts()
{
    [[ ! -f mergeByComponent.C ]] && \
      cp -f $ALICE_ROOT/PWGPP/CalibMacros/CPass0/mergeByComponent.C $PWD && \
      echo "taking the default scripts from $ALICE_ROOT"
    [[ ! -f makeOCDB.C ]] && \
      cp -f $ALICE_ROOT/PWGPP/CalibMacros/CPass0/makeOCDB.C $PWD && \
      echo "taking the default scripts from $ALICE_ROOT"
}

if [ $isLocal -eq 0 ]; then
    #first, make sure we have the scripts
    copyScripts | tee -a copy.log
    ls
    #with alien files copy them first to local
    echo "***********************" 2>&1 | tee -a copy.log
    echo copying files for run $run 2>&1 | tee -a copy.log
    echo from $path 2>&1 | tee -a copy.log
    echo "***********************" 2>&1 | tee -a copy.log
    echo aliroot -b -q "mergeByComponent.C(\"MAKEALIENLIST\",\"$alienFileList\", \"$path\", \"AliESDfriends_v1.root\")" 2>&1 | tee -a copy.log
    aliroot -b -q "mergeByComponent.C(\"MAKEALIENLIST\",\"$alienFileList\", \"$path\", \"AliESDfriends_v1.root\")" 2>&1 | tee -a copy.log
    echo split --numeric-suffixes --suffix-length=4 --lines=20 ${alienFileList} ${partialAlienFileListPrefix} | tee -a copy.log
    split --numeric-suffixes --suffix-length=4 --lines=20 ${alienFileList} ${partialAlienFileListPrefix}
    rm -f $runningMergeByComponentLockFile
    for partialAlienFileList in ${partialAlienFileListPrefix}*
    do
      #copy the alien files to local
      partialAlienFileListPostfix=${partialAlienFileList#$partialAlienFileListPrefix}
      partialLocalFileList=${partialLocalFileListPrefix}${partialAlienFileListPostfix}
      echo aliroot -b -q "mergeByComponent.C(\"COPY\",\"$partialAlienFileList\",\"noPath\",\"noPattern\",10,\"$partialLocalFileList\")" 2>&1 | tee -a copy.log
      aliroot -b -q "mergeByComponent.C(\"COPY\",\"$partialAlienFileList\",\"noPath\",\"noPattern\",10,\"$partialLocalFileList\")" 2>&1 | tee -a copy.log
      
      #handle syswatch
      if [[ -f syswatch_copy.log ]]
      then
        sed '1d' syswatch.log >> syswatch_copy.log
        rm -f syswatch.log
      else 
        mv syswatch.log syswatch_copy.log
      fi

      #merge in parallel, use a simple lockfile
      echo "waitIfLocked $runningMergeByComponentLockFile" | tee -a merge.log
      waitIfLocked $runningMergeByComponentLockFile
      mergeByComponent $partialLocalFileList "doCleanup" 2>&1 | tee -a merge.log &
    done
else
  #locally just use the list directly
  mergeByComponent $path
fi;

#cleanup
rm -f ${partialAlienFileListPrefix}*
rm -f ${partialLocalFileListPrefix}*

# make OCDB
echo "***********************" 2>&1 | tee -a ocdb.log
echo making $det OCDB 2>&1 | tee -a ocdb.log
echo "***********************" 2>&1 | tee -a ocdb.log
echo aliroot -b -q "makeOCDB.C($run, \"$ocdb\", \"$defaultOCDB\")" 2>&1 | tee -a ocdb.log
aliroot -b -q "makeOCDB.C($run, \"$ocdb\", \"$defaultOCDB\")" 2>&1 | tee -a ocdb.log
mv syswatch.log syswatch_makeOCDB.log

# summary
echo "***********************" 2>&1 | tee -a ocdb.log
echo SUMMARY 2>&1 | tee -a ocdb.log
echo "***********************" 2>&1 | tee -a ocdb.log
ls -altr *CalibObjects.root *done 2>&1 | tee -a ocdb.log
