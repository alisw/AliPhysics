#!/bin/bash

# Script to merge objects coming out of the calibration train and export of the OCDB:
# Arguments for local use:
#    1 - directory on which to look for the files to be merged 
#    2 - runNumber
#    3 - OCDB output path
#    optional:
#    4 - input default OCDB, default="raw://"
#    5 - option: "alien_cp" will use alien_cp to copy files from alien
#                "tfilecp" will copy the files first using TFile::Cp()
#                "nocopy" will access files over the network directly TFile::Open()
#                default is tfilecp, use 0
#    6 - number of files to download per iteration, default 20
#
# example local use (files will be downloaded from alien)
#./mergeMakeOCDB.byComponent.sh /alice/data/2012/LHC12g/000188362/cpass0/ 188362 local://./OCDB local:///cvmfs/alice.gsi.de/alice/data/2011/OCDB/ 0 30
#
# on ALIEN just do:
# $1 = directory where to perform the find 
# $2 = runNumber
# $3 = OCDB path

if [[ $# -eq 0 ]]; then
  echo arguments:
  echo  "  1 - directory on which to look for the files to be merged or local file list"
  echo  "  2 - runNumber"
  echo  "  3 - OCDB output path"
  echo  "  optional:"
  echo  "  4 - input default OCDB, default=\"raw://\""
  echo  "  5 - option: \"alien_cp\" will use alien_cp to copy files from alien"
  echo  "              \"tfilecp\" will copy the files first using TFile::Cp()"
  echo  "              \"nocopy\" will access files over the network directly TFile::Open()"
  echo  "              default is tfilecp (use 0)"
  echo  "  6 - number of files to download/merge per iteration, default 20"
  echo
  echo "example:"
  echo  " ./mergeMakeOCDB.byComponent.sh /alice/data/2012/LHC12g/000188362/cpass1/ 188362 local://./OCDB raw:// alien_cp 10"
  echo

  exit 0
fi

# init
path=$1
run=$(echo "$2" | sed 's/^0*//')
ocdb=$3
defaultOCDB="raw://"
[[ -n $4 ]] && defaultOCDB=$4
fileAccessMethod=$5
[[ "$fileAccessMethod" != "alien_cp" && "$fileAccessMethod" != "tfilecp" && "$fileAccessMethod" != "nocopy" ]] && fileAccessMethod="tfilecp"
numberOfFilesInAbunch=50
[[ $6 =~ ^[0-9]+$ ]] && numberOfFilesInAbunch=$6

isLocal=0
[[ -f $path ]] && isLocal=1
cleanup=1
[[ $isLocal -eq 1 ]] && cleanup=0
[[ $isLocal -eq 1 ]] && fileAccessMethod="nocopy"

# setup components to be merged
components="TOF MeanVertex T0 SDD TRD TPCCalib TPCAlign TPCCluster"

#################################################################
echo
echo $0" $*"
echo
echo "***********************" | tee -a merge.log
echo mergeMakeOCDB.byComponent.sh started | tee -a merge.log
echo path = $path | tee -a merge.log
echo run  = $run | tee -a merge.log
echo ocdb = $ocdb | tee -a merge.log
echo defaultOCDB=$defaultOCDB | tee -a merge.log
echo isLocal = $isLocal | tee -a merge.log
echo cleanup = $cleanup | tee -a merge.log
echo fileAccessMethod=$fileAccessMethod | tee -a merge.log
echo numberOfFilesInAbunch = $numberOfFilesInAbunch | tee -a merge.log
echo "***********************" | tee -a merge.log

alienFileList="alien.list"
localFileList="local.list"
partialLocalFileListPrefix=${localFileList}_
partialAlienFileListPrefix=${alienFileList}_
runningMergeByComponentLockFile="runningMergeByComponent.lock"

mergeByComponent()
{
  # process by component
  # first argument is the file list to process

  #lock
  touch $runningMergeByComponentLockFile

  # run inside a dedicated running directory
  # whic means copy the file list to process and prefic each line with ../
  # since the file names have no absolute paths!
  runningDirectory="${runningMergeByComponentLockFile}.${1}.dir"
  fileList="$1"
  cleanup=$2
  
  #sanity checks
  [[ ! -f $fileList ]] && echo "$fileList does not exist" && return
  mkdir -p $runningDirectory
  [[ ! -d $runningDirectory ]] && echo "cannot create the running directory $runningDirectory" && return 1
  
  #add the results of previous iteration to the BEGINNING of the list of files to be merged
  if [[ -f CalibObjects.root ]]; then
    echo "../CalibObjects.root" > $runningDirectory/$fileList
  fi

  #move the to be merged files to the running directory and make a list of available files
  #handle the case of archives (x.zip#y.root) as well
  nFiles=0
  while read entry; do
    if [[ $entry =~ ^.*\.root$ ]]; then
      file=${entry%#*}
      fileContent=${entry##*#}
      [[ "${fileContent}" == "${file}" ]] && fileContent=""
      [[ -n ${fileContent} ]] && fileContent="#${fileContent}"
      if [[ -f ${file} ]]; then
        ((nFiles++))
        if [[ -f "./${file}" ]]; then
          echo "../${file}${fileContent}" >> "${runningDirectory}/$fileList"
        else
          echo "${file}${fileContent}" >> "${runningDirectory}/$fileList"
        fi
      fi
    fi
  done < $fileList
  [[ $nFiles -lt 1 ]] && echo "no new files in $fileList" && rm -rf $runningDirectory && return

  #copy the macro to the running directory
  [[ -f mergeByComponent.C ]] && cp mergeByComponent.C $runningDirectory

  cd $runningDirectory

  echo "processing following files from $fileList:"
  cat $fileList

  for det in $components; do
    # merge
    echo "***********************" 
    echo merging $det data
    echo "***********************"
    echo aliroot -b -q "mergeByComponent.C(\"$det\", \"$fileList\")"
    aliroot -b -q "mergeByComponent.C(\"$det\", \"$fileList\")" 2>&1 | tee -a merge_${det}.log
    if validateMerging ${det}; then
      echo "### merge OK: mv CalibObjects.root CalibObjects_$det.root"
      mv CalibObjects.root CalibObjects_$det.root
    else 
    echo "### merging not validated"
      rm -f CalibObjects.root
    fi
    mv syswatch.log syswatch_merge_$det.log
  done
  
  # global merge
  echo "***********************"
  echo merging ALL data
  echo "***********************"
  partialCalibObjectsList="objects.list.${1}"
  ls -1 CalibObjects_*.root > $partialCalibObjectsList
  echo aliroot -b -q "mergeByComponent.C(\"ALL\", \"$partialCalibObjectsList\")"
  aliroot -b -q "mergeByComponent.C(\"ALL\", \"$partialCalibObjectsList\")" 2>&1 | tee -a merge_ALL.log
  if validateMerging "ALL"; then
    #move the new CalibObjects to parent directory
    echo mv -f CalibObjects.root ..
    mv -f CalibObjects.root ..
  else 
    echo "### merging not validated"
    rm -f CalibObjects.root
  fi

  mv syswatch.log syswatch_ALL.log
  rm $partialCalibObjectsList
  
  #move stuff back to the parent dir and clean up
  #merge the syswatch logs
  for x in syswatch*log; do
    [[ ! -f $x ]] && continue
    if [[ -f ../$x ]]
    then 
      echo "sed '1d' $x  >> ../$x"
      sed '1d' $x >> ../$x
      rm -f $x
    else 
      echo mv -f $x ..
      mv -f $x ..
    fi
  done

  #merge the other logs
  for x in *.log; do
    [[ ! -f $x ]] && continue
    if [[ -f ../$x ]]
    then
      echo "cat $x >> ../$x"
      cat $x >> ../$x
    else
      echo "mv -f $x .."
      mv -f $x ..
    fi
  done
  
  #final cleanup.
  #go to parent dir first to use the original fileList
  #without ../CalibObjects.root
  cd ..
  if [[ $cleanup -eq 1 ]]; then
    echo "cleaning up merged files..."
    while read file; do
      echo rm -rf $file
      rm -rf $file
    done<$fileList
  fi
  rm -rf $runningDirectory

  #unlock
  rm -f $runningMergeByComponentLockFile
  [[ ! -f CalibObjects.root ]] && echo "WARNING: CalibObjects.root not available!"
  echo "***mergeByComponent() DONE"
  echo
  return 0
}

waitIfLocked()
{
  while [ 1 ]; do
    [[ ! -f $1 ]] && break
    sleep 1
  done
  return 0
}

validateMerging()
{
  det=$1
  retCode=0
  [[ ! -f CalibObjects.root ]] && ((retCode++)) && echo "### no CalibObjects.root..."
  [[ ! -f ${det}_merge_done ]] && ((retCode++)) && echo "### no ${det}_merge_done, job finished abnormally..."
  error=$(grep -e "was a crash" *.log)
  [[ -n $error ]] && ((retCode++)) && echo "### error: $error"
  
  return $retCode
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

#first, make sure we have the scripts
copyScripts | tee -a copy.log
ls
#with alien files copy them first to local
echo "***********************" 2>&1 | tee -a copy.log
echo copying files for run $run 2>&1 | tee -a copy.log
echo from $path 2>&1 | tee -a copy.log
echo "***********************" 2>&1 | tee -a copy.log
if [[ $isLocal -eq 0 ]]; then
  if [[ "$fileAccessMethod" == "alien_cp" ]]; then
    echo "alien_find $path AliESDfriends_v1.root | egrep ^/ > $alienFileList" 2>&1 | tee -a copy.log
    alien_find $path "AliESDfriends_v1.root" | egrep "^/" >  $alienFileList
    echo "alien_find done"
    echo
  else 
    echo aliroot -b -q "mergeByComponent.C(\"MAKEALIENLIST\",\"$alienFileList\", \"$path\", \"AliESDfriends_v1.root\")" 2>&1 | tee -a copy.log
    aliroot -b -q "mergeByComponent.C(\"MAKEALIENLIST\",\"$alienFileList\", \"$path\", \"AliESDfriends_v1.root\")" 2>&1 | tee -a copy.log
    echo "MAKEALIENLIST done"
    echo
  fi
else
  cp $path $alienFileList
fi
echo split --numeric-suffixes --suffix-length=6 --lines=$numberOfFilesInAbunch ${alienFileList} ${partialAlienFileListPrefix} | tee -a copy.log
split --numeric-suffixes --suffix-length=6 --lines=$numberOfFilesInAbunch ${alienFileList} ${partialAlienFileListPrefix}
rm -f $runningMergeByComponentLockFile
for partialAlienFileList in ${partialAlienFileListPrefix}*
do
  partialAlienFileListPostfix=${partialAlienFileList#$partialAlienFileListPrefix}
  partialLocalFileList=${partialLocalFileListPrefix}${partialAlienFileListPostfix}

  #copy the files if appropriate and make a list of files to use
  rm -f $partialLocalFileList
  if [[ "$fileAccessMethod" == "alien_cp" ]]; then
    while read x; do
      [[ $x != /* ]] && continue
      localName=${x//"/"/_}
      echo "alien_cp "alien://$x" $localName" 2>&1 | tee -a copy.log
      alien_cp "alien://$x" $localName
      echo $localName>>$partialLocalFileList
    done<$partialAlienFileList
  elif [[ "$fileAccessMethod" == "tfilecp" ]]; then
    echo aliroot -b -q "mergeByComponent.C(\"COPY\",\"$partialAlienFileList\",\"noPath\",\"noPattern\",10,\"$partialLocalFileList\")" 2>&1 | tee -a copy.log
    aliroot -b -q "mergeByComponent.C(\"COPY\",\"$partialAlienFileList\",\"noPath\",\"noPattern\",10,\"$partialLocalFileList\")" 2>&1 | tee -a copy.log
  elif [[ "$fileAccessMethod" == "nocopy" ]]; then
    while read x; do
      [[ $isLocal -eq 0 ]] && echo "alien://$x" >> $partialLocalFileList
      [[ $isLocal -eq 1 ]] && echo "$x" >> $partialLocalFileList
    done<$partialAlienFileList
    cleanup=0
  fi

  #handle syswatch
  if [[ -f syswatch_copy.log ]]
  then
    sed '1d' syswatch.log >> syswatch_copy.log
    rm -f syswatch.log
  else 
    [[ -f syswatch.log ]] && mv syswatch.log syswatch_copy.log
  fi

  #merge in parallel, use a simple lockfile
  #echo "waitIfLocked $runningMergeByComponentLockFile" | tee -a merge.log
  waitIfLocked $runningMergeByComponentLockFile
  if [[ $isLocal -eq 1 ]]; then
    #locally we don't need to run in parallel since we dont download,
    #besides it causes problems on some filesystems (like lustre) with file sync
    mergeByComponent $partialLocalFileList $cleanup 2>&1 | tee -a merge.log
  else
    mergeByComponent $partialLocalFileList $cleanup 2>&1 | tee -a merge.log &
  fi
done

#cleanup
rm -f ${partialAlienFileListPrefix}*
rm -f ${partialLocalFileListPrefix}*
rm -f $alienFileList

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
