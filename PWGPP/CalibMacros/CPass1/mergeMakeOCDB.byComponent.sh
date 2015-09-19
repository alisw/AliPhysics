#!/bin/bash
# example local use (files will be downloaded from alien)
#./mergeMakeOCDB.byComponent.sh /alice/data/2012/LHC12g/000188362/cpass0/ 188362 local://./OCDB defaultOCDB=local:///cvmfs/alice.gsi.de/alice/data/2011/OCDB/ fileAccessMethod=tfilecp numberOfFilesInAbunch=30
#
# on ALIEN just do:
# $1 = directory where to perform the find 
# $2 = runNumber
# $3 = OCDB path

main()
{
  maxChunksTOT=999999

  save_args=("$@") 
  if [[ $# -eq 0 ]]; then
    echo arguments:
    echo  "  1 - directory on which to look for the files to be merged or local file list"
    echo  "  2 - runNumber"
    echo  "  3 - OCDB output path"
    echo
    echo  "  optionally set any of these"
    echo  "  maxChunksTOT=<max number of chunks> (def: ${maxChunksTOT}, max number of chunks to process)"
    echo  "                                       random selection will be done but 1st entr)"
    echo  "  runParallel={0,1}   (1 to run parallel merging and downloading)"
    echo  "  defaultOCDB=raw://  (or any other valid OCDB storage as input database)"
    echo  "  fileAccessMethod={alien_cp,tfilecp,nocopy}   (by default tfilecp)"
    echo  "                    \"alien_cp\" will use alien_cp to copy files from alien"
    echo  "                    \"tfilecp\" will copy the files first using TFile::Cp()"
    echo  "                    \"nocopy\" will access files over the network directly TFile::Open()"
    echo  "  numberOfFilesInAbunch=50                     (how many files get downloaded in one go)"
    echo  "  maxChunksTPC=3000    (max number of chunks to be merged for TPC calibration)"

    echo
    echo "example:"
    echo  " ./mergeMakeOCDB.byComponent.sh /alice/data/2012/LHC12g/000188362/cpass1/ 188362 local://./OCDB runParallel=0 defaultOCDB=raw:// fileAccessMethod=alien_cp"
    echo

    return 0
  fi

  # init
  path=$1
  run=$(echo "$2" | sed 's/^0*//')
  ocdb=$3

  shift 3

  #default values
  runParallel=0
  defaultOCDB="raw://"
  fileAccessMethod="tfilecp"
  numberOfFilesInAbunch=50
  maxChunksTPC=3000
  maxChunksTOT=3000
  [[ -n ${ALIEN_JDL_TTL} ]] && maxTimeToLive=$(( ${ALIEN_JDL_TTL}-2000 ))

  parseConfig "$@"

  #some sanity checks of option values
  [[ ! $numberOfFilesInAbunch =~ ^[0-9]+$ ]] && numberOfFilesInAbunch=50
  [[ "$fileAccessMethod" != "alien_cp" && "$fileAccessMethod" != "tfilecp" && "$fileAccessMethod" != "nocopy" ]] && fileAccessMethod="tfilecp" && echo "using default fileAccessMethod=tfilecp"

  filesAreLocal=0
  [[ -f $path ]] && filesAreLocal=1
  cleanup=1
  [[ $filesAreLocal -eq 1 ]] && cleanup=0
  [[ $filesAreLocal -eq 1 ]] && fileAccessMethod="nocopy"

  # setup components to be merged
  components="TOF MeanVertex T0 SDD TRD TPCCalib TPCCluster TPCAlign"
  #components="TOF MeanVertex T0 SDD TRD TPCCalib"

  # take Data Quality Flags from JDL
  detectorBitsQualityFlag=${ALIEN_JDL_RUNQUALITY--1}

  #################################################################
  echo "" | tee -a merge.log
  echo "$0 ${save_args[*]}" | tee -a merge.log
  echo "" | tee -a merge.log
  echo "***********************" | tee -a merge.log
  echo mergeMakeOCDB.byComponent.sh started | tee -a merge.log
  echo path = $path | tee -a merge.log
  echo run  = $run | tee -a merge.log
  echo ocdb = $ocdb | tee -a merge.log
  echo defaultOCDB = $defaultOCDB | tee -a merge.log
  echo filesAreLocal = $filesAreLocal | tee -a merge.log
  echo cleanup = $cleanup | tee -a merge.log
  echo fileAccessMethod = $fileAccessMethod | tee -a merge.log
  echo numberOfFilesInAbunch = $numberOfFilesInAbunch | tee -a merge.log
  echo runParallel = $runParallel
  echo detectorBitsQualityFlag = $detectorBitsQualityFlag
  echo "***********************" | tee -a merge.log

  alienFileList="alien.list"
  localFileList="local.list"
  filesProcessedTPClist="filesProcessedTPC.log"
  rm -f $filesProcessedTPClist
  touch $filesProcessedTPClist
  partialLocalFileListPrefix=${localFileList}_
  partialAlienFileListPrefix=${alienFileList}_

  #first, make sure we have the scripts
  copyScripts | tee -a copy.log
  ls
  #with alien files copy them first to local
  echo "***********************" 2>&1 | tee -a copy.log
  echo copying files for run $run 2>&1 | tee -a copy.log
  echo from $path 2>&1 | tee -a copy.log
  echo "***********************" 2>&1 | tee -a copy.log
  if [[ $filesAreLocal -eq 0 ]]; then
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
  #randomize the list
  #keep the first line intact (it is the largest file of the entire collection)
  echo "randomizing $maxChunksTOT chunks of the list, preserving the 1st one "
  nchAdd=$[maxChunksTOT-1]
  head -n1 ${alienFileList} > ${alienFileList}.tmp
  tail -n +2 ${alienFileList} | shuf -n $nchAdd >> ${alienFileList}.tmp 
  mv -f ${alienFileList}.tmp ${alienFileList}

  #split into sublists, each to be processed separately
  echo split --numeric-suffixes --suffix-length=6 --lines=$numberOfFilesInAbunch ${alienFileList} ${partialAlienFileListPrefix} | tee -a copy.log
  #split --numeric-suffixes --suffix-length=6 --lines=$numberOfFilesInAbunch ${alienFileList} ${partialAlienFileListPrefix}
  split -a 6 -l $numberOfFilesInAbunch ${alienFileList} ${partialAlienFileListPrefix}

  for partialAlienFileList in ${partialAlienFileListPrefix}*
  do

    #if it takes too long, break
    #[[ ${SECONDS} -gt ${maxTimeToLive} ]] && break

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
        [[ $filesAreLocal -eq 0 ]] && echo "alien://$x" >> $partialLocalFileList
        [[ $filesAreLocal -eq 1 ]] && echo "$x" >> $partialLocalFileList
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

    echo waiting
    wait $!
    if [[ $runParallel -eq 1 ]]; then
      mergeByComponent $partialLocalFileList $cleanup 2>&1 | tee -a merge.log &
    else
      mergeByComponent $partialLocalFileList $cleanup 2>&1 | tee -a merge.log
    fi
  done

  #merge all the subfiles into one, wait for the last one to complete
  echo waiting
  wait $!
  if [[ "$components" =~ ALL && -f CalibObjects_ALL.root ]]; then
    mv -f CalibObjects_ALL.root CalibObjects.root
  else
    echo "***********************"
    echo merging ALL data
    echo "***********************"
    finalCalibObjectsList="finalObjects.list"
    ls -1 CalibObjects_*.root > $finalCalibObjectsList
    echo aliroot -b -q "mergeByComponent.C(\"ALL\", \"$finalCalibObjectsList\")" | tee -a merge.log
    aliroot -b -q "mergeByComponent.C(\"ALL\", \"$finalCalibObjectsList\")" 2>&1 | tee -a merge.log
    mv -f syswatch.log syswatch_merge_ALL.log
  fi

  if ! validateMerging "ALL"; then
    echo final merging not validatet, exiting...
    return 1
  fi
  rm -f CalibObjects_*.root

  #cleanup
  rm -f ${partialAlienFileListPrefix}*
  rm -f ${partialLocalFileListPrefix}*
  rm -f $alienFileList
  rm -f $localFileList

  # make OCDB
  echo "***********************" 2>&1 | tee -a ocdb.log
  echo making ${det} OCDB 2>&1 | tee -a ocdb.log
  echo "***********************" 2>&1 | tee -a ocdb.log
  echo aliroot -b -q "makeOCDB.C($run, \"$ocdb\", \"$defaultOCDB\", $detectorBitsQualityFlag)" 2>&1 | tee -a ocdb.log
  aliroot -b -q "makeOCDB.C($run, \"$ocdb\", \"$defaultOCDB\", $detectorBitsQualityFlag)" 2>&1 | tee -a ocdb.log
  mv syswatch.log syswatch_makeOCDB.log

  # summary
  echo "***********************" 2>&1 | tee -a ocdb.log
  echo SUMMARY 2>&1 | tee -a ocdb.log
  echo "***********************" 2>&1 | tee -a ocdb.log
  ls -altr CalibObjects.root *done 2>&1 | tee -a ocdb.log
}

mergeByComponent()
{
  # process by component
  # first argument is the file list to process

  # run inside a dedicated running directory
  # whic means copy the file list to process and prefic each line with ../
  # since the file names have no absolute paths!
  runningDirectory="${1}.dir"
  fileList="$1"
  cleanup=$2

  #sanity checks
  if [[ ! -f ${fileList} ]]; then
    echo "${fileList} does not exist"
    return 1
  fi
  mkdir -p $runningDirectory
  if [[ ! -d $runningDirectory ]]; then
    echo "cannot create the running directory $runningDirectory"
    return 1
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
          echo "../${file}${fileContent}" >> "${runningDirectory}/${fileList}"
        else
          echo "${file}${fileContent}" >> "${runningDirectory}/${fileList}"
        fi
      fi
    fi
  done < ${fileList}
  if [[ $nFiles -lt 1 ]]; then
    echo "no new files in ${fileList}"
    echo rm -rf $runningDirectory 
    rm -rf $runningDirectory 
    return 1
  fi

  #copy the macro to the running directory
  [[ -f mergeByComponent.C ]] && cp mergeByComponent.C $runningDirectory

  #go to running directory
  cd $runningDirectory

  numberOfChunksTPC=$(cat ../$filesProcessedTPClist 2>/dev/null | wc -l )

  for det in $components; do

    # merge
    echo "***********************" 
    echo merging ${det} data
    echo "***********************"

    #limit the number of chunks processed by TPC
    if [[ "${det}" =~ TPC && $numberOfChunksTPC -ge $maxChunksTPC ]]; then
      echo "Not merging TPC anymore, max number of chunks processed ($maxChunksTPC)"
      continue
    fi

    #add the results of previous iteration to the BEGINNING of the list of files to be merged
    [[ -f ../CalibObjects_${det}.root ]] && echo "../CalibObjects_${det}.root" > ${fileList}_${det}
    cat ${fileList} >> ${fileList}_${det}

    echo "processing following files from ${fileList}_${det}:"
    cat ${fileList}_${det}

    echo aliroot -b -q "mergeByComponent.C(\"${det}\", \"${fileList}_${det}\")"
    aliroot -b -q "mergeByComponent.C(\"${det}\", \"${fileList}_${det}\")" 2>&1 | tee -a merge_${det}.log
    if validateMerging ${det}; then
      echo "### merge OK: mv CalibObjects.root ../CalibObjects_${det}.root"
      mv -f CalibObjects.root ../CalibObjects_${det}.root
      [[ "${det}" =~ TPCCalib ]] && cat ${fileList} >> ../$filesProcessedTPClist
    else 
      echo "### merging not validated"
      rm -f CalibObjects.root
    fi
    mv syswatch.log syswatch_merge_${det}.log

  done


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

  echo "***mergeByComponent() DONE"
  echo
  echo "numberOfChunksTPC=$numberOfChunksTPC"
  return 0
}

parseConfig()
{
  args=("$@")
  for opt in "${args[@]}"; do
    [[ -z ${opt} ]] && continue
    [[ -n ${encodedSpaces} ]] && opt="$(decSpaces ${opt})"
    [[ "${opt}" =~ ^[[:space:]]*$ ]] && continue
    if [[ ! "${opt}" =~ .*=.* ]]; then
      echo "badly formatted option \"${opt}\" should be: option=value, stopping..."
      return 1
    fi
    local var="${opt%%=*}"
    local value="${opt#*=}"
    export ${var}="${value}"
  done

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
    cp -f $ALICE_ROOT/PWGPP/CalibMacros/CPass1/mergeByComponent.C $PWD && \
    echo "taking the default scripts from $ALICE_ROOT"
  [[ ! -f makeOCDB.C ]] && \
    cp -f $ALICE_ROOT/PWGPP/CalibMacros/CPass1/makeOCDB.C $PWD && \
    echo "taking the default scripts from $ALICE_ROOT"
}

#these functions encode strings to and from a space-less form
#use when spaces are not well handled (e.g. in arguments to
#commands in makeflow files, etc.
encSpaces()(echo "${1// /@@@@}")
decSpaces()(echo "${1//@@@@/ }")

main "$@"
