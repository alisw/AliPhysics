#available variables:
#  $dataType     e.g. data or sim
#  $year         e.g. 2011
#  $period       e.g. LHC13g
#  $runNumber     e.g. 169123
#  $pass         e.g. cpass1,pass1,passMC
#  #ocdbStorage  e.g. "raw://", "local://./OCDB"

runLevelQA()
{
  #full path of QAresults.root is provided
  local qaFile=$1

  # This is used only to extract the muon information
  ln -s $ALICE_ROOT/PWGPP/MUON/lite/LoadLibsForMuonQA.C
  ln -s $ALICE_ROOT/PWGPP/MUON/lite/MakeTrend.C
  aliroot -b -l <<EOF
.x LoadLibsForMuonQA.C+("maketrend")
.x MakeTrend.C("${qaFile}",${runNumber})
.q
EOF
  rm LoadLibsForMuonQA.C
  rm MakeTrend.C
  rm *.d *.so
}

periodLevelQA()
{
  #path of the merged period trending.root is provided
  local trendingFile=$1

  local fileList="trendList.txt"
  local tmpFileList="tmp${fileList}"
  local fileNames="QAresults.root QAresults_barrel.root QAresults_outer.root"
  find -L . -name "QAresults*.root" > $fileList

  # If no result is found, it means we have an archive as input
  if [ ! -s $fileList ]; then
    sFileList=$(find -L . -name "QA*archive.zip")
    for ifile in $sFileList; do
      for searchFile in $fileNames; do
        if unzip -l ${ifile} | egrep "$searchFile" &>/dev/null; then
          echo "${ifile}#${searchFile}" >> $fileList
        fi
      done
    done
  fi

  # Assume that "outputDir" is known from the steering runQA.sh
  local cfgFileDir="${outputDir}/configFiles"
  local cfgFileSuffix="${dataType}_${period}.txt"

  #if run list is provided, filter the output limiting to this list
  # FIXME: the code is run in a temporary directory
  # where should we add this file?
  local runList="${cfgFileDir}/runList_${cfgFileSuffix}"
  if [ -e ${runList} ]; then
    sRunList=$(cat ${runList} | xargs)
    tmpFileList="tmp${fileList}"
    mv ${fileList} ${tmpFileList}
    for irun in $sRunList; do
      currFile=$(grep ${irun} $tmpFileList)
      if [ "${currFile}" = "" ]; then
        continue
      fi
      echo "${currFile}" >> ${fileList}
    done
    rm $tmpFileList
  fi

  #if trigger list is provided, filter the tracking output accordngly
  # FIXME: the code is run in a temporary directory
  # where should we add this file?
  local triggerList="$cfgFileDir/trigList_${cfgFileSuffix}"
  if [ -e ${triggerList} ]; then
    triggerList="\"${triggerList}\""
  else
    triggerList="0x0"
  fi

  ln -s $ALICE_ROOT/PWGPP/MUON/lite/LoadLibsForMuonQA.C

  # First run tracker (it merges the QAresults and we need it for
  # scaler trending in trigger
  local mergedQAname="QAresults.root"
  for ifile in $fileNames; do
    hasFile=`grep -c "$ifile" $fileList`
    if [ $hasFile -gt 0 ]; then
      mergedQAname="$ifile"
      break
    fi
  done

  ln -s $ALICE_ROOT/PWGPP/MUON/lite/PlotMuonQA.C
aliroot -b <<EOF
.x LoadLibsForMuonQA.C+("tracktrend");
.x PlotMuonQA.C+(".","${fileList}",${triggerList},kTRUE,"muon_tracker","${mergedQAname}");
.q
EOF
  rm PlotMuonQA.C

  # Then run trigger
  local runScalers="kFALSE"
  if [ "${dataType}" = "data" ]; then
    runScalers="kTRUE";
  fi

  ln -s $ALICE_ROOT/PWGPP/MUON/lite/trigEffQA.C
  aliroot -b <<EOF
.x LoadLibsForMuonQA.C+("trigtrend");
.x trigEffQA.C+("${fileList}","QA_muon_trigger.root","${ocdbStorage}",$runScalers,"${mergedQAname}");
.q
EOF
  rm trigEffQA.C

  rm LoadLibsForMuonQA.C
  rm *.d *.so
  rm ${mergedQAname} # remove merged file, since it is recreated each time

  rm $fileList;
}

#########################################################
#########EXPERTS ONLY####################################
#runLevelHighPtTreeQA()
#{
#  #input is the high pt tree (if available)
#  highPtTree=$1
#}
