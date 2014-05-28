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
  qaFile=$1

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

  #should produce a file trending.root
  #if not, a default one will be provided
}

periodLevelQA()
{
  #path of the merged period trending.root is provided
  trendingFile=$1

  fileList="trendList.txt"
  find -L . -name "QAresults.root" > $fileList

  # Assume that "outputDir" is known from the steering runQA.sh
  cfgFileDir="${outputDir}/configFiles"
  cfgFileSuffix="${dataType}_${period}.txt"

  #if run list is provided, filter the output limiting to this list
  # FIXME: the code is run in a temporary directory
  # where should we add this file?
  runList="${cfgFileDir}/runList_${cfgFileSuffix}"
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
  triggerList="$cfgFileDir/trigList_${cfgFileSuffix}"
  if [ -e ${triggerList} ]; then
    triggerList="\"${triggerList}\""
  else
    triggerList="0x0"
  fi

  ln -s $ALICE_ROOT/PWGPP/MUON/lite/LoadLibsForMuonQA.C

  # First run tracker (we it merges the QAresults and we need it for
  # scaler trending in trigger
  mergedQAname="QAresults.root"

  ln -s $ALICE_ROOT/PWGPP/MUON/lite/PlotMuonQA.C
aliroot -b <<EOF
.x LoadLibsForMuonQA.C+("tracktrend");
.x PlotMuonQA.C+(".","${fileList}",${triggerList},kTRUE,"muon_tracker","${mergedQAname}");
.q
EOF
  rm PlotMuonQA.C

  # Then run trigger
  runScalers="kFALSE"
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
