#available variables:
#  $dataType     e.g. data or sim
#  $year         e.g. 2011
#  $period       e.g. LHC13g
#  $runNumber     e.g. 169123
#  $pass         e.g. cpass1,pass1,passMC
#  #ocdbStorage  e.g. "raw://", "local://./OCDB"

ParseMUconfig()
{
  # Parse muon configuration file

  local cfgFilename="$ALICE_PHYSICS/PWGPP/MUON/lite/MU.config"
  local what="$1"

  awk -v swhat="$what" -v cond="${dataType}:${period}:${pass}" '{
    if ( index($0,cond) ) {
      while ( getline ) {
        if ( index($0,"type:period:pass") ) break;
        if ( index($0,swhat) ) {
          matchLine=$0
          gsub(swhat,"",matchLine)
          gsub("=","",matchLine)
          gsub("{","",matchLine)
          if ( match(matchLine, /[a-zA-Z0-9]/) ) print matchLine
          while ( getline ) {
            if ( index($0,"=") ) break;
            matchLine=$0
            gsub("{","",matchLine)
            gsub("}","",matchLine)
            if ( match(matchLine, /[a-zA-Z0-9]/) ) print matchLine
          }
        }
      }
    }
  }' $cfgFilename

}

GetMUvar()
{

  # Get configuration variable for muon

  local what="$1"
  local filename="${outputDir}/${what}.txt"

  if [ "$what" = "triggerList" ]; then
    ParseMUconfig "$what" > $filename
    if [ -s "$filename" ]; then
      output="\"${filename}\""
    else
      rm $filename
      output="0x0"
    fi
  else
    output="$(ParseMUconfig $what)"
    if [ "$output" = "" ]; then
      # Flag not specified: use default
      if [ "$what" = "isMC" ]; then
        output="0" && [[ "$dataType" = "sim" ]] && output="1"
      elif [ "$what" = "MUenablePhysSel" ]; then
        output="1" && [[ "$dataType" = "sim" ]] && output="0"
      elif [ "$what" = "MUcheckTrigScalers" ]; then
        output="1"
      fi
    fi
  fi
  echo "$output"
}

runLevelQA()
{
  #full path of QAresults.root is provided
  local qaFile=$1
  local isMC
  isMC=$(GetMUvar "isMC")
  local MUenablePhysSel
  MUenablePhysSel=$(GetMUvar "MUenablePhysSel")

  # This is used only to extract the muon information
  ln -s $ALICE_PHYSICS/PWGPP/MUON/lite/MakeTrend.C
  aliroot -b -l <<EOF
gSystem->AddIncludePath("-I$ALICE_ROOT/include -I$ALICE_PHYSICS/include");
.x MakeTrend.C("${qaFile}",${runNumber},$isMC,$MUenablePhysSel)
.q
EOF
  rm MakeTrend.C
}

periodLevelQA()
{
  #path of the merged period trending.root is provided
  # local trendingFile=$1

  local fileList="trendList.txt"
  local tmpFileList="tmp${fileList}"
#  local fileNames="QAresults.root QAresults_barrel.root QAresults_outer.root"
#  find -L . -name "QAresults*.root" > $fileList
  local fileNames="QAresults.root"
  find -L . -name "QAresults.root" > $fileList

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

  local excludeRuns
  excludeRuns=$(GetMUvar "excludeRuns")
  for irun in $excludeRuns; do
    sed -n -i'' '/'"$irun"'/!p' $fileList
  done

  #if run list is provided, filter the output limiting to this list
  local runList
  runList=$(GetMUvar "runList")
  if [ "$runList" != "" ]; then
    tmpFileList="tmp${fileList}"
    mv ${fileList} ${tmpFileList}
    for irun in $runList; do
      currFile=$(grep ${irun} $tmpFileList)
      if [ "${currFile}" = "" ]; then
        continue
      fi
      echo "${currFile}" >> ${fileList}
    done
    rm $tmpFileList
  fi

  #if trigger list is provided, filter the tracking output accordngly
  local triggerList
  triggerList=$(GetMUvar "triggerList")

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

  local MUenablePhysSel
  MUenablePhysSel=$(GetMUvar "MUenablePhysSel")

  ln -s $ALICE_PHYSICS/PWGPP/MUON/lite/PlotMuonQA.C
aliroot -b <<EOF
gSystem->AddIncludePath("-I$ALICE_ROOT/include -I$ALICE_PHYSICS/include");
.x PlotMuonQA.C+(".","${fileList}",${triggerList},${MUenablePhysSel},"muon_tracker","${mergedQAname}");
.q
EOF
  rm PlotMuonQA.C
  if [ -e "${triggerList}" ]; then
    rm "${triggerList}"
  fi

  # Then run trigger
  local MUcheckTrigScalers
  MUcheckTrigScalers=$(GetMUvar "MUcheckTrigScalers")

  ln -s $ALICE_PHYSICS/PWGPP/MUON/lite/trigEffQA.C
  aliroot -b <<EOF
gSystem->AddIncludePath("-I$ALICE_ROOT/include -I$ALICE_PHYSICS/include");
.x trigEffQA.C+("${fileList}","QA_muon_trigger.root","${ocdbStorage}",$MUcheckTrigScalers,"${mergedQAname}");
.q
EOF
  rm trigEffQA.C

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
