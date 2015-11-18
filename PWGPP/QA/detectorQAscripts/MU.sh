#available variables:
#  $dataType     e.g. data or sim
#  $year         e.g. 2011
#  $period       e.g. LHC13g
#  $runNumber     e.g. 169123
#  $pass         e.g. cpass1,pass1,passMC
#  #ocdbStorage  e.g. "raw://", "local://./OCDB"

WriteMuonConfig()
{
  # FIXME
  # This function allows to create special configuration files
  # steer the MU.sh in special cases
  # (e.g. problems with CDB, Physics Selection, etc.)
  # When the macro will run at CERN it will be possible to
  # put custom cfg files when needed
  # and this function can be commented
  local cfgFilename="$1"
  if [ -e $cfgFilename ]; then

    # Change previous configuration if needed
    if [[ "$period" =~ ^LHC15[a-l]$ ]]; then
      local hasFeature=$(grep -c "MUenablePhysSel=0" $cfgFilename)
      if [ $hasFeature -gt 0 ]; then
        local tmpCfgFile=${cfgFilename/".txt"/"_tmp.txt"}
        sed 's/MUenablePhysSel=0/MUenablePhysSel=1/;' $cfgFilename > $tmpCfgFile
        mv $tmpCfgFile $cfgFilename
      fi
    fi
    return
  fi
  local cfgDir=$(dirname "$cfgFilename")
  if [ ! -e $cfgDir ]; then
    mkdir "$cfgDir"
  fi

  if [ "$period" = "LHC15e" ]; then
    echo "MUcheckTrigScalers=0" >> $cfgFilename
  fi
}

GetCfgVarFromFile()
{
  # Read special configuration variables for muon
  # Usually one does not need special configuration
  # However, it might be needed in case of issues
  # (arising e.g. when there are problems with CDB, Physics Selection, etc.)
  local what="$1"
  local filename="$2"
  local fullFlag=""
  if [ -e $filename ]; then
    fullFlag=$(grep "$what" $filename | xargs)
    fullFlag=${fullFlag/"$what"/""}
    fullFlag=${fullFlag/"="/""}
    fullFlag=${fullFlag/" "/""}
  fi
  echo $fullFlag
}

GetMUvar()
{
  # Get configuration variable for muon.

  # Assume that "outputDir" is known from the steering runQA.sh
  local cfgFileDir="${outputDir}/configFiles"
  local cfgFilePrefix="${dataType}_${period}_${pass}"
  local cfgFilename="${cfgFileDir}/${cfgFilePrefix}_cfgMU.txt"

  # FIXME: comment the following line when we will run at CERN
#  WriteMuonConfig "$cfgFilename"

  local what="$1"
  local output=""
  local filename=""

  if [ "$what" = "triggerList" ]; then
    filename="${cfgFileDir}/${cfgFilePrefix}_trigList.txt"
    if [ -e "$filename" ]; then
      output="\"${filename}\""
    else
      output="0x0"
    fi
  elif [ "$what" = "runList" ]; then
    filename="${cfgFileDir}/${cfgFilePrefix}_runList.txt"
    if [ -e "$filename" ]; then
      output="$(cat ${filename} | xargs)"
    else
      output=""
    fi
  elif [ "$what" = "isMC" ]; then
    if [ "$dataType" = "sim" ]; then
      output="1"
    else
      output="0"
    fi
  elif [ "$what" = "MUenablePhysSel" ]; then
    output=$(GetCfgVarFromFile "MUenablePhysSel" "$cfgFilename")
    if [ "$output" = "" ]; then
      if [ "$dataType" = "sim" ]; then
        output="0"
      else
        output="1"
      fi
    fi
  elif [ "$what" = "MUcheckTrigScalers" ]; then
    output=$(GetCfgVarFromFile "MUcheckTrigScalers" "$cfgFilename")
    if [ "$output" = "" ]; then
      output="1"
    fi
  fi
  echo "$output"
}


runLevelQA()
{
  #full path of QAresults.root is provided
  local qaFile=$1
  local isMC=$(GetMUvar "isMC")
  local MUenablePhysSel=$(GetMUvar "MUenablePhysSel")

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
  local trendingFile=$1

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

  #if run list is provided, filter the output limiting to this list
  local runList=$(GetMUvar "runList")
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
  local triggerList=$(GetMUvar "triggerList")

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

  local MUenablePhysSel=$(GetMUvar "MUenablePhysSel")

  ln -s $ALICE_PHYSICS/PWGPP/MUON/lite/PlotMuonQA.C
aliroot -b <<EOF
gSystem->AddIncludePath("-I$ALICE_ROOT/include -I$ALICE_PHYSICS/include");
.x PlotMuonQA.C+(".","${fileList}",${triggerList},${MUenablePhysSel},"muon_tracker","${mergedQAname}");
.q
EOF
  rm PlotMuonQA.C

  # Then run trigger
  local MUcheckTrigScalers=$(GetMUvar "MUcheckTrigScalers")

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
