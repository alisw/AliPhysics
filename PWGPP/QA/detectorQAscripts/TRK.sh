#available variables:
#  $dataType     e.g. data or sim
#  $year         e.g. 2011
#  $period       e.g. LHC13g
#  $runNumer     e.g. 169123
#  $pass         e.g. cpass1,pass1,passMC
#  #ocdbStorage  e.g. "raw://", "local://./OCDB"

runLevelHighPtTreeQA()
{
  highPtTreeFile=${1}
    
  #makePlots needs to compile it locally for now
  cp $ALICE_PHYSICS/PWGPP/QA/Tracking/ExpertQA/AliHighPtTreeAnalysis.h .
  cp $ALICE_PHYSICS/PWGPP/QA/Tracking/ExpertQA/AliHighPtTreeAnalysis.cxx .

  aliroot -l -b -q "$ALICE_PHYSICS/PWGPP/QA/Tracking/ExpertQA/makePlots.C(\"$highPtTreeFile\")"
}

periodLevelQA()
{
  trendingFile=$1

  if ls */genericHistos_Bneg.root&>/dev/null; then
    hadd mergedGenericHistos_Bneg.root */genericHistos_Bneg.root
    mkdir outputBneg
    oldDir=${PWD}
    cd outputBneg
    mv ${oldDir}/mergedGenericHistos_Bneg.root .

    #makePlots needs to compile it locally for now
    cp $ALICE_PHYSICS/PWGPP/QA/Tracking/ExpertQA/AliHighPtTreeAnalysis.h .
    cp $ALICE_PHYSICS/PWGPP/QA/Tracking/ExpertQA/AliHighPtTreeAnalysis.cxx .

    aliroot -l -b -q "$ALICE_PHYSICS/PWGPP/QA/Tracking/ExpertQA/makePeriodPlots.C(\"mergedGenericHistos_Bneg.root\",\"${period}/${pass}\")"
    cd ${oldDir}
  fi
  if ls */genericHistos_Bpos.root&>/dev/null; then
    hadd mergedGenericHistos_Bpos.root */genericHistos_Bpos.root
    mkdir outputBpos
    oldDir=${PWD}
    cd outputBpos

    #makePlots needs to compile it locally for now
    cp $ALICE_PHYSICS/PWGPP/QA/Tracking/ExpertQA/AliHighPtTreeAnalysis.h .
    cp $ALICE_PHYSICS/PWGPP/QA/Tracking/ExpertQA/AliHighPtTreeAnalysis.cxx .

    mv ${oldDir}/mergedGenericHistos_Bpos.root .
    aliroot -l -b -q "$ALICE_PHYSICS/PWGPP/QA/Tracking/ExpertQA/makePeriodPlots.C(\"mergedGenericHistos_Bpos.root\",\"${period}/${pass}\")"
    cd ${oldDir}
  fi

  aliroot -l -b -q "$ALICE_PHYSICS/PWGPP/QA/Tracking/ExpertQA/makeTrendingPlots.C(\"${trendingFile}\")"
}
