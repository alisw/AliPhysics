#available variables:
#  $dataType     e.g. data or sim
#  $year         e.g. 2011
#  $period       e.g. LHC13g
#  $runNumber    e.g. 169123
#  $pass         e.g. cpass1, pass1, passMC

runLevelQA()
{
  qaFile=$1
  cp $ALICE_PHYSICS/PWGPP/HMPID/macros/makeHMPQA.C .
  aliroot -b -q -l "Analyze_QA_run(\"$dataType\", ${year}, \"$period\", \"$pass\",${runNumber})"
}

periodLevelQA()
{
  trendingFile=$1
  cp $ALICE_PHYSICS/PWGPP/HMPID/macros/makeHMPQA.C .
  aliroot -b -q -l "makeHMPQA(\"$dataType\", ${year}, \"$period\", \"$pass\")"
}
