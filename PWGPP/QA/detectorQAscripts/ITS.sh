#available variables:
#  $dataType     e.g. data or sim
#  $year         e.g. 2011
#  $period       e.g. LHC13g
#  $runNumber     e.g. 169123
#  $pass         e.g. cpass1,pass1,passMC
#  #ocdbStorage  e.g. "raw://", "local://./OCDB"

runLevelQA()
{
  qaFile=$1

  cp $ALICE_PHYSICS/PWGPP/ITS/QA/MakeTrendingITSQA.C .
  aliroot -b -q -l "MakeTrendingITSQA.C(\"$qaFile\",${runNumber},kFALSE,kFALSE,kFALSE,\"${ocdbStorage}\")"
}

periodLevelQA()
{
  trendingFile=$1

  cp $ALICE_PHYSICS/PWGPP/ITS/QA/DrawTrendingITSQA.C .
  aliroot -b -q -l "DrawTrendingITSQA.C(\"trending.root\")"
}
