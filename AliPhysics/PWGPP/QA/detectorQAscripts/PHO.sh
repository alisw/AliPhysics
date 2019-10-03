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

  cp $ALICE_PHYSICS/PWGGA/PHOSTasks/AutoTrendQA/MakeTrendingPHOSQA.C .
  aliroot -b -q -l "MakeTrendingPHOSQA.C(\"$qaFile\",${runNumber},kFALSE)" 
}

periodLevelQA()
{
  trendingFile=$1

  cp $ALICE_PHYSICS/PWGGA/PHOSTasks/AutoTrendQA/DrawTrendingPHOSQA.C .
  aliroot -b -q -l "DrawTrendingPHOSQA.C(\"trending.root\")"
}
