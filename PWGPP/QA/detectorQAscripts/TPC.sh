#available variables:
#  $dataType     e.g. data or sim
#  $year         e.g. 2011
#  $period       e.g. LHC13g
#  $runNumer     e.g. 169123
#  $pass         e.g. cpass1,pass1,passMC

runLevelQA()
{
  qaFile=$1

  cp $ALICE_ROOT/PWGPP/TPC/macros/MakeTrend.C .
  aliroot -b -q -l "MakeTrend.C(\"$qaFile\",$runNumber)" 

  cp $ALICE_ROOT/PWGPP/TPC/macros/drawPerformanceTPCQAMatch.C .
  aliroot -b -q -l "drawPerformanceTPCQAMatch.C(\"$qaFile\")"
}

periodLevelQA()
{
  trendingFile=$1

  cp $ALICE_ROOT/PWGPP/TPC/macros/drawPerformanceTPCQAMatchTrends.C .
  cp $ALICE_ROOT/PWGPP/TPC/macros/qaConfig.C .
  aliroot -b -q -l "drawPerformanceTPCQAMatchTrends.C(\"trending.root\",\"PbPb\")"
}
