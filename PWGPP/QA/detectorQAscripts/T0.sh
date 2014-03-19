#available variables:
#  $dataType     e.g. data or sim
#  $year         e.g. 2011
#  $period       e.g. LHC13g
#  $runNumber     e.g. 169123
#  $pass         e.g. cpass1,pass1,passMC

runLevelQA()
{
  qaFile=$1

  cp $ALICE_ROOT/T0/MakeTrendT0.C .
  aliroot -b -q -l "MakeTrendT0.C(\"$qaFile\",${runNumber})" 
}

periodLevelQA()
{
  trendingFile=$1

  cp $ALICE_ROOT/T0/drawPerformanceT0QATrends.C .
  aliroot -b -q -l "drawPerformanceT0QATrends.C(\"$trendingFile\")"
}
