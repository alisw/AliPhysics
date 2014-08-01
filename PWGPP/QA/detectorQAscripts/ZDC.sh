#available variables:
#  $dataType     e.g. data or sim
#  $year         e.g. 2011
#  $period       e.g. LHC13g
#  $runNumber    e.g. 169123
#  $pass         e.g. cpass1,pass1,passMC

runLevelQA()
{
  qaFile=$1

  cp $ALICE_ROOT/PWGPP/ZDC/trending/MakeTrendingZDCQA.C .
  aliroot -b -q -l "MakeTrendingZDCQA.C(\"$qaFile\",${runNumber})" 
}

periodLevelQA()
{
  trendingFile=$1

  cp $ALICE_ROOT/PWGPP/ZDC/trending/DrawTrendingZDCQA.C .
  aliroot -b -q -l "DrawTrendingZDCQA.C(\"trending.root\")"
}
