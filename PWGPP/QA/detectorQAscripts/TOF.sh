#available variables:
#  $dataType     e.g. data or sim
#  $year         e.g. 2011
#  $period       e.g. LHC13g
#  $runNumber     e.g. 169123
#  $pass         e.g. cpass1,pass1,passMC

runLevelQA()
{
  qaFile=$1

  #cp $ALICE_ROOT/TOF/macrosQA/MakeTrendingTOFQA.C .
  #aliroot -b -q -l "MakeTrendingTOFQA.C(\"$qaFile\",${runNumber})" 
}

periodLevelQA()
{
  trendingFile=$1

  #cp $ALICE_ROOT/TOF/macrosQA/DrawTrending.C .
  #rm trending.root
  #ls 000*/trending.root > trending.list
  #aliroot -b -q -l "DrawTrending.C(\"trending.list\")"
}
