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

  cp $ALICE_PHYSICS/PWGPP/TOF/trending/MakeTrendingTOFQAv2.C .
  aliroot -b -q -l "MakeTrendingTOFQAv2.C(\"$qaFile\",${runNumber}, 0, 1, -2., 2., \"${ocdbStorage}\", 0, 1, 1, 1)" 
}

periodLevelQA()
{
  trendingFile=$1

  cp $ALICE_PHYSICS/PWGPP/TOF/trending/DrawTrendingTOFQA.C .
  aliroot -b -q -l "DrawTrendingTOFQA.C(\"trending.root\")"
}
