runLevelEventStatQA()
{
  eventStatFile=$1

  cp $ALICE_ROOT/PWGPP/EVS/runLevelEventStatQA.C .
  cp $ALICE_ROOT/PWGPP/EVS/triggerInfo.C .
  aliroot -b -q -l "runLevelEventStatQA.C(\"$eventStatFile\",${runNumber})" 
}

periodLevelQA()
{
  trendingFile=$1

  cp $ALICE_ROOT/PWGPP/EVS/periodLevelQA.C .
  aliroot -b -q -l "periodLevelQA.C(\"trending.root\")"
}
