runLevelEventStatQA()
{
  eventStatFile=$1

  cp $ALICE_PHYSICS/PWGPP/EVS/runLevelEventStatQA.C .
  cp $ALICE_PHYSICS/PWGPP/EVS/triggerInfo.C .
  aliroot -b -q -l "runLevelEventStatQA.C(\"$eventStatFile\",${runNumber},\"${ocdbStorage}\")" 
}

periodLevelQA()
{
  trendingFile=$1

  cp $ALICE_PHYSICS/PWGPP/EVS/periodLevelQA.C .
  aliroot -b -q -l "periodLevelQA.C(\"trending.root\")"
}
