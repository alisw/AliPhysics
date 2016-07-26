runLevelEventStatQA()
{
  eventStatFile=$1

  cp $ALICE_PHYSICS/PWGPP/EVS/runLevelEventStatQA.C .
  cp $ALICE_PHYSICS/PWGPP/EVS/triggerInfo.C .
  aliroot -b -q -l "runLevelEventStatQA.C(\"$eventStatFile\",${runNumber},\"${ocdbStorage}\")" > runLevel.log
}

periodLevelQA()
{
  trendingFile=$1

  cp $ALICE_PHYSICS/PWGPP/EVS/periodLevelQA.C .
  aliroot -b -q -l "periodLevelQA.C(\"trending.root\")" > periodLevel.log
}
