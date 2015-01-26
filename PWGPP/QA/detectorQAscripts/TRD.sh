#available variables:
#  $dataType     e.g. data or sim
#  $year         e.g. 2011
#  $period       e.g. LHC13g
#  $runNumer     e.g. 169123
#  $pass         e.g. cpass1,pass1,passMC
#  #ocdbStorage  e.g. "raw://", "local://./OCDB"

runLevelQA()
{
  qaFile=$1
  cp $ALICE_PHYSICS/PWGPP/TRD/macros/ProcessTRDRunQA.C .
  aliroot -b -q -l "ProcessTRDRunQA.C(\"${qaFile}\",${runNumber},\"${dataType}\",${year},\"${period}\",\"${pass}\",\"${ocdbStorage}\")"
}

periodLevelQA()
{
  trendingFile=$1

  cp $ALICE_PHYSICS/PWGPP/TRD/macros/DrawTrendingTRDQA.C .
  aliroot -b -q -l "DrawTrendingTRDQA.C(\"${trendingFile}\")"
}
