#available variables:
#  $dataType     e.g. data or sim
#  $year         e.g. 2011
#  $period       e.g. LHC13g
#  $runNumber     e.g. 169123
#  $pass         e.g. cpass1,pass1,passMC
#  #ocdbStorage  e.g. "raw://", "local://./OCDB"

runLevelQA()
{
  qaFile="QAresults_AOD.root"

  cp $ALICE_PHYSICS/PWGPP/macros/PlotAODtrackQA.C .
  root -b -q -l "PlotAODtrackQA.C(\"$qaFile\",\"QA\",${runNumber})"

  cp $ALICE_PHYSICS/PWGPP/macros/PlotAODvertexQA.C .
  root -b -q -l "PlotAODvertexQA.C(\"$qaFile\",\"QA\")"

  cp $ALICE_PHYSICS/PWGPP/macros/MakeTrendingPIDQA.C .
  root -b -q -l "MakeTrendingPIDQA.C(\"$qaFile\",\"PIDqa\",\"all\",${runNumber})"

  cp $ALICE_PHYSICS/PWGHF/vertexingHF/macros/readMCPerform.C .
  root -b -q -l "readMCPerform.C(\"$qaFile\")"

}

periodLevelQA()
{

}