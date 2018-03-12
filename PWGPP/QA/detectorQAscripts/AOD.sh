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
  mkdir -p trackQA
  rm -f trackQA/*.png
  mv *.png trackQA

  cp $ALICE_PHYSICS/PWGPP/macros/PlotAODvertexQA.C .
  root -b -q -l "PlotAODvertexQA.C(\"$qaFile\",\"QA\",${runNumber})"
  mkdir -p vertexQA
  rm -f vertexQA/*.png
  mv *.png vertexQA

  cp $ALICE_PHYSICS/PWGPP/macros/MakeTrendingPIDQA.C .
  root -b -q -l "MakeTrendingPIDQA.C(\"$qaFile\",\"PIDqa\",\"all\",${runNumber})"
  mv trending.root trendingPID.root
  mkdir -p pidQA
  rm -f pidQA/*.png
  mv *.png pidQA

  cp $ALICE_PHYSICS/PWGHF/vertexingHF/macros/readMCPerform.C .
  root -b -q -l "readMCPerform.C(\"$qaFile\",1,${runNumber})"
  mkdir -p hfQA
  rm -f hfQA/*.png
  mv *.png hfQA

  hadd -k trending.root trendingAODtracks.root trendingAODvertex.root trendingPID.root trendingHF.root

  rm -f trendingAODtracks.root
  rm -f trendingAODvertex.root
  rm -f trendingPID.root
  rm -f trendingHF.root

}

periodLevelQA()
{

}