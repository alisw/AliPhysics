
#  Single track efficeincy QA > Standard Variables:
#  $dataType     e.g. data or sim
#  $year         e.g. 2011
#  $period       e.g. LHC13g
#  $runNumber    e.g. 169123
#  $pass         e.g. cpass1, pass1, passMC


runLevelQA()
{
  [[ ${dataType} == "data" ]] && return
  qaFile=$1
  cp $ALICE_PHYSICS/PWGPP/EvTrkSelection/macros/SingleTrackEffTrend.C .
  cp $ALICE_PHYSICS/PWGPP/EvTrkSelection/macros/CalcSingleTrackEffQA.C .
  aliroot -b -q -l "SingleTrackEffTrend.C(\"$qaFile\",$runNumber)"
}


periodLevelQA()
{
  [[ ${dataType} == "data" ]] && return
  trendingFile=$1
  cp $ALICE_PHYSICS/PWGPP/EvTrkSelection/macros/periodLevelQAEff.C .
  root -b -q -l "periodLevelQAEff.C(\"trending.root\")"
}
