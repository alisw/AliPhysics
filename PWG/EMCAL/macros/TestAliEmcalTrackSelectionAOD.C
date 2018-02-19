int TestAliEmcalTrackSelectionAOD(){
  PWG::EMCAL::TestAliEmcalTrackSelectionAOD testrunner;
  testrunner.Init();
  if(testrunner.RunAllTests()) return 0;
  return 1; 
}