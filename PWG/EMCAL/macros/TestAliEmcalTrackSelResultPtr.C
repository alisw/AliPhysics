int TestAliEmcalTrackSelResultPtr() {
  PWG::EMCAL::TestAliEmcalTrackSelResultPtr testrunner;
  if(testrunner.RunAllTests()) return 0;
  return 1; 
}