int TestAliEmcalAODHybridTrackCuts() {
  PWG::EMCAL::TestAliEmcalAODHybridTrackCuts testrunner;
  testrunner.Init();
  if(testrunner.RunAllTests()) return 0;
  return 1; 
}