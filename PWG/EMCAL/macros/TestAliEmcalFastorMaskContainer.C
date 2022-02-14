int TestAliEmcalFastorMaskContainer() {
  PWG::EMCAL::TestAliEmcalFastorMaskContainer testrunner;
  if(testrunner.RunAllTests()) return 0;
  return 1; 
}