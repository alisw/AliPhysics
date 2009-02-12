void SETUP() {
  gSystem->Load("libPWG2flowCommon");
  
  // Set the include paths
  gROOT->ProcessLine(".include PWG2flowCommon/FLOW/AliFlowCommon");

  // Set our location, so that other packages can find us
  gSystem->Setenv("PWG2flowCommon_INCLUDE", "PWG2flowCommon/FLOW/AliFlowCommon");
}

