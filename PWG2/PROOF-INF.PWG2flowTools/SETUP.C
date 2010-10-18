void SETUP() {
  gSystem->Load("libPWG2flowTools");
  
  // Set the include paths
  gROOT->ProcessLine(".include PWG2flowTools/FLOW/AliFlowTools/glauberMC");

  // Set our location, so that other packages can find us
  gSystem->Setenv("PWG2flowTools_INCLUDE", "PWG2flowTools/FLOW/AliFlowTools");
}

