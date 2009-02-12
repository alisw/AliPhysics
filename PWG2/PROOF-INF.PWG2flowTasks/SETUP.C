void SETUP() {
  gSystem->Load("libPWG2flowTasks");
  
  // Set the include paths
  gROOT->ProcessLine(".include PWG2flowTasks/FLOW/AliFlowTasks");

  // Set our location, so that other packages can find us
  gSystem->Setenv("PWG2flowTasks_INCLUDE", "PWG2flowTasks/FLOW/AliFlowTasks");
}

