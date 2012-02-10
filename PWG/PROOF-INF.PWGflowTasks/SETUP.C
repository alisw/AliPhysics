void SETUP() {

  TString dypath = gSystem->GetDynamicPath();
  dypath.Prepend(".:");
  gSystem->SetDynamicPath(dypath);
  gSystem->Load("libPWGflowBase");
  gSystem->Load("libPWGflowTasks");
  
  // Set the include paths
  gROOT->ProcessLine(".include PWGflowTasks/FLOW/Tasks $(ALICE_ROOT)/TOF");

  // Set our location, so that other packages can find us
  gSystem->Setenv("PWGflowTasks_INCLUDE", "PWGflowTasks/FLOW/Tasks");
}

