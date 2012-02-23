void SETUP() {
  // Set the library paths
  TString dypath = gSystem->GetDynamicPath();
  dypath.Prepend("./PWGflowTasks/:");
  gSystem->SetDynamicPath(dypath);
  gSystem->Load("libPWGflowBase");
  gSystem->Load("libPWGflowTasks");
  
  // Set the include paths
  TString incpath = gSystem->GetIncludePath();
  incpath.Prepend("-I./PWGflowTasks/FLOW/Tasks -I$ALICE_ROOT/TOF ");
  gSystem->SetIncludePath(incpath);

  // Set our location, so that other packages can find us
  gSystem->Setenv("PWGflowTasks_INCLUDE", "PWGflowTasks/FLOW/Tasks");
}

