void SETUP() {
  // Set the library paths
  TString dypath = gSystem->GetDynamicPath();
  dypath.Prepend("./PWGflowBase/:");
  gSystem->SetDynamicPath(dypath);
  gSystem->Load("libPWGflowBase");

  // Set the include paths
  TString incpath = gSystem->GetIncludePath();
  incpath.Prepend("-I./PWGflowBase/FLOW/Base ");
  gSystem->SetIncludePath(incpath);

  // Set our location, so that other packages can find us
  gSystem->Setenv("PWGflowBase_INCLUDE", "PWGflowBase/FLOW/Base");
}

