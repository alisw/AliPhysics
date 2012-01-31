void SETUP() {

  TString dypath = gSystem->GetDynamicPath();
  dypath.Prepend(".:");
  gSystem->SetDynamicPath(dypath);
  gSystem->Load("libPWGflowBase");

  // Set the include paths
  gROOT->ProcessLine(".include PWGflowBase/FLOW/Base");

  // Set our location, so that other packages can find us
  gSystem->Setenv("PWGflowBase_INCLUDE", "PWGflowBase/FLOW/Base");
}

