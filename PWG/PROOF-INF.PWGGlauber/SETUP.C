void SETUP() {


  TString dypath = gSystem->GetDynamicPath();
  dypath.Prepend(".:");
  gSystem->SetDynamicPath(dypath);
  gSystem->Load("libPWGGlauber");
  
  // Set the include paths
  gROOT->ProcessLine(".include PWGGlauber/Glauber");

  // Set our location, so that other packages can find us
  gSystem->Setenv("PWGGlauber_INCLUDE", "PWGGlauber/Glauber");
}

