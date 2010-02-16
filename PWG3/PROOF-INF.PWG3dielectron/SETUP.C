{
  gSystem->Load("libDielectron.so");
  // Set the include paths
  gROOT->ProcessLine(".include Dielectron");
  
  // Set our location, so that other packages can find us
  gSystem->Setenv("Dielectron_INCLUDE", "Dielectron");
  
}
