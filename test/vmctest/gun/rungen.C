Int_t n = 100; // This makes n "visible" to the compiled macro
void rungen(Int_t nev, const char* configMacro) {
  // Load libraries
  gSystem->Load("liblhapdf");      // Parton density functions
  gSystem->Load("libEGPythia6");   // TGenerator interface
  gSystem->Load("libpythia6");     // Pythia
  gSystem->Load("libAliPythia6");  // ALICE specific implementations


  // Load configuration macros

  gROOT->LoadMacro("$ALICE_ROOT/test/vmctest/gun/commonConfig.C");
  gROOT->LoadMacro("$ALICE_ROOT/test/vmctest/gun/genGunConfig.C");
  gROOT->LoadMacro(configMacro);  

  n = nev; // Use the requested number of events
  // Execute gen.C
  gROOT->Macro("gen.C(n)");
}
