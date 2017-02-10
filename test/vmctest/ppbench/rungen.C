Int_t n = 100; // This makes n "visible" to the compiled macro
void rungen(Int_t nev) {
  // Load libraries
  gSystem->Load("liblhapdf");      // Parton density functions
  gSystem->Load("libEGPythia6");   // TGenerator interface
  gSystem->Load("libpythia6");     // Pythia
  gSystem->Load("libAliPythia6");  // ALICE specific implementations
  gSystem->Load("libHIJING");
  gSystem->Load("libTHijing");

  // Load configuration macros
  gROOT->LoadMacro("$ALICE_ROOT/test/vmctest/ppbench/genPPbenchConfig.C");


  n = nev; // Use the requested number of events
  // Execute gen.C
  gROOT->Macro("gen.C(n)");
}
