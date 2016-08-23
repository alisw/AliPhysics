Int_t n = 100; // This makes n "visible" to the compiled macro
void runsim(Int_t nev = 1){
  gSystem->Load("liblhapdf");      // Parton density functions
  gSystem->Load("libEGPythia6");   // TGenerator interface
  gSystem->Load("libpythia6");     // Pythia
  gSystem->Load("libAliPythia6");  // ALICE specific implementations
  gSystem->Load("libgeant321");

  n = nev;
  gROOT->Macro("sim.C(n)");
}
