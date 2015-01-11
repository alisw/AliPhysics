void check(){

  gSystem->Load("liblhapdf");      // Parton density functions
  gSystem->Load("libEGPythia6");   // TGenerator interface
  gSystem->Load("libpythia6");     // Pythia
  gSystem->Load("libAliPythia6");  // ALICE specific implementations
  gSystem->Load("libDPMJET");
  gSystem->Load("libTDPMjet");
  gROOT->Macro("$ALICE_ROOT/STEER/macros/CheckESD.C");
}
