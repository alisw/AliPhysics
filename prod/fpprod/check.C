void check(){
  if (!strcmp(gSystem->GetBuildArch(),"win32gcc")) {
    gSystem->Load("libProof");
    gSystem->Load("libGui");
    gROOT->Macro("loadlibs.C");
    new AliRun("gAlice","The ALICE Off-line Simulation Framework");
  }

  gSystem->Load("liblhapdf");      // Parton density functions
  gSystem->Load("libEGPythia6");   // TGenerator interface
  gSystem->Load("libpythia6");     // Pythia
  gSystem->Load("libAliPythia6");  // ALICE specific implementations
  gSystem->Load("libdpmjet");
  gSystem->Load("libTDPMjet");
  gROOT->Macro("$ALICE_ROOT/STEER/CheckESD.C");
}
