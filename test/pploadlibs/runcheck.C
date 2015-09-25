void runcheck(){
  gSystem->Load("libpythia6");     // Pythia
  gROOT->Macro("${ALICE_ROOT}/STEER/macros/CheckESD.C++");
}
