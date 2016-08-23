void runcheck(){
  gSystem->Load("libpythia6");     // Pythia
  gSystem->AddIncludePath("-I${ALICE_ROOT}/include");
  gROOT->Macro("${ALICE_ROOT}/STEER/macros/CheckESD.C++");
}
