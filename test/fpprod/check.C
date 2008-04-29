void check(){

  gSystem->Load("libdpmjet");
  gSystem->Load("libTDPMjet");
  gROOT->Macro("$ALICE_ROOT/STEER/CheckESD.C");
}
