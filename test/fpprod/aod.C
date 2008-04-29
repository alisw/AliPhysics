void aod(){

  gSystem->Load("libdpmjet");
  gSystem->Load("libTDPMjet");
  gROOT->Macro("$ALICE_ROOT/STEER/CreateAODfromESD.C");
}
