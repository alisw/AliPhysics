void aod(){
  if (!strcmp(gSystem->GetBuildArch(),"win32gcc")) {
    gSystem->Load("libProof");
    gSystem->Load("libGui");
    gROOT->Macro("loadlibs.C");
    new AliRun("gAlice","The ALICE Off-line Simulation Framework");
  }

  gSystem->Load("libdpmjet");
  gSystem->Load("libTDPMjet");
  gROOT->Macro("$ALICE_ROOT/STEER/CreateAODfromESD.C");
}
