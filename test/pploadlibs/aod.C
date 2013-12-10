void aod(){

  gROOT->Macro("loadlibs.C");

  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libCORRFW");
  gSystem->Load("libPWGHFbase");
  gSystem->Load("libPWGmuon");

  gROOT->Macro("${ALICE_ROOT}/STEER/CreateAODfromESD.C");
}
