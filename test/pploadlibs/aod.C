void aod(){

  gROOT->Macro("loadlibs.C");

  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libCORRFW");
  gSystem->Load("libPWGHFbase");
  gSystem->Load("libPWGmuon");
  gSystem->Load("libESDfilter");
  gSystem->Load("libTender");
  gSystem->Load("libPWGPP");
  
  gROOT->Macro("${ALICE_ROOT}/STEER/macros/CreateAODfromESD.C(\"AliESDs.root\",\"AliAOD.root\",\"local://$ALICE_ROOT/OCDB\",\"local://.\")");
}
