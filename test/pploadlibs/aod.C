void aod(){

  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libCORRFW");
  gSystem->Load("libPWGHFbase");
  gSystem->Load("libPWGmuon");
  gSystem->Load("libESDfilter");
  gSystem->Load("libTender");
  gSystem->Load("libPWGPP");
  
  gROOT->Macro("${ALICE_ROOT}/STEER/macros/CreateAODfromESD.C(\"AliESDs.root\",\"AliAOD.root\",\"local://$ALIROOT_OCDB_ROOT/OCDB\",\"local://.\")");
}
