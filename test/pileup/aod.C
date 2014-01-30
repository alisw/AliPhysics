void aod(){
  gSystem->Load("liblhapdf");      // Parton density functions
  gSystem->Load("libEGPythia6");   // TGenerator interface
  gSystem->Load("libpythia6");     // Pythia
  gSystem->Load("libAliPythia6");  // ALICE specific implementations
  gSystem->Load("libdpmjet");
  gSystem->Load("libTDPMjet");
 
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libCORRFW");
  gSystem->Load("libPWGHFbase");
  gSystem->Load("libPWGmuon");
  gSystem->Load("libESDfilter");
  gSystem->Load("libTENDER");
  gSystem->Load("libPWGPP");

  gROOT->Macro("${ALICE_ROOT}/STEER/CreateAODfromESD.C(\"AliESDs.root\",\"AliAODs.root\",\"local://$ALICE_ROOT/OCDB\",\"local://.\")");
}
