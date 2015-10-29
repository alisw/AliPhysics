void runReconstruction() {
  
  //AliLog::SetGlobalDebugLevel(1);
  /*
  gSystem->Load("libITSUpgradeBase.so");
  gSystem->Load("libITSUpgradeSim.so");
  gSystem->Load("libITSUpgradeRec.so");

  // Set ITS upgrade reconstructor
  
  gPluginMgr->AddHandler("AliReconstructor", "ITS",
                         "AliITSUReconstructor","ITS", "AliITSUReconstructor()");
  */
  TDatime dt;
  UInt_t seed = dt.Get();
  
  gRandom->SetSeed(seed);
  
  AliReconstruction *reco = new AliReconstruction("galice.root");
  
  // switch off cleanESD
  reco->SetCleanESD(kFALSE);

  reco->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  //reco->SetDefaultStorage("alien://folder=/alice/data/2011/OCDB");
  
  // GRP from local OCDB
  reco->SetSpecificStorage("GRP/GRP/Data",Form("local://%s",gSystem->pwd()));

  // ITS (2 objects)
  //reco->SetSpecificStorage("ITS/Align/Data",      "alien://folder=/alice/simulation/LS1_upgrade/Ideal");
  //reco->SetSpecificStorage("ITS/Calib/RecoParam", "alien://folder=/alice/simulation/LS1_upgrade/Ideal");

  //reco->SetSpecificStorage("ITS/Align/Data", Form("local://%s",gSystem->pwd()));
  //reco->SetSpecificStorage("ITS/Calib/RecoParam", Form("local://%s",gSystem->pwd()));
  
  // MUON Tracker
  
  reco->SetSpecificStorage("MUON/Align/Data",      "alien://folder=/alice/simulation/2008/v4-15-Release/Residual");
  reco->SetSpecificStorage("MUON/Calib/RecoParam", "alien://folder=/alice/cern.ch/user/a/auras/OCDB/");
  reco->SetSpecificStorage("MFT/Align/Data",	   "alien://folder=/alice/cern.ch/user/a/auras/OCDB/");
  reco->SetSpecificStorage("MFT/Calib/RecoParam",  "alien://folder=/alice/cern.ch/user/a/auras/OCDB/");
  /*
  // copied to local
  reco->SetSpecificStorage("MUON/Calib/RecoParam", "local:///users/calcul/vulpescu/alice/Work/MFT/CA-test/CA-ali/feature-itsmft/OCDB");
  reco->SetSpecificStorage("MFT/Align/Data",       "local:///users/calcul/vulpescu/alice/Work/MFT/CA-test/CA-ali/feature-itsmft/OCDB");
  reco->SetSpecificStorage("MFT/Calib/RecoParam",  "local:///users/calcul/vulpescu/alice/Work/MFT/CA-test/CA-ali/feature-itsmft/OCDB");
  */
  reco->SetOption("MUON MFT","SAVEDIGITS");
  reco->SetRunQA(":");
  //reco->SetQAWriteExpert(AliQAv1::kMUON);
  //reco->SetQARefDefaultStorage("local://$ALICE_ROOT/QAref");
  
  reco->SetWriteESDfriend(kFALSE);
  reco->SetStopOnError(kFALSE);
  TStopwatch timer;
  timer.Start();
  reco->Run();
  timer.Stop();
  timer.Print();
  
  delete reco;
  AliCodeTimer::Instance()->Print();

}
