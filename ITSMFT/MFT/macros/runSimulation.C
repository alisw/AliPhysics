void runSimulation(Int_t nevents=1,
                   Int_t runNumber=169099) {
  
  // AliLog::SetGlobalDebugLevel(1);
  
  AliSimulation *simulator = new AliSimulation();
  gSystem->Load("libFITbase.so");
  gSystem->Load("libFITsim.so");
  
  TDatime dt;
  UInt_t seed = dt.Get();
  
  simulator->SetSeed(seed);
  simulator->SetRunNumber(runNumber);
  simulator->SetMakeDigits("MUON MFT");
  simulator->SetMakeSDigits("TRD TOF PHOS HMPID EMCAL MUON ZDC PMD T0 VZERO FMD MFT");
  simulator->SetMakeDigitsFromHits("ITS TPC");
  simulator->SetRunQA(":");
  simulator->SetRunHLT("");
  
  gRandom->SetSeed(seed);
  simulator->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  //simulator->SetDefaultStorage("alien://folder=/alice/data/2011/OCDB");
    
  // MUON Tracker
  
  simulator->SetSpecificStorage("MUON/Align/Data",      "alien://folder=/alice/simulation/2008/v4-15-Release/Ideal");
  simulator->SetSpecificStorage("MUON/Calib/RecoParam", "alien://folder=/alice/cern.ch/user/a/auras/OCDB/");
  simulator->SetSpecificStorage("MFT/Align/Data",       "alien://folder=/alice/cern.ch/user/a/auras/OCDB/");
  simulator->SetSpecificStorage("MFT/Calib/RecoParam",  "alien://folder=/alice/cern.ch/user/a/auras/OCDB/");
  /*
  // copied to local
  simulator->SetSpecificStorage("MUON/Calib/RecoParam", "local:///users/calcul/vulpescu/alice/Work/MFT/CA-test/CA-ali/feature-itsmft/OCDB");
  simulator->SetSpecificStorage("MFT/Align/Data",       "local:///users/calcul/vulpescu/alice/Work/MFT/CA-test/CA-ali/feature-itsmft/OCDB");
  simulator->SetSpecificStorage("MFT/Calib/RecoParam",  "local:///users/calcul/vulpescu/alice/Work/MFT/CA-test/CA-ali/feature-itsmft/OCDB");
  */
  simulator->SetSpecificStorage("GRP/GRP/Data",
                                Form("local://%s",gSystem->pwd()));
  /*
  simulator->SetSpecificStorage("ITS/Align/Data",
                                Form("local://%s",gSystem->pwd()));
  simulator->SetSpecificStorage("ITS/Calib/SimuParam",
                                Form("local://%s",gSystem->pwd()));
  */
  //simulator->UseMagFieldFromGRP();
  
  // The rest
  TStopwatch timer;
  timer.Start();
  simulator->Run(nevents);
  timer.Stop();
  timer.Print();
    AliCodeTimer::Instance()->Print();
}
