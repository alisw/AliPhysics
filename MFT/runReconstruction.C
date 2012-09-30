void runReconstruction(Int_t seed, const Char_t *recOptions) {
  
  gRandom->SetSeed(seed);

  AliReconstruction *reco = new AliReconstruction("galice.root");

  // switch off cleanESD
  reco->SetCleanESD(kFALSE);

  // GRP from local OCDB
  reco->SetSpecificStorage("GRP/GRP/Data",Form("local://%s",gSystem->pwd()));
  
  reco->SetDefaultStorage("alien://folder=/alice/data/2011/OCDB");

  reco->SetSpecificStorage("MUON/Align/Data",                      "alien://folder=/alice/simulation/2008/v4-15-Release/Residual");
  reco->SetSpecificStorage("MUON/Calib/Capacitances",              "alien://folder=/alice/simulation/2008/v4-15-Release/Residual");
  reco->SetSpecificStorage("MUON/Calib/Config",                    "alien://folder=/alice/simulation/2008/v4-15-Release/Residual");
  reco->SetSpecificStorage("MUON/Calib/Gain",                      "alien://folder=/alice/simulation/2008/v4-15-Release/Residual");
  reco->SetSpecificStorage("MUON/Calib/GlobalTriggerBoardMasks",   "alien://folder=/alice/simulation/2008/v4-15-Release/Residual");
  reco->SetSpecificStorage("MUON/Calib/GlobalTriggerCrateConfig",  "alien://folder=/alice/simulation/2008/v4-15-Release/Residual");
  reco->SetSpecificStorage("MUON/Calib/HV",                        "alien://folder=/alice/simulation/2008/v4-15-Release/Residual");
  reco->SetSpecificStorage("MUON/Calib/LocalTriggerBoardMasks",    "alien://folder=/alice/simulation/2008/v4-15-Release/Residual");
  reco->SetSpecificStorage("MUON/Calib/MappingData",               "alien://folder=/alice/simulation/2008/v4-15-Release/Residual");
  reco->SetSpecificStorage("MUON/Calib/MappingRunData",            "alien://folder=/alice/simulation/2008/v4-15-Release/Residual");
  reco->SetSpecificStorage("MUON/Calib/Neighbours",                "alien://folder=/alice/simulation/2008/v4-15-Release/Residual");
  reco->SetSpecificStorage("MUON/Calib/OccupancyMap",              "alien://folder=/alice/simulation/2008/v4-15-Release/Residual");
  reco->SetSpecificStorage("MUON/Calib/Pedestals",                 "alien://folder=/alice/simulation/2008/v4-15-Release/Residual");
  reco->SetSpecificStorage("MUON/Calib/RegionalTriggerBoardMasks", "alien://folder=/alice/simulation/2008/v4-15-Release/Residual");
  reco->SetSpecificStorage("MUON/Calib/RegionalTriggerConfig",     "alien://folder=/alice/simulation/2008/v4-15-Release/Residual");
  reco->SetSpecificStorage("MUON/Calib/RejectList",                "alien://folder=/alice/simulation/2008/v4-15-Release/Residual");
  reco->SetSpecificStorage("MUON/Calib/TriggerDCS",                "alien://folder=/alice/simulation/2008/v4-15-Release/Residual");
  reco->SetSpecificStorage("MUON/Calib/TriggerEfficiency",         "alien://folder=/alice/simulation/2008/v4-15-Release/Residual");
  reco->SetSpecificStorage("MUON/Calib/TriggerLut",                "alien://folder=/alice/simulation/2008/v4-15-Release/Residual");

  reco->SetSpecificStorage("MUON/Calib/RecoParam", "alien://folder=/alice/cern.ch/user/a/auras/OCDB/");
  reco->SetSpecificStorage("MFT/Align/Data",       "alien://folder=/alice/cern.ch/user/a/auras/OCDB/");
  reco->SetSpecificStorage("MFT/Calib/RecoParam",  "alien://folder=/alice/cern.ch/user/a/auras/OCDB/");

  reco->SetRunReconstruction("MUON MFT");
  reco->SetRunLocalReconstruction("MUON MFT");
  reco->SetOption("MUON MFT",recOptions);
  //  reco->SetRunQA("DetectorList:ActionList");
  //  reco->SetQAWriteExpert(AliQAv1::kMUON);

  reco->SetWriteESDfriend(kFALSE);
  reco->SetStopOnError(kFALSE);

  TStopwatch timer;
  timer.Start();
  reco->Run();
  timer.Stop();
  timer.Print();

  delete reco;

}
