void runReconstruction(const Char_t *recOptions) {
  
  TDatime dt;
  UInt_t seed = dt.Get();

  gRandom->SetSeed(seed);

  AliReconstruction *reco = new AliReconstruction("galice.root");

  // switch off cleanESD
  reco->SetCleanESD(kFALSE);

  // GRP from local OCDB
  reco->SetSpecificStorage("GRP/GRP/Data",Form("local://%s",gSystem->pwd()));
  
  // MUON Tracker
  reco->SetSpecificStorage("MUON/Align/Data",     "alien://folder=/alice/simulation/2008/v4-15-Release/Residual");
  reco->SetSpecificStorage("MUON/Calib/RecoParam","alien://folder=/alice/cern.ch/user/a/auras/OCDB/");

  reco->SetOption("MUON MFT",recOptions);
  reco->SetRunQA(":");
//   reco->SetQAWriteExpert(AliQAv1::kMUON);
//   reco->SetQARefDefaultStorage("local://$ALICE_ROOT/QAref");

  reco->SetWriteESDfriend(kFALSE);
  reco->SetStopOnError(kFALSE);

  TStopwatch timer;
  timer.Start();
  reco->Run();
  timer.Stop();
  timer.Print();

  delete reco;

}
