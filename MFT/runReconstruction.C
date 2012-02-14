void runReconstruction(Int_t seed, const Char_t *recOptions) {
  
  gRandom->SetSeed(seed);

  AliReconstruction *reco = new AliReconstruction("galice.root");

  // switch off cleanESD
  reco->SetCleanESD(kFALSE);

  // GRP from local OCDB
  reco->SetSpecificStorage("GRP/GRP/Data",Form("local://%s",gSystem->pwd()));
  
  // MUON Tracker
  // reco->SetSpecificStorage("MUON/Align/Data","local:///$OCDB/simulation/2008/v4-15-Release/Residual");

  reco->SetRunReconstruction("MUON MFT");
  reco->SetRunLocalReconstruction("MUON MFT");
  reco->SetOption("MUON MFT",recOptions);
  reco->SetRunQA("MUON:ALL");
  reco->SetQAWriteExpert(AliQAv1::kMUON);
  reco->SetQARefDefaultStorage("local://$ALICE_ROOT/QAref");

  reco->SetWriteESDfriend(kFALSE);
  reco->SetStopOnError(kFALSE);

  TStopwatch timer;
  timer.Start();
  reco->Run();
  timer.Stop();
  timer.Print();

  delete reco;

}
