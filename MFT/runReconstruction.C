void runReconstruction(Int_t seed, const Char_t *recOptions) {
  
  gRandom->SetSeed(seed);

  AliReconstruction *reco = new AliReconstruction("galice.root");

  // switch off cleanESD
  reco->SetCleanESD(kFALSE);

  // GRP from local OCDB
  reco->SetSpecificStorage("GRP/GRP/Data",Form("local://%s",gSystem->pwd()));
  
  // MUON Tracker -> local:///$OCDB should reflect the content of alien://folder=/alice
  reco->SetDefaultStorage("alien://folder=/alice/simulation/2008/v4-15-Release/Ideal");
  reco->SetSpecificStorage("MUON/Align/Data",     "alien://folder=/alice/simulation/2008/v4-15-Release/Residual");
  reco->SetSpecificStorage("MFT/Align/Data",      "alien://folder=/alice/cern.ch/user/a/auras/OCDB/");
  reco->SetSpecificStorage("MFT/Calib/RecoParam", "alien://folder=/alice/cern.ch/user/a/auras/OCDB/");

  AliMUONRecoParam *param = AliMUONRecoParam::GetLowFluxParam();
  param->SetPadGoodnessMask(0x8080);      
  for (Int_t iCh=0; iCh<10; iCh++) {
    param->SetDefaultNonBendingReso(iCh,0.2);
    param->SetDefaultBendingReso(iCh,0.2);
  }
  param->SetSigmaCutForTracking(5.);
  param->ImproveTracks(kTRUE, 4.);
  param->SetStripCutForTrigger(1.5);
  param->SetSigmaCutForTrigger(4.);
  param->Print("FULL");
  reco->SetRecoParam("MUON", param);

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
