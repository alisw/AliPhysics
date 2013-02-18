Bool_t CheckESD(const char* esdFileName = "AliESDs.root")
{
//nb (commented)
    //AliCDBManager *cdb = AliCDBManager::Instance();
    //cdb->SetDefaultStorage("local://$ALICE_ROOT/OCDB");

  // open the ESD file
  TFile* esdFile = TFile::Open(esdFileName);
  if (!esdFile || !esdFile->IsOpen()) {
    Error("CheckESD", "opening ESD file %s failed", esdFileName);
    return kFALSE;
  }
  AliESDEvent * esd = new AliESDEvent;
  TTree* tree = (TTree*) esdFile->Get("esdTree");
  if (!tree) {
    Error("CheckESD", "no ESD tree found");
    return kFALSE;
  }
  esd->ReadFromTree(tree);

  // loop over events
  for (Int_t iEvent = 0; iEvent < tree->GetEntries(); iEvent++) {

    // get the event summary data
    tree->GetEvent(iEvent);
    if (!esd) {
      Error("CheckESD", "no ESD object found for event %d", iEvent);
      return kFALSE;
    }

    Int_t nTracks = esd->GetNumberOfMuonTracks();
    for(Int_t iTracks = 0; iTracks < nTracks; iTracks++) {
      AliESDMuonTrack* muonTrack = esd->GetMuonTrack(iTracks);
      if (muonTrack->ContainTrackerData()) {
      }
    }

  }

  delete esd;
  esdFile->Close();
  delete esdFile;

  // result of check
  Info("CheckESD", "check of ESD was successfull");
  return kTRUE;

}

