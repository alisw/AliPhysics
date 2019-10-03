Bool_t CheckAOD(const char* aodFileName = "AliAOD.Muons.root")
{
  //nb (commented)
  //AliCDBManager *cdb = AliCDBManager::Instance();
  //cdb->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  
  // open the ESD file
  TFile* aodFile = TFile::Open(aodFileName);
  if (!aodFile || !aodFile->IsOpen()) {
    Error("CheckAOD", "opening AOD file %s failed", aodFileName);
    return kFALSE;
  }
  AliAODEvent *aod = new AliAODEvent();
  TTree* tree = (TTree*) aodFile->Get("aodTree");
  if (!tree) {
    Error("CheckAOD", "no AOD tree found");
    return kFALSE;
  }
  aod->ReadFromTree(tree);
  
  // loop over events
  for (Int_t iEvent = 0; iEvent < tree->GetEntries(); iEvent++) {
    
    // get the event data
    if (tree->GetEvent(iEvent) <= 0) {
      Error("CheckAOD", "no AOD object found for event %d", iEvent);
      return kFALSE;
    }
    
  }
  
  delete aod;
  aodFile->Close();
  delete aodFile;
  
  // result of check
  Info("CheckAOD", "check of AOD was successfull");
  return kTRUE;
  
}

