addObjectDuringAODCreation() {

  // add an object to an aod and write it

  TFile *aodFile = TFile::Open("addAOD.root", "RECREATE");

    // create an AliAOD object 
  AliAODEvent *aod = new AliAODEvent();
  aod->CreateStdContent();
  
  // add new information, we use AliESDtracks for now
  TClonesArray *tracks = new TClonesArray("AliESDtrack", 0);
  aod->AddObject(tracks);

  // go to the file
  aodFile->cd();
  
  // create the tree
  TTree *aodTree = new TTree("aodTree", "AliAOD tree");
  aodTree->Branch(aod->GetList());

  for (Int_t iEvent = 0; iEvent < 10; ++iEvent) {
    // add (part of) standard information
    AliAODHeader *header = aod->GetHeader();

    tracks->Delete(); // delete old objects
    tracks->Expand(iEvent+5/* just to make it a different number each time*/); // expand container (just for speed)
    
    // fill TClonesArray
    TClonesArray &rTracks = *tracks;
    for (Int_t i = 0; i< iEvent+5; i++) {
      new(rTracks[i]) AliESDtrack();
    }

    // fill the tree for this event
    aodTree->Fill();
  } // end of event loop

  aodTree->GetUserInfo()->Add(aod);

  // write the tree to the specified file
  aodFile = aodTree->GetCurrentFile();
  aodFile->cd();
  aodTree->Write();



}
