void cloneAODTreeAndAddObject(const char *newFileName = "AliAOD_new.root", const char *orgFileName = "AliAOD.root") {
  // This little macro takes an already created AOD file and clones it.
  // After adding a new brach, the new TTree is written to a new file.
  //
  // To read back this data (including the newly updated branch), have
  // a look at $ALICE_ROOT/STEER/ReadAOD.C and add one line
  //       TClonesArray *esdTracks = 
  //            (TClonesArray*)ev->GetList()->FindObject("esdTracks");
  // after   
  //       ev->ReadFromTree(aodTree);
  //
  // esdTracks will be 'filled' correctly after each call of
  //       aodTree->GetEvent(nEv);
  //


  // open input file and get the TTree
  TFile orgFile(orgFileName, "READ");
  // get original TTree
  TTree *orgAodTree = (TTree*)orgFile.Get("aodTree");

  // open new output file
  TFile *newFile = new TFile(newFileName, "RECREATE");
  // clone old TTree
  TTree *newAodTree = orgAodTree->CloneTree();

  // do your gymnastics with the new TTree
  AliAODEvent *evNew = new AliAODEvent();
  evNew->ReadFromTree(newAodTree);
    
  // add new information to the list (we use AliESDtracks as an example)
  TClonesArray *tracks = new TClonesArray("AliESDtrack", 0);
  evNew->AddObject(tracks); // add new object to the list

  // define a unique name
  const char *name = "esdTracks";
  tracks->SetName(name); // set this name to be the name of the TClonesArray
  // create the new branch
  TBranch *newBranch = newAodTree->Branch(name, &tracks); // the branch gets the same name

  // loop over events
  Int_t nEvent = newAodTree->GetEntries();

  for(Int_t iEv = 0; iEv < nEvent; iEv++) {

    /*
    // read events (only necessary if you want to access the old data)
    newAodTree->GetEvent(iEv);
    */

    tracks->Delete(); // for each event delete old entries of new TClonesArray

    Int_t nTracks = gRandom->Rndm() * 50; // randomize size of TClonesArray (just for this example)
    tracks->Expand(nTracks); // expand container (just for speed)
    
    // fill TClonesArray
    TClonesArray &rTracks = *tracks;
    for (Int_t iTr = 0; iTr< nTracks; iTr++) {
      new(rTracks[iTr]) AliESDtrack();
    }
    
    // fill the new branch
    newBranch->Fill();
  }

  // delete old and write new UserInfo
  newAodTree->GetUserInfo()->Clear();
  newAodTree->GetUserInfo()->Add(evNew);

  // write new TTree to file
  newAodTree->Write();

  // close files
  newFile->Close();
  delete newFile;

  orgFile.Close();
}
