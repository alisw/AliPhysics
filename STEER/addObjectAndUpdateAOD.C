void addObjectAndUpdateAOD(const char *fileName = "AliAOD.root") {
  // This little macro takes an already created AOD file and adds
  // a new branch to it.
  // You need write access to the AOD file, which is generally not 
  // the case for AOD files available on the GRID/on CAF.
  // Be aware that if something breaks (lost connection, node goes 
  // down, ...) the content of the file will be essentially lost, 
  // meaning: 
  // WITH THIS UPDATE PROCEDURE YOU HAVE THE POWER TO DESTROY DATA!
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
  TFile inFile(fileName, "UPDATE");

  TTree *aodTree = (TTree*)inFile.Get("aodTree");
  AliAODEvent *ev = new AliAODEvent();
  ev->ReadFromTree(aodTree);

  // add new information to the list (we use AliESDtracks as an example)
  TClonesArray *tracks = new TClonesArray("AliESDtrack", 0);
  ev->AddObject(tracks); // add new object to the list

  // define a unique name
  const char *name = "esdTracks";
  tracks->SetName(name); // set this name to be the name of the TClonesArray
  // create the new branch
  TBranch *newBranch = aodTree->Branch(name, &tracks); // the branch gets the same name

  // loop over events
  Int_t nEvent = aodTree->GetEntries();

  for(Int_t iEv = 0; iEv < nEvent; iEv++) {

    /*
    // read events (only necessary if you want to access the old data)
    aodTree->GetEvent(nEv);
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

  // (over)write the tree (to the same file!)

  aodTree->Write("", TObject::kOverwrite);
  inFile.Close();

  return;
}
