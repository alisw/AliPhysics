void ReadAOD(const char *fileName = "AliAOD.root") {

  // open input file and get the TTree
  TFile inFile(fileName, "READ");
  TTree *aodTree = (TTree*)inFile.Get("AOD");

  AliAODEvent *aod = aodTree->GetUserInfo()->FindObject("AliAODEvent");
  TIter next(aod->GetList());
  TNamed *el;
  while(el=(TNamed*)next()) 
    aodTree->SetBranchAddress(el->GetName(),aod->GetList()->GetObjectRef(el));

  // loop over events
  Int_t nEvents = aodTree->GetEntries();
  for (Int_t nEv = 0; nEv < nEvents; nEv++) {
    cout << "Event: " << nEv+1 << "/" << nEvents << endl;

    // read events
    aodTree->GetEvent(nEv);
    
    // set pointers
    aod->GetStdContent();

    //print event info
    aod->GetHeader()->Print();

    // loop over tracks
    Int_t nTracks = aod->GetNTracks();
    for (Int_t nTr = 0; nTr < nTracks; nTr++) {
      
      // print track info
      cout << nTr+1 << "/" << nTracks << ": track pt: " << aod->GetTrack(nTr)->Pt() << ", vertex x of this track: " << aod->GetTrack(nTr)->GetProdVertex()->GetX() << endl;
    }

    // loop over vertices
    Int_t nVtxs = aod->GetNVertices();
    for (Int_t nVtx = 0; nVtx < nVtxs; nVtx++) {
      
      // print track info
      cout << nVtx+1 << "/" << nVtxs << ": vertex z position: " << aod->GetVertex(nVtx)->GetZ() << endl;
    }
  }
  
  return;
}
