void testvert() {

  TFile f("AliESDs.root");

  TTree * tree = (TTree*)f.Get("esdTree");

  AliESDEvent * esd = new AliESDEvent();// The signal ESD object is put here
  esd->ReadFromTree(tree);

  Int_t nev = tree->GetEntriesFast();
  
  for (Int_t iev=0; iev<nev; iev++) {
    cout << "---------- Signal event ----------" << iev << endl;

    // Get ESD
    tree->GetEntry(iev);

    AliESDVertex * vert = esd->GetPrimaryVertex();

    if (vert) cout << vert->GetTitle() << endl;

    if (strstr(vert->GetTitle(),"VertexerTracks")) cout << "OK" << endl;

  }
}
