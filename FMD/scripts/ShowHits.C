//
// Script to dump hit information to std::cout. 
//
ShowHits()
{
  AliRunLoader* runLoader = AliRunLoader::Open("galice.root");
  runLoader->LoadgAlice();
  runLoader->LoadHeader();
  gAlice                   = runLoader->GetAliRun();
  AliFMD*       fmd        = static_cast<AliFMD*>(gAlice->GetDetector("FMD"));
  AliLoader*    fmdLoader  = runLoader->GetLoader("FMDLoader");
  fmdLoader->LoadHits("READ");
  
  Int_t nEvents = runLoader->TreeE()->GetEntries();
  for (Int_t event = 0; event < nEvents; event++) {
    cout << "Event # " << event << endl;
    runLoader->GetEvent(event);
    TClonesArray* hits   = 0;
    TTree*        treeH  = fmdLoader->TreeH();
    TBranch*      branch = treeH->GetBranch("FMD");
    branch->SetAddress(&hits);

    Int_t total = 0;
    Int_t nEntries  = treeH->GetEntries();
    for (Int_t entry = 0; entry < nEntries; entry++) {
      // cout << " Entry # " << entry << endl;
      treeH->GetEntry(entry);

      Int_t nHits = hits->GetEntries();
      for (Int_t i = 0; i < nHits; i++) {
	cout << "  Hit # " << i << "/" << nHits << "\t" << flush;
	AliFMDHit* hit = static_cast<AliFMDHit*>(hits->UncheckedAt(i));
	hit->Print();
	total++;
      }
    }
    cout << "Total number of hits: " << total << endl;
  }
}
