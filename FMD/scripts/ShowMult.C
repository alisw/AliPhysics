//
// Script to dump multiplicity information to std::cout. 
//
void
ShowMult()
{
  AliRunLoader* runLoader = AliRunLoader::Open("galice.root");
  runLoader->LoadgAlice();
  runLoader->LoadHeader();
  gAlice                   = runLoader->GetAliRun();
  AliFMD*       fmd        = static_cast<AliFMD*>(gAlice->GetDetector("FMD"));
  AliLoader*    fmdLoader  = runLoader->GetLoader("FMDLoader");
  fmdLoader->LoadRecPoints("READ");
  
  Int_t nEvents = runLoader->TreeE()->GetEntries();
  for (Int_t event = 0; event < nEvents; event++) {
    cout << "Event # " << event << endl;
    runLoader->GetEvent(event);
    TClonesArray* multStrips  = 0;
    TClonesArray* multRegions = 0;
    TTree*        treeR  = fmdLoader->TreeR();
    TBranch*      branchRegions = treeR->GetBranch("FMDPoisson");
    TBranch*      branchStrips  = treeR->GetBranch("FMDNaiive");
    branchRegions->SetAddress(&multRegions);
    branchStrips->SetAddress(&multStrips);

    Int_t total = 0;
    Int_t nEntries  = treeR->GetEntries();
    for (Int_t entry = 0; entry < nEntries; entry++) {
      cout << " Entry # " << entry << endl;
      treeR->GetEntry(entry);

      Int_t nMults = multStrips->GetLast();
      for (Int_t i = 0; i < nMults; i++) {
	// cout << "  Digit # " << i << endl;
	AliFMDMultStrip* mult = 
	  static_cast<AliFMDMultStrip*>(multStrips->UncheckedAt(i));
	if (mult->Particles() > 0) mult->Print();
      }

      nMults = multRegions->GetLast();
      for (Int_t i = 0; i < nMults; i++) {
	// cout << "  Digit # " << i << endl;
	AliFMDMultStrip* mult = 
	  static_cast<AliFMDMultStrip*>(multRegions->UncheckedAt(i));
	if (mult->Particles() > 0) mult->Print();
      }
    }
  }
}
