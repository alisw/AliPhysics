//
// Script to digit multiplicity information to std::cout. 
//
void
ShowDigits()
{
  AliRunLoader* runLoader = AliRunLoader::Open("galice.root");
  runLoader->LoadgAlice();
  runLoader->LoadHeader();
  gAlice                   = runLoader->GetAliRun();
  AliFMD*       fmd        = static_cast<AliFMD*>(gAlice->GetDetector("FMD"));
  AliLoader*    fmdLoader  = runLoader->GetLoader("FMDLoader");
  fmdLoader->LoadDigits("READ");
  
  Int_t nEvents = runLoader->TreeE()->GetEntries();
  for (Int_t event = 0; event < nEvents; event++) {
    cout << "Event # " << event << endl;
    runLoader->GetEvent(event);
    TClonesArray* digits   = 0;
    TTree*        treeD  = fmdLoader->TreeD();
    TBranch*      branch = treeD->GetBranch("FMD");
    branch->SetAddress(&digits);

    Int_t total = 0;
    Int_t nEntries  = treeD->GetEntries();
    for (Int_t entry = 0; entry < nEntries; entry++) {
      cout << " Entry # " << entry << endl;
      treeD->GetEntry(entry);

      Int_t nDigits = digits->GetLast();
      for (Int_t i = 0; i < nDigits; i++) {
	// cout << "  Digit # " << i << endl;
	AliFMDDigit* digit = static_cast<AliFMDDigit*>(digits->UncheckedAt(i));
	if (digit->Counts() > 12) { 
	  digit->Print();
	  total++;
	}
      }
    }
    cout << "Total number of digits: " << total << endl;
  }
}
