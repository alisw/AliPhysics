// Macro for testing the new C++ reconstruction code

// Arguments:
//   FirstEvent (default 0)
//   LastEvent (default 0)
//   RecGeantHits (1 to reconstruct GEANT hits) (default 0)
//   FileName (for signal) (default "galice.root")
//   BkgGeantFileName (for background),
//      needed only if RecGeantHits = 1 and background to be added
void MUONreco (Int_t FirstEvent = 0, Int_t LastEvent = 0, Int_t RecGeantHits = 0, Text_t *FileName = "galice.root", Text_t *BkgGeantFileName = "")
{
  //
  cout << "MUONtest_reco" << endl;
  cout << "FirstEvent " << FirstEvent << endl;
  cout << "LastEvent " << LastEvent << endl;
  cout << "RecGeantHits " << RecGeantHits << endl;
  cout << "FileName ``" << FileName << "''" << endl;
  cout << "BkgGeantFileName ``" << BkgGeantFileName << "''" << endl;
  // Dynamically link some shared libs                    
  if (gClassTable->GetID("AliRun") < 0) {
    gROOT->LoadMacro("loadlibs.C");
    loadlibs();
  }

  // Connect the Root Galice file containing Geometry, Kine and Hits
  TFile *file = (TFile*) gROOT->GetListOfFiles()->FindObject(FileName);
  if (!file) {
    printf("\n Creating file %s\n", FileName);
    file = new TFile(FileName);
  }
  else printf("\n File %s found in file list\n", FileName);

  // Get AliRun object from file or create it if not on file
  if (!gAlice) {
    gAlice = (AliRun*) file->Get("gAlice");
    if (gAlice) printf("AliRun object found on file\n");
    if (!gAlice) {
      printf("\n Create new gAlice object");
      gAlice = new AliRun("gAlice","Alice test program");
    }
  }

  // Initializations
  // AliMUON *MUON  = (AliMUON*) gAlice->GetModule("MUON"); // necessary ????
  AliMUONEventReconstructor *Reco = new AliMUONEventReconstructor();
  Reco->SetRecGeantHits(RecGeantHits);

  // The right place for changing AliMUONEventReconstructor parameters
  // with respect to the default ones
//   Reco->SetMaxSigma2Distance(100.0);
  Reco->SetPrintLevel(10);
  cout << "AliMUONEventReconstructor: actual parameters" << endl;
  Reco->Dump();

  // Loop over events
  for (Int_t event = FirstEvent; event <= LastEvent; event++) {
    cout << "Event: " << event << endl;
    AliMUON *MUON  = (AliMUON*) gAlice->GetModule("MUON"); // necessary ????
    Int_t nparticles = gAlice->GetEvent(event);
    cout << "nparticles: " << nparticles << endl;
    // prepare background file and/or event if necessary
    if (RecGeantHits == 1) {
      if (event == FirstEvent) Reco->SetBkgGeantFile(BkgGeantFileName);
      if (Reco->GetBkgGeantFile())Reco->NextBkgGeantEvent();
    }
    Reco->EventReconstruct();
    // Dump current event
    Reco->EventDump();
  } // Event loop
}
