void digitsCreate () {

/////////////////////////////////////////////////////////////////////////
//
// Creates the digits from the hit information. An additional hit-tree
// is added to the input file.
//
/////////////////////////////////////////////////////////////////////////

  // Dynamically link some shared libs
  if (gClassTable->GetID("AliRun") < 0) {
    gROOT->LoadMacro("loadlibs.C");
    loadlibs();
  }

  // Input (and output) file name
  Char_t *alifile = "galice_v1.root"; 

  // Event number
  Int_t   nEvent  = 0;

  // Connect the AliRoot file containing Geometry, Kine, and Hits
  TFile *gafl = (TFile*) gROOT->GetListOfFiles()->FindObject(alifile);
  if (!gafl) {
    cout << "Open the ALIROOT-file " << alifile << endl;
    gafl = new TFile(alifile,"UPDATE");
  }
  else {
    cout << alifile << " is already open" << endl;
  }

  // Get AliRun object from file or create it if not on file
  if (!gAlice) {
    gAlice = (AliRun*) gafl->Get("gAlice");
    if (gAlice)  
      cout << "AliRun object found on file" << endl;
    else
      gAlice = new AliRun("gAlice","Alice test program");
  }

  // Import the Trees for the event nEvent in the file
  Int_t nparticles = gAlice->GetEvent(nEvent);
  if (nparticles <= 0) break;

  // Get the pointer to the detector class
  AliTRDv2 *TRD = (AliTRDv2*) gAlice->GetDetector("TRD");

  // Create the digitd and fill the digits-tree
  TRD->Hits2Digits();

  // Write the new tree into the input file
  cout << "Entries in hit tree = " << gAlice->TreeD()->GetEntries()) << endl;
  Char_t treeName[7];
  sprintf(treeName,"TreeD%d",nEvent);
  gAlice->TreeD()->Write(treeName);

}
