void fastClusterCreate() {

///////////////////////////////////////////////////////////////////////// 
//
// Creates cluster from the hit information (fast simulator). 
// An additional hit-tree is added to the input file.
//
///////////////////////////////////////////////////////////////////////// 

  // Dynamically link some shared libs
  if (gClassTable->GetID("AliRun") < 0) {
    gROOT->LoadMacro("loadlibs.C");
    loadlibs();
  }

  // Input (and output) file name
  Char_t *alifile = "galice_c_v0.root";

  // Event number
  Int_t   nEvent  = 0;

  // Connect the AliRoot file containing Geometry, Kine, Hits, and Digits
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

  // Get the pointer to the detector classes
  AliTRDv0 *TRD = (AliTRDv0*) gAlice->GetDetector("TRD");

  // Create the clusters
  TRD->Hits2Clusters();

  // Write the new tree into the input file
  cout << "Entries in digits tree = " << gAlice->TreeD()->GetEntries() << endl;
  Char_t treeName[7];
  sprintf(treeName,"TreeD%d",nEvent);
  gAlice->TreeD()->Write(treeName);

}
