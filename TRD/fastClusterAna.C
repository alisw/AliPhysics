void fastClusterAna() {

/////////////////////////////////////////////////////////////////////////
//
// Example macro for the analysis of the TRD cluster 
//
/////////////////////////////////////////////////////////////////////////

  // Dynamically link some shared libs
  if (gClassTable->GetID("AliRun") < 0) {
    gROOT->LoadMacro("loadlibs.C");
    loadlibs();
  }

  // Input file name
  Char_t *alifile = "galice_c_v0.root"; 

  // Event number
  Int_t   nEvent  = 0;

  // Define the objects
  AliTRDv1      *TRD;
  TClonesArray  *TRDCluster;
  AliTRDcluster *OneTRDcluster;

  TH1F *hZ = new TH1F("hZ","Cluster z-position",700,-350.0,350.0);

  // Connect the AliRoot file containing Geometry, Kine, Hits, and Digits
  TFile *gafl = (TFile*) gROOT->GetListOfFiles()->FindObject(alifile);
  if (!gafl) {
    cout << "Open the ALIROOT-file " << alifile << endl;
    gafl = new TFile(alifile);
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
  cout << "nparticles = " << nparticles << endl;
  if (nparticles <= 0) break;
  
  // Get the pointer to the tree
  TTree *ClusterTree = gAlice->TreeD();

  // Get the pointer to the detector classes
  TRD = (AliTRDv1 *) gAlice->GetDetector("TRD");
  // Get the pointer to the hit container
  if (TRD) TRDCluster = TRD->Cluster();

  // Reconstruct the address
  ClusterTree->GetBranch("TRDcluster")->SetAddress(&TRDCluster);

  Int_t nEntries = ClusterTree->GetEntries();
  cout << "Number of entries in cluster tree = " << nEntries << endl; 

  // Loop through all entries in the tree
  Int_t nbytes;
  for (Int_t iEntry = 0; iEntry < nEntries; iEntry++) {

    cout << "iEntry = " << iEntry << endl;

    // Import the tree
    gAlice->ResetDigits();
    nbytes += ClusterTree->GetEvent(iEntry);

    // Get the number of digits in the detector 
    Int_t nTRDCluster = TRDCluster->GetEntriesFast();
    cout << " nTRDCluster = " << nTRDCluster << endl;    

    // Loop through all TRD digits
    for (Int_t iTRDCluster = 0; iTRDCluster < nTRDCluster; iTRDCluster++) {

      // Get the information for this digit
      OneTRDcluster = (AliTRDcluster*) TRDCluster->UncheckedAt(iTRDCluster);
      hZ->Fill(OneTRDcluster->fZ);

    }

  }

  hZ->Draw();

}
