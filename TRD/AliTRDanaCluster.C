void AliTRDanaCluster() 
{

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
  Char_t *alifile = "galice.root"; 

  // Event number
  Int_t   nEvent  = 0;

  // Define the histograms
  TH1F *hCharge = new TH1F("hCharge","Cluster charge",100,0.0,1000.0);

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
  gAlice = (AliRun*) gafl->Get("gAlice");
  if (gAlice) {
    cout << "AliRun object found on file" << endl;
  }
  else {
    gAlice = new AliRun("gAlice","Alice test program");
  }

  // Import the Trees for the event nEvent in the file
  Int_t nparticles = gAlice->GetEvent(nEvent);
  if (nparticles <= 0) break;
  
  // Get the pointer to the hit-tree
  Char_t treeName[14];
  sprintf(treeName,"TreeR%d_TRD",nEvent);
  TTree          *clusterTree  = gafl->Get(treeName);
  clusterTree->Print();
  // Get the pointer to the detector classes
  AliTRDv1       *trd          = (AliTRDv1*) gAlice->GetDetector("TRD");
  // Get the geometry
  AliTRDgeometry *geo          = trd->GetGeometry();
  // Get the pointer to the hit container
  TObjArray      *clusterArray = trd->RecPoints();
  // Set the branch address
  clusterTree->GetBranch("TRDcluster")->SetAddress(&clusterArray);

  Int_t nEntries = clusterTree->GetEntries();
  cout << "nEntries = " << nEntries << endl;

  // Loop through all entries in the tree
  Int_t nbytes;
  for (Int_t iEntry = 0; iEntry < nEntries; iEntry++) {

    // Import the tree
    nbytes += clusterTree->GetEvent(iEntry);

    // Get the number of points in the detector 
    Int_t ncluster = clusterArray->GetEntriesFast();

    // Loop through all TRD digits
    for (Int_t icluster = 0; icluster < ncluster; icluster++) {

      // Get the information for this digit
      AliTRDcluster *cluster = (AliTRDcluster *) clusterArray->UncheckedAt(icluster);
      Int_t    detector = cluster->GetDetector();      
      Int_t    sector   = geo->GetSector(detector);
      Int_t    plane    = geo->GetPlane(detector);
      Int_t    chamber  = geo->GetChamber(detector);
      Float_t  charge   = cluster->GetQ();

      hCharge->Fill(charge);

    }

  }

  TCanvas *c1 = new TCanvas("c1","Cluster",50,50,600,400);
  hCharge->Draw();

}
