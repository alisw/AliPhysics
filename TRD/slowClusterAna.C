void slowClusterAna() {

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
  Char_t *alifile = "galice_r_v1.root"; 

  // Event number
  Int_t   nEvent  = 0;

  // Define the histograms
  TH1F *hEnergy = new TH1F("hEnergy","Cluster energy",100,0.0,1000.0);

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
  if (gAlice)  
    cout << "AliRun object found on file" << endl;
  else
    gAlice = new AliRun("gAlice","Alice test program");

  // Import the Trees for the event nEvent in the file
  Int_t nparticles = gAlice->GetEvent(nEvent);
  if (nparticles <= 0) break;
  
  // Get the pointer to the hit-tree
  TTree          *RecTree       = gAlice->TreeR();
  RecTree->Print();
  // Get the pointer to the detector classes
  AliTRDv1       *TRD           = (AliTRDv1*) gAlice->GetDetector("TRD");
  // Get the geometry
  AliTRDgeometry *TRDgeometry   = TRD->GetGeometry();
  // Get the pointer to the hit container
  TObjArray      *RecPointArray = TRD->RecPoints();
  // Set the branch address
  RecTree->GetBranch("TRDrecPoints")->SetAddress(&RecPointArray);

  Int_t nEntries = RecTree->GetEntries();
  cout << "nEntries = " << nEntries << endl;

  // Loop through all entries in the tree
  Int_t nbytes;
  for (Int_t iEntry = 0; iEntry < nEntries; iEntry++) {

    cout << "iEntry = " << iEntry << endl;

    // Import the tree
    nbytes += RecTree->GetEvent(iEntry);

    // Get the number of points in the detector 
    Int_t nRecPoint = RecPointArray->GetEntriesFast();
    cout << " nRecPoint = " << nRecPoint << endl;    

    // Loop through all TRD digits
    for (Int_t iRecPoint = 0; iRecPoint < nRecPoint; iRecPoint++) {

      // Get the information for this digit
      AliTRDrecPoint *RecPoint = (AliTRDrecPoint *) RecPointArray->UncheckedAt(iRecPoint);
      Int_t    detector = RecPoint->GetDetector();      
      Int_t    sector   = TRDgeometry->GetSector(detector);
      Int_t    plane    = TRDgeometry->GetPlane(detector);
      Int_t    chamber  = TRDgeometry->GetChamber(detector);
      Int_t    energy   = RecPoint->GetEnergy();
      TVector3 pos;
      RecPoint->GetLocalPosition(pos);

      hEnergy->Fill((Float_t) energy);

    }

  }

  TCanvas *c1 = new TCanvas("c1","Cluster",50,50,600,400);
  hEnergy->Draw();

}
