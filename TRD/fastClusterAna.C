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
  Char_t *alifile = "galice_r_v0.root"; 

  // Event number
  Int_t   nEvent  = 0;

  TH2F *HLocal  = new TH2F("HLocal" ,"rec. points local row/col-position"
                                    ,21,-0.5,20.5,81,-0.5,80.5);
  TH2F *HGlobal = new TH2F("HGlobal","rec. points global x/y-position"   
                                    ,800,-400,400,800,-400,400);

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
  cout << "nparticles = " << nparticles << endl;
  if (nparticles <= 0) break;
  
  // Get the pointer to the tree
  TTree          *RecTree       = gAlice->TreeR();
  RecTree->Print();
  // Get the pointer to the detector classes
  AliTRDv0       *TRD           = (AliTRDv0*) gAlice->GetDetector("TRD");
  // Get the geometry
  AliTRDgeometry *TRDgeometry   = TRD->GetGeometry();
  // Get the pointer to the hit container
  TObjArray      *RecPointArray = TRD->RecPoints();
  // Set the branch address
  RecTree->GetBranch("TRDrecPoints")->SetAddress(&RecPointArray);

  Int_t nEntries = RecTree->GetEntries();
  cout << "Number of entries in reconstruction tree = " << nEntries << endl; 

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
      Float_t  row      = RecPoint->GetLocalRow();
      Float_t  col      = RecPoint->GetLocalCol();

      Int_t    sector   = TRDgeometry->GetSector(detector);
      Int_t    plane    = TRDgeometry->GetPlane(detector);
      Int_t    chamber  = TRDgeometry->GetChamber(detector);

      TVector3 Pos;
      TMatrix  Cov;
      RecPoint->GetGlobalPosition(Pos,Cov);
      HGlobal->Fill(Pos.X(),Pos.Y());
      if ((sector == 17) && (plane == 0) && (chamber == 2)) {
        HLocal->Fill(row,col);
      }

    }

  }

  TCanvas *C = new TCanvas("C","recPoints",10,10,400,600);
  C->Divide(1,2);
  C->cd(1);
  HLocal->Draw("BOX");
  C->cd(2);
  HGlobal->Draw("BOX");

}
