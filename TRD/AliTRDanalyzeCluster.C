Int_t AliTRDanalyzeCluster()
{
  //
  // Analyzes the cluster
  //

  Int_t rc = 0;

  if (!gAlice) {
    cout << "<AliTRDanalyzeCluster> No AliRun object found" << endl;
    rc = 1;
    return rc;
  }
  gAlice->GetEvent(0);

  // Get the pointer to the TRD detector 
  AliTRD *TRD = (AliTRD *) gAlice->GetDetector("TRD");
  if (!TRD) {
    cout << "<AliTRDanalyzeCluster> No TRD detector found" << endl;
    rc = 2;
    return rc;
  }

  // Define the histograms
  TH1F *hClusAll   = new TH1F("hClusAll"  ,"Amplitude of the cluster (all)"     ,501,-0.5,500.5);
  TH1F *hClusNoise = new TH1F("hClusNoise","Amplitude of the cluster (noise)"   ,501,-0.5,500.5);
  TH1F *hClusEl    = new TH1F("hClusEl"   ,"Amplitude of the cluster (electron)",501,-0.5,500.5);
  TH1F *hClusPi    = new TH1F("hClusPi"   ,"Amplitude of the cluster (pion)"    ,501,-0.5,500.5);

  // Get the pointer to the geometry object
  AliTRDgeometry *TRDgeometry;
  if (TRD) {
    TRDgeometry = TRD->GetGeometry();
  }
  else {
    cout << "<AliTRDanalyzeCluster> No TRD geometry found" << endl;
    rc = 3;
    return rc;
  }

  // Get the pointer to the hit-tree
  TFile *file = (TFile *) gROOT->GetListOfFiles()->FindObject("TRD_test.root");
  TTree *RecTree = file->Get("ClusterTree");
  if (!(RecTree)) {
    cout << "<AliTRDanalyzeCluster> No tree with rec. points found" << endl;
    rc = 4;
    return rc;
  }

  // Get the pointer to the hit container
  TObjArray *RecPointArray = TRD->RecPoints();
  if (!(RecPointArray)) {
    cout << "<AliTRDanalyzeCluster> No RecPointArray found" << endl;
    rc = 5;
    return rc;
  }

  // Set the branch address
  RecTree->GetBranch("TRDrecPoints")->SetAddress(&RecPointArray);
  Int_t nEntries = RecTree->GetEntries();
  cout << "Number of entries in the rec. point tree = " << nEntries << endl;

  Int_t countCluster = 0;

  // Loop through all entries in the tree
  Int_t nbytes;
  for (Int_t iEntry = 0; iEntry < nEntries; iEntry++) {

    // Import the tree
    nbytes += RecTree->GetEvent(iEntry);

    // Get the number of points in the detector 
    Int_t nRecPoint = RecPointArray->GetEntriesFast();

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
      Int_t    track0   = RecPoint->GetTrackIndex(0);
      Int_t    track1   = RecPoint->GetTrackIndex(1);
      TParticle *Part = 0;
      if (track0 > -1) {
        Part = gAlice->Particle(track0);
      }

      countCluster++;

      // Total spectrum
      hClusAll->Fill((Float_t) energy);

      // Noise spectrum
      if (track0 < 0) {
        hClusNoise->Fill((Float_t) energy);
      }          

      // Electron cluster
      if ((Part) && (Part->GetPdgCode() ==   11) && (track1 < 0)) {
        hClusEl->Fill((Float_t) energy);
      }

      // Pion cluster
      if ((Part) && (Part->GetPdgCode() == -211) && (track1 < 0)) {
        hClusPi->Fill((Float_t) energy);
      }

    }

  }

  cout << "<AliTRDanalyzeCluster> Found " << countCluster << " cluster in total" << endl;

  TCanvas *cCluster = new TCanvas("cCluster","AliTRDanalyzeCluster",50,50,600,600);
  cCluster->Divide(2,2);

  TF1 *fun;
  cCluster->cd(1);
  gPad->SetLogy();
  hClusAll->Fit("landau","0");
  fun = (TF1 *) hClusAll->GetListOfFunctions()->First();
  Float_t meanAll = fun->GetParameter(1);
  hClusAll->Draw();
  fun->SetLineColor(2);
  fun->Draw("SAME");
  
  cCluster->cd(2);
  gPad->SetLogy();
  gPad->SetLogx();
  Float_t meanNoise = hClusNoise->GetMean();
  hClusNoise->Draw();

  cCluster->cd(3);
  gPad->SetLogy();
  hClusEl->Fit("landau","0");
  fun = (TF1 *) hClusEl->GetListOfFunctions()->First();
  fun->SetLineColor(2);
  Float_t meanEl = fun->GetParameter(1);
  hClusEl->Draw();
  fun->Draw("SAME");

  cCluster->cd(4);
  gPad->SetLogy();
  hClusPi->Fit("landau","0");
  fun = (TF1 *) hClusPi->GetListOfFunctions()->First();
  fun->SetLineColor(2);
  Float_t meanPi = fun->GetParameter(1);
  hClusPi->Draw();
  fun->Draw("SAME");

  cout << endl;
  cout << "##################################################################" << endl;
  cout << "    Mean all       = " << meanAll   << endl;
  cout << "    Mean noise     = " << meanNoise << endl;
  cout << "    Mean electrons = " << meanEl    << endl;
  cout << "    Mean pions     = " << meanPi    << endl;
  cout << "##################################################################" << endl;
  cout << endl;

  return rc;

}
