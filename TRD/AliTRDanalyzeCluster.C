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
  TH1F *hClusAll   = new TH1F("hClusAll"  ,"Amplitude of the cluster (all)"     
                                          ,501,-0.5,500.5);
  TH1F *hClusNoise = new TH1F("hClusNoise","Amplitude of the cluster (noise)"   
                                          ,  5,-0.5,  4.5);
  TH1F *hClusEl    = new TH1F("hClusEl"   ,"Amplitude of the cluster (electron)"
                                          ,501,-0.5,500.5);
  TH1F *hClusPi    = new TH1F("hClusPi"   ,"Amplitude of the cluster (pion)"    
                                          ,501,-0.5,500.5);

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
  TTree *ClusterTree = file->Get("ClusterTree");
  if (!(ClusterTree)) {
    cout << "<AliTRDanalyzeCluster> No tree with clusters found" << endl;
    rc = 4;
    return rc;
  }

  // Get the pointer to the hit container
  TObjArray *ClusterArray = TRD->RecPoints();
  if (!(ClusterArray)) {
    cout << "<AliTRDanalyzeCluster> No ClusterArray found" << endl;
    rc = 5;
    return rc;
  }

  // Set the branch address
  ClusterTree->GetBranch("TRDcluster")->SetAddress(&ClusterArray);
  Int_t nEntries = ClusterTree->GetEntries();
  cout << "<AliTRDanalyzeCluster> Number of entries in the cluster tree = " 
       << nEntries 
       << endl;

  Int_t countCluster = 0;
  Int_t countUnfold  = 0;

  // Loop through all entries in the tree
  Int_t nbytes;
  for (Int_t iEntry = 0; iEntry < nEntries; iEntry++) {

    // Import the tree
    nbytes += ClusterTree->GetEvent(iEntry);

    // Get the number of points in the detector 
    Int_t nCluster = ClusterArray->GetEntriesFast();

    // Loop through all TRD digits
    for (Int_t iCluster = 0; iCluster < nCluster; iCluster++) {

      // Get the information for this digit
      AliTRDcluster *Cluster = (AliTRDcluster *) ClusterArray->UncheckedAt(iCluster);
      Int_t    detector = Cluster->GetDetector();      
      Int_t    sector   = TRDgeometry->GetSector(detector);
      Int_t    plane    = TRDgeometry->GetPlane(detector);
      Int_t    chamber  = TRDgeometry->GetChamber(detector);
      Float_t  energy   = Cluster->GetQ();
      Int_t    track0   = Cluster->GetTrackIndex(0);
      Int_t    track1   = Cluster->GetTrackIndex(1);
      Int_t    track2   = Cluster->GetTrackIndex(2);
      TParticle *Part = 0;
      if (track0 > -1) {
        Part = gAlice->Particle(track0);
      }

      countCluster++;
      if (Cluster->FromUnfolding()) countUnfold++;

      // Total spectrum
      hClusAll->Fill(energy);

      if (!Cluster->FromUnfolding()) {

        // Noise spectrum
        if (track0 < 0) {
          hClusNoise->Fill(energy);
        }          

        // Electron cluster
        if ((Part) && (Part->GetPdgCode() ==   11) && (track1 < 0)) {
          hClusEl->Fill(energy);
        }

        // Pion cluster
        if ((Part) && (Part->GetPdgCode() == -211) && (track1 < 0)) {
          hClusPi->Fill(energy);
        }

      }

    }

  }

  cout << "<AliTRDanalyzeCluster> Found " << countCluster << " cluster in total"       << endl;
  cout << "<AliTRDanalyzeCluster> Found " << countUnfold  << " cluster from unfolding" << endl;
  cout << endl;

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
