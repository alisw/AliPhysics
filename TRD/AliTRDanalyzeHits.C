Int_t AliTRDanalyzeHits()
{
  //
  // Analyzes the hits and fills QA-histograms 
  //

  Int_t rc = 0;

  if (!gAlice) {
    cout << "<AliTRDanalyzeHits> No AliRun object found" << endl;
    rc = 1;
    return rc;
  }
  gAlice->GetEvent(0);

  // Get the pointer to the TRD detector 
  AliTRD *TRD = (AliTRD *) gAlice->GetDetector("TRD");
  if (!TRD) {
    cout << "<AliTRDanalyzeHits> No TRD detector found" << endl;
    rc = 2;
    return rc;
  }

  // Get the pointer to the geometry object
  AliTRDgeometryFull *TRDgeometry;
  if (TRD) {
    TRDgeometry = (AliTRDgeometryFull *) TRD->GetGeometry();
  }
  else {
    cout << "<AliTRDanalyzeHits> No TRD geometry found" << endl;
    rc = 3;
    return rc;
  }

  // Define the histograms
  TH1F *hQdedx  = new TH1F("hQdedx","Charge dedx-hits",100,0.0,1000.0);
  TH1F *hQtr    = new TH1F("hQtr"  ,"Charge TR-hits"  ,100,0.0,1000.0);

  Float_t rmin   = TRDgeometry->Rmin();
  Float_t rmax   = TRDgeometry->Rmax();
  Float_t length = TRDgeometry->GetChamberLengthI(0);
  Float_t width  = TRDgeometry->GetChamberWidth(0);
  Int_t   ncol   = TRDgeometry->GetColMax(0);
  Int_t   nrow   = TRDgeometry->GetRowMax(0,2,13);
  Int_t   ntime  = ((Int_t) (rmax - rmin) / TRDgeometry->GetTimeBinSize());

  TH2F *hZY     = new TH2F("hZY"   ,"Y vs Z (chamber 0)", nrow,-length/2.,length/2.
                                                        ,ntime,      rmin,     rmax);
  TH2F *hXZ     = new TH2F("hXZ"   ,"Z vs X (plane 0)"  , ncol, -width/2., width/2.
                                                        , nrow,-length/2.,length/2.);

  TH1F *hNprimP = new TH1F("hNprimP","Number of primary electrons per cm (MIP pion)"
                                    ,50,0.0,100.0);
  TH1F *hNprimE = new TH1F("hNprimE","Number of primary electrons per cm (3GeV electron)"
                                    ,50,0.0,100.0);
  TH1F *hNtotP  = new TH1F("hNtotP" ,"Number of electrons per cm (MIP pion)"
                                    ,50,0.0,1000.0);
  TH1F *hNtotE  = new TH1F("hNtotE" ,"Number of electrons per cm (3GeV electron)"
                                    ,50,0.0,1000.0);

  // Get the pointer hit tree
  TTree *HitTree = gAlice->TreeH();  
  if (!HitTree) {
    cout << "<AliTRDanalyzeHits> No hit tree found" << endl;
    rc = 4;
    return rc;
  }

  Int_t countHits = 0;
  Int_t nBytes    = 0;

  // Get the number of entries in the hit tree
  // (Number of primary particles creating a hit somewhere)
  Int_t nTrack = (Int_t) HitTree->GetEntries();
  cout << "<AliTRDanalyzeHits> Found " << nTrack 
       << " primary particles with hits" << endl;

  // Loop through all entries in the tree
  for (Int_t iTrack = 0; iTrack < nTrack; iTrack++) {

    gAlice->ResetHits();
    nBytes += HitTree->GetEvent(iTrack);

    Int_t nPrimE = 0;
    Int_t nPrimP = 0;
    Int_t nTotE  = 0;
    Int_t nTotP  = 0;    

    // Get the number of hits in the TRD created by this particle
    Int_t nHit = TRD->Hits()->GetEntriesFast();

    // Loop through the TRD hits  
    for (Int_t iHit = 0; iHit < nHit; iHit++) {

      countHits++;

      AliTRDhit *hit = (AliTRDhit *) TRD->Hits()->UncheckedAt(iHit);

      Float_t x     = hit->X();
      Float_t y     = hit->Y();
      Float_t z     = hit->Z();
      Float_t q     = hit->GetCharge();
      Int_t   track = hit->Track();
      Int_t   det   = hit->GetDetector();
      Int_t   plane = TRDgeometry->GetPlane(det);

      if      (hit->FromDrift()) {
        hQdedx->Fill(q);
        hZY->Fill(z,y);
        if (plane == 0) {
          hXZ->Fill(x,z);
	}
      }
      else if (hit->FromTRphoton()) {
        hQtr->Fill(TMath::Abs(q));
      }

      TParticle *part = gAlice->Particle(track);

      if ((plane == 0) && (hit->FromDrift())) {

        // 3 GeV electrons
        if ((part->GetPdgCode() ==   11) && 
            (part->P()          >   2.9)) {
          nPrimE++;
          nTotE += ((Int_t) q);
        }

        // MIP pions
        if ((part->GetPdgCode() == -211) &&
            (part->P()          >   0.5)) {
          nPrimP++;
          nTotP += ((Int_t) q);
        }

      }

    }

    if (nPrimE > 0) hNprimE->Fill(((Double_t) nPrimE)/3.);
    if (nPrimP > 0) hNprimP->Fill(((Double_t) nPrimP)/3.);
    if (nTotE  > 0) hNtotE->Fill(((Double_t) nTotE)/3.);
    if (nTotP  > 0) hNtotP->Fill(((Double_t) nTotP)/3.);

  }

  cout << "<AliTRDanalyzeHits> Found " << countHits << " hits in total" << endl;

  TCanvas *cHits = new TCanvas("cHits","AliTRDanalyzeHits",50,50,600,600);
  cHits->Divide(2,2);
  cHits->cd(1);
  hXZ->Draw("COL");
  cHits->cd(2);
  hZY->Draw("COL");
  cHits->cd(3);
  gPad->SetLogy();
  hQdedx->Draw();
  cHits->cd(4);
  gPad->SetLogy();
  hQtr->Draw();

  TCanvas *cNel = new TCanvas("cNel","Ionization",50,50,600,600);
  cNel->Divide(2,2);
  cNel->cd(1);
  hNprimE->SetStats();
  hNprimE->Draw();
  cNel->cd(2);
  hNprimP->SetStats();
  hNprimP->Draw();
  cNel->cd(3);
  hNtotE->SetStats();
  hNtotE->Draw();
  cNel->cd(4);
  hNtotP->SetStats();
  hNtotP->Draw();

  TFile *fout = new TFile("TRD_ionization.root","RECREATE");
  hNprimE->Write();
  hNprimP->Write();
  hNtotE->Write();
  hNtotP->Write();
  fout->Close(); 

  return rc;

}
