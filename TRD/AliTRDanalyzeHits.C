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
  AliDetector *TRD = gAlice->GetDetector("TRD");
  if (!TRD) {
    cout << "<AliTRDanalyzeHits> No TRD detector found" << endl;
    rc = 2;
    return rc;
  }

  // Define the histograms
  TH1F *hQdedx = new TH1F("hQdedx","Charge dedx-hits",100,0.0,1000.0);
  TH1F *hQtr   = new TH1F("hQtr"  ,"Charge TR-hits"  ,100,0.0,1000.0);
  TH2F *hZY    = new TH2F("hZY"   ,"Y vs Z",50,-100.0,100.0,40,290.0,370.0);
  TH2F *hXZ    = new TH2F("hXZ"   ,"Z vs X",50,-100.0,100.0,50,-100.0,100.0);

  // Get the pointer hit tree
  TTree *HitTree = gAlice->TreeH();  
  if (!HitTree) {
    cout << "<AliTRDanalyzeHits> No hit tree found" << endl;
    rc = 3;
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

    // Get the number of hits in the TRD created by this particle
    Int_t nHit = TRD->Hits()->GetEntriesFast();
    //cout << "<AliTRDanalyzeHits> Found " << nHit
    //     << " hits for primary particle " << iTrack << endl;

    // Loop through the TRD hits  
    for (Int_t iHit = 0; iHit < nHit; iHit++) {

      countHits++;

      AliTRDhit *hit = (AliTRDhit *) TRD->Hits()->UncheckedAt(iHit);

      Float_t x = hit->X();
      Float_t y = hit->Y();
      Float_t z = hit->Z();
      Float_t q = hit->GetCharge();

      if      (q > 0) { 
        hQdedx->Fill(q);
      }
      else if (q < 0) {
        hQtr->Fill(TMath::Abs(q));
      }

      hZY->Fill(z,y);
      hXZ->Fill(x,z);

    }

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

  return rc;

}
