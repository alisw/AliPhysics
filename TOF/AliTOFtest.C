Int_t AliTOFtest() 
{
  //
  // Test macro for the TOF code
  // report bug to Fabrizio.Pierella@cern.ch
  // Use case:
  // start aliroot
  // root [0] .L AliTOFtest.C
  // root [1] AliTOFtest()

  Int_t rc = 0;

  // Initialize the test setup 

  //gAlice->Init("$(ALICE_ROOT)/TOF/AliTOFconfig.C");
  gAlice->Init("$ALICE_ROOT/TOF/AliTOFconfig.C");

  // Run one central Hijing event and create the hits (time required: 
  // some minuts)

  gAlice->SetDebug(2);
  gAlice->Run(1);
  
  if (gAlice) delete gAlice;
  TFile *file = (TFile *) gROOT->GetListOfFiles()->FindObject("TOF_test.root");
  gAlice = (AliRun *) file->Get("gAlice");

  // Analyze the TOF hits
  if (rc = AliTOFanalyzeHits()) return rc;

  return rc;

}

//_____________________________________________________________________________
Int_t AliTOFanalyzeHits()
{
  //
  // Analyzes the hits and fills QA-histograms 
  //

  Int_t rc = 0;

  if (!gAlice) {
    cout << "<AliTOFanalyzeHits> No AliRun object found" << endl;
    rc = 1;
    return rc;
  }
  gAlice->GetEvent(0);

  // Get the pointer to the TOF detector 
  AliTOF *tof = (AliTOF *) gAlice->GetDetector("TOF");
  if (!tof) {
    cout << "<AliTOFanalyzeHits> No TOF detector found" << endl;
    rc = 2;
    return rc;
  }

  // Define the histograms
  // x,y,z, rho, tof, padx, padz, sector, plate, strip, (x vs y)

  // hit-map in a plane
  TH2F *h2hitMap = new TH2F("h2hitMap","Hit Map (projection on the plane)",2500,-12500.,12500.,800,-400.,400.);
  // time of flight distribution for primaries and secondaries
  TH1F *htofp    = new TH1F("htofp","Time of Flight (primaries)",800,0.,80.);
  TH1F *htofs    = new TH1F("htofs","Time of Flight (secondaries)",800,0.,80.);
  // momentum when striking the TOF for primaries and secondaries
  TH1F *htofmomp = new TH1F("htofmomp","Momentum at TOF (primaries)",100,0.,10.);
  TH1F *htofmoms = new TH1F("htofmoms","Momentum at TOF (secondaries)",100,0.,10.);
  // TOF hit volumes
  TH1F *hsector  = new TH1F("hsector","Sector",20,0.,20.);
  TH1F *hplate   = new TH1F("hplate","Plate ", 6,0., 6.);
  TH1F *hstrip   = new TH1F("hstrip","Strip ",25,0.,25.);
  TH1F *hpadz    = new TH1F("hpadz","Pad along z ",3,0.,3.);
  TH1F *hpadx    = new TH1F("hpadx","Pad along x",50,0.,50.);
  // track length when striking the TOF (used by AliTOFT0)
  TH1F *htrackLenp= new TH1F("htrackLenp","Track Length on TOF for Primaries",800,0.,800.);

  // Get the pointer hit tree
  TTree *hitTree = gAlice->TreeH();  
  if (!hitTree) {
    cout << "<AliTOFanalyzeHits> No hit tree found" << endl;
    rc = 4;
    return rc;
  }

  Int_t countHits = 0;
  Int_t nBytes    = 0;

  // Get the number of entries in the hit tree
  // (Number of primary particles creating a hit somewhere)
  Int_t nTrack = (Int_t) hitTree->GetEntries();
  cout << "<AliTOFanalyzeHits> Found " << nTrack 
       << " primary particles with hits" << endl;

  Int_t nPrimaryOnTof = 0;
  Int_t nSecondaryOnTof = 0;
  Int_t nelectron  = 0;
  Int_t npion      = 0;
  Int_t nkaon      = 0;
  Int_t nproton    = 0;    
  Int_t nmuon      = 0;
  
  // Loop through all entries in the tree
  for (Int_t iTrack = 0; iTrack < nTrack; iTrack++) {

    gAlice->ResetHits();
    nBytes += hitTree->GetEvent(iTrack);


    // Loop through the TOF hits  
    Int_t iHit = 0;
    AliTOFhitT0 *hit = (AliTOFhitT0 *) tof->FirstHit(-1);
    while (hit) {

      countHits++;
      iHit++;

      Float_t x     = hit->X();
      Float_t y     = hit->Y();
      Float_t z     = hit->Z();
      Float_t circleLen=TMath::Sqrt(x*x+y*y)*TMath::ATan2(y,x);
      h2hitMap->Fill(circleLen,z);

      Float_t flightTime = hit->GetTof(); // [s]
      flightTime*=1.e+09; // convert in [ns]
      Float_t angle = hit->GetIncA();
      Float_t tofmom= hit->GetMom(); // [GeV/c]
      Float_t trackLen= hit->GetLen(); // [cm]

      // TOF hit volumes
      Int_t sector    = hit->GetSector(); // range [1-18]
      Int_t plate     = hit->GetPlate();  // range [1- 5]
      Int_t strip     = hit->GetStrip();  // range [1-20]
      Int_t padz      = hit->GetPadz();   // range [1- 2]
      Int_t padx      = hit->GetPadx();   // range [1-48]
      // it is QA, then I perform QA!
      Bool_t isHitBad = (sector<1 || sector>18 || plate<1 || plate >5 || padz<1 || padz>2 || padx<1 || padx>48);

      if (isHitBad) {
	cout << "<AliTOFanalyzeHits>  strange hit found" << endl;
	rc = 3;
	return rc;
      }
      // filling hit volume histos
      hsector->Fill(sector);
      hplate->Fill(plate);
      hstrip->Fill(strip);
      hpadx->Fill(padx);
      hpadz->Fill(padz);


      Int_t track     = hit->Track();
      TParticle *part = gAlice->Particle(track);

      // getting MC info for the current track
      if (part->GetFirstMother()<0){
	Int_t icod = TMath::Abs(part->GetPdgCode());
	switch (icod) {
	case 211:
	  npion++;
	  break ;
	case 321:
	  nkaon++;
	  break ;
	case 2212:
	  nproton++;
	  break ;
	case 11:
	  nelectron++;
	  break ;
	case 13:
	  nmuon++;
	  break ;
	}
	htofp->Fill(flightTime);
	htofmomp->Fill(tofmom);
	htrackLenp->Fill(trackLen);
      } else {
	htofs->Fill(flightTime);
	htofmoms->Fill(tofmom);
      }

      // go to next hit
      hit = (AliTOFhitT0 *) tof->NextHit();         

    }

  }

  cout << "<AliTOFanalyzeHits> Found " << countHits << " hits in total" << endl;
  cout << npion     << " primary pions reached the TOF detector"     << endl;
  cout << nkaon     << " primary kaons reached the TOF detector"     << endl;
  cout << nproton   << " primary protons reached the TOF detector"   << endl;
  cout << nelectron << " primary electrons reached the TOF detector" << endl;
  cout << nmuon     << " primary muons reached the TOF detector"     << endl;

  
  TCanvas *cHits = new TCanvas("cHits","AliTOFanalyzeHits hit volumes",50,50,900,900);
  cHits->Divide(3,2);
  cHits->cd(1);
  hsector->Draw();
  cHits->cd(2);
  hplate->Draw();
  cHits->cd(3);
  hstrip->Draw();
  cHits->cd(4);
  hpadz->Draw();
  cHits->cd(5);
  hpadx->Draw();

  TCanvas *chitmap = new TCanvas("chitmap","AliTOFanalyzeHits Hit Map",50,50,600,600);
  chitmap->cd();
  h2hitMap->Draw();

  TCanvas *ctrackLen = new TCanvas("ctrackLen","AliTOFanalyzeHits Track Length for primaries on TOF",50,50,400,400);
  ctrackLen->cd();
  htrackLenp->Draw();

  TCanvas *ctofmom = new TCanvas("ctofmom","AliTOFanalyzeHits flight times",50,50,700,700);
  ctofmom->Divide(2,2);
  ctofmom->cd(1);
  gPad->SetLogy();
  htofp->Draw();
  ctofmom->cd(2);
  gPad->SetLogy();
  htofs->Draw();
  ctofmom->cd(3);
  gPad->SetLogy();
  htofmomp->Draw();
  ctofmom->cd(4);
  gPad->SetLogy();
  htofmoms->Draw();
  

  // save histos into file TOF_hitsQA.root
  TFile *fout = new TFile("TOF_hitsQA.root","RECREATE");
  h2hitMap->Write();
  htofp->Write();
  htofs->Write();
  htofmomp->Write();
  htofmoms->Write();
  hsector->Write();
  hplate->Write();
  hstrip->Write();
  hpadz->Write();
  hpadx->Write();
  htrackLenp->Write();
  fout->Close(); 

  return rc;

}
