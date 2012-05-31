void bcmhits(Int_t events  = 10)
{
    TH1::AddDirectory(kFALSE);
    zH = new TH1F("zH", "z coordinate of hit", 210., -2100., 2100.);
    rH = new TH1F("rH", "r coordinate of hit", 200.,     0.,   20.);
    tH = new TH1F("tH", "time of hits",        100.,     0.,  100.);
    
  // Creating Run Loader and openning file containing Hits
    TClonesArray* hits   = new TClonesArray ("AliBCMHit", 1000);    
    AliRunLoader* rl = AliRunLoader::Open("galice.root");
    if (rl ==0x0) {
	printf(">>> Error : Error Opening file \n");
	return;
    }
    rl->LoadgAlice();
    rl->LoadHeader();
    rl->LoadKinematics();
    gAlice = rl->GetAliRun();
    AliBCM*  pBCM  = (AliBCM*) gAlice->GetDetector("BCM");

// Histos 
    for (Int_t nev = 0; nev < events; nev++) {
	// Load event
	rl->GetEvent(nev);
	AliLoader* bcmL = rl->GetLoader("BCMLoader");
//	bcmL->InitDefaults();
	bcmL->LoadHits();
	TTree* treeH = bcmL->TreeH();
	Int_t ntracks = (Int_t) treeH->GetEntries();
	printf("Number of particles %6d %6d \n", nev, ntracks);
	for (Int_t itr = 0; itr < ntracks; itr++) {
	    gAlice->ResetHits ();
	    treeH->GetEvent(itr);
	    hits = pBCM->Hits();
	    Int_t nhits = hits->GetEntriesFast();
	    for (Int_t hit = 0; hit < nhits; hit++)
	    {
		
		AliBCMHit*  mHit = (AliBCMHit *) hits->UncheckedAt(hit);
		Float_t x = mHit->X();
		Float_t y = mHit->Y();
		Float_t z = mHit->Z();
		Float_t r = TMath::Sqrt(x * x + y * y);
		
		zH->Fill(z);
		tH->Fill(mHit->Time() * 1.e9);
		rH->Fill(r);
		
		
	    } // hits
	} // tracks
    } // events 
    new TCanvas("c1");
    zH->Draw();
   
    new TCanvas("c2");
    tH->Draw();

    new TCanvas("c3");
    rH->Draw();
    
}
