void testEMCALGeom(const char *filename="galice.root",){
/////////////////////////////////////////////////////////////////////////
//   This macro computes the sampling fraction in the EMCAL
//   
//     Root > .L ComputeSamplingFraction.C //this loads the macro in memory
//     Root > ComputeSamplingFraction();   //by default process first event   
//     Root > ComputeSamplingFraction("galice2.root"); //process from the file 
//                                                       galice2.root .
//Begin_Html
/*
<img src="picts/ComputeSamplingFraction.gif">
*/
//End_Html
/////////////////////////////////////////////////////////////////////////
    if(gAlice){
	delete gAlice;
	gAlice=0;
    }else{
	// Dynamically link some shared libs
	if(gClassTable->GetID("AliRun") < 0) {
	    gROOT->LoadMacro("loadlibs.C");
	    loadlibs();
	} // end if
    } // end if gAlice
    // Connect the Root Galice file containing Geometry, Kine and Hits
    TFile *file = (TFile*)gROOT->GetListOfFiles()->FindObject(filename);
    if(!file) file = new TFile(filename);

    // Get AliRun object from file or create it if not on file
    if(!gAlice) {
	gAlice = (AliRun*)file->Get("gAlice");
	if(gAlice) printf("AliRun object found on file\n");
	if(!gAlice) gAlice = new AliRun("gAlice","Alice test program");
    } // end if !gAlice
      
    // Set event pointer to this event
    Int_t nparticles = gAlice->GetEvent(0);
    if (nparticles <= 0){
	cout << "No particles found for event " << evNumber;
	cout << " in file " << filename << endl;
	return;
    } // end if nparticles <=0

    // Pointer to specific detector hits.
    AliEMCALHit  *emcalHit;

    // Get pointers to ALL Alice detectors and Hits containers
    AliEMCALv0  *EMCAL  = (AliEMCALv0*)  gAlice->GetDetector("EMCAL");
    if(!EMCAL) {
	cout << "EMCAL not Found. Exiting." << endl;
	return;
    } // end if !EMCAL

    AliEMCALGeometry *g = EMCAL->GetGeometry();
    
    Int_t i,in,ix,ieta,iphi,ipre;
    Float_t eta,phi;
    cout << "i    ieta  iphi  ipre     in       eta          phi   ix" << endl;
    for(i=0;i<=2*g->GetNEta()*g->GetNPhi();i++){
	g->TowerIndexes(i,ieta,iphi,ipre);
	in = g->TowerIndex(ieta,iphi,ipre);
	g->EtaPhiFromIndex(i,eta,phi);
	ix = g->TowerIndexFromEtaPhi(eta,phi);
	cout << i << "     " << ieta << "     " << iphi << "     " << ipre;
	cout << "     " << in << "     " << eta << "     " << phi << "     " << ix << endl;
    } //
    return;
    // Get pointer to the particle
    TParticle    *part;
    Int_t         ipart;

    // Create histograms
    TH1F *hEMCAL = new TH1F("hEMCAL" ,"Energy",100,0.,2.);
    TH2F *hSfs   = new TH2F("hSfs","Sampling Fraction",200,-0.7,0.7,
			    200,0.0,120.0);
    TH2F *hSfi   = new TH2F("hSfi","Sampling Fraction Normilization",
			    200,-0.7,0.7,200,0.0,120.0);

    Int_t track,ntracks = gAlice->TreeH()->GetEntries();
    Float_t etot = 0.0;
    Float_t egam = 0.0;
    Bool_t first = kTRUE;
    Float_t eta,phi,x,y,z;
    // Start loop on tracks in the hits containers
    for(track=0; track<ntracks;track++){
	//MI change
	gAlice->ResetHits();
	gAlice->TreeH()->GetEvent(track);
	for(emcalHit=(AliEMCALHit*)EMCAL->FirstHit(-1);emcalHit;
	    emcalHit=(AliEMCALHit*)EMCAL->NextHit()) {
	    ipart = emcalHit->GetTrack();
	    part  = gAlice->Particle(ipart);
	    if(part->GetFirstMother()!=0) continue;
	    // only those that don't interact.
	    hEMCAL->Fill((Float_t)(emcalHit->GetEnergy()));
	    x = emcalHit->X();
	    y = emcalHit->Y();
	    z = emcalHit->Z();
	    findetaphi(x,y,z,eta,phi);
	    if(first){
		egam += part->Energy();
		hSfi->Fill(eta,phi,part->Energy());
	    } // end if first.
	    first = kFALSE;
	    etot += emcalHit->GetEnergy();
	    hSfs->Fill(eta,phi,emcalHit->GetEnergy());
	} // end for tpcHit
	first = kTRUE;
    } // end for track
    hSfs->Divide(hSfi);
    cout << "Energy of gamma =" << egam << "GeV. Energy in EMCAL =" << etot;
    cout << "GeV. Sampling Fractions is therefore f=" << etot/egam << endl;

    //Create a canvas, set the view range, show histograms
    TCanvas *c0 = new TCanvas("c0","Alice Detectors",400,10,600,700);
    c0->Divide(1,3);
    c0->cd(1);
    hEMCAL->SetFillColor(42);
    hEMCAL->Draw();
    c0->cd(2);
    hSfi->Draw("lego");
    c0->cd(3);
    hSfs->Draw("lego");
//    c0->SaveAs("analHitsEMCAL.eps");
}
void findetaphi(Float_t x,Float_t y,Float_t z,Float_t &eta,Float_t &phi){
    // Compute Eta and Phi from x, y, z. Assumes vertex is a zero.
    Float_t r;

    r = x*x+y*y;
    r = TMath::Sqrt(r);
    phi = 180.*TMath::ATan2(y,x)/TMath::Pi();
    if(phi<0.0) phi += 180.;
    r = TMath::ATan2(r,z);
    r = TMath::Tan(0.5*r);
    eta = -TMath::Exp(r);
    return;
}
