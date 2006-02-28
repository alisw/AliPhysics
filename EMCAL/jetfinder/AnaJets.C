void AnaJets(Int_t evNumber1=0, Int_t evNumber2=0) 
{
//*-- Author: Andreas Morsch (CERN)

    if (gClassTable->GetID("AliRun") < 0) {
	gROOT->LoadMacro("loadlibs.C");
	loadlibs();
    }
// Connect the Root Galice file containing Geometry, Kine and Hits
    
    TFile *file = (TFile*)gROOT->GetListOfFiles()->FindObject("jets.root");
    TFile *source = (TFile*)gROOT->GetListOfFiles()->FindObject("galice.root");

    if (!file) {
	printf("\n Creating galice.root \n");
	file = new TFile("jets.root");
    } else {
	printf("\n galice.root found in file list");
    }

 if (!source) {
	printf("\n Creating galice.root \n");
	source = new TFile("galice.root");
    } else {
	printf("\n galice.root found in file list");
    }
// Get AliRun object from file or create it if not on file
 //   if (!gAlice) {
	gAlice = (AliRun*)(source->Get("gAlice"));
	if (gAlice) printf("AliRun object found on file\n");
	if (!gAlice) {
	    printf("\n create new gAlice object");
	    gAlice = new AliRun("gAlice","Alice test program");
	}
	// }
// Book histos    
    TH1F *eH   = new TH1F("eH","Energy",    200,  0.0, 200.);
    TH1F *etaH = new TH1F("eEta","Eta",     180, -0.9,  0.9);
    TH1F *phiH = new TH1F("ePhi","Phi",      62, -3.1,  3.1);
    TH1F *tH   = new TH1F("tH","n tracks",   30,  0.5, 29.5);
    TH1F *ptH  = new TH1F("ptH","Track pT", 100., 0., 100.);
    TH1F *drH  = new TH1F("drH","Track dR", 120., 0.,   6.);    
    
    Float_t phiT[50], etaT[50], ptT[50];
    

    TClonesArray* jets = new TClonesArray("AliEMCALJet",10000);
    
    for (int nev=0; nev<= evNumber2; nev++) {
	printf("\n Event .............%d", nev);
	Int_t nparticles = gAlice->GetEvent(nev);
	Int_t nbytes     = 0;
	AliEMCAL *pEMCAL  = (AliEMCAL*) gAlice->GetModule("EMCAL");
	if (pEMCAL) {
	  TTree *TR =(TTree *)(file->Get("TreeR0"));
	 
	    Int_t nent=TR->GetEntries();
	    TR->SetBranchAddress("EMCALJets", &jets);
	    nbytes += TR->GetEntry(0);
	    Int_t nJet = jets->GetEntries();
	    printf("\n Number of Jets %d", nJet);
	    AliEMCALJet  *mJet;
	    for (Int_t ij=0; ij < nJet; ij++) {
		mJet = (AliEMCALJet*)jets->UncheckedAt(ij);
		Float_t eta = mJet->Eta();
		
		printf("\n Jet:%d E %f phi %f eta %f tracks %d\n", ij, 
		       mJet->Energy(), mJet->Phi(), eta,
		       mJet->NTracks());
		etaH->Fill(mJet->Eta());
		phiH->Fill(mJet->Phi());
		if (TMath::Abs(eta) < 0.4) 
		    eH->Fill(mJet->Energy());
		tH  ->Fill((Float_t)mJet->NTracks());
		
		mJet->TrackList(ptT, etaT, phiT);
		for (Int_t it = 0; it < mJet->NTracks(); it++)
		{
		    printf(" Track: %5d pT %8.3f eta %8.3f phi %8.3f \n",
			   it, ptT[it], etaT[it], phiT[it]);
		    ptH->Fill(ptT[it]);
		    Float_t dPhi = phiT[it]-mJet->Phi();
		    Float_t dEta = etaT[it]-mJet->Eta();
		    Float_t dr = TMath::Sqrt(dPhi*dPhi+dEta*dEta);
		    drH->Fill(dr);
		}
	   } // jet
       } // ?EMCAL
   } // event
    TCanvas *c1 = new TCanvas("c1","Canvas 1",400,10,600,700);
    c1->Divide(2,2);
    c1->cd(1); eH->Draw();
    c1->cd(2); etaH->Draw();
    c1->cd(3); phiH->Draw();
    c1->cd(4); tH->Draw();

    TCanvas *c2 = new TCanvas("c2","Canvas 2",400,10,600,700);
    c2->Divide(2,2);
    c2->cd(1); ptH->Draw();
    c2->cd(2); drH->Draw();
}

