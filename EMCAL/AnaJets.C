void AnaJets(Int_t evNumber1=0, Int_t evNumber2=0) 
{
//*-- Author: Andreas Morsch (CERN)

    if (gClassTable->GetID("AliRun") < 0) {
	gROOT->LoadMacro("../macros/loadlibs.C");
	loadlibs();
    }
// Connect the Root Galice file containing Geometry, Kine and Hits

    TFile *file = (TFile*)gROOT->GetListOfFiles()->FindObject("galice.root");
    
    if (!file) {
	printf("\n Creating galice.root \n");
	file = new TFile("galice.root");
    } else {
	printf("\n galice.root found in file list");
    }
// Get AliRun object from file or create it if not on file
    if (!gAlice) {
	gAlice = (AliRun*)(file->Get("gAlice"));
	if (gAlice) printf("AliRun object found on file\n");
	if (!gAlice) {
	    printf("\n create new gAlice object");
	    gAlice = new AliRun("gAlice","Alice test program");
	}
    }
// Book histos    
    TH1F *eH   = new TH1F("eH","Energy", 150,  0.0, 150.);
    TH1F *etaH = new TH1F("eEta","Eta",  180, -0.9, 0.9);
    TH1F *phiH = new TH1F("ePhi","Phi",   62, -3.1, 3.1);
 

    TClonesArray* jets = new TClonesArray("AliEMCALJet",10000);

    for (int nev=0; nev<= evNumber2; nev++) {
	printf("\n Event .............%d", nev);
	Int_t nparticles = gAlice->GetEvent(nev);
	Int_t nbytes     = 0;
	AliEMCAL *pEMCAL  = (AliEMCAL*) gAlice->GetModule("EMCAL");
	if (pEMCAL) {
	    TTree *TR = gAlice->TreeR();
	    Int_t nent=TR->GetEntries();
	    TR->SetBranchAddress("Jets", &jets);
	    nbytes += TR->GetEntry(0);
	    Int_t nJet = jets->GetEntries();
	    printf("\n Number of Jets %d", nJet);
	    AliEMCALJet  *mJet;
	    for (Int_t ij=0; ij < nJet; ij++) {
		mJet = (AliEMCALJet*)jets->UncheckedAt(ij);
		printf("\n Jet:%d E %f phi %f eta %f\n", ij, 
		       mJet->Energy(), mJet->Phi(), mJet->Eta());
		etaH->Fill(mJet->Eta());
		phiH->Fill(mJet->Phi());
		eH  ->Fill(mJet->Energy());		
	   } // jet
       } // ?EMCAL
   } // event
    TCanvas *c1 = new TCanvas("c1","Canvas 1",400,10,600,700);
    c1->Divide(2,2);
    c1->cd(1); eH->Draw();
    c1->cd(2); etaH->Draw();
    c1->cd(3); phiH->Draw();
}
