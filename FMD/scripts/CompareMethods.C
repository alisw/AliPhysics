#ifndef __CINT__
#include <TH1.h>
#include <TTree.h>
#include <TClonesArray.h>
#include <AliFMD.h>
#include <AliFMDHit.h>
#include <AliFMDGeometry.h>
#include <AliFMDMultStrip.h>
#include <AliFMDMultRegion.h>
#include <AliLoader.h>
#include <AliRunLoader.h>
#include <AliRun.h>
#include <AliStack.h>
#include <AliHeader.h>
#include <AliGenEventHeader.h>
#include <TCanvas.h>
#include <TVector3.h>
#include <TLegend.h>
#include <TClassTable.h>
#include <TParticle.h>
#include <TLegend.h>
#endif
/* Script to compare the Naiive and Poisson method, to the number of
   hits and primaries. */

Double_t Strip2Eta(UShort_t detector, Char_t ring, UShort_t sector, 
		   UShort_t strip, Double_t vz)
{
  //  cout<<" Strip2Eta "<<detector<<" "<<ring<<" "<<strip;
  AliFMDGeometry* fmdgeo = AliFMDGeometry::Instance();
  Double_t x, y, z;
  fmdgeo->Detector2XYZ(detector, ring, sector, strip, x, y, z);
  Double_t realZ = z - vz;
  Double_t r     = TMath::Sqrt(x * x + y * y);
  Double_t theta = TMath::ATan2(r, realZ);
  Double_t eta   = - TMath::Log(TMath::Tan(theta / 2));
  return eta;		       
}


void CompareMethods()
{
  // Comparison of reconstruction Poisson and Naiive methods and real multiplicity

 // Dynamically link some shared libs
#ifdef __CINT__
  if (gClassTable->GetID("AliRun") < 0) {
    gROOT->LoadMacro("$ALICE_ROOT/macros/loadlibs.C");
    loadlibs();
  }
#endif

  Double_t bins1I[15]={3.61728,
		       3.7029,
		       3.8029,
		       3.9029,
		       4.0029,
		       4.1029,
		       4.2029,
		       4.3029,
		       4.4029,
		       4.5029,
		       4.6029,
		       4.7029,
		       4.8029,
		       4.9029,
		       5.0029}; 

  Double_t bins2I[15]={2.28235,
		       2.35884,
		       2.45884,
		       2.55884,
		       2.65884,
		       2.75884,
		       2.85884,
		       2.95884,
		       3.05884,
		       3.15884,
		       3.25884,
		       3.35884,
		       3.45884,
		       3.55884,
		       3.65884};

  Double_t bins3I[15]={-3.37566,
		       -3.27566,
		       -3.17566,
		       -3.07566,
		       -2.97566,
		       -2.87566,
		       -2.77566,
		       -2.67566,
		       -2.57566,
		       -2.47566,
		       -2.37566,
		       -2.27566,
		       -2.17566,
		       -2.07566,
		       -2.00644};

  Double_t bins2O[7]={1.71408,
		      1.77662,
		      1.87662,
		      1.97662,
		      2.07662,
		      2.17662,
		      2.27662};
  Double_t bins3O[7]={-2.27662,
		      -2.17662,
		      -2.07662,
		      -1.97662,
		      -1.87662,
		      -1.77662,
		      -1.71408};
  //__________________________________________________________________
  TH1F* hHits    = new TH1F("hits",    "Hits",    100, -5, 5);
  TH1F* hPrimary = new TH1F("primary", "Primary", 100, -5, 5);
  TH1F* hAll     = new TH1F("all",     "All",     100, -5, 5);
  TH1F* hNaiive  = new TH1F("naiive",  "Naiive",  100, -5, 5);
  TH1F* hPoisson = new TH1F("poisson", "Poisson", 100, -5, 5);
  hNaiive->SetLineStyle(1); hNaiive->SetLineColor(2); 
  hNaiive->SetFillStyle(3004); hNaiive->SetFillColor(2);
  hPoisson->SetLineStyle(1); hPoisson->SetLineColor(3); 
  hPoisson->SetFillStyle(3005); hPoisson->SetFillColor(3);
  hHits->SetLineStyle(1); hHits->SetLineColor(1); 
  hHits->SetFillStyle(3001); hHits->SetFillColor(1);
  hPrimary->SetLineStyle(1); hPrimary->SetLineColor(4); 
  hPrimary->SetFillStyle(3001); hPrimary->SetFillColor(4);
  hAll->SetLineStyle(1); hAll->SetLineColor(6); 
  hAll->SetFillStyle(3001); hAll->SetFillColor(6);
  

  //__________________________________________________________________
  TH1F *h1IHits = new TH1F ("h1IHits", "Hits in FMD1I",
			    14,bins1I[0],bins1I[14]);
  TH1F *h2IHits = new TH1F ("h2IHits", "Hits in FMD2I",
			    14,bins2I[0],bins2I[14]);
  TH1F *h3IHits = new TH1F ("h3IHits", "Hits in FMD3I",
			    14,bins3I[0],bins3I[14]);
  TH1F *h2OHits = new TH1F ("h2OHits", "Hits in FMD2O",
			    6,bins2O[0],bins2O[6]);
  TH1F *h3OHits = new TH1F ("h3OHits", "Hits in FMD3O",
			    6,bins3O[0],bins3O[6]);
  h1IHits->SetLineColor(1); // h1IHits->SetFillColor(1);
  h1IHits->SetLineStyle(1); // h1IHits->SetFillStyle(3000 + 1); 
  h2IHits->SetLineColor(1); // h2IHits->SetFillColor(1);
  h2IHits->SetLineStyle(1); // h2IHits->SetFillStyle(3000 + 1); 
  h2OHits->SetLineColor(1); // h2OHits->SetFillColor(1);
  h2OHits->SetLineStyle(1); // h2OHits->SetFillStyle(3000 + 1); 
  h3IHits->SetLineColor(1); // h3IHits->SetFillColor(1);
  h3IHits->SetLineStyle(1); // h3IHits->SetFillStyle(3000 + 1); 
  h3OHits->SetLineColor(1); // h3OHits->SetFillColor(1);
  h3OHits->SetLineStyle(1); // h3OHits->SetFillStyle(3000 + 1); 
  
  //__________________________________________________________________
  TH1F *h1INaiive = new TH1F ("h1INaiive", "Naiive in FMD1I",
				 14,bins1I[0],bins1I[14]);
  TH1F *h2INaiive = new TH1F ("h2INaiive", "Naiive in FMD2I",
				 14,bins2I[0],bins2I[14]);
  TH1F *h3INaiive = new TH1F ("h3INaiive", "Naiive in FMD3I",
				 14,bins3I[0],bins3I[14]);
  TH1F *h2ONaiive = new TH1F ("h2ONaiive", "Naiive in FMD2O",
				 6,bins2O[0],bins2O[6]);
  TH1F *h3ONaiive = new TH1F ("h3ONaiive", "Naiive in FMD3O",
				 6,bins3O[0],bins3O[6]);
  h1INaiive->SetLineColor(2); // h1INaiive->SetFillColor(2);
  h1INaiive->SetLineStyle(2); // h1INaiive->SetFillStyle(3000 + 2); 
  h2INaiive->SetLineColor(2); // h2INaiive->SetFillColor(2);
  h2INaiive->SetLineStyle(2); // h2INaiive->SetFillStyle(3000 + 2); 
  h2ONaiive->SetLineColor(2); // h2ONaiive->SetFillColor(2);
  h2ONaiive->SetLineStyle(2); // h2ONaiive->SetFillStyle(3000 + 2); 
  h3INaiive->SetLineColor(2); // h3INaiive->SetFillColor(2);
  h3INaiive->SetLineStyle(2); // h3INaiive->SetFillStyle(3000 + 2); 
  h3ONaiive->SetLineColor(2); // h3ONaiive->SetFillColor(2);
  h3ONaiive->SetLineStyle(2); // h3ONaiive->SetFillStyle(3000 + 2); 

  //__________________________________________________________________
  TH1F *h1IPoisson = new TH1F ("h1IPoisson", "Poisson in FMD1I",
				  14,bins1I[0],bins1I[14]);
  TH1F *h2IPoisson = new TH1F ("h2IPoisson", "Poisson in FMD2I",
				  14,bins2I[0],bins2I[14]);
  TH1F *h3IPoisson = new TH1F ("h3IPoisson", "Poisson in FMD3I",
				  14,bins3I[0],bins3I[14]);
  TH1F *h2OPoisson = new TH1F ("h2OPoisson", "Poisson in FMD2O",
				  6,bins2O[0],bins2O[6]);
  TH1F *h3OPoisson = new TH1F ("h3OPoisson", "Poisson in FMD3O",
				  6,bins3O[0],bins3O[6]);
  h1IPoisson->SetLineColor(3); // h1IPoisson->SetFillColor(3);
  h1IPoisson->SetLineStyle(3); // h1IPoisson->SetFillStyle(3000 + 3); 
  h2IPoisson->SetLineColor(3); // h2IPoisson->SetFillColor(3);
  h2IPoisson->SetLineStyle(3); // h2IPoisson->SetFillStyle(3000 + 3); 
  h2OPoisson->SetLineColor(3); // h2OPoisson->SetFillColor(3);
  h2OPoisson->SetLineStyle(3); // h2OPoisson->SetFillStyle(3000 + 3); 
  h3IPoisson->SetLineColor(3); // h3IPoisson->SetFillColor(3);
  h3IPoisson->SetLineStyle(3); // h3IPoisson->SetFillStyle(3000 + 3); 
  h3OPoisson->SetLineColor(3); // h3OPoisson->SetFillColor(3);
  h3OPoisson->SetLineStyle(3); // h3OPoisson->SetFillStyle(3000 + 3); 
  
  //__________________________________________________________________
  TH1F *h1IPrimary = new TH1F ("h1IPrimary", "Primary in FMD1I",
				  14,bins1I[0],bins1I[14]);
  TH1F *h2IPrimary = new TH1F ("h2IPrimary", "Primary in FMD2I",
				  14,bins2I[0],bins2I[14]);
  TH1F *h3IPrimary = new TH1F ("h3IPrimary", "Primary in FMD3I",
				  14,bins3I[0],bins3I[14]);
  TH1F *h2OPrimary = new TH1F ("h2OPrimary", "Primary in FMD2O",
				  6,bins2O[0],bins2O[6]);
  TH1F *h3OPrimary = new TH1F ("h3OPrimary", "Primary in FMD3O",
				  6,bins3O[0],bins3O[6]);
  TH1F *hEdep = new TH1F("hEdep"," edep",100,0,1);
  h1IPrimary->SetLineColor(4); // h1IPrimary->SetFillColor(4);
  h1IPrimary->SetLineStyle(4); // h1IPrimary->SetFillStyle(3000 + 4); 
  h2IPrimary->SetLineColor(4); // h2IPrimary->SetFillColor(4);
  h2IPrimary->SetLineStyle(4); // h2IPrimary->SetFillStyle(3000 + 4); 
  h2OPrimary->SetLineColor(4); // h2OPrimary->SetFillColor(4);
  h2OPrimary->SetLineStyle(4); // h2OPrimary->SetFillStyle(3000 + 4); 
  h3IPrimary->SetLineColor(4); // h3IPrimary->SetFillColor(4);
  h3IPrimary->SetLineStyle(4); // h3IPrimary->SetFillStyle(3000 + 4); 
  h3OPrimary->SetLineColor(4); // h3OPrimary->SetFillColor(4);
  h3OPrimary->SetLineStyle(4); // h3OPrimary->SetFillStyle(3000 + 4); 

  //__________________________________________________________________
  AliRunLoader* runLoader = AliRunLoader::Open("galice.root");
  runLoader->LoadgAlice();
  runLoader->LoadHeader();
  runLoader->LoadKinematics();
  gAlice                   = runLoader->GetAliRun();
  // AliFMD*       fmd     = static_cast<AliFMD*>(gAlice->GetDetector("FMD"));
  AliLoader*    fmdLoader  = runLoader->GetLoader("FMDLoader");
  fmdLoader->LoadHits("READ");
  fmdLoader->LoadRecPoints("READ");
  TArrayF vtx;
  
  Int_t nEvents = runLoader->TreeE()->GetEntries();
  for (Int_t event = 0; event < nEvents && event < 1 ; event++) {
    cout << "Event # " << event << flush;
    runLoader->GetEvent(event);
    AliHeader*         header      = runLoader->GetHeader();
    AliGenEventHeader* eventHeader = header->GenEventHeader();
    eventHeader->PrimaryVertex(vtx);
    Double_t vz = vtx[2];
    cout<<" Vz= "<< vz << flush;

    // Fill primaries 
    Int_t nparticles= runLoader->Stack()->GetNtrack();
    cout << " # particles=" << nparticles << flush;
    for(Int_t ipart=0; ipart<nparticles; ipart++) {
      TParticle *p=runLoader->Stack()->Particle(ipart);
      Float_t etaKine=p->Eta();
      if (p->GetMother(0) == -1) {
	hPrimary->Fill(etaKine);
	h1IPrimary->Fill(etaKine);
	h2IPrimary->Fill(etaKine);
	h3IPrimary->Fill(etaKine);
	h2OPrimary->Fill(etaKine);
	h3OPrimary->Fill(etaKine);
      }
      hAll->Fill(etaKine);
    }
    cout << endl;
    
    // Get the hits and histogram them 
    TClonesArray* hits   = 0;
    TTree*        treeH  = fmdLoader->TreeH();
    TBranch*      branch = treeH->GetBranch("FMD");
    branch->SetAddress(&hits);
    Int_t totalHits   = 0;
    Int_t nHitEntries = treeH->GetEntries();
    cout << "  # hit entries=" << nHitEntries << " " << flush;
    for (Int_t iHitEntry = 0; iHitEntry <  nHitEntries ; iHitEntry++) {
      treeH->GetEntry(iHitEntry);
      
      if (iHitEntry % 1000 == 0) cout << "." << flush;
      Int_t nHits = hits->GetEntries();
      for (Int_t ihit = 0; ihit < nHits; ihit++) {
	AliFMDHit* hit = static_cast<AliFMDHit*>(hits->UncheckedAt(ihit));
	UShort_t detector = hit->Detector();
	Char_t   ring     = hit->Ring();
	UShort_t sector   = hit->Sector();
	UShort_t strip    = hit->Strip();
	Float_t  edep     = hit->Edep();
	if (edep >0 )     hEdep->Fill(edep); 
	Float_t  eta      =  Strip2Eta(detector, ring, sector, strip, vz);
	if (edep > 0.0136) {
	  totalHits++;
	  hHits->Fill(eta);
	  switch (detector) {
	  case 1: h1IHits->Fill(eta); break;
	  case 2: 
	    switch (ring){
	    case 'I': h2IHits->Fill(eta); break;
	    case 'O': h2OHits->Fill(eta); break;
	    }
	    break;
	  case 3: 
	    switch (ring){
	    case 'I': h3IHits->Fill(eta); break;
	    case 'O': h3OHits->Fill(eta); break;
	    }
	    break;
	  }
	} // if (edep > 0.0136)
      } // for (Int_t i = 0; i < nHits; i++)
    } // for (Int_t entry = 0; entry < nEntry; entry++)
    cout << endl;

    // Get the recPoints, and histogram them 
    TClonesArray* multNaiive  = 0;
    TClonesArray* multPoisson = 0;
    TTree*        treeR  = fmdLoader->TreeR();
    TBranch*      branchPoisson = treeR->GetBranch("FMDPoisson");
    TBranch*      branchNaiive  = treeR->GetBranch("FMDNaiive");
    branchPoisson->SetAddress(&multPoisson);
    branchNaiive->SetAddress(&multNaiive);
    
    Int_t nRecEntries  = treeR->GetEntries();
    cout << "  # rec entries=" << nRecEntries << endl;
    for (Int_t iRecEntry = 0; iRecEntry < nRecEntries; iRecEntry++) {
      treeR->GetEntry(iRecEntry);

      // Naiive reconstrution 
      Int_t nNaiives = multNaiive->GetLast();
      cout << "      # naiive=" << nNaiives << " " << flush;
      for (Int_t inaiive = 0; inaiive < nNaiives; inaiive++) {
	if (inaiive % 1000 == 0) cout << "." << flush;
	AliFMDMultStrip* naiive = 
	  static_cast<AliFMDMultStrip*>(multNaiive->UncheckedAt(inaiive));
	Float_t  nParticles = naiive->Particles();
	Float_t  eta        = naiive->Eta();
	Char_t   ring       = naiive->Ring();
	UShort_t detector   = naiive->Detector();
	hNaiive->Fill(eta, nParticles);
	switch (detector) {
	case 1: h1INaiive->Fill(eta,nParticles); break;
	case 2: 
	  switch (ring){
	  case 'i':
	  case 'I': h2INaiive->Fill(eta,nParticles); break;
	  case 'o':
	  case 'O': h2ONaiive->Fill(eta,nParticles); break;
	  }
	  break;
	case 3: 
	  switch (ring){
	  case 'i':
	  case 'I': h3INaiive->Fill(eta,nParticles); break;
	  case 'o':
	  case 'O': h3ONaiive->Fill(eta,nParticles); break;
	  }
	  break;
	}
      } // for (inrec = 0; ...)
      cout << endl;

      // Poisson reconstruction 
      Int_t nPoissons = multPoisson->GetLast();
      cout << "      # poisson=" << nPoissons << " " << flush;
      for (Int_t ipoisson = 0; ipoisson <= nPoissons; ipoisson++) {
	cout << "." << flush;
	AliFMDMultRegion* poisson = 
	  static_cast<AliFMDMultRegion*>(multPoisson->UncheckedAt(ipoisson));

	Float_t    nParticles = poisson->Particles();
	UShort_t   detector   = poisson->Detector();
	Char_t     ring       = poisson->Ring();
	Float_t    eta        = (poisson->MaxEta() + poisson->MinEta()) / 2;
	// poisson->Print("EPTD");
	// cout << "  Eta returned: " << eta << endl;
	hPoisson->Fill(eta, nParticles);
	switch (detector) {
	case 1: h1IPoisson->Fill(eta,nParticles); break;
	case 2: 
	  switch (ring){
	  case 'i':
	  case 'I': h2IPoisson->Fill(eta,nParticles); break;
	  case 'o':
	  case 'O': h2OPoisson->Fill(eta,nParticles); break;
	  }
	  break;
	case 3: 
	  switch (ring){
	  case 'i':
	  case 'I': h3IPoisson->Fill(eta,nParticles); break;
	  case 'o':
	  case 'O': h3OPoisson->Fill(eta,nParticles); break;
	  }
	  break;
	}
      } // for (ipoisson = 0; ...)
      cout << endl;
    }
  }
  cout << "All done, drawing ... "  << endl;
  
  TCanvas* nullc = new TCanvas("nullc", "Digit Data");
  TH1* null1 = new TH1D("null1", "null", 10, 0, 1);
  null1->Draw();
  
  TCanvas* c2 = new TCanvas("hit", "Hit Data");
  c2->cd();
  hAll->Draw();
  hPrimary->Draw("same");
  hHits->Draw("same");
  hNaiive->Draw("same");
  hPoisson->Draw("same");
  TLegend* l = new TLegend(.7, .7, 1, 1);
  l->AddEntry(hAll,  "All", "L");
  l->AddEntry(hPrimary, "Primary", "L");
  l->AddEntry(hHits,    "Hits", "L");
  l->AddEntry(hNaiive,  "Naiive", "L");
  l->AddEntry(hPoisson, "Poisson", "L");
  l->Draw();
  

  TCanvas* c1 = new TCanvas("perDetector", "Per detector", 1200, 800);
  c1->Divide(3,2);
  c1->cd(6);
  TH1* null = new TH1D("null", "null", 10, 0, 1);
  null->Draw("A");
  TLegend* l2 = new TLegend(0, 0, 1, 1);
  l2->AddEntry(h1INaiive,  "Naiive", "L");
  l2->AddEntry(h1IPoisson, "Poission", "L");
  l2->AddEntry(h1IHits,    "Hits", "L");
  l2->AddEntry(h1IPrimary, "Primary", "L");
  l2->Draw();
    
  c1->cd(3);
  h1IPrimary->SetMinimum(0); h1IPrimary->SetMaximum(300);
  h1IPrimary->Draw();;
  h1IHits->Draw("same");
  h1IPoisson->Draw("same");
  h1INaiive->Draw("same");

  c1->cd(5);
  h2OPrimary->SetMinimum(0); h2OPrimary->SetMaximum(300);
  h2OPrimary->Draw();;
  h2OHits->Draw("same");
  h2OPoisson->Draw("same");
  h2ONaiive->Draw("same");
  
  c1->cd(2);
  h2IPrimary->SetMinimum(0); h2IPrimary->SetMaximum(300);
  h2IPrimary->Draw();;
  h2IHits->Draw("same");
  h2IPoisson->Draw("same");
  h2INaiive->Draw("same");


  c1->cd(4);
  h3OPrimary->SetMinimum(0); h3OPrimary->SetMaximum(300);
  h3OPrimary->Draw();;
  h3OHits->Draw("same");
  h3ONaiive->Draw("same");
  h3OPoisson->Draw("same");

  c1->cd(1);
  h3IPrimary->SetMinimum(0); h3IPrimary->SetMaximum(300);
  h3IPrimary->Draw();;
  h3IHits->Draw("same");
  h3INaiive->Draw("same");
  h3IPoisson->Draw("same");


  
  TFile *fileHist = new TFile("compare2000.root","RECREATE");
  hHits->Write();
  hPrimary->Write();
  hNaiive->Write();
  hPoisson->Write();
  hAll->Write();
  
  h1INaiive->Write();
  h2INaiive->Write();
  h3INaiive->Write();
  h2ONaiive->Write();
  h3ONaiive->Write();
  h1IPoisson->Write();
  h2IPoisson->Write();
  h3IPoisson->Write();
  h2OPoisson->Write();
  h3OPoisson->Write();
  h1IHits->Write();
  h2IHits->Write();
  h3IHits->Write();
  h2OHits->Write();
  h3OHits->Write();
  h1IPrimary->Write();
  h2IPrimary->Write();
  h3IPrimary->Write();
  h2OPrimary->Write();
  h3OPrimary->Write();
  hEdep->Write();
  nullc->Close();

}
