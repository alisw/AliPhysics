/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
 
//---------------------------------------------------------------------
// JetAnalysis class 
// Analyse Jets
// Author: andreas.morsch@cern.ch
//---------------------------------------------------------------------
 
#include "AliJetAnalysis.h"
ClassImp(AliJetAnalysis)
 
 
////////////////////////////////////////////////////////////////////////
#include <TH1F.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TLorentzVector.h>

#include "AliJetProductionDataPDC2004.h"
#include "AliJet.h"
#include "AliJetESDReaderHeader.h"
#include "AliUA1JetHeader.h"
#include "AliLeading.h"

    AliJetAnalysis::AliJetAnalysis()
{
    // Constructor
    fDirectory    = 0x0;   
    fEventMin     =  0;
    fEventMax     = -1;
    fRunMin       =  0;
    fRunMax       = 11;
}

void AliJetAnalysis::Analyse() 
{
    //
    // Some histos
    //
    TH1F::AddDirectory(0);
    TProfile::AddDirectory(0);
    
    TH1F* e0H    = new TH1F("e0H"  ,"Jet Energy (reconstructed)",      40,  0., 200.);
    TH1F* e1H    = new TH1F("e1H" , "Jet Energy (generated)",          40,  0., 200.);
    TH1F* e2H    = new TH1F("e2H" , "Jet Energy (generated, nrec = 0", 40,  0., 200.);
    TH1F* e3H    = new TH1F("e3H" , "Jet Energy (leading)",            40,  0., 200.);
    TH1F* e4H    = new TH1F("e4H" , "Jet Energy (reconstructed: 105 < Egen < 125", 40,  0., 200.);

    TH1F* e5H    = new TH1F("e5H" , "Jet Energy (generated)", 40,  0., 200.);
    TH1F* e6H    = new TH1F("e6H" , "Jet Energy (generated)", 40,  0., 200.);
    TH1F* e7H    = new TH1F("e7H" , "Jet Energy (generated)", 40,  0., 200.);
    TH1F* e8H    = new TH1F("e8H" , "Jet Energy (generated)", 40,  0., 200.);

    TProfile* r5H    = new TProfile("r5H" , "rec/generated", 20,  0., 200, 0., 1., "S");
    TProfile* r6H    = new TProfile("r6H" , "rec/generated", 20,  0., 200, 0., 1., "S");

    TProfile* r7H    = new TProfile("r7H" , "rec/generated", 20,  0., 200, 0., 1., "S");
    TProfile* r8H    = new TProfile("r8H" , "rec/generated", 20,  0., 200, 0., 1., "S");


    TH1F* dr1H = new TH1F("dr1H", "delta R",  160, 0.,   2.);
    TH1F* dr2H = new TH1F("dr2H", "delta R",  160, 0.,   2.);
    TH1F* dr3H = new TH1F("dr4H", "delta R",  160, 0.,   2.);

    TH1F* etaH  = new TH1F("etaH",  "eta",  160, -2.,   2.);
    TH1F* eta1H = new TH1F("eta1H", "eta",  160, -2.,   2.);
    TH1F* eta2H = new TH1F("eta2H", "eta",  160, -2.,   2.);

    TH1F* phiH   = new TH1F("phiH",   "phi",  160, -3.,   3.);
    TH1F* dphiH  = new TH1F("dphiH",  "phi",  160,  0.,   3.1415);
    TH1F* phi1H  = new TH1F("phi1H", "phi",   160,  0.,   6.28);
    TH1F* phi2H  = new TH1F("phi2H", "phi",   160,  0.,   6.28);
    

    TProfile* drP1    = new TProfile("drP1" , "Delta_R", 20,  0., 200, -1., 1., "S");
    TProfile* drP2    = new TProfile("drP2" , "Delta_R", 20,  0., 200, -1., 1., "S");

  // Run data 
    AliJetProductionDataPDC2004* runData = new AliJetProductionDataPDC2004();
    

    // Loop over runs
    TFile* jFile = 0x0;
  
    for (Int_t iRun = fRunMin; iRun <= fRunMax; iRun++)
    {

	// Open file
	char fn[20];
	sprintf(fn, "%s/%s.root", fDirectory, (runData->GetRunTitle(iRun)).Data());
	
	
	jFile = new TFile(fn);

	printf("  Analyzing run: %d %s\n", iRun,fn);	
	// Get jet header and display parameters
	AliUA1JetHeader* jHeader = 
	    (AliUA1JetHeader*) (jFile->Get("AliUA1JetHeader"));
	// jHeader->PrintParameters();
	
	// Get reader header and events to be looped over
	AliJetESDReaderHeader *jReaderH = 
	    (AliJetESDReaderHeader*)(jFile->Get("AliJetKineReaderHeader"));

	if (fEventMin == -1) fEventMin =  jReaderH->GetFirstEvent();
	if (fEventMax == -1) {
	    fEventMax =  jReaderH->GetLastEvent();
	} else {
	    fEventMax = TMath::Min(fEventMax, jReaderH->GetLastEvent());
	}
	
	
	// Calculate weight
	Float_t wgt = runData->GetWeight(iRun) / Float_t(fEventMax - fEventMin + 1);
	Float_t ptmin, ptmax;
	runData->GetPtHardLimits(iRun, ptmin, ptmax);
	
	
	// Loop over events
	AliJet *jets  = 0x0; 
	AliJet *gjets = 0x0;
	AliLeading *leading = 0x0;
	Float_t egen, etag, econ, erec;
	
	
	for (Int_t i = fEventMin; i < fEventMax; i++) {
	    printf("  Analyzing run: %d  Event %d / %d \n", 
		   iRun, i, fEventMax);
	    
	    // Het next tree with AliJet
	    char nameT[100];
	    sprintf(nameT, "TreeJ%d",i);
	    TTree *jetT =(TTree *)(jFile->Get(nameT));
	    jetT->SetBranchAddress("FoundJet",    &jets);
	    jetT->SetBranchAddress("GenJet",      &gjets);
	    jetT->SetBranchAddress("LeadingPart", &leading);
	    jetT->GetEntry(0);
	    
//
//    Find the jet with the highest E_T within fiducial region
//
	    Int_t njets = jets->GetNJets();
	    Int_t imax = -1;
	    Int_t jmax = -1;
	    Float_t emax = 0.;
	    
	    for (Int_t ij = 0; ij < njets; ij++) {
		if (jets->GetPt(ij) > emax && 
		    TMath::Abs(jets->GetEta(ij)) < 0.60) {
		    emax = jets->GetPt(ij);
		    jmax = imax;
		    imax = ij;
		}
	    }
	    
	    
	    if (imax == -1) {
		Int_t ngen = gjets->GetNJets();
		if(ngen>0) e2H->Fill(gjets->GetPt(0), wgt);
	    } else {
//	      printf("Reconstructed Jet %5d %13.3f %13.3f %13.3f\n", 
//                    imax, emax, jets->GetEta(imax), jets->GetPhi(imax));
		econ = jets->GetPt(imax);
		erec = jets->GetPt(imax) / 0.65;
		dr1H->Fill(jets->GetEta(imax) - gjets->GetEta(0));
//
//    Find the generated jet closest to the reconstructed 
//
		Float_t rmin;
		Float_t etamin = 1.e6;
		
		Int_t   igen;
		Float_t etaj = jets->GetEta(imax);
		Float_t phij = jets->GetPhi(imax);
		
		Int_t ngen = gjets->GetNJets();
		
		if (ngen != 0) {
		    rmin = 1.e6;
		    igen = -1;
		    for (Int_t j = 0; j < ngen; j++) {
			etag = gjets->GetEta(j);
			Float_t phig = gjets->GetPhi(j);
			Float_t deta = etag - etaj;
			Float_t dphi = TMath::Abs(phig - phij);
			if (dphi > TMath::Pi()) dphi = 2. * TMath::Pi() - dphi;
			Float_t r = TMath::Sqrt(deta * deta + dphi * dphi);
			if (r  < rmin) {
			    rmin = r;
			    etamin = deta;
			    igen = j;
			}
		    }
		    
		    egen = gjets->GetPt(igen);
		    e1H->Fill(gjets->GetPt(igen), wgt);
		    etag = gjets->GetEta(igen);
		    Float_t phig = gjets->GetPhi(igen);
		    Float_t dphi = phig - phij;
		    
//		  if (econ < ptmax) {
		    e0H->Fill(erec, wgt);
//		  } else {
//		      e0H->Fill(erec, 6.7e-6);
//		  }
		    
		    
		    
		    if (egen > 20. && egen < 40.) {
			phiH->Fill((dphi));
			etaH->Fill(etag - etaj, wgt);
			phi1H->Fill(phig);
			phi2H->Fill(phij);		  
			e4H->Fill(jets->GetPt(imax));
		    }
		    
		    if (erec > 90. && erec < 110. && rmin < 0.1) {
			e5H->Fill(egen, wgt);
			dr2H->Fill(rmin);
			if (egen < 30.) {
			    printf("Strange jet %6d %13.3f %13.3f %13.3f \n", 
				   imax, etaj, phij, erec);
			    for (Int_t j = 0; j < ngen; j++) {
				printf("Generated %6d %13.3f %13.3f %13.3f\n", 
				       j, gjets->GetEta(j), 
				       gjets->GetPhi(j), gjets->GetPt(j));
				
			    }
			}
		    }
		    
		    if (rmin < 0.1) {
			r5H->Fill(egen, jets->GetPt(imax)/egen, wgt);
			r6H->Fill(jets->GetPt(imax) 
				  / 0.4, jets->GetPt(imax)/egen, wgt);
			e8H->Fill(erec, wgt);
		    }
		    
		    if (rmin < 0.1) {
			drP1->Fill(egen, etamin, wgt);
		    }
		} // ngen !=0
	    } // has reconstructed jet
	    
//
//   Leading particle
//
	    if (leading->LeadingFound()) {
		Float_t etal = leading->GetLeading()->Eta();
		Float_t phil = leading->GetLeading()->Phi();
		Float_t el   = leading->GetLeading()->E();
//	  printf("Leading %f %f %f \n", etal, phil, el);
	  
		e3H->Fill(el, wgt);
	  
		Float_t rmin = 1.e6;
		Float_t etamin = 1.e6;
	  
		Int_t igen = -1;
		Int_t ngen = gjets->GetNJets();
		for (Int_t j = 0; j < ngen; j++) {
		    etag = gjets->GetEta(j);
		    Float_t phig = gjets->GetPhi(j);
		    Float_t deta = etag-etal;
		    Float_t dphi = TMath::Abs(phig - phil);
		    if (dphi > TMath::Pi()) dphi = 2. * TMath::Pi() - dphi;
		    
		    Float_t r = TMath::Sqrt(deta * deta + dphi * dphi);
	      
		    if (r  < rmin) {
			rmin = r;
			igen = j;
			etamin = deta;
		    }
		}
		if (egen > 20. && egen < 40.) {
		    dr3H->Fill(rmin);
		    eta1H->Fill(etag-etal, wgt);
		    
		}
	  
		if (el > 54. && el  < 66.) e6H->Fill(egen, wgt);
		if (rmin < 0.3) {
		    r7H->Fill(egen, el/egen, wgt);
		    r8H->Fill(el / 0.2, el/egen, wgt);
		}
		
		if (rmin < 0.1) {
		    drP2->Fill(egen, etamin, wgt);
		}
	    } // if leading particle
	    


//
// Generated Jet
//
	    Int_t ngen = gjets->GetNJets();
	    emax = 0.;
	    imax = -1;
	    
	    for (Int_t j = 0; j < ngen; j++) {
		if (gjets->GetPt(j) > 
		    emax && TMath::Abs(gjets->GetEta(j)) < 0.5) {
		    emax = gjets->GetPt(j);
		    imax = j;
		}
	    }
	    if (imax != -1) e7H->Fill(emax, wgt);
	    
	    
	    delete jetT;
	    
	} // events
	if (jFile) jFile->Close();
	delete jFile;
	
    } // runs
    
//  delete jFile;
//  if (jFile) jFile->Close();  
/*
  TFile* f = new TFile("j.root", "recreate");
  e0H->Write();
  e1H->Write();
  e2H->Write();
  e3H->Write();
  e4H->Write();
  e7H->Write();
  e8H->Write();
  f->Close();
*/

    // Get Control Plots
//  gStyle->SetOptStat(0);
    
    TCanvas* c1 = new TCanvas("c1");
  e0H->Draw();
  e1H->SetLineColor(2);
  e2H->SetLineColor(4);
  e3H->SetLineColor(5);
  e1H->Draw("same");
  e3H->Draw("same");


  TCanvas* c2 = new TCanvas("c2");
//  dr1H->Draw();
  dr2H->SetLineColor(2);
  dr2H->Draw("");
  
  TCanvas* c3 = new TCanvas("c3");
  dr2H->Draw();
  dr3H->Draw("same");

  TCanvas* c4 = new TCanvas("c4");
  e0H->Draw();

  TCanvas* c5 = new TCanvas("c5");
  etaH->SetXTitle("#eta_{rec} - #eta_{gen}");
  
  etaH->Draw();
  eta1H->SetLineColor(2);
  
  eta1H->Draw("same");
  
  TCanvas* c5a = new TCanvas("c5a");
  eta1H->Draw();

  TCanvas* c5b = new TCanvas("c5b");
  eta2H->Draw();

  TCanvas* c6 = new TCanvas("c6");
  e4H->Draw();
  TCanvas* c7 = new TCanvas("c7");
  phiH->Draw();

  TCanvas* c7a = new TCanvas("c7a");
  phi1H->Draw();
  TCanvas* c7b = new TCanvas("c7b");
  phi2H->Draw();

  TCanvas* c8 = new TCanvas("c8");
  e5H->SetXTitle("E_{gen} (GeV)");
  e5H->Draw();
  e6H->SetLineColor(2);
  e6H->Draw("same");
  
  TCanvas* c9 = new TCanvas("c9");
  e6H->Draw();

  
  
  TCanvas* c10 = new TCanvas("c10");

  
  r5H->SetMaximum(1.);
  r5H->Draw();
  r5H->SetXTitle("E_{gen} (GeV)");
  r5H->SetYTitle("E_{leading} / E_{gen}");
  r6H->SetLineColor(2);
  r6H->Draw("same");

  TCanvas* c11 = new TCanvas("c11");
//
// Leading particle rec/gen
//
  r7H->SetMaximum(1.);
  
  r7H->Draw();
  r7H->SetXTitle("E_{gen} (GeV)");
  r7H->SetYTitle("E_{leading} / E_{gen}");
  
  r8H->SetLineColor(2);
  r8H->Draw("same");
  
  TCanvas* c12 = new TCanvas("c12");
  drP1->SetXTitle("E_{gen} (GeV)");
  drP1->SetYTitle("#eta_{rec} - #eta_{gen}");
  drP1->Draw();

  TCanvas* c13 = new TCanvas("c13");
  drP2->SetXTitle("E_{gen} (GeV)");
  drP2->SetYTitle("#eta_{leading} - #eta_{gen}");
  
  drP2->Draw();
  
  TCanvas* c14 = new TCanvas("c14");
  dphiH->Draw();

/*
  e1H->Write();
  e2H->Write();
  e3H->Write();
  e4H->Write();
*/
} 

