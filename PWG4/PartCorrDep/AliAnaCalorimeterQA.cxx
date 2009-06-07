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
/* $Id: $ */

//_________________________________________________________________________
// Class to check results from simulations or reconstructed real data. 
// Fill few histograms and do some checking plots
//
//-- Author: Gustavo Conesa (INFN-LNF)
//_________________________________________________________________________


// --- ROOT system ---
#include "Riostream.h"
#include "TRefArray.h"
#include "TParticle.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TROOT.h"
//#include "TH3F.h"
#include "TH2F.h"
#include "TLegend.h"

//---- AliRoot system ----
#include "AliAnaCalorimeterQA.h"
#include "AliCaloTrackReader.h"
#include "AliStack.h"
#include "AliAODCaloCells.h"
#include "AliAODCaloCluster.h"
#include "AliFidutialCut.h"
#include "AliESDtrack.h"
#include "AliESDCaloCluster.h"
#include "AliESDEvent.h"

ClassImp(AliAnaCalorimeterQA)
  
//____________________________________________________________________________
  AliAnaCalorimeterQA::AliAnaCalorimeterQA() : 
    AliAnaPartCorrBaseClass(), fCalorimeter(""), fStyleMacro(""), fhE(0),fhPt(0),fhPhi(0),fhEta(0), fhEtaPhi(0), 
    fhECharged(0),fhPtCharged(0),fhPhiCharged(0),fhEtaCharged(0), fhEtaPhiCharged(0), 
    fhEChargedNoOut(0),fhPtChargedNoOut(0),fhPhiChargedNoOut(0),fhEtaChargedNoOut(0), fhEtaPhiChargedNoOut(0), 
    fhDeltaE(0), fhDeltaPt(0),fhDeltaPhi(0),fhDeltaEta(0), fhRatioE(0), fhRatioPt(0),fhRatioPhi(0),fhRatioEta(0),
    fh2E(0),fh2Pt(0),fh2Phi(0),fh2Eta(0), fhIM(0), fhAsym(0), fhNCellsPerCluster(0), fhNClusters(0), fhNCells(0), fhAmplitude(0), 
    fhGenGamPt(0),fhGenGamEta(0),fhGenGamPhi(0),fhGenPi0Pt(0),fhGenPi0Eta(0),fhGenPi0Phi(0),
    fhGenEtaPt(0),fhGenEtaEta(0),fhGenEtaPhi(0),fhGenOmegaPt(0),fhGenOmegaEta(0),fhGenOmegaPhi(0),
    fhGenElePt(0),fhGenEleEta(0),fhGenElePhi(0), fhEMVxyz(0),  fhEMR(0), fhHaVxyz(0),  fhHaR(0),
    fhGamE(0),fhGamPt(0),fhGamPhi(0),fhGamEta(0), 
    fhGamDeltaE(0), fhGamDeltaPt(0),fhGamDeltaPhi(0),fhGamDeltaEta(0), fhGamRatioE(0), fhGamRatioPt(0),fhGamRatioPhi(0),fhGamRatioEta(0),
    fhGam2E(0),fhGam2Pt(0),fhGam2Phi(0),fhGam2Eta(0),
    fhEleE(0),fhElePt(0),fhElePhi(0),fhEleEta(0),
    fhPi0E(0),fhPi0Pt(0),fhPi0Phi(0),fhPi0Eta(0), fhNeHadE(0),fhNeHadPt(0),fhNeHadPhi(0),fhNeHadEta(0), 
    fhChHadE(0),fhChHadPt(0),fhChHadPhi(0),fhChHadEta(0),
    fhGenGamAccE(0),fhGenGamAccPt(0),fhGenGamAccEta(0),fhGenGamAccPhi(0),
    fhGenPi0AccE(0),fhGenPi0AccPt(0),fhGenPi0AccEta(0),fhGenPi0AccPhi(0),
    fh1pOverE(0),fh1dR(0),fh2EledEdx(0), fh2MatchdEdx(0)
{
  //Default Ctor

  //Initialize parameters
  InitParameters();
}

//____________________________________________________________________________
AliAnaCalorimeterQA::AliAnaCalorimeterQA(const AliAnaCalorimeterQA & qa) :   
  AliAnaPartCorrBaseClass(qa), fCalorimeter(qa.fCalorimeter), fStyleMacro(qa.fStyleMacro),
  fhE(qa.fhE),fhPt(qa.fhPt), fhPhi(qa.fhPhi), fhEta(qa.fhEta), 
  fhEtaPhi(qa.fhEtaPhi), fhECharged(qa.fhECharged),fhPtCharged(qa.fhPtCharged),fhPhiCharged(qa.fhPhiCharged),
  fhEtaCharged(qa.fhEtaCharged), fhEtaPhiCharged(qa.fhEtaPhiCharged), 
  fhEChargedNoOut(qa.fhEChargedNoOut),fhPtChargedNoOut(qa.fhPtChargedNoOut),fhPhiChargedNoOut(qa.fhPhiChargedNoOut),
  fhEtaChargedNoOut(qa.fhEtaChargedNoOut), fhEtaPhiChargedNoOut(qa.fhEtaPhiChargedNoOut), 
  fhDeltaE(qa.fhDeltaE), fhDeltaPt(qa.fhDeltaPt), fhDeltaPhi(qa.fhDeltaPhi), fhDeltaEta(qa.fhDeltaEta),
  fhRatioE(qa.fhRatioE), fhRatioPt(qa.fhRatioPt), fhRatioPhi(qa.fhRatioPhi), fhRatioEta(qa.fhRatioEta),
  fh2E(qa.fh2E), fh2Pt(qa.fh2Pt), fh2Phi(qa.fh2Phi),fh2Eta(qa.fh2Eta), 
  fhIM(qa.fhIM), fhAsym(qa.fhAsym), fhNCellsPerCluster(qa.fhNCellsPerCluster), fhNClusters(qa.fhNClusters), 
  fhNCells(qa.fhNCells), fhAmplitude(qa.fhAmplitude),
  fhGenGamPt(qa.fhGenGamPt), fhGenGamEta(qa.fhGenGamEta), fhGenGamPhi(qa.fhGenGamPhi),
  fhGenPi0Pt(qa.fhGenPi0Pt), fhGenPi0Eta(qa.fhGenPi0Eta), fhGenPi0Phi(qa.fhGenPi0Phi),
  fhGenEtaPt(qa.fhGenEtaPt), fhGenEtaEta(qa.fhGenEtaEta), fhGenEtaPhi(qa.fhGenEtaPhi),
  fhGenOmegaPt(qa.fhGenOmegaPt), fhGenOmegaEta(qa.fhGenOmegaEta), fhGenOmegaPhi(qa.fhGenOmegaPhi),
  fhGenElePt(qa.fhGenElePt), fhGenEleEta(qa.fhGenEleEta), fhGenElePhi(qa.fhGenElePhi), 
  fhEMVxyz(qa.fhEMVxyz),  fhEMR(qa.fhEMR), fhHaVxyz(qa.fhHaVxyz),  fhHaR(qa.fhHaR),
  fhGamE(qa.fhGamE),fhGamPt(qa.fhGamPt),fhGamPhi(qa.fhGamPhi),fhGamEta(qa.fhGamEta), 
  fhGamDeltaE(qa.fhGamDeltaE), fhGamDeltaPt(qa.fhGamDeltaPt), fhGamDeltaPhi(qa.fhGamDeltaPhi), fhGamDeltaEta(qa.fhGamDeltaEta),
  fhGamRatioE(qa.fhGamRatioE), fhGamRatioPt(qa.fhGamRatioPt), fhGamRatioPhi(qa.fhGamRatioPhi), fhGamRatioEta(qa.fhGamRatioEta),
  fhGam2E(qa.fhGam2E), fhGam2Pt(qa.fhGam2Pt), fhGam2Phi(qa.fhGam2Phi),fhGam2Eta(qa.fhGam2Eta), 
  fhEleE(qa.fhEleE),fhElePt(qa.fhElePt),fhElePhi(qa.fhElePhi),fhEleEta(qa.fhEleEta),
  fhPi0E(qa.fhPi0E),fhPi0Pt(qa.fhPi0Pt),fhPi0Phi(qa.fhPi0Phi),fhPi0Eta(qa.fhPi0Eta), 
  fhNeHadE(qa.fhNeHadE),fhNeHadPt(qa.fhNeHadPt),fhNeHadPhi(qa.fhNeHadPhi),fhNeHadEta(qa.fhNeHadEta), 
  fhChHadE(qa.fhChHadE),fhChHadPt(qa.fhChHadPt),fhChHadPhi(qa.fhChHadPhi),fhChHadEta(qa.fhChHadEta),
  fhGenGamAccE(qa.fhGenGamAccE),fhGenGamAccPt(qa.fhGenGamAccPt),fhGenGamAccEta(qa.fhGenGamAccEta),fhGenGamAccPhi(qa.fhGenGamAccPhi),
  fhGenPi0AccE(qa.fhGenPi0AccE),fhGenPi0AccPt(qa.fhGenPi0AccPt),fhGenPi0AccEta(qa.fhGenPi0AccEta),fhGenPi0AccPhi(qa.fhGenPi0AccPhi),
  fh1pOverE(qa.fh1pOverE),fh1dR(qa.fh1dR),fh2EledEdx(qa.fh2EledEdx), fh2MatchdEdx(qa.fh2MatchdEdx)
{
  // cpy ctor
  
}

//_________________________________________________________________________
AliAnaCalorimeterQA & AliAnaCalorimeterQA::operator = (const AliAnaCalorimeterQA & qa)
{
  // assignment operator

  if(this == &qa)return *this;
  ((AliAnaPartCorrBaseClass *)this)->operator=(qa);

  fCalorimeter  = qa.fCalorimeter;
  fStyleMacro   = qa.fStyleMacro;	
	
  fhE      = qa.fhE;
  fhPt     = qa.fhPt;
  fhPhi    = qa.fhPhi;
  fhEta    = qa.fhEta;
  fhEtaPhi = qa.fhEtaPhi;

  fhECharged      = qa.fhECharged;
  fhPtCharged     = qa.fhPtCharged;
  fhPhiCharged    = qa.fhPhiCharged;
  fhEtaCharged    = qa.fhEtaCharged;
  fhEtaPhiCharged = qa.fhEtaPhiCharged;

  fhEChargedNoOut      = qa.fhEChargedNoOut;
  fhPtChargedNoOut     = qa.fhPtChargedNoOut;
  fhPhiChargedNoOut    = qa.fhPhiChargedNoOut;
  fhEtaChargedNoOut    = qa.fhEtaChargedNoOut;
  fhEtaPhiChargedNoOut = qa.fhEtaPhiChargedNoOut;
	
  fhIM   = qa.fhIM;
  fhAsym = qa.fhAsym;
	
  fhNCellsPerCluster = qa.fhNCellsPerCluster;
  fhNClusters        = qa.fhNClusters;
	
  fhDeltaE   = qa.fhDeltaE;	
  fhDeltaPt  = qa.fhDeltaPt;
  fhDeltaPhi = qa.fhDeltaPhi;
  fhDeltaEta = qa.fhDeltaEta;
	
  fhRatioE   = qa.fhRatioE;	
  fhRatioPt  = qa.fhRatioPt;
  fhRatioPhi = qa.fhRatioPhi;
  fhRatioEta = qa.fhRatioEta;
	
	
  fh2E   = qa.fh2E;	
  fh2Pt  = qa.fh2Pt;
  fh2Phi = qa.fh2Phi;
  fh2Eta = qa.fh2Eta;
	
  fhNCells    = qa.fhNCells;
  fhAmplitude = qa.fhAmplitude;

  fhGenGamPt = qa.fhGenGamPt  ; fhGenGamEta = qa.fhGenGamEta  ; fhGenGamPhi = qa.fhGenGamPhi  ;
  fhGenPi0Pt = qa.fhGenPi0Pt  ; fhGenPi0Eta = qa.fhGenPi0Eta  ; fhGenPi0Phi = qa.fhGenPi0Phi  ;
  fhGenEtaPt = qa.fhGenEtaPt  ; fhGenEtaEta = qa.fhGenEtaEta  ; fhGenEtaPhi = qa.fhGenEtaPhi  ;
  fhGenOmegaPt = qa.fhGenOmegaPt  ; fhGenOmegaEta = qa.fhGenOmegaEta  ; fhGenOmegaPhi = qa.fhGenOmegaPhi  ;
  fhGenElePt = qa.fhGenElePt  ; fhGenEleEta = qa.fhGenEleEta  ; fhGenElePhi = qa.fhGenElePhi ;	

  fhEMVxyz = qa.fhEMVxyz ;  fhEMR = qa.fhEMR ; fhHaVxyz = qa.fhHaVxyz ;  fhHaR = qa.fhHaR ;
  fhGamE = qa.fhGamE ;fhGamPt = qa.fhGamPt ;fhGamPhi = qa.fhGamPhi ;fhGamEta = qa.fhGamEta ; 
  fhGamDeltaE   = qa.fhDeltaE; fhGamDeltaPt  = qa.fhDeltaPt; fhGamDeltaPhi = qa.fhDeltaPhi; fhGamDeltaEta = qa.fhDeltaEta;
	
  fhGamRatioE   = qa.fhGamRatioE;  fhGamRatioPt  = qa.fhGamRatioPt;  fhGamRatioPhi = qa.fhGamRatioPhi;  fhGamRatioEta = qa.fhGamRatioEta;
  fhGam2E   = qa.fhGam2E;	 fhGam2Pt  = qa.fhGam2Pt; fhGam2Phi = qa.fhGam2Phi; fhGam2Eta = qa.fhGam2Eta;

  fhEleE = qa.fhEleE ;fhElePt = qa.fhElePt ;fhElePhi = qa.fhElePhi ;fhEleEta = qa.fhEleEta ;
  fhPi0E = qa.fhPi0E ;fhPi0Pt = qa.fhPi0Pt ;fhPi0Phi = qa.fhPi0Phi ;fhPi0Eta = qa.fhPi0Eta ; 
  fhNeHadE = qa.fhNeHadE ;fhNeHadPt = qa.fhNeHadPt ;fhNeHadPhi = qa.fhNeHadPhi ;fhNeHadEta = qa.fhNeHadEta ; 
  fhChHadE = qa.fhChHadE ;fhChHadPt = qa.fhChHadPt ;fhChHadPhi = qa.fhChHadPhi ;fhChHadEta = qa.fhChHadEta ;
  fhGenGamAccE = qa.fhGenGamAccE ;  fhGenPi0AccE = qa.fhGenPi0AccE ;
  fhGenGamAccPt = qa.fhGenGamAccPt ;fhGenGamAccEta = qa.fhGenGamAccEta ;fhGenGamAccPhi = qa.fhGenGamAccPhi ;
  fhGenPi0AccPt = qa.fhGenPi0AccPt ;fhGenPi0AccEta = qa.fhGenPi0AccEta; fhGenPi0AccPhi = qa.fhGenPi0AccPhi ;	
		
  fh1pOverE = qa.fh1pOverE;
  fh1dR = qa.fh1dR;
  fh2MatchdEdx = qa.fh2MatchdEdx;
  fh2EledEdx = qa.fh2EledEdx;	
	
  return *this;

}

//________________________________________________________________________
TList *  AliAnaCalorimeterQA::GetCreateOutputObjects()
{  
	// Create histograms to be saved in output file and 
	// store them in fOutputContainer
    
	TList * outputContainer = new TList() ; 
	outputContainer->SetName("ExampleHistos") ; 
	
	Int_t nptbins  = GetHistoNPtBins();
	Int_t nphibins = GetHistoNPhiBins();
	Int_t netabins = GetHistoNEtaBins();
	Float_t ptmax  = GetHistoPtMax();
	Float_t phimax = GetHistoPhiMax();
	Float_t etamax = GetHistoEtaMax();
	Float_t ptmin  = GetHistoPtMin();
	Float_t phimin = GetHistoPhiMin();
	Float_t etamin = GetHistoEtaMin();	
	
	fhE  = new TH1F ("hE","E reconstructed", nptbins,ptmin,ptmax); 
	fhE->SetXTitle("E (GeV)");
	outputContainer->Add(fhE);
	
	fhPt  = new TH1F ("hPt","p_{T} reconstructed", nptbins,ptmin,ptmax); 
	fhPt->SetXTitle("p_{T} (GeV/c)");
	outputContainer->Add(fhPt);
	
	fhPhi  = new TH1F ("hPhi","#phi reconstructed",nphibins,phimin,phimax); 
	fhPhi->SetXTitle("#phi (rad)");
	outputContainer->Add(fhPhi);
	
	fhEta  = new TH1F ("hEta","#eta reconstructed",netabins,etamin,etamax); 
	fhEta->SetXTitle("#eta ");
	outputContainer->Add(fhEta);
	
	fhEtaPhi  = new TH2F ("hEtaPhi","#eta vs #phi, reconstructed",netabins,etamin,etamax,nphibins,phimin,phimax); 
	fhEtaPhi->SetXTitle("#eta ");
	fhEtaPhi->SetYTitle("#phi ");
	outputContainer->Add(fhEtaPhi);
	
	fhECharged  = new TH1F ("hECharged","E reconstructed, matched with track", nptbins,ptmin,ptmax); 
	fhECharged->SetXTitle("E (GeV)");
	outputContainer->Add(fhECharged);
	
	fhPtCharged  = new TH1F ("hPtCharged","p_{T} reconstructed, matched with track", nptbins,ptmin,ptmax); 
	fhPtCharged->SetXTitle("p_{T} (GeV/c)");
	outputContainer->Add(fhPtCharged);
	
	fhPhiCharged  = new TH1F ("hPhiCharged","#phi reconstructed, matched with track",nphibins,phimin,phimax); 
	fhPhiCharged->SetXTitle("#phi (rad)");
	outputContainer->Add(fhPhiCharged);
	
	fhEtaCharged  = new TH1F ("hEtaCharged","#eta reconstructed, matched with track",netabins,etamin,etamax); 
	fhEtaCharged->SetXTitle("#eta ");
	outputContainer->Add(fhEtaCharged);
	
	fhEtaPhiCharged  = new TH2F ("hEtaPhiCharged","#eta vs #phi, reconstructed , matched with track",netabins,etamin,etamax,nphibins,phimin,phimax); 
	fhEtaPhiCharged->SetXTitle("#eta ");
	fhEtaPhiCharged->SetYTitle("#phi ");
	outputContainer->Add(fhEtaPhiCharged);	
	

	fhEChargedNoOut  = new TH1F ("hEChargedNoOut","E reconstructed, matched with track, no output track params", nptbins,ptmin,ptmax); 
	fhEChargedNoOut->SetXTitle("E (GeV)");
	outputContainer->Add(fhEChargedNoOut);
	
	fhPtChargedNoOut  = new TH1F ("hPtChargedNoOut","p_{T} reconstructed, matched with track, no output track params", nptbins,ptmin,ptmax); 
	fhPtChargedNoOut->SetXTitle("p_{T} (GeV/c)");
	outputContainer->Add(fhPtChargedNoOut);
	
	fhPhiChargedNoOut  = new TH1F ("hPhiChargedNoOut","#phi reconstructed, matched with track, no output track params",nphibins,phimin,phimax); 
	fhPhiChargedNoOut->SetXTitle("#phi (rad)");
	outputContainer->Add(fhPhiChargedNoOut);
	
	fhEtaChargedNoOut  = new TH1F ("hEtaChargedNoOut","#eta reconstructed, matched with track, no output track params",netabins,etamin,etamax); 
	fhEtaChargedNoOut->SetXTitle("#eta ");
	outputContainer->Add(fhEtaChargedNoOut);
	
	fhEtaPhiChargedNoOut  = new TH2F ("hEtaPhiChargedNoOut","#eta vs #phi, reconstructed , matched with track, no output track params",netabins,etamin,etamax,nphibins,phimin,phimax); 
	fhEtaPhiChargedNoOut->SetXTitle("#eta ");
	fhEtaPhiChargedNoOut->SetYTitle("#phi ");
	outputContainer->Add(fhEtaPhiChargedNoOut);	

	fhIM  = new TH2F ("hIM","Cluster pairs Invariant mass vs reconstructed pair energy",nptbins,ptmin,ptmax,200,0,1); 
	fhIM->SetXTitle("E_{cluster pairs} (GeV) ");
	fhIM->SetYTitle("M_{cluster pairs} (GeV/c^{2})");
	outputContainer->Add(fhIM);
	
	fhAsym  = new TH2F ("hAssym","Cluster pairs Asymmetry vs reconstructed pair energy",nptbins,ptmin,ptmax,100,0,1); 
	fhAsym->SetXTitle("E_{cluster pairs} (GeV) ");
	fhAsym->SetYTitle("Asymmetry");
	outputContainer->Add(fhAsym);	
	
	fhNCellsPerCluster  = new TH2F ("hNCellsPerCluster","# cells per cluster vs cluster energy", nptbins,ptmin,ptmax, 100,0,1000); 
	fhNCellsPerCluster->SetXTitle("E (GeV)");
	fhNCellsPerCluster->SetYTitle("n cells");
	outputContainer->Add(fhNCellsPerCluster);
	
	fhNClusters  = new TH1F ("hNClusters","# clusters", 300,0,300); 
	fhNClusters->SetXTitle("number of clusters");
	outputContainer->Add(fhNClusters);
	
	//Calo cells
	fhNCells  = new TH1F ("hNCells","# cells", 13000,0,13000); 
	fhNCells->SetXTitle("n cells");
	outputContainer->Add(fhNCells);
    
	fhAmplitude  = new TH1F ("hAmplitude","Amplitude", 100,0,1000); 
	fhAmplitude->SetXTitle("Amplitude ");
	outputContainer->Add(fhAmplitude);
	
	if(IsDataMC()){
		
		fhDeltaE  = new TH1F ("hDeltaE","MC - Reco E ", 200,-50,50); 
		fhDeltaE->SetXTitle("#Delta E (GeV)");
		outputContainer->Add(fhDeltaE);
		
		fhDeltaPt  = new TH1F ("hDeltaPt","MC - Reco p_{T} ", 200,-50,50); 
		fhDeltaPt->SetXTitle("#Delta p_{T} (GeV/c)");
		outputContainer->Add(fhDeltaPt);
		
		fhDeltaPhi  = new TH1F ("hDeltaPhi","MC - Reco #phi ",100,-2,2); 
		fhDeltaPhi->SetXTitle("#Delta #phi (rad)");
		outputContainer->Add(fhDeltaPhi);
		
		fhDeltaEta  = new TH1F ("hDeltaEta","MC- Reco #eta",100,-1,1); 
		fhDeltaEta->SetXTitle("#Delta #eta ");
		outputContainer->Add(fhDeltaEta);
		
		fhRatioE  = new TH1F ("hRatioE","Reco/MC E ", 200,0,2); 
		fhRatioE->SetXTitle("E_{reco}/E_{gen}");
		outputContainer->Add(fhRatioE);
		
		fhRatioPt  = new TH1F ("hRatioPt","Reco/MC p_{T} ", 200,0,2); 
		fhRatioPt->SetXTitle("p_{T, reco}/p_{T, gen}");
		outputContainer->Add(fhRatioPt);
		
		fhRatioPhi  = new TH1F ("hRatioPhi","Reco/MC #phi ",200,0,2); 
		fhRatioPhi->SetXTitle("#phi_{reco}/#phi_{gen}");
		outputContainer->Add(fhRatioPhi);
		
		fhRatioEta  = new TH1F ("hRatioEta","Reco/MC #eta",200,0,2); 
		fhRatioEta->SetXTitle("#eta_{reco}/#eta_{gen} ");
		outputContainer->Add(fhRatioEta);
		
		fh2E  = new TH2F ("h2E","E distribution, reconstructed vs generated", nptbins,ptmin,ptmax,nptbins,ptmin,ptmax); 
		fh2E->SetXTitle("E_{rec} (GeV)");
		fh2E->SetYTitle("E_{gen} (GeV)");
		outputContainer->Add(fh2E);	  
		
		fh2Pt  = new TH2F ("h2Pt","p_T distribution, reconstructed vs generated", nptbins,ptmin,ptmax,nptbins,ptmin,ptmax); 
		fh2Pt->SetXTitle("p_{T,rec} (GeV/c)");
		fh2Pt->SetYTitle("p_{T,gen} (GeV/c)");
		outputContainer->Add(fh2Pt);
		
		fh2Phi  = new TH2F ("h2Phi","#phi distribution, reconstructed vs generated", nphibins,phimin,phimax, nphibins,phimin,phimax); 
		fh2Phi->SetXTitle("#phi_{rec} (rad)");
		fh2Phi->SetYTitle("#phi_{gen} (rad)");
		outputContainer->Add(fh2Phi);
		
		fh2Eta  = new TH2F ("h2Eta","#eta distribution, reconstructed vs generated", netabins,etamin,etamax,netabins,etamin,etamax); 
		fh2Eta->SetXTitle("#eta_{rec} ");
		fh2Eta->SetYTitle("#eta_{gen} ");
		outputContainer->Add(fh2Eta);
		
		//Fill histos depending on origin of cluster
		fhGamE  = new TH2F ("hGamE","E reconstructed vs E generated from #gamma", nptbins,ptmin,ptmax, nptbins,ptmin,ptmax); 
		fhGamE->SetXTitle("E (GeV)");
		outputContainer->Add(fhGamE);
		
		fhGamPt  = new TH2F ("hGamPt","p_{T} reconstructed vs E generated from #gamma", nptbins,ptmin,ptmax, nptbins,ptmin,ptmax); 
		fhGamPt->SetXTitle("p_{T} (GeV/c)");
		outputContainer->Add(fhGamPt);
		
		fhGamPhi  = new TH2F ("hGamPhi","#phi reconstructed vs E generated from #gamma",nphibins,phimin,phimax,nphibins,phimin,phimax); 
		fhGamPhi->SetXTitle("#phi (rad)");
		outputContainer->Add(fhGamPhi);
		
		fhGamEta  = new TH2F ("hGamEta","#eta reconstructed vs E generated from #gamma",netabins,etamin,etamax,netabins,etamin,etamax); 
		fhGamEta->SetXTitle("#eta ");
		outputContainer->Add(fhGamEta);
		
		fhGamDeltaE  = new TH1F ("hGamDeltaE","#gamma MC - Reco E ", 200,-50,50); 
		fhGamDeltaE->SetXTitle("#Delta E (GeV)");
		outputContainer->Add(fhGamDeltaE);
		
		fhGamDeltaPt  = new TH1F ("hGamDeltaPt","#gamma MC - Reco p_{T} ", 200,-50,50); 
		fhGamDeltaPt->SetXTitle("#Delta p_{T} (GeV/c)");
		outputContainer->Add(fhGamDeltaPt);
		
		fhGamDeltaPhi  = new TH1F ("hGamDeltaPhi","#gamma MC - Reco #phi ",100,-2,2); 
		fhGamDeltaPhi->SetXTitle("#Delta #phi (rad)");
		outputContainer->Add(fhGamDeltaPhi);
		
		fhGamDeltaEta  = new TH1F ("hGamDeltaEta","#gamma MC- Reco #eta",100,-1,1); 
		fhGamDeltaEta->SetXTitle("#Delta #eta ");
		outputContainer->Add(fhGamDeltaEta);
		
		fhGamRatioE  = new TH1F ("hGamRatioE","#gamma Reco/MC E ", 200,0,2); 
		fhGamRatioE->SetXTitle("E_{reco}/E_{gen}");
		outputContainer->Add(fhGamRatioE);
		
		fhGamRatioPt  = new TH1F ("hGamRatioPt","#gamma Reco/MC p_{T} ", 200,0,2); 
		fhGamRatioPt->SetXTitle("p_{T, reco}/p_{T, gen}");
		outputContainer->Add(fhGamRatioPt);
		
		fhGamRatioPhi  = new TH1F ("hGamRatioPhi","#gamma Reco/MC #phi ",200,0,2); 
		fhGamRatioPhi->SetXTitle("#phi_{reco}/#phi_{gen}");
		outputContainer->Add(fhGamRatioPhi);
		
		fhGamRatioEta  = new TH1F ("hGamRatioEta","#gamma Reco/MC #eta",200,0,2); 
		fhGamRatioEta->SetXTitle("#eta_{reco}/#eta_{gen} ");
		outputContainer->Add(fhGamRatioEta);

		fhGam2E  = new TH2F ("hGam2E","#gamma E distribution, reconstructed vs generated", nptbins,ptmin,ptmax,nptbins,ptmin,ptmax); 
		fhGam2E->SetXTitle("E_{rec} (GeV)");
		fhGam2E->SetYTitle("E_{gen} (GeV)");
		outputContainer->Add(fhGam2E);	  
		
		fhGam2Pt  = new TH2F ("hGam2Pt","#gamma p_T distribution, reconstructed vs generated", nptbins,ptmin,ptmax,nptbins,ptmin,ptmax); 
		fhGam2Pt->SetXTitle("p_{T,rec} (GeV/c)");
		fhGam2Pt->SetYTitle("p_{T,gen} (GeV/c)");
		outputContainer->Add(fhGam2Pt);
		
		fhGam2Phi  = new TH2F ("hGam2Phi","#gamma #phi distribution, reconstructed vs generated", nphibins,phimin,phimax, nphibins,phimin,phimax); 
		fhGam2Phi->SetXTitle("#phi_{rec} (rad)");
		fhGam2Phi->SetYTitle("#phi_{gen} (rad)");
		outputContainer->Add(fhGam2Phi);
		
		fhGam2Eta  = new TH2F ("hGam2Eta","#gamma #eta distribution, reconstructed vs generated", netabins,etamin,etamax,netabins,etamin,etamax); 
		fhGam2Eta->SetXTitle("#eta_{rec} ");
		fhGam2Eta->SetYTitle("#eta_{gen} ");
		outputContainer->Add(fhGam2Eta);

		fhPi0E  = new TH2F ("hPi0E","E reconstructed vs E generated from #pi^0", nptbins,ptmin,ptmax, nptbins,ptmin,ptmax); 
		fhPi0E->SetXTitle("E (GeV)");
		outputContainer->Add(fhPi0E);
		
		fhPi0Pt  = new TH2F ("hPi0Pt","p_{T} reconstructed vs E generated from #pi^{0}", nptbins,ptmin,ptmax, nptbins,ptmin,ptmax); 
		fhPi0Pt->SetXTitle("p_{T} (GeV/c)");
		outputContainer->Add(fhPi0Pt);
		
		fhPi0Phi  = new TH2F ("hPi0Phi","#phi reconstructed vs E generated from #pi^{0}",nphibins,phimin,phimax,nphibins,phimin,phimax); 
		fhPi0Phi->SetXTitle("#phi (rad)");
		outputContainer->Add(fhPi0Phi);
		
		fhPi0Eta  = new TH2F ("hPi0Eta","#eta reconstructed vs E generated from #pi^{0}",netabins,etamin,etamax,netabins,etamin,etamax); 
		fhPi0Eta->SetXTitle("#eta ");
		outputContainer->Add(fhPi0Eta);
		
		fhEleE  = new TH2F ("hEleE","E reconstructed vs E generated from e^{#pm}", nptbins,ptmin,ptmax, nptbins,ptmin,ptmax); 
		fhEleE->SetXTitle("E (GeV)");
		outputContainer->Add(fhEleE);		
		
		fhElePt  = new TH2F ("hElePt","p_{T} reconstructed vs E generated from e^{#pm}", nptbins,ptmin,ptmax, nptbins,ptmin,ptmax); 
		fhElePt->SetXTitle("p_{T} (GeV/c)");
		outputContainer->Add(fhElePt);
		
		fhElePhi  = new TH2F ("hElePhi","#phi reconstructed vs E generated frome^{#pm}",nphibins,phimin,phimax,nphibins,phimin,phimax); 
		fhElePhi->SetXTitle("#phi (rad)");
		outputContainer->Add(fhElePhi);
		
		fhEleEta  = new TH2F ("hEleEta","#eta reconstructed vs E generated from e^{#pm}",netabins,etamin,etamax,netabins,etamin,etamax); 
		fhEleEta->SetXTitle("#eta ");
		outputContainer->Add(fhEleEta);
		
		fhNeHadE  = new TH2F ("hNeHadE","E reconstructed vs E generated from #pi^0", nptbins,ptmin,ptmax, nptbins,ptmin,ptmax); 
		fhNeHadE->SetXTitle("E (GeV)");
		outputContainer->Add(fhNeHadE);
		
		fhNeHadPt  = new TH2F ("hNeHadPt","p_{T} reconstructed vs E generated from neutral hadron", nptbins,ptmin,ptmax, nptbins,ptmin,ptmax); 
		fhNeHadPt->SetXTitle("p_{T} (GeV/c)");
		outputContainer->Add(fhNeHadPt);
		
		fhNeHadPhi  = new TH2F ("hNeHadPhi","#phi reconstructed vs E generated from neutral hadron",nphibins,phimin,phimax,nphibins,phimin,phimax); 
		fhNeHadPhi->SetXTitle("#phi (rad)");
		outputContainer->Add(fhNeHadPhi);
		
		fhNeHadEta  = new TH2F ("hNeHadEta","#eta reconstructed vs E generated from neutral hadron",netabins,etamin,etamax,netabins,etamin,etamax); 
		fhNeHadEta->SetXTitle("#eta ");
		outputContainer->Add(fhNeHadEta);
		
		fhChHadE  = new TH2F ("hChHadE","E reconstructed vs E generated from #pi^0", nptbins,ptmin,ptmax, nptbins,ptmin,ptmax); 
		fhChHadE->SetXTitle("E (GeV)");
		outputContainer->Add(fhChHadE);
		
		fhChHadPt  = new TH2F ("hChHadPt","p_{T} reconstructed vs E generated from charged hadron", nptbins,ptmin,ptmax, nptbins,ptmin,ptmax); 
		fhChHadPt->SetXTitle("p_{T} (GeV/c)");
		outputContainer->Add(fhChHadPt);
		
		fhChHadPhi  = new TH2F ("hChHadPhi","#phi reconstructed vs E generated from charged hadron",nphibins,phimin,phimax,nphibins,phimin,phimax); 
		fhChHadPhi->SetXTitle("#phi (rad)");
		outputContainer->Add(fhChHadPhi);
		
		fhChHadEta  = new TH2F ("hChHadEta","#eta reconstructed vs E generated from charged hadron",netabins,etamin,etamax,netabins,etamin,etamax); 
		fhChHadEta->SetXTitle("#eta ");
		outputContainer->Add(fhChHadEta);
		
		//Vertex of generated particles 
		
		fhEMVxyz  = new TH2F ("hEMVxyz","Production vertex of reconstructed ElectroMagnetic particles",100,0,500,100,0,500);//,100,0,500); 
		fhEMVxyz->SetXTitle("v_{x}");
		fhEMVxyz->SetYTitle("v_{y}");
		//fhEMVxyz->SetZTitle("v_{z}");
		outputContainer->Add(fhEMVxyz);
		
		fhHaVxyz  = new TH2F ("hHaVxyz","Production vertex of reconstructed hadrons",100,0,500,100,0,500);//,100,0,500); 
		fhHaVxyz->SetXTitle("v_{x}");
		fhHaVxyz->SetYTitle("v_{y}");
		//fhHaVxyz->SetZTitle("v_{z}");
		outputContainer->Add(fhHaVxyz);
		
		fhEMR  = new TH2F ("hEMR","Distance to production vertex of reconstructed ElectroMagnetic particles vs E rec",nptbins,ptmin,ptmax,100,0,500); 
		fhEMR->SetXTitle("E (GeV)");
		fhEMR->SetYTitle("TMath::Sqrt(v_{x}^{2}+v_{y}^{2})");
		outputContainer->Add(fhEMR);
		
		fhHaR  = new TH2F ("hHaR","Distance to production vertex of reconstructed Hadrons vs E rec",nptbins,ptmin,ptmax,100,0,500); 
		fhHaR->SetXTitle("E (GeV)");
		fhHaR->SetYTitle("TMath::Sqrt(v_{x}^{2}+v_{y}^{2})");
		outputContainer->Add(fhHaR);
		
		
		//Pure MC
		fhGenGamPt  = new TH1F("hGenGamPt" ,"p_{T} of generated #gamma",nptbins,ptmin,ptmax);
		fhGenGamEta = new TH1F("hGenGamEta","Y of generated #gamma",netabins,etamin,etamax);
		fhGenGamPhi = new TH1F("hGenGamPhi","#phi of generated #gamma",nphibins,phimin,phimax);
		
		fhGenPi0Pt  = new TH1F("hGenPi0Pt" ,"p_{T} of generated #pi^{0}",nptbins,ptmin,ptmax);
		fhGenPi0Eta = new TH1F("hGenPi0Eta","Y of generated #pi^{0}",netabins,etamin,etamax);
		fhGenPi0Phi = new TH1F("hGenPi0Phi","#phi of generated #pi^{0}",nphibins,phimin,phimax);
		
		fhGenEtaPt  = new TH1F("hGenEtaPt" ,"p_{T} of generated #eta",nptbins,ptmin,ptmax);
		fhGenEtaEta = new TH1F("hGenEtaEta","Y of generated #eta",netabins,etamin,etamax);
		fhGenEtaPhi = new TH1F("hGenEtaPhi","#phi of generated #eta",nphibins,phimin,phimax);
		
		fhGenOmegaPt  = new TH1F("hGenOmegaPt" ,"p_{T} of generated #omega",nptbins,ptmin,ptmax);
		fhGenOmegaEta = new TH1F("hGenOmegaEta","Y of generated #omega",netabins,etamin,etamax);
		fhGenOmegaPhi = new TH1F("hGenOmegaPhi","#phi of generated #omega",nphibins,phimin,phimax);		
		
		fhGenElePt  = new TH1F("hGenElePt" ,"p_{T} of generated e^{#pm}",nptbins,ptmin,ptmax);
		fhGenEleEta = new TH1F("hGenEleEta","Y of generated  e^{#pm}",netabins,etamin,etamax);
		fhGenElePhi = new TH1F("hGenElePhi","#phi of generated  e^{#pm}",nphibins,phimin,phimax);		
		
		fhGenGamPt->SetXTitle("p_{T} (GeV/c)");
		fhGenGamEta->SetXTitle("#eta");
		fhGenGamPhi->SetXTitle("#phi (rad)");
		outputContainer->Add(fhGenGamPt);
		outputContainer->Add(fhGenGamEta);
		outputContainer->Add(fhGenGamPhi);

		fhGenPi0Pt->SetXTitle("p_{T} (GeV/c)");
		fhGenPi0Eta->SetXTitle("#eta");
		fhGenPi0Phi->SetXTitle("#phi (rad)");
		outputContainer->Add(fhGenPi0Pt);
		outputContainer->Add(fhGenPi0Eta);
		outputContainer->Add(fhGenPi0Phi);
		
		fhGenEtaPt->SetXTitle("p_{T} (GeV/c)");
		fhGenEtaEta->SetXTitle("#eta");
		fhGenEtaPhi->SetXTitle("#phi (rad)");
		outputContainer->Add(fhGenEtaPt);
		outputContainer->Add(fhGenEtaEta);
		outputContainer->Add(fhGenEtaPhi);
				
		fhGenOmegaPt->SetXTitle("p_{T} (GeV/c)");
		fhGenOmegaEta->SetXTitle("#eta");
		fhGenOmegaPhi->SetXTitle("#phi (rad)");
		outputContainer->Add(fhGenOmegaPt);
		outputContainer->Add(fhGenOmegaEta);
		outputContainer->Add(fhGenOmegaPhi);
		
		fhGenElePt->SetXTitle("p_{T} (GeV/c)");
		fhGenEleEta->SetXTitle("#eta");
		fhGenElePhi->SetXTitle("#phi (rad)");
		outputContainer->Add(fhGenElePt);
		outputContainer->Add(fhGenEleEta);
		outputContainer->Add(fhGenElePhi);
		
		fhGenGamAccE   = new TH1F("hGenGamAccE" ,"E of generated #gamma in calorimeter acceptance",nptbins,ptmin,ptmax);
		fhGenGamAccPt  = new TH1F("hGenGamAccPt" ,"p_{T} of generated #gamma in calorimeter acceptance",nptbins,ptmin,ptmax);
		fhGenGamAccEta = new TH1F("hGenGamAccEta","Y of generated #gamma in calorimeter acceptance",netabins,etamin,etamax);
		fhGenGamAccPhi = new TH1F("hGenGamAccPhi","#phi of generated #gamma  in calorimeter acceptance",nphibins,phimin,phimax);
		
		fhGenPi0AccE  = new TH1F("hGenPi0AccE" ,"E of generated #pi^{0} in calorimeter acceptance",nptbins,ptmin,ptmax);
		fhGenPi0AccPt  = new TH1F("hGenPi0AccPt" ,"p_{T} of generated #pi^{0} in calorimeter acceptance",nptbins,ptmin,ptmax);
		fhGenPi0AccEta = new TH1F("hGenPi0AccEta","Y of generated #pi^{0} in calorimeter acceptance",netabins,etamin,etamax);
		fhGenPi0AccPhi = new TH1F("hGenPi0AccPhi","#phi of generated #pi^{0} in calorimeter acceptance",nphibins,phimin,phimax);
		
		fhGenGamAccE  ->SetXTitle("E (GeV)");
		fhGenGamAccPt ->SetXTitle("p_{T} (GeV/c)");
		fhGenGamAccEta->SetXTitle("#eta");
		fhGenGamAccPhi->SetXTitle("#phi (rad)");
		outputContainer->Add(fhGenGamAccE);		
		outputContainer->Add(fhGenGamAccPt);
		outputContainer->Add(fhGenGamAccEta);
		outputContainer->Add(fhGenGamAccPhi);
		
		fhGenPi0AccE  ->SetXTitle("E (GeV)");		
		fhGenPi0AccPt ->SetXTitle("p_{T} (GeV/c)");
		fhGenPi0AccEta->SetXTitle("#eta");
		fhGenPi0AccPhi->SetXTitle("#phi (rad)");
		outputContainer->Add(fhGenPi0AccE);		
		outputContainer->Add(fhGenPi0AccPt);
		outputContainer->Add(fhGenPi0AccEta);
		outputContainer->Add(fhGenPi0AccPhi);
		
	}
	
	
	fh1pOverE = new TH1F("h1pOverE","EMCAL-TRACK matches p/E",100,0.,10.);
	fh1pOverE->SetXTitle("p/E");
	outputContainer->Add(fh1pOverE);
	
	fh1dR = new TH1F("h1dR","EMCAL-TRACK matches dR",300, 0.,TMath::Pi());
	fh1dR->SetXTitle("#Delta R (rad)");
	outputContainer->Add(fh1dR) ;
	
	fh2MatchdEdx = new TH2F("h2MatchdEdx","dE/dx vs. p for all matches",200,0.,50.,200,0.,400.);
	fh2MatchdEdx->SetXTitle("p (GeV/c)");
	fh2MatchdEdx->SetYTitle("<dE/dx>");
	outputContainer->Add(fh2MatchdEdx);
	
	fh2EledEdx = new TH2F("h2EledEdx","dE/dx vs. p for electrons",200,0.,50.,200,0.,400.);
	fh2EledEdx->SetXTitle("p (GeV/c)");
	fh2EledEdx->SetYTitle("<dE/dx>");
	outputContainer->Add(fh2EledEdx) ;
	
	return outputContainer;
}

//__________________________________________________
void AliAnaCalorimeterQA::Init()
{ 
	//Check if the data or settings are ok
	if(fCalorimeter != "PHOS" && fCalorimeter !="EMCAL"){
		printf("AliAnaCalorimeterQA::Init() - Wrong calorimeter name <%s>, END\n", fCalorimeter.Data());
		abort();
	}	
	
	if(GetReader()->GetDataType()== AliCaloTrackReader::kMC){
		printf("AliAnaCalorimeterQA::Init() - Analysis of reconstructed data, MC reader not aplicable\n");
		abort();
	}	
	
}


//__________________________________________________
void AliAnaCalorimeterQA::InitParameters()
{ 
  //Initialize the parameters of the analysis.
  AddToHistogramsName("AnaCaloQA_");

  fCalorimeter = "EMCAL"; //or PHOS
  fStyleMacro = "" ;
}

//__________________________________________________________________
void AliAnaCalorimeterQA::Print(const Option_t * opt) const
{
  //Print some relevant parameters set for the analysis
  if(! opt)
    return;
  
  printf("**** Print %s %s ****\n", GetName(), GetTitle() ) ;
  AliAnaPartCorrBaseClass::Print(" ");

  printf("Select Calorimeter %s \n",fCalorimeter.Data());
  printf("Plots style macro %s \n",fStyleMacro.Data()); 
} 

//__________________________________________________________________
void  AliAnaCalorimeterQA::MakeAnalysisFillHistograms() 
{
	//Do analysis and fill histograms
	TLorentzVector mom ;
	TLorentzVector mom2 ;
	//Play with the MC stack if available
	AliStack * stack = 0x0;
	TParticle * primary = 0x0;
	if(IsDataMC()) stack =  GetMCStack() ;
	
	//Get List with clusters  
	TRefArray * partList = new TRefArray();
	if(fCalorimeter == "EMCAL") partList = GetAODEMCAL();
	else if(fCalorimeter == "PHOS") partList = GetAODPHOS();
	else {
		printf("AliAnaCalorimeterQA::MakeAnalysisFillHistograms() - Wrong calorimeter name <%s>, END\n", fCalorimeter.Data());
		abort();
	}
	
	if(!partList || partList->GetEntriesFast() == 0) return ;
	
	Int_t nclusters = partList->GetEntriesFast() ; 
	fhNClusters->Fill(nclusters);
	
	if(GetDebug() > 0)
		printf("AliAnaCalorimeterQA::MakeAnalysisFillHistograms() - In %s there are %d clusters \n", fCalorimeter.Data(), nclusters);
	
	//Get vertex for photon momentum calculation
	Double_t v[3] ; //vertex ;
	GetReader()->GetVertex(v);
	
	for(Int_t iclus = 0; iclus < partList->GetEntriesFast(); iclus++){
		
		AliAODCaloCluster * calo =  (AliAODCaloCluster*) (partList->At(iclus));
		//if(calo->GetNCells() <= 2) continue;
		//Get cluster kinematics
		calo->GetMomentum(mom,v);
		Float_t e   = mom.E();
		//if(e < 0.5) continue;
		//printf("e %2.2f, n %f\n",e, calo->GetNCells());
		Float_t pt  = mom.Pt();
		Float_t eta = mom.Eta();
		Float_t phi = mom.Phi();
		if(phi < 0) phi +=TMath::TwoPi();
		
		fhE     ->Fill(e);		
		fhPt    ->Fill(pt);
		fhPhi   ->Fill(phi);
		fhEta   ->Fill(eta);
		fhEtaPhi->Fill(eta,phi);

		//matched cluster with tracks
		if(calo->GetNTracksMatched() > 0){
			fhECharged     ->Fill(e);		
			fhPtCharged    ->Fill(pt);
			fhPhiCharged   ->Fill(phi);
			fhEtaCharged   ->Fill(eta);
			fhEtaPhiCharged->Fill(eta,phi);				
			if((!strcmp(GetReader()->GetInputEvent()->GetName(),"AliESDEvent"))) {
				AliESDEvent *esd = (AliESDEvent*) GetReader()->GetInputEvent();
				AliESDCaloCluster * esdcalo = (AliESDCaloCluster*) esd->GetCaloCluster(calo->GetID());
				Int_t trackIndex = esdcalo->GetTrackMatched();
				//printf("track index %d ntracks %d\n", trackIndex, esd->GetNumberOfTracks());
				if(trackIndex >= 0){
					AliESDtrack* track = (AliESDtrack*) esd->GetTrack(trackIndex);
					if (track && track->GetOuterParam() ) {
				
						Double_t tphi = track->GetOuterParam()->Phi();
						Double_t teta = track->GetOuterParam()->Eta();
						Double_t tmom = track->GetOuterParam()->P();
						Double_t deta = teta - eta;
						Double_t dphi = tphi - phi;
						if(dphi > TMath::Pi()) dphi -= 2*TMath::Pi();
						if(dphi < -TMath::Pi()) dphi += 2*TMath::Pi();
						Double_t dR = sqrt(dphi*dphi + deta*deta);
						
						Double_t pOverE = tmom/e;
						
						fh1pOverE->Fill(pOverE);
						fh1dR->Fill(dR);
						fh2MatchdEdx->Fill(track->P(),track->GetTPCsignal());
						
						int nITS = track->GetNcls(0);
						int nTPC = track->GetNcls(1);
						if(dR < 0.02 && pOverE > 0.5 && pOverE < 1.5
						   && calo->GetNCells() > 1 && nITS > 3 && nTPC > 20) {
							fh2EledEdx->Fill(track->P(),track->GetTPCsignal());
							//cout<<track->P()<<" "<<track->P()<<endl;
						}
					} 
					else if(!track->GetOuterParam()){
					  ULong_t status=AliESDtrack::kTPCrefit;
					  status|=AliESDtrack::kITSrefit;
					  //printf("track status %d\n", track->GetStatus() );
					  fhEChargedNoOut     ->Fill(e);		
					  fhPtChargedNoOut     ->Fill(pt);
					  fhPhiChargedNoOut    ->Fill(phi);
					  fhEtaChargedNoOut    ->Fill(eta);
					  fhEtaPhiChargedNoOut ->Fill(eta,phi);	
					  if ( (track->GetStatus() & status) == status) printf("ITS+TPC\n");
					}
					else {
						Printf("ERROR: Could not receive track %d", trackIndex);
					}
				}// non negative track index
			}//do only if input are ESDs
		}// at least one track matched
		
		if(IsDataMC()){
	        //Play with the MC stack if available
		    Int_t label = calo->GetLabel(0);
			
	        if(label < 0 || !stack) {
				printf("AliAnaCalorimeterQA::MakeAnalysisFillHistograms() *** bad label or no stack ***:  label %d \n", label);
				continue;
	        }
	        
	        if(label >=  stack->GetNtrack()) {
				printf("AliAnaCalorimeterQA::MakeAnalysisFillHistograms() *** large label ***:  label %d, n tracks %d \n", label, stack->GetNtrack());
				continue ;
	        }
	        //cout<<"LABEL > "<<label<<endl;
	        primary = GetMCStack()->Particle(label);
		    Int_t pdg = primary->GetPdgCode();
		    Float_t vx = primary->Vx();
		    Float_t vy = primary->Vy();
		    //Float_t vz = primary->Vz();
		    Float_t r = TMath::Sqrt(vx*vx + vy*vy);
		    if((pdg == 22 || TMath::Abs(pdg)==11) && primary->GetStatusCode()!=1) {
		      fhEMVxyz   ->Fill(vx,vy);//,vz);
		      fhEMR      ->Fill(e,r);
		    }
		    
		    //			printf("reco e %f, pt %f, phi %f, eta %f \n", e, pt, phi, eta);
		    //			printf("prim e %f, pt %f, phi %f, eta %f \n", primary->Energy(),primary->Pt() ,primary->Phi() ,primary->Eta() );
		    //			printf("vertex: vx %f, vy %f, vz %f, r %f \n", vx, vy, vz, r);
		    
		    //Get final particle, no conversion products
		    Int_t status =  primary->GetStatusCode();
		    Int_t mother= primary->GetFirstMother();
		    Int_t finallabel = -1;
		    if(status == 0){
		      while(mother >= 0){
			primary = GetMCStack()->Particle(mother);
			finallabel = mother;
			status = primary->GetStatusCode();
			mother= primary->GetFirstMother();
			if(status == 1) break;		   
			//printf("mother %d\n",mother);
		      }	
		    }
		    
			fh2E      ->Fill(e, primary->Energy());
			fh2Pt     ->Fill(pt, primary->Pt());
			fh2Phi    ->Fill(phi, primary->Phi());
			fh2Eta    ->Fill(eta, primary->Eta());
			fhDeltaE  ->Fill(primary->Energy()-e);
			fhDeltaPt ->Fill(primary->Pt()-pt);
			fhDeltaPhi->Fill(primary->Phi()-phi);
			fhDeltaEta->Fill(primary->Eta()-eta);
			if(primary->Energy() > 0) fhRatioE  ->Fill(e/primary->Energy());
		    if(primary->Pt()     > 0) fhRatioPt ->Fill(pt/primary->Pt());
		    if(primary->Phi()    > 0) fhRatioPhi->Fill(phi/primary->Phi());
		    if(primary->Eta()    > 0) fhRatioEta->Fill(eta/primary->Eta());			
			
		    //cout<<"Final Label "<<finallabel<<" mother "<<mother<<endl;
		    pdg = primary->GetPdgCode();
			Double_t  charge = TDatabasePDG::Instance()->GetParticle(pdg)->Charge();

		    if(pdg == 22){
		      //cout<<"pdg "<<pdg<<" status "<<status<<" "<<primary->GetStatusCode()<<endl;	
		      TParticle *pi0 = 0x0;
		      TParticle *p1 = 0x0;
		      TParticle *p2 = 0x0;
		      TParticle * tmp = 0x0;
		      Int_t pdgpi0 = 0;
		      Int_t mothertmp = 0;
		      //Check if it is a decay photon from pi0, in this case, check if both decay contribute cluster
		      Bool_t ok1 = kFALSE;
		      Bool_t ok2 = kFALSE;
		      
		      if(calo->GetNLabel() > 1){
			if(pdg !=111){
			  while(mother >= 0){
			    pi0 = GetMCStack()->Particle(mother);
			    pdgpi0 = pi0->GetPdgCode();
			    if(pdgpi0 == 111) break;
			    mother= pi0->GetFirstMother();
			    //printf("mother %d\n",mother);
			  }
			}
			else   pi0 = primary;
			//cout<<"MOTHER PI0 LABEL "<<mother<<" pt" << pi0->Pt()<<endl;
			if(pi0->GetNDaughters() == 2){
			  //  cout<<"pi0, 2 daughters "<<endl;
			  Int_t id1 = pi0->GetFirstDaughter();
			  Int_t id2 = pi0->GetFirstDaughter()+1;
			  p1=GetMCStack()->Particle(id1);
			  p2=GetMCStack()->Particle(id2);
			  
			  if(p1->GetFirstMother()!=p2->GetFirstMother()) cout <<"Decay photon mothers are not the same!!"<<endl;
			  if(p1->GetPdgCode()==22 && p2->GetPdgCode()==22){
			    // cout<<"2 photons, labels "<< id1<<" "<<id2<<endl;
			    for(UInt_t ilabel = 0; ilabel < calo->GetNLabel(); ilabel++){
			      Int_t iprim = calo->GetLabel(ilabel);
			      //cout<<"iprim "<<iprim<<endl;
			      if (iprim == id1) ok1 = kTRUE;
			      else if (iprim == id2) ok2 = kTRUE;
			      mothertmp = iprim;
			      while(mothertmp >= 0){
				tmp = GetMCStack()->Particle(mothertmp);				    
				mothertmp= tmp->GetFirstMother();
				//	cout<<"mothertmp "<<mothertmp<<" "<<tmp->GetName()<< " pt "<<tmp->Pt()<<endl;
				if (mothertmp == id1) ok1 = kTRUE;
				else if (mothertmp == id2) ok2 = kTRUE;
				if(ok1 && ok2) break;
				//printf("mother %d\n",mother);
			      }
			    } 
			  }//2 photon daughters
			}////mother pi0, 2 daughers
		      }//more than one contribution to clust
		      
		      if(ok1 && ok2){
			//cout<<"Fill pi0"<< "E  "<< e <<" prim E "<<primary->Energy()<<endl;
			fhPi0E     ->Fill(e,primary->Energy());	
			fhPi0Pt    ->Fill(pt,primary->Pt());
			fhPi0Eta   ->Fill(eta,primary->Eta());	
			fhPi0Phi   ->Fill(phi,primary->Phi());					
		      }
		      else{
			fhGamE     ->Fill(e,primary->Energy());	
			fhGamPt    ->Fill(pt,primary->Pt());
			fhGamEta   ->Fill(eta,primary->Eta());	
			fhGamPhi   ->Fill(phi,primary->Phi());
			fhGamDeltaE  ->Fill(primary->Energy()-e);
			fhGamDeltaPt ->Fill(primary->Pt()-pt);
			fhGamDeltaPhi->Fill(primary->Phi()-phi);
			fhGamDeltaEta->Fill(primary->Eta()-eta);
			if(primary->Energy() > 0) fhGamRatioE  ->Fill(e/primary->Energy());
			if(primary->Pt()     > 0) fhGamRatioPt ->Fill(pt/primary->Pt());
			if(primary->Phi()    > 0) fhGamRatioPhi->Fill(phi/primary->Phi());
			if(primary->Eta()    > 0) fhGamRatioEta->Fill(eta/primary->Eta());
			fhGam2E      ->Fill(e, primary->Energy());
			fhGam2Pt     ->Fill(pt, primary->Pt());
			fhGam2Phi    ->Fill(phi, primary->Phi());
			fhGam2Eta    ->Fill(eta, primary->Eta());
		      }
		    }//pdg == 22
		    else if(TMath::Abs(pdg) == 11) {
		      fhEleE     ->Fill(e,primary->Energy());	
		      fhElePt    ->Fill(pt,primary->Pt());
		      fhEleEta   ->Fill(eta,primary->Eta());	
		      fhElePhi   ->Fill(phi,primary->Phi());
		      fhEMVxyz   ->Fill(vx,vy);//,vz);
		      fhEMR      ->Fill(e,r);
		    }
		    else if(pdg == 111) {
		      fhPi0E     ->Fill(e,primary->Energy());	
		      fhPi0Pt    ->Fill(pt,primary->Pt());
		      fhPi0Eta   ->Fill(eta,primary->Eta());	
		      fhPi0Phi   ->Fill(phi,primary->Phi());
		    }
		    else if(charge == 0){
		      fhNeHadE     ->Fill(e,primary->Energy());	
		      fhNeHadPt    ->Fill(pt,primary->Pt());
		      fhNeHadEta   ->Fill(eta,primary->Eta());	
		      fhNeHadPhi   ->Fill(phi,primary->Phi());	
		      fhHaVxyz     ->Fill(vx,vy);//,vz);
		      fhHaR        ->Fill(e,r);
		    }
		    else if(charge!=0){
		      fhChHadE     ->Fill(e,primary->Energy());	
		      fhChHadPt    ->Fill(pt,primary->Pt());
		      fhChHadEta   ->Fill(eta,primary->Eta());	
		      fhChHadPhi   ->Fill(phi,primary->Phi());	
		      fhHaVxyz     ->Fill(vx,vy);//,vz);
		      fhHaR        ->Fill(e,r);
		    }
		}//Work with stack also
		
		//Invariant mass
		if (nclusters > 1 ) {
			for(Int_t jclus = iclus + 1 ; jclus < nclusters ; jclus++) {
				AliAODCaloCluster * calo2 =  (AliAODCaloCluster*) (partList->At(jclus));
				//Get cluster kinematics
				calo2->GetMomentum(mom2,v);
				fhIM  ->Fill((mom+mom2).E(),(mom+mom2).M());
				fhAsym->Fill((mom+mom2).E(),TMath::Abs((e-mom2.E())/(e+mom2.E())));
			}// 2nd cluster loop
		}////more than 1 cluster in calorimeter    
		
		//Cells per cluster
		fhNCellsPerCluster->Fill(e,calo->GetNCells());
		
	}// 1st cluster loop
	
	//CaloCells
	AliAODCaloCells * cell = new AliAODCaloCells ;
	if(fCalorimeter == "PHOS") 
		cell = (AliAODCaloCells *) GetPHOSCells();
	else  
		cell = (AliAODCaloCells *) GetEMCALCells();	
	
	if(!cell) {
		printf("AliAnaCalorimeterQA::MakeAnalysisFillHistograms() - ABORT: No CELLS available for analysis");
		abort();
	}
	
	//Some prints
	if(GetDebug() > 0 && cell )
		printf("AliAnaCalorimeterQA::MakeAnalysisFillHistograms() - In ESD %s cell entries %d\n", fCalorimeter.Data(), cell->GetNumberOfCells());    
	
	Int_t ncells = cell->GetNumberOfCells() ;
	fhNCells->Fill(ncells) ;
	
	for (Int_t iCell = 0; iCell < ncells; iCell++) {      
		if(GetDebug() > 2)  printf("AliAnaCalorimeterQA::MakeAnalysisFillHistograms() - Cell : amp %f, absId %d \n", cell->GetAmplitude(iCell), cell->GetCellNumber(iCell));
		fhAmplitude->Fill(cell->GetAmplitude(iCell));
	}
	
	if(GetDebug() > 0) 
		printf("AliAnaCalorimeterQA::MakeAnalysisFillHistograms() - end \n");  
	
	//Monte Carlo
	if(IsDataMC()){
		//Play with the MC stack if available
		for(Int_t i=8 ; i<stack->GetNprimary(); i++){
			primary = stack->Particle(i) ;
			
			//if (!primary->IsPrimary()) continue;
			if (TMath::Abs(primary->Eta()) > 1) continue;
			
			Int_t kf = primary->GetPdgCode();
			//printf("kf %d\n",kf);
			
			Bool_t in = kTRUE;
			primary->Momentum(mom);
			if(IsFidutialCutOn()) in =  GetFidutialCut()->IsInFidutialCut(mom,fCalorimeter) ;
			
			if (kf==22) {
				fhGenGamPt ->Fill(primary->Pt());
				fhGenGamEta->Fill(primary->Eta());
				fhGenGamPhi->Fill(primary->Phi());
				if(in){
					fhGenGamAccE  ->Fill(primary->Energy());
					fhGenGamAccPt ->Fill(primary->Pt());
					fhGenGamAccEta->Fill(primary->Eta());
					fhGenGamAccPhi->Fill(primary->Phi());					
				}
			}
			else if (kf==111) {
				fhGenPi0Pt ->Fill(primary->Pt());
				fhGenPi0Eta->Fill(primary->Eta());
				fhGenPi0Phi->Fill(primary->Phi());
				if(in){
					fhGenPi0AccE  ->Fill(primary->Energy());					
					fhGenPi0AccPt ->Fill(primary->Pt());
					fhGenPi0AccEta->Fill(primary->Eta());
					fhGenPi0AccPhi->Fill(primary->Phi());					
				}
			}
			else if (kf==221) {
				fhGenEtaPt ->Fill(primary->Pt());
				fhGenEtaEta->Fill(primary->Eta());
				fhGenEtaPhi->Fill(primary->Phi());
			}
			else if (kf==223) {
				fhGenOmegaPt ->Fill(primary->Pt());
				fhGenOmegaEta->Fill(primary->Eta());
				fhGenOmegaPhi->Fill(primary->Phi());
			}
			else if (TMath::Abs(kf)==11) {
				fhGenElePt ->Fill(primary->Pt());
				fhGenEleEta->Fill(primary->Eta());
				fhGenElePhi->Fill(primary->Phi());
			}
		} //primary loop
	} //Is data MC
	
}

//________________________________________________________________________
void AliAnaCalorimeterQA::ReadHistograms(TList* outputList)
{
	// Needed when Terminate is executed in distributed environment
	// Refill analysis histograms of this class with corresponding histograms in output list. 
	
	// Histograms of this analsys are kept in the same list as other analysis, recover the position of
	// the first one and then add the next 
	Int_t index = outputList->IndexOf(outputList->FindObject(GetAddedHistogramsStringToName()+"hE"));
	//printf("Calo: %s, index: %d\n",fCalorimeter.Data(),index);
	//Read histograms, must be in the same order as in GetCreateOutputObject.
	fhE      = (TH1F *) outputList->At(index++); 	
	fhPt     = (TH1F *) outputList->At(index++); 
	fhPhi    = (TH1F *) outputList->At(index++); 
	fhEta    = (TH1F *) outputList->At(index++);
	fhEtaPhi = (TH2F *) outputList->At(index++);
	
	fhECharged      = (TH1F *) outputList->At(index++); 	
	fhPtCharged     = (TH1F *) outputList->At(index++); 
	fhPhiCharged    = (TH1F *) outputList->At(index++); 
	fhEtaCharged    = (TH1F *) outputList->At(index++);
	fhEtaPhiCharged = (TH2F *) outputList->At(index++);
	
	fhEChargedNoOut      = (TH1F *) outputList->At(index++); 	
	fhPtChargedNoOut     = (TH1F *) outputList->At(index++); 
	fhPhiChargedNoOut    = (TH1F *) outputList->At(index++); 
	fhEtaChargedNoOut    = (TH1F *) outputList->At(index++);
	fhEtaPhiChargedNoOut = (TH2F *) outputList->At(index++);

	fhIM     = (TH2F *) outputList->At(index++);
	fhAsym   = (TH2F *) outputList->At(index++);
	
	fhNCellsPerCluster = (TH2F *) outputList->At(index++);
	fhNClusters  = (TH1F *) outputList->At(index++); 
	fhNCells     = (TH1F *) outputList->At(index++); 
	fhAmplitude  = (TH1F *) outputList->At(index++); 
	
	if(IsDataMC()){
		fhDeltaE   = (TH1F *) outputList->At(index++); 
		fhDeltaPt  = (TH1F *) outputList->At(index++); 
		fhDeltaPhi = (TH1F *) outputList->At(index++); 
		fhDeltaEta = (TH1F *) outputList->At(index++); 
		
		fhRatioE   = (TH1F *) outputList->At(index++); 
		fhRatioPt  = (TH1F *) outputList->At(index++); 
		fhRatioPhi = (TH1F *) outputList->At(index++); 
		fhRatioEta = (TH1F *) outputList->At(index++); 
		
		fh2E       = (TH2F *) outputList->At(index++); 
		fh2Pt      = (TH2F *) outputList->At(index++); 
		fh2Phi     = (TH2F *) outputList->At(index++); 
		fh2Eta     = (TH2F *) outputList->At(index++); 
		
		fhGamE     = (TH2F *) outputList->At(index++); 
		fhGamPt    = (TH2F *) outputList->At(index++); 
		fhGamPhi   = (TH2F *) outputList->At(index++); 
		fhGamEta   = (TH2F *) outputList->At(index++); 
		
		fhGamDeltaE   = (TH1F *) outputList->At(index++); 
		fhGamDeltaPt  = (TH1F *) outputList->At(index++); 
		fhGamDeltaPhi = (TH1F *) outputList->At(index++); 
		fhGamDeltaEta = (TH1F *) outputList->At(index++); 
		
		fhGamRatioE   = (TH1F *) outputList->At(index++); 
		fhGamRatioPt  = (TH1F *) outputList->At(index++); 
		fhGamRatioPhi = (TH1F *) outputList->At(index++); 
		fhGamRatioEta = (TH1F *) outputList->At(index++); 

		fhGam2E       = (TH2F *) outputList->At(index++); 
		fhGam2Pt      = (TH2F *) outputList->At(index++); 
		fhGam2Phi     = (TH2F *) outputList->At(index++); 
		fhGam2Eta     = (TH2F *) outputList->At(index++); 

		fhPi0E     = (TH2F *) outputList->At(index++); 
		fhPi0Pt    = (TH2F *) outputList->At(index++); 
		fhPi0Phi   = (TH2F *) outputList->At(index++); 
		fhPi0Eta   = (TH2F *) outputList->At(index++); 		
		
		fhEleE     = (TH2F *) outputList->At(index++); 
		fhElePt    = (TH2F *) outputList->At(index++); 
		fhElePhi   = (TH2F *) outputList->At(index++); 
		fhEleEta   = (TH2F *) outputList->At(index++); 		
		
		fhNeHadE     = (TH2F *) outputList->At(index++); 
		fhNeHadPt    = (TH2F *) outputList->At(index++); 
		fhNeHadPhi   = (TH2F *) outputList->At(index++); 
		fhNeHadEta   = (TH2F *) outputList->At(index++); 		
		
		fhChHadE     = (TH2F *) outputList->At(index++); 
		fhChHadPt    = (TH2F *) outputList->At(index++); 
		fhChHadPhi   = (TH2F *) outputList->At(index++); 
		fhChHadEta   = (TH2F *) outputList->At(index++); 				
		
//		fhEMVxyz     = (TH3F *) outputList->At(index++); 
//		fhHaVxyz     = (TH3F *) outputList->At(index++); 
		
		fhEMVxyz     = (TH2F *) outputList->At(index++); 
		fhHaVxyz     = (TH2F *) outputList->At(index++); 
		fhEMR        = (TH2F *) outputList->At(index++); 
		fhHaR        = (TH2F *) outputList->At(index++); 
		
		fhGenGamPt    = (TH1F *) outputList->At(index++); 
		fhGenGamEta   = (TH1F *) outputList->At(index++); 
		fhGenGamPhi   = (TH1F *) outputList->At(index++); 
		
		fhGenPi0Pt    = (TH1F *) outputList->At(index++); 
		fhGenPi0Eta   = (TH1F *) outputList->At(index++); 
		fhGenPi0Phi   = (TH1F *) outputList->At(index++); 
		
		fhGenEtaPt    = (TH1F *) outputList->At(index++); 
		fhGenEtaEta   = (TH1F *) outputList->At(index++); 
		fhGenEtaPhi   = (TH1F *) outputList->At(index++); 
		
		fhGenOmegaPt  = (TH1F *) outputList->At(index++); 
		fhGenOmegaEta = (TH1F *) outputList->At(index++); 
		fhGenOmegaPhi = (TH1F *) outputList->At(index++); 
		
		fhGenElePt    = (TH1F *) outputList->At(index++); 
		fhGenEleEta   = (TH1F *) outputList->At(index++); 
		fhGenElePhi   = (TH1F *) outputList->At(index++); 
		
		fhGenGamAccE   = (TH1F *) outputList->At(index++); 		
		fhGenGamAccPt  = (TH1F *) outputList->At(index++); 
		fhGenGamAccEta = (TH1F *) outputList->At(index++); 
		fhGenGamAccPhi = (TH1F *) outputList->At(index++); 
		
		fhGenPi0AccE   = (TH1F *) outputList->At(index++); 		
		fhGenPi0AccPt  = (TH1F *) outputList->At(index++); 
		fhGenPi0AccEta = (TH1F *) outputList->At(index++); 
		fhGenPi0AccPhi = (TH1F *) outputList->At(index++); 
		
	}//Is data MC	
	
	fh1pOverE =    (TH1F *) outputList->At(index++);
	fh1dR =        (TH1F *) outputList->At(index++);
	fh2MatchdEdx = (TH2F *) outputList->At(index++);
	fh2EledEdx =   (TH2F *) outputList->At(index++);
	
//	for(Int_t i = 0;  i<index ; i++) cout<<outputList->At(i)->GetName()<<endl;
}

//__________________________________________________________________
void  AliAnaCalorimeterQA::Terminate(TList* outputList) 
{
	
	//Do some plots to end
	 if(fStyleMacro!="")gROOT->Macro(fStyleMacro); 
	//Recover histograms from output histograms list, needed for distributed analysis.	
	ReadHistograms(outputList);
	
	//printf(" AliAnaCalorimeterQA::Terminate()  *** %s Report:", GetName()) ; 
	//printf(" AliAnaCalorimeterQA::Terminate()        pt         : %5.3f , RMS : %5.3f \n", fhPt->GetMean(),   fhPt->GetRMS() ) ;

	char name[128];
	char cname[128];
	
	//Reconstructed distributions
	//printf("c\n");
	sprintf(cname,"QA_%s_rec",fCalorimeter.Data());
	TCanvas  * c = new TCanvas(cname, "Reconstructed distributions", 400, 400) ;
	c->Divide(2, 2);
	
	c->cd(1) ; 
	gPad->SetLogy();
	fhE->SetLineColor(4);
	fhE->Draw();
	
	c->cd(2) ; 
	gPad->SetLogy();
	fhPt->SetLineColor(4);
	fhPt->Draw();
	
	c->cd(3) ; 
	fhPhi->SetLineColor(4);
	fhPhi->Draw();
	
	c->cd(4) ; 
	fhEta->SetLineColor(4);
	fhEta->Draw();
	
	sprintf(name,"QA_%s_ReconstructedDistributions.eps",fCalorimeter.Data());
	c->Print(name);
	
	//Reconstructed distributions, matched with tracks
	//printf("c2\n");
	sprintf(cname,"QA_%s_rectrackmatch",fCalorimeter.Data());
	TCanvas  * c2 = new TCanvas(cname, "Reconstructed distributions, matched with tracks", 400, 400) ;
	c2->Divide(2, 2);
	
	c2->cd(1) ; 
	gPad->SetLogy();
	fhECharged->SetLineColor(4);
	fhECharged->Draw();
	
	c2->cd(2) ; 
	gPad->SetLogy();
	fhPtCharged->SetLineColor(4);
	fhPtCharged->Draw();
	
	c2->cd(3) ; 
	fhPhiCharged->SetLineColor(4);
	fhPhiCharged->Draw();
	
	c2->cd(4) ; 
	fhEtaCharged->SetLineColor(4);
	fhEtaCharged->Draw();
	
	sprintf(name,"QA_%s_ReconstructedDistributions_TrackMatched.eps",fCalorimeter.Data());
	c2->Print(name);
	
	TH1F *	hEChargedClone   = (TH1F*)   fhECharged->Clone("EChargedClone");
	TH1F *	hPtChargedClone  = (TH1F*)   fhPtCharged->Clone("PtChargedClone");
	TH1F *	hEtaChargedClone = (TH1F*)   fhEtaCharged->Clone("EtaChargedClone");
	TH1F *	hPhiChargedClone = (TH1F*)   fhPhiCharged->Clone("PhiChargedClone");

	TH1F *	hEChargedClone2   = (TH1F*)   fhECharged->Clone("EChargedClone2");
	TH1F *	hPtChargedClone2  = (TH1F*)   fhPtCharged->Clone("PtChargedClone2");
	TH1F *	hEtaChargedClone2 = (TH1F*)   fhEtaCharged->Clone("EtaChargedClone2");
	TH1F *	hPhiChargedClone2 = (TH1F*)   fhPhiCharged->Clone("PhiChargedClone2");

	//Ratio: reconstructed track matched/ all reconstructed
	//printf("c3\n");
	sprintf(cname,"QA_%s_rectrackmatchrat",fCalorimeter.Data());
	TCanvas  * c3 = new TCanvas(cname, "Ratio: reconstructed track matched/ all reconstructed", 400, 400) ;
	c3->Divide(2, 2);
	
	c3->cd(1) ;
	hEChargedClone->SetYTitle("track matched / all");
	hEChargedClone->Divide(fhE);
	hEChargedClone->Draw();
	
	c3->cd(2) ; 
	hPtChargedClone->SetYTitle("track matched / all");
	hPtChargedClone->Divide(fhPt);
	hPtChargedClone->Draw();
	
	c3->cd(3) ;
	hPhiChargedClone->SetYTitle("track matched / all");
	hPhiChargedClone->Divide(fhPhi);
	hPhiChargedClone->Draw();
	
	c3->cd(4) ; 
	hEtaChargedClone->SetYTitle("track matched / all");
	hEtaChargedClone->Divide(fhEta);
	hEtaChargedClone->Draw();
	
	sprintf(name,"QA_%s_RatioReconstructedMatchedDistributions.eps",fCalorimeter.Data());
	c3->Print(name);
	
	//Ratio: reconstructed track matched (minus no track param) / all
	//printf("c333\n");
	sprintf(cname,"QA_%s_rectrackmatchratout",fCalorimeter.Data());
	TCanvas  * c333 = new TCanvas(cname, "Ratio: reconstructed track matched (with outer track param)/ all", 400, 400) ;
	c333->Divide(2, 2);

	c333->cd(1) ;
	hEChargedClone2->Add(fhEChargedNoOut,-1);
	hEChargedClone2->SetYTitle("track matched / all");
	hEChargedClone2->Divide(fhE);
	hEChargedClone2->Draw();
	
	c333->cd(2) ; 
	hPtChargedClone2->Add(fhPtChargedNoOut,-1);
	hPtChargedClone2->SetYTitle("track matched / all");
	hPtChargedClone2->Divide(fhPt);
	hPtChargedClone2->Draw();
	
	c333->cd(3) ;
	hPhiChargedClone2->Add(fhPhiChargedNoOut,-1);
	hPhiChargedClone2->SetYTitle("track matched / all");
	hPhiChargedClone2->Divide(fhPhi);
	hPhiChargedClone2->Draw();
	
	c333->cd(4) ; 
	hEtaChargedClone2->Add(fhEtaChargedNoOut,-1);
	hEtaChargedClone2->SetYTitle("track matched / all");
	hEtaChargedClone2->Divide(fhEta);
	hEtaChargedClone2->Draw();
	
	sprintf(name,"QA_%s_RatioReconstructedMatchedDistributionsOuter.eps",fCalorimeter.Data());
	c333->Print(name);

	//Reconstructed distributions, matched with tracks
	//printf("c2\n");
	sprintf(cname,"QA_%s_rectrackmatch_noout",fCalorimeter.Data());
	TCanvas  * c22 = new TCanvas(cname, "Reconstructed distributions, matched with tracks, no outer track param", 400, 400) ;
	c22->Divide(2, 2);
	
	c22->cd(1) ; 
	gPad->SetLogy();
	fhEChargedNoOut->SetLineColor(4);
	fhEChargedNoOut->Draw();
	
	c22->cd(2) ; 
	gPad->SetLogy();
	fhPtChargedNoOut->SetLineColor(4);
	fhPtChargedNoOut->Draw();
	
	c22->cd(3) ; 
	fhPhiChargedNoOut->SetLineColor(4);
	fhPhiChargedNoOut->Draw();
	
	c22->cd(4) ; 
	fhEtaChargedNoOut->SetLineColor(4);
	fhEtaChargedNoOut->Draw();
	
	sprintf(name,"QA_%s_ReconstructedDistributions_TrackMatched_NoOutParam.eps",fCalorimeter.Data());
	c22->Print(name);
	
	//Ratio: reconstructed track matched/ all reconstructed
	//printf("c33\n");
	TH1F *	hEChargedNoOutClone   = (TH1F*)   fhEChargedNoOut->Clone("EChargedNoOutClone");
	TH1F *	hPtChargedNoOutClone  = (TH1F*)   fhPtChargedNoOut->Clone("PtChargedNoOutClone");
	TH1F *	hEtaChargedNoOutClone = (TH1F*)   fhEtaChargedNoOut->Clone("EtaChargedNoOutClone");
	TH1F *	hPhiChargedNoOutClone = (TH1F*)   fhPhiChargedNoOut->Clone("PhiChargedNoOutClone");	
	
	sprintf(cname,"QA_%s_rectrackmatchratnoout",fCalorimeter.Data());
	TCanvas  * c33 = new TCanvas(cname, "Ratio: reconstructed track matched/ all reconstructed", 400, 400) ;
	c33->Divide(2, 2);
	
	c33->cd(1) ;
	hEChargedNoOutClone->SetYTitle("track matched no out/ all matched");
	hEChargedNoOutClone->Divide(fhECharged);
	hEChargedNoOutClone->Draw();
	
	c33->cd(2) ; 
	hPtChargedNoOutClone->SetYTitle("track matched no out / all matched");
	hPtChargedNoOutClone->Divide(fhPtCharged);
	hPtChargedNoOutClone->Draw();
	
	c33->cd(3) ;
	hPhiChargedNoOutClone->SetYTitle("track matched no out/ all matched");
	hPhiChargedNoOutClone->Divide(fhPhiCharged);
	hPhiChargedNoOutClone->Draw();
	
	c33->cd(4) ; 
	hEtaChargedNoOutClone->SetYTitle("track matched no out/ all matched");
	hEtaChargedNoOutClone->Divide(fhEtaCharged);
	hEtaChargedNoOutClone->Draw();
	
	sprintf(name,"QA_%s_RatioMatchedDistributionsAllToNoOut.eps",fCalorimeter.Data());
	c33->Print(name);

	//eta vs phi
	//printf("c4\n");
	sprintf(cname,"QA_%s_etavsphi",fCalorimeter.Data());
	TCanvas  * c4 = new TCanvas(cname, "reconstructed #eta vs #phi", 200, 600) ;
	c4->Divide(1, 3);
		
	c4->cd(1) ;
	fhEtaPhi->Draw("cont");
	
	c4->cd(2) ; 
	fhEtaPhiCharged->Draw("cont");
	
	c4->cd(3) ; 
	fhEtaPhiChargedNoOut->Draw("cont");

	sprintf(name,"QA_%s_ReconstructedEtaVsPhi.eps",fCalorimeter.Data());
	c4->Print(name);
	
	//Invariant mass
	Int_t binmin = -1;
	Int_t binmax = -1;

	if(fhIM->GetEntries() > 1){
		Int_t nebins  = fhIM->GetNbinsX();
		Int_t emax = (Int_t) fhIM->GetXaxis()->GetXmax();
		Int_t emin = (Int_t) fhIM->GetXaxis()->GetXmin();
		if (emin != 0 ) printf("emin != 0 \n");
		//printf("IM: nBinsX %d, emin %2.2f, emax %2.2f\n",nebins,emin,emax);
		
		sprintf(cname,"QA_%s_IM",fCalorimeter.Data());
	//	printf("c5\n");
		TCanvas  * c5 = new TCanvas(cname, "Invariant mass", 400, 400) ;
		c5->Divide(2, 2);
		
		c5->cd(1) ; 
		fhIM->SetLineColor(4);
		fhIM->Draw();
		
		c5->cd(2) ; 
		binmin = 0;
		binmax =  (Int_t) (5-emin)*nebins/emax;
		TH1D *pyim5 = fhIM->ProjectionY("pyim5",binmin,binmax);
		pyim5->SetTitle("E_{pair} < 5 GeV");
		pyim5->SetLineColor(4);
		pyim5->Draw();
		
		c5->cd(3) ; 
		binmin =  (Int_t) (5-emin)*nebins/emax;
		binmax =  (Int_t) (10-emin)*nebins/emax;
		TH1D *pyim510 = fhIM->ProjectionY("pyim5_10",binmin,binmax);
		pyim510->SetTitle("5 < E_{pair} < 10 GeV");
		pyim510->SetLineColor(4);
		pyim510->Draw();
		
		c5->cd(4) ;
		binmin =  (Int_t) (10-emin)*nebins/emax;
		binmax = -1;
		TH1D *pyim10 = fhIM->ProjectionY("pyim10",binmin,binmax);
		pyim10->SetTitle("E_{pair} > 10 GeV");
		pyim10->SetLineColor(4);
		pyim10->Draw();
		
		sprintf(name,"QA_%s_InvariantMass.eps",fCalorimeter.Data());
		c5->Print(name);
	}
	
	//Asymmetry
	if(fhAsym->GetEntries() > 1){
		Int_t nebins  = fhAsym->GetNbinsX();
		Int_t emax = (Int_t) fhAsym->GetXaxis()->GetXmax();
		Int_t emin = (Int_t) fhAsym->GetXaxis()->GetXmin();
		if (emin != 0 ) printf("emin != 0 \n");
		//printf("Asym: nBinsX %d, emin %2.2f, emax %2.2f\n",nebins,emin,emax);
		
		sprintf(cname,"QA_%s_Asym",fCalorimeter.Data());
		//	printf("c5\n");
		TCanvas  * c5b = new TCanvas(cname, "Asymmetry", 400, 400) ;
		c5b->Divide(2, 2);
		
		c5b->cd(1) ; 
		fhAsym->SetLineColor(4);
		fhAsym->Draw();
		
		c5b->cd(2) ; 
		binmin = 0;
		binmax = (Int_t) (5-emin)*nebins/emax;
		TH1D *pyAsym5 = fhAsym->ProjectionY("pyAsym5",binmin,binmax);
		pyAsym5->SetTitle("E_{pair} < 5 GeV");
		pyAsym5->SetLineColor(4);
		pyAsym5->Draw();
		
		c5b->cd(3) ; 
		binmin = (Int_t) (5-emin)*nebins/emax;
		binmax = (Int_t) (10-emin)*nebins/emax;
		TH1D *pyAsym510 = fhAsym->ProjectionY("pyAsym5_10",binmin,binmax);
		pyAsym510->SetTitle("5 < E_{pair} < 10 GeV");
		pyAsym510->SetLineColor(4);
		pyAsym510->Draw();
		
		c5b->cd(4) ;
		binmin = (Int_t) (10-emin)*nebins/emax;
		binmax = -1;
		TH1D *pyAsym10 = fhAsym->ProjectionY("pyAsym10",binmin,binmax);
		pyAsym10->SetTitle("E_{pair} > 10 GeV");
		pyAsym10->SetLineColor(4);
		pyAsym10->Draw();
		
		sprintf(name,"QA_%s_Asymmetry.eps",fCalorimeter.Data());
		c5b->Print(name);
	}
	
	
	//Reconstructed vs MC distributions
	//printf("c6\n");
	sprintf(cname,"QA_%s_recvsmc",fCalorimeter.Data());
	TCanvas  * c6 = new TCanvas(cname, "Reconstructed vs MC distributions", 400, 400) ;
	c6->Divide(2, 2);
	
	c6->cd(1) ; 
	fh2E->SetLineColor(4);
	fh2E->Draw();
	
	c6->cd(2) ; 
	fh2Pt->SetLineColor(4);
	fh2Pt->Draw();
	
	c6->cd(3) ; 
	fh2Phi->SetLineColor(4);
	fh2Phi->Draw();
	
	c6->cd(4) ; 
	fh2Eta->SetLineColor(4);
	fh2Eta->Draw();
	
	sprintf(name,"QA_%s_ReconstructedVSMCDistributions.eps",fCalorimeter.Data());
	c6->Print(name);	
	
	
	//Reconstructed vs MC distributions
	//printf("c6\n");
	sprintf(cname,"QA_%s_gamrecvsmc",fCalorimeter.Data());
	TCanvas  * c6Gam = new TCanvas(cname, "Reconstructed vs MC distributions", 400, 400) ;
	c6Gam->Divide(2, 2);
	
	c6Gam->cd(1) ; 
	fhGam2E->SetLineColor(4);
	fhGam2E->Draw();
	
	c6Gam->cd(2) ; 
	fhGam2Pt->SetLineColor(4);
	fhGam2Pt->Draw();
	
	c6Gam->cd(3) ; 
	fhGam2Phi->SetLineColor(4);
	fhGam2Phi->Draw();
	
	c6Gam->cd(4) ; 
	fhGam2Eta->SetLineColor(4);
	fhGam2Eta->Draw();
	
	sprintf(name,"QA_%s_GammaReconstructedVSMCDistributions.eps",fCalorimeter.Data());
	c6->Print(name);	

	//Generated - reconstructed  
	//printf("c7\n");
	sprintf(cname,"QA_%s_diffgenrec",fCalorimeter.Data());
	TCanvas  * c7 = new TCanvas(cname, "generated - reconstructed", 400, 400) ;
	c7->Divide(2, 2);
	
	c7->cd(1) ; 
	gPad->SetLogy();
	fhGamDeltaE->SetLineColor(4);
	fhDeltaE->Draw();
	fhGamDeltaE->Draw("same");

	TLegend pLegendd(0.65,0.55,0.9,0.8);
	pLegendd.SetTextSize(0.06);
	pLegendd.AddEntry(fhDeltaE,"all","L");
	pLegendd.AddEntry(fhGamDeltaE,"from  #gamma","L");
	pLegendd.SetFillColor(10);
	pLegendd.SetBorderSize(1);
	pLegendd.Draw();

	c7->cd(2) ; 
	gPad->SetLogy();
	fhGamDeltaPt->SetLineColor(4);
	fhDeltaPt->Draw();
       	fhGamDeltaPt->Draw("same");

	c7->cd(3) ; 
	fhGamDeltaPhi->SetLineColor(4);
	fhDeltaPhi->Draw();
       	fhGamDeltaPhi->Draw("same");
	
	c7->cd(4) ; 
	fhGamDeltaEta->SetLineColor(4);
	fhDeltaEta->Draw();
	fhGamDeltaEta->Draw("same");

	sprintf(name,"QA_%s_DiffGeneratedReconstructed.eps",fCalorimeter.Data());
	c7->Print(name);
	
	
	// Reconstructed / Generated 
	//printf("c8\n");
	sprintf(cname,"QA_%s_ratiorecgen",fCalorimeter.Data());
	TCanvas  * c8 = new TCanvas(cname, " reconstructed / generated", 400, 400) ;
	c8->Divide(2, 2);
	
	c8->cd(1) ; 
	gPad->SetLogy();
	fhGamRatioE->SetLineColor(4);
	fhRatioE->Draw();
	fhGamRatioE->Draw("same");

	TLegend pLegendr(0.65,0.55,0.9,0.8);
	pLegendr.SetTextSize(0.06);
	pLegendr.AddEntry(fhRatioE,"all","L");
	pLegendr.AddEntry(fhGamRatioE,"from  #gamma","L");
	pLegendr.SetFillColor(10);
	pLegendr.SetBorderSize(1);
	pLegendr.Draw();

	c8->cd(2) ; 
	gPad->SetLogy();
	fhGamRatioPt->SetLineColor(4);
	fhRatioPt->Draw();
       	fhGamRatioPt->Draw("same");

	c8->cd(3) ; 
	fhGamRatioPhi->SetLineColor(4);
	fhRatioPhi->Draw();
       	fhGamRatioPhi->Draw("same");

	c8->cd(4) ; 
	fhGamRatioEta->SetLineColor(4);
	fhRatioEta->Draw();
       	fhGamRatioEta->Draw("same");

	sprintf(name,"QA_%s_ReconstructedDivGenerated.eps",fCalorimeter.Data());
	c8->Print(name);
	
	// CaloClusters CaloCells
	//printf("c9\n");
	sprintf(cname,"QA_%s_nclustercellsamp",fCalorimeter.Data());
	TCanvas  * c9 = new TCanvas(cname, " CaloClusters and CaloCells", 400, 400) ;
	c9->Divide(2, 2);
	
	c9->cd(1) ; 
	gPad->SetLogy();
	gPad->SetLogx();
	fhNClusters->SetLineColor(4);
	fhNClusters->Draw();
	
	c9->cd(2) ; 
	gPad->SetLogy();
	gPad->SetLogx();
	fhNCells->SetLineColor(4);
	fhNCells->Draw();
	
	c9->cd(3) ; 
	gPad->SetLogy();
	gPad->SetLogx();
	fhNCellsPerCluster->SetLineColor(4);
	fhNCellsPerCluster->Draw();
	
	c9->cd(4) ; 
	gPad->SetLogy();
	//gPad->SetLogx();
	fhAmplitude->SetLineColor(4);
	fhAmplitude->Draw();
	
	sprintf(name,"QA_%s_CaloClustersAndCaloCells.eps",fCalorimeter.Data());
	c9->Print(name);
	
	
	//MC
	
	//Generated distributions
	//printf("c1\n");
	sprintf(cname,"QA_%s_gen",fCalorimeter.Data());
	TCanvas  * c10 = new TCanvas(cname, "Generated distributions", 600, 200) ;
	c10->Divide(3, 1);
	
	c10->cd(1) ; 
	gPad->SetLogy();
	TH1F * haxispt  = (TH1F*) fhGenPi0Pt->Clone("axispt");  
	haxispt->SetTitle("Generated Particles p_{T}, |#eta| < 1");
	fhGenPi0Pt->SetLineColor(1);
	fhGenGamPt->SetLineColor(4);
	fhGenEtaPt->SetLineColor(2);
	fhGenOmegaPt->SetLineColor(7);
	fhGenElePt->SetLineColor(6);
	
	//Select the maximum of the histogram to show all lines.
	if(fhGenPi0Pt->GetMaximum() >= fhGenGamPt->GetMaximum() && fhGenPi0Pt->GetMaximum() >= fhGenEtaPt->GetMaximum() && 
	   fhGenPi0Pt->GetMaximum() >= fhGenOmegaPt->GetMaximum() && fhGenPi0Pt->GetMaximum() >= fhGenElePt->GetMaximum())
		haxispt->SetMaximum(fhGenPi0Pt->GetMaximum());
	else if(fhGenGamPt->GetMaximum() >= fhGenPi0Pt->GetMaximum() && fhGenGamPt->GetMaximum() >= fhGenEtaPt->GetMaximum() && 
	   fhGenGamPt->GetMaximum() >= fhGenOmegaPt->GetMaximum() && fhGenGamPt->GetMaximum() >= fhGenElePt->GetMaximum())
		haxispt->SetMaximum(fhGenGamPt->GetMaximum());
	else if(fhGenEtaPt->GetMaximum() >= fhGenPi0Pt->GetMaximum() && fhGenEtaPt->GetMaximum() >= fhGenGamPt->GetMaximum() && 
			fhGenEtaPt->GetMaximum() >= fhGenOmegaPt->GetMaximum() && fhGenEtaPt->GetMaximum() >= fhGenElePt->GetMaximum())
		haxispt->SetMaximum(fhGenEtaPt->GetMaximum());	
	else if(fhGenOmegaPt->GetMaximum() >= fhGenPi0Pt->GetMaximum() && fhGenOmegaPt->GetMaximum() >= fhGenEtaPt->GetMaximum() && 
			fhGenOmegaPt->GetMaximum() >= fhGenGamPt->GetMaximum() && fhGenOmegaPt->GetMaximum() >= fhGenElePt->GetMaximum())
		haxispt->SetMaximum(fhGenOmegaPt->GetMaximum());
	else if(fhGenElePt->GetMaximum() >= fhGenPi0Pt->GetMaximum() && fhGenElePt->GetMaximum() >= fhGenEtaPt->GetMaximum() && 
			fhGenElePt->GetMaximum() >= fhGenOmegaPt->GetMaximum() && fhGenElePt->GetMaximum() >= fhGenGamPt->GetMaximum())
		haxispt->SetMaximum(fhGenElePt->GetMaximum());
	
	haxispt->Draw("axis");
	fhGenPi0Pt->Draw("same");
	fhGenGamPt->Draw("same");
	fhGenEtaPt->Draw("same");
	fhGenOmegaPt->Draw("same");
	fhGenElePt->Draw("same");
	
	TLegend pLegend(0.75,0.45,0.9,0.8);
	pLegend.SetTextSize(0.06);
	pLegend.AddEntry(fhGenPi0Pt,"#pi^{0}","L");
	pLegend.AddEntry(fhGenGamPt,"#gamma","L");
	pLegend.AddEntry(fhGenEtaPt,"#eta","L");
	pLegend.AddEntry(fhGenOmegaPt,"#omega","L");
	pLegend.AddEntry(fhGenElePt,"e^{#pm}","L");
	pLegend.SetFillColor(10);
	pLegend.SetBorderSize(1);
	pLegend.Draw();
	
	c10->cd(2) ;
	gPad->SetLogy();
	TH1F * haxiseta  = (TH1F*) fhGenPi0Eta->Clone("axiseta");  
	haxiseta->SetTitle("Generated Particles #eta, |#eta| < 1");
	fhGenPi0Eta->SetLineColor(1);
	fhGenGamEta->SetLineColor(4);
	fhGenEtaEta->SetLineColor(2);
	fhGenOmegaEta->SetLineColor(7);
	fhGenEleEta->SetLineColor(6);
	//Select the maximum of the histogram to show all lines.
	if(fhGenPi0Eta->GetMaximum() >= fhGenGamEta->GetMaximum() && fhGenPi0Eta->GetMaximum() >= fhGenEtaEta->GetMaximum() && 
	   fhGenPi0Eta->GetMaximum() >= fhGenOmegaEta->GetMaximum() && fhGenPi0Eta->GetMaximum() >= fhGenEleEta->GetMaximum())
		haxiseta->SetMaximum(fhGenPi0Eta->GetMaximum());
	else if(fhGenGamEta->GetMaximum() >= fhGenPi0Eta->GetMaximum() && fhGenGamEta->GetMaximum() >= fhGenEtaEta->GetMaximum() && 
			fhGenGamEta->GetMaximum() >= fhGenOmegaEta->GetMaximum() && fhGenGamEta->GetMaximum() >= fhGenEleEta->GetMaximum())
		haxiseta->SetMaximum(fhGenGamEta->GetMaximum());
	else if(fhGenEtaEta->GetMaximum() >= fhGenPi0Eta->GetMaximum() && fhGenEtaEta->GetMaximum() >= fhGenGamEta->GetMaximum() && 
			fhGenEtaEta->GetMaximum() >= fhGenOmegaEta->GetMaximum() && fhGenEtaEta->GetMaximum() >= fhGenEleEta->GetMaximum())
		haxiseta->SetMaximum(fhGenEtaEta->GetMaximum());	
	else if(fhGenOmegaEta->GetMaximum() >= fhGenPi0Eta->GetMaximum() && fhGenOmegaEta->GetMaximum() >= fhGenEtaEta->GetMaximum() && 
			fhGenOmegaEta->GetMaximum() >= fhGenGamEta->GetMaximum() && fhGenOmegaEta->GetMaximum() >= fhGenEleEta->GetMaximum())
		haxiseta->SetMaximum(fhGenOmegaEta->GetMaximum());
	else if(fhGenEleEta->GetMaximum() >= fhGenPi0Eta->GetMaximum() && fhGenEleEta->GetMaximum() >= fhGenEtaEta->GetMaximum() && 
			fhGenEleEta->GetMaximum() >= fhGenOmegaEta->GetMaximum() && fhGenEleEta->GetMaximum() >= fhGenGamEta->GetMaximum())
		haxiseta->SetMaximum(fhGenEleEta->GetMaximum());
	
	haxiseta->Draw("axis");
	fhGenPi0Eta->Draw("same");
	fhGenGamEta->Draw("same");
	fhGenEtaEta->Draw("same");
	fhGenOmegaEta->Draw("same");
	fhGenEleEta->Draw("same");
	
	
	c10->cd(3) ; 
	gPad->SetLogy();
	TH1F * haxisphi  = (TH1F*) fhGenPi0Phi->Clone("axisphi");  
	haxisphi->SetTitle("Generated Particles #phi, |#eta| < 1");
	fhGenPi0Phi->SetLineColor(1);
	fhGenGamPhi->SetLineColor(4);
	fhGenEtaPhi->SetLineColor(2);
	fhGenOmegaPhi->SetLineColor(7);
	fhGenElePhi->SetLineColor(6);
	//Select the maximum of the histogram to show all lines.
	if(fhGenPi0Phi->GetMaximum() >= fhGenGamPhi->GetMaximum() && fhGenPi0Phi->GetMaximum() >= fhGenEtaPhi->GetMaximum() && 
	   fhGenPi0Phi->GetMaximum() >= fhGenOmegaPhi->GetMaximum() && fhGenPi0Phi->GetMaximum() >= fhGenElePhi->GetMaximum())
		haxisphi->SetMaximum(fhGenPi0Phi->GetMaximum());
	else if(fhGenGamPhi->GetMaximum() >= fhGenPi0Phi->GetMaximum() && fhGenGamPhi->GetMaximum() >= fhGenEtaPhi->GetMaximum() && 
			fhGenGamPhi->GetMaximum() >= fhGenOmegaPhi->GetMaximum() && fhGenGamPhi->GetMaximum() >= fhGenElePhi->GetMaximum())
		haxisphi->SetMaximum(fhGenGamPhi->GetMaximum());
	else if(fhGenEtaPhi->GetMaximum() >= fhGenPi0Phi->GetMaximum() && fhGenEtaPhi->GetMaximum() >= fhGenGamPhi->GetMaximum() && 
			fhGenEtaPhi->GetMaximum() >= fhGenOmegaPhi->GetMaximum() && fhGenEtaPhi->GetMaximum() >= fhGenElePhi->GetMaximum())
		haxisphi->SetMaximum(fhGenEtaPhi->GetMaximum());	
	else if(fhGenOmegaPhi->GetMaximum() >= fhGenPi0Phi->GetMaximum() && fhGenOmegaPhi->GetMaximum() >= fhGenEtaPhi->GetMaximum() && 
			fhGenOmegaPhi->GetMaximum() >= fhGenGamPhi->GetMaximum() && fhGenOmegaPhi->GetMaximum() >= fhGenElePhi->GetMaximum())
		haxisphi->SetMaximum(fhGenOmegaPhi->GetMaximum());
	else if(fhGenElePhi->GetMaximum() >= fhGenPi0Phi->GetMaximum() && fhGenElePhi->GetMaximum() >= fhGenEtaPhi->GetMaximum() && 
			fhGenElePhi->GetMaximum() >= fhGenOmegaPhi->GetMaximum() && fhGenElePhi->GetMaximum() >= fhGenGamPhi->GetMaximum())
		haxisphi->SetMaximum(fhGenElePhi->GetMaximum());
	
	haxisphi->Draw("axis");
	fhGenPi0Phi->Draw("same");
	fhGenGamPhi->Draw("same");
	fhGenEtaPhi->Draw("same");
	fhGenOmegaPhi->Draw("same");
	fhGenElePhi->Draw("same");
	
	sprintf(name,"QA_%s_GeneratedDistributions.eps",fCalorimeter.Data());
	c10->Print(name);
	
	
	//Reconstructed clusters depending on its original particle.
	//printf("c1\n");
	sprintf(cname,"QA_%s_recgenid",fCalorimeter.Data());
	TCanvas  * c11 = new TCanvas(cname, "Reconstructed particles, function of their original particle ID", 400, 400) ;
	c11->Divide(2, 2);
	
	
	c11->cd(1) ; 
	gPad->SetLogy();
	TH1F * hGamE   = (TH1F*) fhGamE->ProjectionY("hGamE",-1,-1);
	TH1F * hPi0E   = (TH1F*) fhPi0E->ProjectionY("hPi0E",-1,-1);
	TH1F * hEleE   = (TH1F*) fhEleE->ProjectionY("hEleE",-1,-1);
	TH1F * hNeHadE = (TH1F*) fhNeHadE->ProjectionY("hNeHadE",-1,-1);
	TH1F * hChHadE = (TH1F*) fhNeHadE->ProjectionY("hChHadE",-1,-1);
	TH1F * haxisE  = (TH1F*) hPi0E->Clone("axisE");  
	haxisE->SetTitle("Reconstructed particles E, function of their original particle ID");
	hPi0E->SetLineColor(1);
	hGamE->SetLineColor(4);
	hNeHadE->SetLineColor(2);
	hChHadE->SetLineColor(7);
	hEleE->SetLineColor(6);
	
	//Select the maximum of the histogram to show all lines.
	if(hPi0E->GetMaximum() >= hGamE->GetMaximum() && hPi0E->GetMaximum() >= hNeHadE->GetMaximum() && 
	   hPi0E->GetMaximum() >= hChHadE->GetMaximum() && hPi0E->GetMaximum() >= hEleE->GetMaximum())
		haxisE->SetMaximum(hPi0E->GetMaximum());
	else if(hGamE->GetMaximum() >= hPi0E->GetMaximum() && hGamE->GetMaximum() >= hNeHadE->GetMaximum() && 
			hGamE->GetMaximum() >= hChHadE->GetMaximum() && hGamE->GetMaximum() >= hEleE->GetMaximum())
		haxisE->SetMaximum(hGamE->GetMaximum());
	else if(hNeHadE->GetMaximum() >= hPi0E->GetMaximum() && hNeHadE->GetMaximum() >= hGamE->GetMaximum() && 
			hNeHadE->GetMaximum() >= hChHadE->GetMaximum() && hNeHadE->GetMaximum() >= hEleE->GetMaximum())
		haxisE->SetMaximum(hNeHadE->GetMaximum());	
	else if(hChHadE->GetMaximum() >= hPi0E->GetMaximum() && hChHadE->GetMaximum() >= hNeHadE->GetMaximum() && 
			hChHadE->GetMaximum() >= hGamE->GetMaximum() && hChHadE->GetMaximum() >= hEleE->GetMaximum())
		haxisE->SetMaximum(hChHadE->GetMaximum());
	else if(hEleE->GetMaximum() >= hPi0E->GetMaximum() && hEleE->GetMaximum() >= hNeHadE->GetMaximum() && 
			hEleE->GetMaximum() >= hChHadE->GetMaximum() && hEleE->GetMaximum() >= hGamE->GetMaximum())
		haxisE->SetMaximum(hEleE->GetMaximum());
	
	haxisE->Draw("axis");
	hPi0E->Draw("same");
	hGamE->Draw("same");
	hNeHadE->Draw("same");
	hChHadE->Draw("same");
	hEleE->Draw("same");
	
	TLegend pLegend2(0.75,0.45,0.9,0.8);
	pLegend2.SetTextSize(0.06);
	pLegend2.AddEntry(hPi0E,"#pi^{0}","L");
	pLegend2.AddEntry(hGamE,"#gamma","L");
	pLegend2.AddEntry(hEleE,"e^{#pm}","L");
	pLegend2.AddEntry(hChHadE,"h^{#pm}","L");
	pLegend2.AddEntry(hNeHadE,"h^{0}","L");
	pLegend2.SetFillColor(10);
	pLegend2.SetBorderSize(1);
	pLegend2.Draw();
	
	
	c11->cd(2) ; 
	gPad->SetLogy();
	//printf("%s, %s, %s, %s, %s\n",fhGamPt->GetName(),fhPi0Pt->GetName(),fhElePt->GetName(),fhNeHadPt->GetName(), fhChHadPt->GetName());
	TH1F * hGamPt   = (TH1F*) fhGamPt->ProjectionY("hGamPt",-1,-1);
	TH1F * hPi0Pt   = (TH1F*) fhPi0Pt->ProjectionY("hPi0Pt",-1,-1);
	TH1F * hElePt   = (TH1F*) fhElePt->ProjectionY("hElePt",-1,-1);
	TH1F * hNeHadPt = (TH1F*) fhNeHadPt->ProjectionY("hNeHadPt",-1,-1);
	TH1F * hChHadPt = (TH1F*) fhNeHadPt->ProjectionY("hChHadPt",-1,-1);
	haxispt  = (TH1F*) hPi0Pt->Clone("axispt");  
	haxispt->SetTitle("Reconstructed particles p_{T}, function of their original particle ID");
	hPi0Pt->SetLineColor(1);
	hGamPt->SetLineColor(4);
	hNeHadPt->SetLineColor(2);
	hChHadPt->SetLineColor(7);
	hElePt->SetLineColor(6);
	
	//Select the maximum of the histogram to show all lines.
	if(hPi0Pt->GetMaximum() >= hGamPt->GetMaximum() && hPi0Pt->GetMaximum() >= hNeHadPt->GetMaximum() && 
	   hPi0Pt->GetMaximum() >= hChHadPt->GetMaximum() && hPi0Pt->GetMaximum() >= hElePt->GetMaximum())
		haxispt->SetMaximum(hPi0Pt->GetMaximum());
	else if(hGamPt->GetMaximum() >= hPi0Pt->GetMaximum() && hGamPt->GetMaximum() >= hNeHadPt->GetMaximum() && 
			hGamPt->GetMaximum() >= hChHadPt->GetMaximum() && hGamPt->GetMaximum() >= hElePt->GetMaximum())
		haxispt->SetMaximum(hGamPt->GetMaximum());
	else if(hNeHadPt->GetMaximum() >= hPi0Pt->GetMaximum() && hNeHadPt->GetMaximum() >= hGamPt->GetMaximum() && 
			hNeHadPt->GetMaximum() >= hChHadPt->GetMaximum() && hNeHadPt->GetMaximum() >= hElePt->GetMaximum())
		haxispt->SetMaximum(hNeHadPt->GetMaximum());	
	else if(hChHadPt->GetMaximum() >= hPi0Pt->GetMaximum() && hChHadPt->GetMaximum() >= hNeHadPt->GetMaximum() && 
			hChHadPt->GetMaximum() >= hGamPt->GetMaximum() && hChHadPt->GetMaximum() >= hElePt->GetMaximum())
		haxispt->SetMaximum(hChHadPt->GetMaximum());
	else if(hElePt->GetMaximum() >= hPi0Pt->GetMaximum() && hElePt->GetMaximum() >= hNeHadPt->GetMaximum() && 
			hElePt->GetMaximum() >= hChHadPt->GetMaximum() && hElePt->GetMaximum() >= hGamPt->GetMaximum())
		haxispt->SetMaximum(hElePt->GetMaximum());
	
	haxispt->Draw("axis");
	hPi0Pt->Draw("same");
	hGamPt->Draw("same");
	hNeHadPt->Draw("same");
	hChHadPt->Draw("same");
	hElePt->Draw("same");
	
	
	c11->cd(3) ;
	gPad->SetLogy();

	TH1F * hGamEta   = (TH1F*) fhGamEta->ProjectionY("hGamEta",-1,-1);
	TH1F * hPi0Eta   = (TH1F*) fhPi0Eta->ProjectionY("hPi0Eta",-1,-1);
	TH1F * hEleEta   = (TH1F*) fhEleEta->ProjectionY("hEleEta",-1,-1);
	TH1F * hNeHadEta = (TH1F*) fhNeHadEta->ProjectionY("hNeHadEta",-1,-1);
	TH1F * hChHadEta = (TH1F*) fhNeHadEta->ProjectionY("hChHadEta",-1,-1);
	haxiseta  = (TH1F*) hPi0Eta->Clone("axiseta");  
	haxiseta->SetTitle("Reconstructed particles #eta, function of their original particle ID");
	hPi0Eta->SetLineColor(1);
	hGamEta->SetLineColor(4);
	hNeHadEta->SetLineColor(2);
	hChHadEta->SetLineColor(7);
	hEleEta->SetLineColor(6);
	//Select the maximum of the histogram to show all lines.
	if(hPi0Eta->GetMaximum() >= hGamEta->GetMaximum() && hPi0Eta->GetMaximum() >= hNeHadEta->GetMaximum() && 
	   hPi0Eta->GetMaximum() >= hChHadEta->GetMaximum() && hPi0Eta->GetMaximum() >= hEleEta->GetMaximum())
		haxiseta->SetMaximum(hPi0Eta->GetMaximum());
	else if(hGamEta->GetMaximum() >= hPi0Eta->GetMaximum() && hGamEta->GetMaximum() >= hNeHadEta->GetMaximum() && 
			hGamEta->GetMaximum() >= hChHadEta->GetMaximum() && hGamEta->GetMaximum() >= hEleEta->GetMaximum())
		haxiseta->SetMaximum(hGamEta->GetMaximum());
	else if(hNeHadEta->GetMaximum() >= hPi0Eta->GetMaximum() && hNeHadEta->GetMaximum() >= hGamEta->GetMaximum() && 
			hNeHadEta->GetMaximum() >= hChHadEta->GetMaximum() && hNeHadEta->GetMaximum() >= hEleEta->GetMaximum())
		haxiseta->SetMaximum(hNeHadEta->GetMaximum());	
	else if(hChHadEta->GetMaximum() >= hPi0Eta->GetMaximum() && hChHadEta->GetMaximum() >= hNeHadEta->GetMaximum() && 
			hChHadEta->GetMaximum() >= hGamEta->GetMaximum() && hChHadEta->GetMaximum() >= hEleEta->GetMaximum())
		haxiseta->SetMaximum(hChHadEta->GetMaximum());
	else if(hEleEta->GetMaximum() >= hPi0Eta->GetMaximum() && hEleEta->GetMaximum() >= hNeHadEta->GetMaximum() && 
			hEleEta->GetMaximum() >= hChHadEta->GetMaximum() && hEleEta->GetMaximum() >= hGamEta->GetMaximum())
		haxiseta->SetMaximum(hEleEta->GetMaximum());
	
	haxiseta->Draw("axis");
	hPi0Eta->Draw("same");
	hGamEta->Draw("same");
	hNeHadEta->Draw("same");
	hChHadEta->Draw("same");
	hEleEta->Draw("same");
	
	
	c11->cd(4) ; 
	gPad->SetLogy();
	TH1F * hGamPhi   = (TH1F*) fhGamPhi->ProjectionY("hGamPhi",-1,-1);
	TH1F * hPi0Phi   = (TH1F*) fhPi0Phi->ProjectionY("hPi0Phi",-1,-1);
	TH1F * hElePhi   = (TH1F*) fhElePhi->ProjectionY("hElePhi",-1,-1);
	TH1F * hNeHadPhi = (TH1F*) fhNeHadPhi->ProjectionY("hNeHadPhi",-1,-1);
	TH1F * hChHadPhi = (TH1F*) fhNeHadPhi->ProjectionY("hChHadPhi",-1,-1);
	haxisphi  = (TH1F*) hPi0Phi->Clone("axisphi");  
	haxisphi->SetTitle("Reconstructed particles #phi, function of their original particle ID");

	hPi0Phi->SetLineColor(1);
	hGamPhi->SetLineColor(4);
	hNeHadPhi->SetLineColor(2);
	hChHadPhi->SetLineColor(7);
	hElePhi->SetLineColor(6);
	//Select the maximum of the histogram to show all lines.
	if(hPi0Phi->GetMaximum() >= hGamPhi->GetMaximum() && hPi0Phi->GetMaximum() >= hNeHadPhi->GetMaximum() && 
	   hPi0Phi->GetMaximum() >= hChHadPhi->GetMaximum() && hPi0Phi->GetMaximum() >= hElePhi->GetMaximum())
		haxisphi->SetMaximum(hPi0Phi->GetMaximum());
	else if(hGamPhi->GetMaximum() >= hPi0Phi->GetMaximum() && hGamPhi->GetMaximum() >= hNeHadPhi->GetMaximum() && 
			hGamPhi->GetMaximum() >= hChHadPhi->GetMaximum() && hGamPhi->GetMaximum() >= hElePhi->GetMaximum())
		haxisphi->SetMaximum(hGamPhi->GetMaximum());
	else if(hNeHadPhi->GetMaximum() >= hPi0Phi->GetMaximum() && hNeHadPhi->GetMaximum() >= hGamPhi->GetMaximum() && 
			hNeHadPhi->GetMaximum() >= hChHadPhi->GetMaximum() && hNeHadPhi->GetMaximum() >= hElePhi->GetMaximum())
		haxisphi->SetMaximum(hNeHadPhi->GetMaximum());	
	else if(hChHadPhi->GetMaximum() >= hPi0Phi->GetMaximum() && hChHadPhi->GetMaximum() >= hNeHadPhi->GetMaximum() && 
			hChHadPhi->GetMaximum() >= hGamPhi->GetMaximum() && hChHadPhi->GetMaximum() >= hElePhi->GetMaximum())
		haxisphi->SetMaximum(hChHadPhi->GetMaximum());
	else if(hElePhi->GetMaximum() >= hPi0Phi->GetMaximum() && hElePhi->GetMaximum() >= hNeHadPhi->GetMaximum() && 
			hElePhi->GetMaximum() >= hChHadPhi->GetMaximum() && hElePhi->GetMaximum() >= hGamPhi->GetMaximum())
		haxisphi->SetMaximum(hElePhi->GetMaximum());
	
	haxisphi->Draw("axis");
	hPi0Phi->Draw("same");
	hGamPhi->Draw("same");
	hNeHadPhi->Draw("same");
	hChHadPhi->Draw("same");
	hElePhi->Draw("same");
	
	sprintf(name,"QA_%s_RecDistributionsGenID.eps",fCalorimeter.Data());
	c11->Print(name);
	
	
	//Ratio reconstructed clusters / generated particles in acceptance, for different particle ID
	//printf("c12\n");
	
	TH1F *	hPi0EClone   = (TH1F*)   hPi0E  ->Clone("hPi0EClone");
	TH1F *	hGamEClone   = (TH1F*)   hGamE  ->Clone("hGamEClone");
	TH1F *	hPi0PtClone  = (TH1F*)   hPi0Pt ->Clone("hPi0PtClone");
	TH1F *	hGamPtClone  = (TH1F*)   hGamPt ->Clone("hGamPtClone");	
	TH1F *	hPi0EtaClone = (TH1F*)   hPi0Eta->Clone("hPi0EtaClone");
	TH1F *	hGamEtaClone = (TH1F*)   hGamEta->Clone("hGamEtaClone");	
	TH1F *	hPi0PhiClone = (TH1F*)   hPi0Phi->Clone("hPi0PhiClone");
	TH1F *	hGamPhiClone = (TH1F*)   hGamPhi->Clone("hGamPhiClone");
	sprintf(cname,"QA_%s_recgenidratio",fCalorimeter.Data());
	TCanvas  * c12 = new TCanvas(cname, "Ratio reconstructed clusters / generated particles in acceptance, for different particle ID", 400, 400) ;
	c12->Divide(2, 2);
	
	c12->cd(1) ; 
	gPad->SetLogy();
	haxisE->SetTitle("Ratio reconstructed clusters / generated particles in acceptance, for different particle ID");
	hPi0EClone->Divide(fhGenPi0AccE);
	hGamEClone->Divide(fhGenGamAccE);
	haxisE->SetMaximum(1.2);
	haxisE->SetYTitle("ratio = rec/gen");
	haxisE->Draw("axis");
	hPi0EClone->SetLineColor(1);
	hGamEClone->SetLineColor(4);
	hPi0EClone->Draw("same");
	hGamEClone->Draw("same");
	
	TLegend pLegend3(0.75,0.45,0.9,0.8);
	pLegend3.SetTextSize(0.06);
	pLegend3.AddEntry(hPi0EClone,"#pi^{0}","L");
	pLegend3.AddEntry(hGamEClone,"#gamma","L");
	pLegend3.SetFillColor(10);
	pLegend3.SetBorderSize(1);
	pLegend3.Draw();
	
	c12->cd(2) ; 
	gPad->SetLogy();
	haxispt->SetTitle("Ratio reconstructed clusters / generated particles in acceptance, for different particle ID");
	hPi0PtClone->Divide(fhGenPi0AccPt);
	hGamPtClone->Divide(fhGenGamAccPt);
	haxispt->SetMaximum(1.2);
	haxispt->SetYTitle("ratio = rec/gen");
	haxispt->Draw("axis");
	hPi0PtClone->SetLineColor(1);
	hGamPtClone->SetLineColor(4);
	hPi0PtClone->Draw("same");
	hGamPtClone->Draw("same");
	
	c12->cd(3) ;
	gPad->SetLogy();
	
	haxiseta->SetTitle("Ratio reconstructed clusters / generated particles in acceptance, for different particle ID");
	hPi0EtaClone->Divide(fhGenPi0AccEta);
	hGamEtaClone->Divide(fhGenGamAccEta);
	haxiseta->SetMaximum(1.2);
	haxiseta->SetYTitle("ratio = rec/gen");
	haxiseta->Draw("axis");
	hPi0EtaClone->SetLineColor(1);
	hGamEtaClone->SetLineColor(4);
	hPi0EtaClone->Draw("same");
	hGamEtaClone->Draw("same");
	
	
	c12->cd(4) ; 
	gPad->SetLogy();
	haxisphi->SetTitle("Ratio reconstructed clusters / generated particles in acceptance, for different particle ID");
	hPi0PhiClone->Divide(fhGenPi0AccPhi);
	hGamPhiClone->Divide(fhGenGamAccPhi);
	haxisphi->SetYTitle("ratio = rec/gen");
	haxisphi->SetMaximum(1.2);
	haxisphi->Draw("axis");
	hPi0PhiClone->SetLineColor(1);
	hGamPhiClone->SetLineColor(4);
	hPi0PhiClone->Draw("same");
	hGamPhiClone->Draw("same");
	
	
	sprintf(name,"QA_%s_EfficiencyGenID.eps",fCalorimeter.Data());
	c12->Print(name);
	
	
	//Reconstructed distributions
	//printf("c1\n");
	sprintf(cname,"QA_%s_vertex",fCalorimeter.Data());
	TCanvas  * c13 = new TCanvas(cname, "Particle vertex", 400, 400) ;
	c13->Divide(2, 2);
	
	c13->cd(1) ; 
	//gPad->SetLogy();
	fhEMVxyz->Draw();
	
	c13->cd(2) ; 
	//gPad->SetLogy();
	fhHaVxyz->Draw();
	
	c13->cd(3) ;
	gPad->SetLogy();
	TH1F * hEMR = (TH1F*) fhEMR->ProjectionY("hEM",-1,-1); 
	hEMR->SetLineColor(4);
	hEMR->Draw();
	
	c13->cd(4) ; 
	gPad->SetLogy();
	TH1F * hHaR = (TH1F*) fhHaR->ProjectionY("hHa",-1,-1); 
	hHaR->SetLineColor(4);
	hHaR->Draw();
	
	
	sprintf(name,"QA_%s_ParticleVertex.eps",fCalorimeter.Data());
	c13->Print(name);
	
	
	//Track-matching distributions
	if(!strcmp(GetReader()->GetInputEvent()->GetName(),"AliESDEvent")){
		sprintf(cname,"QA_%s_trkmatch",fCalorimeter.Data());
		TCanvas *cme = new TCanvas(cname,"Track-matching distributions", 400, 400);
		cme->Divide(2,2);
		
		cme->cd(1);
		fh1pOverE->Draw();
		
		cme->cd(2);
		fh1dR->Draw();
		
		cme->cd(3);
		fh2MatchdEdx->Draw();
		
		cme->cd(4);
		fh2EledEdx->Draw();
		
		sprintf(name,"QA_%s_TrackMatchingEleDist.eps",fCalorimeter.Data());
		cme->Print(name);       
	}
	char line[1024] ; 
	sprintf(line, ".!tar -zcf QA_%s_%s.tar.gz *%s*.eps", fCalorimeter.Data(), GetName(),fCalorimeter.Data()) ; 
	gROOT->ProcessLine(line);
	sprintf(line, ".!rm -fR *.eps"); 
	gROOT->ProcessLine(line);
	
	printf("AliAnaCalorimeterQA::Terminate() - !! All the eps files are in QA_%s_%s.tar.gz !!!\n",  fCalorimeter.Data(), GetName());
	
}
