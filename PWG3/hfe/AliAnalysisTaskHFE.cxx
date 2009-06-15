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
/*
 * The analysis task:
 * Filling an AliCFContainer with the quantities pt, eta and phi
 * for tracks which survivied the particle cuts (MC resp. ESD tracks)
 * Track selection is done using the AliHFE package
 * 
 * Author:
 *  Raphaelle Bailhache <R.Bailhache@gsi.de>
 *  Markus Fasel <M.Fasel@gsi.de>
 *  MinJung Kweon <minjung@physi.uni-heidelberg.de>
 */
#include <TAxis.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TDirectory.h>
#include <TH1F.h>
#include <TH1I.h>
#include <TH2F.h>
#include <THnSparse.h>
#include <TIterator.h>
#include <TList.h>
#include <TLegend.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TParticle.h>
#include <TProfile.h>
#include <TTree.h>

#include "AliCFContainer.h"
#include "AliCFManager.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliLog.h"
#include "AliAnalysisManager.h"
#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
#include "AliMCParticle.h"
#include "AliPID.h"
#include "AliStack.h"

#include "AliHFEpid.h"
#include "AliHFEcuts.h"
#include "AliHFEmcQA.h"
#include "AliHFEsecVtx.h"
#include "AliAnalysisTaskHFE.h"

//____________________________________________________________
AliAnalysisTaskHFE::AliAnalysisTaskHFE():
  AliAnalysisTask("PID efficiency Analysis", "")
  , fESD(0x0)
  , fMC(0x0)
  , fCFM(0x0)
  , fCorrelation(0x0)
  , fPID(0x0)
  , fCuts(0x0)
  , fSecVtx(0x0)
  , fAnalysisMCQA(0x0)
  , fNEvents(0x0)
  , fQA(0x0)
  , fOutput(0x0)
  , fHistMCQA(0x0)
  , fHistSECVTX(0x0)
{
	DefineInput(0, TChain::Class());
	DefineOutput(0, TH1I::Class());
	DefineOutput(1, TList::Class());
  DefineOutput(2, TList::Class());

  // Initialize cuts
  fCuts = new AliHFEcuts;
  fPID = new AliHFEpid;
}

//____________________________________________________________
AliAnalysisTaskHFE::~AliAnalysisTaskHFE(){
	if(fESD) delete fESD;
	if(fMC) delete fMC;
	if(fPID) delete fPID;
	if(fQA){
    fQA->Clear();
    delete fQA;
  }
  if(fOutput){ 
    fOutput->Clear();
    delete fOutput;
  }
  if(fHistMCQA){
    fHistMCQA->Clear();
    delete fHistMCQA;
  }
  if(fHistSECVTX){
    fHistSECVTX->Clear();
    delete fHistSECVTX;
  }
  if(fCuts) delete fCuts;
  if(fSecVtx) delete fSecVtx;
  if(fAnalysisMCQA) delete fAnalysisMCQA;
	if(fNEvents) delete fNEvents;
  if(fCorrelation) delete fCorrelation;
  if(fFakeElectrons) delete fFakeElectrons;
}

//____________________________________________________________
void AliAnalysisTaskHFE::ConnectInputData(Option_t *){
	TTree *esdchain = dynamic_cast<TChain *>(GetInputData(0));
	if(!esdchain){
		AliError("ESD chain empty");
		return;
	} else {
		esdchain->SetBranchStatus("Tracks", 1);
	}
	AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler *>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
	if(!esdH){	
		AliError("No ESD input handler");
		return;
	} else {
		fESD = esdH->GetEvent();
	}
	AliMCEventHandler *mcH = dynamic_cast<AliMCEventHandler *>(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
	if(!mcH){	
		AliError("No MC truth handler");
		return;
	} else {
		fMC = mcH->MCEvent();
	}
}

//____________________________________________________________
void AliAnalysisTaskHFE::CreateOutputObjects(){
	fNEvents = new TH1I("nEvents", "Number of Events in the Analysis", 2, 0, 2); // Number of Events neccessary for the analysis and not a QA histogram
	// First Step: TRD alone
	if(!fQA) fQA = new TList;
	fQA->AddAt(new TProfile("conr", "Electron PID contamination", 20, 0, 20), 0);
	fQA->AddAt(new TH1F("alpha_rec", "Alpha from reconstructed tracks with TRD hits", 36, -TMath::Pi(), TMath::Pi()), 1);
	fQA->AddAt(new TH1F("alpha_sim", "Alpha from simulated electron tracks", 36, -TMath::Pi(), TMath::Pi()), 2);
	fQA->AddAt(new TH1F("nElectron", "Number of electrons", 100, 0, 100), 3);
	fQA->AddAt(new TProfile("pidquality", "TRD PID quality as function of momentum", 20, 0, 20), 4);
	fQA->AddAt(new TProfile("ntrdclusters", "Number of TRD clusters as function of momentum", 20, 0, 20), 5);
	fQA->AddAt(new TH1F("chi2TRD","#chi2 per TRD cluster", 20, 0, 20), 6);

  if(!fOutput) fOutput = new TList;
  // Initialize correction Framework and Cuts
  fCFM = new AliCFManager;
  MakeParticleContainer();
  // Temporary fix: Initialize particle cuts with 0x0
  for(Int_t istep = 0; istep < fCFM->GetParticleContainer()->GetNStep(); istep++)
    fCFM->SetParticleCutsList(istep, 0x0);
  if(IsQAOn()){
    fCuts->SetDebugMode();
    fQA->Add(fCuts->GetQAhistograms());
  }
  fCuts->CreateStandardCuts();
  fCuts->Initialize(fCFM);
  // add output objects to the List
  fOutput->AddAt(fCFM->GetParticleContainer(), 0);
  fOutput->AddAt(fCorrelation, 1);
  fOutput->AddAt(fFakeElectrons, 2);

  // Initialize PID
  if(IsQAOn()){
    fPID->SetQAOn();
    fQA->Add(fPID->GetQAhistograms());
  }
  fPID->SetHasMCData(kTRUE);
  fPID->InitializePID("TRD");

  // [mj] mcQA----------------------------------
  if (IsMCQAOn()) {
    if(!fHistMCQA) fHistMCQA = new TList();
    fQA->Add(fHistMCQA);
    fAnalysisMCQA->CreateHistograms(AliHFEmcQA::fkCharm,"mcqa_");               // create histograms for charm
    fAnalysisMCQA->CreateHistograms(AliHFEmcQA::fkBeauty,"mcqa_");              // create histograms for beauty
    //Rossella
//     fAnalysisMCQA->CreateHistosHadrons();                   // create histograms for hadrons

  }
  // [mj] secvtx----------------------------
  if (IsSecVtxOn()) {
    fOutput->Add(fHistSECVTX);
    if(!fHistSECVTX) fHistSECVTX = new TList();
    fSecVtx->CreateHistograms("secvtx_");
  }

  TIter next_(gDirectory->GetList());
  TObject *obj_;
  int counter_ = 0;
  int counter__ = 0;
  TString objname;
  while ((obj_ = next_.Next())) {
    objname = obj_->GetName();
    TObjArray *toks = objname.Tokenize("_");
    if (toks->GetEntriesFast()){
      TObjString *fpart = (TObjString *)(toks->UncheckedAt(0));
      if ((fpart->String()).CompareTo("secvtx") == 0){
        fHistSECVTX->AddAt(obj_, counter_++);
      }
      else if ((fpart->String()).CompareTo("mcqa") == 0){
        fHistMCQA->AddAt(obj_, counter__++);
      }
    }
  }
  //--------------------------------------- 
}

//____________________________________________________________
void AliAnalysisTaskHFE::Exec(Option_t *){
	//
	// Run the analysis
	//
	if(!fESD){
		AliError("No ESD Event");
		return;
	}
	if(!fMC){
		AliError("No MC Event");
		return;
	}
	fCFM->SetEventInfo(fMC);
  fPID->SetMCEvent(fMC);

	//fCFM->CheckEventCuts(AliCFManager::kEvtGenCuts, fMC);

	Double_t container[6];

	// Loop over the Monte Carlo tracks to see whether we have overlooked any track
	AliMCParticle *mctrack = 0x0;
	Int_t nElectrons = 0;

  // [mj] run MC QA ------------------------------------------------
  if (IsMCQAOn()) {

    fAnalysisMCQA->SetStack(fMC->Stack());
    fAnalysisMCQA->Init();

    Int_t nPrims = fMC->Stack()->GetNprimary();
    Int_t nMCTracks = fMC->Stack()->GetNtrack();

    // loop over primary particles for quark and heavy hadrons
    for (Int_t igen = 0; igen < nPrims; igen++){
      fAnalysisMCQA->GetQuarkKine(igen, AliHFEmcQA::fkCharm);
      fAnalysisMCQA->GetQuarkKine(igen, AliHFEmcQA::fkBeauty);
    }
    fAnalysisMCQA->EndOfEventAna(AliHFEmcQA::fkCharm);
    fAnalysisMCQA->EndOfEventAna(AliHFEmcQA::fkBeauty);

    // loop over all tracks for decayed electrons
    for (Int_t igen = 0; igen < nMCTracks; igen++){
      fAnalysisMCQA->GetDecayedKine(igen, AliHFEmcQA::fkCharm, AliHFEmcQA::fkElectron);
      fAnalysisMCQA->GetDecayedKine(igen, AliHFEmcQA::fkBeauty, AliHFEmcQA::fkElectron);
    }

  } // end of MC QA loop
  // -----------------------------------------------------------------

	for(Int_t imc = fMC->GetNumberOfTracks(); imc--;){
		mctrack = fMC->GetTrack(imc);
		container[0] = mctrack->Pt();
		container[1] = mctrack->Eta();
    container[2] = mctrack->Phi();

		if(!fCFM->CheckParticleCuts(AliHFEcuts::kStepMCGenerated, mctrack)) continue;
		fCFM->GetParticleContainer()->Fill(container, AliHFEcuts::kStepMCGenerated);
		(dynamic_cast<TH1F *>(fQA->At(2)))->Fill(mctrack->Phi() - TMath::Pi());
		if(!fCFM->CheckParticleCuts(AliHFEcuts::kStepMCInAcceptance, mctrack)) continue;
		// find the label in the vector
		fCFM->GetParticleContainer()->Fill(container, AliHFEcuts::kStepMCInAcceptance);
		nElectrons++;
	}
	(dynamic_cast<TH1F *>(fQA->At(3)))->Fill(nElectrons);

	// fCFM->CheckEventCuts(AliCFManager::kEvtRecCuts, fESD);
	AliESDtrack *track = 0x0, *htrack = 0x0;
  Int_t pid = 0;
	for(Int_t itrack = 0; itrack < fESD->GetNumberOfTracks(); itrack++){
		track = fESD->GetTrack(itrack);
		container[0] = track->Pt();
		container[1] = track->Eta();
    container[2] = track->Phi();
		if(!fCFM->CheckParticleCuts(AliHFEcuts::kStepRecKine, track)) continue;
		fCFM->GetParticleContainer()->Fill(container, AliHFEcuts::kStepRecKine);

		// Check TRD criterions (outside the correction framework)
		if(track->GetTRDncls()){
			(dynamic_cast<TH1F *>(fQA->At(6)))->Fill(track->GetTRDchi2()/track->GetTRDncls());
			(dynamic_cast<TH1F *>(fQA->At(1)))->Fill(track->GetAlpha());	// Check the acceptance without tight cuts
			(dynamic_cast<TProfile *>(fQA->At(4)))->Fill(container[0], track->GetTRDpidQuality());
			(dynamic_cast<TProfile *>(fQA->At(5)))->Fill(container[0], track->GetTRDncls());
		}
		if(!fCFM->CheckParticleCuts(AliHFEcuts::kStepRecPrim, track)) continue;
    fCFM->GetParticleContainer()->Fill(container, AliHFEcuts::kStepRecPrim);
    if(fCFM->CheckParticleCuts(AliHFEcuts::kStepHFEcuts, track)) continue;
    fCFM->GetParticleContainer()->Fill(container, AliHFEcuts::kStepHFEcuts);
    // track accepted, do PID
		if(!fPID->IsSelected(track)) continue;
    // Track selected: distinguish between true and fake
    if(!(mctrack = fMC->GetTrack(TMath::Abs(track->GetLabel())))) continue;
    if((pid = TMath::Abs(mctrack->Particle()->GetPdgCode())) == 11){
  	  fCFM->GetParticleContainer()->Fill(container, AliHFEcuts::kStepHFEcuts + 1);
      // dimensions 3&4&5 : pt,eta,phi (MC)
      container[3] = mctrack->Pt();
      container[4] = mctrack->Eta();
      container[5] = mctrack->Phi();
      fCorrelation->Fill(container);
      // pair analysis [mj]
      if (IsSecVtxOn()) {
        for(Int_t jtrack = 0; jtrack < fESD->GetNumberOfTracks(); jtrack++){
          htrack = fESD->GetTrack(jtrack);
          if ( itrack == jtrack ) continue;
          //if( fPID->IsSelected(htrack) && (itrack < jtrack)) continue;
          if( abs(fSecVtx->GetMCPID(track)) == 11 && (itrack < jtrack)) continue;
          fSecVtx->AnaPair(track, htrack, fSecVtx->PairCode(track,htrack), jtrack);
        }
      }
    } else {
      // Fill THnSparse with the information for Fake Electrons
      switch(pid){
        case 13:    container[3] = AliPID::kMuon;
        case 211:   container[3] = AliPID::kPion;
        case 321:   container[3] = AliPID::kKaon;
        case 2212:  container[3] = AliPID::kProton;
      }
      fFakeElectrons->Fill(container);
    }
	}
	fNEvents->Fill(1);

	// Done!!!
	PostData(0, fNEvents);
	PostData(1, fOutput);
	PostData(2, fQA);
}

//____________________________________________________________
void AliAnalysisTaskHFE::Terminate(Option_t *){
	//
	// Terminate not implemented at the moment
	//
}

//____________________________________________________________
void AliAnalysisTaskHFE::MakeParticleContainer(){
  //
  // Create the particle container for the correction framework manager and 
  // link it
  //
  const Int_t kNvar   = 3 ; //number of variables on the grid:pt,eta, phi
  const Double_t kPtmin = 0., kPtmax = 10.;
  const Double_t kEtamin = -0.9, kEtamax = 0.9;
  const Double_t kPhimin = 0., kPhimax = 2. * TMath::Pi();

  //arrays for the number of bins in each dimension
  Int_t iBin[kNvar];
  iBin[0] = 20; //bins in pt
  iBin[1] =  8; //bins in eta 
  iBin[2] = 18; // bins in phi

  //arrays for lower bounds :
  Double_t* binEdges[kNvar];
  for(Int_t ivar = 0; ivar < kNvar; ivar++)
    binEdges[ivar] = new Double_t[iBin[ivar] + 1];

  //values for bin lower bounds
  for(Int_t i=0; i<=iBin[0]; i++) binEdges[0][i]=(Double_t)kPtmin + (kPtmax-kPtmin)/iBin[0]*(Double_t)i; 
  for(Int_t i=0; i<=iBin[1]; i++) binEdges[1][i]=(Double_t)kEtamin  + (kEtamax-kEtamin)/iBin[1]*(Double_t)i;
  for(Int_t i=0; i<=iBin[2]; i++) binEdges[2][i]=(Double_t)kPhimin  + (kPhimax-kPhimin)/iBin[2]*(Double_t)i;

  //one "container" for MC
  AliCFContainer* container = new AliCFContainer("container","container for tracks", AliHFEcuts::kNcutSteps + 1, kNvar, iBin);
  //setting the bin limits
  for(Int_t ivar = 0; ivar < kNvar; ivar++)  container -> SetBinLimits(0, binEdges[ivar]);
  fCFM->SetParticleContainer(container);

  //create correlation matrix for unfolding
  Int_t thnDim[2*kNvar];
  for (int k=0; k<kNvar; k++) {
    //first half  : reconstructed 
    //second half : MC
    thnDim[k]      = iBin[k];
    thnDim[k+kNvar] = iBin[k];
  }

  fCorrelation = new THnSparseF("correlation","THnSparse with correlations",2*kNvar,thnDim);
  for (int k=0; k<kNvar; k++) {
    fCorrelation->SetBinEdges(k,binEdges[k]);
    fCorrelation->SetBinEdges(k+kNvar,binEdges[k]);
  }
  fCorrelation->Sumw2();

  // Add a histogram for Fake electrons
  thnDim[3] = AliPID::kSPECIES;
  fFakeElectrons = new THnSparseF("fakeEkectrons", "Output for Fake Electrons", kNvar + 1, thnDim);
  for(Int_t idim = 0; idim < kNvar; idim++)
    fFakeElectrons->SetBinEdges(idim, binEdges[idim]);
}
