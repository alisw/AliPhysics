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
 *  Markus Fasel <M.Fasel@gsi.de>
 */
#include <TAxis.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TH1F.h>
#include <TH1I.h>
#include <TH2F.h>
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

#include "AliHFEpid.h"
#include "AliHFEcuts.h"
#include "AliAnalysisElectronTask.h"

//____________________________________________________________
AliAnalysisElectronTask::AliAnalysisElectronTask():
	AliAnalysisTask("PID efficiency Analysis", "")
	, fESD(0x0)
	, fMC(0x0)
	, fCFM(0x0)
	, fPID(0x0)
  , fCuts(0x0)
	, fNEvents(0x0)
	, fQA(0x0)
{
	DefineInput(0, TChain::Class());
	DefineOutput(0, TH1I::Class());
	DefineOutput(1, AliCFContainer::Class());
	DefineOutput(2, TList::Class());

  // Initialize cuts
  fCuts = new AliHFEcuts;
  fPID = new AliHFEpid;
}

//____________________________________________________________
AliAnalysisElectronTask::~AliAnalysisElectronTask(){
	if(fESD) delete fESD;
	if(fMC) delete fMC;
	if(fPID) delete fPID;
	if(fQA) delete fQA;
  if(fCuts) delete fCuts;
	if(fNEvents) delete fNEvents;
	fQA = 0x0; fMC = 0x0; fESD = 0x0; fCFM = 0x0; fPID = 0x0; fCuts = 0x0;
}

//____________________________________________________________
void AliAnalysisElectronTask::ConnectInputData(Option_t *){
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
void AliAnalysisElectronTask::CreateOutputObjects(){
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

  // Initialize correction Framework and Cuts
  fCFM = new AliCFManager;
  MakeParticleContainer();
  // Temporary fix: Initialize particle cuts with 0x0
  for(Int_t istep = 0; istep < fCFM->GetParticleContainer()->GetNStep(); istep++)
    fCFM->SetParticleCutsList(istep, 0x0);
  if(IsQAOn()){
    fCuts->SetDebugMode();
    fQA->AddAt(fCuts->GetQAhistograms(), 7);
  }
  fCuts->CreateStandardCuts();
  fCuts->Initialize(fCFM);

  // Initialize PID
  if(IsQAOn()){
    fPID->SetQAOn();
    fQA->AddAt(fPID->GetQAhistograms(), 8);
  }
  fPID->SetHasMCData(kTRUE);
  fPID->InitializePID("TRD");
}

//____________________________________________________________
void AliAnalysisElectronTask::Exec(Option_t *){
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

	Double_t pt = 0;
	Double_t container[3];

	// Loop over the Monte Carlo tracks to see whether we have overlooked any track
	AliMCParticle *mctrack = 0x0;
	Int_t nElectrons = 0;
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
	AliESDtrack *track = 0x0;
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
  	fCFM->GetParticleContainer()->Fill(container, AliHFEcuts::kStepHFEcuts + 1);
	}
	fNEvents->Fill(1);

	// Done!!!
	PostData(0, fNEvents);
	PostData(1, fCFM->GetParticleContainer());
	PostData(2, fQA);
}

//____________________________________________________________
void AliAnalysisElectronTask::Terminate(Option_t *){
	//
	// Terminate not implemented at the moment
	//
}

//____________________________________________________________
void AliAnalysisElectronTask::MakeParticleContainer(){
  //
  // Create the particle container for the correction framework manager and 
  // link it
  //
  const Int_t nvar   = 3 ; //number of variables on the grid:pt,eta, phi
  const Double_t ptmin = 0., ptmax = 10.;
  const Double_t etamin = -0.9, etamax = 0.9;
  const Double_t phimin = 0., phimax = 2. * TMath::Pi();
  

  //arrays for the number of bins in each dimension
  Int_t iBin[nvar];
  iBin[0] = 20; //bins in pt
  iBin[1] =  8; //bins in eta 
  iBin[2] = 18; // bins in phi

  //arrays for lower bounds :
  Double_t *binLim1 = new Double_t[iBin[0] + 1];
  Double_t *binLim2 = new Double_t[iBin[1] + 1];
  Double_t *binLim3 = new Double_t[iBin[2] + 1];

  //values for bin lower bounds
  for(Int_t i=0; i<=iBin[0]; i++) binLim1[i]=(Double_t)ptmin + (ptmax-ptmin)/iBin[0]*(Double_t)i; 
  for(Int_t i=0; i<=iBin[1]; i++) binLim2[i]=(Double_t)etamin  + (etamax-etamin)/iBin[1]*(Double_t)i;
  for(Int_t i=0; i<=iBin[2]; i++) binLim3[i]=(Double_t)phimin  + (phimax-phimin)/iBin[2]*(Double_t)i;

  //one "container" for MC
  AliCFContainer* container = new AliCFContainer("container","container for tracks", AliHFEcuts::kNcutSteps + 1, nvar, iBin);
  //setting the bin limits
  container -> SetBinLimits(0,binLim1);
  container -> SetBinLimits(1,binLim2);
  container -> SetBinLimits(2,binLim3);
  fCFM->SetParticleContainer(container);
}
