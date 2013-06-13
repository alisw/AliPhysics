/************************************************************************* 
* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. * 
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

// -----------------------------------------------------------------------
//  In this analysis task we compare multiple AOD filtermasks.
// -----------------------------------------------------------------------
//  Author: Misha Veldhoen (misha.veldhoen@cern.ch)

#include <iostream>

// Basic Includes
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TFile.h"
#include "TFormula.h"
#include "TChain.h"
#include "TObject.h"
#include "TRandom3.h"

// Analysis Includes
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSE.h"
#include "AliInputEventHandler.h"

// PID includes.
#include "AliPIDResponse.h"

// AOD includes.
#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliAODHandler.h"
#include "AliAODVertex.h"
#include "AliAODInputHandler.h"
#include "AliAODMCParticle.h"
#include "AliAODMCHeader.h"

// Includes from own library
#include "AliTrackDiHadronPID.h"
#include "AliAODTrackCutsDiHadronPID.h"
#include "AliAODEventCutsDiHadronPID.h"

// AnalysisTask Header
#include "AliAnalysisTaskCompareAODTrackCuts.h"

using namespace std;

ClassImp(AliAnalysisTaskCompareAODTrackCuts);

// -----------------------------------------------------------------------
AliAnalysisTaskCompareAODTrackCuts::AliAnalysisTaskCompareAODTrackCuts():
	AliAnalysisTaskSE(),
	fPIDResponse(0x0),
	fOutputList(0x0),
	fIsMC(kFALSE),
	fVerbose(kFALSE),
	fCalculateTOFMismatch(kFALSE),	
	fMismatchMethod(0),	
	fEventCuts(0x0),
	fTrackCuts(0x0),
	fInclusiveTimes(0x0),
	fExternalTOFfile(0x0),
	fInclusiveTimesIn(0x0),
	fT0Fill(0x0),
	fLvsEta(0x0),
	fCurrentAODEvent(0x0),
	fCurrentAODTrack(0x0),
	fCurrentDiHadronPIDTrack(0x0),
	fGlobalTracksArray(0x0)

{

	if (fDebug > 0) {cout << Form("File: %s, Line: %i, Function: %s",__FILE__,__LINE__,__func__) << endl;}

	// Default Constructor.
	for (Int_t ii = 0; ii<200; ii++) {fLvsEtaProj[ii] = 0x0;}
	for (Int_t ii = 0; ii<100; ii++) {fInclusiveTimesInProj[ii] = 0x0;}

}

// -----------------------------------------------------------------------
AliAnalysisTaskCompareAODTrackCuts::AliAnalysisTaskCompareAODTrackCuts(const char* name):
	AliAnalysisTaskSE(name),
	fPIDResponse(0x0),
	fOutputList(0x0),
	fIsMC(kFALSE),
	fVerbose(kFALSE),
	fCalculateTOFMismatch(kFALSE),	
	fMismatchMethod(0),	
	fEventCuts(0x0),
	fTrackCuts(0x0),
	fInclusiveTimes(0x0),
	fExternalTOFfile(0x0),
	fInclusiveTimesIn(0x0),
	fT0Fill(0x0),
	fLvsEta(0x0),
	fCurrentAODEvent(0x0),
	fCurrentAODTrack(0x0),
	fCurrentDiHadronPIDTrack(0x0),
	fGlobalTracksArray(0x0)

{

	if (fDebug > 0) {cout << Form("File: %s, Line: %i, Function: %s",__FILE__,__LINE__,__func__) << endl;}

	// Named Constructor. 
	fTrackCuts = new TObjArray();
	fTrackCuts->SetName("fTrackCuts");

	for (Int_t ii = 0; ii<200; ii++) {fLvsEtaProj[ii] = 0x0;}
	for (Int_t ii = 0; ii<100; ii++) {fInclusiveTimesInProj[ii] = 0x0;}

	DefineInput(0,TChain::Class());
	DefineOutput(1, TList::Class());

}

// -----------------------------------------------------------------------
AliAnalysisTaskCompareAODTrackCuts::~AliAnalysisTaskCompareAODTrackCuts() {

	if (fDebug > 0) {cout << Form("File: %s, Line: %i, Function: %s",__FILE__,__LINE__,__func__) << endl;}

}

// -----------------------------------------------------------------------
void AliAnalysisTaskCompareAODTrackCuts::UserCreateOutputObjects() {

	if (fDebug > 0) {cout << Form("File: %s, Line: %i, Function: %s",__FILE__,__LINE__,__func__) << endl;}

	// Getting a pointer to the analysis manager and input handler.
	AliAnalysisManager* manager = AliAnalysisManager::GetAnalysisManager();
	AliInputEventHandler* inputHandler = dynamic_cast<AliInputEventHandler*> (manager->GetInputEventHandler());
	
	// Getting the pointer to the PID response object.
	fPIDResponse = inputHandler->GetPIDResponse();

	// Create the output list.
	fOutputList = new TList();
	fOutputList->SetOwner(kTRUE);

	// Adding Event Cuts to the output.
	fEventCuts->CreateHistos();				// Generating all requested histograms.	
	fOutputList->Add(fEventCuts);			

	// Adding Track Cuts to the output.
	for (Int_t iCuts = 0; iCuts < fTrackCuts->GetEntries(); iCuts++) {
		AliAODTrackCutsDiHadronPID* currentcuts = (AliAODTrackCutsDiHadronPID*)fTrackCuts->At(iCuts);
		currentcuts->CreateHistos(); 		// Generating all requested histograms.
		fOutputList->Add(currentcuts);		
	}

	// Creating inclusive times histogram.
	fInclusiveTimes = new TH2F("fInclusiveTimes","Inclusive Times;#eta;t (ps)",100,0,1.,1500,10000.,25000.);
	fOutputList->Add(fInclusiveTimes);

	// Creating Global Tracks Array
	fGlobalTracksArray = new TObjArray();

	// Loading the appropriate external mismatch histograms.
	if (fCalculateTOFMismatch) LoadExternalMismatchHistos(fMismatchMethod);

	PostData(1,fOutputList);

}

// -----------------------------------------------------------------------
void AliAnalysisTaskCompareAODTrackCuts::FillGlobalTracksArray() {

	if (fDebug > 0) {cout << Form("File: %s, Line: %i, Function: %s",__FILE__,__LINE__,__func__) << endl;}

	// Initialize the mapping for corresponding PID tracks.
	// See AliAnalysisTaskESDFilter.cxx for explanation about the mapping
	// between TPC-Only and Global tracks
		
	// Create TObjArray for Global Tracks.

	//cout<<"Global Tracks Array: "<<fGlobalTracksArray<<endl;

	// Clear previous tracks...
	fGlobalTracksArray->Clear();

	//if (fGlobalTracksArray) delete fGlobalTracksArray;
	//fGlobalTracksArray = new TObjArray();

	AliAODTrack* track = 0x0;
		
	for (Int_t iTrack = 0; iTrack < fCurrentAODEvent->GetNumberOfTracks(); iTrack++) {	

		track = fCurrentAODEvent->GetTrack(iTrack);
        if (track->GetID()>-1) fGlobalTracksArray->AddAtAndExpand(track,track->GetID());
	}

}

// -----------------------------------------------------------------------
AliAODTrack* AliAnalysisTaskCompareAODTrackCuts::GetGlobalTrack(AliAODTrack* track) {

	// Returns Global Track corresponding to a TPC-Only Track.
	if (fDebug > 0) {cout << Form("File: %s, Line: %i, Function: %s",__FILE__,__LINE__,__func__) << endl;}

	if (track->GetID() >= 0) {
		cout<<"AliAnalysisTaskCompareAODTrackCuts::GetGlobalTrack -> Input Track is not TPC-Only."<<endl;
		return 0x0;
	}

	if (!fGlobalTracksArray) {
		cout<<"AliAnalysisTaskCompareAODTrackCuts::GetGlobalTrack -> Global Tracks Array Does not Exist."<<endl;
		return 0x0;
	}

    AliAODTrack* globaltrack = (AliAODTrack*)(fGlobalTracksArray->At(-track->GetID()-1));
	if (!globaltrack) {
		cout<<"AliAnalysisTaskCompareAODTrackCuts::GetGlobalTrack -> No Global Track Found."<<endl;
		return 0x0;
	}

	return globaltrack;
	
}

// -----------------------------------------------------------------------
void AliAnalysisTaskCompareAODTrackCuts::UserExec(Option_t*) {

	if (fDebug > 0) {cout << Form("File: %s, Line: %i, Function: %s",__FILE__,__LINE__,__func__) << endl;}

	// Input Current Event.
	fCurrentAODEvent = dynamic_cast<AliAODEvent*>(InputEvent());
	if (!fCurrentAODEvent) {
		cout<<"AliAnalysisTaskCompareAODTrackCuts::UserExec -> AOD Event not found."<<endl;
		return;
	}

	// Check Event Cuts Object.
	if (!fEventCuts) AliFatal("No Event Cuts Object Found!");

	// Perform Event Cuts.
	if (!fEventCuts->IsSelected(fCurrentAODEvent)) {return;}

	// Check Track Cuts Array.
	if (!fTrackCuts) AliFatal("No Track Cuts Array Found!");
	if (fTrackCuts->GetEntries() == 0) AliFatal("Track Cuts Array is Empty!");

	// If MC, then fill MC reconstructed QA histograms.
	TClonesArray* mcArray = 0x0;
	TObjArray* mcArrayLabel = 0x0; 
	if (fIsMC) {
		
		mcArray = dynamic_cast<TClonesArray*>(fCurrentAODEvent->FindListObject(AliAODMCParticle::StdBranchName()));
    	if (!mcArray) {
        	AliFatal("No MC array found in the AOD.");
        	return;
    	}

    	mcArrayLabel = new TObjArray(10000);

    	// Loop over MC particles.
    	for (Int_t iParticle = 0; iParticle < mcArray->GetEntriesFast(); iParticle++) {
		
			// Put the MC Particle in the Label array.
			AliAODMCParticle* CurrentAODMCParticle = (AliAODMCParticle*) mcArray->At(iParticle);
			mcArrayLabel->AddAtAndExpand(CurrentAODMCParticle,CurrentAODMCParticle->Label());
    		//cout<<"Index: "<<iParticle<<" Label: "<<CurrentAODMCParticle->Label()<<endl;

			//if (CurrentAODMCParticle->Label()!=iParticle) cout<<"Index unequal to particle's label!"<<endl;

			// Loop over all Track Cuts.
			for (Int_t iCuts = 0; iCuts < fTrackCuts->GetEntries(); iCuts++) {
				((AliAODTrackCutsDiHadronPID*)fTrackCuts->At(iCuts))->IsSelectedGeneratedMC(CurrentAODMCParticle);
			}
		}
	}

	//for (Int_t ii=0; ii<200; ii++) cout<<fLvsEtaProj[ii]<<endl;		

	// Create mapping to Global Tracks.
	FillGlobalTracksArray();

	// Tell the track cuts object that a new event has started.
	for (Int_t iCuts = 0; iCuts < fTrackCuts->GetEntries(); iCuts++) {
		AliAODTrackCutsDiHadronPID* currentcuts = (AliAODTrackCutsDiHadronPID*)fTrackCuts->At(iCuts);
		if (fVerbose) AliInfo(Form("Starting new event for cuts: %s",currentcuts->GetName()));
		currentcuts->StartNewEvent();
	}

	// Loop over all reconstructed tracks.
	for (Int_t iTrack = 0; iTrack < fCurrentAODEvent->GetNumberOfTracks(); iTrack++) {

		// Get the Current Track.
		fCurrentAODTrack = (AliAODTrack*)fCurrentAODEvent->GetTrack(iTrack);
		if (!fCurrentAODTrack) {
			cout<<"AliAnalysisTaskCompareAODTrackCuts::UserExec -> AOD Track not found."<<endl;
			continue;
		}

		// Ignore muon tracks.
		if (fCurrentAODTrack->IsMuonTrack()) {continue;}

		// Copy Track info into AliTrackDiHadronPID object.
		fCurrentDiHadronPIDTrack = 0x0;

		// Check if we can indeed find a MC particle corresponding to the reconstructed track.
		/*
		AliAODMCParticle* MCparticleCheck = (AliAODMCParticle*)mcArrayLabel->At(TMath::Abs(fCurrentAODTrack->GetLabel()));
		Double_t MCPt1 = -999.;
		if (MCparticleCheck) MCPt1 = MCparticleCheck->Pt();
		AliAODMCParticle* MCparticleCheck2 = (AliAODMCParticle*)mcArray->At(TMath::Abs(fCurrentAODTrack->GetLabel()));	
		Double_t MCPt2 = -999.;
		if (MCparticleCheck2) MCPt2 = MCparticleCheck2->Pt();

		cout<<"AOD track abs label: "<<TMath::Abs(fCurrentAODTrack->GetLabel())<<" MC particle with same label found: "<<(Bool_t)MCparticleCheck<<" and with same index: "<<(Bool_t)MCparticleCheck2<<endl;
		cout<<"AOD track pt: "<<fCurrentAODTrack->Pt()<<" same label pt: "<<MCPt1<<" same index pt: "<<MCPt2<<endl;
		//if (!MCparticleCheck) cout<<"MC Particle for a reconstructed track not found."<<endl;
*/
		// Check if it's a TPC-Only or Global Track.
		if (fCurrentAODTrack->GetID() < 0) {
			// Q: Do we really need to create the arraylabel object or is the original MC array already nicely ordered. (NOTE THAT WE TAKE THE MC PARTICLES FROM THE mcArray)
			if (fIsMC) fCurrentDiHadronPIDTrack = new AliTrackDiHadronPID(fCurrentAODTrack,GetGlobalTrack(fCurrentAODTrack),(AliAODMCParticle*)mcArray->At(TMath::Abs(fCurrentAODTrack->GetLabel())),fPIDResponse);
			else fCurrentDiHadronPIDTrack = new AliTrackDiHadronPID(fCurrentAODTrack,GetGlobalTrack(fCurrentAODTrack),0x0,fPIDResponse);
		} else {
			if (fIsMC) fCurrentDiHadronPIDTrack = new AliTrackDiHadronPID(fCurrentAODTrack,0x0,(AliAODMCParticle*)mcArray->At(TMath::Abs(fCurrentAODTrack->GetLabel())),fPIDResponse);
			else fCurrentDiHadronPIDTrack = new AliTrackDiHadronPID(fCurrentAODTrack,0x0,0x0,fPIDResponse);
		}

		if (!fCurrentDiHadronPIDTrack) {
			cout<<"AliAnalysisTaskCompareAODTrackCuts::UserExec -> Copying to DiHadronPIDTrack failed."<<endl;
			continue;
		}

		// Filling random times histogram:
		ULong_t requestedflags = (UInt_t)(AliAODTrack::kTOFout)|(UInt_t)(AliAODTrack::kTIME);
		ULong_t trackflags = fCurrentDiHadronPIDTrack->GetFlags();
		if (requestedflags&trackflags) {fInclusiveTimes->Fill(TMath::Abs(fCurrentDiHadronPIDTrack->Eta()),fCurrentDiHadronPIDTrack->GetTOFsignal());}

		Double_t rndhittime = -1.e21;
		if (fCalculateTOFMismatch) rndhittime = GenerateRandomHit(fMismatchMethod,fCurrentDiHadronPIDTrack->Eta());

		// Loop over all Track Cuts.
		for (Int_t iCuts = 0; iCuts < fTrackCuts->GetEntries(); iCuts++) {
			if (fIsMC) ((AliAODTrackCutsDiHadronPID*)fTrackCuts->At(iCuts))->IsSelectedReconstructedMC(fCurrentDiHadronPIDTrack);
			else ((AliAODTrackCutsDiHadronPID*)fTrackCuts->At(iCuts))->IsSelectedData(fCurrentDiHadronPIDTrack,rndhittime);
		}

		// Delete Current DiHadronPIDTrack.
		delete fCurrentDiHadronPIDTrack;
		fCurrentDiHadronPIDTrack = 0x0;

	}

	// Tell the track cuts object that the event has been processed.
	for (Int_t iCuts = 0; iCuts < fTrackCuts->GetEntries(); iCuts++) {
		((AliAODTrackCutsDiHadronPID*)fTrackCuts->At(iCuts))->EventIsDone(fIsMC);
	}
	
	PostData(1,fOutputList);

} 

// -----------------------------------------------------------------------
Bool_t AliAnalysisTaskCompareAODTrackCuts::LoadExternalMismatchHistos(Int_t mismatchmethod) {

	if (fDebug > 0) {cout << Form("File: %s, Line: %i, Function: %s",__FILE__,__LINE__,__func__) << endl;}

	// Inputted histograms are loaded and projected properly.
	if ((mismatchmethod < 0)||(mismatchmethod > 1)) {return 0;}

	if (mismatchmethod == 0) {

		AliInfo(Form("Loading histograms in file TOFmismatchHistos.root..."));

		//fExternalTOFfile = TFile::Open("TOFmismatchHistos.root");
		fExternalTOFfile = new TFile("TOFmismatchHistos.root","read");
		fLvsEta = (TH2F*)fExternalTOFfile->Get("hLvsEta");
		fT0Fill = (TH1F*)fExternalTOFfile->Get("hNewT0Fill");

		// Make projections on the L axis.
		for (Int_t ii = 0; ii < fLvsEta->GetNbinsX(); ii++) {

			fLvsEtaProj[ii] = (TH1F*)fLvsEta->ProjectionY(Form("fLvsEtaProj_%i",ii),ii+1,ii+1);
			fLvsEtaProj[ii]->SetDirectory(0);

		}

	}

	if (mismatchmethod == 1) {

		AliInfo(Form("Loading histograms in file InclusiveTimes.root..."));

		//fExternalTOFfile = TFile::Open("InclusiveTimes.root");
		fExternalTOFfile = new TFile("InclusiveTimes.root","read");
		fInclusiveTimesIn = (TH2F*)fExternalTOFfile->Get("fInclusiveTimesIn");

		// Make projections on the time axis.
		for (Int_t ii = 0; ii < fInclusiveTimesIn->GetNbinsX(); ii++) {

			fInclusiveTimesInProj[ii] = (TH1F*)fInclusiveTimesIn->ProjectionY(Form("fInclusiveTimesInProj_%i",ii),ii+1,ii+1);
			fInclusiveTimesInProj[ii]->SetDirectory(0);

		}
	}

	return 0;

}

// -----------------------------------------------------------------------
Bool_t AliAnalysisTaskCompareAODTrackCuts::UnLoadExternalMismatchHistos() {

	if (fDebug > 0) {cout << Form("File: %s, Line: %i, Function: %s",__FILE__,__LINE__,__func__) << endl;}

	AliInfo("Unloading external histo's");

	for (Int_t ii = 0; ii < 100; ii++){
		if (fLvsEtaProj[ii]) {
			delete fLvsEtaProj[ii]; 
			fLvsEtaProj[ii] = 0x0;
		}
	} 
	for (Int_t ii = 0; ii < 200; ii++) {
		if (fInclusiveTimesInProj[ii]) {
			delete fInclusiveTimesInProj[ii]; 
			fInclusiveTimesInProj[ii] = 0x0;
			}
		}

	if (fExternalTOFfile) {fExternalTOFfile->Close();}
	
	return 0;

}

// -----------------------------------------------------------------------
Double_t AliAnalysisTaskCompareAODTrackCuts::GenerateRandomHit(Int_t mismatchmethod, Double_t eta) {

	if (fDebug > 0) {cout << Form("File: %s, Line: %i, Function: %s",__FILE__,__LINE__,__func__) << endl;}

	if ((mismatchmethod < 0)||(mismatchmethod > 1)) {return -1.e21;}
	if (TMath::Abs(eta) > 0.8) {return -1.e21;}

	Double_t rndhittime = -1.e21;

	// Method as in Roberto's code.
	if (mismatchmethod == 0) {

		Int_t etabin = (Int_t)((eta+1.)*100.);

		Double_t currentRndLength = fLvsEtaProj[etabin]->GetRandom(); // in cm.

		Double_t currentRndTime = currentRndLength / (TMath::C() * 1.e2 / 1.e12);
		Double_t t0fill = -1.26416e+04;

		rndhittime = fT0Fill->GetRandom() - t0fill + currentRndTime;

	}

	// Using length inclusive time distributions.
	if (fMismatchMethod == 1) {
		
		Double_t abseta = TMath::Abs(eta);
		Int_t etabin = (Int_t)(abseta*100.);
		rndhittime = fInclusiveTimesInProj[etabin]->GetRandom();
		
	}

	return rndhittime;

}

// -----------------------------------------------------------------------
void AliAnalysisTaskCompareAODTrackCuts::Terminate(Option_t*) {

	if (fDebug > 0) {cout << Form("File: %s, Line: %i, Function: %s",__FILE__,__LINE__,__func__) << endl;}

	if (fCalculateTOFMismatch) UnLoadExternalMismatchHistos();

}

// -----------------------------------------------------------------------