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
	fUseMismatchFileFromHomeDir(kTRUE),
	fUseNSigmaOnPIDAxes(kFALSE),	
	fEventCuts(0x0),
	fTrackCuts(0x0),
	fInclusiveTimes(0x0),
	fT0Fill(0x0),
	fLvsEta(0x0),
	fGlobalPtvsTPCPt(0x0),
	fLvsEtaProjections(0x0),
	fCurrentAODEvent(0x0),
	fCurrentAODTrack(0x0),
	fCurrentDiHadronPIDTrack(0x0),
	fGlobalTracksArray(0x0)

{

	if (fDebug > 0) {cout << Form("File: %s, Line: %i, Function: %s",__FILE__,__LINE__,__func__) << endl;}

}

// -----------------------------------------------------------------------
AliAnalysisTaskCompareAODTrackCuts::AliAnalysisTaskCompareAODTrackCuts(const char* name):
	AliAnalysisTaskSE(name),
	fPIDResponse(0x0),
	fOutputList(0x0),
	fIsMC(kFALSE),
	fVerbose(kFALSE),
	fCalculateTOFMismatch(kFALSE),	
	fUseMismatchFileFromHomeDir(kTRUE),	
	fUseNSigmaOnPIDAxes(kFALSE),	
	fEventCuts(0x0),
	fTrackCuts(0x0),
	fInclusiveTimes(0x0),
	fT0Fill(0x0),
	fLvsEta(0x0),
	fGlobalPtvsTPCPt(0x0),
	fLvsEtaProjections(0x0),
	fCurrentAODEvent(0x0),
	fCurrentAODTrack(0x0),
	fCurrentDiHadronPIDTrack(0x0),
	fGlobalTracksArray(0x0)

{

	if (fDebug > 0) {cout << Form("File: %s, Line: %i, Function: %s",__FILE__,__LINE__,__func__) << endl;}

	// Named Constructor. 
	fTrackCuts = new TObjArray();
	fTrackCuts->SetName("fTrackCuts");

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
	if (!manager) {AliFatal("Could not obtain analysis manager.");}
	AliInputEventHandler* inputHandler = dynamic_cast<AliInputEventHandler*> (manager->GetInputEventHandler());
	if (!inputHandler) {AliFatal("Could not obtain input handler."); return;}	

	// Getting the pointer to the PID response object.
	fPIDResponse = inputHandler->GetPIDResponse();
	if (!fPIDResponse) {AliFatal("Could not obtain PID response."); return;}

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

	// Create the diagram correlating TPC and Global transverse momenta.
	fGlobalPtvsTPCPt = new TH2F("fGlobalPtvsTPCPt","Global p_{T} vs TPC p_{T}; Global p_{T}; TPC p_{T}",100,0,10,100,0,10);
	fOutputList->Add(fGlobalPtvsTPCPt);

	// Creating Global Tracks Array
	fGlobalTracksArray = new TObjArray();

	// Loading the appropriate external mismatch histograms.
	if (fCalculateTOFMismatch) LoadExternalMismatchHistos();

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

		track = dynamic_cast<AliAODTrack*>(fCurrentAODEvent->GetTrack(iTrack));
		if(!track) AliFatal("Not a standard AOD");
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
	Int_t nmismatched = 0;
	Int_t nmatched = 0;
	Int_t nnotof = 0;
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
			//cout << "PR: " << CurrentAODMCParticle->IsPhysicalPrimary() << " SW: " << CurrentAODMCParticle->IsSecondaryFromWeakDecay() << " SM: " << CurrentAODMCParticle->IsSecondaryFromMaterial() << endl;
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
		fCurrentAODTrack = dynamic_cast<AliAODTrack*>(fCurrentAODEvent->GetTrack(iTrack));
		if(!fCurrentAODTrack) AliFatal("Not a standard AOD");
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
		
			// Fill histogram of Global p_T vs TPC-only p_T
			fGlobalPtvsTPCPt->Fill(GetGlobalTrack(fCurrentAODTrack)->Pt(), fCurrentAODTrack->Pt());

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
		if (fCalculateTOFMismatch) rndhittime = GenerateRandomHit(fCurrentDiHadronPIDTrack->Eta());

		// Loop over all Track Cuts.
		for (Int_t iCuts = 0; iCuts < fTrackCuts->GetEntries(); iCuts++) {
			if (fIsMC) ((AliAODTrackCutsDiHadronPID*)fTrackCuts->At(iCuts))->IsSelectedReconstructedMC(fCurrentDiHadronPIDTrack);
			else ((AliAODTrackCutsDiHadronPID*)fTrackCuts->At(iCuts))->IsSelectedData(fCurrentDiHadronPIDTrack,rndhittime);
		}

		Int_t tofmatchstat = fCurrentDiHadronPIDTrack->GetTOFMatchingStatus();
		if (tofmatchstat==0) {nmatched++;}
		if (tofmatchstat==1) {nmismatched++;}
		if (tofmatchstat==2) {nnotof++;}

		// Delete Current DiHadronPIDTrack.
		delete fCurrentDiHadronPIDTrack;
		fCurrentDiHadronPIDTrack = 0x0;

	}

	// Tell the track cuts object that the event has been processed.
	for (Int_t iCuts = 0; iCuts < fTrackCuts->GetEntries(); iCuts++) {
		((AliAODTrackCutsDiHadronPID*)fTrackCuts->At(iCuts))->EventIsDone(fIsMC);
	}
	
	//cout << "Matched: "<<nmatched<<" Mismatched: "<<nmismatched<<" No TOF hit: "<<nnotof<<endl;

	PostData(1,fOutputList);

} 

// -----------------------------------------------------------------------
Bool_t AliAnalysisTaskCompareAODTrackCuts::LoadExternalMismatchHistos() {

	//
	// Attempting to load a root file containing information needed
	// to generate random TOF hits.
 	//

	if (fDebug > 0) {cout << Form("File: %s, Line: %i, Function: %s",__FILE__,__LINE__,__func__) << endl;}

	// Opening external TOF file.
	TFile* fin = 0x0;
	
	// The default is that the file TOFmismatchHistos.root is taken from the /rootfiles/ directory, 
	// and this works fine when running on the train. When the user submits the jobs himself, he can
	// choose to not take the file from the home dir, but upload one along with the source code.
	// If in this case the file is not found, the program tries to get the file from the root directory
	// anyway.
	if (fUseMismatchFileFromHomeDir == kFALSE) {
		fin = TFile::Open("TOFmismatchHistos.root");
		if (!fin) {AliWarning("Tried to open uploaded TOFmismatchHistos.root, but failed");}
	}

	if (!fin) {fin = TFile::Open("alien:///alice/cern.ch/user/m/mveldhoe/rootfiles/TOFmismatchHistos.root");}
	
	if (!fin) {
		AliWarning("Couln't open TOFmismatchHistos.root, will not calculate mismatches...");
		fCalculateTOFMismatch = kFALSE;
		return kFALSE;
	} else {
		AliInfo("Sucessfully loaded TOFmismatchHistos.root");
	}

	// Check if the required histograms are present.
	TH1F* tmp1 = (TH1F*)fin->Get("hNewT0Fill");
	if (!tmp1) {
		AliWarning("Couln't find hNewT0Fill, will not calculate mismatches...");
		fCalculateTOFMismatch = kFALSE;
		return kFALSE;	
	}
	TH2F* tmp2 = (TH2F*)fin->Get("hLvsEta");
	if (!tmp2) {
		AliWarning("Couln't find hLvsEta, will not calculate mismatches...");
		fCalculateTOFMismatch = kFALSE;
		return kFALSE;	
	}	

	// Make a deep copy of the files in the histogram.
	fT0Fill = (TH1F*)tmp1->Clone("fT0Fill");
	fLvsEta = (TH2F*)tmp2->Clone("fLvsEta");

	// Close the external file.
	AliInfo("Closing external file.");
	fin->Close();

	// Creating a TObjArray for LvsEta projections.
	const Int_t nbinseta = fLvsEta->GetNbinsX();
	fLvsEtaProjections = new TObjArray(nbinseta);
	fLvsEtaProjections->SetOwner(kTRUE);

	// Making the projections needed (excluding underflow/ overflow).
	for (Int_t iEtaBin = 1; iEtaBin < (nbinseta + 1); iEtaBin++) {
		TH1F* tmp = (TH1F*)fLvsEta->ProjectionY(Form("LvsEtaProjection_%i",iEtaBin),iEtaBin,iEtaBin);
		tmp->SetDirectory(0);
		fLvsEtaProjections->AddAt(tmp,iEtaBin - 1);
	}

	return kTRUE;

}

// -----------------------------------------------------------------------
Double_t AliAnalysisTaskCompareAODTrackCuts::GenerateRandomHit(Double_t eta) {

	//
	// Returns a random TOF time.
	//

	if (fDebug > 0) {cout << Form("File: %s, Line: %i, Function: %s",__FILE__,__LINE__,__func__) << endl;}

	// Default (error) value:
	Double_t rndhittime = -1.e21;

	// TOF mismatch flag is not turned on.
	if (!fCalculateTOFMismatch) {
		AliFatal("Called GenerateRandomHit() method, but flag fCalculateTOFMismatch not set.");
		return rndhittime;
	}

	// TOF doesn't extend much further than 0.8.
	if (TMath::Abs(eta) > 0.8) {
		if (fDebug) {AliInfo("Tried to get a random hit for a track with eta > 0.8.");}
		return rndhittime;
	}

	// Finding the bin of the eta.
	TAxis* etaAxis = fLvsEta->GetXaxis();
	Int_t etaBin = etaAxis->FindBin(eta);
	if (etaBin == 0 || (etaBin == etaAxis->GetNbins() + 1)) {return rndhittime;}

	const TH1F* lengthDistribution = (const TH1F*)fLvsEtaProjections->At(etaBin - 1);

	if (!lengthDistribution) {
		AliFatal("length Distribution not found.");
		return rndhittime;
	}

	Double_t currentRndLength = lengthDistribution->GetRandom(); // in cm.

	// Similar to Roberto's code.
	Double_t currentRndTime = currentRndLength / (TMath::C() * 1.e2 / 1.e12);
	Double_t t0fill = -1.26416e+04;
	rndhittime = fT0Fill->GetRandom() - t0fill + currentRndTime;

	return rndhittime;

}

// -----------------------------------------------------------------------
void AliAnalysisTaskCompareAODTrackCuts::SetUseNSigmaOnPIDAxes(Bool_t UseNSigma) {

	// Will use NSigma on all PID axes. Will also change all track cuts objects
	// owned by this task.
	if (fDebug > 0) {cout << Form("File: %s, Line: %i, Function: %s",__FILE__,__LINE__,__func__) << endl;}

	fUseNSigmaOnPIDAxes = UseNSigma;

	if (fTrackCuts) {
		for (Int_t iCut = 0; iCut < fTrackCuts->GetSize(); ++iCut) {
			AliAODTrackCutsDiHadronPID* cutstmp = (AliAODTrackCutsDiHadronPID*)(fTrackCuts->At(iCut));
			if (cutstmp) {cutstmp->SetUseNSigmaOnPIDAxes(UseNSigma);}
			else {cout << Form("%s -> WARNING: Found an empty spot in the track cuts array...",__func__) << endl;}
		}
	}
}

// -----------------------------------------------------------------------
void AliAnalysisTaskCompareAODTrackCuts::SetEventCuts(AliAODEventCutsDiHadronPID* eventcuts) {

	if (fDebug > 0) {cout << Form("File: %s, Line: %i, Function: %s",__FILE__,__LINE__,__func__) << endl;}
	if (!eventcuts) {cout << Form("%s -> ERROR: No Event Cuts Object provided.",__func__) << endl; return;}

	fEventCuts = eventcuts;

}

// -----------------------------------------------------------------------
void AliAnalysisTaskCompareAODTrackCuts::AddTrackCuts(AliAODTrackCutsDiHadronPID* trackcuts) {

	if (fDebug > 0) {cout << Form("File: %s, Line: %i, Function: %s",__FILE__,__LINE__,__func__) << endl;}
	if (!trackcuts) {cout << Form("%s -> ERROR: No Track Cuts Object provided.",__func__) << endl; return;}
	if (!fTrackCuts) {cout << Form("%s -> ERROR: No Track Cuts array available.",__func__) << endl; return;}

	// The setting of the task propagates to the imported track cuts object.
	trackcuts->SetUseNSigmaOnPIDAxes(fUseNSigmaOnPIDAxes);

	fTrackCuts->AddLast(trackcuts);

}

// -----------------------------------------------------------------------
void AliAnalysisTaskCompareAODTrackCuts::SetDebugLevel(Int_t debuglvl) {

	// Sets debug level to a certain value, as well as the debug level of the
	// track cuts objects and event cut object.
	if (fDebug > 0) {cout << Form("File: %s, Line: %i, Function: %s",__FILE__,__LINE__,__func__) << endl;}

	fDebug = debuglvl;

	if (fEventCuts) {fEventCuts->SetDebugLevel(debuglvl);}

	if (fTrackCuts) {
		for (Int_t iTrackCutObj = 0; iTrackCutObj < fTrackCuts->GetEntriesFast(); ++iTrackCutObj) {
			((AliTrackDiHadronPID*)fTrackCuts->At(iTrackCutObj))->SetDebugLevel(debuglvl);
		}
	}

}

// -----------------------------------------------------------------------
void AliAnalysisTaskCompareAODTrackCuts::Terminate(Option_t*) {

	if (fDebug > 0) {cout << Form("File: %s, Line: %i, Function: %s",__FILE__,__LINE__,__func__) << endl;}

	if (fCalculateTOFMismatch) {
		delete fT0Fill;
		fT0Fill = 0x0;
		delete fLvsEta;
		fLvsEta = 0x0;
		delete fLvsEtaProjections;
		fLvsEtaProjections = 0x0;
	}

}
