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
 
//-------------------------------------------------------------------------
//  Class AliRsnSelectorRL
// ------------------------
// Reader for conversion of ESD output into the internal format
// used for resonance study.
// ---
// author: A. Pulvirenti             (email: alberto.pulvirenti@ct.infn.it)
// ---
// adapted for TSelector compliance
// by    : R. Vernet                 (email: renaud.vernet@cern.ch)
//-------------------------------------------------------------------------

#include <Riostream.h>

#include <TFile.h>
#include <TChain.h>
#include <TParticle.h>
#include <TRandom.h>
#include <TObjString.h>
#include <TObjectTable.h>
#include <TOrdCollection.h>
#include "AliRun.h"
#include "AliESD.h"
#include "AliStack.h"
#include "AliHeader.h"
#include "AliESDtrack.h"
#include "AliRunLoader.h"
#include "AliRsnDaughter.h"
#include "AliRsnEvent.h"
#include "AliRsnSelectorRL.h"
#include "TAlienFile.h"

ClassImp(AliRsnSelectorRL)

//--------------------------------------------------------------------------------------------------------
AliRsnSelectorRL::AliRsnSelectorRL(TTree*) :
  AliSelectorRL(),
  fOutputPath(0),
  fDebugFlag(0),
  fStoreKineInfo(0),
  fCheckITSRefit(0),
  fRejectFakes(0),
  fCopyMomentum(0),
  fIsRunLoaderOpen(0),
  fRsnEventTree(0),
  fRsnEvent(0),
  fRsnEventBranch(0)
{
//
// Constructor.
// Initializes all working parameters to default values:
// - ESD particle identification
// - rejection of non-ITS-refitted tracks
// - maximum distance allowed from primary vertex = 3.0 cm (beam pipe)
//
	fPIDMethod = kESDPID;
	for (Int_t i = 0; i < AliPID::kSPECIES; i++) fPrior[i] = 1.0;
	fPtLimit4PID = 4.0;
	fProbThreshold = 0.0;
	fMaxRadius = 3.0;      // the beam pipe
}

AliRsnSelectorRL::~AliRsnSelectorRL() {
  Clear();
}
//--------------------------------------------------------------------------------------------------------
AliRsnSelectorRL::AliRsnSelectorRL(const AliRsnSelectorRL& obj) :
  AliSelectorRL(), // not implemented a copy constructor for AliRsnSelectorRL
  fOutputPath(obj.fOutputPath),
  fDebugFlag(obj.fDebugFlag),
  fStoreKineInfo(obj.fStoreKineInfo),
  fCheckITSRefit(obj.fCheckITSRefit),
  fRejectFakes(obj.fRejectFakes),
  fCopyMomentum(obj.fCopyMomentum),
  fIsRunLoaderOpen(obj.fIsRunLoaderOpen),
  fRsnEventTree(obj.fRsnEventTree),
  fRsnEvent(obj.fRsnEvent),
  fRsnEventBranch(obj.fRsnEventBranch)
{
//
// Copy constructor
//
	fPIDMethod = obj.fPIDMethod;
	for (Int_t i = 0; i < AliPID::kSPECIES; i++) fPrior[i] = obj.fPrior[i];
	fPtLimit4PID = obj.fPtLimit4PID;
	fProbThreshold = obj.fProbThreshold;
	fMaxRadius = obj.fMaxRadius;
}
//--------------------------------------------------------------------------------------------------------
AliRsnSelectorRL& AliRsnSelectorRL::operator=(const AliRsnSelectorRL& obj) 
{
//
// Assignment operator
// works in the same way as the copy constructor
//
	if (this!=&obj) {
		fDebugFlag = obj.fDebugFlag;
		fOutputPath = obj.fOutputPath;
		fStoreKineInfo = obj.fStoreKineInfo;
		fCheckITSRefit = obj.fCheckITSRefit;
		fRejectFakes = obj.fRejectFakes;
		fCopyMomentum = obj.fCopyMomentum;
		fIsRunLoaderOpen = obj.fIsRunLoaderOpen;
		fRsnEventTree = obj.fRsnEventTree;
		fRsnEvent = obj.fRsnEvent;
		fRsnEventBranch = obj.fRsnEventBranch;
		fPIDMethod=obj.fPIDMethod;
		for (Int_t i = 0; i < AliPID::kSPECIES; i++) fPrior[i] = obj.fPrior[i];
		fPtLimit4PID = obj.fPtLimit4PID;
		fProbThreshold = obj.fProbThreshold;
		fMaxRadius = obj.fMaxRadius;
	}
	
	return *this;
}
//--------------------------------------------------------------------------------------------------------
void AliRsnSelectorRL::Clear(Option_t * /*option*/)
{
  //
  // Does nothing.
  //
}
//--------------------------------------------------------------------------------------------------------
Double_t* AliRsnSelectorRL::GetPIDprobabilities(AliRsnDaughter track)
{
//
// Computes PID probabilites starting from priors and weights
//
	Int_t i;
	Double_t *prob = new Double_t[AliPID::kSPECIES];
	
	// step 1 - compute the normalization factor
	Double_t sum = 0.0;
	for (i = 0; i < AliPID::kSPECIES; i++) {
		prob[i] = fPrior[i] * track.GetPIDweight(i);
		sum += prob[i];
	}
	if (sum <= 0.0) return 0;
	
	// step 2 - normalize PID weights by the relative prior probability
	for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
		prob[i] /= sum;
	}
	
	return prob;
}
//--------------------------------------------------------------------------------------------------------
void AliRsnSelectorRL::Identify(AliRsnDaughter &track)
{
//
// Identifies a track according to method selected
//
	Bool_t doESDPID = (fPIDMethod == kESDPID);
	Bool_t keepRecSign = (fPIDMethod == kPerfectPIDwithRecSign);
	
	if (doESDPID) {
		// when doing ESD PID it is supposed that charge sign
		// comes from ESD track and it is not modified
		Double_t pt = track.GetPt();
		if (pt <= fPtLimit4PID) {
	      Double_t *prob = GetPIDprobabilities(track);
    	  if (!prob) track.SetPDG(0);
	      Int_t idx[AliPID::kSPECIES];
    	  TMath::Sort(AliPID::kSPECIES, prob, idx);
	      Int_t imax = idx[0];
    	  Double_t maxprob = prob[imax];
	      if (maxprob >= fProbThreshold) {
			track.SetPDG((UShort_t)AliPID::ParticleCode(imax));
	      }
     	 delete [] prob;
    	}
	    else {
    	  track.SetPDG(0);
    	}
	}
	else {
	    Short_t truePDG = track.GetTruePDG();
    	track.SetPDG((UShort_t)TMath::Abs(truePDG));
		if (!keepRecSign) {
			if (TMath::Abs(truePDG) <= 13) {
				if (truePDG > 0) track.SetSign(-1); else track.SetSign(1);
			}
			else {
				if (truePDG > 0) track.SetSign(1); else track.SetSign(-1);
			}
		}
	}
}
//--------------------------------------------------------------------------------------------------------
void AliRsnSelectorRL::SetPriorProbabilities(Double_t *prior)
{
//
// Set prior probabilities to be used in case of ESD PID.
//
	if (!prior) return;
	Int_t i = 0;
	for (i = 0; i < AliPID::kSPECIES; i++) fPrior[i] = prior[i];
}
//--------------------------------------------------------------------------------------------------------
void AliRsnSelectorRL::SetPriorProbability(AliPID::EParticleType type, Double_t value)
{
//
// Sets prior probability referred to a single particle species.
//
	if (type >= AliPID::kElectron && type < AliPID::kPhoton) fPrior[type] = value;
}
//--------------------------------------------------------------------------------------------------------
AliPID::EParticleType AliRsnSelectorRL::FindType(Int_t pdg)
{
//
// Finds the enum value corresponding to a PDG code
//
	pdg = TMath::Abs(pdg);
	switch (pdg) {
		case   11: return AliPID::kElectron; break;
		case   13: return AliPID::kMuon; break;
		case  211: return AliPID::kPion; break;
		case  321: return AliPID::kKaon; break;
		case 2212: return AliPID::kProton; break;
		default  : return AliPID::kPhoton;
	}
}
//--------------------------------------------------------------------------------------------------------

//--------------------------------------------------------
//             The following is Selector stuff
//--------------------------------------------------------

//--------------------------------------------------------------------------------------------------------
void AliRsnSelectorRL::Begin(TTree *)
{
//
// Implementation of BEGIN method
//
	Info("Begin", "");
	TString option = GetOption();
}
//--------------------------------------------------------------------------------------------------------
void AliRsnSelectorRL::SlaveBegin(TTree * tree)
{
//
// Implementation of secondary BEGIN which 
// is called by separate process managers
//
	Info("SlaveBegin", "");
	Init(tree);
	TString option = GetOption();
}
//--------------------------------------------------------------------------------------------------------
void AliRsnSelectorRL::Init(TTree *tree)
{
//
// Initializer
// Connects the selector to a TTree and links the branch
// which is used to make translation ESD --> Rsn
//
	Info("Init","");
	
	if (!tree) return;
	fTree = tree;
	if (fDebugFlag) {
		Info("Init", "fTree=%p   fTree->GetCurrentFile()=%p", fTree, fTree->GetCurrentFile());
	}
	fTree->SetBranchAddress("ESD", &fESD);
	fRsnEventTree = new TTree("selection", "AliRsnEvents");
	TTree::SetBranchStyle(1);
	fRsnEvent = new AliRsnEvent;
	fRsnEventBranch = fRsnEventTree->Branch("events", "AliRsnEvent", &fRsnEvent, 32000, 1);
	fRsnEventBranch->SetAutoDelete(kFALSE);
}
//--------------------------------------------------------------------------------------------------------
Bool_t AliRsnSelectorRL::Process(Long64_t entry)
{
//
// Main core of the Selector processing
// Reads the ESD input and creates a TTree or AliRsnEvent objects
//
	if (fDebugFlag) Info("Process", "Processing event %d", entry);
	if (!AliSelectorRL::Process(entry)) return kFALSE;
	
	AliStack* stack = GetStack();
	if (!stack) {
		Warning("Process", "NULL stack: cannot get kinematics info");
		fStoreKineInfo = kFALSE;
	}
		
	Int_t ntracks, nSelTracks = 0;
		
	// primary vertex
	Double_t vertex[3];
	vertex[0] = fESD->GetVertex()->GetXv();
	vertex[1] = fESD->GetVertex()->GetYv();
	vertex[2] = fESD->GetVertex()->GetZv();
	
	// multiplicity
	ntracks = fESD->GetNumberOfTracks();
	if (!ntracks) {
		Warning("Process", "Event %d has no tracks!", entry);
		return kFALSE;
	}
	if (fDebugFlag) Info("Process", "Number of ESD tracks : %d", ntracks);
	
	// create new AliRsnEvent object
	fRsnEvent->Clear("DELETE");
	fRsnEvent->Init();
	
	// store tracks from ESD
	Int_t index, label;
	Double_t vtot, v[3];
	AliRsnDaughter track;
	for (index = 0; index < ntracks; index++) {
		if (fDebugFlag) Info("Process","Track # %d", index);
		AliESDtrack *esdTrack = fESD->GetTrack(index);
		
		// check against vertex constraint
		esdTrack->GetXYZ(v);
		vtot  = (v[0] - vertex[0])*(v[0] - vertex[0]);
		vtot += (v[1] - vertex[1])*(v[1] - vertex[1]);
		vtot += (v[2] - vertex[2])*(v[2] - vertex[2]);
		vtot  = TMath::Sqrt(vtot);
		if (vtot > fMaxRadius) continue;
		
		// check for fakes
		label = esdTrack->GetLabel();
		if (fRejectFakes && (label < 0)) continue;
		
		// copy ESDtrack data into RsnDaughter (and make Bayesian PID)
		if (!track.Adopt(esdTrack, fCheckITSRefit)) continue;
		track.SetIndex(index);
		
		// if requested, get Kine info
		if (fStoreKineInfo) {
			if (fDebugFlag) Info("Process", "Getting part label=%d stack=%p", label, stack);
			TParticle *part = stack->Particle(TMath::Abs(label));
			if (fDebugFlag) Info("Process", "part=%p", part);
			track.SetTruePDG(part->GetPdgCode());
			Int_t mother = part->GetFirstMother();
			track.SetMother(mother);
			if (mother >= 0) {
				TParticle *mum = stack->Particle(mother);
				track.SetMotherPDG(mum->GetPdgCode());
			}
			// if requested, reconstructed momentum is replaced with true momentum
			if (fCopyMomentum) track.SetPxPyPz(part->Px(), part->Py(), part->Pz());
		}
		
		// identification
		Identify(track);
		
		// store in TClonesArray
		if (fDebugFlag) track.Print("SLITVPMNW");
		fRsnEvent->AddTrack(track);
		nSelTracks++;
	}
	
	// compute total multiplicity
	fRsnEvent->ComputeMultiplicity();
	
	// link to events tree and fill
	fRsnEventTree->Fill();
	
	return kTRUE;
}
//--------------------------------------------------------------------------------------------------------
void AliRsnSelectorRL::SlaveTerminate()
{
//
// SlaveTerminate
// Partial termination method
//
	Info("SlaveTerminate", "");
	
	// Add the histograms to the output on each slave server
	// fOutput->Add(fRsnEventTree);
}
//--------------------------------------------------------------------------------------------------------
void AliRsnSelectorRL::Terminate()
{
//
// Global termination method
//
	Info("Terminate","");
	// fRsnEventTree = dynamic_cast<TTree*>(fOutput->FindObject("fRsnEventTree"));
	
	//AliSelector::Terminate();
	cout << fOutputPath << endl;
	Info("Terminate", Form("Saving in: %s", fOutputPath->Data()));
	TFile* file = TFile::Open(fOutputPath->Data(), "RECREATE");
	fRsnEventTree->Write();
	file->Close();
	
	delete fRsnEventTree;
	delete file;
	delete fOutputPath;
}
