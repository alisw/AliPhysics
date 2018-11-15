#ifndef ALITRACKDIHADRONPID_H
#define ALITRACKDIHADRONPID_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. * 
* See cxx source for full Copyright notice */ 
/* $Id$ */

#include <iostream>
using namespace std;

#include "AliAODTrack.h"
#include "AliAODEvent.h"
#include "AliAODMCParticle.h"
#include "AliPIDResponse.h"

class AliTrackDiHadronPID : public TObject {

public:
	AliTrackDiHadronPID();
	AliTrackDiHadronPID(AliAODTrack* track, AliAODTrack* globaltrack = 0x0, AliAODMCParticle* mcparticle = 0x0, AliPIDResponse* pidresponse = 0x0);
	virtual ~AliTrackDiHadronPID() {}

private:
	AliTrackDiHadronPID(const AliTrackDiHadronPID&);
	AliTrackDiHadronPID& operator=(const AliTrackDiHadronPID&);

// Internal copy functions.
private:
	Bool_t CopyBasicTrackInfo();
	Bool_t CopyFlags();
	Bool_t CopyDCAInfo();
	Bool_t CopyITSInfo();
	Bool_t CopyTPCInfo();
	Bool_t CopyTOFInfo();
	Bool_t CopyMCInfo();

// Check functions.
public:
	Bool_t IsBasicTrackInfoAvailable() const {return fBasicInfoAvailable;}
	Bool_t IsFlagInfoAvailable() const {return fFlagsAvailable;}
	Bool_t IsDCAInfoAvailable() const {return fDCAInfoAvailable;}
	Bool_t IsITSInfoAvailable() const {return fITSInfoAvailable;}
	Bool_t IsTPCInfoAvailable() const {return fTPCInfoAvailable;}
	Bool_t IsTOFInfoAvailable() const {return fTOFInfoAvailable;}
	Bool_t IsMCInfoAvailable() const {return fMCInfoAvailable;}			

public:
// Getting Track Parameters. Functionality is the same as AOD track,
// unless stated otherwise.
	Bool_t UnknownSpecies(Int_t species) const;
	Double_t Pt() const {return fPt;}
	Double_t Eta() const {return fEta;}
	Double_t Phi() const {return fPhi;}
	Double_t Y(Int_t species) {
		if (UnknownSpecies(species)) return fY[species];
		else return -999.;
	}

	ULong64_t GetFlags() const {return fFlags;}
	ULong64_t GetStatus() const {return GetFlags();}
	UInt_t GetFilterMap() const {return fFilterMap;}
	Bool_t TestFilterMask(UInt_t filterMask) const {return (Bool_t)((filterMask & fFilterMap) == filterMask);}

	Int_t GetID() const {return fID;}
	Int_t GetLabel() const {return fLabel;}
	void GetTOFLabel(Int_t* label) const {label[0] = fTOFLabel[0]; label[1] = fTOFLabel[1]; label[2] = fTOFLabel[2];}
	Short_t Charge() const {return fCharge;}
	Int_t GetNclsTPC() const {return fNclsTPC;}

	Double_t GetZAtDCA() const {return fDCAz;}
	Double_t GetXYAtDCA() const {return fDCAxy;}

	// TOF related Getters.
	Double_t GetTOFsignal() const {return fTOFsignal;}
	Double_t GetTOFsignalMinusExpected(Int_t species) const {
		if (UnknownSpecies(species)) {return -10e10;}
		return fTOFsignalMinusExpected[species];
	}
	Double_t GetTOFsignalExpected(Int_t species) const {
		if (UnknownSpecies(species)) {return -10e10;}
		return (fTOFsignal - fTOFsignalMinusExpected[species]);
	}
	Double_t GetTOFsigmaExpected(Int_t species) {
		if (UnknownSpecies(species)) {return -10e10;}
		if (GetNumberOfSigmasTOF(species) < 10e-10) {return -10e30;}
		return (GetTOFsignalMinusExpected(species)/GetNumberOfSigmasTOF(species));
	}
	Double_t GetNumberOfSigmasTOF(Int_t species) const {
		if (UnknownSpecies(species)) {return -10e10;}
		return fTOFNsigma[species];
	}
	Int_t GetTOFMatchingStatus() const {return fTOFMatchingStatus;}
	Bool_t IsTOFMismatch() const {
		if (fTOFMatchingStatus==1) {return kTRUE;}
		else {return kFALSE;}
	}

	Double_t GetTPCsignal() const {return fTPCsignal;}
	Double_t GetTPCsignalMinusExpected(Int_t species) const {return fTPCsignalMinusExpected[species];}
	Double_t GetTPCmomentum() const {return fTPCmomentum;}
	Double_t GetNumberOfSigmasTPC(Int_t species) const {
		if (UnknownSpecies(species)) {return -10e10;}
		return fTPCNsigma[species];
	}

	Char_t GetITSClusterMap() const {return fITSClusterMap;}
	Bool_t HasPointOnITSLayer(Int_t layer) const {
		if ((layer > -1) && (layer < 6)) return fITSHits[layer];
		else {
			cout<<"ITS has only 6 layers."<<endl;
			return kFALSE;
		}
	}

	Double_t MCPt() const {return fMCPt;}
	Double_t MCEta() const {return fMCEta;}
	Double_t MCPhi() const {return fMCPhi;}
	Double_t MCY() const {return fMCY;}
	Int_t GetMCSpecies() const {
		Int_t abspdg = TMath::Abs(GetPdgCode());
		if (abspdg == 211) return 0;
		else if (abspdg == 321) return 1;
		else if (abspdg == 2212) return 2;
		else return -999;
	}
	
	Int_t GetPdgCode() const {return fPdgCode;}
	Bool_t IsPhysicalPrimary() const {return fIsPhysicalPrimary;}
	Bool_t IsSecondaryFromWeakDecay() const {return fIsSecondaryFromWeakDecay;}
	Bool_t IsSecondaryFromMaterial() const {return fIsSecondaryFromMaterial;}

	void ForgetAboutPointers() {
		// AOD tracks are usually deleted, while the 
		// AliTrackDiHadronPID is not. This method ensures that
		// the pointers to the track/event objects etc. don't
		// point to deleted objects.
		fAODTrack = 0x0;
		fAODGlobalTrack = 0x0;
		fAODEvent = 0x0;
		fAODMCParticle = 0x0;
		fPIDResponse = 0x0;
	}

	void SetDebugLevel(Int_t debuglevel) {fDebug = debuglevel;}
	Int_t GetDebugLevel() const {return fDebug;}

private:
// Pointers to original Tracks etc (cannot be returned, will not be streamed/saved).
	AliAODTrack*		fAODTrack;			//! Original AOD Track.
	AliAODTrack*		fAODGlobalTrack;	//! Corresponding Global AOD Track.
	AliAODEvent*		fAODEvent;			//!	Original AOD Event.
	AliAODMCParticle*	fAODMCParticle;		//! Original MC Particle.
	AliPIDResponse*		fPIDResponse;		//!	Original PID Response.

// Flags if a certain type of information is available for this track.
	Bool_t				fBasicInfoAvailable;// Basic Track Info.
	Bool_t				fFlagsAvailable;
	Bool_t				fDCAInfoAvailable;	// DCA Info.
	Bool_t				fITSInfoAvailable;	// ITS Info.
	Bool_t				fTPCInfoAvailable;	// TPC Info.
	Bool_t 				fTOFInfoAvailable;	// TOF Info.
	Bool_t				fMCInfoAvailable;	// MC Info.

// Basic Track Info.
	Double_t 			fPt;				// Reconstructed Pt.
	Double_t			fEta;				// Reconstructed Eta.
	Double_t			fY[3];				// Reconstructed Rapidity for pi,K,p.
	Double_t			fPhi;				// Reconstructed Phi.

	ULong64_t			fFlags;				// Reconstruction Flags.
	UInt_t				fFilterMap;			// FilterMap.

	Short_t				fID;				// Track ID.
	Int_t 				fLabel;				// Track Label.
	Int_t				fTOFLabel[3];		// Track TOF label.

	Short_t				fCharge;			// Charge (is a Char_t in AliAODTrack)
	Int_t				fNclsTPC;			// Number of clusters in TPC.
	
// DCA Info.
	Double_t			fDCAz;				// z at DCA.
	Double_t			fDCAxy;				// xy at DCA.

// PID Info.
	Double_t			fTOFsignal;			
	Double_t			fTOFsignalMinusExpected[3];
	Double_t			fTOFNsigma[3];
	Int_t				fTOFMatchingStatus;		// 0 -> match, 1 -> mismatch, 2 -> no TOF hit.
	Double_t			fTPCsignal;			
	Double_t			fTPCsignalMinusExpected[3];
	Double_t			fTPCNsigma[3];
	Double_t			fTPCmomentum;		
	UChar_t				fITSClusterMap;
	Bool_t				fITSHits[6];

// MC Info.
	Double_t			fMCPt;
	Double_t			fMCEta;
	Double_t			fMCPhi;
	Double_t			fMCY;
	Int_t				fPdgCode;
	Bool_t				fIsPhysicalPrimary;
	Bool_t				fIsSecondaryFromWeakDecay;
	Bool_t				fIsSecondaryFromMaterial;

// Static variables.
public:
	static Double_t 	fSigmaTOFStd;
	static Double_t 	fSigmaTPCStd;

// Debug.
private:
	Int_t				fDebug;				// Debug flag.

	ClassDef(AliTrackDiHadronPID,2);

};

#endif
