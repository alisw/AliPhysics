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
	Bool_t IsBasicTrackInfoAvailable() const { 
		if (!fBasicInfoAvailable) cout<<"Basic Track Info not available."<<endl; 
		return fBasicInfoAvailable;
	}
	Bool_t IsFlagInfoAvailable() const { 
		if (!fFlagsAvailable) cout<<"Flag Info not available."<<endl; 
		return fFlagsAvailable;
	}
	Bool_t IsDCAInfoAvailable() const { 
		if (!fDCAInfoAvailable) cout<<"DCA Info not available."<<endl; 
		return fDCAInfoAvailable;
	}
	Bool_t IsITSInfoAvailable() const { 
		if (!fITSInfoAvailable) cout<<"ITS Info not available."<<endl; 
		return fITSInfoAvailable;
	}
	Bool_t IsTPCInfoAvailable() const { 
		if (!fTPCInfoAvailable) cout<<"TPC Info not available."<<endl; 
		return fTPCInfoAvailable;
	}
	Bool_t IsTOFInfoAvailable() const { 
		if (!fTOFInfoAvailable) cout<<"TOF Info not available."<<endl; 
		return fTOFInfoAvailable;
	}
	Bool_t IsMCInfoAvailable() const { 
		if (!fMCInfoAvailable) cout<<"MC Info not available."<<endl; 
		return fMCInfoAvailable;
	}			

public:
// Getting Track Parameters. Functionality is the same as AOD track,
// unless stated otherwise.
	Double_t Pt() const {return fPt;}
	Double_t Eta() const {return fEta;}
	Double_t Phi() const {return fPhi;}
	Double_t Y(Int_t species) {
		if (species >= 0 && species < 3) return fY[species];
		else return -999.;
	}

	ULong_t GetFlags() const {return fFlags;}
	ULong_t GetStatus() const {return GetFlags();}
	UInt_t GetFilterMap() const {return fFilterMap;}
	Bool_t TestFilterMask(UInt_t filterMask) const {return (Bool_t)((filterMask & fFilterMap) == filterMask);}

	Int_t GetID() const {return fID;}
	Int_t GetLabel() const {return fLabel;}
	Short_t Charge() const {return fCharge;}

	Double_t GetZAtDCA() const {return fDCAz;}
	Double_t GetXYAtDCA() const {return fDCAxy;}

	Double_t GetTOFsignal() const {return fTOFsignal;}
	Double_t GetTOFsignalMinusExpected(Int_t species) const {
		if (species < 0 || species > 2) {
			cout<<"ERROR: Unknown species"<<endl;
			return -10e10;
		}
		return fTOFsignalMinusExpected[species];
	}
	Double_t GetTOFsignalExpected(Int_t species) const {
		if (species < 0 || species > 2) {
			cout<<"ERROR: Unknown species"<<endl;
			return -10e10;
		}
		return (fTOFsignal - fTOFsignalMinusExpected[species]);
	}
	Double_t GetNumberOfSigmasTOF(Int_t species) const {
		if (species < 0 || species > 2) {
			cout<<"ERROR: Unknown species"<<endl;
			return -10e10;
		}
		return fTOFNsigma[species];
	}
	Double_t GetNumberOfSigmasTPC(Int_t species) const {
		if (species < 0 || species > 2) {
			cout<<"ERROR: Unknown species"<<endl;
			return -10e10;
		}
		return fTPCNsigma[species];
	}
	Bool_t IsTOFmismatch() const {return fIsTOFmismatch;}

	Double_t GetTPCsignal() const {return fTPCsignal;}
	Double_t GetTPCsignalMinusExpected(Int_t species) const {return fTPCsignalMinusExpected[species];}
	Double_t GetTPCmomentum() const {return fTPCmomentum;}

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
	const AliAODEvent*	fAODEvent;			//!	Original AOD Event.
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

	ULong_t				fFlags;				// Reconstruction Flags.
	UInt_t				fFilterMap;			// FilterMap.

	Short_t				fID;				// Track ID.
	Int_t 				fLabel;				// Track Label.

	Short_t				fCharge;			// Charge (is a Char_t in AliAODTrack)
	
// DCA Info.
	Double_t			fDCAz;				// z at DCA.
	Double_t			fDCAxy;				// xy at DCA.

// PID Info.
	Double_t			fTOFsignal;
	Double_t			fTOFsignalMinusExpected[3];
	Double_t			fTOFNsigma[3];
	Bool_t				fIsTOFmismatch;
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

// Debug.
	Int_t				fDebug;				// Debug flag.

	ClassDef(AliTrackDiHadronPID,1);

};

#endif
