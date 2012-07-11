#ifndef AliHFAssociatedTrackCuts_H
#define AliHFAssociatedTrackCuts_H
/**************************************************************************
 * Copyright(c) 1998-2010, ALICE Experiment at CERN, All rights reserved. *
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

/* $Id$ */

////////////////////////////////////////////////////////////////////////
//
// Base class for cuts on Associated tracks for HF Correlation analysis
//
// Author: S.Bjelogrlic (Utrecht) sandro.bjelogrlic@cern.ch
////////////////////////////////////////////////////////////////////////

#include <TString.h>
#include "AliAnalysisCuts.h"
#include "AliESDtrackCuts.h"
#include "AliAODPidHF.h"
#include "AliAODEvent.h"
#include "AliAODRecoDecayHF2Prong.h"
#include <TClonesArray.h>

class AliAODTrack;
class AliAODEvent;


//
class AliHFAssociatedTrackCuts : public AliAnalysisCuts
{
	public:	
	AliHFAssociatedTrackCuts();
	AliHFAssociatedTrackCuts(const char* name, const char* title);
	
	
	AliHFAssociatedTrackCuts(const AliHFAssociatedTrackCuts& source);
	AliHFAssociatedTrackCuts& operator=(const AliHFAssociatedTrackCuts& source);
	
	virtual ~AliHFAssociatedTrackCuts(); // destructor
	Bool_t IsSelected(TList*  list) {if(list) return kTRUE; return kFALSE;};
	Bool_t IsSelected(TObject*  obj) {if(obj) return kTRUE; return kFALSE;};
	Bool_t IsInAcceptance();
	Bool_t IsHadronSelected(AliAODTrack * track);
	Bool_t CheckHadronKinematic(Double_t pt, Double_t d0); 
	Bool_t Charge(Short_t charge, AliAODTrack* track);
	Bool_t CheckKaonCompatibility(AliAODTrack * track, Bool_t useMc, TClonesArray* mcArray, Int_t method=1);
	Bool_t IsKZeroSelected(AliAODv0 *vzero, AliAODVertex *vtx1);
	Int_t IsMCpartFromHF(Int_t label, TClonesArray*mcArray);
	Bool_t InvMassDstarRejection(AliAODRecoDecayHF2Prong* d, AliAODTrack *track, Int_t hypD0) const;	
	
	//getters
	Int_t GetMaxNEventsInPool() {return fPoolMaxNEvents;}
	Int_t GetMinNTracksInPool() {return fPoolMinNTracks;}
	Int_t GetMinEventsToMix(){return fMinEventsToMix;}
	Int_t GetNZvtxPoolBins() {return fNzVtxBins;}
	Double_t *GetZvtxPoolBins(){return fZvtxBins;}
	Int_t GetNCentPoolBins() {return fNCentBins;}
	Double_t *GetCentPoolBins(){return fCentBins;}
	
	
	
	void AddTrackCuts(const AliESDtrackCuts *cuts) {
		delete fESDTrackCuts; 
		fESDTrackCuts=new AliESDtrackCuts(*cuts); 
		return;
	}
	//setters
	//event pool settings
	void SetMaxNEventsInPool(Int_t events){fPoolMaxNEvents=events;}
	void SetMinNTracksInPool(Int_t tracks){fPoolMinNTracks=tracks;}
	void SetMinEventsToMix(Int_t events){fMinEventsToMix=events;}
	
	void SetNofPoolBins(Int_t Nzvtxbins, Int_t Ncentbins){
		fNzVtxBins=Nzvtxbins;
		fNzVtxBinsDim=Nzvtxbins+1;
		
	    fNCentBins=Ncentbins;
		fNCentBinsDim=Ncentbins+1;
	}
	
	void SetPoolBins(Double_t *ZvtxBins, Double_t* CentBins){
		fZvtxBins=ZvtxBins; 
		fCentBins=CentBins;
	}
	//cut settings
	void SetAODTrackCuts(Float_t *cutsarray);
	void SetTrackCutsNames(/*TString *namearray*/);
	void SetAODvZeroCuts(Float_t *cutsarray);
	void SetvZeroCutsNames(/*TString *namearray*/);
	void SetPidHF(AliAODPidHF* pid) {fPidObj = pid; return;}
	void SetCharge(Short_t charge) {fCharge = charge;}
	void SetFilterBit(Int_t bit) {fBit = bit;}
	virtual void PrintAll();
	virtual void PrintPoolParameters();
	Int_t GetNVarsTrack(){return fNTrackCuts;}
	
	
	
	void SetNVarsTrack(Int_t nVars){fNTrackCuts=nVars;}
	void SetNVarsVzero(Int_t nVars){fNvZeroCuts=nVars;}
	
private:
	//AliESDtrackCuts *fTrackCuts;
	AliESDtrackCuts *fESDTrackCuts; // track cut object
	AliAODPidHF * fPidObj;     /// PID object
	
	Int_t fPoolMaxNEvents; // set maximum number of events in the pool
	Int_t fPoolMinNTracks; // se minimum number of tracks in the pool
	Int_t fMinEventsToMix; // set the minimum number of events you wanna mix
	
	Int_t fNzVtxBins; // number of z vrtx bins
	Int_t fNzVtxBinsDim; // number of z vrtx bins +1 : necessary to initialize correctly the array
	Double_t* fZvtxBins; // [fNzVtxBinsDim]
	
	
	Int_t fNCentBins; //number of centrality bins
	Int_t fNCentBinsDim; //number of centrality bins bins +1 : necessary to initialize correctly the array
	Double_t* fCentBins; // [fNCentBinsDim]
	
	Int_t fNTrackCuts;     // array dimension
	Float_t* fAODTrackCuts;//[fNTrackCuts]
	TString * fTrackCutsNames;//[fNTrackCuts]
	Int_t fNvZeroCuts;// array dimension
	Float_t *fAODvZeroCuts;//[fNvZeroCuts]
	TString * fvZeroCutsNames;//[fNvZeroCuts]
	Int_t fBit; // filterBit
	Short_t fCharge; // charge (+1 or -1)
	
	
	ClassDef(AliHFAssociatedTrackCuts,3);
};


#endif
