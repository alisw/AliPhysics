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
#include "AliESDVertex.h"
#include "AliAODPidHF.h"
#include "AliAODEvent.h"
#include "AliAODRecoDecayHF2Prong.h"
#include <TClonesArray.h>
#include <TH3D.h>


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
	Bool_t IsHadronSelected(AliAODTrack * track,const AliESDVertex *primary=0x0,const Double_t magfield=0);
	Bool_t CheckHadronKinematic(Double_t pt, Double_t d0); 
	Bool_t Charge(Short_t charge, AliAODTrack* track);
	Bool_t CheckKaonCompatibility(AliAODTrack * track, Bool_t useMc, TClonesArray* mcArray, Int_t method=1);
	Bool_t IsKZeroSelected(AliAODv0 *vzero, AliAODVertex *vtx1);
	Bool_t *IsMCpartFromHF(Int_t label, TClonesArray*mcArray);
	Bool_t InvMassDstarRejection(AliAODRecoDecayHF2Prong* d, AliAODTrack *track, Int_t hypD0) const;
	void SetPidAssociated();	
	
	// getters
	AliESDtrackCuts * GetESDTrackCuts() const {return fESDTrackCuts;}
	AliAODPidHF * GetPIDObject() const {return fPidObj;}
	TH3D * GetEfficiencyWeight() const {return fEffWeights;}
	
	Int_t GetMaxNEventsInPool() const {return fPoolMaxNEvents;}
	Int_t GetMinNTracksInPool() const {return fPoolMinNTracks;}
	Int_t GetMinEventsToMix() const {return fMinEventsToMix;}
	Int_t GetNZvtxPoolBins() const {return fNzVtxBins;}
	Double_t *GetZvtxPoolBins() const {return fZvtxBins;}
	Int_t GetNCentPoolBins() const {return fNCentBins;}
	Double_t *GetCentPoolBins() const {return fCentBins;}
	
	Int_t GetNofMCEventType() const {return fNofMCEventType;}
	Int_t *GetMCEventType() const {return fMCEventType;}
	
	Int_t GetNTrackCuts() const {return fNTrackCuts;}
	Float_t* GetAODTrackCuts() const {return fAODTrackCuts;}
	TString * GetTrackCutNames() const {return fTrackCutsNames;}
	Int_t GetNvZeroCuts() const {return fNvZeroCuts;}
	Float_t * GetAODvZeroCuts() const {return fAODvZeroCuts;}
	TString * GetvZeroCutNames() const {return fvZeroCutsNames;}
	Int_t GetFilterBit() const {return fBit;}
	Short_t GetCharge() const {return fCharge;}
	TString GetDescription() const {return fDescription;}
	
	
	
	void AddTrackCuts(const AliESDtrackCuts *cuts) {
		delete fESDTrackCuts; 
		fESDTrackCuts=new AliESDtrackCuts(*cuts); 
		return;
	}
	
	void AddDescription(TString description){fDescription=description;}
	
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
	
	// set MC events to process
	
	void SetNofMCEventTypes(Int_t k) {fNofMCEventType=k;}
	void SetMCEventTypes(Int_t *MCEventTypeArray);
	
	//cut settings
	void SetAODTrackCuts(Float_t *cutsarray);
	void SetTrackCutsNames(/*TString *namearray*/);
	void SetAODvZeroCuts(Float_t *cutsarray);
	void SetvZeroCutsNames(/*TString *namearray*/);
	void SetPidHF(AliAODPidHF* pid) {fPidObj = pid; return;}
	void SetCharge(Short_t charge) {fCharge = charge;}
	void SetFilterBit(Int_t bit) {fBit = bit;}
	void SetEfficiencyWeightMap(TH3D *hMap){if(fEffWeights)delete fEffWeights;fEffWeights=(TH3D*)hMap->Clone();}
	Double_t GetTrackWeight(Double_t pt, Double_t eta,Double_t zvtx);
	void Print(Option_t *option) const;
	virtual void PrintAll() const;
	virtual void PrintPoolParameters() const;
	virtual void PrintSelectedMCevents() const;

	
	
	
	void SetNVarsTrack(Int_t nVars){fNTrackCuts=nVars;}
	void SetNVarsVzero(Int_t nVars){fNvZeroCuts=nVars;}
	
	
	
private:
	AliESDtrackCuts *fESDTrackCuts; // track cut object
	AliAODPidHF * fPidObj;     /// PID object
	TH3D *fEffWeights;     // weight map (pt,eta,zvtx) to account for single track efficiency  
	Int_t fPoolMaxNEvents; // set maximum number of events in the pool
	Int_t fPoolMinNTracks; // se minimum number of tracks in the pool
	Int_t fMinEventsToMix; // set the minimum number of events you wanna mix
	
	Int_t fNzVtxBins; // number of z vrtx bins
	Int_t fNzVtxBinsDim; // number of z vrtx bins +1 : necessary to initialize correctly the array
	Double_t* fZvtxBins; // [fNzVtxBinsDim]
	
	
	Int_t fNCentBins; //number of centrality bins
	Int_t fNCentBinsDim; //number of centrality bins bins +1 : necessary to initialize correctly the array
	Double_t* fCentBins; // [fNCentBinsDim]
	
	Int_t fNofMCEventType;// number of event types to be selected in MC simultaneously;
	Int_t *fMCEventType;//[fNofMCEventType]
	
	Int_t fNTrackCuts;     // array dimension
	Float_t* fAODTrackCuts;//[fNTrackCuts]
	TString * fTrackCutsNames;//[fNTrackCuts]
	Int_t fNvZeroCuts;// array dimension
	Float_t *fAODvZeroCuts;//[fNvZeroCuts]
	TString * fvZeroCutsNames;//[fNvZeroCuts]
	Int_t fBit; // filterBit
	Short_t fCharge; // charge (+1 or -1)
	TString fDescription; // additional description to the cuts
	
	
	ClassDef(AliHFAssociatedTrackCuts,5);
};


#endif
