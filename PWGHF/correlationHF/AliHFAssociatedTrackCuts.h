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
	Bool_t IsHadronSelected(AliAODTrack * track, AliAODVertex *vtx1, Double_t bz);
	Bool_t CheckKaonCompatibility(AliAODTrack * track, Bool_t useMc, TClonesArray* mcArray);
	Bool_t IsKZeroSelected(AliAODv0 *vzero, AliAODVertex *vtx1);
	Int_t IsMCpartFromHF(Int_t label, TClonesArray*mcArray);
	
	
	
	void AddTrackCuts(const AliESDtrackCuts *cuts) {
		delete fESDTrackCuts; 
		fESDTrackCuts=new AliESDtrackCuts(*cuts); 
		return;
	}
	
	void SetAODTrackCuts(Float_t *cutsarray);
	void SetTrackCutsNames(/*TString *namearray*/);
	void SetAODvZeroCuts(Float_t *cutsarray);
	void SetvZeroCutsNames(/*TString *namearray*/);
	void SetPidHF(AliAODPidHF* pid) {fPidObj = pid; return;}
	virtual void PrintAll();
	Int_t GetNVarsTrack(){return fNTrackCuts;}
	
	
	
	void SetNVarsTrack(Int_t nVars){fNTrackCuts=nVars;}
	void SetNVarsVzero(Int_t nVars){fNvZeroCuts=nVars;}
	
private:
	//AliESDtrackCuts *fTrackCuts;
	AliESDtrackCuts *fESDTrackCuts; // track cut object
	AliAODPidHF * fPidObj;     /// PID object
	Int_t fNTrackCuts;     // array dimension
	Float_t* fAODTrackCuts;//[fNTrackCuts]
	TString * fTrackCutsNames;//[fNTrackCuts]
	Int_t fNvZeroCuts;// array dimension
	Float_t *fAODvZeroCuts;//[fNvZeroCuts]
	TString * fvZeroCutsNames;//[fNvZeroCuts]
	
	
	ClassDef(AliHFAssociatedTrackCuts,1);
};


#endif
