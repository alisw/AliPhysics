/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               *
 **************************************************************************/

//-------------------------------------------------------------------------
//                      Class AliRsnEvent
//  Simple collection of reconstructed tracks, selected from an ESD event
// 
// author: A. Pulvirenti             (email: alberto.pulvirenti@ct.infn.it)
//-------------------------------------------------------------------------

#ifndef ALIRSNEVENT_H
#define ALIRSNEVENT_H

#include "AliPID.h"

class AliRsnDaughter;

class AliRsnEvent : public TObject
{
public:
	               AliRsnEvent();
	               AliRsnEvent(const AliRsnEvent& copy);
				   AliRsnEvent& operator=(const AliRsnEvent& copy);
			 
	virtual       ~AliRsnEvent() {Clear("DELETE");}
	
	void           AddTrack(AliRsnDaughter track);
	void           Clear(Option_t *option = "");
	Int_t          GetMultiplicity(Bool_t recalc = kFALSE);
	const char*    GetOriginFileName() const;
	Double_t       GetPrimaryVertexX() const {return fPVx;}
	Double_t       GetPrimaryVertexY() const {return fPVy;}
	Double_t       GetPrimaryVertexZ() const {return fPVz;}
	void           GetPrimaryVertex(Double_t &x, Double_t &y, Double_t &z) const {x=fPVx;y=fPVy;z=fPVz;}
	TClonesArray*  GetTracks(Char_t sign, AliPID::EParticleType type);
	void           Init();
	Int_t          PDG2Enum(Int_t pdgcode);
	void           PrintTracks();
	void           SetESD(Bool_t yesno = kTRUE) {fIsESD=yesno;}
	void           SetPath(TString path) {fPath=path;}
	void           SetPrimaryVertex(Double_t x, Double_t y, Double_t z) {fPVx=x;fPVy=y;fPVz=z;}

private:

	Bool_t         fIsESD;                    // if true, it is ESD event, otherwise it comes from Kine
	TString        fPath;                     // complete path where input event file is stored
	
	Double_t       fPVx;  			          // position of
	Double_t       fPVy;  			          // primary
	Double_t       fPVz;  			          // vertex
	
	Int_t          fMultiplicity;             // global event multiplicity

	TClonesArray  *fPos[AliPID::kSPECIES];    // collections of positive particles
	TClonesArray  *fNeg[AliPID::kSPECIES];    // collections of negative particles
	
	ClassDef(AliRsnEvent,1);
};

#endif
