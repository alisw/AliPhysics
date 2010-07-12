#ifndef _ALIGENDEUTERON_H_
#define _ALIGENDEUTERON_H_

/* Copyright(c) 2009-2010, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// Afterburner to simulate coalescence of nucleons
// Author: Eulogio Serradilla <eulogio.serradilla@ciemat.es>
//         Arturo Menchaca <menchaca@fisica.unam.mx>

#include "AliGenerator.h"

class AliGenDeuteron: public AliGenerator
{

 public:

	AliGenDeuteron(Int_t sign=1, Double_t pmax=0.3, Double_t rmax=2.1, Double_t rsrc=1.5, Int_t model=1);
	virtual ~AliGenDeuteron();

	virtual void Init();
	virtual void Generate();
	
	Double_t GetCoalescenceDistance() const { return fRmax*1.e+13; }
	Double_t GetCoalescenceMomentum() const { return fPmax; }
	Double_t GetSourceRadius() const { return fRsrc*1.e+13; }
	Double_t GetSpinProbability() const { return fSpinProb; }
	Int_t GetFreezeOutModel() const { return fModel; }
	Int_t GetSign() const { return fSign;}
	
	void SetCoalescenceDistance(Double_t r=2.1) { fRmax = r*1.e-13; }
	void SetCoalescenceMomentum(Double_t p=0.3) { fPmax = p; }
	void SetSourceRadius(Double_t r=1.5) { fRsrc = r*1.e-13; }
	void SetSpinProbability(Double_t s=0.75) { fSpinProb = s; }
	void SetFreezeOutModel(Int_t model = kThermal) { fModel = model; }
	void SetSign(Int_t sign) {  fSign = sign > 0 ? 1 : -1;}
	
 public:
	enum { kExpansion, kThermal };
	
 private:
 
	AliGenDeuteron(const AliGenDeuteron &other);
	AliGenDeuteron& operator=(const AliGenDeuteron &other);
	
	void FixProductionVertex(class TParticle* i);
	
 private:
	
	const Double_t fDeuteronMass;
	Double_t fPmax; // Maximum p-n momentum difference (GeV/c)
	Double_t fRmax; // Maximum p-n distance (cm)
	Double_t fRsrc; // Emitting source radius (cm)
	Double_t fSpinProb; // cluster formation probability due to spin
	Int_t fModel;   // How to place the nucleons at freeze-out stage
	Int_t fSign; // +1 for deuteron, -1 for antideuterons
	
	ClassDef(AliGenDeuteron,1)
};

#endif // _ALIGENDEUTERON_H_
