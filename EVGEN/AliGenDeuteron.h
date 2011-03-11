#ifndef ALIGENDEUTERON_H
#define ALIGENDEUTERON_H


/* Copyright(c) 2009-2010, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// Afterburner to simulate coalescence of (anti)nucleons
// Author: Eulogio Serradilla <eulogio.serradilla@ciemat.es>
//         Arturo Menchaca <menchaca@fisica.unam.mx>

#include "AliGenerator.h"


class AliGenDeuteron: public AliGenerator
{

 public:

	AliGenDeuteron(Int_t sign=1, Double_t pmax=0.3, Double_t rmax=2.1, Int_t cluster=0 );
	virtual ~AliGenDeuteron();

	virtual void Init();
	virtual void Generate();
	
	Int_t GetSign() const { return fSign;}
	Double_t GetCoalescenceMomentum() const { return fPmax; }
	Double_t GetCoalescenceDistance() const { return fRmax; }
	Double_t GetSpinProbability() const { return fSpinProb; }
	Double_t GetFreezeOutRadius() const { return fMaxRadius; }
	Double_t GetCoalescenceProbability(const TParticle* nucleon1, const TParticle* nucleon2) const;
	Int_t GetClusterType() const { return fClusterType; }
	Int_t GetFreezeOutModel() const { return fModel; }
	Double_t GetFreezeOutTime() const { return fTimeLength; }
	
	void SetSign(Int_t sign) {  fSign = sign > 0 ? 1 : -1;}
	void SetCoalescenceMomentum(Double_t p) { fPmax = p; }
	void SetCoalescenceDistance(Double_t r) { fRmax = r; }
	void SetSpinProbability(Double_t s) { fSpinProb = s; }
	void SetFreezeOutRadius(Double_t r) { fMaxRadius = r; }
	void SetClusterType(Int_t cluster) { fClusterType = cluster; }
	void SetFreezeOutModel(Int_t model, Double_t timeLength=2.5) { fModel = model; fTimeLength=timeLength;}
	
 public:

	enum { kFirstPartner, kLowestMomentum, kLowestDistance, kBoth };
	enum { kNone, kThermal, kExpansion };
	
 private:
 
	AliGenDeuteron(const AliGenDeuteron &other);
	AliGenDeuteron& operator=(const AliGenDeuteron &other);
	
	void FixProductionVertex(class TParticle* i);
	void FirstPartner(const TList* protons, TList* neutrons);
	void WeightMatrix(const TList* protons, const TList* neutrons);
	void PushDeuteron(TParticle* parent1, TParticle* parent2);
	
 private:
	
	Int_t fSign;          // +1 for deuterons, -1 for antideuterons
	Double_t fPmax;       // Maximum p-n momentum difference (GeV/c)
	Double_t fRmax;       // Maximum p-n distance (fm)
	Double_t fSpinProb;   // cluster formation probability due to spin
	Double_t fMaxRadius;  // Maximum freeze-out radius (fm)
	Int_t fClusterType;   // Probability criteria to find clusters
	Int_t fModel;         // Model to override generator source
	Double_t fTimeLength; // Thermal and chemical freeze-out time (fm/c)
	Double_t fB;          // Impact parameter (fm)
	Double_t fR;          // Projectile/Target nuclear radius (fm)
	Double_t fPsiR;       // Reaction plane angle
	AliStack* fCurStack;  //! current event stack
	Int_t fNtrk;          //! number of the stored track
	
	ClassDef(AliGenDeuteron,1)
};

#endif // ALIGENDEUTERON_H
