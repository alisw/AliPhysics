#ifndef ALIGENLIGHTNUCLEI_H
#define ALIGENLIGHTNUCLEI_H

/* Copyright(c) 2009-2010, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// afterburner to generate light nuclei
// Author: Eulogio Serradilla <eulogio.serradilla@cern.h>

#include "AliGenCocktail.h"

class TParticle;
class AliStack;
class TLorentzVector;

class AliGenLightNuclei: public AliGenCocktail
{

 public:

	AliGenLightNuclei();
	virtual ~AliGenLightNuclei();

	virtual void Generate();
	
	Double_t GetCoalescenceMomentum() const { return fP0; }
	
	void SetCoalescenceMomentum(Double_t p0) { fP0 = p0; }
	
	void SetNucleusPdgCode(Int_t pdg) { fPdg = pdg; }
	
	enum {kDeuteron=1000010020, kAntiDeuteron=-1000010020, kTriton=1000010030, kAntiTriton=-1000010030, kHe3Nucleus=1000020030, kAntiHe3Nucleus=-1000020030, kAlpha=1000020040, kAntiAlpha=-1000020040, kHyperTriton=1010010030, kAntiHyperTriton=-1010010030 };

	enum {kCluster=77};
	
 private:
 
	AliGenLightNuclei(const AliGenLightNuclei &other);
	AliGenLightNuclei& operator=(const AliGenLightNuclei &other);
	
	Bool_t Coalescence(const TLorentzVector& p1, const TLorentzVector& p2) const;
	Bool_t Coalescence(const TLorentzVector& p1, const TLorentzVector& p2, const TLorentzVector& p3) const;
	Bool_t Coalescence(const TLorentzVector& p1, const TLorentzVector& p2, const TLorentzVector& p3, const TLorentzVector& p4) const;
	
	Int_t GenerateNuclei(Int_t pdg, Double_t mass, const TList* l1, const TList* l2);
	Int_t GenerateNuclei(Int_t pdg, Double_t mass, const TList* l1, const TList* l2, const TList* l3);
	Int_t GenerateNuclei(Int_t pdg, Double_t mass, const TList* l1, const TList* l2, const TList* l3, const TList* l4);
	
	void PushNucleus(Int_t pdg, Double_t mass, TParticle* parent1, TParticle* parent2, TParticle* parent3=0, TParticle* parent4=0);
	
 private:
	
	Double_t fP0; // coalescence momentum (radius of the sphere)
	Int_t fPdg;   // nucleus PDG code
	
	ClassDef(AliGenLightNuclei,2)
};

#endif // ALIGENLIGHTNUCLEI_H
