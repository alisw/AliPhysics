#ifndef ALIGENLIGHTNUCLEI_H
#define ALIGENLIGHTNUCLEI_H

/* Copyright(c) 2009-2016, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// afterburner to generate light nuclei
// Author: Eulogio Serradilla <eulogio.serradilla@cern.h>

#include "AliGenerator.h"

class TParticle;
class AliStack;
class TLorentzVector;

class AliGenLightNuclei: public AliGenerator
{

 public:

	AliGenLightNuclei();
	AliGenLightNuclei(Int_t pdg, Double_t p0) { fPdg = pdg; fP0 = p0; }
	
	virtual ~AliGenLightNuclei();
	
	virtual void Generate();
	
	Double_t GetCoalescenceMomentum() const { return fP0; }
	
	void SetNucleusPdgCode(Int_t pdg) { fPdg = pdg; }
	
	void SetCoalescenceMomentum(Double_t p0) { fP0 = p0; }
	
	enum {   kDeuteron    = 1000010020 // pn
	       , kLambdaN     = 1010000020 // ln
	       , kHDibarion   = 1020000020 // ll
	       , kTriton      = 1000010030 // pnn
	       , kHyperTriton = 1010010030 // lpn
	       , kLambdaNN    = 1010000030 // lnn
	       , kLambdaLN    = 1020000030 // lln
	       , kLambdaLP    = 1020010030 // llp
	       , kHe3Nucleus  = 1000020030 // pnp
	       , kAlpha       = 1000020040 // pnpn
	};
	
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
	
	Int_t fPdg;          // nucleus PDG code
	Double_t fP0;        // coalescence momentum
	
	ClassDef(AliGenLightNuclei, 7)
};

#endif // ALIGENLIGHTNUCLEI_H
