#ifndef ALIGENLIGHTNUCLEI_H
#define ALIGENLIGHTNUCLEI_H

/* Copyright(c) 2009-2010, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// afterburner to generate light nuclei
// Author: Eulogio Serradilla <eulogio.serradilla@cern.h>

#include "AliGenCocktail.h"

class TParticle;
class AliStack;

class AliGenLightNuclei: public AliGenCocktail
{

 public:

	AliGenLightNuclei();
	virtual ~AliGenLightNuclei();

	virtual void Generate();
	
	Double_t GetCoalescenceMomentum() const { return fP0; }
	
	void SetCoalescenceMomentum(Double_t p0) { fP0 = p0; }
	
	void SetDeuterons(Bool_t flag = kTRUE) { fGenDeuterons = flag; }
	void SetTritons(Bool_t flag = kTRUE)   { fGenTritons = flag; }
	void SetHe3Nuclei(Bool_t flag = kTRUE) { fGenHe3Nuclei = flag; }
	void SetHe4Nuclei(Bool_t flag = kTRUE) { fGenHe4Nuclei = flag; }
	
	enum {kDeuteron=1000010020, kAntiDeuteron=-1000010020, kTriton=1000010030, kAntiTriton=-1000010030, kHe3Nucleus=1000020030, kAntiHe3Nucleus=-1000020030, kAlpha=1000020040, kAntiAlpha=-1000020040 };

	enum {kCluster=77};
	
 private:
 
	AliGenLightNuclei(const AliGenLightNuclei &other);
	AliGenLightNuclei& operator=(const AliGenLightNuclei &other);
	
	Int_t GenerateDeuterons(const TList* protons, const TList* neutrons);
	Int_t GenerateTritons(const TList* protons, const TList* neutrons);
	Int_t GenerateHe3Nuclei(const TList* protons, const TList* neutrons);
	Int_t GenerateHe4Nuclei(const TList* protons, const TList* neutrons);
	
	void PushDeuteron(TParticle* parent1, TParticle* parent2);
	void PushTriton(TParticle* parent1, TParticle* parent2, TParticle* parent3);
	void PushHe3Nucleus(TParticle* parent1, TParticle* parent2, TParticle* parent3);
	void PushHe4Nucleus(TParticle* parent1, TParticle* parent2, TParticle* parent3, TParticle* parent4);
	
	void PushNucleus(Int_t pdg, Double_t mass, TParticle* parent1, TParticle* parent2, TParticle* parent3=0, TParticle* parent4=0);
	
	Bool_t Coalescence(const TParticle* n1, const TParticle* n2) const;
	TParticle* FindPartner(const TParticle* n0, const TList* nucleons, const TParticle* nx=0) const;
	
	Double_t GetS(Double_t p1x, Double_t p1y, Double_t p1z, Double_t m1, Double_t p2x, Double_t p2y, Double_t p2z, Double_t m2) const;
	Double_t GetPcm(Double_t p1x, Double_t p1y, Double_t p1z, Double_t m1, Double_t p2x, Double_t p2y, Double_t p2z, Double_t m2) const;
	
 private:
	
	Double_t fP0;         // coalescence momentum (radius of the sphere)
	Bool_t fGenDeuterons; // generate deuterons and anti-deuterons
	Bool_t fGenTritons;   // generate tritons and anti-tritons
	Bool_t fGenHe3Nuclei; // generate He3 and anti-He3 nuclei
	Bool_t fGenHe4Nuclei; // generate He4 and anti-He4 nuclei
	
	ClassDef(AliGenLightNuclei,1)
};

#endif // ALIGENLIGHTNUCLEI_H
