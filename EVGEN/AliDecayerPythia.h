#ifndef ALIDECAYERPYTHIA_H
#define ALIDECAYERPYTHIA_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "AliDecayer.h"
#include "AliPythia.h"

class AliDecayerPythia :
public AliDecayer
{
 public:
//
    AliDecayerPythia();
    virtual void    Init();
    virtual void    Decay(Int_t idpart, TLorentzVector *p);
    virtual Int_t   ImportParticles(TClonesArray *particles)
	{return fPythia->ImportParticles(particles, "All");}
    virtual void    SetForceDecay(Decay_t decay) {fDecay=decay;}
    virtual void    ForceDecay();
    
    virtual Float_t GetPartialBranchingRatio(Int_t ipart);
 private:
    void     DefineParticles();
    void     Lu1Ent(Int_t flag, Int_t idpart, 
		    Double_t mom, Double_t theta, Double_t phi);
    Int_t    CountProducts(Int_t channel, Int_t particle);
    void     ForceParticleDecay(Int_t particle, Int_t product, Int_t mult);
    void     ForceHadronicD();    
    void     AllowAllDecays();
    Float_t  GetBraPart(Int_t kf);
    
    
 private:
    AliPythia*  fPythia;          // Pointer to AliPythia
    Decay_t     fDecay;           // Forced decay mode
    Float_t     fBraPart[501];    // Branching ratios

    ClassDef(AliDecayerPythia,1)  // AliDecayer implementation using Pythia  
};
#endif







