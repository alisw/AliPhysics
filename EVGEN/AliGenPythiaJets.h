#ifndef ALIGENPYTHIAJETS_H
#define ALIGENPYTHIAJETS_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//
// Generator using the TPythia interface (via AliPythia)
// to generate jets in pp collisions.
//
// andreas.morsch@cern.ch
//

#include "AliGenMC.h"
#include "AliGenPythia.h"

class TParticle;

class AliGenPythiaJets : public AliGenPythia
{
 public:
    AliGenPythiaJets();
    AliGenPythiaJets(Int_t npart);
    AliGenPythiaJets(const AliGenPythiaJets &Pythia);
    virtual ~AliGenPythiaJets();
    virtual void    Init();
    virtual void    Generate();
    virtual void    TransformEvent(Float_t beta, Float_t gamma);
    virtual Bool_t  CheckTrigger();
    virtual void    SetQuenchingFactor(Float_t quench = -1) {fQuench = quench;}
    
    // Assignment Operator
    AliGenPythiaJets & operator=(const AliGenPythiaJets & rhs);
 protected:
    Float_t fQuench;                // Quench factor
    Float_t fEtMinJetQ[2];          // Minimum et of triggered Jet
    Float_t fEtMaxJetQ[2];          // Maximum et of triggered Jet
    Float_t fPtHardMinQ[2];         // Lower pT-hard cut 
    Float_t fPtHardMaxQ[2];         // Higher pT-hard cut
    ClassDef(AliGenPythiaJets,1)    // AliGenerator Interface to Pythia Jet Production
};
#endif





