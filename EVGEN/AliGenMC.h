#ifndef ALIGENMC_H
#define ALIGENMC_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */


#include "AliGenerator.h"
#include "AliDecayer.h"
#include <TArrayI.h>    

class TParticle;

class AliGenMC : public AliGenerator
{
 public:
    AliGenMC();
    AliGenMC(Int_t npart);
    AliGenMC(const AliGenMC &MC);
    virtual ~AliGenMC();
    virtual void Init();
    virtual void SetForceDecay(Decay_t decay = kAll) {fForceDecay = decay;}
    AliGenMC & operator=(const AliGenMC & rhs);
        virtual void SetCutOnChild(Int_t flag = 0) {fCutOnChild = flag;}
    virtual void SetChildMomentumRange(Float_t pmin = 0, Float_t pmax = 1.e10)
	{fChildPMin = pmin; fChildPMax = pmax;}
    virtual void SetChildPtRange(Float_t ptmin = 0, Float_t ptmax = 20.)
	{fChildPtMin = ptmin; fChildPtMax = ptmax;}
    virtual void SetChildPhiRange(Float_t phimin = -180., Float_t phimax = 180)
	{fChildPhiMin = TMath::Pi()*phimin/180;
	fChildPhiMax  = TMath::Pi()*phimax/180;}
    virtual void SetChildThetaRange(Float_t thetamin = 0, Float_t thetamax = 180)
	{fChildThetaMin = TMath::Pi()*thetamin/180;
	fChildThetaMax  = TMath::Pi()*thetamax/180;}
    virtual void SetChildYRange(Float_t ymin = -12, Float_t ymax = 12)
	{fChildYMin = ymin;
	fChildYMax  = ymax;}
  virtual void SetMaximumLifetime(Float_t time = 1.e-15) {fMaxLifeTime = time;}
 protected:
    // check if particle is selected as parent particle
    Bool_t ParentSelected(Int_t ip);
    // check if particle is selected as child particle
    Bool_t ChildSelected(Int_t ip);
    // all kinematic selection cuts go here 
    Bool_t KinematicSelection(TParticle *particle, Int_t flag);
    Int_t  CheckPDGCode(Int_t pdgcode);

 protected:
    TArrayI     fParentSelect;  //!Parent particles to be selected 
    TArrayI     fChildSelect;   //!Decay products to be selected
    Int_t       fCutOnChild;    // Cuts on decay products (children)  are enabled/disabled
    Float_t     fChildPtMin;    // Children minimum pT
    Float_t     fChildPtMax;    // Children maximum pT
    Float_t     fChildPMin;     // Children minimum p
    Float_t     fChildPMax;     // Children maximum p
    Float_t     fChildPhiMin;   // Children minimum phi
    Float_t     fChildPhiMax;   // Children maximum phi
    Float_t     fChildThetaMin; // Children minimum theta
    Float_t     fChildThetaMax; // Children maximum theta
    Float_t     fChildYMin;     // Children minimum y
    Float_t     fChildYMax;     // Children maximum y
    Decay_t     fForceDecay;    // Decay channel forced
    Float_t     fMaxLifeTime;   // Maximum lifetime for unstable particles
    
    ClassDef(AliGenMC,4)       // AliGenerator implementation for generators using MC methods
};
#endif





