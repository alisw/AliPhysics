#ifndef ALIGENHIJING_H
#define ALIGENHIJING_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */


#include "AliGenerator.h"
#include "GenTypeDefs.h"
#include <THijing.h>
#include <TString.h>
#include <TArrayI.h>

class TArrayI;
class TParticle;

class AliGenHijing : public AliGenerator
{
    enum {kNoTrigger, kHardProcesses, kDirectPhotons};

 public:
    AliGenHijing();
    AliGenHijing(Int_t npart);
    AliGenHijing(const AliGenHijing &Hijing);
    virtual ~AliGenHijing();
    virtual void    Generate();
    virtual void    Init();
    // set centre of mass energy
    virtual void    SetEnergyCMS(Float_t energy=5500) {fEnergyCMS=energy;}
    virtual void    SetReferenceFrame(TString frame="CMS")
	{fFrame=frame;}
    virtual void    SetProjectile(TString proj="A",
				  Int_t a=208, Int_t z=82)
	{fProjectile = proj;
	fAProjectile = a;
	fZProjectile = z;}    
    virtual void    SetTarget(TString tar="A",
			      Int_t a=208, Int_t z=82)
	{fTarget = tar;
	fATarget = a;
	fZTarget =z;}    
    virtual void    SetImpactParameterRange(Float_t bmin = 0, Float_t bmax =15.)
	{ fMinImpactParam=bmin;
	fMaxImpactParam=bmax;
	}
    virtual void    KeepFullEvent();
    virtual void    SetJetQuenching(Int_t flag=1)     {fQuench    = flag;}
    virtual void    SetShadowing(Int_t flag=1)        {fShadowing = flag;}
    virtual void    SetDecaysOff(Int_t flag=1)        {fDecaysOff = flag;}
    virtual void    SetTrigger(Int_t flag=kNoTrigger) {fTrigger   = flag;}
    virtual void    SetFlavor(Int_t flag=0)           {fFlavor    = flag;}    
    virtual void    SetEvaluate(Int_t flag=0)         {fEvaluate  = flag;}
    virtual void    SetSelectAll(Int_t flag=0)        {fSelectAll = flag;}    
    AliGenHijing &  operator=(const AliGenHijing & rhs);
// Physics Routines	    
    virtual void EvaluateCrossSections();
 protected:
    Bool_t SelectFlavor(Int_t pid);
    void   MakeHeader();

 protected:
    TString     fFrame;         // Reference frame 
    TString     fProjectile;    // Projectile
    TString     fTarget;        // Target
    Int_t       fAProjectile;    // Projectile A
    Int_t       fZProjectile;    // Projectile Z
    Int_t       fATarget;        // Target A
    Int_t       fZTarget;        // Target Z
    Float_t     fMinImpactParam; // minimum impact parameter
    Float_t     fMaxImpactParam; // maximum impact parameter	
    Int_t       fKeep;           // Flag to keep full event information
    Int_t       fQuench;         // Flag to switch on jet quenching
    Int_t       fShadowing;      // Flag to switch on voclear effects on parton distribution function
    Int_t       fDecaysOff;      // Flag to turn off decays of pi0, K_s, D, Lambda, sigma
    Int_t       fTrigger;        // Trigger type
    Int_t       fEvaluate;       // Evaluate total and partial cross-sections
    Int_t       fSelectAll;      // Flag to write the full event
    Int_t       fFlavor;         // Selected particle flavor 4: charm+beauty 5: beauty
    Decay_t     fForceDecay;     // Decay channel  are forced
    Float_t     fEnergyCMS;      // Centre of mass energy
    Float_t     fKineBias;       // Bias from kinematic selection
    Int_t       fTrials;         // Number of trials
    TArrayI     fParentSelect;   // Parent particles to be selected 
    TArrayI     fChildSelect;    // Decay products to be selected
    Float_t     fXsection;       // Cross-section
    THijing    *fHijing;         // Hijing
    Float_t     fPtHardMin;      // lower pT-hard cut 
    Float_t     fPtHardMax;      // higher pT-hard cut

 private:
    // check if particle is selected as parent particle
    Bool_t ParentSelected(Int_t ip);
    // check if particle is selected as child particle
    Bool_t ChildSelected(Int_t ip);
    // all kinematic selection cuts go here 
    Bool_t KinematicSelection(TParticle *particle);
    // adjust the weight from kinematic cuts
    void   AdjustWeights();
    // check seleted daughters
    Bool_t DaughtersSelection(TParticle* iparticle, TClonesArray* particles);
    ClassDef(AliGenHijing,1) // AliGenerator interface to Hijing
};
#endif





