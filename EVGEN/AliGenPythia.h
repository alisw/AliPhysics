#ifndef ALIGENPYTHIA_H
#define ALIGENPYTHIA_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */


#include "AliGenerator.h"
#include "GenTypeDefs.h"
#include <TArrayI.h>    

class AliPythia;
class TParticle;

class AliGenPythia : public AliGenerator
{
 public:
    AliGenPythia();
    AliGenPythia(Int_t npart);
    AliGenPythia(const AliGenPythia &Pythia);
    virtual ~AliGenPythia();
    virtual void    Generate();
    virtual void    Init();
    // select process type
    virtual void    SetProcess(Process_t proc=charm) {fProcess=proc;}
    // select structure function
    virtual void    SetStrucFunc(StrucFunc_t func=GRV_HO) {fStrucFunc=func;}
    // select pt of hard scattering 
    virtual void    SetPtHard(Float_t ptmin=0, Float_t ptmax=1.e10)
	{fPtHardMin=ptmin; fPtHardMax=ptmax; }
    // set centre of mass energy
    virtual void    SetEnergyCMS(Float_t energy=5500) {fEnergyCMS=energy;}
    // force decay type
    virtual void    SetForceDecay(Decay_t decay=semimuonic) {fForceDecay=decay;}
    // get cross section of process
    virtual Float_t GetXsection() {return fXsection;}
    // Check PDG code
    virtual Int_t CheckPDGCode(Int_t pdgcode);
    // Assignment Operator
    AliGenPythia & operator=(const AliGenPythia & rhs);
 protected:
    Process_t   fProcess;       // Process type
    StrucFunc_t fStrucFunc;     // Structure Function
    Decay_t     fForceDecay;    // Decay channel  are forced
    Float_t     fEnergyCMS;     // Centre of mass energy
    Float_t     fKineBias;      // Bias from kinematic selection
    Int_t       fTrials;        // Number of trials
    TArrayI     fParentSelect;  // Parent particles to be selected 
    TArrayI     fChildSelect;   // Decay products to be selected
    Float_t     fXsection;      // Cross-section
    AliPythia   *fPythia;       // Pythia 
    Float_t     fPtHardMin;     // lower pT-hard cut 
    Float_t     fPtHardMax;     // higher pT-hard cut

 private:
    // check if particle is selected as parent particle
    Bool_t ParentSelected(Int_t ip);
    // check if particle is selected as child particle
    Bool_t ChildSelected(Int_t ip);
    // all kinematic selection cuts go here 
    Bool_t KinematicSelection(TParticle *particle);
    // adjust the weight from kinematic cuts
    void   AdjustWeights();

    ClassDef(AliGenPythia,1) // AliGenerator interface to Pythia
};
#endif





