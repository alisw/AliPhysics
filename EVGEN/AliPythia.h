#ifndef ALIPYTHIA_H
#define ALIPYTHIA_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include <TPythia.h>
#include "GenTypeDefs.h"

class AliPythia:public TPythia
{
 public:
    AliPythia();
    virtual ~AliPythia(){;}
    // convert to compressed code and print result (for debugging only)
    virtual Int_t CheckedLuComp(Int_t kf);
    // entry to the corresponding lujet function
    virtual void Lu1Ent(int flag, int idpart, 
			float mom, float theta,float phi);
    // Decay a Particle
    virtual void DecayParticle
	(Int_t idpart, Float_t mom, Float_t theta,Float_t phi);
    //
    // Pythia initialisation for selected processes
    virtual void ProcInit
	(Process_t process, Float_t energy, StrucFunc_t strucfunc);
    //
    // Count decay products
    virtual Int_t CountProducts(Int_t channel, Int_t particle);
    
    //
    // Force decay modes
    //
    // select type and multiplicity of the decay product
    virtual void ForceParticleDecay(Int_t particle,Int_t product,Int_t mult);
    //
    // force a particular decy-type
    virtual void ForceDecay(Decay_t  decay);
    // don't force any decay
    virtual void AllowAllDecays();
    
    //
    // Define heavy mesons to GEANT and make correspondance
    // GEANT - Pythia particle code
    virtual void DefineParticles();
    //
    // Convert from Pythia to Geant particle code
    virtual Int_t  GetGeantCode(Int_t kf);
    //
    // Get sum of branching ratios for forced decays
    virtual Float_t GetBraPart(Int_t kf);
 protected:
    Process_t     fProcess;         // Process type
    Decay_t       fDecay;           // Decay channel 
    Float_t       fEcms;            // Centre of mass energy
    StrucFunc_t   fStrucFunc;       // Structure function
    Int_t         fGPCode[501][2];  // Particle codes
    Float_t       fBraPart[501];    // Branching ratios
    static        Int_t fgInit;     // Make sure that only one object of type
                                    // AliPythia is ceated

    ClassDef(AliPythia,1) //ALICE UI to PYTHIA
};

#endif



