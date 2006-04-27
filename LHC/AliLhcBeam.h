#ifndef ALILHCBEAM_H
#define ALILHCBEAM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//
// Class that holds all parameters about an LHC beam.
// The parameters can change with time.
// A monitor can be set that stores the time distribution of 
// emittance and number of particles per bunch.
// Author: Andreas Morsch
// andreas.morsch@cern.ch
//

#include <TNamed.h>
#include "AliLhcMonitor.h"

class AliLHC;

class AliLhcBeam : public TNamed, public AliLhcMonitor
{
 public:
    AliLhcBeam(AliLHC* lhc);
    AliLhcBeam(const AliLhcBeam &beam);
    virtual ~AliLhcBeam();
    
    virtual void Init();
    
    virtual Float_t N() const {return fN;}
    virtual Float_t A() const {return fA;}
    virtual Float_t Z() const {return fZ;}
    virtual Float_t Emittance() const   {return fEmittance;}
    virtual Float_t Energy()    const   {return fEnergy;}
    virtual Float_t Gamma()     const   {return fGamma;}
    virtual Float_t LongEmittance() const {return fEmittanceL;}
    virtual Float_t EnergySpread()  const {return fEnergySpread;}

    virtual void  SetParticle(Float_t a, Float_t z) {fA = a; fZ = z;}
    virtual void  SetN(Float_t n) {fN = n;}	    
    virtual void  SetAccelerator(AliLHC* acc) {fAccelerator = acc;}
    virtual void  SetEnergy(Float_t  e) {fEnergy = e;}
    virtual void  SetNEmittance(Float_t  e) {fNEmittance = e;}
    virtual void  SetLongEmittance(Float_t  e) {fEmittanceL = e;}
    virtual void  SetEnergySpread(Float_t b) {fEnergySpread = b;}    
    
    virtual void  RemoveParticles(Float_t loss);
    virtual void  IncreaseEmittance(Float_t de, Float_t del);
    virtual void  SetMonitor(Int_t n);
    virtual void  Record();
    virtual void  DrawPlots();


    AliLhcBeam & operator=(const AliLhcBeam & rhs);
    
 protected:
    AliLHC* fAccelerator;         // Accelerator
    Float_t fN;                   // Number of Particles
    Float_t fN0;                  // Initial Number of Particles
    Float_t fNEmittance;          // Normalized Emittance
    Float_t fEmittance;           // Emittance
    Float_t fEmittance0;          // Initial Emittance
    Float_t fEmittanceL;          // Longitudinal Emittance
    Float_t fEmittanceL0;         // Longitudinal Emittance
    Float_t fEnergySpread;        // Energy Spread

    Float_t fA;                   // Atomic Number
    Float_t fZ;                   // Charge Number 
    Float_t fEnergy;              // Energy   
    Float_t fGamma;               // relativistic gamma 
    //
    Float_t* fTimeArray;          // [fNmax] Time array
    Float_t* fEmittanceArray;     // [fNmax] Emittance array
    Float_t* fEmittanceLArray;    // [fNmax] Long. Emittance array      
//
    ClassDef(AliLhcBeam,1) // LHC Beam
};

#endif





