#ifndef ALIEMCALRECPARTICLE_H
#define ALIEMCALRECPARTICLE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//_________________________________________________________________________
//  A Reconstructed Particle in EMCAL    
//  To become a general class of AliRoot ?        
//  why not      
//*-- Author: Yves Schutz (SUBATECH)

// --- ROOT system ---

// --- Standard library ---

// --- AliRoot header files ---

#include "AliEMCALFastRecParticle.h"
class TParticle ;

class AliEMCALRecParticle : public AliEMCALFastRecParticle {

 public:
  
  AliEMCALRecParticle() { fEMCALTrackSegment = 0 ; fDebug = kFALSE ; } 
  AliEMCALRecParticle(const AliEMCALRecParticle & rp) ;  // ctor
  virtual ~AliEMCALRecParticle(){  }

  Int_t   GetEMCALTSIndex()const {    return fEMCALTrackSegment ;  }
  virtual const Int_t GetNPrimariesToRecParticles() const ;
  virtual const Int_t GetNPrimaries() const ;
  virtual const TParticle * GetPrimary(Int_t index) const ;
  void    SetDebug() { fDebug = kTRUE ; } 
  void    UnsetDebug() { fDebug = kFALSE ; }
  void    SetTrackSegment(Int_t index){fEMCALTrackSegment = index; }

  typedef TClonesArray RecParticlesList ; 
  
 private:

  Int_t fEMCALTrackSegment ; // pointer to the associated track segment in EMCAL  
  Bool_t fDebug ; // to steer debug output 

  ClassDef(AliEMCALRecParticle,2)  // Reconstructed Particle
};

#endif // AliEMCALRECPARTICLE_H
