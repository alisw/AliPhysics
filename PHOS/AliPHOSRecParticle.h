#ifndef ALIPHOSRECPARTICLE_H
#define ALIPHOSRECPARTICLE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//_________________________________________________________________________
//  A Reconstructed Particle in PHOS    
//  To become a general class of AliRoot ?        
//  why not      
//*-- Author: Yves Schutz (SUBATECH)

// --- ROOT system ---

// --- Standard library ---

// --- AliRoot header files ---

#include "AliPHOSFastRecParticle.h"
class TParticle ;

class AliPHOSRecParticle : public AliPHOSFastRecParticle {

 public:
  
  AliPHOSRecParticle() { fPHOSTrackSegment = 0 ; fDebug = kFALSE ; } 
  AliPHOSRecParticle(const AliPHOSRecParticle & rp) ;  // ctor
  virtual ~AliPHOSRecParticle(){  }

  Int_t   GetPHOSTSIndex()const {    return fPHOSTrackSegment ;  }
  virtual const Int_t GetNPrimariesToRecParticles() const ;
  virtual const Int_t GetNPrimaries() const ;
  virtual const TParticle * GetPrimary(Int_t index) const ;
  void    SetDebug() { fDebug = kTRUE ; } 
  void    UnsetDebug() { fDebug = kFALSE ; }
  void    SetTrackSegment(Int_t index){fPHOSTrackSegment = index; }

  typedef TClonesArray RecParticlesList ; 
  
 private:

  Int_t fPHOSTrackSegment ; // pointer to the associated track segment in PHOS  
  Bool_t fDebug ; // to steer debug output 

  ClassDef(AliPHOSRecParticle,2)  // Reconstructed Particle
};

#endif // AliPHOSRECPARTICLE_H
