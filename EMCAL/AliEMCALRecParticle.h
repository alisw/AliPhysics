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

#include "AliESDtrack.h" 
#include "AliEMCALFastRecParticle.h"
class TParticle ;
#include  "TVector3.h"  

class AliEMCALRecParticle : public AliEMCALFastRecParticle {

 public:
  
  AliEMCALRecParticle() ;  
  AliEMCALRecParticle(const AliEMCALRecParticle & rp) ;  // ctor
  virtual ~AliEMCALRecParticle(){  }

  Int_t   GetEMCALRPIndex()const {    return fEMCALRecPoint ;  }
  virtual const Int_t GetNPrimariesToRecParticles() const ;
  virtual const Int_t GetNPrimaries() const ;
  TVector3 GetPos() const { return fPos ; } 
  virtual const TParticle * GetPrimary(Int_t index) const ;
  const Double_t *GetPID();
  void    SetDebug() { fDebug = kTRUE ; } 
  void    SetPos(TVector3 pos) { fPos.SetXYZ( pos.X(), pos.Y(), pos.Z() ); } 
  void    UnsetDebug() { fDebug = kFALSE ; }
  void    SetRecPoint(Int_t index){fEMCALRecPoint = index; }

  typedef TClonesArray RecParticlesList ; 
  
 private:

  Int_t fEMCALRecPoint ; // pointer to the associated track segment in EMCAL  
  Bool_t fDebug ; // to steer debug output 
  TVector3 fPos ; // position in the global alice coordinate system 
  Double_t fPID[AliESDtrack::kSPECIESN] ; // PID probability densities

  ClassDef(AliEMCALRecParticle,3)  // Reconstructed Particle
};

#endif // AliEMCALRECPARTICLE_H
