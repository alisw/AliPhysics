#ifndef ALIPHOSFASTRECPARTICLE_H
#define ALIPHOSFASTRECPARTICLE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//_________________________________________________________________________
//  A  Particle modified by PHOS response and produced by AliPHOSvFast
//  To become a general class of AliRoot ?    
//               
//*-- Author: Yves Schutz (SUBATECH)

// --- ROOT system ---

#include "TParticle.h"
#include "TVector3.h"

// --- Standard library ---

// --- AliRoot header files ---

const static Int_t kUNDEFINED =-1 ;
const static Int_t kNEUTRALEM = 0 ;
const static Int_t kNEUTRALHA = 1 ;
const static Int_t kGAMMA     = 2 ; 
const static Int_t kGAMMAHA   = 3 ; 
const static Int_t kABSURDEM  = 4 ;
const static Int_t kABSURDHA  = 5 ;
const static Int_t kELECTRON  = 6 ;
const static Int_t kCHARGEDHA = 7 ;


class AliPHOSFastRecParticle : public TParticle {

 public:
  
  AliPHOSFastRecParticle() {
    // ctor 
  };         
  AliPHOSFastRecParticle(const AliPHOSFastRecParticle & rp) ;  // ctor
  AliPHOSFastRecParticle(const TParticle & p) ;  // ctor
  virtual ~AliPHOSFastRecParticle(){
    // dtor
  }
  virtual Int_t DistancetoPrimitive(Int_t px, Int_t py) ; 
  virtual void Draw(Option_t *option) ;  
  virtual void ExecuteEvent(Int_t event, Int_t px, Int_t py) ;
  Int_t GetIndexInList() const { 
    // returns the index of this in the list
    return fIndexInList ; 
  } 
  virtual Int_t * GetPrimaries(Int_t & number) ;
  Int_t GetType() { 
    // returns the type of the particle
    return fType ; 
  } 
  TString Name() ; 
  virtual void Paint(Option_t * option="");
  void Print() ; 
  void SetPrimary(Int_t index) { 
    // sets the primary particle index
    fPrimary = index ; 
  }
  void SetType(Int_t type) { 
    // sets the particle type 
    fType = type ; 
  } 
  void SetIndexInList(Int_t val) { 
    // sets the value of the index in the list 
    fIndexInList = val ; 
  } 

 protected:

  Int_t fIndexInList ; // the index of this RecParticle in the list stored in TreeR (to be set by analysis)
  Int_t fPrimary ;     // (unique) primary particle index 
  Int_t fType ;        // particle type obtained by "virtual" reconstruction

 private:

  typedef TClonesArray  FastRecParticlesList ; 

  ClassDef(AliPHOSFastRecParticle,1)  // Reconstructed Particle produced by the fast simulation 

};

#endif // AliPHOSFASTRECPARTICLE_H
