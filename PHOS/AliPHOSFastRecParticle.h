#ifndef ALIPHOSFASTRECPARTICLE_H
#define ALIPHOSFASTRECPARTICLE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////
//  A  Particle modified by PHOS response     //
//  Yves Schutz SUBATECH                      //
//  To become a general class of AliRoot ?    //  
//                                            //
////////////////////////////////////////////////

// --- ROOT system ---

#include "TParticle.h"
#include "TVector3.h"

// --- Standard library ---

// --- AliRoot header files ---

#include "AliPHOSRecParticle.h" 

typedef TClonesArray  FastRecParticlesList ; 


class AliPHOSFastRecParticle : public TParticle {

public:
  
  AliPHOSFastRecParticle() {};          // ctor
  AliPHOSFastRecParticle(const AliPHOSFastRecParticle & rp) ;  // ctor
  AliPHOSFastRecParticle(const TParticle & p) ;  // ctor
  virtual ~AliPHOSFastRecParticle() ; // dtor

  virtual Int_t DistancetoPrimitive(Int_t px, Int_t py) ; 
  virtual void Draw(Option_t *option) ;  
  virtual void ExecuteEvent(Int_t event, Int_t px, Int_t py) ; 
  Int_t GetType() { return fType ; } 
  TString Name() ; 
  virtual void Paint(Option_t * option="");
  void Print() ; 
  void SetType(Int_t type) { fType = type ; } 

private:

  Int_t fType ;   // identified particle type

  ClassDef(AliPHOSFastRecParticle,1)  // Reconstructed Particle, version 1

};

#endif // AliPHOSFASTRECPARTICLE_H
