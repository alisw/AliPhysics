#ifndef ALIPHOSPIDV1_H
#define ALIPHOSPIDV1_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////
//  Particle Identifier   class for PHOS         //
//  Version SUBATECH                             //
//  Author Yves Schutz                           //
//     comment: identify the type of particle    //  
//              PHOS SubTrack alone              //
///////////////////////////////////////////////////

// --- ROOT system ---

// --- Standard library ---

// --- AliRoot header files ---

#include "AliPHOSPID.h"

class  AliPHOSPIDv1 : public AliPHOSPID {

public:

  AliPHOSPIDv1() ;                     
  virtual ~ AliPHOSPIDv1() ; // dtor

  virtual void GetParticleType(TrackSegmentsList * trsl, RecParticlesList * rpl ) ; // does the job
  void Print() ; 
  virtual void SetDispersionCutOff(Float_t Dcut) {fCutOnDispersion = Dcut ; }    
  virtual void SetShowerProfileCuts(Float_t l1m, Float_t l1M, Float_t l2m, Float_t l2M) ; 

 private:

  // cuts on the shower profile 
  Float_t fLambda1m ; // minimum value for first elips axis
  Float_t fLambda1M ; // maximum value for first elips axis
  Float_t fLambda2m ; // minimum value for second elips axis
  Float_t fLambda2M ; // maximum value for second elips axis

  Float_t fCutOnDispersion ; // cut on the shower dispersion to distinguish hadronic from EM showers

  ClassDef( AliPHOSPIDv1,1)  // particle identifier implementation , version 1

};

#endif // AliPHOSPIDV1_H
