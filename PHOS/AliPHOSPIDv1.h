#ifndef ALIPHOSPIDV1_H
#define ALIPHOSPIDV1_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */


//_________________________________________________________________________
// Implementation version v1 of the PHOS particle identifier 
// Identification is based on information from PPSD and EMC
//                  
//*-- Author: Yves Schutz (SUBATECH)

// --- ROOT system ---

// --- Standard library ---

// --- AliRoot header files ---

#include "AliPHOSPID.h"

class  AliPHOSPIDv1 : public AliPHOSPID {

public:

  AliPHOSPIDv1() 
  { 
    fCutOnDispersion = 1.5; 
    fCutOnRelativeDistance = 3.0 ;
  }
                     
  virtual ~ AliPHOSPIDv1(){} ; // dtor


  Float_t GetDistanceInPHOSPlane(AliPHOSEmcRecPoint * emcclu, AliPHOSPpsdRecPoint * PpsdClu, Bool_t &toofar, Option_t * Axis) ; // Relative Distance PPSD-EMC
  virtual void MakeParticles(AliPHOSTrackSegment::TrackSegmentsList * trsl, 
			     AliPHOSRecParticle::RecParticlesList * rpl ) ; // does the job
  virtual void Print(const char *) ; 
  virtual void SetDispersionCutOff(Float_t Dcut) {fCutOnDispersion = Dcut ; }    
  virtual void SetShowerProfileCuts(Float_t l1m, Float_t l1M, Float_t l2m, Float_t l2M) ; 
  virtual void SetRelativeDistanceCut(Float_t CutOnRelativeDistance) ;
 

 private:

  // cuts on the shower profile 
  Float_t fLambda1m ;        // minimum value for first elips axis
  Float_t fLambda1M ;        // maximum value for first elips axis
  Float_t fLambda2m ;        // minimum value for second elips axis
  Float_t fLambda2M ;        // maximum value for second elips axis
  Float_t fCutOnDispersion ; // cut on the shower dispersion to distinguish hadronic from EM showers
  Float_t fCutOnRelativeDistance; //Cut on the relative distance between PPSD and EMC

  ClassDef( AliPHOSPIDv1,1)  // Particle identifier implementation version 1

};

#endif // AliPHOSPIDV1_H
