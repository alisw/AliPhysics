#ifndef ALIPHOSPID_H
#define ALIPHOSPID_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
                            
/* $Id$ */

//_________________________________________________________________________
//  Algorithm class for the identification of particles detected in PHOS        
//  base  class                             
//                  
//*-- Author: Yves Schutz (SUBATECH)

// --- ROOT system ---

#include "TObject.h" 
#include "TClonesArray.h"

// --- Standard library ---

// --- AliRoot header files ---

#include "AliPHOSTrackSegment.h"
#include "AliPHOSRecParticle.h"



class AliPHOSPID : public TObject {

public:

  AliPHOSPID() ;          // ctor            
  virtual ~AliPHOSPID() ; // dtor

  virtual void MakeParticles(AliPHOSTrackSegment::TrackSegmentsList * trsl, 
			     AliPHOSRecParticle::RecParticlesList * rpl) {} ; 
  virtual void SetShowerProfileCuts(Float_t, Float_t, Float_t, Float_t) {} ; 
  virtual void SetDispersionCutOff(Float_t ) {}    

  ClassDef(AliPHOSPID,1)  // Particle Identifier algorithm (base class)

} ;

#endif // ALIPHOSPID_H
