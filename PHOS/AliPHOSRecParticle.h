#ifndef ALIPHOSRECPARTICLE_H
#define ALIPHOSRECPARTICLE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////
//  A Reconstructed Particle in PHOS          //
//  Yves Schutz SUBATECH                      //
//  To become a general class of AliRoot ?    //  
//                                            //
////////////////////////////////////////////////

// --- ROOT system ---

#include "TParticle.h"
#include "TVector3.h"

// --- Standard library ---

// --- AliRoot header files ---

#include "AliPHOSTrackSegment.h"
#include "AliPHOSFastRecParticle.h"

class AliPHOSRecParticle : public AliPHOSFastRecParticle {

public:
  
  AliPHOSRecParticle() {};          // ctor
  AliPHOSRecParticle(AliPHOSTrackSegment * ts) ;  // ctor
  AliPHOSRecParticle(const AliPHOSRecParticle & rp) ;  // ctor
  virtual ~AliPHOSRecParticle() ; // dtor

  AliPHOSTrackSegment * GetPHOSTrackSegment() const { return fPHOSTrackSegment ; } 

private:

  AliPHOSTrackSegment * fPHOSTrackSegment ; // pointer to the associated track segment in PHOS  

  ClassDef(AliPHOSRecParticle,1)  // Reconstructed Particle, version 1

};

#endif // AliPHOSRECPARTICLE_H
