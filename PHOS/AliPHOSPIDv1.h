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

  void    GetParticleType(TrackSegmentsList * trsl, RecParticlesList * rpl ) ; // does the job

  ClassDef( AliPHOSPIDv1,1)  // particle identifier implementation , version 1

};

#endif // AliPHOSPIDV1_H
