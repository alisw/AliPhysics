#ifndef ALIPHOSPARTICLEGUESSERV1_H
#define ALIPHOSPARTICLEGUESSERV1_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////
//  Particle Guesser    class for PHOS           //
//  Version SUBATECH                             //
//  Author Yves Schutz                           //
//     comment: guess the type of particle       //  
//              PHOS SubTrack alone              //
///////////////////////////////////////////////////

// --- ROOT system ---

// --- Standard library ---

// --- AliRoot header files ---

#include "AliPHOSParticleGuesser.h"

class  AliPHOSParticleGuesserv1 : public AliPHOSParticleGuesser {

public:

  AliPHOSParticleGuesserv1() ;                     
  virtual ~ AliPHOSParticleGuesserv1() ; // dtor

  void    GuessParticleType(TrackSegmentsList * trsl, RecParticlesList * rpl ) ; // does the job

  ClassDef( AliPHOSParticleGuesserv1,1)  // particle guesser implementation , version 1

};

#endif // AliPHOSPARTICLEGUESSERV1_H
