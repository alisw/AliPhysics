#ifndef ALIPHOSPARTICLEGUESSER_H
#define ALIPHOSPARTICLEGUESSER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
                            
/* $Id$ */

////////////////////////////////////////////////
//  Algorithme class for the guess of         //
//          particles detected in PHOS        //
//  interface class                           //
//  Version SUBATECH                          //
//  Author Yves Schutz     SUBATECH           //
//                                            //  
//   pABC                                     //
////////////////////////////////////////////////

// --- ROOT system ---

#include "TObject.h" 
#include "TClonesArray.h"

// --- Standard library ---

// --- AliRoot header files ---

#include "AliPHOSTrackSegmentMaker.h"


typedef TClonesArray RecParticlesList ; 

class AliPHOSParticleGuesser : public TObject {

public:

  AliPHOSParticleGuesser() ;          // ctor            
  virtual ~AliPHOSParticleGuesser() ; // dtor

  virtual void GuessParticleType(TrackSegmentsList * trsl, RecParticlesList * rpl) {} ; 

  ClassDef(AliPHOSParticleGuesser,1)  // Particle Guesser interface, version 1

} ;

#endif // ALIPHOSPARTICLEGUESSER_H
