#ifndef ALIPHOSPID_H
#define ALIPHOSPID_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
                            
/* $Id$ */

////////////////////////////////////////////////
//  Algorithme class for the identification of//
//          particles detected in PHOS        //
//  base  class                               //
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

class AliPHOSPID : public TObject {

public:

  AliPHOSPID() ;          // ctor            
  virtual ~AliPHOSPID() ; // dtor

  virtual void GetParticleType(TrackSegmentsList * trsl, RecParticlesList * rpl) {} ; 

  ClassDef(AliPHOSPID,1)  // Particle Identifier interface, version 1

} ;

#endif // ALIPHOSPID_H
