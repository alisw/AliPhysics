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
#include <cassert>

const static Int_t kUNDEFINED     = -1; 
const static Int_t kGAMMA         = 0 ; 
const static Int_t kELECTRON      = 1 ;
const static Int_t kNEUTRAL       = 2 ;  
const static Int_t kCHARGED       = 3 ;  
const static Int_t kCHARGEDHADRON = 4 ;  
const static Int_t kNEUTRALHADRON = 5 ;  
const static Int_t kNEUTRALEM     = 6 ;  
const static Int_t kGAMMAHADRON   = 7 ; 

class AliPHOSRecParticle : public TParticle {

public:
  
  AliPHOSRecParticle() {};          // ctor
  AliPHOSRecParticle(AliPHOSTrackSegment * ts) ;  // ctor
  AliPHOSRecParticle(const AliPHOSRecParticle & rp) ;  // ctor
  AliPHOSRecParticle& operator= (const AliPHOSRecParticle& rp) { assert(0==1) ; } // forbidden

  virtual ~AliPHOSRecParticle() ; // dtor

  virtual Int_t DistancetoPrimitive(Int_t px, Int_t py) ; 
  virtual void Draw(Option_t *option) ;  
  virtual void ExecuteEvent(Int_t event, Int_t px, Int_t py) ; 
  AliPHOSTrackSegment * GetPHOSTrackSegment() const { return fPHOSTrackSegment ; } 
  Int_t GetType() { return fType ; } 
  TString Name() ; 
  virtual void Paint(Option_t * option="");
  void Print() ; 
  void SetType(Int_t type) { fType = type ; } 

private:

  AliPHOSTrackSegment * fPHOSTrackSegment ; // pointer to the associated track segment in PHOS  
  Int_t fType ;                             // identified particle type

  ClassDef(AliPHOSRecParticle,1)  // Reconstructed Particle, version 1

};

#endif // AliPHOSRECPARTICLE_H
