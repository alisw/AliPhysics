#ifndef ALIPHOSINDEXTOOBJECT_H
#define ALIPHOSINDEXTOOBJECT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//_________________________________________________________________________
//  A singleton that retrieves objets from an array stored in a Tree on a disk file
//    1. AliPHOSDigit from TreeD     
//                  
//*-- Author: Yves Schutz (SUBATECH)


// --- ROOT system ---

#include "TFile.h"
#include "TString.h"
#include "TParticle.h"

// --- Standard library ---

#include "assert.h"

// --- AliRoot header files ---

#include "AliPHOS.h" 
#include "AliRun.h" 
#include "AliPHOSDigit.h" 
#include "AliPHOSEmcRecPoint.h"
#include "AliPHOSPpsdRecPoint.h"
#include "AliPHOSTrackSegment.h"
#include "AliPHOSRecParticle.h"

class AliPHOSIndexToObject : public TObject {

public:

  AliPHOSIndexToObject(){ 
    // ctor: this is a singleton, the ctor should never be called but cint needs it as publiv
    assert(0==1) ; 
  } 
  virtual ~AliPHOSIndexToObject(){
    // dtor
  }

  static AliPHOSIndexToObject * GetInstance(AliPHOS * det) ; 
  static AliPHOSIndexToObject * GetInstance() ; 
  
  AliPHOSDigit *        GimeDigit(Int_t index) ; 
  TParticle *           GimePrimaryParticle(Int_t index) ;
  AliPHOSRecParticle *  GimeRecParticle(Int_t index) ; 
  AliRecPoint *         GimeRecPoint(Int_t index, TString s) ; 
  AliPHOSTrackSegment * GimeTrackSegment(Int_t index) ;
  
 private:
  
  AliPHOSIndexToObject(AliPHOS * det) ; 

  AliPHOS * fDetector ;                    // the detector 
  TTree * fReconstruct ;                   // the reconstruction tree  

  static AliPHOSIndexToObject * fgObjGetter ; // pointer to the unique instance of the singleton 

  ClassDef(AliPHOSIndexToObject,1)  // Algorithm class that provides methods to retrieve objects from a list knowing the index 

};

#endif // AliPHOSINDEXTOOBJECT_H
