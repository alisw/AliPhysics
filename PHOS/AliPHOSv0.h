#ifndef ALIPHOSV0_H
#define ALIPHOSV0_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//_________________________________________________________________________
// Implementation version v0 of PHOS Manager class 
// Layout EMC + PPSD has name GPS2  
//                  
//*-- Author: Yves Schutz (SUBATECH)

// --- ROOT system ---
#include "TClonesArray.h"

// --- AliRoot header files ---
#include "AliPHOS.h"
#include "AliPHOSGeometry.h"
#include "AliPHOSReconstructioner.h"
#include "AliPHOSTrackSegmentMaker.h"
#include "AliPHOSPID.h"

class AliPHOSv0 : public AliPHOS {

public:

  AliPHOSv0(void){
    // ctor
  } 
  AliPHOSv0(const char *name, const char *title="") ;
  AliPHOSv0(AliPHOSReconstructioner * Reconstructioner, const char *name, const char *title="") ;
  AliPHOSv0(const AliPHOSv0 & phos) {
    // cpy ctor: no implementation yet
    // requested by the Coding Convention
    assert(0==1) ; 
  } 
  virtual ~AliPHOSv0(void){
    // dtor
  } 

  virtual void   BuildGeometry(void) ;                              // creates the geometry for the ROOT display
  void           BuildGeometryforPHOS(void) ;                       // creates the PHOS geometry for the ROOT display
  void           BuildGeometryforPPSD(void) ;                       // creates the PPSD geometry for the ROOT display
  virtual void   CreateGeometry(void) ;                             // creates the geometry for GEANT
  void           CreateGeometryforPHOS(void) ;                      // creates the PHOS geometry for GEANT
  void           CreateGeometryforPPSD(void) ;                      // creates the PPSD geometry for GEANT
  virtual AliPHOSGeometry * GetGeometry() { return fGeom ; }  
  virtual void   Init(void) ;                                       // does nothing
  Int_t IsVersion(void) const { return 0 ; }
  virtual TString Version(void){ return TString("v0"); }
  
  AliPHOSv0 & operator = (const AliPHOSv0 & rvalue)  {
    // assignement operator requested by coding convention
    // but not needed
    assert(0==1) ;
    return *this ; 
  }
  
 protected:
  
  AliPHOSGeometry * fGeom ;                       // Geometry definition
  
  ClassDef(AliPHOSv0,1)  // Implementation of PHOS manager class for layout EMC+PPSD
    
    };
    
#endif // AliPHOSV0_H
