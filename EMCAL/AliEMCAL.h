#ifndef ALIEMCAL_H
#define ALIEMCAL_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

/* $Id$ */

//_________________________________________________________________________
//  Base Class for EMCAL     
//                  
//*-- Author: Yves Schutz (SUBATECH)

#include <assert.h>

// --- ROOT system ---
class TString ;
class TClonesArray ;

// --- AliRoot header files ---
#include <stdlib.h>
#include "AliDetector.h"
#include "AliEMCALGeometry.h"
class AliEMCALGeometry ; 

class AliEMCAL : public AliDetector {

 public:

  AliEMCAL(); 
  AliEMCAL(const char* name, const char* title="");
  AliEMCAL(const AliEMCAL & emcal) {
    // cpy ctor: no implementation yet
    // requested by the Coding Convention
    abort() ; 
  }
  virtual ~AliEMCAL() ; 
  virtual void   AddHit(Int_t, Int_t*, Float_t *) {
    // do not use this definition but the one below
    abort() ;
  }
  virtual void   AddHit( Int_t shunt, Int_t primary, Int_t track, 
			 Int_t id, Float_t *hits ) = 0 ;


  virtual void   CreateMaterials() ;                     
  virtual AliEMCALGeometry * GetGeometry()  = 0 ;   
  Int_t   IsVersion(void) const { return -1 ; } 
  virtual void  SetTreeAddress() ;               
  virtual TString Version() {return TString(" ") ; }  
  AliEMCAL & operator = (const AliEMCAL & rvalue)  {
    // assignement operator requested by coding convention
    // but not needed
    abort() ;
    return *this ; 
  }
 
protected:

  AliEMCALGeometry * fGeom ;                       // Geometry definition

  ClassDef(AliEMCAL,1) // Electromagnetic calorimeter (base class)

} ;

#endif // ALIEMCAL_H
