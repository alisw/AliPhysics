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

#include "AliDetector.h"
class AliEMCALGeometry ; 

class AliEMCAL : public AliDetector {

 public:

  AliEMCAL() ;
  AliEMCAL(const char* name, const char* title="");
  AliEMCAL(const AliEMCAL & emcal) {
    // cpy ctor: no implementation yet
    // requested by the Coding Convention
    assert(0==1) ; 
  }
  virtual ~AliEMCAL() ; 
  virtual void   CreateMaterials() ;                     
  virtual AliEMCALGeometry * GetGeometry() { return fGeom ; }  
  Int_t   IsVersion(void) const { return -1 ; } 
  virtual void  SetTreeAddress() ;               
  TClonesArray *SDigits() const {return fSDigits;}
  virtual TString Version() {return TString(" ") ; }  
  AliEMCAL & operator = (const AliEMCAL & rvalue)  {
    // assignement operator requested by coding convention
    // but not needed
    assert(0==1) ;
    return *this ; 
  }
 
protected:

  AliEMCALGeometry * fGeom ;                       // Geometry definition
  TClonesArray   *fSDigits      ; // List of summable digits
  TClonesArray   *fDigits      ;  // List of digits

  ClassDef(AliEMCAL,2) // Electromagnetic calorimeter (base class)

} ;

#endif // ALIEMCAL_H
