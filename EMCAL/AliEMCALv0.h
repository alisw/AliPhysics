#ifndef ALIEMCALV0_H
#define ALIEMCALV0_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                      
         */
/* $Id$ */

//_________________________________________________________________________
// Implementation version v0 of EMCAL Manager class 
//*--                  
//*-- Author: Yves Schutz (SUBATECH)

#include <assert.h>

// --- ROOT system ---

// --- AliRoot header files ---
#include "AliEMCAL.h"
class AliEMCALGeometry ; 

class AliEMCALv0 : public AliEMCAL {

 public:

  AliEMCALv0() {fGeom=0;}
  AliEMCALv0(const char *name, const char *title="") ;
  AliEMCALv0(const AliEMCALv0 & emcal) {
    // cpy ctor: no implementation yet
    // requested by the Coding Convention
    assert(0==1) ; 
  } 
  virtual ~AliEMCALv0(void){} 

  virtual void   BuildGeometry(void) ;             // creates the geometry for the ROOT display
  virtual void   CreateGeometry(void) ;            // creates the geometry for GEANT

  virtual void   Init(void) ;                                       // does nothing
  virtual Int_t  IsVersion(void) const { 
    // Gives the version number 
    return 0 ; 
  }
  virtual TString Version(void){ 
    // As above
    return TString("v0") ; 
  }
  
  AliEMCALv0 & operator = (const AliEMCALv0 & rvalue)  {
    // assignement operator requested by coding convention but not needed
    assert(0==1) ;
    return *this ; 
  }
  
 protected:
    
  ClassDef(AliEMCALv0,1)  // Implementation of EMCAL manager class for layout EMC+PPSD
    
    };
    
#endif // AliEMCALV0_H
