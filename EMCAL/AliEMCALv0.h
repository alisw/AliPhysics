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
//*-- and   : Sahal Yacoob (LBL / UCT) 
//#include <assert.h>

// --- ROOT system ---

class TFile;

// --- AliRoot header files ---
#include "AliEMCAL.h"

//class AliEMCALGeometry ; 

class AliEMCALv0 : public AliEMCAL {

 public:

  AliEMCALv0():AliEMCAL() {}
  AliEMCALv0(const char *name, const char *title="") ;
  AliEMCALv0(const AliEMCALv0 & emcal):AliEMCAL(emcal) {
    // cpy ctor: no implementation yet
    // requested by the Coding Convention
    Fatal("cpy ctor", "not implemented") ;  
  } 
  virtual ~AliEMCALv0(){} 

  virtual void  AddHit( Int_t shunt, Int_t primary, Int_t track, 
			Int_t id, Float_t *hits ) {
    // no hits - useless
  }

  virtual void BuildGeometry();// creates the geometry for the ROOT display
  virtual void CreateGeometry() ;// creates the geometry for GEANT
  virtual void   Init(void) ;                                       // does nothing
  virtual Int_t  IsVersion(void) const { 
    // Gives the version number 
    return 0 ; 
  }
  virtual const TString Version(void) const{ 
    // As above
    return TString("v0") ; 
  }
  
  AliEMCALv0 & operator = (const AliEMCALv0 & rvalue)  {
    // assignement operator requested by coding convention but not needed
    Fatal("operator =", "not implemented") ;  
    return *this ; 
  }
  
 protected:

  ClassDef(AliEMCALv0,2)  // Implementation of EMCAL manager class for midrapidity barrel layout between 0 and 120 degrees 
    
    };
    
#endif // AliEMCALV0_H
