#ifndef ALIPHOS_H
#define ALIPHOS_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
/* $Id$ */

//_________________________________________________________________________
//  Base Class for PHOS     
//                  
//*-- Author: Laurent Aphecetche & Yves Schutz (SUBATECH)

// --- ROOT system ---

class TString ; 

// --- AliRoot header files ---

#include "AliDetector.h" 
class AliPHOSGeometry ; 

class AliPHOS : public AliDetector {

 public:

  AliPHOS() {}
  AliPHOS(const char* name, const char* title="") {}
  AliPHOS(const AliPHOS & phos) {
    // cpy ctor: no implementation yet
    // requested by the Coding Convention
    abort() ; 
  }
  virtual ~AliPHOS() {}
  virtual void   AddHit(Int_t, Int_t*, Float_t *) {
    // do not use this definition but the one below
    abort() ; 
  }
  virtual void   AddHit( Int_t shunt, Int_t primary, Int_t track, Int_t id, Float_t *hits ) = 0 ;   
  virtual void   CreateMaterials() ;                     
  virtual  AliPHOSGeometry * GetGeometry() = 0 ;

  Int_t   IsVersion(void) const { return -1 ; } 
  virtual void    SetTreeAddress();                
  virtual TString Version() {return TString(" ") ; } 
 
  AliPHOS & operator = (const AliPHOS & rvalue)  {
    // assignement operator requested by coding convention
    // but not needed
    abort() ;
    return *this ; 
  }
 
protected:

  ClassDef(AliPHOS,2) // Photon Spectrometer Detector (base class)

} ;

#endif // ALIPHOS_H
