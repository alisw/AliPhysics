#ifndef ALIEMCAL_H
#define ALIEMCAL_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

/* $Id$ */

//_________________________________________________________________________
//  Base Class for EMCAL     
//                  
//*-- Author: Yves Schutz (SUBATECH)

// --- ROOT system ---

class TString ;
class TTask ;
class TFolder ;

// --- AliRoot header files ---

#include "AliDetector.h"
class AliEMCALGeometry ; 
class AliEMCAL : public AliDetector {

 public:

  AliEMCAL(); 
  AliEMCAL(const char* name, const char* title="");
  AliEMCAL(const AliEMCAL& emcal) : AliDetector(emcal) {
    // cpy ctor: no implementation yet
    // requested by the Coding Convention
    Fatal("cpy ctor", "not implemented") ;  
  }
  virtual ~AliEMCAL() ; 
  virtual void   AddHit(Int_t, Int_t*, Float_t *) const{
    Fatal("AddHit(Int_t, Int_t*, Float_t *", "not to be used: use AddHit( Int_t shunt, Int_t primary, Int_t track,Int_t id, Float_t *hits )") ;  
  }
  virtual void  CreateMaterials() ;   
  virtual void  FinishRun() {}                  
  virtual AliEMCALGeometry * GetGeometry() const ;   
  virtual Int_t   IsVersion(void) const = 0 ;   
  virtual void  SetTreeAddress() ;              
  virtual const TString Version() const {return TString(" ") ; }   
  AliEMCAL & operator = (const AliEMCAL & /*rvalue*/)  {
    Fatal("operator =", "not implemented") ;  return *this ; }
 
  virtual AliLoader* MakeLoader(const char* topfoldername);
  
  virtual void    Hits2SDigits();
  virtual AliDigitizer* CreateDigitizer(AliRunDigitizer* manager) const;

protected:
  AliEMCALGeometry * fGeom ;   // the geometry object

  ClassDef(AliEMCAL,5) // Electromagnetic calorimeter (base class)

} ;

#endif // ALIEMCAL_H
