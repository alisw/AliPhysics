#ifndef ALIEMCAL_H
#define ALIEMCAL_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

/* $Id$ */

//_________________________________________________________________________
//  Base Class for EMCAL     
//                  
//*-- Author: Yves Schutz (SUBATECH)

#include <stdlib.h>

// --- ROOT system ---

class TString ;
class TTask ;
class TFolder ;

// --- AliRoot header files ---

#include "AliDetector.h"
//class AliDetector;
class AliEMCALGeometry ; 
//class AliEMCALQAChecker ;

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
  virtual void   AddHit(Int_t, Int_t*, Float_t *) {
    // do not use this definition but the one below
    Fatal("AddHit(Int_t, Int_t*, Float_t *", 
	  "not to be used: use AddHit( Int_t shunt, Int_t primary, Int_t track,Int_t id, Float_t *hits )") ;  

  }
  virtual void  AddHit( Int_t shunt, Int_t primary, Int_t track, 
			 Int_t id, Float_t *hits ) = 0 ;
  virtual void  CreateMaterials() ;   
  virtual void  FinishRun() {WriteQA();}                  
  virtual AliEMCALGeometry * GetGeometry() const ;   
  virtual Int_t   IsVersion(void) const = 0 ; 
  //AliEMCALQAChecker * QAChecker() const {return fQATask;}  
  virtual void  SetTreeAddress() ;
  virtual TTree * TreeQA() const {return fTreeQA; }                
  virtual const TString Version() const {return TString(" ") ; }  
  virtual void WriteQA() ; 
  AliEMCAL & operator = (const AliEMCAL & rvalue)  {
    // assignement operator requested by coding convention
    // but not needed
    Fatal("operator =", "not implemented") ;  
    return *this ; 
  }
 
  virtual AliLoader* MakeLoader(const char* topfoldername);
  
protected:

  //AliEMCALQAChecker * fQATask ; //! PHOS checkers container
  TTree * fTreeQA ;            // the QA tree that contains the alarms
  AliEMCALGeometry * fGeom ;   // the geometry object

  ClassDef(AliEMCAL,4) // Electromagnetic calorimeter (base class)

} ;

#endif // ALIEMCAL_H
