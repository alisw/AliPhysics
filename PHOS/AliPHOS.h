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
class TTask ;
class TFolder ;
class TTree ; 

// --- AliRoot header files ---
#include "AliDetector.h" 
class AliPHOSGeometry ; 
class AliPHOSQAChecker ;

class AliPHOS : public AliDetector {

 public:

  AliPHOS() ;
  AliPHOS(const char* name, const char* title="") ;  
  AliPHOS(const AliPHOS & phos) : AliDetector(phos) {
    // cpy ctor: no implementation yet
    // requested by the Coding Convention
    Fatal("cpy ctor", "not implemented") ;
  }
  virtual ~AliPHOS() ; 
  virtual void   AddHit(Int_t, Int_t*, Float_t *) {
    // do not use this definition but the one below
    Fatal("AddHit(Int_t, Int_t*, Float_t *)", "do not use") ;
    
  }
  virtual void   AddHit( Int_t shunt, Int_t primary, Int_t track, 
			 Int_t id, Float_t *hits ) = 0 ;   
  virtual void   CreateMaterials() ;                     
  virtual void  FinishRun() {WriteQA();}
  virtual AliPHOSGeometry * GetGeometry() const ;
  virtual Int_t   IsVersion(void) const = 0 ;  
  AliPHOSQAChecker * QAChecker() {return fQATask;}  
  virtual void    SetTreeAddress();   
  virtual TTree * TreeQA() const {return fTreeQA; } 
  virtual const TString Version() const {return TString(" ") ; } 
  virtual void WriteQA() ; 
  AliPHOS & operator = (const AliPHOS & /*rvalue*/)  {
    Fatal("operator =", "not implemented") ; return *this ; }

  virtual AliLoader* MakeLoader(const char* topfoldername);
 
  virtual void    Hits2SDigits();
  virtual AliDigitizer* CreateDigitizer(AliRunDigitizer* manager);

protected:
  
  AliPHOSQAChecker * fQATask ; //! PHOS checkers container
  TTree * fTreeQA ;            // the QA tree that contains the alarms

  ClassDef(AliPHOS,2) // Photon Spectrometer Detector (base class)

} ;

#endif // ALIPHOS_H
