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
#include "AliPHOSGeometry.h" 
class AliPHOSQAChecker ;

class AliPHOS : public AliDetector {

 public:

  AliPHOS() ;
  AliPHOS(const char* name, const char* title="") ;  
  AliPHOS(AliPHOS & phos) : AliDetector(phos) {
    Copy(*this) ; 
  }
  virtual ~AliPHOS() ; 
  virtual void Copy(AliPHOS & phos) ; 
  virtual void   AddHit(Int_t, Int_t*, Float_t *) {
    // do not use this definition but the one below
    Fatal("AddHit(Int_t, Int_t*, Float_t *)", "do not use") ;
    
  }
  virtual void   AddHit( Int_t shunt, Int_t primary, Int_t track, 
			 Int_t id, Float_t *hits ) = 0 ;   
  virtual AliDigitizer* CreateDigitizer(AliRunDigitizer* manager) const;
  virtual void  CreateMaterials() ;            
  virtual void  FinishRun() {WriteQA();}
  virtual AliPHOSGeometry * GetGeometry() const 
  {return AliPHOSGeometry::GetInstance(GetTitle(),"") ;  }
  virtual void    Hits2SDigits();
  virtual Int_t   IsVersion(void) const = 0 ;  
  virtual AliLoader* MakeLoader(const char* topfoldername);
  AliPHOSQAChecker * QAChecker() {return fQATask;}  
  void    SetDebug() { fDebug = kTRUE ; }
  void    ResetDebug() { fDebug = kFALSE ; }
  Bool_t  Debug() const { return fDebug ; } 
  virtual void    SetTreeAddress();   
  virtual TTree * TreeQA() const {return fTreeQA; } 
  virtual const TString Version() const {return TString(" ") ; } 
  virtual void WriteQA() ; 
  AliPHOS & operator = (const AliPHOS & /*rvalue*/)  {
    Fatal("operator =", "not implemented") ; return *this ; }

protected:
  
  AliPHOSQAChecker * fQATask ; //! PHOS checkers container
  TTree * fTreeQA ;            // the QA tree that contains the alarms
  ClassDef(AliPHOS,2) // Photon Spectrometer Detector (base class)

} ;

#endif // ALIPHOS_H
