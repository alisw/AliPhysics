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
class TRandom ; 

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
  virtual void   AddHit(Int_t, Int_t*, Float_t *) {
    // do not use this definition but the one below
    Fatal("AddHit(Int_t, Int_t*, Float_t *)", "do not use") ;
    
  }
  virtual void   AddHit( Int_t shunt, Int_t primary, Int_t track, 
			 Int_t id, Float_t *hits ) = 0 ;   
  virtual void Copy(AliPHOS & phos) ; 
  virtual AliDigitizer* CreateDigitizer(AliRunDigitizer* manager) const;
  virtual void  CreateMaterials() ;            
  virtual void  Digits2Raw();
  virtual void  FinishRun() {WriteQA();}
  virtual AliPHOSGeometry * GetGeometry() const 
  {return AliPHOSGeometry::GetInstance(GetTitle(),"") ;  }
  virtual void    Hits2SDigits();
  virtual Int_t   IsVersion(void) const = 0 ;  
  Int_t GetRawFormatHighGainFactor() const { return fHighGainFactor ; }  
  Int_t GetRawFormatHighGainOffset() const { return fHighGainOffset ; }  
  Int_t GetRawFormatTimeBins() const { return fkTimeBins ; }    
  Double_t GetRawFormatTimeMax() const { return fTimeMax ; }   
  Double_t GetRawFormatTimePeak() const { return fTimePeak ; }    
  Double_t GetRawFormatTimeRes() const { return fTimeRes ; }   
  virtual AliLoader* MakeLoader(const char* topfoldername);
  AliPHOSQAChecker * QAChecker() {return fQATask;}  
  static Double_t  RawResponseFunction(Double_t *x, Double_t *par) ; 
  virtual void    SetTreeAddress();   
  virtual TTree * TreeQA() const {return fTreeQA; } 
  virtual const TString Version() const {return TString(" ") ; } 
  virtual void WriteQA() ; 
  AliPHOS & operator = (const AliPHOS & /*rvalue*/)  {
    Fatal("operator =", "not implemented") ; return *this ; }


protected:

  Bool_t  RawSampledResponse(const Float_t dtime, const Int_t damp, Int_t * adcH, Int_t * adcL) const ; 


  AliPHOSQAChecker * fQATask ; //! PHOS checkers container
  TTree * fTreeQA ;            // the QA tree that contains the alarms
  Int_t    fHighGainFactor ;   // High gain attenuation factor of the raw RO signal
  Int_t    fHighGainOffset ;   // offset added to the module id to distinguish high and low gain data
  static const Int_t fkTimeBins = 256 ;     // number of sampling bins of the raw RO signal  
  Double_t fTimeMax ;          // maximum sampled time of the raw RO signal
  Double_t fTimePeak ;         // peaking time of the raw RO signal
  Double_t fTimeRes ;          // decay rime width of the raw RO signal 

  ClassDef(AliPHOS,3) // Photon Spectrometer Detector (base class)

} ;

#endif // ALIPHOS_H
