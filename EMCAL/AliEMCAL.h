#ifndef ALIEMCAL_H
#define ALIEMCAL_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

/* $Id$ */
/* History of cvs commits:
 *
 * $Log$
 *
 */
//_________________________________________________________________________
//  Base Class for EMCAL     
//  holds all geant information of
//  materials, etc.
//                  
//*-- Author: Yves Schutz (SUBATECH)

// --- ROOT system ---

class TString ;
class TTask ;
class TFolder ;
class TRandom ; 
class TGraph;
class TF1;

// --- AliRoot header files ---
class AliRawReader;
#include "AliDetector.h"
#include "AliEMCALGeometry.h" 
#include "AliEMCALTrigger.h" 

class AliEMCAL : public AliDetector {

 public:
  
  AliEMCAL(); 
  AliEMCAL(const char* name, const char* title="");

  virtual ~AliEMCAL() ; 
  virtual void   AddHit(Int_t, Int_t*, Float_t *) {
    Fatal("AddHit(Int_t, Int_t*, Float_t *", "not to be used: use AddHit( Int_t shunt, Int_t primary, Int_t track,Int_t id, Float_t *hits )") ;  
  }
  virtual AliDigitizer* CreateDigitizer(AliRunDigitizer* manager) const;
  virtual void  CreateMaterials() ;   
  //  virtual void  
  virtual void  Digits2Raw();
  
  using AliDetector::Raw2Digits;
  virtual void  Raw2Digits(AliRawReader *reader);
  
  virtual void  FinishRun() {}                  
  virtual AliEMCALGeometry * GetGeometry() const 
    {return AliEMCALGeometry::GetInstance(GetTitle(),"") ;  }   
  virtual void    Hits2SDigits();
  virtual Int_t   IsVersion(void) const = 0 ;   
  
  virtual AliTriggerDetector* CreateTriggerDetector() const 
    { return new AliEMCALTrigger(); }

  // Raw Read Out
  Double_t GetRawFormatCapa() const { return fgCapa ; }   
  Double_t GetRawFormatHighCharge() const { return fHighCharge ; }  
  Double_t GetRawFormatHighGain() const { return fHighGain ; }  
  Double_t GetRawFormatHighLowGainFactor() const { return fHighLowGainFactor ; }  
  Double_t GetRawFormatLowCharge() const { return ( fHighCharge *  fHighLowGainFactor ) ; }  
  Double_t GetRawFormatLowGain() const { return ( fHighGain / fHighLowGainFactor ) ; }  
  Int_t GetRawFormatLowGainOffset() const { return fLowGainOffset ; }  
  Int_t GetRawFormatOrder() const { return fgOrder ; }   
  Int_t GetRawFormatTimeBins() const { return fgkTimeBins ; }    
  Double_t GetRawFormatTimeMax() const { return fgTimeMax ; }   
  Double_t GetRawFormatTimePeak() const { return fgTimePeak ; }    
  Double_t GetRawFormatTimeTrigger() const { return fgTimeTrigger ; }
  Int_t GetRawFormatThreshold() const { return fgThreshold ; }       
  Int_t GetRawFormatDDLPerSuperModule() const { return fgDDLPerSuperModule ; }       
  static Double_t RawResponseFunctionMax(Double_t charge, Double_t gain) ;
  Bool_t   RawSampledResponse(Double_t dtime, Double_t damp, Int_t * adcH, Int_t * adcL) const ; 
  //  
  virtual AliLoader* MakeLoader(const char* topfoldername);
  virtual const TString Version() const {return TString(" ") ; }   

protected:
  
  static Double_t RawResponseFunction(Double_t *x, Double_t *par) ; 
  void FitRaw(Bool_t lowGainFlag, TGraph * gLowGain, TGraph * gHighGain, TF1* signalF, Double_t & energy, Double_t & time) ;

  void Init(void);  //initializes some params

  Int_t fBirkC0;    // constants for Birk's Law implementation
  Double_t fBirkC1; // constants for Birk's Law implementation
  Double_t fBirkC2; // constants for Birk's Law implementation

  static Double_t fgCapa ;              // capacitor of the preamplifier for the raw RO signal
  Double_t fHighCharge ;                // high charge (to convert energy to charge) for the raw RO signal
  Double_t fHighGain ;                  // high gain for the raw RO signal
  Double_t fHighLowGainFactor ;         // high to low gain factor for the raw RO signal
  Int_t    fLowGainOffset ;             // to separate high from low gain in the DDL
  static Int_t fgOrder ;                // order of the gamma function for the RO signal
  static const Int_t fgkTimeBins = 256 ; // number of sampling bins of the raw RO signal  
  static Double_t fgTimeMax ;           // maximum sampled time of the raw RO signal                             
  static Double_t fgTimePeak ;          // peaking time of the raw RO signal                                    
  static Double_t fgTimeTrigger ;       // time of the trigger for the RO signal 
  static Int_t fgThreshold;             // threshold
  static Int_t fgDDLPerSuperModule;        // number of DDL per SuperModule
 
private:
  AliEMCAL(const AliEMCAL& emcal);
  AliEMCAL & operator = (const AliEMCAL & /*rvalue*/);

  ClassDef(AliEMCAL,9) // Electromagnetic calorimeter (base class)
    
    } ;

#endif // ALIEMCAL_H
