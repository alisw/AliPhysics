#ifndef ALIEMCALDIGIT_H
#define ALIEMCALDIGIT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//_________________________________________________________________________
//  EMCAL digit: 
// 
//  A  Digit is the sum of energy in a Tower (Hit sum) and stores information, about primaries
//  and entering particle contributing to a Digit
//
//*-- Author: Sahal Yacoob (LBL)
// based on : AliPHOSDigit
//___________________________________________________________________________

// --- ROOT system ---

#include "TObject.h" 

// --- Standard library ---

// --- AliRoot header files ---
#include "AliDigitNew.h"

class AliEMCALDigit : public AliDigitNew {

  friend ostream& operator << ( ostream& , const AliEMCALDigit&) ;

 public:
  
  AliEMCALDigit() ;
  AliEMCALDigit(Int_t primary, Int_t iparent, Int_t id, Float_t digEnergy, Float_t time, Int_t type,Int_t index = -1, Float_t chi2=0, Int_t ndf=0, Float_t dE = 0) ;
  AliEMCALDigit(const AliEMCALDigit & digit) ;
  virtual ~AliEMCALDigit() ;

  Bool_t operator==(const AliEMCALDigit &rValue) const;
  AliEMCALDigit operator+(const AliEMCALDigit &rValue) ;
  AliEMCALDigit operator*(Float_t factor) ; 
  AliEMCALDigit& operator = (const AliEMCALDigit & digit) ;
  
  enum  digitType{kUnknown=-1, kHG=0, kLG=1, kLGnoHG=2, kTrigger=3, kEmbedded = 4};

  void     Clear(const Option_t*) ;	
  Int_t    Compare(const TObject * obj) const ;
  Float_t  GetAmplitude()   const { if(!fAmp)return fAmpFloat ; else return fAmp ;}//Keep backward compatibility.
  Float_t  GetEta()         const ; 
  Int_t    GetNprimary()    const { return fNprimary ;}
  Int_t    GetPrimary(Int_t index)   const ; 
  Float_t  GetDEPrimary(Int_t index) const ; 
  Int_t    GetNiparent()    const { return fNiparent ;}
  Int_t    GetIparent(Int_t index)   const ;
  Float_t  GetDEParent(Int_t index)  const ; 
  Float_t  GetPhi()         const ;
  Float_t  GetTime(void)    const { return fTime  ;}
  Float_t  GetTimeR(void)   const { return fTimeR ;}
  Float_t  GetChi2(void)    const { return fChi2  ;}
  Int_t    GetNDF(void)     const { return fNDF   ;}
  Bool_t   IsSortable()     const { return kTRUE  ;}
  Int_t    GetType()        const { return fDigitType ;}
	
  void     SetAmp(Int_t amp)         { fAmp       = amp  ; } //old
  void     SetAmplitude(Float_t amp) { fAmpFloat  = amp  ; }
  void     SetId(Int_t idt)          { fId        = idt  ; }
  void     SetTime(Float_t time)     { fTime      = time ; }
  void     SetTimeR(Float_t time)    { fTimeR     = time ; }
  void     SetChi2(Float_t chi)      { fChi2      = chi  ; }
  void     SetNDF(Int_t ndf)         { fNDF       = ndf  ; }
  void     SetType(Int_t t)          { fDigitType = t    ; }
  void     ShiftPrimary(Int_t shift); // shift to separate different TreeK in merging

  //Raw time sample
  //ALTRO
  Int_t    GetNALTROSamplesLG() const {if(fDigitType==kLG)      return fNSamples;   else return 0 ; }
  Bool_t   GetALTROSampleLG(const Int_t iSample, Int_t& timeBin, Int_t& amp) const;
  Int_t    GetNALTROSamplesHG() const {if(fDigitType==kHG)      return fNSamplesHG; else return 0 ; }
  Bool_t   GetALTROSampleHG(const Int_t iSample, Int_t& timeBin, Int_t& amp) const;
  //FALTRO, trigger. Same data members as Low Gain	
  Int_t    GetNFALTROSamples()  const {if(fDigitType==kTrigger) return fNSamples;   else return 0 ; }
  Bool_t   GetFALTROSample(const Int_t iSample, Int_t& timeBin, Int_t& amp)  const ;
	
  void     SetALTROSamplesHG (const Int_t nSamplesHG, Int_t *samplesHG);
  void     SetALTROSamplesLG (const Int_t nSamplesLG, Int_t *samplesLG);
  void     SetFALTROSamples  (const Int_t nSamples,   Int_t *samples) 
  { if(fDigitType==kTrigger) SetALTROSamplesLG(nSamples, samples) ; } 

  void     SetCalibAmp(Float_t amp) { fAmpCalib = amp  ; }
  Double_t GetCalibAmp()   const    { return fAmpCalib ; }

  void     Print(const Option_t* /*opt*/) const;
	
 private: 
	
  Float_t  fAmpFloat;     // Cell amplitude, float
  Int_t    fNSamples;     // Number of time samples, Low Gain for ALTRO, used also for FALTRO 
  Int_t   *fSamples;	  //[fNSamples], list of time bin constents, Low Gain for ALTRO, used also for FALTRO 
  Int_t    fNSamplesHG;   // Number of time samples, High Gain for ALTRO
  Int_t   *fSamplesHG;	  //[fNSamples], list of time bin constents, High Gain for ALTRO, used also for FALTRO 
	
  Int_t    fNprimary ;    // Number of primaries
  Int_t    fNMaxPrimary ; // Max Number of primaries
  Int_t   *fPrimary ;     //[fNMaxPrimary]  Array of primaries       
  Float_t *fDEPrimary ;   //[fNMaxPrimary]  Array of primary energy contributions
    
  Int_t    fNiparent ;    // Number of initial parents 
  Int_t    fNMaxiparent ; // Max Number of parents 
  Int_t   *fIparent ;     //[fNMaxiparent] Array of parents       
  Float_t *fDEParent;     //[fNMaxiparent]  Array of parent energy contributions
  Int_t    fMaxIter  ;    // Number to Increment Maxiparent, and MaxPrimary if default is not sufficient
  Float_t  fTime ;        // Calculated time  
  Float_t  fTimeR ;       // Earliest time: to be used by Digits2Raw
  //Fit quality parameters
  Float_t  fChi2;         // Fit Chi aquare	
  Int_t    fNDF;          // Fit Number of Degrees of Freedom
	
  Int_t    fDigitType;    // This is a trigger digit(0), HG (1) or LG (3)
  Float_t  fAmpCalib;     //! Calibrated energy

  ClassDef(AliEMCALDigit,6)   // Digit in EMCAL 
} ;

#endif //  ALIEMCALDIGIT_H
