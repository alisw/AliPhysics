#ifndef ALIEMCALDIGIT_H
#define ALIEMCALDIGIT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////////////////////////////////////
///  
/// \class AliEMCALDigit
/// \ingroup EMCALbase
/// \brief EMCal digits object  
///
///  A Digit is the sum of the energy lost in an EMCAL Tower
///  It also stores information on Primary, and enterring particle
///  tracknumbers Digits are created using AliEMCALSDigitizer, followed
///  by AliEMCALDigitizer. Based on : AliPHOSDigit
///
/// \author Sahal Yacoob (LBL)
/// \author Gustavo Conesa Balbastre <Gustavo.Conesa.Balbastre@cern.ch>, LPSC-IN2P3-CNRS
///
////////////////////////////////////////////////////////////////////////////////

// --- ROOT system ---
#include "TObject.h" 

// --- AliRoot header files ---
#include "AliDigitNew.h"

using std::ostream;

class AliEMCALDigit : public AliDigitNew 
{

  friend ostream& operator << ( ostream& , const AliEMCALDigit&) ;

 public:
  
  AliEMCALDigit() ;
  AliEMCALDigit(Int_t primary, Int_t iparent, Int_t id, Float_t digEnergy, 
                Float_t time, Int_t type,Int_t index = -1, 
                Float_t chi2=0, Int_t ndf=0, Float_t dE = 0) ;
  
  AliEMCALDigit(const AliEMCALDigit & digit) ;
  
  virtual ~AliEMCALDigit() ;

  Bool_t operator==(const AliEMCALDigit &rValue) const;
  AliEMCALDigit operator+(const AliEMCALDigit &rValue) ;
  AliEMCALDigit operator*(Float_t factor) ; 
  AliEMCALDigit& operator = (const AliEMCALDigit & digit) ;
  
  enum  digitType{kUnknown=-1, kHG=0, kLG=1, kLGnoHG=2, kTrigger=3, kEmbedded = 4};

  void     Clear(Option_t*) ;	
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

  //
  // Raw time sample
  //
  
  // ALTRO
  Int_t    GetNALTROSamplesLG() const {if(fDigitType==kLG)      return fNSamples;   else return 0 ; }
  Bool_t   GetALTROSampleLG(const Int_t iSample, Int_t& timeBin, Int_t& amp) const;
  Int_t    GetNALTROSamplesHG() const {if(fDigitType==kHG)      return fNSamplesHG; else return 0 ; }
  Bool_t   GetALTROSampleHG(const Int_t iSample, Int_t& timeBin, Int_t& amp) const;
  
  // FALTRO, trigger. Same data members as Low Gain	
  Int_t    GetNFALTROSamples()  const {if(fDigitType==kTrigger) return fNSamples;   else return 0 ; }
  Bool_t   GetFALTROSample(const Int_t iSample, Int_t& timeBin, Int_t& amp)  const ;
	
  void     SetALTROSamplesHG (const Int_t nSamplesHG, Int_t *samplesHG);
  void     SetALTROSamplesLG (const Int_t nSamplesLG, Int_t *samplesLG);
  void     SetFALTROSamples  (const Int_t nSamples,   Int_t *samples) 
  { if(fDigitType==kTrigger) SetALTROSamplesLG(nSamples, samples) ; } 

  //
  // Primary/Parents array creation
  // Used at analysis level while reclusterizing
  //
  void SetListOfPrimaries(Int_t npri, Int_t * prilist, Float_t * edepList) ;  
  void SetListOfParents  (Int_t npar, Int_t * parlist, Float_t * edepList) ;  
  
  //
  // Other
  //
  void     SetCalibAmp(Float_t amp) { fAmpCalib = amp  ; }
  Double_t GetCalibAmp()   const    { return fAmpCalib ; }

  void     Print(const Option_t* /*opt*/) const;
	
 private: 
	
  Float_t  fAmpFloat;     ///< Cell amplitude, float
  
  Int_t    fNSamples;     ///< Number of time samples, Low Gain for ALTRO, used also for FALTRO 
  
  /// List of time bin constents, Low Gain for ALTRO, used also for FALTRO
  Int_t   *fSamples;	    //[fNSamples]  
  
  Int_t    fNSamplesHG;   ///< Number of time samples, High Gain for ALTRO
  
  /// List of time bin constents, High Gain for ALTRO, used also for FALTRO 
  Int_t   *fSamplesHG;	  //[fNSamples] 
	
  Int_t    fNprimary ;    ///< Number of primaries
  
  Int_t    fNMaxPrimary ; ///< Max Number of primaries
  
  /// Array of primary labels
  Int_t   *fPrimary ;     //[fNMaxPrimary]       
  
  /// Array of primary energy contributions
  Float_t *fDEPrimary ;   //[fNMaxPrimary]  
    
  Int_t    fNiparent ;    ///< Number of initial parents 
  
  Int_t    fNMaxiparent ; ///< Max Number of parents 
  
  /// Array of parents labels
  Int_t   *fIparent ;     //[fNMaxiparent]       
  
  /// Array of parent energy contributions
  Float_t *fDEParent;     //[fNMaxiparent] 
  
  Int_t    fMaxIter  ;    ///< Number to Increment Maxiparent, and MaxPrimary if default is not sufficient
  
  Float_t  fTime ;        ///< Calculated time  
  
  Float_t  fTimeR ;       ///< Earliest time: to be used by Digits2Raw
  
  Float_t  fChi2;         ///< Fit quality parameter, chi square	
  
  Int_t    fNDF;          ///< Fit quality parameter, number of Degrees of Freedom
	
  Int_t    fDigitType;    ///< This is a trigger digit(0), HG (1) or LG (3)
  
  Float_t  fAmpCalib;     ///< Calibrated energy

  /// \cond CLASSIMP
  ClassDef(AliEMCALDigit,7)  ;
  /// \endcond

} ;

// inline definitions

///
/// Compares two digits with respect to its Id
/// to sort according increasing Id
//____________________________________________________________________________
inline Int_t AliEMCALDigit::Compare(const TObject * obj) const
{
  AliEMCALDigit * digit = (AliEMCALDigit *)obj ;

  Int_t iddiff = fId - digit->GetId() ;

  if ( iddiff > 0 )
    return 1 ;
  else if ( iddiff < 0 )
    return -1 ;
  else
    return 0 ;
}

#endif //  ALIEMCALDIGIT_H
