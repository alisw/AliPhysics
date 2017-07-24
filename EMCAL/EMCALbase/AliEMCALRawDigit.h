#ifndef ALIEMCALRAWDIGIT_H
#define ALIEMCALRAWDIGIT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////////////////////////////////////
///  
/// \class AliEMCALRawDigit
/// \ingroup EMCALbase
/// \brief EMCal raw digits object  
///
///  A Digit containing raw samples, for trigger.
///
/// \author R. Guernane, LPSC-IN2P3-CNRS
///
////////////////////////////////////////////////////////////////////////////////

#include "TObject.h" 

class AliEMCALRawDigit : public TObject 
{
  
public:
  
  AliEMCALRawDigit();
  AliEMCALRawDigit(Int_t id, Int_t timeSamples[], Int_t nSamples);
  
  virtual ~AliEMCALRawDigit();
  
  void Clear(Option_t *);
  
  Bool_t  IsSortable() const { return kTRUE;}
  
  Int_t   Compare(const TObject* obj) const;
  
  void    SetId(Int_t id)           { fId        = id ; }
  void    SetAmplitude(Float_t amp) { fAmplitude = amp  ; }
  void    SetTime(Float_t time)     { fTime      = time ; }
  
  void    SetTimeSamples(const Int_t timeSamples[], const Int_t nSamples);
  
  Int_t   GetId()        const { return fId        ; }	
  Float_t GetAmplitude() const { return fAmplitude ; }
  Float_t GetTime()      const { return fTime      ; }
  Int_t   GetNSamples()  const { return fNSamples  ; }
  
  Bool_t  GetTimeSample(const Int_t iSample, Int_t& timeBin, Int_t& amp) const;
  Bool_t  GetMaximum(Int_t& amplitude, Int_t& time) const;
  
  virtual void Print(const Option_t* opt) const;
  
protected: 
  
  AliEMCALRawDigit           (const AliEMCALRawDigit &cd); // Not implemented
  AliEMCALRawDigit &operator=(const AliEMCALRawDigit &cd); // Not implemented
  
  Int_t   fId;            ///< Absolute id
  
  Int_t   fNSamples;      ///< Number of time samples
  
  /// Samples stored here 
  Int_t*  fSamples;	      //[fNSamples]
  
  Float_t fAmplitude;     ///< digit amplitude
  
  Float_t fTime;          ///< digit time 
  
  /// \cond CLASSIMP
  ClassDef(AliEMCALRawDigit,1)  ;
  /// \endcond
  
};

#endif // ALIEMCALRAWDIGIT_H

