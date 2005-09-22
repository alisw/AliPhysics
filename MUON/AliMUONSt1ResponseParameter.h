#ifndef ALI_MUON_ST1_RESPONSE_PARAMETER_H
#define ALI_MUON_ST1_RESPONSE_PARAMETER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
// Revision of includes 07/05/2004

/// \ingroup sim
/// \class AliMUONSt1ResponseParameter
/// \brief Describes a set of filters to be applied to a digital value
///
/// Describes a set of filters to be applied to a digital value
/// in order to simulate electronics characteristics 
/// (pedestal, noise, sticky bits, etc....)
/// Threshold levels for the MANU zero supression algorithm are included.
///
/// Authors: David Guez, Ivana Hrivnacova, Marion MacCormick; IPN Orsay

#include <TNamed.h>

class TString;

class AliMUONSt1ResponseParameter : public TNamed 
{
  public:
   typedef enum {kNone,kValue,kGauss,kFile} TMode;
   typedef struct {Double_t mean; Double_t sigma;} TGaussParam;

  public:
    AliMUONSt1ResponseParameter();
    AliMUONSt1ResponseParameter(const TString& name,const TString& title);
    virtual ~AliMUONSt1ResponseParameter();
    
    void SetState(Bool_t state) ;
    void SetPedestal(Double_t val);
    void SetPedestal(Double_t mean,Double_t sigma);
    void SetPedestal(const TString& fileName);
    void UnSetPedestal();
    void SetNoise(Double_t val);
    void SetNoise(Double_t mean,Double_t sigma);
    void SetNoise(const TString& fileName);
    void SetNofSigma(Int_t nofSigma);
    void SetStickyBitOn (Int_t bit,Int_t val=1);
    void SetStickyBitOff(Int_t bit,Int_t val=1);
    Int_t ApplyPedestal(Int_t base,Int_t GC) const;
    Int_t ApplyStickyBits(Int_t base) const;
    Bool_t HasPedestal() const {return fPedestalMode != kNone;}
    Bool_t GetState() const {return fState;}
        
 private:
    static const Int_t fgkNofChannels=64;  // number of channels
    typedef union {
      //Double_t values[fgkNofChannels];
      Double_t values[64];
      Double_t value;
      TGaussParam gauss;
    } TParam;

    Double_t Choose(TMode mode,TParam param,Int_t GC) const;
    TMode  fPedestalMode;  // mode for pedestal values
    TParam fPedestalParam; // pedestal access parameters
    TMode  fNoiseMode;     // mode for noise values
    TParam fNoiseParam;    // noise access parameters 
    Int_t  fNofSigma;      // No of sigma for threshold (zero supression)
    Bool_t fState;         // is the element on/off
    Int_t  fStickyOn;      // which bits are always on (mask) .
    Int_t  fStickyOff;     // which bits are always off (mask).

  ClassDef(AliMUONSt1ResponseParameter,1) // electronics parmeters for Response
};
#endif //ALI_MUON_ST1_RESPONSE_PARAMETER_H
