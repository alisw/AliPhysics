#ifndef AliTRDCALGLOBALS_H
#define AliTRDCALGLOBALS_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TRD calibration class for global TRD parameters        //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TNamed.h"

class AliTRDCalGlobals : public TNamed {
  public:
    AliTRDCalGlobals();
    AliTRDCalGlobals(const Text_t* name, const Text_t* title);
    virtual ~AliTRDCalGlobals() {};
    
    void SetSamplingFrequency(Float_t freq)                { fSamplingFrequency = freq; };
    Float_t GetSamplingFrequency()                   const { return fSamplingFrequency; };
    
    void SetNumberOfTimeBins(Int_t value)     { fNumberOfTimeBins = value; };
    Int_t GetNumberOfTimeBins()         const { return fNumberOfTimeBins; };
  
  protected:
    Float_t fSamplingFrequency;                  // Sampling Frequency in MHz
    Int_t fNumberOfTimeBins;                     // Number of timebins  
    
    void Init();
    
    ClassDef(AliTRDCalGlobals,1)                      //  TRD calibration class for global TRD parameters
};

#endif
