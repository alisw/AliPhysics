#ifndef AliTRDCALMONITORING_H
#define AliTRDCALMONITORING_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TRD calibration class for monitoring data                                //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TNamed.h"

class AliTRDCalMonitoring : public TNamed {
  public:
    AliTRDCalMonitoring();
    AliTRDCalMonitoring(const Text_t* name, const Text_t* title);
    virtual ~AliTRDCalMonitoring() {};

    //void SetSamplingFrequency(Float_t freq)                { fSamplingFrequency = freq; };
    //Float_t GetSamplingFrequency()                   const { return fSamplingFrequency; };

  protected:
    Int_t fADCTresholds[6700];
    Float_t fDriftVelocity;                      // Drift velocity from the monitor
    TString fGasComposition;			 // Gas composition
    Float_t fEnvironmentTemperature;

    //Float_t fMCMTemperature[6700];

    Float_t fAnodeCurrentsMin[540];
    Float_t fAnodeCurrentsMax[540];
    Float_t fDriftCurrentsMin[540];
    Float_t fDriftCurrentsMax[540];
    Float_t fAnodeVoltagesMin[540];
    Float_t fAnodeVoltagesMax[540];
    Float_t fDriftVoltagesMin[540];
    Float_t fDriftVoltagesMax[540];

    Float_t fLVVoltage[360];
    Float_t fLVCurrent[360];

    void Init();

    ClassDef(AliTRDCalMonitoring,1)                      //  TRD calibration class for global TRD parameters
};

#endif
