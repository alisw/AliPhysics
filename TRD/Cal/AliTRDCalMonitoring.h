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

 protected:

  Int_t   fADCTresholds[6700];                 //  Threshold voltage for ADCs
  Float_t fDriftVelocity;                      //  Drift velocity from the monitor
  TString fGasComposition;		       //  Gas composition
  Float_t fEnvironmentTemperature;             //  Environment temperature

  Float_t fAnodeCurrentsMin[540];              //  Minimum anode current
  Float_t fAnodeCurrentsMax[540];              //  Maximum anode current
  Float_t fDriftCurrentsMin[540];              //  Minimum drift current
  Float_t fDriftCurrentsMax[540];              //  Maximum drift current
  Float_t fAnodeVoltagesMin[540];              //  Minimum anode voltage
  Float_t fAnodeVoltagesMax[540];              //  Maximum anode voltage
  Float_t fDriftVoltagesMin[540];              //  Minimum drift voltage
  Float_t fDriftVoltagesMax[540];              //  Maximum drift voltage

  Float_t fLVVoltage[360];                     //  Low voltage
  Float_t fLVCurrent[360];                     //  Low voltage current

  ClassDef(AliTRDCalMonitoring,1)              //  TRD calibration class for global TRD parameters

};

#endif
