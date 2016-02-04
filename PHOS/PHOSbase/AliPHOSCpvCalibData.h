#ifndef ALIPHOSCPVCALIBDATA_H
#define ALIPHOSCPVCALIBDATA_H

/* Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////
//  class for CPV calibration                 //
////////////////////////////////////////////////

#include "TNamed.h"
#include "AliPHOSCpvParam.h"

class AliPHOSCpvCalibData: public TNamed {

 public:
  AliPHOSCpvCalibData();
  AliPHOSCpvCalibData(const char* name);
  AliPHOSCpvCalibData(const AliPHOSCpvCalibData &calibda);
  AliPHOSCpvCalibData& operator= (const AliPHOSCpvCalibData &calibda);
  virtual ~AliPHOSCpvCalibData();
  void Reset();
  virtual void Print(Option_t *option = "") const; 
  //
  Float_t GetADCchannelCpv(Int_t module, Int_t column, Int_t row) const;
  Float_t GetADCpedestalCpv(Int_t module, Int_t column, Int_t row) const;
  //
  void SetADCchannelCpv(Int_t module, Int_t column, Int_t row, Float_t value);
  void SetADCpedestalCpv(Int_t module, Int_t column, Int_t row, Float_t value);

 protected:
  Float_t  fADCchannelCpv[AliPHOSCpvParam::kNModules][AliPHOSCpvParam::kPadPcX][AliPHOSCpvParam::kPadPcY];  // width of one CPV ADC channel ([mod][col][row])
  Float_t  fADCpedestalCpv[AliPHOSCpvParam::kNModules][AliPHOSCpvParam::kPadPcX][AliPHOSCpvParam::kPadPcY]; // value of the CPV ADC pedestal ([mod][col][row])
  //
  ClassDef(AliPHOSCpvCalibData,2)    // CPV Calibration data

};

#endif
