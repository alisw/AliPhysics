#ifndef ALIPHOSCPVCALIBDATA_H
#define ALIPHOSCPVCALIBDATA_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////
//  class for CPV calibration                 //
////////////////////////////////////////////////

#include "TNamed.h"

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
  Float_t  fADCchannelCpv[5][56][128];  // width of one CPV ADC channel ([mod][col][row])
  Float_t  fADCpedestalCpv[5][56][128]; // value of the CPV ADC pedestal ([mod][col][row])
  //
  ClassDef(AliPHOSCpvCalibData,1)    // CPV Calibration data

};

#endif
