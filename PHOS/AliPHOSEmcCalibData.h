#ifndef ALIPHOSEMCCALIBDATA_H
#define ALIPHOSEMCCALIBDATA_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////
//  class for EMC calibration                 //
////////////////////////////////////////////////

#include "TNamed.h"

class AliPHOSEmcCalibData: public TNamed {

 public:
  AliPHOSEmcCalibData();
  AliPHOSEmcCalibData(const char* name);
  AliPHOSEmcCalibData(const AliPHOSEmcCalibData &calibda);
  AliPHOSEmcCalibData& operator= (const AliPHOSEmcCalibData &calibda);
  virtual ~AliPHOSEmcCalibData();
  void Reset();
  virtual void Print(Option_t *option = "") const; 
  //
  Float_t GetADCchannelEmc(Int_t module, Int_t column, Int_t row) const;
  Float_t GetADCpedestalEmc(Int_t module, Int_t column, Int_t row) const;
  //
  void SetADCchannelEmc(Int_t module, Int_t column, Int_t row, Float_t value);
  void SetADCpedestalEmc(Int_t module, Int_t column, Int_t row, Float_t value);

 protected:
  Float_t  fADCchannelEmc[5][56][64] ;  // width of one EMC ADC channel in GeV ([mod][col][row])
  Float_t  fADCpedestalEmc[5][56][64] ; // value of the EMC ADC pedestal ([mod][col][row])
  //
  ClassDef(AliPHOSEmcCalibData,1)    // PHOS EMC calibration data
};

#endif
