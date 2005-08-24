#ifndef ALIPHOSCALIBDATA_H
#define ALIPHOSCALIBDATA_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////
//  class for PHOS calibration                 //
////////////////////////////////////////////////

#include "TNamed.h"
#include "AliPHOS.h"

class AliPHOSCalibData: public TNamed {

 public:
  AliPHOSCalibData();
  AliPHOSCalibData(const char* name);
  AliPHOSCalibData(const AliPHOSCalibData &calibda);
  AliPHOSCalibData& operator= (const AliPHOSCalibData &calibda);
  virtual ~AliPHOSCalibData();
  void Reset();
  virtual void Print(Option_t *option = "") const; 
  //
  Float_t GetADCchannelEmc(Int_t module, Int_t column, Int_t row) const;
  Float_t GetADCpedestalEmc(Int_t module, Int_t column, Int_t row) const;
  //
  void SetADCchannelEmc(Int_t module, Int_t column, Int_t row, Float_t value);
  void SetADCpedestalEmc(Int_t module, Int_t column, Int_t row, Float_t value);

 protected:
  Float_t  fADCchannelEmc[5][56][64] ;  // width of one ADC channel in GeV ([mod][col][row])
  Float_t  fADCpedestalEmc[5][56][64] ; // value of the EMC ADC pedestal ([mod][col][row])
  //
  ClassDef(AliPHOSCalibData,1)    // PHOS Calibration data
};

#endif
