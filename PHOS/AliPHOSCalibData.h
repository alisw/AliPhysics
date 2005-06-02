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
  Float_t GetADCchannelEmc(Int_t module, Int_t column, Int_t row) const {return fADCchannelEmc[module][column][row];}
  Float_t GetADCpedestalEmc(Int_t module, Int_t column, Int_t row) const {return fADCpedestalEmc[module][column][row];}
  //
  void SetADCchannelEmc(Int_t module, Int_t column, Int_t row, Float_t value)  {fADCchannelEmc[module][column][row] = value;}
  void SetADCpedestalEmc(Int_t module, Int_t column, Int_t row, Float_t value) {fADCpedestalEmc[module][column][row] = value;}

 protected:
  Float_t  fADCchannelEmc[5][64][56] ;           // width of one ADC channel in GeV
  Float_t  fADCpedestalEmc[5][64][56] ;          // value of the EMC ADC pedestal
  //
  ClassDef(AliPHOSCalibData,1)    // PHOS Calibration data
};

#endif
