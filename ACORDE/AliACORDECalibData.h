#ifndef ALIACORDECALIBDATA_H
#define ALIACORDECALIBDATA_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "TNamed.h"
class AliACORDECalibData: public TNamed {

 public:
  AliACORDECalibData();
  AliACORDECalibData(const char* name);
  AliACORDECalibData(const AliACORDECalibData &calibda);
  AliACORDECalibData& operator= (const AliACORDECalibData &calibda);
  virtual ~AliACORDECalibData();
  void Reset();

  Float_t* GetEfficiencies() const { return (float*)fEfficiencies; }
  Float_t  GetEfficiency(Int_t i) const { return fEfficiencies[i-1];}
  Float_t* GetRates() const {return (float*)fRates;}
  Float_t GetRate(Int_t i) const {return fRates[i-1];}
  void SetRates(Float_t* Rt);
  void SetRate(Float_t rate, Int_t mod){fRates[mod-1]=rate;}
  void SetEfficiencies(Float_t* Eff);
  void SetEfficiency(Float_t eff, Int_t mod) {fEfficiencies[mod-1]=eff;}

 protected:
  Float_t fEfficiencies[60];
  Float_t fRates[60];

  ClassDef(AliACORDECalibData,1)    // ACORDE Calibration data
};

#endif

