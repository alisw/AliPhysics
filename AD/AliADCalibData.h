#ifndef ALIADCALIBDATA_H
#define ALIADCALIBDATA_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "TNamed.h"
#include "TH1D.h"
class AliADCalibData: public TNamed {

 public:
  AliADCalibData();
  AliADCalibData(const char* name);
  
  AliADCalibData(const AliADCalibData &calibda);
  AliADCalibData& operator= (const AliADCalibData &calibda);
  virtual ~AliADCalibData();
  void Reset();

  Float_t* GetEfficiencies() const { return (float*)fEfficiencies; }
  Float_t  GetEfficiency(Int_t i) const { return fEfficiencies[i-1];}
  Float_t* GetRates() const {return (float*)fRates;}
  Float_t GetRate(Int_t i) const {return fRates[i-1];}
  Float_t* GetModulesActivity() const {return (float*)fModulesActivity;}
  Float_t GetModuleActivity(Int_t i) const {return fModulesActivity[i-1];}
 // TList*  GetHistos()const {return Hist;} 
  void SetRates(Float_t* Rt);
  void SetRate(Float_t rate, Int_t mod){fRates[mod-1]=rate;}
  void SetEfficiencies(Float_t* Eff);
  void SetEfficiency(Float_t eff, Int_t mod) {fEfficiencies[mod-1]=eff;}
  void Draw(Option_t *option="");
  void SetModulesActivity(Float_t* Mac);
  void SetModuleActivity(Float_t mac,Int_t mod){fModulesActivity[mod-1]=mac;}


 protected:
  Float_t fEfficiencies[60];
  Float_t fRates[60];
  Float_t fModulesActivity[60];

  ClassDef(AliADCalibData,1)    // AD Calibration data
};

#endif

