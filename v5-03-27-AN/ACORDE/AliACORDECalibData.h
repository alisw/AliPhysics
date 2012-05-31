#ifndef ALIACORDECALIBDATA_H
#define ALIACORDECALIBDATA_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "TNamed.h"
#include "TH1D.h"
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
  Float_t* GetModulesActivity() const {return (float*)fModulesActivity;}
  Float_t GetModuleActivity(Int_t i) const {return fModulesActivity[i-1];}
 // TList*  GetHistos()const {return Hist;} 
  void SetRates(Float_t* Rt);
  void SetRate(Float_t rate, Int_t mod){fRates[mod-1]=rate;}
  void SetEfficiencies(Float_t* Eff);
  void SetEfficiency(Float_t eff, Int_t mod) {fEfficiencies[mod-1]=eff;}
  void AddHHits(TH1D  *Histo){fHits=(TH1D*)Histo->Clone("Hits");}// Hits
  void AddHTHits(TH1D *Histo){fTHits=(TH1D*)Histo->Clone("Total Hits");}//Total Hits 
  void AddHMultiHits(TH1D  *Histo){fMultiHits=(TH1D*)Histo->Clone("MultiHits");}//
  void AddHTMultiHits(TH1D *Histo){fTMultiHits=(TH1D*)Histo->Clone("Total Multi Hits");}
  void Draw(Option_t *option="");
  void SetModulesActivity(Float_t* Mac);
  void SetModuleActivity(Float_t mac,Int_t mod){fModulesActivity[mod-1]=mac;}


 protected:
  Float_t fEfficiencies[60];
  Float_t fRates[60];
  Float_t fModulesActivity[60];
  TH1D *fHits;
  TH1D *fTHits;
  TH1D *fMultiHits;
  TH1D *fTMultiHits;

  ClassDef(AliACORDECalibData,3)    // ACORDE Calibration data
};

#endif

