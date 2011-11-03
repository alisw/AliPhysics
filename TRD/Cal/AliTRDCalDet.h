#ifndef ALITRDCALDET_H
#define ALITRDCALDET_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TRD calibration class for parameters which are saved per detector        //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TNamed.h"

#include "../AliTRDgeometry.h"

class TH1F;
class TH2F;

class AliTRDCalDet : public TNamed {

 public:
 
  enum { kNplan = 6, kNcham = 5, kNsect = 18, kNdet = 540 };

  AliTRDCalDet();
  AliTRDCalDet(const Text_t* name, const Text_t* title);
  AliTRDCalDet(const AliTRDCalDet &c);   
  virtual      ~AliTRDCalDet();
  AliTRDCalDet &operator=(const AliTRDCalDet &c);

  virtual void  Copy(TObject &c) const;

  Float_t       GetValue(Int_t d) const          { return fData[d];  };
  Float_t       GetValue(Int_t p, Int_t c, Int_t s) const 
                                                 { return fData[AliTRDgeometry::GetDetector(p,c,s)];  };

  void          SetValue(Int_t d, Float_t value) { fData[d] = value; };
  void          SetValue(Int_t p, Int_t c, Int_t s, Float_t value) 
                                                 { fData[AliTRDgeometry::GetDetector(p,c,s)] = value; };

  // statistic
  Double_t GetMean(AliTRDCalDet * const outlierDet=0) const;
  Double_t GetRMS(AliTRDCalDet * const outlierDet=0) const;
  Double_t GetRMSRobust(Double_t robust=0.92) const;
  Double_t GetMedian(AliTRDCalDet * const outlierDet=0) const;
  Double_t GetLTM(Double_t * sigma=0, Double_t fraction=0.9, AliTRDCalDet * const outlierDet=0);     
  Double_t CalcMean(Bool_t wghtPads=kFALSE);
  Double_t CalcMean(Bool_t wghtPads, Int_t &calib);
  Double_t CalcRMS(Bool_t wghtPads=kFALSE);
  Double_t CalcRMS(Bool_t wghtPads, Int_t &calib);

  // Plot functions
  TH1F * MakeHisto1Distribution(Float_t min=4, Float_t max=-4, Int_t type=0);     
  TH1F * MakeHisto1DAsFunctionOfDet(Float_t min=4, Float_t max=-4, Int_t type=0);
  TH2F * MakeHisto2DCh(Int_t ch, Float_t min=4, Float_t max=-4, Int_t type=0);  
  TH2F * MakeHisto2DSmPl(Int_t sm, Int_t pl, Float_t min=4, Float_t max=-4, Int_t type=0); 

  // algebra functions
  void Add(Float_t c1);
  void Multiply(Float_t c1);
  void Add(const AliTRDCalDet * calDet, Double_t c1 = 1);
  void Multiply(const AliTRDCalDet * calDet);
  void Divide(const AliTRDCalDet * calDet);
    
 protected:

  Float_t  fData[kNdet];       //[kNdet] Data

  ClassDef(AliTRDCalDet,1)     //  TRD calibration class for parameters which are saved per detector

};

#endif
