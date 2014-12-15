#ifndef ALITRDCALPAD_H
#define ALITRDCALPAD_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TRD calibration class for parameters which are saved per pad             //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TNamed.h"

class AliTRDCalROC;
class AliTRDCalDet;
class TH2F;
class TH1F;

class AliTRDCalPad : public TNamed 
{

 public:
 
  enum { kNplan = 6, kNcham = 5, kNsect = 18, kNdet = 540 };

  AliTRDCalPad();
  AliTRDCalPad(const Text_t* name, const Text_t* title);
  AliTRDCalPad(const AliTRDCalPad &c);   
  virtual            ~AliTRDCalPad();
  AliTRDCalPad        &operator=(const AliTRDCalPad &c);

  virtual void        Copy(TObject &c) const;

  static  Int_t       GetDet(Int_t p, Int_t c, Int_t s) { return p+c*kNplan+s*kNplan*kNcham; };

  AliTRDCalROC       *GetCalROC(Int_t d) const          { return fROC[d]; };
  AliTRDCalROC       *GetCalROC(Int_t p, Int_t c, Int_t s) const
                                                        { return fROC[GetDet(p,c,s)]; };
  
  Bool_t              ScaleROCs(const AliTRDCalDet* values);

  void                SetCalROC(Int_t det, AliTRDCalROC *calroc);

  // Statistic
  Double_t GetMeanRMS(Double_t &rms, const AliTRDCalDet *calDet = 0, Int_t type = 0);
  Double_t GetMean(const AliTRDCalDet *calDet = 0, Int_t type = 0, AliTRDCalPad* const outlierPad = 0);
  Double_t GetRMS(const AliTRDCalDet *calDet = 0, Int_t type = 0, AliTRDCalPad* const outlierPad = 0);
  Double_t GetMedian(const AliTRDCalDet *calDet = 0, Int_t type = 0, AliTRDCalPad* const outlierPad = 0);
  Double_t GetLTM(Double_t *sigma=0, Double_t fraction=0.9
                , const AliTRDCalDet *calDet = 0, Int_t type = 0, AliTRDCalPad* const outlierPad = 0);

  // Plot functions
  TH1F    *MakeHisto1D(const AliTRDCalDet *calDet = 0, Int_t typedet=0, Float_t min=4, Float_t max=-4,Int_t type=0);
  TH2F    *MakeHisto2DSmPl(Int_t sm, Int_t pl, const AliTRDCalDet *calDet = 0, Int_t typedet=0, Float_t min=4, Float_t max=-4,Int_t type=0);
  TH2F    *MakeHisto2DCh(Int_t ch, const AliTRDCalDet *calDet = 0, Int_t typedet=0, Float_t min=4, Float_t max=-4,Int_t type=0);

  // Algebra functions
  Bool_t Add(Float_t c1);
  Bool_t Multiply(Float_t c1);
  Bool_t Add(const AliTRDCalPad * pad, Double_t c1 = 1, const AliTRDCalDet * calDet1 = 0, const AliTRDCalDet *calDet2 = 0, Int_t type = 0);
  Bool_t Multiply(const AliTRDCalPad * pad, const AliTRDCalDet * calDet1 = 0, const AliTRDCalDet *calDet2 = 0, Int_t type = 0);
  Bool_t Divide(const AliTRDCalPad * pad, const AliTRDCalDet * calDet1 = 0, const AliTRDCalDet *calDet2 = 0, Int_t type = 0);

 protected:

  AliTRDCalROC *fROC[kNdet];  //  Array of ROC objects which contain the values per pad

  ClassDef(AliTRDCalPad,1)    //  TRD calibration class for parameters which are saved per pad

};

#endif
