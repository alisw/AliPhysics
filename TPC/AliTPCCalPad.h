#ifndef ALITPCCALPAD_H
#define ALITPCCALPAD_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TPC calibration class for parameters which are saved per pad                 //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TNamed.h"

class AliTPCCalROC;
class AliTPCCalDet;
class TObjArray;
class TGraph;
class TH2F;
class TH1F;

class AliTPCCalPad : public TNamed {
 public:
  enum { kNsec = 72 };
  AliTPCCalPad();
  AliTPCCalPad(const Text_t* name, const Text_t* title);
  AliTPCCalPad(const AliTPCCalPad &c);   
  AliTPCCalPad(TObjArray *arrayROC);
  virtual ~AliTPCCalPad();
  AliTPCCalPad &operator=(const AliTPCCalPad &c);
  virtual void     Copy(TObject &c) const;
  AliTPCCalROC *GetCalROC(Int_t sector) const { return fROC[sector]; };  
  //
  // algebra
  void Add(Float_t c1);
  void Multiply(Float_t c1);
  void Add(const AliTPCCalPad * roc, Double_t c1 = 1);
  void Multiply(const AliTPCCalPad * pad);
  void Divide(const AliTPCCalPad * pad);
  //
  Double_t GetMeanRMS(Double_t &rms);
  Double_t GetMean();
  Double_t GetRMS() ;
  Double_t GetMedian() ;
  Double_t GetLTM(Double_t *sigma=0, Double_t fraction=0.9);
  TGraph  *MakeGraph(Int_t type=0, Float_t ratio=0.7);
  TH2F    *MakeHisto2D(Int_t side=0);
  TH1F    *MakeHisto1D(Float_t min=4, Float_t max=-4, Int_t type=0);
 protected:
  AliTPCCalROC *fROC[kNsec];                    //  Array of ROC objects which contain the values per pad
  ClassDef(AliTPCCalPad,1)                      //  TPC calibration class for parameters which are saved per pad
};

#endif
