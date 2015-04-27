#ifndef ALITRDCALROC_H
#define ALITRDCALROC_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//////////////////////////////////////////////////
//                                              //
//  TRD calibration base class for one ROC      //
//                                              //
//////////////////////////////////////////////////

#include <TObject.h>

class TArrayI;
class TArrayF;
class TH2F;
class TH1F;
class TGraph2D;
class TH2D;

//_____________________________________________________________________________
class AliTRDCalROC : public TObject 
{

 public:

  AliTRDCalROC();
  AliTRDCalROC(Int_t p, Int_t c);
  AliTRDCalROC(const AliTRDCalROC &c);
  virtual      ~AliTRDCalROC();
  AliTRDCalROC &operator=(const AliTRDCalROC &c);
  virtual void  Copy(TObject &c) const;

  Int_t         GetNrows() const                   { return fNrows; };
  Int_t         GetNcols() const                   { return fNcols; };

  Int_t         GetChannel(Int_t c, Int_t r) const { return r+c*fNrows; };
  Int_t         GetNchannels() const               { return fNchannels; };

  Float_t       GetValue(Int_t ich) const          { return (Float_t) fData[ich] / 10000;  };
  Float_t       GetValue(Int_t col, Int_t row)     { return GetValue(GetChannel(col,row)); };

  void          SetValue(Int_t ich, Float_t value) { fData[ich] = (UShort_t) (value * 10000); };
  void          SetValue(Int_t col, Int_t row, Float_t value) 
                                                   { SetValue(GetChannel(col,row), value);    };

  // statistic
  //
  Double_t GetMean(AliTRDCalROC * const outlierROC=0) const;
  Double_t GetMeanNotNull() const;
  Double_t GetRMS(AliTRDCalROC * const outlierROC=0) const;
  Double_t GetRMSNotNull() const;
  Double_t GetMedian(AliTRDCalROC * const outlierROC=0) const;
  Double_t GetLTM(Double_t *sigma=0, Double_t fraction=0.9, AliTRDCalROC * const outlierROC=0);

  // algebra
  Bool_t Add(Float_t c1);
  Bool_t Multiply(Float_t c1);
  Bool_t Add(const AliTRDCalROC * roc, Double_t c1 = 1);
  Bool_t Multiply(const AliTRDCalROC * roc);   
  Bool_t Divide(const AliTRDCalROC * roc);   

  // noise
  Bool_t Unfold();
  
  //Plots
  TH2F *   MakeHisto2D(Float_t min, Float_t max,Int_t type, Float_t mu = 1.0);
  TH1F *   MakeHisto1D(Float_t min, Float_t max,Int_t type, Float_t mu = 1.0);
  
 protected:

  Int_t     fPla;              //  Plane number
  Int_t     fCha;              //  Chamber number

  Int_t     fNrows;            //  Number of rows
  Int_t     fNcols;            //  Number of columns

  Int_t     fNchannels;        //  Number of channels
  UShort_t *fData;             //[fNchannels] Data

  ClassDef(AliTRDCalROC, 2)    //  TRD ROC calibration class

};

#endif
