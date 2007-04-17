#ifndef ALITPCCALROC_H
#define ALITPCCALROC_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTPCCalROC.h,v */

//////////////////////////////////////////////////
//                                              //
//  TPC calibration base class for one ROC      //
//                                              //
//////////////////////////////////////////////////

#include <TObject.h>
#include <TMath.h>
#include <AliTPCROC.h>
class TH1F;
class TH2F;
class TArrayI;
class TLinearFitter;
//_____________________________________________________________________________
class AliTPCCalROC : public TObject {

 public:
  
  AliTPCCalROC();
  AliTPCCalROC(UInt_t sector);
  AliTPCCalROC(const AliTPCCalROC &c);
 AliTPCCalROC &operator = (const AliTPCCalROC & param);
  virtual           ~AliTPCCalROC();  
  UInt_t        GetNrows() const               { return fNRows;};
  UInt_t        GetNchannels()       const     { return fNChannels;};
  UInt_t        GetNPads(UInt_t row)  const     { return (row<fNRows)? AliTPCROC::Instance()->GetNPads(fSector,row):0;};
  Float_t      GetValue(UInt_t row, UInt_t pad) const { return ( (row<fNRows) && (fIndexes[row]+pad)<fNChannels)? fData[fIndexes[row]+pad]: 0; };
  Float_t      GetValue(UInt_t channel) const { return  fData[channel]; };
  void         SetValue(UInt_t row, UInt_t pad, Float_t vd) { if ( row<fNRows && (fIndexes[row]+pad)<fNChannels)fData[fIndexes[row]+pad]= vd; };
  void         SetValue(UInt_t channel, Float_t vd) {fData[channel]= vd; };
  virtual void Draw(Option_t* option = "");
  //
  // algebra
  void Add(Float_t c1);
  void Multiply(Float_t c1);
  void Add(const AliTPCCalROC * roc, Double_t c1 = 1);
  void Multiply(const AliTPCCalROC * roc);   
  void Divide(const AliTPCCalROC * roc);   
  // statistic
  //
  Double_t GetMean(){return TMath::Mean(fNChannels, fData);}
  Double_t GetRMS() {return TMath::RMS(fNChannels, fData);}
  Double_t GetMedian() {return TMath::Median(fNChannels, fData);}
  Double_t GetLTM(Double_t *sigma=0, Double_t fraction=0.9);
  TH1F * MakeHisto1D(Float_t min=4, Float_t max=-4, Int_t type=0);     
  TH2F * MakeHisto2D(Float_t min=4, Float_t max=-4, Int_t type=0);   
  TH2F * MakeHistoOutliers(Float_t delta=4, Float_t fraction=0.7, Int_t mode=0);

  AliTPCCalROC * LocalFit(Int_t rowRadius, Int_t padRadius, AliTPCCalROC* ROCoutliers = 0, Bool_t robust = kFALSE);
  
  
  static void Test();
 protected:
  
  Double_t GetNeighbourhoodValue(TLinearFitter* fitterQ, Int_t row, Int_t pad, Int_t rRadius, Int_t pRadius, AliTPCCalROC* ROCoutliers, Bool_t robust);
  void GetNeighbourhood(TArrayI* &rowArray, TArrayI* &padArray, Int_t row, Int_t pad, Int_t rRadius, Int_t pRadius);
  
  
  UInt_t     fSector;          // sector number
  UInt_t     fNChannels;       // number of channels
  UInt_t     fNRows;           // number of rows
  const UInt_t* fIndexes;      //!indexes
  Float_t  *fData;            //[fNChannels] Data
  ClassDef(AliTPCCalROC,1)    //  TPC ROC calibration class

};

#endif
