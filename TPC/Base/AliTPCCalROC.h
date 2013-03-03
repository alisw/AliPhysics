#ifndef ALITPCCALROC_H
#define ALITPCCALROC_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTPCCalROC.h -1   */

//////////////////////////////////////////////////
//                                              //
//  TPC calibration base class for one ROC      //
//                                              //
//////////////////////////////////////////////////
#include <TObject.h>
#include <TMath.h>
#include <AliTPCROC.h>
#include "TLinearFitter.h"

class TH1F;
class TH2F;
class TArrayI;
//_____________________________________________________________________________
class AliTPCCalROC : public TNamed {

 public:
  
  AliTPCCalROC();
  AliTPCCalROC(UInt_t sector);
  AliTPCCalROC(const AliTPCCalROC &c);
 AliTPCCalROC &operator = (const AliTPCCalROC & param);
  virtual           ~AliTPCCalROC();  
  UInt_t        GetSector() const { return fSector;}
  UInt_t        GetNrows() const               { return fNRows;};
  UInt_t        GetNchannels()       const     { return fNChannels;};
  UInt_t        GetNPads(UInt_t row)  const     { return (row<fNRows)? AliTPCROC::Instance()->GetNPads(fSector,row):0;};
  Float_t      GetValue(UInt_t row, UInt_t pad) const { return ( (row<fNRows) && (fkIndexes[row]+pad)<fNChannels)? fData[fkIndexes[row]+pad]: 0; };
  Float_t      GetValue(UInt_t channel) const { return  fData[channel]; };
  void         SetValue(UInt_t row, UInt_t pad, Float_t vd) { if ( row<fNRows && (fkIndexes[row]+pad)<fNChannels)fData[fkIndexes[row]+pad]= vd; };
  void         SetValue(UInt_t channel, Float_t vd) {fData[channel]= vd; };
  virtual void Draw(Option_t* option = "");
  //
  // algebra
  void Add(Float_t c1); // add c1 to each channel of the ROC
  void Multiply(Float_t c1); // multiply each channel of the ROC with c1
  void Add(const AliTPCCalROC * roc, Double_t c1 = 1);  // multiply AliTPCCalROC roc by c1 and add each channel to the coresponing channel in the ROC
  void Multiply(const AliTPCCalROC * roc);   // multiply each channel of the ROC with the coresponding channel of 'roc'
  void Divide(const AliTPCCalROC * roc);   // divide each channel of the ROC by the coresponding value of 'roc'
  // statistic
  //
  Double_t GetMean(AliTPCCalROC *const outlierROC = 0) const;
  Double_t GetRMS(AliTPCCalROC *const outlierROC = 0) const;
  Double_t GetMedian(AliTPCCalROC *const outlierROC = 0) const;
  Double_t GetLTM(Double_t *const sigma=0, Double_t fraction=0.9, AliTPCCalROC *const outlierROC = 0);
  TH1F * MakeHisto1D(Float_t min=4, Float_t max=-4, Int_t type=0);     
  TH2F * MakeHisto2D(Float_t min=4, Float_t max=-4, Int_t type=0);   
  TH2F * MakeHistoOutliers(Float_t delta=4, Float_t fraction=0.7, Int_t mode=0);

  AliTPCCalROC * LocalFit(Int_t rowRadius, Int_t padRadius, AliTPCCalROC* ROCoutliers = 0, Bool_t robust = kFALSE, Double_t chi2Threshold = 5, Double_t robustFraction = 0.7);
  void GlobalFit(const AliTPCCalROC* ROCoutliers, Bool_t robust, TVectorD &fitParam, TMatrixD &covMatrix, Float_t & chi2, Int_t fitType = 1, Double_t chi2Threshold = 5, Double_t robustFraction = 0.7, Double_t err=1);
  
  static AliTPCCalROC* CreateGlobalFitCalROC(TVectorD &fitParam, Int_t sector);
  
  static void Test();
 protected:
  
  Double_t GetNeighbourhoodValue(TLinearFitter* fitterQ, Int_t row, Int_t pad, Int_t rRadius, Int_t pRadius, AliTPCCalROC *const ROCoutliers, Bool_t robust, Double_t chi2Threshold, Double_t robustFraction);
  void GetNeighbourhood(TArrayI* &rowArray, TArrayI* &padArray, Int_t row, Int_t pad, Int_t rRadius, Int_t pRadius);
  
  UInt_t     fSector;          // sector number
  UInt_t     fNChannels;       // number of channels
  UInt_t     fNRows;           // number of rows
  const UInt_t* fkIndexes;      //!indexes
  Float_t  *fData;            //[fNChannels] Data
  ClassDef(AliTPCCalROC,2)    //  TPC ROC calibration class

};

#endif
