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
//#include <TMath.h>
//#include <AliTPCROC.h>
#include "TLinearFitter.h"
#include "TVectorD.h"
//#include <iostream>



class AliTPCCalROC;
class AliTPCCalDet;
class TObjArray;
class TGraph;
class TH2F;
class TH1F;
class TCanvas;
class TTree;

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
  AliTPCCalROC *GetCalROC(Int_t sector) const {return fROC[sector]; };  
  void SetCalROC(AliTPCCalROC* roc, Int_t sector = -1);  
  virtual void	Draw(Option_t* option = "");
  //
  // algebra
  void Add(Float_t c1);   // add constant c1 to all channels of all ROCs
  void Multiply(Float_t c1);   // multiply each channel of all ROCs with c1
  void Add(const AliTPCCalPad * roc, Double_t c1 = 1);   // multiply AliTPCCalPad 'pad' by c1 and add each channel to the coresponing channel in all ROCs
  void Multiply(const AliTPCCalPad * pad);  // multiply each channel of all ROCs with the coresponding channel of 'pad'
  void Divide(const AliTPCCalPad * pad);    // divide each channel of all ROCs by the coresponding channel of 'pad'
  //
  Double_t GetMeanRMS(Double_t &rms);   // Calculates mean and RMS of all ROCs
  Double_t GetMean(AliTPCCalPad* outlierPad = 0);   // return mean of the mean of all ROCs
  Double_t GetRMS(AliTPCCalPad* outlierPad = 0) ;   // return mean of the RMS of all ROCs
  Double_t GetMedian(AliTPCCalPad* outlierPad = 0) ;   // return mean of the median of all ROCs
  Double_t GetLTM(Double_t *sigma=0, Double_t fraction=0.9, AliTPCCalPad* outlierPad = 0);   // return mean of the LTM and sigma of all ROCs
  TGraph  *MakeGraph(Int_t type=0, Float_t ratio=0.7);
  TH2F    *MakeHisto2D(Int_t side=0);
  TH1F    *MakeHisto1D(Float_t min=4, Float_t max=-4, Int_t type=0, Int_t side=0);  

  AliTPCCalPad* LocalFit(const char* padName, Int_t rowRadius, Int_t padRadius, AliTPCCalPad* Padoutliers = 0, Bool_t robust = kFALSE, Double_t chi2Threshold = 5, Double_t robustFraction = 0.7, Bool_t printCurrentSector = kFALSE) const;
  AliTPCCalPad* GlobalFit(const char* padName, AliTPCCalPad* Padoutliers = 0, Bool_t robust = kFALSE, Int_t fitType = 1, Double_t chi2Threshold = 5, Double_t robustFraction = 0.7, Double_t err=1, TObjArray *fitParArr=0x0, TObjArray *fitCovArr=0x0);
  
  void GlobalSidesFit(const AliTPCCalPad* PadOutliers, const char* fitFormula, TVectorD &fitParamSideA, TVectorD &fitParamSideC, TMatrixD &covMatrixSideA, TMatrixD &covMatrixSideC, Float_t &chi2SideA, Float_t &chi2SideC, AliTPCCalPad *pointError=0, Bool_t robust = kFALSE, Double_t robustFraction = 0.7);
  
  static AliTPCCalPad* CreateCalPadFit(const char* fitFormula, const TVectorD &fitParamSideA, const TVectorD &fitParamSideC);
    
  static TObjArray *CreateFormulaArray(const char *fitFormula);
  static void EvalFormulaArray(const TObjArray &arrFitFormulas, TVectorD &results,
                               const Int_t sec, const Int_t row, const Int_t pad);
  //
  // default report
  //
  static TCanvas * MakeReportPadSector(TTree *chain, const char* varName, const char*varTitle, const char *axisTitle, Float_t min, Float_t max, const char * cutUser="");
  static TCanvas * MakeReportPadSector2D(TTree *chain, const char* varName, const char*varTitle, const char *axisTitle, Float_t min, Float_t max, const char *cutUser="");

 protected:
  AliTPCCalROC *fROC[kNsec];                    //  Array of ROC objects which contain the values per pad
  ClassDef(AliTPCCalPad,1)                      //  TPC calibration class for parameters which are saved per pad
};

#endif
