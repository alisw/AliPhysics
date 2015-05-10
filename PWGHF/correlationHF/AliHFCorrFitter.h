#ifndef ALIHFCORRFITTER_H
#define ALIHFCORRFITTER_H

/* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
/*____________________________________________________________
| Class for performing the fit of azimuthal correlations           |      
| Example of its usage in the macro PWGHF/correlationHF/FitPlots.C |
|                                                                  |
|  Author: S. Kar (somnath.kar@cern.ch),                           |
|          A. Rossi (andrea.rossi@cern.ch)                         |
|                                                                  |
|_____________________________________________________________*/


#include <TFile.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLatex.h>

class TH1F;
class TCanvas;
class TF1;


class AliHFCorrFitter{

 public:
    
  // enums
  enum FunctionType{kConstwoGaus = 1, kTwoGausPeriodicity =2, kConstThreeGausPeriodicity = 3, ConstThreeGausPeriodicityAS =4, kv2Modulation =5, kTwoGausPeriodicityWithv2Modulation = 6};  

  // constructors
  AliHFCorrFitter();
  AliHFCorrFitter(const TH1F* histoToFit, Double_t min, Double_t max);
  virtual ~AliHFCorrFitter();
  AliHFCorrFitter(const AliHFCorrFitter &source);
  AliHFCorrFitter& operator=(const AliHFCorrFitter &cfit);


  //------------------------Setters---------->
  void SetHisto(const TH1F *histoToFit);
  void SetFunction();
  void SetRange(Double_t a,Double_t b){fMin=a;fMax=b;}
  void SetFuncType(FunctionType fittype){fTypeOfFitfunc = fittype;}
  void SetFixBasetype(Int_t fixbasetype){fFixBase=fixbasetype;}
  void SetFixMeanType(Int_t fixmeantype){fFixMean=fixmeantype;}
  void SetBaselineExt(Double_t baseval){
    fBaseline=baseval;
    fFixBase=6;
  }
  void SetBaselineEstimationRange(Double_t min, Double_t max){
        fMinBaselineRange = min; // minimum for estimation of the baseline
        fMaxBaselineRange = max; // max for estimation of the baseline
  }
  //---------------------Getters----------->
  Double_t GetNSSigma(){return fFit->GetParameter("NS #sigma");}
  Double_t GetNSYield(){return fFit->GetParameter("NS Y");}
  Double_t GetASSigma(){return fFit->GetParameter("AS #sigma");}
  Double_t GetASYield(){return fFit->GetParameter("AS Y");}
  Double_t GetPedestal(){return fBaseline;}
  Double_t Getv2hadron(){return fFit->GetParameter("v_{2} hadron");}
  Double_t Getv2Dmeson(){return fFit->GetParameter("v_{2} D meson");}
  
  Double_t GetNSSigmaError(){return fFit->GetParError(fFit->GetParNumber("NS #sigma"));}
  Double_t GetNSYieldError(){return fFit->GetParError(fFit->GetParNumber("NS Y"));}
  Double_t GetASSigmaError(){return fFit->GetParError(fFit->GetParNumber("AS #sigma"));}
  Double_t GetASYieldError(){return fFit->GetParError(fFit->GetParNumber("AS Y"));}
  Double_t GetPedestalError(){return fErrbaseline;}
  Double_t Getv2hadronError(){return fFit->GetParError(fFit->GetParNumber("v_{2} hadron"));}
  Double_t Getv2DmesonError(){return fFit->GetParError(fFit->GetParNumber("v_{2} D meson"));}

  //------------------Functions---------->
  TF1 * GetFitFunction(){
    if(!fFit) {printf("AliHFCorrFitter::GetFitFunction NoFitFunc - error"); return NULL;}
    return fFit;
  }

  Double_t FindBaseline();
  void Fitting(Bool_t drawSplitTerm=kTRUE);
  void YieldErrAboveBaseline();
  TH1F* SubtractBaseline();
  void DrawLegendWithParameters();
  void SetSingleTermsForDrawing(Bool_t draw);


 private:
  //---------Data Members------------->
  TH1F      *fHist;    //histogram to fit
  TF1       *fFit;
  TF1       *fGausNS;
  TF1       *fGausNS2;
  TF1       *fGausAS;
  TF1       *fGausAS2;
  TF1       *fPed;
  Double_t            fMin;         //histo min range
  Double_t            fMax;         //histo max range
  Int_t               fFixBase;      //to fix the baseline
  Int_t               fFixMean;      //to fix the mean
  Double_t            fBaseline;     //
  Double_t            fErrbaseline;
  Double_t            fNsybc;
  Double_t            fEnsybc;
  Double_t            fAsybc;
  Double_t            fEasybc;
  Double_t fMinBaselineRange; // minimum for estimation of the baseline
  Double_t fMaxBaselineRange; // maximum for estimation of the baseline
  Int_t        fTypeOfFitfunc;     //type of functions
  Int_t           fDmesonType;
  ClassDef(AliHFCorrFitter,1);
};
#endif
