#ifndef ALIHFDMESONCORRAVERAGE_H
#define ALIHFDMESONCORRAVERAGE_H
/* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: $ */

/////////////////////////////////////////////////////////////
// class to average D meson -hadron correlations
//
// Author: A. Rossi, andrea.rossi@cern.ch
/////////////////////////////////////////////////////////////
#include "AliHFDhadronCorrSystUnc.h"
#include <TMath.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TGraphAsymmErrors.h>

class AliHFDmesonCorrAverage : public TNamed {
  

 public:
  
  AliHFDmesonCorrAverage();
  AliHFDmesonCorrAverage(const char* name);
  ~AliHFDmesonCorrAverage();
  
  
  Bool_t InitSystematicUncertainty(Int_t system=-1,Int_t year=-1,Int_t centbin=0);
  void SetIncludeDzero(Bool_t inclDzero){fincludeDzero=inclDzero;}
  void SetIncludeDstar(Bool_t inclDstar){fincludeDstar=inclDstar;}
  void SetIncludeDplus(Bool_t inclDplus){fincludeDplus=inclDplus;}
  
  void SetDzeroHisto(TH1D *h){fhDzero=(TH1D*)h->Clone("hInputDzero");}
  void SetDplusHisto(TH1D *h){fhDplus=(TH1D*)h->Clone("hInputDplus");}
  void SetDstarHisto(TH1D *h){fhDstar=(TH1D*)h->Clone("hInputDstar");}
  void SetArithmeticAverage(Bool_t averType){fArithmeticAverage=averType;}  
  void SetMethod(Int_t method){fmethod=method;}
  void SetSystem(Int_t sys,Int_t year){fsys=sys;fyear=year;}
  void SetCentBin(Int_t centbin){fCentBin=centbin;}
  void SetUncertaintyFromOnlyOneMeson(Bool_t unc){fAverUncOnlyOneMeson=unc;}
  void SetMomentumRanges(Double_t minptD,Double_t maxptD,Double_t minptAsso,Double_t maxptAsso){fptminD=minptD;fptmaxD=maxptD;fptminAsso=minptAsso;fptmaxAsso=maxptAsso;}
  void CalculateAverage();  
  TH1D *GetAverageHisto(){return fhDaverage;}
  void InitAverageHisto(TH1D *h);
  TH1D *ReflectHisto(TH1D *h);
  TH1D *GetWeightsUsedDzero(){
    return fhUsedWeightsDzero;
  }
  TH1D *GetWeightsUsedDstar(){
    return fhUsedWeightsDstar;
  }
  TH1D *GetWeightsUsedDplus(){
    return fhUsedWeightsDplus;
  }
  void SetDzeroSystUnc(AliHFDhadronCorrSystUnc *sys){fSystDzero=(AliHFDhadronCorrSystUnc*)sys->Clone("fSystDzero");}
  void SetDplusSystUnc(AliHFDhadronCorrSystUnc *sys){fSystDplus=(AliHFDhadronCorrSystUnc*)sys->Clone("fSystDplus");}
  void SetDstarSystUnc(AliHFDhadronCorrSystUnc *sys){fSystDstar=(AliHFDhadronCorrSystUnc*)sys->Clone("fSystDstar");}

  void SetSystAreAlreadySet(Bool_t syst){fSystAlreadySet=syst;}
  AliHFDhadronCorrSystUnc* GetAverageSystUncertainty(){return fSystDaverage;}
 private:
  void SetWeights();


  AliHFDhadronCorrSystUnc *fSystDzero;     // Dzero syst unc
  AliHFDhadronCorrSystUnc *fSystDstar;     // Dzero syst unc
  AliHFDhadronCorrSystUnc *fSystDplus;     // Dzero syst unc
  AliHFDhadronCorrSystUnc *fSystDaverage;   // Average syst unc
  Bool_t fincludeDzero;                    //flag to include Dzero
  Bool_t fincludeDstar;                    // flag to include Dstar 
  Bool_t fincludeDplus;                    // flag to include Dplus
  Int_t fmethod;                     // flag to use abs uncertainty (first digit (unit)=0) or rel uncertainty (first digit=1) and stat only (second digit=0) or stat+uncorr syst (2nd digit =1)
  Double_t fptminD;                        // min D pt
  Double_t fptmaxD;                        // max D pt
  Double_t fptminAsso;                        // min associated track pt
  Double_t fptmaxAsso;                        // max associated pt
  TH1D *fhDzero;                           //  Dzero input histo
  TH1D *fhDstar;                           //  Dstar input histo
  TH1D *fhDplus;                           //  Dplus input histo
  TH1D *fhDaverage;                         // D average
  TGraphAsymmErrors *fgrTotSystAverage;        // 
  TGraphAsymmErrors *fgrFDSystAverage;        // 
  TGraphAsymmErrors *fgrNonFDSystAverage;        // 
  Double_t *fweightsDzeroStat;                       // Dzero weights used
  Double_t *fweightsDstarStat;                       // Dstar weights used
  Double_t *fweightsDplusStat;                       // Dplus weights used
  Double_t *fweightsDzeroSystYield;                       // Dzero weights used
  Double_t *fweightsDstarSystYield;                       // Dstar weights used
  Double_t *fweightsDplusSystYield;                       // Dplus weights used
  Double_t *fweightsDzeroSystBkg;                       // Dzero weights used
  Double_t *fweightsDstarSystBkg;                       // Dstar weights used
  Double_t *fweightsDplusSystBkg;                       // Dplus weights used
  Int_t   fnbinsphi;                        // nbins phi
  Int_t	  fsys;					//system (0=pp, 1=pPb)
  Int_t	  fyear;				// year  (2010 for pp@7 TeV, 2013 and 2016 for pPb@5.02 TeV, 2017 for pp 5 TeV, 2018 for pp 13 TeV)
  Int_t   fCentBin;                             // centrality bin
  Bool_t fSystAlreadySet;                       // Set it to kTRUE when systematic uncertainties from external files are set
  Bool_t fArithmeticAverage;                   // flag to perform arithmetic average
  Bool_t fAverUncOnlyOneMeson;                  // flag to use as uncertainty that of single meson
  TH1D *fhUsedWeightsDzero;                    // histogram with final weights used for Dzero
  TH1D *fhUsedWeightsDstar;                   // histogram with final weights used for Dstar
  TH1D *fhUsedWeightsDplus;                   // histogram with final weights used for Dplus

  ClassDef(AliHFDmesonCorrAverage,4); //class for averaging D meson -hadron correlations
};


#endif
