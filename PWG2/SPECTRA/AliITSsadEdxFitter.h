#ifndef ALIITSSADEDXFITTER_H
#define ALIITSSADEDXFITTER_H
/* Copyright(c) 2007-2011, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: $ */

////////////////////////////////////////////////////
//class to perform different gaussian fits to the //
//dEdx distribution, using different approach,    //
//for the ITS stand-alone track spectra analysis  //
//E. Biolcati, F. Prino                           //
////////////////////////////////////////////////////

#include <TObject.h>
class TGraph;

class AliITSsadEdxFitter  : public TObject {

 public:
  AliITSsadEdxFitter();  
  virtual ~AliITSsadEdxFitter(){};

  static Double_t CalcSigma(Int_t code,Float_t x, Bool_t mc);
  static Int_t CalcMean(Int_t code,Float_t x, Float_t mean0, Float_t &mean1, Float_t &mean2);

  void GetFitPar(Double_t *fitpar, Double_t *fitparerr) const;
  void DoFitTail(TH1F *h, Int_t bin, Int_t code);
  void DoFit(TH1F *h, Int_t bin, Int_t code, Bool_t mc, TGraph *gres);
  void DoFitProton(TH1F *h, Int_t bin, Int_t code, Bool_t mc, TGraph *gres);
  void FillHisto(TH1F *hsps, Int_t bin, Float_t binsize, Int_t code);
  void FillHistoMC(TH1F *hsps, Int_t bin, Int_t code, TH1F *h);
  Bool_t IsGoodBin(Int_t bin,Int_t code);

  void SetRangeStep1(Double_t dxlow=-0.2, Double_t dxup=0.3){
    fRangeStep1[0]=dxlow;
    fRangeStep1[1]=dxup;
  }
  void SetRangeStep2(Double_t dxlow=-0.1, Double_t dxup=0.3){
    fRangeStep2[0]=dxlow;
    fRangeStep2[1]=dxup;
  }
  void SetRangeStep3(Double_t dxlow=-0.1, Double_t dxup=2.5){
    fRangeStep3[0]=dxlow;
    fRangeStep3[1]=dxup;
  }
  void SetRangeFinalStep(Double_t dxlow=-3.5, Double_t dxup=3.5){
    fRangeStep4[0]=dxlow;
    fRangeStep4[1]=dxup;
  }
  void SetLimitsOnSigmaPion(Double_t smin=0.98, Double_t smax=1.02){
    fLimitsOnSigmaPion[0]=smin;
    fLimitsOnSigmaPion[1]=smax;
  }
  void SetLimitsOnSigmaKaon(Double_t smin=0.98, Double_t smax=1.02){
    fLimitsOnSigmaKaon[0]=smin;
    fLimitsOnSigmaKaon[1]=smax;
  }
  void SetLimitsOnSigmaProton(Double_t smin=0.98, Double_t smax=1.02){
    fLimitsOnSigmaProton[0]=smin;
    fLimitsOnSigmaProton[1]=smax;
  }

  void PrintAll() const;
  void CalcResidual(TH1F *h,TF1 *fun,TGraph *gres) const;
  Double_t GausPlusTail(const Double_t x, const Double_t mean, Double_t rms, Double_t c, Double_t slope, Double_t cut ) const;  
  Double_t GausOnBackground(const Double_t* x, const Double_t *par) const;
  void DrawFitFunction(TF1 *fun) const;

 private:
  Double_t fFitpar[5];     // array with fit parameters
  Double_t fFitparErr[5];  // array with fit parameter errors 
  Double_t fRangeStep1[2]; // Range for Step1 (w.r.t pion peak)
  Double_t fRangeStep2[2]; // Range for Step2 (w.r.t kaon/proton peak)
  Double_t fRangeStep3[2]; // Range for Step3 (w.r.t proton/kaon peak)
  Double_t fRangeStep4[2]; // Range for Last Fit
  Double_t fLimitsOnSigmaPion[2]; // limits on sigma pions
  Double_t fLimitsOnSigmaKaon[2]; // limits on sigma pions
  Double_t fLimitsOnSigmaProton[2]; // limits on sigma protons

  ClassDef(AliITSsadEdxFitter,1);
};

#endif

