#ifndef ALIITSSADEDXFITTER_H
#define ALIITSSADEDXFITTER_H
/* Copyright(c) 2007-2011, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////
// Class with the fits algorithms to be used in the identified       //
// spectra analysis using the ITS in stand-alone mode                //
// Author: E.Biolcati, biolcati@to.infn.it                           //
//         F.Prino, prino@to.infn.it                                 //
///////////////////////////////////////////////////////////////////////

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
	void DoFitOnePeak(TH1F *h, Int_t bin, Int_t signedcode, Bool_t mc);
	void DoFitProtonFirst(TH1F *h, Int_t bin, Int_t signedcode, Bool_t mc, TGraph *gres);
	void GetInitialParam(TH1F* h,Bool_t mc,Int_t code,Int_t bin, Float_t &pt, Float_t &ampl, Float_t &mean1, Float_t &mean2, Float_t &mean3, Float_t &sigma1, Float_t &sigma2, Float_t &sigma3);
  void FillHisto(TH1F *hsps, Int_t bin, Float_t binsize, Int_t code);
  void FillHistoMC(TH1F *hsps, Int_t bin, Int_t code, TH1F *h);
  Bool_t IsGoodBin(Int_t bin,Int_t code);

  void SetRangeStep1(Double_t dxlow=-0.2, Double_t dxup=0.3){
    fRangeStep1[0]=dxlow;
    fRangeStep1[1]=dxup;
  }
  void SetRangeStep2(Double_t dxlow=-0.1, Double_t dxup=0.4){
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
  void SetLimitsOnSigmaPion(Double_t smin=0.95, Double_t smax=1.05){
    fLimitsOnSigmaPion[0]=smin;
    fLimitsOnSigmaPion[1]=smax;
  }
  void SetLimitsOnSigmaKaon(Double_t smin=0.95, Double_t smax=1.05){
    fLimitsOnSigmaKaon[0]=smin;
    fLimitsOnSigmaKaon[1]=smax;
  }
  void SetLimitsOnSigmaProton(Double_t smin=0.95, Double_t smax=1.05){
    fLimitsOnSigmaProton[0]=smin;
    fLimitsOnSigmaProton[1]=smax;
  }
	void SetBinsUsedPion(Int_t bmin=2, Int_t bmax=14){
		fBinsUsedPion[0]=bmin;
		fBinsUsedPion[1]=bmax;
	}
	void SetBinsUsedKaon(Int_t bmin=7, Int_t bmax=12){
		fBinsUsedKaon[0]=bmin;
		fBinsUsedKaon[1]=bmax;
	}
	void SetBinsUsedProton(Int_t bmin=8, Int_t bmax=17){
		fBinsUsedProton[0]=bmin;
		fBinsUsedProton[1]=bmax;
	}

  void PrintAll() const;
  void CalcResidual(TH1F *h,TF1 *fun,TGraph *gres) const;
  Double_t GausPlusTail(const Double_t x, const Double_t mean, Double_t rms, Double_t c, Double_t slope, Double_t cut ) const;  
  Double_t GausOnBackground(const Double_t* x, const Double_t *par) const;
  void DrawFitFunction(TF1 *fun) const;

 private:
  Double_t fFitpar[3];     // array with fit parameters
  Double_t fFitparErr[3];  // array with fit parameter errors 
  Double_t fRangeStep1[2]; // Range for Step1 (w.r.t pion peak)
  Double_t fRangeStep2[2]; // Range for Step2 (w.r.t kaon/proton peak)
  Double_t fRangeStep3[2]; // Range for Step3 (w.r.t proton/kaon peak)
  Double_t fRangeStep4[2]; // Range for Last Fit
  Double_t fLimitsOnSigmaPion[2]; // limits on sigma pions
  Double_t fLimitsOnSigmaKaon[2]; // limits on sigma kaons
  Double_t fLimitsOnSigmaProton[2]; // limits on sigma protons
  Int_t fBinsUsedPion[2];   // limits on bins used pions
  Int_t fBinsUsedKaon[2];   // limits on bins used kaons
  Int_t fBinsUsedProton[2]; // limits on bins used protons

  ClassDef(AliITSsadEdxFitter,2);
};

#endif

