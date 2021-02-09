#ifndef ALIHFINVMASSFITTER_H
#define ALIHFINVMASSFITTER_H

/* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


/////////////////////////////////////////////////////////////
///
/// \class AliHFInvMassFitter
/// \brief AliHFInvMassFitter class for the fit of 
///  invariant mass distribution of charm hadrons
///
/// \author Author: F.Prino, A. Rossi, C. Bianchin
/////////////////////////////////////////////////////////////


#include <TNamed.h>
#include "AliLog.h"

class TF1;
class TH1F;

class AliHFInvMassFitter : public TNamed {
 public:

  enum ETypeOfBkg{ kExpo=0, kLin=1, kPol2=2, kNoBk=3, kPow=4, kPowEx=5};
  enum ETypeOfSgn{ kGaus=0, k2Gaus=1, k2GausSigmaRatioPar=2 };
  AliHFInvMassFitter();
  AliHFInvMassFitter(const TH1F* histoToFit, Double_t minvalue, Double_t maxvalue, Int_t fittypeb=kExpo, Int_t fittypes=kGaus);
  ~AliHFInvMassFitter();

  void     SetHistogramFit(const TH1F* histoToFit){
    if (fHistoInvMass) delete fHistoInvMass;
    fHistoInvMass=(TH1F*)histoToFit->Clone("fHistoInvMass");
    fHistoInvMass->SetDirectory(0);
  }
  void     SetRangeFit(Double_t minvalue, Double_t maxvalue){
    fMinMass=minvalue; fMaxMass=maxvalue;
  }   
  void     SetFitFunctions(Int_t fittypeb, Int_t fittypes){
    fTypeOfFit4Bkg=fittypeb; fTypeOfFit4Sgn=fittypes;
    SetNumberOfParams();
  }
  void     SetSigmaLimit(Double_t sigmaVar, Double_t sigmalimit){
    fSigmaVar=sigmaVar; fParSig=sigmalimit;
  } 

  void SetUseLikelihoodFit(){fFitOption="L,E";}
  void SetUseLikelihoodWithWeightsFit(){fFitOption="WL,E";}
  void SetUseChi2Fit(){fFitOption="E";}
  void SetFitOption(TString opt){fFitOption=opt.Data();};
  void SetParticlePdgMass(Double_t mass){fMassParticle=mass;}
  Double_t GetParticlePdgMass(){return fMassParticle;}
  void SetPolDegreeForBackgroundFit(Int_t deg){
    if(fTypeOfFit4Bkg!=6) AliFatal("fTypeOfFit4Bkg should be set to 6 to use higher order polynomials\n");
    fPolDegreeBkg=deg;
    SetNumberOfParams();
  }
  void SetInitialGaussianMean(Double_t mean) {fMass=mean;} 
  void SetInitialGaussianSigma(Double_t sigma) {fSigmaSgn=sigma;}
  void SetInitialSecondGaussianSigma(Double_t sigma) {fSigmaSgn2Gaus=sigma;}
  void SetInitialFrac2Gaus(Double_t frac) {fFrac2Gaus=frac;}
  void SetInitialRatio2GausSigma(Double_t fracsigma) {fRatio2GausSigma=fracsigma;}
  void SetFixGaussianMean(Double_t mean){
    SetInitialGaussianMean(mean); 
    fFixedMean=kTRUE;
  }
  void SetFixGaussianSigma(Double_t sigma){
    SetInitialGaussianSigma(sigma);
    fFixedSigma=kTRUE;
  }
  void SetBoundGaussianSigma(Double_t sigma, Double_t sigmalimit){
    SetInitialGaussianSigma(sigma);
    SetSigmaLimit(sigma, sigmalimit);
    fBoundSigma=kTRUE;
    fFitOption="L,E,B";
  }
  void SetFixSecondGaussianSigma(Double_t sigma){
    if(fTypeOfFit4Sgn!=k2Gaus) AliFatal("fTypeOfFit4Sgn should be set to k2Gaus to fix ratio between gaussians\n");
    SetInitialSecondGaussianSigma(sigma);
    fFixedSigma2Gaus=kTRUE;
  }
  void SetFixFrac2Gaus(Double_t frac){
    if(fTypeOfFit4Sgn!=k2Gaus && fTypeOfFit4Sgn!=k2GausSigmaRatioPar) AliFatal("fTypeOfFit4Sgn should be set to k2Gaus to fix relative integral of two gaussians\n");
    SetInitialFrac2Gaus(frac);
    fFixedFrac2Gaus=kTRUE;
  }
  void SetFixRatio2GausSigma(Double_t sigmafrac){
    if(fTypeOfFit4Sgn!=k2GausSigmaRatioPar) AliFatal("fTypeOfFit4Sgn should be set to k2GausSigmaRatioPar to fix ratio between gaussian sigmas\n");
    SetInitialRatio2GausSigma(sigmafrac);
    fFixedRatio2GausSigma=kTRUE;
  }
  void SetFixSignalYield(Double_t yield){
    fFixedRawYield=yield;
  }
  void SetNSigma4SideBands(Double_t ns=4.){
    fNSigma4SideBands=ns;
  }
  TH1F* SetTemplateReflections(const TH1 *h, TString opt,Double_t minRange,Double_t maxRange);
  void SetInitialReflOverS(Double_t rovers){fRflOverSig=rovers;}
  void     SetFixReflOverS(Double_t rovers){
    SetInitialReflOverS(rovers);
    fFixRflOverSig=kTRUE;
  }
  void SetSmoothReflectionTemplate(Bool_t opt){fSmoothRfl=opt;}

  void IncludeSecondGausPeak(Double_t mass, Bool_t fixm, Double_t width, Bool_t fixw){
    fSecondPeak=kTRUE; fSecMass=mass; fSecWidth=width;
    fFixSecMass=fixm;  fFixSecWidth=fixw;
  }
  void SetCheckSignalCountsAfterFirstFit(Bool_t opt){fCheckSignalCountsAfterFirstFit=opt;}
  
  Double_t GetRawYield()const {return fRawYield;}
  Double_t GetRawYieldError()const {return fRawYieldErr;}
  Double_t GetMean() const {return fMass;}
  Double_t GetMeanUncertainty() const {return fMassErr;}
  Double_t GetSigma()const {return fSigmaSgn;}
  Double_t GetSigmaUncertainty()const { return fSigmaSgnErr;}
  Double_t GetReflOverSig()const{
    if(fRflFunc) return fRflFunc->GetParameter(0);
    else return 0;
  }
  Double_t GetReflOverSigUncertainty()const{
    if(fRflFunc) return fRflFunc->GetParError(0);
    else return 0;
  }
  TF1*     GetBackgroundFullRangeFunc(){return fBkgFunc;}
  TF1*     GetBackgroundRecalcFunc(){return fBkgFuncRefit;}
  TF1*     GetBkgPlusReflFunc(){return fBkRFunc;}
  TF1*     GetSignalFunc(){return fSigFunc;}
  TF1*     GetMassFunc(){return fTotFunc;}
  TF1*     GetSecondPeakFunc(){return fSecFunc;}
  TF1*     GetReflFunc(){return fRflFunc;}
  Double_t GetChiSquare() const{
    if(fTotFunc) return fTotFunc->GetChisquare();
    else return -1;
  }
  Double_t GetReducedChiSquare() const{
    if(fTotFunc) return fTotFunc->GetChisquare()/fTotFunc->GetNDF();
    else return -1;
  }
  Double_t GetFitProbability() const{
    if(fTotFunc) return fTotFunc->GetProb();
    else return -1;
  }
  TH1F*    GetHistoClone() const{
    TH1F* hout=(TH1F*)fHistoInvMass->Clone(Form("%scloned",fHistoInvMass->GetName()));
    return hout;
  }
  Double_t GetRawYieldBinCounting(Double_t& errRyBC, Double_t nSigma=3., Int_t option=0, Int_t pdgCode=0) const;
  Double_t GetRawYieldBinCounting(Double_t& errRyBC, Double_t minMass, Double_t maxMass, Int_t option=0) const;
  Int_t    MassFitter(Bool_t draw=kTRUE);
  Double_t FitFunction4Sgn (Double_t* x, Double_t* par);
  Double_t FitFunction4Bkg (Double_t* x, Double_t* par);
  Double_t FitFunction4Refl(Double_t *x,Double_t *par);
  Double_t FitFunction4BkgAndRefl(Double_t *x,Double_t *par);
  Double_t FitFunction4SecPeak (Double_t* x, Double_t* par);
  Double_t FitFunction4Mass (Double_t* x, Double_t* par);
  virtual  void     Signal(Double_t nOfSigma,Double_t &signal,Double_t &errsignal) const;
  virtual  void     Signal(Double_t min,Double_t max,Double_t &signal,Double_t &errsignal) const;
  void Background(Double_t nOfSigma, Double_t &background,Double_t &errbackground) const;
  void Background(Double_t min, Double_t max, Double_t &background,Double_t &errbackground) const;
  void DrawHere(TVirtualPad* c, Double_t nsigma=3,Int_t writeFitInfo=2);
  void DrawHistoMinusFit(TVirtualPad* c,Int_t writeFitInfo=1);
  void Significance(Double_t nOfSigma, Double_t &significance,Double_t &errsignificance) const;
  void Significance(Double_t min, Double_t max, Double_t &significance,Double_t &errsignificance) const;
  // Function to compute the residuals histo-fit
  //   option=0 -> use the total fit function (signal+background+reflections)
  //   option=1 -> use only the background fit function (to isolate signal+reflections)
  //   option=2 -> use background+reflection fit function (to isolate the signal line shape)
  TH1F* GetResidualsAndPulls(TH1 *hPulls=0x0,TH1 *hResidualTrend=0x0,TH1 *hPullsTrend=0x0,Double_t minrange=0,Double_t maxrange=-1, Int_t option=0);
  // Interface method to compute the residuals histo-fit function of the background
  TH1F* GetOverBackgroundResidualsAndPulls(TH1 *hPulls=0x0,TH1 *hResidualTrend=0x0,TH1 *hPullsTrend=0x0, Double_t minrange=0,Double_t maxrange=-1);
  // Interface method to compute the residuals histo-(fit function of the background+reflections)
  TH1F* GetOverBackgroundPlusReflResidualsAndPulls(TH1 *hPulls=0x0,TH1 *hResidualTrend=0x0,TH1 *hPullsTrend=0x0, Double_t minrange=0,Double_t maxrange=-1);
  void PrintFunctions();

 private:
  AliHFInvMassFitter(const AliHFInvMassFitter &source);
  AliHFInvMassFitter& operator=(const AliHFInvMassFitter& source); 

  void  SetNumberOfParams();
  Double_t CheckForSignal(Double_t mean, Double_t sigma);
  TF1*	   CreateBackgroundFitFunction(TString fname, Double_t integral);
  TF1*	   CreateSignalFitFunction(TString fname, Double_t integral);
  TF1*	   CreateSecondPeakFunction(TString fname, Double_t integral);
  TF1*	   CreateReflectionFunction(TString fname);
  TF1*	   CreateBackgroundPlusReflectionFunction(TString fname);
  TF1*	   CreateTotalFitFunction(TString fname);
  Bool_t   PrepareHighPolFit(TF1 *fback);
  Double_t BackFitFuncPolHelper(Double_t *x,Double_t *par);

  void DrawFit();

  TH1F*     fHistoInvMass;         /// histogram to fit
  Double_t  fMinMass;              /// lower mass limit
  Double_t  fMaxMass;              /// upper mass limit
  Int_t     fTypeOfFit4Bkg;        /// background fit func
  Int_t     fPolDegreeBkg;         /// degree of polynomial expansion for back fit (option 6 for back)
  Int_t     fCurPolDegreeBkg;      /// help variable
  Double_t  fMassParticle;         /// pdg value of particle mass
  Int_t     fTypeOfFit4Sgn;        /// signal fit func
  Double_t  fMass;                 /// signal gaussian mean value
  Double_t  fMassErr;              /// unc on signal gaussian mean value
  Double_t  fSigmaSgn;             /// signal gaussian sigma
  Double_t  fSigmaSgnErr;          /// unc on signal gaussian sigma
  Double_t  fSigmaSgn2Gaus;        /// signal second gaussian sigma in case of k2Gaus
  Bool_t    fFixedMean;            /// switch for fix mean of gaussian
  Bool_t    fFixedSigma;           /// switch for fix Sigma of gaussian
  Bool_t    fBoundSigma;           /// switch for bound Sigma of gaussian
  Double_t  fSigmaVar;             /// value of bound Sigma of gaussian
  Double_t  fParSig;               /// +/- range variation of bound Sigma of gaussian in %
  Bool_t    fFixedSigma2Gaus;      /// switch for fix Sigma of second gaussian in case of k2Gaus
  Double_t  fFixedRawYield;        /// initialization for wa yield
  Double_t  fFrac2Gaus;            /// initialization for fraction of 2nd gaussian in case of k2Gaus or k2GausSigmaRatioPar
  Bool_t    fFixedFrac2Gaus;       /// switch for fixed fraction of 2nd gaussian in case of k2Gaus or k2GausSigmaRatioPar
  Double_t  fRatio2GausSigma;      /// initialization for ratio between two gaussian sigmas in case of k2GausSigmaRatioPar
  Bool_t    fFixedRatio2GausSigma; /// switch for fixed ratio between two gaussian sigmas in case of k2GausSigmaRatioPar
  Int_t     fNParsSig;             /// fit parameters in signal fit function
  Int_t     fNParsBkg;             /// fit parameters in background fit function
  Bool_t    fOnlySideBands;        /// kTRUE = only side bands considered
  Double_t  fNSigma4SideBands;     /// number of sigmas to veto the signal peak
  Bool_t    fCheckSignalCountsAfterFirstFit; /// switch for check after first fit 
  TString   fFitOption;            /// L, LW or Chi2
  Double_t  fRawYield;             /// signal gaussian integral
  Double_t  fRawYieldErr;          /// err on signal gaussian integral
  TF1*      fSigFunc;              /// Signal fit function
  TF1*      fBkgFuncSb;            /// background fit function (1st step, side bands only)
  TF1*      fBkgFunc;              /// background fit function (1st step, extended in peak region)
  TF1*      fBkgFuncRefit;         /// background fit function (2nd step)
  Bool_t    fReflections;          /// flag use/not use reflections
  Int_t     fNParsRfl;             /// fit parameters in reflection fit function
  Double_t  fRflOverSig;           /// reflection/signal
  Bool_t    fFixRflOverSig;        /// switch for fix refl/signal
  TH1F*     fHistoTemplRfl;        /// histogram with reflection template
  Bool_t    fSmoothRfl;            /// switch for smoothing of reflection template
  Double_t  fRawYieldHelp;         /// internal variable for fit with reflections
  TF1*      fRflFunc;              /// fit function for reflections
  TF1*      fBkRFunc;              /// fit function for reflections
  Bool_t    fSecondPeak;           /// switch off/on second peak (for D+->KKpi in Ds)
  Int_t     fNParsSec;             /// fit parameters in 2nd peak fit function
  Double_t  fSecMass;              /// position of the 2nd peak
  Double_t  fSecWidth;             /// width of the 2nd peak
  Bool_t    fFixSecMass;           /// flag to fix the position of the 2nd peak
  Bool_t    fFixSecWidth;          /// flag to fix the width of the 2nd peak
  TF1*      fSecFunc;              /// fit function for second peak
  TF1*      fTotFunc;              /// total fit function

  /// \cond CLASSIMP     
  ClassDef(AliHFInvMassFitter,7); /// class for invariant mass fit
  /// \endcond
};

#endif
