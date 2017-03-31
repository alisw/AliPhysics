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
  enum ETypeOfSgn{ kGaus=0, k2Gaus=1  };
  AliHFInvMassFitter();
  AliHFInvMassFitter(const TH1F* histoToFit, Double_t minvalue, Double_t maxvalue, Int_t fittypeb=kExpo, Int_t fittypes=kGaus);
  ~AliHFInvMassFitter();

  void     SetRangeFit(Double_t minvalue, Double_t maxvalue){
    fMinMass=minvalue; fMaxMass=maxvalue;
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
  void SetFixGaussianMean(Double_t mean){
    SetInitialGaussianMean(mean); 
    fFixedMean=kTRUE;
  }
  void SetFixGaussianSigma(Double_t sigma){
    SetInitialGaussianSigma(sigma); 
    fFixedSigma=kTRUE;
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
  Double_t GetRawYield()const {return fRawYield;}
  Double_t GetRawYieldError()const {return fRawYieldErr;}
  Double_t GetMean() const {return fMass;}
  Double_t GetMeanUncertainty() const {return fMassErr;}
  Double_t GetSigma()const {return fSigmaSgn;}
  Double_t GetSigmaUncertainty()const { return fSigmaSgnErr;}
  TF1*     GetBackgroundFullRangeFunc(){return fBkgFunc;}
  TF1*     GetBackgroundRecalcFunc(){return fBkgFuncRefit;}
  TF1*     GetMassFunc(){return fTotFunc;}
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
  Double_t GetRawYieldBinCounting(Double_t& errRyBC, Double_t nSigma=3., Int_t option=0, Int_t pdgCode=0) const;
  Double_t GetRawYieldBinCounting(Double_t& errRyBC, Double_t minMass, Double_t maxMass, Int_t option=0) const;

  Bool_t   MassFitter(Bool_t draw=kTRUE);
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
  void DrawHere(TVirtualPad* c, Double_t nsigma=3,Int_t writeFitInfo=1);
  void Significance(Double_t nOfSigma, Double_t &significance,Double_t &errsignificance) const;
  void Significance(Double_t min, Double_t max, Double_t &significance,Double_t &errsignificance) const;
  TH1F* GetResidualsAndPulls(TH1 *hPulls=0x0,TH1 *hResidualTrend=0x0,TH1 *hPullsTrend=0x0,Double_t minrange=0,Double_t maxrange=-1);
  void PrintFunctions();

 private:
  AliHFInvMassFitter(const AliHFInvMassFitter &source);
  AliHFInvMassFitter& operator=(const AliHFInvMassFitter& source); 

  void  SetNumberOfParams();
  Double_t  CheckForSignal(Double_t mean, Double_t sigma);
  TF1*  CreateBackgroundFitFunction(TString fname, Double_t integral);
  TF1*  CreateSignalFitFunction(TString fname, Double_t integral);
  TF1*  CreateSecondPeakFunction(TString fname, Double_t integral);
  TF1* CreateReflectionFunction(TString fname);
  TF1* CreateBackgroundPlusReflectionFunction(TString fname);
  TF1* CreateTotalFitFunction(TString fname);
  Bool_t PrepareHighPolFit(TF1 *fback);
  Double_t BackFitFuncPolHelper(Double_t *x,Double_t *par);

  void DrawFit();

  TH1F*     fHistoInvMass;     /// histogram to fit
  Double_t  fMinMass;          /// lower mass limit
  Double_t  fMaxMass;          /// upper mass limit
  Int_t     fTypeOfFit4Bkg;    /// background fit func
  Int_t     fPolDegreeBkg;     /// degree of polynomial expansion for back fit (option 6 for back)
  Int_t     fCurPolDegreeBkg;  /// help variable
  Double_t  fMassParticle;     /// pdg value of particle mass
  Int_t     fTypeOfFit4Sgn;    /// signal fit func
  Double_t  fMass;             /// signal gaussian mean value
  Double_t  fMassErr;          /// unc on signal gaussian mean value  
  Double_t  fSigmaSgn;         /// signal gaussian sigma
  Double_t  fSigmaSgnErr;      /// unc on signal gaussian sigma
  Bool_t    fFixedMean;        /// switch for fix mean of gaussian 
  Bool_t    fFixedSigma;       /// switch for fix Sigma of gaussian 
  Double_t  fFixedRawYield;    /// initialization for wa yield
  Int_t     fNParsSig;         /// fit parameters in signal fit function
  Int_t     fNParsBkg;         /// fit parameters in background fit function
  Bool_t    fOnlySideBands;    /// kTRUE = only side bands considered
  Double_t  fNSigma4SideBands; /// number of sigmas to veto the signal peak
  TString   fFitOption;        /// L, LW or Chi2
  Double_t  fRawYield;         /// signal gaussian integral
  Double_t  fRawYieldErr;      /// err on signal gaussian integral
  TF1*      fSigFunc;          /// Signal fit function 
  TF1*      fBkgFuncSb;        /// background fit function (1st step, side bands only)
  TF1*      fBkgFunc;          /// background fit function (1st step, extended in peak region) 
  TF1*      fBkgFuncRefit;     /// background fit function (2nd step)
  Bool_t    fReflections;      /// flag use/not use reflections
  Int_t     fNParsRfl;         /// fit parameters in reflection fit function
  Double_t  fRflOverSig;       /// reflection/signal
  Bool_t    fFixRflOverSig;    /// switch for fix refl/signal
  TH1F*     fHistoTemplRfl;    /// histogram with reflection template
  Bool_t    fSmoothRfl;        /// switch for smoothing of reflection template
  Double_t  fRawYieldHelp;     /// internal variable for fit with reflections
  TF1*      fRflFunc;          /// fit function for reflections
  TF1*      fBkRFunc;          /// fit function for reflections
  Bool_t    fSecondPeak;       /// switch off/on second peak (for D+->KKpi in Ds)
  Int_t     fNParsSec          ;/// fit parameters in 2nd peak fit function
  Double_t  fSecMass;          /// position of the 2nd peak
  Double_t  fSecWidth;         /// width of the 2nd peak
  Bool_t    fFixSecMass;       /// flag to fix the position of the 2nd peak
  Bool_t    fFixSecWidth;      /// flag to fix the width of the 2nd peak
  TF1*      fSecFunc;          /// fit function for second peak
  TF1*      fTotFunc;          /// total fit function

  /// \cond CLASSIMP     
  ClassDef(AliHFInvMassFitter,3); /// class for invariant mass fit
  /// \endcond
};

#endif
