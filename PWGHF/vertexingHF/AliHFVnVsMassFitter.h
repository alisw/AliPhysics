#ifndef ALIHFVNVSMASSFITTER_H
#define ALIHFVNVSMASSFITTER_H
  /// \class AliHfVnVsMassFitter
  /// \class that performs the vn vs mass simultaneus fit for D mesons
  /// \author Fabrizio Grosa <grosa@to.infn.it>, INFN Torino

#include <TObject.h>
#include <Riostream.h>
#include <TVirtualPad.h>
#include <TH1F.h>
#include "Fit/Fitter.h"
#include "Fit/Chi2FCN.h"
#include "Math/WrappedMultiTF1.h"
#include "AliHFInvMassFitter.h"

class AliHFVnVsMassFitter : public TObject {

public:
  AliHFVnVsMassFitter();
  AliHFVnVsMassFitter(TH1F* hMass, TH1F* hvn, Double_t min, Double_t max, Int_t funcMassBkg, Int_t funcMassSgn, Int_t funcvnBkg);
  ~AliHFVnVsMassFitter();

  enum ETypeOfBkg{kExpo=0, kLin=1, kPol2=2, kNoBk=3, kPow=4, kPowEx=5, kPoln=6};
  enum ETypeOfSgn{kGaus=0, k2Gaus=1};

  Bool_t SimultaneusFit(Bool_t drawFit=kTRUE);
  void DrawHere(TVirtualPad* c);

  //setters
  void SetInitialGaussianSigma(Double_t sigma) {fSigmaInit=sigma; fSigmaFixed=1;}
  void SetInitialGaussianMean(Double_t mean) {fMeanInit=mean; fMeanFixed=1;}
  void SetParticlePdgMass(Double_t mass){fMassParticle=mass;}
  void SetMassSgnFunc(Int_t functype) {fMassSgnFuncType=functype;}
  void SetMassBkgFunc(Int_t functype) {fMassBkgFuncType=functype;}
  void SetVnBkgFunc(Int_t functype) {fVnBkgFuncType=functype;}
  void FixSigmaFromMassFit() {fSigmaFixedFromMassFit=kTRUE;}
  void FixMeanFromMassFit() {fMeanFixedFromMassFit=kTRUE;}
  void SetNSigmaForVnSB(Int_t nsigma=4) {fNSigmaForSB=nsigma;}
  void SetPolDegreeForBackgroundFit(Int_t deg){
    if(fMassBkgFuncType!=6) AliFatal("fMassBkgFuncType should be set to 6 to use higher order polynomials\n");
    fPolDegreeBkg=deg;
  }
  void SetTemplateReflections(const TH1 *h, TString opt, Double_t minRange, Double_t maxRange) {
    fHistoTemplRflInit=(TH1F*)h->Clone();
    /// option could be:
    ///    "template"                use MC histograms
    ///    "1gaus" ot "singlegaus"   single gaussian function fit to MC templates
    ///    "2gaus" ot "doublegaus"   double gaussian function fit to MC templates
    ///    "pol3"                    3rd order polynomial fit to MC templates
    ///    "pol6"                    6th order polynomial fit to MC templates
    fRflOpt=opt;
    fMinRefl=minRange;
    fMaxRefl=maxRange;
    fReflections=kTRUE;
  }
  void SetInitialReflOverS(Double_t rovers){fRflOverSig=rovers;}
  void SetFixReflOverS(Double_t rovers){
    SetInitialReflOverS(rovers);
    fFixRflOverSig=kTRUE;
  }
  void IncludeSecondGausPeak(Double_t mass, Bool_t fixm, Double_t width, Bool_t fixw, Bool_t doVn){
    fSecondPeak=kTRUE; fSecMass=mass; fSecWidth=width;
    fFixSecMass=fixm;  fFixSecWidth=fixw;
    fDoSecondPeakVn=doVn;
  }
  void SetHarmonic(Int_t harmonic=2) {fHarmonic=harmonic;}

  //getters
  Double_t GetVn() const {return fVn;}
  Double_t GetVnUncertainty() const {return fVnUncertainty;}
  Double_t GetMean() const {return fMean;}
  Double_t GetMeanUncertainty() const {return fMeanUncertainty;}
  Double_t GetSigma() const {return fSigma;}
  Double_t GetSigmaUncertainty() const {return fSigmaUncertainty;}
  Double_t GetRawYield() const {return fRawYield;}
  Double_t GetRawYieldUncertainty() const {return fRawYieldUncertainty;}
  Double_t GetChiSquare() const {return fChiSquare;}
  Int_t GetNDF() const {return fNDF;}
  Double_t GetReducedChiSquare() const {return fChiSquare/fNDF;}
  Double_t GetFitProbability() const {return fProb;}
  Double_t GetParticlePdgMass() const {return fMassParticle;}
  TH1F* GetTemplateReflections() {
    if(fHistoTemplRfl) {return (TH1F*)fHistoTemplRfl->Clone();}
    else if(fHistoTemplRflInit) {return (TH1F*)fHistoTemplRflInit->Clone();}
    else {return 0;}
  }
  void Signal(Double_t nOfSigma,Double_t &signal,Double_t &errsignal) const;
  void Signal(Double_t min,Double_t max,Double_t &signal,Double_t &errsignal) const;
  void Background(Double_t nOfSigma, Double_t &background,Double_t &errbackground) const;
  void Background(Double_t min, Double_t max, Double_t &background,Double_t &errbackground) const;
  void Significance(Double_t nOfSigma, Double_t &significance,Double_t &errsignificance) const;
  void Significance(Double_t min, Double_t max, Double_t &significance,Double_t &errsignificance) const;
  TF1* GetMassTotFitFunc() {
    if(fMassTotFunc) return fMassTotFunc;
    else return 0x0;
  }
  TF1* GetVnVsMassTotFitFunc() {
    if(fVnTotFunc) return fVnTotFunc;
    else return 0x0;
  }

  //struct for global chi2 (for simultaneus fit)
  struct AliHFGlobalChi2 {
    AliHFGlobalChi2(ROOT::Math::IMultiGenFunction & f1,ROOT::Math::IMultiGenFunction & f2) : fChi2_1(&f1), fChi2_2(&f2) {}

    double operator() (const double *par) const {
        return (*fChi2_1)(par) + (*fChi2_2)(par);
    }
    const  ROOT::Math::IMultiGenFunction * fChi2_1;
    const  ROOT::Math::IMultiGenFunction * fChi2_2;
  };

private:

    ///fit functions
  Double_t GetGausPDF(Double_t x, Double_t mean, Double_t sigma);
  Double_t GetExpoPDF(Double_t x, Double_t slope);
  Double_t GetPolPDF(Double_t x, Double_t *pars, Int_t order, Bool_t isnorm=kTRUE);
  Double_t GetPowerFuncPDF(Double_t x, Double_t *pars);
  Double_t GetPowerExpoPDF(Double_t x, Double_t *pars);
  Double_t GetHigherPolFuncPDF(Double_t x, Double_t *pars);
  Double_t MassSignal(Double_t *m, Double_t *pars);
  Double_t MassBkg(Double_t *m, Double_t *pars);
  Double_t MassRfl(Double_t *m,Double_t *par);
  Double_t MassBkgRfl(Double_t *m,Double_t *par);
  Double_t MassSecondPeak(Double_t *m,Double_t *par);
  Double_t vnBkgFunc(Double_t *m, Double_t *pars);
  Double_t MassFunc(Double_t *m, Double_t *pars);
  Double_t vnFunc(Double_t *m, Double_t *pars);

    ///private methods
  void DefineNumberOfParameters();
  Bool_t MassPrefit();
  Bool_t VnSBPrefit();
  void DrawFit();
  void SetParNames();

    ///data members
  TH1F*               fMassHisto;             /// mass histogram to fit
  TH1F*               fVnVsMassHisto;         /// vn vs. mass histogram to fit
  Int_t               fMassSgnFuncType;       /// type of mass signal fit function
  Int_t               fMassBkgFuncType;       /// type of mass bkg fit function
  Int_t               fVnBkgFuncType;         /// type of vn bkg fit function
  TF1*                fMassFuncFromPrefit;    /// mass fit function (1st step, from prefit)
  TF1*                fMassBkgFunc;           /// mass bkg fit function (final, after simultaneus fit)
  TF1*                fMassSgnFunc;           /// mass signal fit function (final, after simultaneus fit)
  TF1*                fMassTotFunc;           /// mass fit function (final, after simultaneus fit)
  TF1*                fVnBkgFuncSb;           /// vn bkg fit function (1st step from SB prefit)
  TF1*                fVnBkgFunc;             /// vn bkg fit function (final, after simultaneus fit)
  TF1*                fVnTotFunc;             /// vn fit function (final, after simultaneus fit)
  AliHFInvMassFitter* fMassFitter;            /// mass fitter for mass prefit
  Double_t            fMassMin;               /// upper mass limit
  Double_t            fMassMax;               /// lower mass limit
  Double_t            fVn;                    /// vn of the signal from fit
  Double_t            fVnUncertainty;         /// uncertainty on vn of the signal from simultaneus fit
  Double_t            fSigma;                 /// mass peak width from simultaneus fit
  Double_t            fSigmaUncertainty;      /// uncertainty on mass peak width from simultaneus fit
  Double_t            fMean;                  /// mass peak position from simultaneus fit
  Double_t            fMeanUncertainty;       /// uncertainty on mass peak position from simultaneus fit
  Double_t            fRawYield;              /// raw yield from simultaneus fit
  Double_t            fRawYieldUncertainty;   /// uncertainty raw yield from simultaneus fit
  Double_t            fChiSquare;             /// simultaneus fit chi square
  Int_t               fNDF;                   /// simultaneus fit number of degree of freedom
  Double_t            fProb;                  /// simultaneus fit probability
  Int_t               fNSigmaForSB;           /// number of sigma for sidebands region (vn bkg prefit)
  Double_t            fSigmaInit;             /// initialization for peak width
  Double_t            fMeanInit;              /// initialization for peak position
  Bool_t              fSigmaFixedFromMassFit; /// flag to fix peak width from mass prefit
  Bool_t              fMeanFixedFromMassFit;  /// flag to fix peak position from mass prefit
  Double_t            fMassParticle;          /// mass of selected particle
  Int_t               fNParsMassSgn;          /// number of parameters in mass signal fit function
  Int_t               fNParsMassBkg;          /// number of parameters in mass bkg fit function
  Int_t               fNParsVnBkg;            /// number of parameters in vn bkg fit function
  Int_t               fSigmaFixed;            /// flag to fix peak width
  Int_t               fMeanFixed;             /// flag to fix peak position
  Int_t               fPolDegreeBkg;          /// degree of polynomial expansion for back fit (option 6 for back)
  Bool_t              fReflections;           /// flag use/not use reflections
  Int_t               fNParsRfl;              /// fit parameters in reflection fit function
  Double_t            fRflOverSig;            /// reflection/signal
  Bool_t              fFixRflOverSig;         /// switch for fix refl/signal
  TH1F*               fHistoTemplRfl;         /// histogram with reflection template
  TH1F*               fHistoTemplRflInit;     /// initial histogram with reflection template
  TF1*                fMassRflFunc;               /// fit function for reflections
  TF1*                fMassBkgRflFunc;       /// mass bkg fit function plus reflections (final, after simultaneus fit)
  TString             fRflOpt;                /// refelction option
  Double_t            fMinRefl;               /// minimum for refelction histo
  Double_t            fMaxRefl;               /// maximum for refelction histo
  Bool_t              fSmoothRfl;             /// switch for smoothing of reflection template
  Double_t            fRawYieldHelp;          /// internal variable for fit with reflections
  Bool_t              fSecondPeak;            /// switch off/on second peak (for D+->KKpi in Ds)
  TF1*                fMassSecPeakFunc;           /// fit function for second peak
  Int_t               fNParsSec;              /// number of parameters in second peak fit function
  Double_t            fSecMass;               /// position of the 2nd peak
  Double_t            fSecWidth;              /// width of the 2nd peak
  Bool_t              fFixSecMass;            /// flag to fix the position of the 2nd peak
  Bool_t              fFixSecWidth;           /// flag to fix the width of the 2nd peak
  Double_t            fVnSecPeak;             /// vn of second peak from fit
  Bool_t              fDoSecondPeakVn;        /// flag to introduce second peak vn in the vn vs. mass fit
  Double_t            fVnSecPeakUncertainty;  /// vn uncertainty of second peak from fit
  Int_t               fHarmonic;              /// harmonic number for drawing

    /// \cond CLASSDEF
  ClassDef(AliHFVnVsMassFitter,1);
    /// \endcond
};
#endif //ALIHFVNVSMASSFITTER
