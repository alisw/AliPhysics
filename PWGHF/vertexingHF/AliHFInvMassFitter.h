#ifndef ALIHFINVMASSFITTER_H
#define ALIHFINVMASSFITTER_H

#include <TNamed.h>

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
  TF1*     GetBackgroundRecalcFunc(){return fBkgFuncRef;}
  TF1*     GetMassFunc(){return fFuncTot;}
  Double_t GetChiSquare() const{
    if(fFuncTot) return fFuncTot->GetChisquare();
    else return -1;
  }
  Double_t GetReducedChiSquare() const{
    if(fFuncTot) return fFuncTot->GetChisquare()/fFuncTot->GetNDF();
    else return -1;
  }
  Double_t GetFitProbability() const{
    if(fFuncTot) return fFuncTot->GetProb();  
    else return -1;
  }
    
  Bool_t   MassFitter(Bool_t draw=kTRUE);
  Double_t FitFunction4Sgn (Double_t* x, Double_t* par);
  Double_t FitFunction4Bkg (Double_t* x, Double_t* par);
  Double_t FitFunction4SecPeak (Double_t* x, Double_t* par);
  Double_t FitFunction4Mass (Double_t* x, Double_t* par);
  virtual  void     Signal(Double_t nOfSigma,Double_t &signal,Double_t &errsignal) const;
  virtual  void     Signal(Double_t min,Double_t max,Double_t &signal,Double_t &errsignal) const;
  void Background(Double_t nOfSigma, Double_t &background,Double_t &errbackground) const;
  void Background(Double_t min, Double_t max, Double_t &background,Double_t &errbackground) const;
  void DrawHere(TVirtualPad* c);
  void Significance(Double_t nOfSigma, Double_t &significance,Double_t &errsignificance) const;
  void Significance(Double_t min, Double_t max, Double_t &significance,Double_t &errsignificance) const;

 private:
  AliHFInvMassFitter(const AliHFInvMassFitter &source);
  AliHFInvMassFitter& operator=(const AliHFInvMassFitter& source); 

  void  SetNumberOfParams();
  Double_t  CheckForSignal(Double_t mean, Double_t sigma);
  TF1*  CreateBackgroundFitFunction(TString fname, Double_t integral);
  TF1*  CreateSignalFitFunction(TString fname, Double_t integral);
  TF1*  CreateSecondPeakFunction(TString fname, Double_t integral);
  TF1* CreateTotalFitFunction(TString fname);

  void DrawFit();

  TH1F*     fHistoInvMass;     /// histogram to fit
  Double_t  fMinMass;          /// lower mass limit
  Double_t  fMaxMass;          /// upper mass limit
  Int_t     fTypeOfFit4Bkg;    /// background fit func
  Int_t     fTypeOfFit4Sgn;    /// signal fit func
  Double_t  fMass;             /// signal gaussian mean value
  Double_t  fMassErr;          /// unc on signal gaussian mean value  
  Double_t  fSigmaSgn;         /// signal gaussian sigma
  Double_t  fSigmaSgnErr;      /// unc on signal gaussian sigma
  Bool_t    fFixedMean;        /// switch for fix mean of gaussian 
  Bool_t    fFixedSigma;       /// switch for fix Sigma of gaussian 
  Double_t  fFixedRawYield;    /// initialization for wa yield
  Int_t     fNSigPars;         /// fit parameters in signal fit function
  Int_t     fNBkgPars;         /// fit parameters in background fit function
  Bool_t    fOnlySideBands;    /// kTRUE = only side bands considered
  TString   fFitOption;        /// L, LW or Chi2
  Double_t  fRawYield;         /// signal gaussian integral
  Double_t  fRawYieldErr;      /// err on signal gaussian integral
  TF1*      fSigFunc;          /// Signal fit function 
  TF1*      fBkgFuncSb;        /// background fit function (1st step)
  TF1*      fBkgFunc;          /// background fit function (1st step)
  TF1*      fBkgFuncRef;       /// background fit function (2nd step)
  Bool_t fSecondPeak;          /// swicth off/on second peak (for D+->KKpi in Ds)
  Double_t fSecMass;           /// position of the 2nd peak
  Double_t fSecWidth;          /// width of the 2nd peak
  Bool_t fFixSecMass;          /// flag to fix the position of the 2nd peak
  Bool_t fFixSecWidth;         /// flag to fix the width of the 2nd peak
  TF1*      fSecFunc;          /// fit function for second peak
  TF1*      fFuncTot;          /// total fit function 

  /// \cond CLASSIMP     
  ClassDef(AliHFInvMassFitter,1); /// class for invariant mass fit
  /// \endcond
};

#endif
