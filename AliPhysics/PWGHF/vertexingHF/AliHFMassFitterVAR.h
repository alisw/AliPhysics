#ifndef ALIHFMASSFITTERVAR_H
#define ALIHFMASSFITTERVAR_H
/* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */ 

/////////////////////////////////////////////////////////////
///
/// \class AliHFMassFitterVAR
/// \brief AliHFMassFitterVAR for the fit of invariant mass distribution
/// of charmed mesons
///
/// \author Author: C.Bianchin, chiara.bianchin@pd.infn.it
/////////////////////////////////////////////////////////////

#include <TNamed.h>
#include <TString.h>
#include "AliHFMassFitter.h"
class TF1;
class TNtuple;
class TFile;
class TList;
class TH1F;
class TVirtualPad;
class TPaveText;

class AliHFMassFitterVAR : public AliHFMassFitter {

 public:
  AliHFMassFitterVAR();
  AliHFMassFitterVAR(const TH1F* histoToFit, Double_t minvalue, Double_t maxvalue, Int_t rebin=1,Int_t fittypeb=0,Int_t fittypes=0);
  virtual ~AliHFMassFitterVAR();

  AliHFMassFitterVAR(const AliHFMassFitterVAR &mfit);
  AliHFMassFitterVAR& operator=(const AliHFMassFitterVAR &mfit);

  /// setters
/*   void     SetHisto(const TH1F *histoToFit); */
/*   void     SetRangeFit(Double_t minvalue, Double_t maxvalue){fminMass=minvalue; fmaxMass=maxvalue; CheckRangeFit();} */
/*   void     SetMinRangeFit(Double_t minvalue){fminMass=minvalue;printf("CheckRangeFit after SetMaxRangeFit is also set\n");} */
/*   void     SetMaxRangeFit(Double_t maxvalue){fmaxMass=maxvalue;printf("CheckRangeFit after SetMinRangeFit is also set\n");} */
/*   void     SetBinN(Int_t newbinN){fNbin=newbinN;} */
/*   void     SetType(Int_t fittypeb, Int_t fittypes); */
/*   void     SetReflectionSigmaFactor(Int_t constant) {ffactor=constant;} */
/*   void     SetInitialGaussianMean(Double_t mean) {fMass=mean;} // change the default value of the mean */
/*   void     SetInitialGaussianSigma(Double_t sigma) {fSigmaSgn=sigma;} // change the default value of the sigma */
/*   void     SetSideBands(Bool_t onlysidebands=kTRUE) {fSideBands=onlysidebands;} // consider only side bands */
/*   void     SetFixParam(Bool_t *fixpar){fFixPar=fixpar;} */
/*   void     SetDefaultFixParam(); */
/*   Bool_t   SetFixThisParam(Int_t thispar,Bool_t fixpar); */

/*   //getters */
/*   TH1F*    GetHistoClone() const; //return the histogram */
/*   void     GetRangeFit(Double_t &minvalue, Double_t &maxvalue) const {minvalue=fminMass; maxvalue=fmaxMass;} */
/*   Double_t GetMinRangeFit()const {return fminMass;} */
/*   Double_t GetMaxRangeFit()const {return fmaxMass;} */
/*   Int_t    GetBinN()       const {return fNbin;} */
/*   void     GetFitPars(Float_t* pars) const; */
/*   Int_t    GetNFinalPars() const {return fNFinalPars;} */
/*   void     GetTypeOfFit(Bool_t &background, Int_t &typeb) const {background = fWithBkg; typeb = ftypeOfFit4Bkg;} */
/*   Int_t    GetReflectionSigmaFactor() const {return ffactor;}  */
/*   Double_t GetMean() const {return fMass;} */
/*   Double_t GetMeanUncertainty() const {return fMassErr;} */
/*   Double_t GetSigma()const {return fSigmaSgn;} */
/*   Double_t GetSigmaUncertainty()const { return fSigmaSgnErr;} */
/*   Double_t GetRawYield()const {return fRawYield;} */
/*   Double_t GetRawYieldError()const {return fRawYieldErr;} */
/*   Double_t GetChiSquare() const; */
/*   Double_t GetReducedChiSquare() const; */
/*   void     GetSideBandsBounds(Int_t& lb, Int_t& hb) const; */
/*   Bool_t*  GetFixParam()const {return fFixPar;} */
/*   Bool_t   GetFixThisParam(Int_t thispar)const; */
/*   TVirtualPad* GetPad(Double_t nsigma=3,Int_t writeFitInfo=1)const; */

/*   void     PrintParTitles() const; */

/*   void     InitNtuParam(TString ntuname="ntupar"); // initialize TNtuple to store the parameters */
/*   void     FillNtuParam(); //Fill the TNtuple with the current parameters */
/*   TNtuple* GetNtuParam() const {return fntuParam;} // return the TNtuple */
/*   TNtuple* NtuParamOneShot(TString ntuname="ntupar"); // the three functions above all together */
/*   void     WriteHisto(TString path="./") const; // write the histogram */
/*   void     WriteNtuple(TString path="./") const; // write the TNtuple */
  void     WriteCanvas(TString userIDstring="",TString path="./",Double_t nsigma=3,Int_t writeFitInfo=1,Bool_t draw=kFALSE) const; //write the canvas in a root file
  void     DrawHere(TVirtualPad* pd,Double_t nsigma=3,Int_t writeFitInfo=1); 
/*   void     DrawFit(Double_t nsigma=3) const; */
/*   void     Reset(); */
  
  void     IntS(Float_t *valuewitherror) const;    /// integral of signal given my the fit with error
  Double_t IntTot() const {return fhistoInvMass->Integral("width");}  /// return total integral of the histogram
  void     Signal(Double_t nOfSigma,Double_t &signal,Double_t &errsignal) const; /// signal in nsigma with error
  void     Signal(Double_t min,Double_t max,Double_t &signal,Double_t &errsignal) const; /// signal in (min, max) with error
  void     Background(Double_t nOfSigma,Double_t &background,Double_t &errbackground) const; /// backgournd in nsigma with error
  void     Background(Double_t min,Double_t max,Double_t &background,Double_t &errbackground) const; /// backgournd in (min, max) with error
  void     Significance(Double_t nOfSigma,Double_t &significance,Double_t &errsignificance) const; /// significance in nsigma with error
  void     Significance(Double_t min,Double_t max,Double_t &significance,Double_t &errsignificance) const; /// significance in (min, max) with error
  Double_t GetReflOverSignal(Double_t &err)const;
  TH1F*    SetTemplateReflections(const TH1F *h,TString option="templ",Double_t minRange=1.72,Double_t maxRange=2.05);
  void     SetFixReflOverS(Double_t val,Bool_t fixpar=kTRUE){fReflInit=val;
    if(ftypeOfFit4Sgn!=2)Printf("Set or Fixing Refl/Signal but fit with template not foreseen");
    fFixParReflExternalValue[0]=fixpar;
    fparReflFixExt[0]=val;
  }
  void     SetFixGaussianMean(Double_t mean=1.865,Bool_t fixpar=kTRUE){fFixParSignExternalValue[1]=fixpar;     fparSignFixExt[1]=mean;} 
  void     SetFixGaussianSigma(Double_t sigma=0.012, Bool_t fixpar=kTRUE){fFixParSignExternalValue[2]=fixpar;     fparSignFixExt[2]=sigma;} 
  void SetBackHighPolDegree(Int_t deg);
  Double_t BackFitFuncPolHelper(Double_t *x,Double_t *par);
  Bool_t PrepareHighPolFit(TF1 *fback);
  void SetParticlePdgMass(Double_t mass){fMassParticle=mass;}
  Double_t GetParticlePdgMass(){return fMassParticle;}
  Double_t FitFunction4MassDistr (Double_t* x, Double_t* par);
  Double_t FitFunction4Sgn (Double_t* x, Double_t* par);
  Double_t FitFunction4Bkg (Double_t* x, Double_t* par);
  Double_t FitFunction4Refl (Double_t* x, Double_t* par);
  Double_t FitFunction4BkgAndReflDraw(Double_t *x, Double_t *par);
  Bool_t   MassFitter(Bool_t draw=kTRUE);
  Bool_t   RefitWithBkgOnly(Bool_t draw=kTRUE);
  TH1F* GetOverBackgroundResidualsAndPulls(Double_t minrange=0,Double_t maxrange=-1,TH1 *hPulls=0x0,TH1 *hResidualTrend=0x0,TH1 *hPullsTrend=0x0);
  TH1F* GetAllRangeResidualsAndPulls(Double_t minrange=0,Double_t maxrange=-1,TH1 *hPulls=0x0,TH1 *hResidualTrend=0x0,TH1 *hPullsTrend=0x0);
  TPaveText* GetYieldBox(Double_t nsigma=3.);
  TPaveText* GetFitParametersBox(Double_t nsigma=3.,Int_t mode=0);

/*   void     RebinMass(Int_t bingroup=1); */
/*   TF1*     GetBackgroundFullRangeFunc(){ */
/*     return fhistoInvMass->GetFunction("funcbkgFullRange"); */
/*   } */
/*   TF1*     GetBackgroundRecalcFunc(){ */
/*     return fhistoInvMass->GetFunction("funcbkgRecalc"); */
/*   } */
/*   TF1*     GetMassFunc(){ */
/*     return fhistoInvMass->GetFunction("funcmass"); */
/*   } */


 private:

  void     PlotFitVAR(TVirtualPad* pd,Double_t nsigma=3,Int_t writeFitInfo=1);

  void     ComputeParSize();
  void     ComputeNFinalPars();
  Bool_t   SideBandsBounds();
  Bool_t   CheckRangeFit();
  void     AddFunctionsToHisto();
  Int_t fNparSignal;           /// number of signal parameters
  Int_t fNparBack;           /// number of bkg parameters
  Int_t fNparRefl;          /// number of reflection parameters
  TH1F *fhTemplRefl;        /// template of reflection contribution
  TString *fSignParNames;   /// signal parameter names
  TString *fBackParNames;   /// back parameter names
  TString *fReflParNames;   /// refl parameter names
  Bool_t fSmoothRefl;      /// smoothing refl template
  Double_t fReflInit;       /// initial value of Refl/Signal
  Bool_t *fFixParSign;  /// fix signal parameter from ext value
  Bool_t *fFixParBack;  /// fix signal parameter from ext value
  Bool_t *fFixParRefl;  /// fix signal parameter from ext value
  Bool_t *fFixParSignExternalValue;  /// fix signal parameter from ext value
  Bool_t *fFixParBackExternalValue;  /// fix signal parameter from ext value
  Bool_t *fFixParReflExternalValue;  /// fix signal parameter from ext value
  Double_t *fparSignFixExt; /// external values to fix signal parameters
  Double_t *fparBackFixExt; /// external values to fix back parameters
  Double_t *fparReflFixExt; /// external values to fix refl parameters
  Double_t fRawYieldHelp;   /// internal variable used when fitting with reflections
  Int_t fpolbackdegreeTay; /// degree of polynomial expansion for back fit (option 6 for back)
  Int_t   fpolbackdegreeTayHelp; /// help variable
  Double_t fMassParticle;       /// pdg value of particle mass
/*   TH1F*     fhistoInvMass;     // histogram to fit */
/*   Double_t  fminMass;          // lower mass limit */
/*   Double_t  fmaxMass;          // upper mass limit */
/*   Int_t     fminBinMass;       // bin corresponding to fminMass */
/*   Int_t     fmaxBinMass;       // bin corresponding to fmaxMass */
/*   Int_t     fNbin;             // number of bins */
/*   Int_t     fParsSize;         // size of fFitPars array */
/*   Int_t     fNFinalPars;       // number of parameters of the final function */
/*   Float_t*  fFitPars;          //[fParsSize] array of fit parameters */
/*   Bool_t    fWithBkg;          // signal+background (kTRUE) or signal only (kFALSE) */
/*   Int_t     ftypeOfFit4Bkg;    // 0 = exponential; 1 = linear; 2 = pol2 */
/*   Int_t     ftypeOfFit4Sgn;    // 0 = gaus; 1 = gaus+gaus broadened */
/*   Int_t     ffactor;           // number to multiply to the sigma of the signal to obtain the reflected gaussian */
/*   TNtuple*  fntuParam;         // contains fit parameters */
/*   Double_t  fMass;             // signal gaussian mean value */
/*   Double_t  fMassErr;          // err signal gaussian mean value */
/*   Double_t  fSigmaSgn;         // signal gaussian sigma */
/*   Double_t  fSigmaSgnErr;      // err signal gaussian sigma */
/*   Double_t  fRawYield;         // signal gaussian integral */
/*   Double_t  fRawYieldErr;      // err on signal gaussian integral */
/*   Bool_t    fSideBands;        // kTRUE = only side bands considered */
/*   Bool_t*   fFixPar;           //[fNFinalPars] for each par if kTRUE it is fixed in fit */
/*   Int_t     fSideBandl;        // left side band limit (bin number) */
/*   Int_t     fSideBandr;        // right side band limit (bin number) */
/*   Int_t     fcounter;          // internal counter */
/*   TList*    fContourGraph;     // TList of TGraph containing contour plots */

  /// \cond CLASSIMP
  ClassDef(AliHFMassFitterVAR,2); /// class for invariant mass fit
  /// \endcond
};

#endif


