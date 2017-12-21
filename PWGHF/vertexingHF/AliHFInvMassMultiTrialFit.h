#ifndef ALIHFINVMASSMULTITRIALFIT_H
#define ALIHFINVMASSMULTITRIALFIT_H
/* Copyright(c) 2008-2019, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <TNamed.h>
#include <TString.h>
#include <TPad.h>
#include <set>
#include <vector>

class TNtuple;
class AliHFInvMassFitter;

/// \class AliHFInvMassMultiTrialFit

class AliHFInvMassMultiTrialFit : public TNamed {

 public:
  AliHFInvMassMultiTrialFit();
  virtual ~AliHFInvMassMultiTrialFit();

  void ConfigureRebinSteps(Int_t nSteps, Int_t*  values){
    if(fRebinSteps) delete [] fRebinSteps;
    fNumOfRebinSteps=nSteps;
    fRebinSteps = new Int_t[fNumOfRebinSteps];
    for(Int_t ib=0; ib<fNumOfRebinSteps; ib++) fRebinSteps[ib]=values[ib];
  }

  void SetNumOfFirstBinSteps(Int_t nfst){
    fNumOfFirstBinSteps=nfst;
  }

  void ConfigureLowLimFitSteps(Int_t nSteps, Double_t*  values){
    if(fLowLimFitSteps) delete [] fLowLimFitSteps;
    fNumOfLowLimFitSteps=nSteps;
    fLowLimFitSteps = new Double_t[fNumOfLowLimFitSteps];
    for(Int_t ib=0; ib<fNumOfLowLimFitSteps; ib++) fLowLimFitSteps[ib]=values[ib];
  }

  void ConfigureUpLimFitSteps(Int_t nSteps, Double_t*  values){
    if(fUpLimFitSteps) delete [] fUpLimFitSteps;
    fNumOfUpLimFitSteps=nSteps;
    fUpLimFitSteps = new Double_t[fNumOfUpLimFitSteps];
    for(Int_t ib=0; ib<fNumOfUpLimFitSteps; ib++) fUpLimFitSteps[ib]=values[ib];
  }

  void ConfigurenSigmaBinCSteps(Int_t nSteps, Double_t*  values){
    if(fnSigmaBinCSteps) delete [] fnSigmaBinCSteps;
    fNumOfnSigmaBinCSteps=nSteps;
    fnSigmaBinCSteps = new Double_t[fNumOfnSigmaBinCSteps];
    for(Int_t ib=0; ib<fNumOfnSigmaBinCSteps; ib++) fnSigmaBinCSteps[ib]=values[ib];
  }

  void GetGlobalMinMaxYield(Double_t& min, Double_t& max) {
    min = fMinYieldGlob;
    max = fMaxYieldGlob;
  }

  void SetMass(Double_t mass){fMassD=mass;}
  void SetSigmaGaussMC(Double_t sig){fSigmaGausMC=sig;}
  void SetSigmaMCVariation(Double_t var=0.15){fSigmaMCVariation=var;}

  void SetSuffixForHistoNames(const Char_t* name){
    fSuffix=name;
  }
  void SetUseChi2Fit(){fFitOption=1;}
  void SetUseLogLikelihoodFit(){fFitOption=0;}


  void SetUseExpoBackground(Bool_t opt=kTRUE){fUseExpoBkg=opt;}
  void SetUseLinBackground(Bool_t opt=kTRUE){fUseLinBkg=opt;}
  void SetUsePol2Background(Bool_t opt=kTRUE){fUsePol2Bkg=opt;}
  void SetUsePol3Background(Bool_t opt=kTRUE){fUsePol3Bkg=opt;}
  void SetUsePol4Background(Bool_t opt=kTRUE){fUsePol4Bkg=opt;}
  void SetUsePol5Background(Bool_t opt=kTRUE){fUsePol5Bkg=opt;}
  void SetUsePowerLawBackground(Bool_t opt=kTRUE){fUsePowLawBkg=opt;}
  void SetUsePowerLawTimesExpoBackground(Bool_t opt=kTRUE){fUsePowLawTimesExpoBkg=opt;}

  void SetUseFixSigUpFreeMean(Bool_t opt=kTRUE) {fUseFixSigUpFreeMean=opt;}
  void SetUseFixSigDownFreeMean(Bool_t opt=kTRUE) {fUseFixSigDownFreeMean=opt;}
  void SetUseFreeS(Bool_t opt=kTRUE) {fUseFreeS=opt;}
  void SetUseFixedMeanFreeS(Bool_t opt=kTRUE) {fUseFixedMeanFreeS=opt;}
  void SetUseFixSigFreeMean(Bool_t opt=kTRUE) {fUseFixSigFreeMean=opt;}
  void SetUseFixSigFixMean(Bool_t opt=kTRUE) {fUseFixSigFixMean=opt;}

  void SetSaveBkgValue(Bool_t opt=kTRUE, Double_t nsigma=3) {fSaveBkgVal=opt; fnSigmaForBkgEval=nsigma;}

  void SetDrawIndividualFits(Bool_t opt=kTRUE){fDrawIndividualFits=opt;}

  Bool_t DoMultiTrials(TH1D* hInvMassHisto, TPad* thePad=0x0);
  void SaveToRoot(TString fileName, TString option="recreate") const;
  void DrawHistos(TCanvas* cry) const;
  void SetTemplatesForReflections(const TH1F *hTemplRefl, const TH1F *hTemplSig);

  void SetFixRefoS(Float_t refloS){fFixRefloS=refloS;}

  void IncludeSecondGausPeak(Double_t mass, Bool_t fixm, Double_t width, Bool_t fixw){
    fUseSecondPeak=kTRUE;
    fMassSecondPeak=mass;   fFixMassSecondPeak=fixm;
    fSigmaSecondPeak=width; fFixSigmaSecondPeak=fixw;
  }

  void AddInvMassFitSaveAsFormat(std::string format) { fInvMassFitSaveAsFormats.insert(format); }
  void DisableInvMassFitSaveAs() { fInvMassFitSaveAsFormats.clear(); }


  enum EBkgFuncCases{ kExpoBkg, kLinBkg, kPol2Bkg, kPol3Bkg, kPol4Bkg, kPol5Bkg, kPowBkg, kPowTimesExpoBkg, kNBkgFuncCases };
  enum EFitParamCases{ kFixSigFreeMean, kFixSigUpFreeMean, kFixSigDownFreeMean, kFreeSigFreeMean, kFixSigFixMean, kFreeSigFixMean, kNFitConfCases};

 private:

  Bool_t CreateHistos();
  Bool_t DoFitWithPol3Bkg(TH1F* histoToFit, Double_t  hmin, Double_t  hmax,
			  Int_t theCase);

  AliHFInvMassMultiTrialFit(const AliHFInvMassMultiTrialFit &source);
  AliHFInvMassMultiTrialFit& operator=(const AliHFInvMassMultiTrialFit& source);

  std::set<std::string> fInvMassFitSaveAsFormats; /// saves the invariant mass fit canvases in the file formats listed in this vector (if empty, does nothing)
  Int_t fNumOfRebinSteps; /// number of rebin steps
  Int_t* fRebinSteps;     //[fNumOfRebinSteps] values of rebin
  Int_t fNumOfFirstBinSteps; /// number of steps in the first bin for rebin
  Int_t fNumOfLowLimFitSteps; /// number of steps on the min. mass for fit
  Double_t* fLowLimFitSteps; //[fNumOfLowLimFitSteps] values of low limits for fit
  Int_t fNumOfUpLimFitSteps; /// number of steps on the max. mass for fit
  Double_t* fUpLimFitSteps; //[fNumOfUpLimFitSteps] values of up limits for fit
  Int_t fNumOfnSigmaBinCSteps; /// number of steps on the bin counting
  Double_t* fnSigmaBinCSteps; //[fNumOfnSigmaBinCSteps] values of nsigma for bin count
  Double_t fnSigmaForBkgEval; //value of sigma in which to extract bkg value

  Double_t fSigmaGausMC; /// sigma of D meson peak from MC
  Double_t fSigmaMCVariation; /// relative variation of the sigma
  Double_t fMassD;       /// mass of D meson
  TString fSuffix;       /// name to characterize analysis case
  Int_t fFitOption;      /// LL or chi2 fit
  Bool_t fUseExpoBkg;    /// switch for exponential background
  Bool_t fUseLinBkg;    /// switch for linear background
  Bool_t fUsePol2Bkg;    /// switch for pol2 background
  Bool_t fUsePol3Bkg;    /// switch for pol3 background
  Bool_t fUsePol4Bkg;    /// switch for pol4 background
  Bool_t fUsePol5Bkg;    /// switch for pol5 background
  Bool_t fUsePowLawBkg;  /// switch for power law background
  Bool_t fUsePowLawTimesExpoBkg;  /// switch for power law background
  Bool_t fUseFixSigUpFreeMean;    /// switch for FixSigUpFreeMean
  Bool_t fUseFixSigDownFreeMean;  /// switch for FixSigDownFreeMean
  Bool_t fUseFreeS;              /// switch for FreeSigma
  Bool_t fUseFixedMeanFreeS;     ///  switch for FixedMeanFreeS
  Bool_t fUseFixSigFreeMean;     ///  switch for FixSigFreeMean
  Bool_t fUseFixSigFixMean;      ///  switch for FixSigFixMean

  Bool_t   fUseSecondPeak;      /// switch off/on second peak (for D+->KKpi in Ds)
  Double_t fMassSecondPeak;     /// position of the 2nd peak
  Double_t fSigmaSecondPeak;    /// width of the 2nd peak
  Bool_t   fFixMassSecondPeak;  /// flag to fix the position of the 2nd peak
  Bool_t   fFixSigmaSecondPeak; /// flag to fix the width of the 2nd peak

  Bool_t fSaveBkgVal;		/// switch for saving bkg values in nsigma

  Bool_t fDrawIndividualFits; /// flag for drawing fits

  TH1F* fHistoRawYieldDistAll;  /// histo with yield from all trials
  TH1F* fHistoRawYieldTrialAll; /// histo with yield from all trials
  TH1F* fHistoSigmaTrialAll;    /// histo with gauss sigma from all trials
  TH1F* fHistoMeanTrialAll;     /// histo with gauss mean from all trials
  TH1F* fHistoChi2TrialAll;     /// histo with chi2 from all trials
  TH1F* fHistoSignifTrialAll;     /// histo with chi2 from all trials
  TH1F* fHistoBkgTrialAll;     /// histo with bkg from all trials
  TH1F* fHistoBkgInBinEdgesTrialAll;    /// histo with bkg in mass bin edges from all trials

  TH1F* fHistoRawYieldDistBinC0All; /// histo with bin counts from all trials
  TH2F* fHistoRawYieldTrialBinC0All; /// histo with bin counts from all trials
  TH1F* fHistoRawYieldDistBinC1All; /// histo with bin counts from all trials
  TH2F* fHistoRawYieldTrialBinC1All; /// histo with bin counts from all trials

  TH1F** fHistoRawYieldDist;  /// histo with yield from subsamples of trials
  TH1F** fHistoRawYieldTrial; /// histo with yield from subsamples of trials
  TH1F** fHistoSigmaTrial;    /// histo with gauss sigma from subsamples of trials
  TH1F** fHistoMeanTrial;     /// histo with gauss mean from subsamples of trials
  TH1F** fHistoChi2Trial;     /// histo with chi2 from subsamples of trials
  TH1F** fHistoSignifTrial;     /// histo with chi2 from subsamples of trials
  TH1F** fHistoBkgTrial;     /// histo with bkg from subsamples of trials
  TH1F** fHistoBkgInBinEdgesTrial;    /// histo with bkg in mass bin edges from subsamples of trials

  TH1F** fHistoRawYieldDistBinC0;  /// histo with bin counts from subsamples of trials
  TH2F** fHistoRawYieldTrialBinC0; /// histo with bin counts from subsamples of trials
  TH1F** fHistoRawYieldDistBinC1;  /// histo with bin counts from subsamples of trials
  TH2F** fHistoRawYieldTrialBinC1; /// histo with bin counts from subsamples of trials

  TH1F *fhTemplRefl;        /// template of reflection contribution
  TH1F *fhTemplSign;        /// template of signal contribution
  Float_t fFixRefloS;
  TNtuple* fNtupleMultiTrials; /// tree

  Double_t fMinYieldGlob;   /// minimum yield
  Double_t fMaxYieldGlob;   /// maximum yield

  std::vector<AliHFInvMassFitter*> fMassFitters; //!<! Mass fitters

  /// \cond CLASSIMP
  ClassDef(AliHFInvMassMultiTrialFit,2); /// class for multiple trials of invariant mass fit
  /// \endcond
};

#endif

