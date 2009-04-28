#ifndef ALIHFMASSFITTER_H
#define ALIHFMASSFITTER_H
/* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/////////////////////////////////////////////////////////////
//
// AliHFMassFitter for the fit of invariant mass distribution
// of charmed mesons
//
// Author: C.Bianchin, chiara.bianchin@pd.infn.it
/////////////////////////////////////////////////////////////

#include <TNamed.h>
#include <Riostream.h>
#include <TH1F.h>
#include <TF1.h>
#include <TMath.h>
#include <TNtuple.h>
#include <TFile.h>
#include <TString.h>


class AliHFMassFitter : public TNamed {

 public:
  AliHFMassFitter();
  AliHFMassFitter(TH1F* histoToFit, Double_t minvalue, Double_t maxvalue, Int_t rebin=1,Int_t fittypeb=0,Int_t fittypes=0);
  virtual ~AliHFMassFitter(){if(fhistoInvMass) delete fhistoInvMass; if(fntuParam) delete fntuParam;}

  AliHFMassFitter(const AliHFMassFitter &mfit);
  AliHFMassFitter& operator=(const AliHFMassFitter &mfit);

  //setters
  void     SetHisto(TH1F *histoToFit) {fhistoInvMass=histoToFit;}
  void     SetRangeFit(Double_t minvalue, Double_t maxvalue){fminMass=minvalue; fmaxMass=maxvalue;};
  void     SetMinRangeFit(Double_t minvalue){fminMass=minvalue;};
  void     SetMaxRangeFit(Double_t maxvalue){fmaxMass=maxvalue;};
  void     SetBinN(Int_t newbinN){fNbin=newbinN;};
  void     SetType(Int_t fittypeb, Int_t fittypes) {ftypeOfFit4Bkg = fittypeb; ftypeOfFit4Sgn = fittypes; };
  void     SetReflectionSigmaFactor(Int_t constant) {ffactor=constant;}
  void     SetInitialGaussianMean(Double_t mean) {fMass=mean;}
  void     SetInitialGaussianSigma(Double_t sigma) {fSigmaSgn=sigma;}
  void     SetSideBands(Bool_t onlysidebands=kTRUE) {fSideBands=onlysidebands;}

  //getters
  TH1F*    GetHisto()      const {return fhistoInvMass;};
  Double_t GetMinRangeFit()const {return fminMass;};
  Double_t GetMaxRangeFit()const {return fmaxMass;};
  Int_t    GetBinN()       const {return fNbin;};
  void     GetFitPars(Float_t*) const;
  void     GetTypeOfFit(Bool_t background, Int_t typeb) const {background = fWithBkg; typeb = ftypeOfFit4Bkg;}
  Int_t    GetReflectionSigmaFactor() const {return ffactor;} 

  void     InitNtuParam(char *ntuname="ntupar");
  void     FillNtuParam();
  TNtuple* GetNtuParam() const {return fntuParam;}
  TNtuple* NtuParamOneShot(char *ntuname="ntupar");

  void     Signal(Double_t nOfSigma,Double_t &signal,Double_t &errsignal) const; 
  void     Background(Double_t nOfSigma,Double_t &background,Double_t &errbackground) const; 
  void     Significance(Double_t nOfSigma,Double_t &significance,Double_t &errsignificance) const;

  Double_t FitFunction4MassDistr (Double_t*, Double_t*);
  Double_t FitFunction4Sgn (Double_t*, Double_t*);
  Double_t FitFunction4Bkg (Double_t*, Double_t*);
  void     MassFitter(Bool_t draw=kTRUE);
  void     RebinMass(Int_t binground=1);
  

 private:

  TH1F    *fhistoInvMass;  // histogram to fit
  Double_t fminMass;       // lower mass limit
  Double_t fmaxMass;       // upper mass limit
  Int_t    fNbin;          // number of bins
  Float_t  fFitPars[30];   // array of fit parameters
  Bool_t   fWithBkg;       // signal+background (kTRUE) or signal only (kFALSE)
  Int_t    ftypeOfFit4Bkg; // 0 = exponential; 1 = linear; 2 = pol2
  Int_t    ftypeOfFit4Sgn; // 0 = gaus; 1 = gaus+gaus broadened
  Int_t    ffactor;         // number to multiply to the sigma of the signal to obtain the reflected gaussian
  TNtuple *fntuParam;      // contains fit parameters
  Double_t fMass;          // signal gaussian mean value
  Double_t fSigmaSgn;      // signal gaussian sigma
  Bool_t   fSideBands;     // kTRUE = only side bands considered

  ClassDef(AliHFMassFitter,0); // class for invariant mass fit
};

#endif


