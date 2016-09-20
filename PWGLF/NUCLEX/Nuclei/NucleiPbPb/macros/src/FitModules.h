//
//  FitModules.hpp
//  TrackerCA
//
//  Created by Maximiliano Puccio on 15/10/15.
//  Copyright © 2015 ALICE. All rights reserved.
//

#ifndef FitModules_hpp
#define FitModules_hpp

#include <RooGlobalFunc.h>
#include <RooRealVar.h>
#include <RooAbsPdf.h>
#include <RooAbsData.h>
#include <RooPlot.h>
#include <TString.h>

using namespace RooFit;

class RooPlot;
class RooAddPdf;
class RooExponential;
class RooChebychev;

class FitModule {
public:
  FitModule(RooRealVar *xm)
  : mX(xm)
  , mBkgCounts(new RooRealVar("mBkgCounts","Bkg counts",1000.,0.,1.e8))
  , mSigCounts(new RooRealVar("mSigCounts","Sig counts",1000.,0.,1.e8))
  , mMu(new RooRealVar("mMu","Mu",-0.1,0.1))
  , mSigma(new RooRealVar("mSigma","Sigma",0.02,.61))
    {}

  virtual void FitData(TH1* h, TString &name, TString &title, TString range = "");
  RooPlot* GetPlot() { return mPlot; }
  static RooPlot* FitAndPlot(RooRealVar &x, RooAbsData &data, RooAbsPdf &model, RooAbsPdf &sig,
                             RooAbsPdf &bkg,TString range) {
    RooPlot *plot = x.frame();
    RooFitResult *res = model.fitTo(data,Extended(),Range(range),Verbose(kFALSE),PrintLevel(-1));
    data.plotOn(plot,Name("data"));
    model.plotOn(plot,Name("model"),Range(range),NormRange(range));
    model.plotOn(plot,Components(bkg),LineStyle(kDashed),Range(range),NormRange(range));
    model.plotOn(plot,Components(sig),LineStyle(kDotted),Range(range),NormRange(range));
    model.paramOn(plot,Label(Form("#chi^{2}/NDF = %2.4f",plot->chiSquare("model","data"))));
    return plot;
  }

  virtual void SetHighPt(bool) {/* flag) {
    if (flag) {
      mSigma->setVal(0.63); // Set looking at sigma trending
      mSigma->setConstant(true);
    }
    else {
      mSigma->setRange(0.01,0.63);
      mSigma->setConstant(false);
    }*/
  }

  RooRealVar  *mX;
  RooRealVar  *mBkgCounts;
  RooRealVar  *mSigCounts;
  RooRealVar  *mMu;
  RooRealVar  *mSigma;
  RooAddPdf   *mTemplate;
  RooPlot     *mPlot;
  RooAbsPdf   *mSignal;
  RooAbsPdf   *mBackground;
  float        mChi2;
};

class FitExpGaus : public FitModule {
public:
  FitExpGaus(RooRealVar *xm);
  virtual ~FitExpGaus();

  RooRealVar  *mTau;
  RooRealVar  *mA1;
  RooRealVar  *mA2;
  RooRealVar  *mA3;

};

class FitExpCB : public FitModule {
public:
  FitExpCB(RooRealVar *xm);
  virtual ~FitExpCB();

  RooRealVar  *mTau;
  RooRealVar  *mAlpha;
  RooRealVar  *mN;
};

class FitExpExpCB : public FitModule {
public:
  FitExpExpCB(RooRealVar *xm);
  virtual ~FitExpExpCB();

  RooExponential *mBkg0;
  RooExponential *mBkg1;
  RooRealVar  *mTau0;
  RooRealVar  *mTau1;
  RooRealVar  *mKbkg;
  RooRealVar  *mAlpha;
  RooRealVar  *mN;
};

class FitExpExpGaus : public FitModule {
public:
  FitExpExpGaus(RooRealVar *xm);
  virtual ~FitExpExpGaus();

  RooExponential *mBkg0;
  RooExponential *mBkg1;
  RooRealVar  *mTau0;
  RooRealVar  *mTau1;
  RooRealVar  *mKbkg;
};

class FitExpExpTailGaus : public FitModule {
public:
  FitExpExpTailGaus(RooRealVar *xm);
  virtual ~FitExpExpTailGaus();
  virtual void SetHighPt(bool);

  RooExponential *mBkg0;
  RooExponential *mBkg1;
  RooRealVar  *mTau0;
  RooRealVar  *mTau1;
  RooRealVar  *mKbkg;
  RooRealVar  *mAlpha0;
};

class FitExpPolTailGaus : public FitModule {
public:
  FitExpPolTailGaus(RooRealVar *xm);
  virtual ~FitExpPolTailGaus();
  virtual void SetHighPt(bool);

  RooExponential *mBkg0;
  RooChebychev *mBkg1;
  RooRealVar  *mTau0;
  RooRealVar  *mA0;
  RooRealVar  *mA1;
  RooRealVar  *mKbkg;
  RooRealVar  *mAlpha0;
};

#endif /* FitModules_hpp */
