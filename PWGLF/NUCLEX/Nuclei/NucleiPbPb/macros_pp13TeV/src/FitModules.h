#ifndef FitModules_hpp
#define FitModules_hpp

#include <RooAbsData.h>
#include <RooAbsPdf.h>
#include <RooAddPdf.h>
#include <RooChebychev.h>
#include <RooExponential.h>
#include <RooGaussian.h>
#include <RooPlot.h>
#include <RooRealVar.h>

#include <TString.h>

#include <memory>

using namespace RooFit;

class FitModule {
public:
  FitModule(RooRealVar *xm, bool extended_likelihood=true)
  : mX(xm)
  , mBkgCounts(new RooRealVar("N_{bkg}","Bkg counts",1000.,0.,1.e10))
  , mSigCounts(new RooRealVar("N_{sig}","Sig counts",1000.,0.,1.e10))
  , mFraction(new RooRealVar("f_{sig}","Signal fraction",0.5,0.,1.))
  , mMu(new RooRealVar("#mu","Mu",0.,-0.15,0.15))
  , mSigma(new RooRealVar("#sigma","Sigma",0.2,0.02,.61))
    {
      mExtended = extended_likelihood;
      mNentries = 0.;
    }

  virtual RooPlot* FitData(TH1* h, TString name, TString title, TString range = "", TString plotrange="",bool isColored = false);

  static RooPlot* FitAndPlot(RooRealVar &x, RooAbsData &data, RooAbsPdf &model, RooAbsPdf &sig,
                             RooAbsPdf &bkg,TString range,TString plotrange) {
    RooPlot *plot = x.frame();
    RooFitResult *res = model.fitTo(data,Extended(),Range(range),Verbose(kFALSE),PrintLevel(-1));
    data.plotOn(plot,Name("data"));
    model.plotOn(plot,Name("model"),Range(plotrange),NormRange(plotrange));
    model.plotOn(plot,Components(bkg),LineStyle(kDashed),Range(plotrange),NormRange(plotrange));
    model.plotOn(plot,Components(sig),LineStyle(kDotted),Range(plotrange),NormRange(plotrange));
    model.paramOn(plot,Label(Form("#chi^{2}/NDF = %2.4f",plot->chiSquare("model","data"))));
    return plot;
  }

  void UseBackground(bool);
  void UseSignal(bool);

  RooRealVar  *mX;
  unique_ptr<RooRealVar> mFraction;
  unique_ptr<RooRealVar> mBkgCounts;
  unique_ptr<RooRealVar> mSigCounts;
  unique_ptr<RooRealVar> mMu;
  unique_ptr<RooRealVar> mSigma;
  unique_ptr<RooAddPdf>  mTemplate;
  unique_ptr<RooAbsPdf>  mSignal;
  unique_ptr<RooAbsPdf>  mBackground;
  float        mChi2;
  bool         mExtended;
  double       mNentries;
};

class FitGausGaus : public FitModule {
public:
  FitGausGaus(RooRealVar *xm, bool extended_likelihood=true);

  unique_ptr<RooRealVar> mMuBkg;
  unique_ptr<RooRealVar> mSigmaBkg;

};

class FitExpTailGaus : public FitModule {
public:
  FitExpTailGaus(RooRealVar *xm, bool extended_likelihood=true);

  unique_ptr<RooRealVar> mTau0;
  unique_ptr<RooRealVar> mAlpha0;
};

class FitExpTailTailGaus : public FitModule {
public:
  FitExpTailTailGaus(RooRealVar *xm, bool extended_likelihood=true);

  unique_ptr<RooRealVar> mTau0;
  unique_ptr<RooRealVar> mAlpha0;
  unique_ptr<RooRealVar> mAlpha1;
};

class FitExpGaus : public FitModule {
public:
  FitExpGaus(RooRealVar *xm, bool extended_likelihood=true);

  unique_ptr<RooRealVar> mTau;
  unique_ptr<RooRealVar> mA1;
  unique_ptr<RooRealVar> mA2;
  unique_ptr<RooRealVar> mA3;

};

class FitExpCB : public FitModule {
public:
  FitExpCB(RooRealVar *xm, bool extended_likelihood=true);

  unique_ptr<RooRealVar> mTau;
  unique_ptr<RooRealVar> mAlpha;
  unique_ptr<RooRealVar> mN;
};

class FitExpExpCB : public FitModule {
public:
  FitExpExpCB(RooRealVar *xm, bool extended_likelihood=true);

  unique_ptr<RooExponential> mBkg0;
  unique_ptr<RooExponential> mBkg1;
  unique_ptr<RooRealVar> mTau0;
  unique_ptr<RooRealVar> mTau1;
  unique_ptr<RooRealVar> mKbkg;
  unique_ptr<RooRealVar> mAlpha;
  unique_ptr<RooRealVar> mN;
};

class FitExpExpGaus : public FitModule {
public:
  FitExpExpGaus(RooRealVar *xm, bool extended_likelihood=true);

  unique_ptr<RooExponential> mBkg0;
  unique_ptr<RooExponential> mBkg1;
  unique_ptr<RooRealVar> mTau0;
  unique_ptr<RooRealVar> mTau1;
  unique_ptr<RooRealVar> mKbkg;
};

class FitExpExpTailGaus : public FitModule {
public:
  FitExpExpTailGaus(RooRealVar *xm, bool extended_likelihood=true);

  unique_ptr<RooExponential> mBkg0;
  unique_ptr<RooExponential> mBkg1;
  unique_ptr<RooRealVar> mTau0;
  unique_ptr<RooRealVar> mTau1;
  unique_ptr<RooRealVar> mKbkg;
  unique_ptr<RooRealVar> mAlpha0;
};

class FitExpExpTailTailGaus : public FitModule {
public:
  FitExpExpTailTailGaus(RooRealVar *xm, bool extended_likelihood=true);

  unique_ptr<RooExponential> mBkg0;
  unique_ptr<RooExponential> mBkg1;
  unique_ptr<RooRealVar> mTau0;
  unique_ptr<RooRealVar> mTau1;
  unique_ptr<RooRealVar> mKbkg;
  unique_ptr<RooRealVar> mAlpha0;
  unique_ptr<RooRealVar> mAlpha1;
};


class FitExpPolTailGaus : public FitModule {
public:
  FitExpPolTailGaus(RooRealVar *xm, bool extended_likelihood=true);

  unique_ptr<RooExponential> mBkg0;
  unique_ptr<RooChebychev> mBkg1;
  unique_ptr<RooRealVar> mTau0;
  unique_ptr<RooRealVar> mA0;
  unique_ptr<RooRealVar> mA1;
  unique_ptr<RooRealVar> mKbkg;
  unique_ptr<RooRealVar> mAlpha0;
};

class FitGausExpTailGaus : public FitModule {
public:
  FitGausExpTailGaus(RooRealVar *xm, bool extended_likelihood=true);

  unique_ptr<RooGaussian> mBkg0;
  unique_ptr<RooExponential> mBkg1;
  unique_ptr<RooRealVar> mMuBkg;
  unique_ptr<RooRealVar> mSigmaBkg;
  unique_ptr<RooRealVar> mTau0;
  unique_ptr<RooRealVar> mKbkg;
  unique_ptr<RooRealVar> mAlpha0;
};

class FitTailGausTailGaus : public FitModule {
public:
  FitTailGausTailGaus(RooRealVar *xm, bool extended_likelihood=true);

  unique_ptr<RooRealVar> mMuBkg;
  unique_ptr<RooRealVar> mSigmaBkg;
  unique_ptr<RooRealVar> mAlpha0Bkg;
  unique_ptr<RooRealVar> mAlpha0;
};

#endif /* FitModules_hpp */
