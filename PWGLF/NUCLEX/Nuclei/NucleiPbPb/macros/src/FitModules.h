#ifndef FitModules_hpp
#define FitModules_hpp

#include <RooAbsData.h>
#include <RooAbsPdf.h>
#include <RooAddPdf.h>
#include <RooChebychev.h>
#include <RooExponential.h>
#include <RooPlot.h>
#include <RooRealVar.h>

#include <TString.h>

#include <memory>

using namespace RooFit;

class FitModule {
public:
  FitModule(RooRealVar *xm)
  : mX(xm)
  , mBkgCounts(new RooRealVar("mBkgCounts","Bkg counts",1000.,0.,1.e8))
  , mSigCounts(new RooRealVar("mSigCounts","Sig counts",1000.,0.,1.e8))
  , mMu(new RooRealVar("mMu","Mu",-0.1,0.1))
  , mSigma(new RooRealVar("mSigma","Sigma",0.02,.61))
    {}

  virtual RooPlot* FitData(TH1* h, TString name, TString title, TString range = "", TString plotrange="", bool change_range = false, float low_x = -2., float high_x = 2.);
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
  std::unique_ptr<RooRealVar> mBkgCounts;
  std::unique_ptr<RooRealVar> mSigCounts;
  std::unique_ptr<RooRealVar> mMu;
  std::unique_ptr<RooRealVar> mSigma;
  std::unique_ptr<RooAddPdf>  mTemplate;
  std::unique_ptr<RooAbsPdf>  mSignal;
  std::unique_ptr<RooAbsPdf>  mBackground;
  float        mChi2;
};

class FitGausGaus : public FitModule {
public:
  FitGausGaus(RooRealVar *xm);

  std::unique_ptr<RooRealVar> mMuBkg;
  std::unique_ptr<RooRealVar> mSigmaBkg;

};

class FitExpTailGaus : public FitModule {
public:
  FitExpTailGaus(RooRealVar *xm);

  std::unique_ptr<RooRealVar> mTau0;
  std::unique_ptr<RooRealVar> mAlpha0;
};

class FitExpGaus : public FitModule {
public:
  FitExpGaus(RooRealVar *xm);

  std::unique_ptr<RooRealVar> mTau;
  std::unique_ptr<RooRealVar> mA1;
  std::unique_ptr<RooRealVar> mA2;
  std::unique_ptr<RooRealVar> mA3;

};

class FitExpCB : public FitModule {
public:
  FitExpCB(RooRealVar *xm);

  std::unique_ptr<RooRealVar> mTau;
  std::unique_ptr<RooRealVar> mAlpha;
  std::unique_ptr<RooRealVar> mN;
};

class FitExpExpCB : public FitModule {
public:
  FitExpExpCB(RooRealVar *xm);

  std::unique_ptr<RooExponential> mBkg0;
  std::unique_ptr<RooExponential> mBkg1;
  std::unique_ptr<RooRealVar> mTau0;
  std::unique_ptr<RooRealVar> mTau1;
  std::unique_ptr<RooRealVar> mKbkg;
  std::unique_ptr<RooRealVar> mAlpha;
  std::unique_ptr<RooRealVar> mN;
};

class FitExpExpGaus : public FitModule {
public:
  FitExpExpGaus(RooRealVar *xm);

  std::unique_ptr<RooExponential> mBkg0;
  std::unique_ptr<RooExponential> mBkg1;
  std::unique_ptr<RooRealVar> mTau0;
  std::unique_ptr<RooRealVar> mTau1;
  std::unique_ptr<RooRealVar> mKbkg;
};

class FitExpExpTailGaus : public FitModule {
public:
  FitExpExpTailGaus(RooRealVar *xm);

  std::unique_ptr<RooExponential> mBkg0;
  std::unique_ptr<RooExponential> mBkg1;
  std::unique_ptr<RooRealVar> mTau0;
  std::unique_ptr<RooRealVar> mTau1;
  std::unique_ptr<RooRealVar> mKbkg;
  std::unique_ptr<RooRealVar> mAlpha0;
};

class FitExpExpTailTailGaus : public FitModule {
public:
  FitExpExpTailTailGaus(RooRealVar *xm);

  std::unique_ptr<RooExponential> mBkg0;
  std::unique_ptr<RooExponential> mBkg1;
  std::unique_ptr<RooRealVar> mTau0;
  std::unique_ptr<RooRealVar> mTau1;
  std::unique_ptr<RooRealVar> mKbkg;
  std::unique_ptr<RooRealVar> mAlpha0;
  std::unique_ptr<RooRealVar> mAlpha1;
};


class FitExpPolTailGaus : public FitModule {
public:
  FitExpPolTailGaus(RooRealVar *xm);

  std::unique_ptr<RooExponential> mBkg0;
  std::unique_ptr<RooChebychev> mBkg1;
  std::unique_ptr<RooRealVar> mTau0;
  std::unique_ptr<RooRealVar> mA0;
  std::unique_ptr<RooRealVar> mA1;
  std::unique_ptr<RooRealVar> mKbkg;
  std::unique_ptr<RooRealVar> mAlpha0;
};

#endif /* FitModules_hpp */
