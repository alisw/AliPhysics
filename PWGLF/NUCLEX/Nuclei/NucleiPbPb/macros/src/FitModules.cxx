//
//  FitModules.cpp
//
//  Created by Maximiliano Puccio on 15/10/15.
//  Copyright © 2015 ALICE. All rights reserved.
//

#include "FitModules.h"
#include <RooDataHist.h>
#include <RooArgList.h>
#include <RooArgSet.h>
#include <RooExponential.h>
#include <RooRealVar.h>
#include <RooAddPdf.h>
#include <RooGaussian.h>
#include <RooGenericPdf.h>
#include <RooCBShape.h>
#include <RooPolynomial.h>
#include <RooFitResult.h>
#include <RooFormulaVar.h>
#include <RooArgSet.h>
#include <RooDataHist.h>
#include <RooHistPdf.h>
#include <RooMCStudy.h>
#include <RooExtendPdf.h>
#include <RooCurve.h>
#include <RooChebychev.h>
#include <RooBreitWigner.h>
#include <RooBernstein.h>
#include "RooGausDExp.h"


void FitModule::FitData(TH1* dat,TString &name, TString &title, TString range) {
  RooDataHist data("data","data",RooArgList(*mX),Import(*dat));
  mPlot = mX->frame();
  mPlot->SetTitle(title.Data());
  mPlot->SetName(name.Data());
  RooFitResult *res = mTemplate->fitTo(data,Extended(),Verbose(kFALSE),PrintLevel(-1),Range(range));
  data.plotOn(mPlot,Name("data"));
  mTemplate->plotOn(mPlot,Name("model"),Range(range),NormRange(range));
  mTemplate->plotOn(mPlot,Components(*mBackground),LineStyle(kDashed),Range(range),NormRange(range));
  mTemplate->plotOn(mPlot,Components(*mSignal),LineStyle(kDotted),Range(range),NormRange(range));
  mTemplate->paramOn(mPlot,Label(Form("#chi^{2}/NDF = %2.4f",mPlot->chiSquare("model","data"))));
  data.removeSelfFromDir();
}

FitExpGaus::FitExpGaus(RooRealVar *x) : FitModule(x) {
  mTau = new RooRealVar("mTau","tau bkg",-5.,0.);
  mBackground = new RooExponential("mBackground","Background",*mX,*mTau);
  mMu = new RooRealVar("mMu","Mu",-0.1,0.1);
  mSigma = new RooRealVar("mSigma","Sigma",0.01,0.6);
  mSignal = new RooGaussian("mSignal","Signal",*mX,*mMu,*mSigma);
  mTemplate = new RooAddPdf("mTemplate","Template",RooArgList(*mSignal,*mBackground),
                            RooArgList(*mSigCounts,*mBkgCounts));
}

FitExpGaus::~FitExpGaus() {
  delete mBackground;
  delete mMu;
  delete mSigma;
  delete mSignal;
  delete mBkgCounts;
  delete mSigCounts;
  delete mTemplate;
}

FitExpCB::FitExpCB(RooRealVar *x) : FitModule(x) {
  mTau = new RooRealVar("mTau","tau bkg",-5.,0.);
  mBackground = new RooExponential("mBackground","Background",*mX,*mTau);
  mAlpha = new RooRealVar("mAlpha","Alpha",-4.,-1.75);
  mMu = new RooRealVar("mMu","Mu",-0.1,0.1);
  mSigma = new RooRealVar("mSigma","Sigma",0.01,0.6);
  mN = new RooRealVar("mN","n",3.,10.);
  mSignal = new RooCBShape("mSignal","Signal",*mX,*mMu,*mSigma,*mAlpha,*mN);
  mTemplate = new RooAddPdf("mTemplate","Template",RooArgList(*mSignal,*mBackground),
                            RooArgList(*mSigCounts,*mBkgCounts));
  
}

FitExpCB::~FitExpCB() {
  delete mTau;
  delete mBackground;
  delete mAlpha;
  delete mMu;
  delete mSigma;
  delete mN;
  delete mSignal;
  delete mBkgCounts;
  delete mSigCounts;
  delete mTemplate;
}

FitExpExpCB::FitExpExpCB(RooRealVar *x) : FitModule(x) {
  mTau0 = new RooRealVar("mTau0","#tau_{0}",-6.5,-1.5);
  mTau1 = new RooRealVar("mTau1","#tau_{1}",-1.5,0.);
  mKbkg = new RooRealVar("mKbkg","mKbkg",0.,0.,1.);
  mKbkg->setConstant(false);
  mBkg0 = new RooExponential("mBkg0","background1",*mX,*mTau0);
  mBkg1 = new RooExponential("mBkg1","background2",*mX,*mTau1);
  mBackground = new RooAddPdf("mBackground","Background",RooArgList(*mBkg0,*mBkg1),RooArgList(*mKbkg));
  mAlpha = new RooRealVar("mAlpha","Alpha",-3.,-1.8);
  mMu = new RooRealVar("mMu","Mu",0.,-0.2,0.2);
  mMu->setConstant(false);
  mSigma = new RooRealVar("mSigma","Sigma",0.025,0.6);
  mN = new RooRealVar("mN","n",3.,40.);
  mSignal = new RooCBShape("mSignal","Signal",*mX,*mMu,*mSigma,*mAlpha,*mN);
  mTemplate = new RooAddPdf("mTemplate","Template",RooArgList(*mSignal,*mBackground),
                            RooArgList(*mSigCounts,*mBkgCounts));
}

FitExpExpCB::~FitExpExpCB() {
  delete mBackground;
  delete mTau0;
  delete mTau1;
  delete mKbkg;
  delete mMu;
  delete mSigma;
  delete mSignal;
  delete mAlpha;
  delete mN;
  delete mBkgCounts;
  delete mSigCounts;
  delete mTemplate;
}

FitExpExpGaus::FitExpExpGaus(RooRealVar *x) : FitModule(x) {
  mTau0 = new RooRealVar("mTau0","#tau_{0}",-6.5,-1.5);
  mTau1 = new RooRealVar("mTau1","#tau_{1}",-1.5,0.);
  mKbkg = new RooRealVar("mKbkg","mKbkg",0.,0.,1.);
  mKbkg->setConstant(false);
  mBkg0 = new RooExponential("mBkg0","background1",*mX,*mTau0);
  mBkg1 = new RooExponential("mBkg1","background2",*mX,*mTau1);
  mBackground = new RooAddPdf("mBackground","Background",RooArgList(*mBkg0,*mBkg1),RooArgList(*mKbkg));
  mMu = new RooRealVar("mMu","Mu",0.,-0.2,0.2);
  mMu->setConstant(false);
  mSigma = new RooRealVar("mSigma","Sigma",0.025,0.9);
  mSignal = new RooGaussian("mSignal","Signal",*mX,*mMu,*mSigma);
  mTemplate = new RooAddPdf("mTemplate","Template",RooArgList(*mSignal,*mBackground),
                            RooArgList(*mSigCounts,*mBkgCounts));
}

FitExpExpGaus::~FitExpExpGaus() {
  delete mBackground;
  delete mTau0;
  delete mTau1;
  delete mKbkg;
  delete mMu;
  delete mSigma;
  delete mSignal;
  delete mBkgCounts;
  delete mSigCounts;
  delete mTemplate;
}

FitExpExpTailGaus::FitExpExpTailGaus(RooRealVar *x) : FitModule(x) {
  mTau0 = new RooRealVar("mTau0","#tau_{0}",-6.5,-1.5);
  mTau1 = new RooRealVar("mTau1","#tau_{1}",-1.5,0.);
  mKbkg = new RooRealVar("mKbkg","mKbkg",0.,0.,1.);
  mKbkg->setConstant(false);
  mBkg0 = new RooExponential("mBkg0","background1",*mX,*mTau0);
  mBkg1 = new RooExponential("mBkg1","background2",*mX,*mTau1);
  mBackground = new RooAddPdf("mBackground","Background",RooArgList(*mBkg0,*mBkg1),RooArgList(*mKbkg));
  mAlpha0 = new RooRealVar("mAlpha0","Alpha0",-1.25,-1.25);
  mAlpha1 = new RooRealVar("mAlpha1","Alpha1",1.2,5.);
  mMu = new RooRealVar("mMu","Mu",0.,-0.2,0.2);
  mMu->setConstant(false);
  mSigma = new RooRealVar("mSigma","Sigma",0.025,0.6);
  mSignal = new RooGausDExp("mSignal","Signal",*mX,*mMu,*mSigma,*mAlpha0,*mAlpha1);
  mTemplate = new RooAddPdf("mTemplate","Template",RooArgList(*mSignal,*mBackground),
                            RooArgList(*mSigCounts,*mBkgCounts));
}

void FitExpExpTailGaus::SetHighPt(bool flag) {
  if (flag) {
    mAlpha0->setVal(-5.);
    mAlpha0->setConstant(true);
  } else {
    mAlpha0->setConstant(false);
    mAlpha0->setRange(-5.,-1.5);
    mAlpha0->setVal(-1.8);
  }
}


FitExpExpTailGaus::~FitExpExpTailGaus() {
  delete mBackground;
  delete mTau0;
  delete mTau1;
  delete mKbkg;
  delete mMu;
  delete mSigma;
  delete mSignal;
  delete mAlpha0;
  delete mAlpha1;
  delete mBkgCounts;
  delete mSigCounts;
  delete mTemplate;
}
