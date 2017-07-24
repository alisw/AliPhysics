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
#include "RooGausExp.h"
#include "RooGausDExp.h"


void FitModule::FitData(TH1* dat,TString name, TString title, TString range, TString plotrange, bool change_range, float low_x, float high_x) {
  if (plotrange == "") plotrange = range;
  RooDataHist data("data","data",RooArgList(*mX),Import(*dat));
  mPlot = mX->frame();
  mPlot->SetTitle(title.Data());
  mPlot->SetName(name.Data());
  mPlot->GetYaxis()->SetTitle(Form("Counts / (%.2f Gev/#it{c}^{2})",mPlot->GetXaxis()->GetBinWidth(1)));
  for (int i = 2; i--;) RooFitResult *res = mTemplate->fitTo(data,Extended(),Verbose(kFALSE),PrintLevel(-1),Range(range));
  if(change_range) mPlot->GetXaxis()->SetRangeUser(low_x,high_x);
  data.plotOn(mPlot,Name("data"),DrawOption("pz"));
  mTemplate->plotOn(mPlot,Name("model"),Range(range),NormRange(range));
  mTemplate->plotOn(mPlot,Name("bkg"),Components(*mBackground),LineStyle(kDashed),LineColor(kRed),Range(range),NormRange(range));
  //mTemplate->plotOn(mPlot,Components(*mSignal),LineStyle(kDotted),LineColor(kGreen+3),Range(range),NormRange(range));
  mChi2 = mPlot->chiSquare("model","data");
  mPlot->remove("model",false);
  mPlot->remove("bkg",false);
  mTemplate->plotOn(mPlot,Name("model"),Range(plotrange),NormRange(range));
  mTemplate->plotOn(mPlot,Name("bkg"),Components(*mBackground),LineStyle(kDashed),LineColor(kRed),Range(plotrange),NormRange(range));
  //mTemplate->plotOn(mPlot,Components(*mSignal),LineStyle(kDotted),LineColor(kGreen+3),Range(plotrange),NormRange(range));
  mTemplate->paramOn(mPlot,Label(Form("#chi^{2}/NDF = %2.4f",mChi2)));
  data.removeSelfFromDir();
}

void FitModule::UseBackground(bool useBkg){
  if(useBkg==false){
    mBkgCounts->setVal(0);
    mBkgCounts->setConstant(true);
  }
  else{
    mBkgCounts->setVal(1000);
    mBkgCounts->setConstant(false);
  }
  RooArgSet* bkg_params = mBackground->getVariables();
  auto iter = bkg_params->createIterator();
  RooRealVar* var = (RooRealVar*)iter->Next();
  while(var){
    var->setConstant(!useBkg);
    var = (RooRealVar*)iter->Next();
  }
}

void FitModule::UseSignal(bool useSig){
  if(useSig==false){
    mSigCounts->setVal(0);
    mSigCounts->setConstant(true);
  }
  else{
    mSigCounts->setVal(1000);
    mSigCounts->setConstant(false);
  }
  RooArgSet* sig_params = mSignal->getVariables();
  auto iter = sig_params->createIterator();
  RooRealVar* var = (RooRealVar*)iter->Next();
  while(var){
    var->setConstant(!useSig);
    var = (RooRealVar*)iter->Next();
  }
}

FitGausGaus::FitGausGaus(RooRealVar *x) : FitModule(x){
  mMuBkg = new RooRealVar("mMuBkg","mu bkg",-6,6); // Better to set the range
  mSigmaBkg = new RooRealVar("mSigmaBkg","sigma bkg",0.01,2.);
  mBackground = new RooGaussian("mBackground","Background",*mX,*mMuBkg,*mSigmaBkg);
  mSignal = new RooGaussian("mSignal","Signal",*mX,*mMu,*mSigma);
  mTemplate = new RooAddPdf("mTemplate","Template",RooArgList(*mSignal,*mBackground),
                            RooArgList(*mSigCounts,*mBkgCounts));
}

FitGausGaus::~FitGausGaus() {
  delete mMuBkg;
  delete mSigmaBkg;
  delete mBackground;
  delete mMu;
  delete mSigma;
  delete mSignal;
  delete mBkgCounts;
  delete mSigCounts;
  delete mTemplate;
}

FitExpTailGaus::FitExpTailGaus(RooRealVar *x) : FitModule(x){
  mTau0 = new RooRealVar("mTau0","tau",-10.,-0.00001);
  mBackground = new RooExponential("mBackground","Background",*mX,*mTau0);
  mAlpha0 = new RooRealVar("mAlpha0","Alpha0",1.6,3.);
  mSignal = new RooGausExp("mSignal","Signal",*mX,*mMu,*mSigma,*mAlpha0);
  mTemplate = new RooAddPdf("mTemplate","Template",RooArgList(*mSignal,*mBackground),
                            RooArgList(*mSigCounts,*mBkgCounts));
}

FitExpTailGaus::~FitExpTailGaus() {
  delete mTau0;
  delete mBackground;
  delete mMu;
  delete mSigma;
  delete mAlpha0;
  delete mSignal;
  delete mBkgCounts;
  delete mSigCounts;
  delete mTemplate;
}

FitExpGaus::FitExpGaus(RooRealVar *x) : FitModule(x) {
  mTau = new RooRealVar("mTau","tau bkg",-5.,0.);
  mBackground = new RooExponential("mBackground","Background",*mX,*mTau);
  mSignal = new RooGaussian("mSignal","Signal",*mX,*mMu,*mSigma);
  mTemplate = new RooAddPdf("mTemplate","Template",RooArgList(*mSignal,*mBackground),
                            RooArgList(*mSigCounts,*mBkgCounts));
}

FitExpGaus::~FitExpGaus() {
  delete mTau;
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
  mTau0 = new RooRealVar("mTau0","#tau_{0}",-6.5,-2.);
  mTau1 = new RooRealVar("mTau1","#tau_{1}",-2.,-0.2);
  mKbkg = new RooRealVar("mKbkg","mKbkg",0.,0.,1.);
  mKbkg->setConstant(false);
  mBkg0 = new RooExponential("mBkg0","background1",*mX,*mTau0);
  mBkg1 = new RooExponential("mBkg1","background2",*mX,*mTau1);
  mBackground = new RooAddPdf("mBackground","Background",RooArgList(*mBkg0,*mBkg1),RooArgList(*mKbkg));
  mAlpha = new RooRealVar("mAlpha","Alpha",-3.,-0.8);
  mMu->setConstant(false);
  mN = new RooRealVar("mN","n",3.5,40.);
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
  mTau0 = new RooRealVar("mTau0","#tau_{0}",-6.5,-2.);
  mTau1 = new RooRealVar("mTau1","#tau_{1}",-2.,-0.2);
  mKbkg = new RooRealVar("mKbkg","mKbkg",0.,0.,1.);
  mKbkg->setConstant(false);
  mBkg0 = new RooExponential("mBkg0","background1",*mX,*mTau0);
  mBkg1 = new RooExponential("mBkg1","background2",*mX,*mTau1);
  mBackground = new RooAddPdf("mBackground","Background",RooArgList(*mBkg0,*mBkg1),RooArgList(*mKbkg));
  mMu->setConstant(false);
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
  mTau0 = new RooRealVar("mTau0","#tau_{0}",-6.5,-0.7);
  mTau1 = new RooRealVar("mTau1","#tau_{1}",-0.7,-0.01);
  mKbkg = new RooRealVar("mKbkg","mKbkg",0.,0.,1.);
  mKbkg->setConstant(false);
  mBkg0 = new RooExponential("mBkg0","background1",*mX,*mTau0);
  mBkg1 = new RooExponential("mBkg1","background2",*mX,*mTau1);
  mBackground = new RooAddPdf("mBackground","Background",RooArgList(*mBkg0,*mBkg1),RooArgList(*mKbkg));
  mAlpha0 = new RooRealVar("mAlpha0","Alpha0",1.6,3.); // tight range based on low pT
  mAlpha0->setConstant(false);
  mMu->setConstant(false);
  mSignal = new RooGausExp("mSignal","Signal",*mX,*mMu,*mSigma,*mAlpha0);
  mTemplate = new RooAddPdf("mTemplate","Template",RooArgList(*mSignal,*mBackground),
                            RooArgList(*mSigCounts,*mBkgCounts));
}

void FitExpExpTailGaus::SetHighPt(bool flag) {
  FitModule::SetHighPt(flag);
  /*if (flag) {
    mAlpha0->setVal(-5.);
    mAlpha0->setConstant(true);
  } else {
    mAlpha0->setConstant(false);
    mAlpha0->setRange(-5.,-1.5);
    mAlpha0->setVal(-1.8);
  }*/
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
  delete mBkgCounts;
  delete mSigCounts;
  delete mTemplate;
}

FitExpExpTailTailGaus::FitExpExpTailTailGaus(RooRealVar *x) : FitModule(x) {
  mTau0 = new RooRealVar("mTau0","#tau_{0}",-6.5,-2.);
  mTau1 = new RooRealVar("mTau1","#tau_{1}",-2.,-0.001);
  mKbkg = new RooRealVar("mKbkg","mKbkg",0.,0.,1.);
  mKbkg->setConstant(false);
  mBkg0 = new RooExponential("mBkg0","background1",*mX,*mTau0);
  mBkg1 = new RooExponential("mBkg1","background2",*mX,*mTau1);
  mBackground = new RooAddPdf("mBackground","Background",RooArgList(*mBkg0,*mBkg1),RooArgList(*mKbkg));
  mAlpha0 = new RooRealVar("mAlpha0","Alpha0",-3.,1.); // tight range based on low pT
  mAlpha1 = new RooRealVar("mAlpha1","Alpha1",1.,3.); // tight range based on low pT
  mSignal = new RooGausDExp("mSignal","Signal",*mX,*mMu,*mSigma,*mAlpha0,*mAlpha1);
  mTemplate = new RooAddPdf("mTemplate","Template",RooArgList(*mSignal,*mBackground),
                            RooArgList(*mSigCounts,*mBkgCounts));
}

void FitExpExpTailTailGaus::SetHighPt(bool flag) {
  FitModule::SetHighPt(flag);
}


FitExpExpTailTailGaus::~FitExpExpTailTailGaus() {
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

FitExpPolTailGaus::FitExpPolTailGaus(RooRealVar *x) : FitModule(x) {
  mTau0 = new RooRealVar("mTau0","#tau_{0}",-6.5,-2.);
  mA0 = new RooRealVar("mA0","#a_{0}",-0.15,-2.,-0.1);
  mA1 = new RooRealVar("mA1","#a_{1}",-2.,-0.1);
  mKbkg = new RooRealVar("mKbkg","mKbkg",0.,0.,1.);
  mKbkg->setConstant(false);
  mBkg0 = new RooExponential("mBkg0","background1",*mX,*mTau0);
  mBkg1 = new RooChebychev("mBkg1","background2",*mX,RooArgList(*mA0));
  mBackground = new RooAddPdf("mBackground","Background",RooArgList(*mBkg0,*mBkg1),RooArgList(*mKbkg));
  mAlpha0 = new RooRealVar("mAlpha0","Alpha0",1.6,3.); // tight range based on low pT
  mAlpha0->setConstant(false);
  mMu->setConstant(false);
  mSignal = new RooGausExp("mSignal","Signal",*mX,*mMu,*mSigma,*mAlpha0);
  mTemplate = new RooAddPdf("mTemplate","Template",RooArgList(*mSignal,*mBackground),RooArgList(*mSigCounts,*mBkgCounts));
}

void FitExpPolTailGaus::SetHighPt(bool flag) {
  FitModule::SetHighPt(flag);
}


FitExpPolTailGaus::~FitExpPolTailGaus() {
  delete mBackground;
  delete mTau0;
  delete mA0;
  delete mKbkg;
  delete mMu;
  delete mSigma;
  delete mSignal;
  //delete mAlpha0;
  delete mBkgCounts;
  delete mSigCounts;
  delete mTemplate;
}
