#include "FitModules.h"
#include "Common.h"
#include "RooGausExp.h"
#include "RooGausDExp.h"
#include "Utils.h"
using namespace utils;

#include <TLegend.h>

#include <RooArgList.h>
#include <RooArgSet.h>
#include <RooCBShape.h>
#include <RooChebychev.h>
#include <RooCurve.h>
#include <RooDataHist.h>
#include <RooFitResult.h>
#include <RooGaussian.h>


RooPlot* FitModule::FitData(TH1* dat,TString name, TString title, TString range, TString plotrange, bool change_range, float low_x, float high_x) {
  if (plotrange == "") plotrange = range;
  RooDataHist data("data","data",RooArgList(*mX),Import(*dat));
  mNentries = data.sum(false);
  RooPlot* plot = mX->frame();
  plot->SetTitle(title.Data());
  plot->SetName(name.Data());
  plot->GetYaxis()->SetTitle(Form("Counts / (%.2f %s)",plot->GetXaxis()->GetBinWidth(1),mX->getUnit()));
  for (int i = 2; i--;){
    if(mExtended) RooFitResult *res = mTemplate->fitTo(data,Extended(),Verbose(kFALSE),PrintEvalErrors(-1),PrintLevel(-1),Range(range));
    else RooFitResult *res = mTemplate->fitTo(data,Verbose(kFALSE),PrintEvalErrors(-1),PrintLevel(-1),Range(range));
  }
  if(change_range) plot->GetXaxis()->SetRangeUser(low_x,high_x);
  data.plotOn(plot,Name("data"),DrawOption("pz"));
  mTemplate->plotOn(plot,Name("model"),Range(range),NormRange(range));
  mTemplate->plotOn(plot,Name("bkg"),Components(*mBackground),LineStyle(kDashed),LineColor(kRed),Range(range),NormRange(range));
  mTemplate->plotOn(plot,Name("sig"),Components(*mSignal),LineStyle(kDashed),LineColor(kGreen+3),Range(range),NormRange(range));
  mChi2 = plot->chiSquare("model","data");
  plot->remove("model",false);
  plot->remove("bkg",false);
  plot->remove("sig",false);
  mTemplate->plotOn(plot,Name("model"),Range(plotrange),NormRange(range));
  mTemplate->plotOn(plot,Name("bkg"),Components(*mBackground),LineStyle(kDashed),LineColor(kRed),Range(plotrange),NormRange(range));
  mTemplate->plotOn(plot,Name("sig"),Components(*mSignal),LineStyle(kDashed),LineColor(kGreen+3),Range(plotrange),NormRange(range));
  mTemplate->paramOn(plot,Label(Form("#chi^{2}/NDF = %2.4f",mChi2)),Layout(0.64,0.92,0.86));
  plot->getAttLine()->SetLineWidth(0);
  data.removeSelfFromDir();
  return plot;
}

void FitModule::UseBackground(bool useBkg){
  if(mExtended){
    mBkgCounts->setConstant(!useBkg);
    mBkgCounts->setVal(1000 * int(useBkg));
  }
  else{
    mFraction->setConstant(!useBkg);
    if(!useBkg) mFraction->setVal(1.);
  }
  auto iter = mBackground->getVariables()->createIterator();
  RooRealVar* var = (RooRealVar*)iter->Next();
  while(var){
    var->setConstant(!useBkg);
    var = (RooRealVar*)iter->Next();
  }
}

void FitModule::UseSignal(bool useSig){
  if(mExtended){
    mSigCounts->setConstant(!useSig);
    mSigCounts->setVal(1000 * int(useSig));
  }
  else{
    mFraction->setConstant(!useSig);
    if(!useSig) mFraction->setVal(0.);
  }
  auto iter = mSignal->getVariables()->createIterator();
  RooRealVar* var = (RooRealVar*)iter->Next();
  while(var){
    var->setConstant(!useSig);
    var = (RooRealVar*)iter->Next();
  }
}

FitGausGaus::FitGausGaus(RooRealVar *x, bool extended_likelihood) : FitModule(x, extended_likelihood){
  mSigmaBkg = make_unique<RooRealVar>("#sigma_{bkg}","#sigma_{Bkg}",0.01,2.);
  mMuBkg = make_unique<RooRealVar>("#mu_{bkg}","#mu_{Bkg}",-5.,-3.);
  mBackground = make_unique<RooGaussian>("mBackground","Background",*mX,*mMuBkg,*mSigmaBkg);
  mSignal = make_unique<RooGaussian>("mSignal","Signal",*mX,*mMu,*mSigma);
  mTemplate = (mExtended) ? make_unique<RooAddPdf>("mTemplate","Template",RooArgList(*mSignal,*mBackground),
                                     RooArgList(*mSigCounts,*mBkgCounts)) : make_unique<RooAddPdf>("mTemplate","Template",RooArgList(*mSignal,*mBackground),
                                                                        RooArgList(*mFraction));
}

FitExpTailGaus::FitExpTailGaus(RooRealVar *x, bool extended_likelihood) : FitModule(x, extended_likelihood){
  mTau0 = make_unique<RooRealVar>("#tau_{0}","tau",-10.,-0.00001);
  mBackground = make_unique<RooExponential>("mBackground","Background",*mX,*mTau0);
  mAlpha0 = make_unique<RooRealVar>("#alpha_{0}","Alpha0",1.6,3.);
  mSignal = make_unique<RooGausExp>("mSignal","Signal",*mX,*mMu,*mSigma,*mAlpha0);
  mTemplate = (mExtended) ? make_unique<RooAddPdf>("mTemplate","Template",RooArgList(*mSignal,*mBackground),
                                     RooArgList(*mSigCounts,*mBkgCounts)) : make_unique<RooAddPdf>("mTemplate","Template",RooArgList(*mSignal,*mBackground),
                                                                        RooArgList(*mFraction));
}

FitExpGaus::FitExpGaus(RooRealVar *x, bool extended_likelihood) : FitModule(x, extended_likelihood) {
  mTau = make_unique<RooRealVar>("mTau","tau bkg",-5.,0.);
  mBackground = make_unique<RooExponential>("mBackground","Background",*mX,*mTau);
  mSignal = make_unique<RooGaussian>("mSignal","Signal",*mX,*mMu,*mSigma);
  mTemplate = (mExtended) ? make_unique<RooAddPdf>("mTemplate","Template",RooArgList(*mSignal,*mBackground),
                                     RooArgList(*mSigCounts,*mBkgCounts)) : make_unique<RooAddPdf>("mTemplate","Template",RooArgList(*mSignal,*mBackground),
                                                                        RooArgList(*mFraction));
}

FitExpCB::FitExpCB(RooRealVar *x, bool extended_likelihood) : FitModule(x, extended_likelihood) {
  mTau = make_unique<RooRealVar>("mTau","tau bkg",-5.,0.);
  mBackground = make_unique<RooExponential>("mBackground","Background",*mX,*mTau);
  mAlpha = make_unique<RooRealVar>("mAlpha","Alpha",-4.,-1.75);
  mN = make_unique<RooRealVar>("mN","n",3.,10.);
  mSignal = make_unique<RooCBShape>("mSignal","Signal",*mX,*mMu,*mSigma,*mAlpha,*mN);
  mTemplate = (mExtended) ? make_unique<RooAddPdf>("mTemplate","Template",RooArgList(*mSignal,*mBackground),
                                     RooArgList(*mSigCounts,*mBkgCounts)) : make_unique<RooAddPdf>("mTemplate","Template",RooArgList(*mSignal,*mBackground),
                                                                        RooArgList(*mFraction));

}

FitExpExpCB::FitExpExpCB(RooRealVar *x, bool extended_likelihood) : FitModule(x, extended_likelihood) {
  mTau0 = make_unique<RooRealVar>("#tau_{0}","#tau_{0}",-6.5,-2.);
  mTau1 = make_unique<RooRealVar>("#tau_{1}","#tau_{1}",-2.,-0.2);
  mKbkg = make_unique<RooRealVar>("K_{bkg}","K_{bkg}",0.,0.,1.);
  mBkg0 = make_unique<RooExponential>("mBkg0","background1",*mX,*mTau0);
  mBkg1 = make_unique<RooExponential>("mBkg1","background2",*mX,*mTau1);
  mBackground = make_unique<RooAddPdf>("mBackground","Background",RooArgList(*mBkg0,*mBkg1),
                                       RooArgList(*mKbkg));
  mAlpha = make_unique<RooRealVar>("mAlpha","Alpha",-3.,-0.8);
  mN = make_unique<RooRealVar>("mN","n",3.5,40.);
  mSignal = make_unique<RooCBShape>("mSignal","Signal",*mX,*mMu,*mSigma,*mAlpha,*mN);
  mTemplate = (mExtended) ? make_unique<RooAddPdf>("mTemplate","Template",RooArgList(*mSignal,*mBackground),
                                     RooArgList(*mSigCounts,*mBkgCounts)) : make_unique<RooAddPdf>("mTemplate","Template",RooArgList(*mSignal,*mBackground),
                                                                        RooArgList(*mFraction));
}

FitExpExpGaus::FitExpExpGaus(RooRealVar *x, bool extended_likelihood) : FitModule(x, extended_likelihood) {
  mTau0 = make_unique<RooRealVar>("#tau_{0}","#tau_{0}",-6.5,-2.);
  mTau1 = make_unique<RooRealVar>("#tau_{1}","#tau_{1}",-2.,-0.2);
  mKbkg = make_unique<RooRealVar>("K_{bkg}","K_{bkg}",0.,0.,1.);
  mBkg0 = make_unique<RooExponential>("mBkg0","background1",*mX,*mTau0);
  mBkg1 = make_unique<RooExponential>("mBkg1","background2",*mX,*mTau1);
  mBackground = make_unique<RooAddPdf>("mBackground","Background",RooArgList(*mBkg0,*mBkg1),
                                       RooArgList(*mKbkg));
  mSignal = make_unique<RooGaussian>("mSignal","Signal",*mX,*mMu,*mSigma);
  mTemplate = (mExtended) ? make_unique<RooAddPdf>("mTemplate","Template",RooArgList(*mSignal,*mBackground),
                                     RooArgList(*mSigCounts,*mBkgCounts)) : make_unique<RooAddPdf>("mTemplate","Template",RooArgList(*mSignal,*mBackground),
                                                                        RooArgList(*mFraction));
}

FitExpExpTailGaus::FitExpExpTailGaus(RooRealVar *x, bool extended_likelihood) : FitModule(x, extended_likelihood) {
  mTau0 = make_unique<RooRealVar>("#tau_{0}","#tau_{0}",-10.,-0.5);
  mTau1 = make_unique<RooRealVar>("#tau_{1}","#tau_{1}",-0.5,-0.01);
  mKbkg = make_unique<RooRealVar>("K_{bkg}","K_{bkg}",0.5,0.,1.);
  mBkg0 = make_unique<RooExponential>("mBkg0","background1",*mX,*mTau0);
  mBkg1 = make_unique<RooExponential>("mBkg1","background2",*mX,*mTau1);
  mBackground = make_unique<RooAddPdf>("mBackground","Background",RooArgList(*mBkg0,*mBkg1),
                                       RooArgList(*mKbkg));
  mAlpha0 = make_unique<RooRealVar>("#alpha","Alpha0",1.6,3.); // tight range based on low pT
  mSignal = make_unique<RooGausExp>("mSignal","Signal",*mX,*mMu,*mSigma,*mAlpha0);
  mTemplate = (mExtended) ? make_unique<RooAddPdf>("mTemplate","Template",RooArgList(*mSignal,*mBackground),
                                     RooArgList(*mSigCounts,*mBkgCounts)) : make_unique<RooAddPdf>("mTemplate","Template",RooArgList(*mSignal,*mBackground),
                                                                        RooArgList(*mFraction));
}

FitExpExpTailTailGaus::FitExpExpTailTailGaus(RooRealVar *x, bool extended_likelihood) : FitModule(x, extended_likelihood) {
  mTau0 = make_unique<RooRealVar>("#tau_{0}","#tau_{0}",-6.5,-2.);
  mTau1 = make_unique<RooRealVar>("#tau_{1}","#tau_{1}",-2.,-0.001);
  mKbkg = make_unique<RooRealVar>("K_{bkg}","K_{bkg}",0.5,0.,1.);
  mBkg0 = make_unique<RooExponential>("mBkg0","background1",*mX,*mTau0);
  mBkg1 = make_unique<RooExponential>("mBkg1","background2",*mX,*mTau1);
  mBackground = make_unique<RooAddPdf>("mBackground","Background",RooArgList(*mBkg0,*mBkg1),
                                       RooArgList(*mKbkg));
  mAlpha0 = make_unique<RooRealVar>("#alpha_{0}","Alpha0",-3.,-1.); // tight range based on low pT
  mAlpha1 = make_unique<RooRealVar>("#alpha_{1}","Alpha1",1.,3.); // tight range based on low pT
  mSignal = make_unique<RooGausDExp>("mSignal","Signal",*mX,*mMu,*mSigma,*mAlpha0,*mAlpha1);
  mTemplate = (mExtended) ? make_unique<RooAddPdf>("mTemplate","Template",RooArgList(*mSignal,*mBackground),
                                     RooArgList(*mSigCounts,*mBkgCounts)) : make_unique<RooAddPdf>("mTemplate","Template",RooArgList(*mSignal,*mBackground),
                                                                        RooArgList(*mFraction));
}

FitExpPolTailGaus::FitExpPolTailGaus(RooRealVar *x, bool extended_likelihood) : FitModule(x, extended_likelihood) {
  mTau0 = make_unique<RooRealVar>("#tau_{0}","#tau_{0}",-6.5,-2.);
  mA0 = make_unique<RooRealVar>("mA0","#a_{0}",-0.15,-2.,-0.1);
  mA1 = make_unique<RooRealVar>("mA1","#a_{1}",-2.,-0.1);
  mKbkg = make_unique<RooRealVar>("K_{bkg}","K_{bkg}",0.,0.,1.);
  mBkg0 = make_unique<RooExponential>("mBkg0","background1",*mX,*mTau0);
  mBkg1 = make_unique<RooChebychev>("mBkg1","background2",*mX,RooArgList(*mA0));
  mBackground = make_unique<RooAddPdf>("mBackground","Background",RooArgList(*mBkg0,*mBkg1),
                                       RooArgList(*mKbkg));
  mAlpha0 = make_unique<RooRealVar>("#alpha_{0}","Alpha0",1.6,3.); // tight range based on low pT
  mSignal = make_unique<RooGausExp>("mSignal","Signal",*mX,*mMu,*mSigma,*mAlpha0);
  mTemplate = (mExtended) ? make_unique<RooAddPdf>("mTemplate","Template",RooArgList(*mSignal,*mBackground),
                                     RooArgList(*mSigCounts,*mBkgCounts)) : make_unique<RooAddPdf>("mTemplate","Template",RooArgList(*mSignal,*mBackground),
                                                                        RooArgList(*mFraction));
}

FitGausExpTailGaus::FitGausExpTailGaus(RooRealVar *x, bool extended_likelihood) : FitModule(x, extended_likelihood) {
  mSigmaBkg = make_unique<RooRealVar>("#sigma_{bkg}","#sigma_{Bkg}",0.05,2.);
  mMuBkg = make_unique<RooRealVar>("#mu_{bkg}","#mu_{Bkg}",-5.,-1.);
  mKbkg = make_unique<RooRealVar>("K_{bkg}","K_{bkg}",0.5,0.,1.);
  mBkg0 = make_unique<RooGaussian>("mBkg0","background1",*mX,*mMuBkg,*mSigmaBkg);
  mTau1 = make_unique<RooRealVar>("#tau_{1}","#tau_{1}",-10,-0.01);
  mBkg1 = make_unique<RooExponential>("mBkg1","background2",*mX,*mTau1);
  mBackground = make_unique<RooAddPdf>("mBackground","Background",RooArgList(*mBkg0,*mBkg1),
                                       RooArgList(*mKbkg));
  mAlpha0 = make_unique<RooRealVar>("#alpha","Alpha0",1.6,3.); // tight range based on low pT
  mSignal = make_unique<RooGausExp>("mSignal","Signal",*mX,*mMu,*mSigma,*mAlpha0);
  mTemplate = (mExtended) ? make_unique<RooAddPdf>("mTemplate","Template",RooArgList(*mSignal,*mBackground),
                                     RooArgList(*mSigCounts,*mBkgCounts)) : make_unique<RooAddPdf>("mTemplate","Template",RooArgList(*mSignal,*mBackground),
                                                                        RooArgList(*mFraction));
}
