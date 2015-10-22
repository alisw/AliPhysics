/**************************************************************************
 * Copyright(c) 2004, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
//__________________________________________________________________
//
// Yves?
// What 
// is 
// this 
// class 
// supposed 
// to
// do?
//__________________________________________________________________
//
// --- ROOT system ---
#include <TClass.h>
#include <TH1F.h> 
#include <TF1.h> 
#include <TH2.h> 
#include <TH1I.h> 
#include <TIterator.h> 
#include <TKey.h> 
#include <TFile.h> 
#include <iostream>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TFitResult.h>
#include <TParameter.h>
#include <TMacro.h>
#include <TPaveText.h>
#include <TVirtualFitter.h>

// --- AliRoot header files ---
#include "AliLog.h"
#include "AliQAv1.h"
#include "AliQAChecker.h"
#include "AliFMDQAChecker.h"
#include "AliFMDQADataMakerRec.h"
#include "AliRecoParam.h"
#include <AliCDBManager.h>
#include <AliCDBEntry.h>
#include <AliCDBId.h>
#include <AliQAThresholds.h>

ClassImp(AliFMDQAChecker)
#if 0
; // This is for Emacs! - do not delete
#endif

namespace {
  void addFitsMacro(TList* l) {
    TMacro* m = new TMacro("fits");
    m->AddLine("void fits() {");
    m->AddLine("  if (!gPad) { Printf(\"No gPad\"); return; }");
    m->AddLine("  TList* lp = gPad->GetListOfPrimitives();");
    m->AddLine("  if (!lp) return;");
    m->AddLine("  TObject* po = 0;");
    m->AddLine("  TIter next(lp);");
    m->AddLine("  while ((po = next())) {");
    m->AddLine("    if (!po->IsA()->InheritsFrom(TH1::Class())) continue;");
    m->AddLine("    TH1*   htmp = dynamic_cast<TH1*>(po);");
    m->AddLine("    TList* lf   = htmp->GetListOfFunctions();");
    m->AddLine("    TObject* pso = (lf ? lf->FindObject(\"stats\") : 0);");
    m->AddLine("    if (!pso) continue;");
    m->AddLine("    TPaveStats* ps = static_cast<TPaveStats*>(pso);");
    m->AddLine("    ps->SetOptFit(111);");
    m->AddLine("    UShort_t qual = htmp->GetUniqueID();");
    m->AddLine("    ps->SetFillColor(qual >= 3 ? kRed-4 : qual >= 2 ? kOrange-4 : qual >= 1 ? kYellow-4 : kGreen-4);");
    // m->AddLine("    lf->Remove(lf->FindObject(\"fits\"));");
    // m->AddLine("    ps->Paint();");
    m->AddLine("    break;");
    m->AddLine("  }");
    // m->AddLine("  gPad->Modified(); gPad->Update(); gPad->cd();");
    m->AddLine("}");

    TObject* old = l->FindObject(m->GetName());
    if (old) l->Remove(old);
    l->Add(m);
  }
  
  const Double_t kROErrorsLabelY    = .30;
  
  const Int_t    kConvolutionSteps  = 100;
  const Double_t kConvolutionNSigma = 5;

  // 
  // The shift of the most probable value for the ROOT function TMath::Landau 
  //
  const Double_t  kMpShift  = -0.22278298;
  // 
  // Integration normalisation 
  //
  const Double_t  kInvSq2Pi = 1. / TMath::Sqrt(2*TMath::Pi());

  Double_t landau(Double_t x, Double_t delta, Double_t xi)
  {
    // 
    // Calculate the shifted Landau
    // @f[
    //    f'_{L}(x;\Delta,\xi) = f_L(x;\Delta+0.22278298\xi)
    // @f]
    //
    // where @f$ f_{L}@f$ is the ROOT implementation of the Landau
    // distribution (known to have @f$ \Delta_{p}=-0.22278298@f$ for
    // @f$\Delta=0,\xi=1@f$.
    // 
    // Parameters:
    //    x      Where to evaluate @f$ f'_{L}@f$ 
    //    delta  Most probable value 
    //    xi     The 'width' of the distribution 
    //
    // Return:
    //    @f$ f'_{L}(x;\Delta,\xi) @f$
    //
    return TMath::Landau(x, delta - xi * kMpShift, xi);
  }
  Double_t landauGaus(Double_t x, Double_t delta, Double_t xi,
		      Double_t sigma, Double_t sigmaN)
  {
    // 
    // Calculate the value of a Landau convolved with a Gaussian 
    // 
    // @f[ 
    // f(x;\Delta,\xi,\sigma') = \frac{1}{\sigma' \sqrt{2 \pi}}
    //    \int_{-\infty}^{+\infty} d\Delta' f'_{L}(x;\Delta',\xi)
    //    \exp{-\frac{(\Delta-\Delta')^2}{2\sigma'^2}}
    // @f]
    // 
    // where @f$ f'_{L}@f$ is the Landau distribution, @f$ \Delta@f$
    // the energy loss, @f$ \xi@f$ the width of the Landau, and @f$
    // \sigma'^2=\sigma^2-\sigma_n^2 @f$.  Here, @f$\sigma@f$ is the
    // variance of the Gaussian, and @f$\sigma_n@f$ is a parameter
    // modelling noise in the detector.
    //
    // Note that this function uses the constants kConvolutionSteps and
    // kConvolutionNSigma
    // 
    // References: 
    //  - <a href="http://dx.doi.org/10.1016/0168-583X(84)90472-5">Nucl.Instrum.Meth.B1:16</a>
    //  - <a href="http://dx.doi.org/10.1103/PhysRevA.28.615">Phys.Rev.A28:615</a>
    //  - <a href="http://root.cern.ch/root/htmldoc/tutorials/fit/langaus.C.html">ROOT implementation</a>
    // 
    // Parameters:
    //    x         where to evaluate @f$ f@f$
    //    delta     @f$ \Delta@f$ of @f$ f(x;\Delta,\xi,\sigma')@f$
    //    xi        @f$ \xi@f$ of @f$ f(x;\Delta,\xi,\sigma')@f$
    //    sigma     @f$ \sigma@f$ of @f$\sigma'^2=\sigma^2-\sigma_n^2 @f$
    //    sigma_n   @f$ \sigma_n@f$ of @f$\sigma'^2=\sigma^2-\sigma_n^2 @f$
    // 
    // Return:
    //    @f$ f@f$ evaluated at @f$ x@f$.  
    //
    Double_t deltap = delta - xi * kMpShift;
    Double_t sigma2 = sigmaN*sigmaN + sigma*sigma;
    Double_t sigma1 = sigmaN == 0 ? sigma : TMath::Sqrt(sigma2);
    Double_t xlow   = x - kConvolutionNSigma * sigma1;
    Double_t xhigh  = x + kConvolutionNSigma * sigma1;
    Double_t step   = (xhigh - xlow) / kConvolutionSteps;
    Double_t sum    = 0;
  
    for (Int_t i = 0; i <= kConvolutionSteps/2; i++) { 
      Double_t x1 = xlow  + (i - .5) * step;
      Double_t x2 = xhigh - (i - .5) * step;
    
      sum += TMath::Landau(x1, deltap, xi, kTRUE) * TMath::Gaus(x, x1, sigma1);
      sum += TMath::Landau(x2, deltap, xi, kTRUE) * TMath::Gaus(x, x2, sigma1);
    }
    return step * sum * kInvSq2Pi / sigma1;
  }

  // 
  // Utility function to use in TF1 defintition 
  //
  Double_t landauGaus1(Double_t* xp, Double_t* pp) 
  {
    Double_t x        = xp[0];
    Double_t constant = pp[0];
    Double_t delta    = pp[1];
    Double_t xi       = pp[2];
    Double_t sigma    = pp[3];
    Double_t sigmaN   = pp[4];

    return constant * landauGaus(x, delta, xi, sigma, sigmaN);
  }

  //____________________________________________________________________
  TF1* makeLandauGaus(const char* , 
		      Double_t  c=1, 
		      Double_t  delta=.5,  Double_t xi=0.07, 
		      Double_t  sigma=.1,  Double_t sigmaN=-1, 
		      Double_t  xmin=0,    Double_t xmax=15)
  {
    // 
    // Generate a TF1 object of @f$ f_I@f$ 
    // 
    // Parameters:
    //    c        Constant
    //    delta    @f$ \Delta@f$ 
    //    xi       @f$ \xi_1@f$	       
    //    sigma    @f$ \sigma_1@f$ 	       
    //    sigma_n  @f$ \sigma_n@f$ 	       
    //    xmin     Least value of range
    //    xmax     Largest value of range
    // 
    // Return:
    //    Newly allocated TF1 object
    //
    Int_t npar     = 5;
    TF1*  func     = new TF1("landauGaus", 
			     &landauGaus1,xmin,xmax,npar);
    // func->SetLineStyle(((i-2) % 10)+2); // start at dashed
    func->SetLineColor(kBlack); 
    func->SetLineWidth(2);
    func->SetNpx(500);
    func->SetParNames("C","#Delta_{p}","#xi", "#sigma", "#sigma_{n}");

    // Set the initial parameters from the seed fit 
    func->SetParameter(0,      c);       
    func->SetParameter(1,  delta);   
    func->SetParameter(2,     xi);      
    func->SetParameter(3,  sigma);   
    func->SetParameter(4, sigmaN); 

    func->SetParLimits(1, 0,    xmax);
    func->SetParLimits(2, 0,    xmax);
    func->SetParLimits(3, 0.01, 1);
    
    if (sigmaN < 0) func->FixParameter(4, 0);
    else            func->SetParLimits(4, 0, xmax);

    return func;
  }
}
  
//__________________________________________________________________
AliFMDQAChecker::AliFMDQAChecker() 
  : AliQACheckerBase("FMD","FMD Quality Assurance Checker") ,
    fDoScale(false),
    fDidExternal(false), 
    fShowFitResults(true),
    fELossLowCut(0.2), 
    fELossNRMS(3), 
    fELossBadChi2Nu(10), 
    fELossFkupChi2Nu(100), 
    fELossMinEntries(1000),
    fELossMaxEntries(-1),
    fELossGoodParError(0.1),
    fELossMinSharing(0.1),
    fROErrorsBad(0.3), 
    fROErrorsFkup(0.5),
    fMaxNProblem(10),
    fMaxNBad(10),
    fMinMPV(.2),
    fMaxXi(.3),
    fMaxSigma(.6),
    fNoFits(false)
{
}          

//__________________________________________________________________
void
AliFMDQAChecker::ProcessExternalParams()
{
  ProcessExternalParam("ELossLowCut",		fELossLowCut);
  ProcessExternalParam("ELossNRMS",		fELossNRMS);
  ProcessExternalParam("ELossBadChi2Nu",	fELossBadChi2Nu);
  ProcessExternalParam("ELossFkupChi2Nu",	fELossFkupChi2Nu);
  ProcessExternalParam("ELossGoodParError",	fELossGoodParError);
  ProcessExternalParam("ROErrorsBad",		fROErrorsBad);
  ProcessExternalParam("ROErrorsFkup",		fROErrorsFkup);
  ProcessExternalParam("ELossMinSharing",       fELossMinSharing);
  Double_t tmp = 0;
  ProcessExternalParam("CommonScale", tmp);
  fDoScale = tmp > 0; tmp = fELossMinEntries;
  ProcessExternalParam("ELossMinEntries", tmp);
  fELossMinEntries = tmp; tmp = fELossMaxEntries;
  ProcessExternalParam("ELossMaxEntries", tmp);
  fELossMaxEntries = tmp; tmp = fMaxNProblem;
  ProcessExternalParam("MaxNProblem", tmp);
  fMaxNProblem = tmp; tmp = 0;
  fELossMaxEntries = tmp; tmp = fMaxNBad;
  ProcessExternalParam("MaxNBad", tmp);
  fMaxNBad = tmp; tmp = 0;
  ProcessExternalParam("NoFits", tmp);
  fNoFits = tmp > 0; tmp = 0;

  ProcessExternalParam("ELossMinMPV",   fMinMPV);
  ProcessExternalParam("ELossMaxXi",    fMaxXi);
  ProcessExternalParam("ELossMaxSigma", fMaxSigma);
  
  GetThresholds();

  fDidExternal = true;
}
//__________________________________________________________________
void
AliFMDQAChecker::ProcessExternalParam(const char* name, Double_t& v)
{
  TObject* o = fExternParamList->FindObject(name);
  if (!o) return; 
  TParameter<double>* p = static_cast<TParameter<double>*>(o);
  v = p->GetVal();
  AliDebugF(3, "External parameter: %-20s=%lf", name, v);
}

//__________________________________________________________________
void
AliFMDQAChecker::GetThresholds()
{
  const char*    path   = "GRP/Calib/QAThresholds";
  AliCDBManager* cdbMan = AliCDBManager::Instance();
  AliCDBEntry*   cdbEnt = cdbMan->Get(path);
  if (!cdbEnt) { 
    AliWarningF("Failed to get CDB entry at %s", path);
    return;
  }
  
  TObjArray*     cdbObj = static_cast<TObjArray*>(cdbEnt->GetObject());
  if (!cdbObj) { 
    AliWarningF("Failed to get CDB object at %s", path);
    return;
  }
  
  TObject*       fmdObj = cdbObj->FindObject("FMD");
  if (!fmdObj) { 
    AliWarningF("Failed to get FMD object at from CDB %s", path);
    return;
  }
  
  AliQAThresholds* qaThr = static_cast<AliQAThresholds*>(fmdObj);
  Int_t nThr = qaThr->GetSize();
  for (Int_t i = 0; i < nThr; i++) { 
    TObject* thr = qaThr->GetThreshold(i);
    if (!thr) continue; 

    TParameter<double>* d = dynamic_cast<TParameter<double>*>(thr);
    if (!d) { 
      AliWarningF("Parameter %s not of type double", thr->GetName());
      continue;
    }
    Double_t val = d->GetVal();
    TString  name(thr->GetName());
    if      (name.EqualTo("ELossBadChi2Nu"))     fELossBadChi2Nu    = val;
    else if (name.EqualTo("ELossFkupChi2Nu"))    fELossFkupChi2Nu   = val;
    else if (name.EqualTo("ELossGoodParError"))  fELossGoodParError = val;
    else if (name.EqualTo("ROErrorsBad"))        fROErrorsBad       = val;
    else if (name.EqualTo("ROErrorsFkup"))       fROErrorsFkup      = val;    
    else if (name.EqualTo("MaxNProblem"))        fMaxNProblem       = val;
    else if (name.EqualTo("MaxNBad"))            fMaxNBad           = val;
    else if (name.EqualTo("ELossMinMPV"))	 fMinMPV	    = val;
    else if (name.EqualTo("ELossMaxXi"))	 fMaxXi		    = val;
    else if (name.EqualTo("ELossMaxSigma"))	 fMaxSigma	    = val;

    AliDebugF(3, "Threshold %s=%f", name.Data(), val);
  }
}

//__________________________________________________________________
AliQAv1::QABIT_t
AliFMDQAChecker::Quality2Bit(UShort_t value) const
{
  AliQAv1::QABIT_t  ret   = AliQAv1::kINFO; // Assume success 
  if      (value >= kWhatTheFk) ret = AliQAv1::kFATAL;
  else if (value >= kBad)       ret = AliQAv1::kERROR;
  else if (value >= kProblem)   ret = AliQAv1::kWARNING;
  
  return ret;
}

//__________________________________________________________________
void
AliFMDQAChecker::SetQA(AliQAv1::ALITASK_t index, Double_t* values) const
{
  AliQAv1 * qa = AliQAv1::Instance(index) ;

  for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) {
    // Check if specie is defined 
    if (!qa->IsEventSpecieSet(AliRecoParam::ConvertIndex(specie)))
      continue ;
    
    // No checker is implemented, set all QA to Fatal
    if (!values) { 
      qa->Set(AliQAv1::kFATAL, specie) ; 
      continue;
    }
    
    UShort_t          value = values[specie];
    AliQAv1::QABIT_t  ret   = Quality2Bit(value);
    
    qa->Set(ret, AliRecoParam::ConvertIndex(specie));
    AliDebugF(3, "Quality of %s: %d -> %d", 
	      AliRecoParam::GetEventSpecieName(specie), value, ret);
  }
}

//__________________________________________________________________
UShort_t
AliFMDQAChecker::BasicCheck(TH1* hist) const
{
  if (hist->GetEntries() <= 0) return kOK;
  return (hist->GetMean() > 0 ? kOK : kProblem);
}

//__________________________________________________________________
UShort_t
AliFMDQAChecker::CheckOne(AliQAv1::ALITASK_t          what,
			  AliRecoParam::EventSpecie_t specie, 
			  TH1*                        hist) const
{
  if(what == AliQAv1::kESD) return CheckESD(specie, hist);
  if(what == AliQAv1::kRAW) return CheckRaw(specie, hist);
  if(what == AliQAv1::kSIM) return CheckSim(specie, hist);
  if(what == AliQAv1::kREC) return CheckRec(specie, hist);
  return 0;
}
//__________________________________________________________________
UShort_t 
AliFMDQAChecker::CheckESD(AliRecoParam::EventSpecie_t /* specie*/, 
			  TH1*                        hist) const
{
  return BasicCheck(hist);
}
namespace {
  Double_t Chi2Scale(TH1* h, Double_t base=10000) 
  {
    return 1. / TMath::Max(1., h->GetEntries() / base);
  }
  void AddLine(TObjArray* lines, 
	       Double_t x1, Double_t x2, Double_t x3, 
	       Double_t y, Double_t dy,
	       const char* name, Double_t val, Double_t lim, 
	       Bool_t ok, Int_t color)
  {
    TString n; n.Form("%s:", name);
    TLatex* ltx = new TLatex(x1, y, n);
    ltx->SetNDC(true);
    ltx->SetTextSize(dy-0.01);
    ltx->SetTextColor(color);
    lines->Add(ltx);
    
    n.Form("%7.3f", val);
    ltx = new TLatex(x2, y, n);
    ltx->SetNDC(true);
    ltx->SetTextSize(dy-0.01);
    ltx->SetTextColor(color);
    lines->Add(ltx);
    
    if (lim < 0) n = "(ignored)";
    else  n.Form("%c %4.2f", ok ? '<' : '>', lim);
    ltx = new TLatex(x3, y, n);
    ltx->SetNDC(true);
    ltx->SetTextSize(dy-0.01);
    ltx->SetTextColor(color);
    lines->Add(ltx);
  }
}

//__________________________________________________________________
UShort_t
AliFMDQAChecker::CheckFit(TH1* hist, const TFitResultPtr& res, 
			  Double_t low, Double_t high, Int_t& color) const
{
  color = kGreen+4;

  // Check if there's indeed a result - if not, flag as OK
  if (!res.Get()) return 0;

  UShort_t   ret   = 0;
  Int_t      nPar  = res->NPar();
  Double_t   dy    = .06;
  Double_t   x     = .2;
  Double_t   x2    = .3;
  Double_t   x3    = .4;
  Double_t   y     = .9-dy;
  Double_t   chi2  = res->Chi2();
  Int_t      nu    = res->Ndf();
  Double_t   s     = Chi2Scale(hist,fELossMinEntries);
  Double_t   red   = (nu == 0 ? fELossFkupChi2Nu : chi2 / nu);
  TObjArray* lines = 0;
  // TLatex*    lRed  = 0;
  TLatex*    ltx   = 0;
  Int_t      chi2Check = 0;
  Double_t   chi2Lim   = fELossBadChi2Nu;
  if (AliDebugLevel() > 0) 
    printf("FIT: %s, 1, %d, %f, %f\n", hist->GetName(),
	    Int_t(hist->GetEntries()), red, s * red);
  red *= s;
  if (red > fELossBadChi2Nu) { // || res->Prob() < .01) { 
    // AliWarningF("Fit gave chi^2/nu=%f/%d=%f>%f (%f)", 
    //             res->Chi2(), res->Ndf(), red, fELossBadChi2Nu, 
    //             fELossFkupChi2Nu);
    // res->Print();
    chi2Check++;
    if (red > fELossFkupChi2Nu) { 
      chi2Check++;
      chi2Lim = fELossFkupChi2Nu;
    }
  }
  ret += chi2Check;

  if (fShowFitResults) { 
    lines = new TObjArray(nPar+3);
    lines->SetName("lines");
    lines->SetOwner(true);

    AddLine(lines, x, x2, x3, y, dy, "#chi^{2}/#nu", red, chi2Lim,
	    chi2Check < 1, chi2Check < 1 ? color : 
	    chi2Check < 2 ? kOrange+2 : kRed+2);
    
    Double_t x1 = .85;
    Double_t y1 = .5;

    // y1 -= dy;
    ltx = new TLatex(x1, y1, Form("Fit range: [%6.2f,%6.2f]", low, high));
    ltx->SetTextColor(kGray+3);
    ltx->SetTextSize(dy-.01);
    ltx->SetTextAlign(31);
    ltx->SetNDC(true);
    lines->Add(ltx);

    y1 -= dy;
    ltx = new TLatex(x1, y1, Form("Entries: %d (%d)", 
				  Int_t(hist->GetEffectiveEntries()),
				  fELossMaxEntries));
    ltx->SetTextColor(kGray+3);
    ltx->SetTextSize(dy-.01);
    ltx->SetTextAlign(31);
    ltx->SetNDC(true);
    lines->Add(ltx);

    y1 -= dy;
    ltx = new TLatex(x1, y1, Form("%s: %f #pm %f", 
				  res->ParName(1).c_str(),
				  res->Parameter(1),
				  res->ParError(1)));
    ltx->SetTextColor(kGray+3);
    ltx->SetTextSize(dy-.01);
    ltx->SetTextAlign(31);
    ltx->SetNDC(true);
    lines->Add(ltx);
  }
  
  // Now check the relative error on the fit parameters 
  Int_t parsOk = 0;
  for (Int_t i = 0; i < nPar; i++) { 
    if (res->IsParameterFixed(i)) continue; 
    Double_t thr = fELossGoodParError;
    Double_t pv  = res->Parameter(i);
    Double_t pe  = res->ParError(i);
    Double_t rel = (pv == 0 ? 100 : pe / pv);
    Bool_t   ok  = (i == 3) || (rel < thr);
    if (i == 1 && pv < fMinMPV) { // Low peak
      ok = false; ret++; }   // Double penalty 
    if (i == 2 && pv > fMaxXi) { // Large xi
      ok = false; ret++; }   // Double penalty 
    if (i == 3 && pv > fMaxSigma) { // Large sigma
      ok = false; ret++; }   // Double penalty 
    if (lines) {
      y -= dy;
      AddLine(lines, x, x2, x3, y, dy,Form("#delta%s/%s", 
					   res->ParName(i).c_str(),
					   res->ParName(i).c_str()),
	      rel, (i == 3 ? -1 : thr), ok, ok ? color : kOrange+2);
    }
    if (i == 3) continue; // Skip sigma 
    if (ok) parsOk++;
  }
  if (parsOk > 0) 
    ret = TMath::Max(ret-(parsOk-1),0);
  if (ret > 1) color = kRed+2;
  if (ret > 0) color = kOrange+2;

  if (lines) {
    TList*   lf  = hist->GetListOfFunctions();
    TObject* old = lf->FindObject(lines->GetName());
    if (old) {
      lf->Remove(old);
      delete old;
    }
    lf->Add(lines);
  }
  hist->SetStats(false);
    
  return ret;
}

//__________________________________________________________________
UShort_t
AliFMDQAChecker::CheckRaw(AliRecoParam::EventSpecie_t specie, 
			  TH1*                        hist) const
{
  Int_t ret = BasicCheck(hist);
  TString name(hist->GetName());
  if (name.Contains("readouterrors", TString::kIgnoreCase)) { 
    // Check the mean number of errors per event 
    TH2*  roErrors = static_cast<TH2*>(hist);
    Int_t nY       = roErrors->GetNbinsY();

    TLatex* ltx = new TLatex(.15, .9, Form("Thresholds: %5.2f,%5.2f",
					   fROErrorsBad, fROErrorsFkup));
    ltx->SetName("thresholds");
    ltx->SetTextColor(kGray+3);
    ltx->SetNDC();

    TList* ll = hist->GetListOfFunctions(); 
    TObject* old = ll->FindObject(ltx->GetName());
    if (old) { 
      ll->Remove(old);
      delete old;
    }
    ll->Add(ltx);
    
    for (Int_t i = 1; i <= 3; i++) {
      Double_t sum = 0;
      Int_t    cnt = 0;
      for (Int_t j = 1; j <= nY; j++) {
	Int_t n =  roErrors->GetBinContent(i, j);
	sum     += n * roErrors->GetYaxis()->GetBinCenter(j);
	cnt     += n;
      }
      Double_t mean = (cnt <= 0 ? 0 : sum / cnt);
      Double_t x    = ((i-.5) * (1-0.1-0.1) / 3 + 0.1);
      
      ltx = new TLatex(x, kROErrorsLabelY, Form("Mean: %6.3f", mean));
      ltx->SetName(Form("FMD%d", i));
      ltx->SetNDC();
      ltx->SetTextAngle(90);
      ltx->SetTextColor(kGreen+4);
      old = ll->FindObject(ltx->GetName());
      if (old) { 
	ll->Remove(old);
	delete old;
      }
      ll->Add(ltx);

      if (mean > fROErrorsBad) { 
	AliWarningF("Mean of readout errors for FMD%d = %f > %f (%f)", 
		    i, mean, fROErrorsBad, fROErrorsFkup);
	ret++;
	ltx->SetTextColor(kOrange+2);
	if (mean > fROErrorsFkup) {
	  ret++;
	  ltx->SetTextColor(kRed+2);
	}
      }
    }
  }
  else if (name.Contains("eloss",TString::kIgnoreCase)) { 
    // If we' asked to not fit the data, return immediately
    if (fNoFits) return ret;
    // Do not fit cosmic or calibration data 
    if (specie == AliRecoParam::kCosmic || 
	specie == AliRecoParam::kCalib) return ret;
    // Do not fit `expert' histograms 
    if (hist->TestBit(AliQAv1::GetExpertBit())) return ret;
    // Do not fit histograms with too little data 
    if (hist->GetEntries() < fELossMinEntries) return ret;

    // Try to fit a function to the histogram 
    Double_t xMin  = hist->GetXaxis()->GetXmin();
    Double_t xMax  = hist->GetXaxis()->GetXmax();

    hist->GetXaxis()->SetRangeUser(fELossLowCut, xMax);
    Int_t    bMaxY = hist->GetMaximumBin();
    Double_t xMaxY = hist->GetXaxis()->GetBinCenter(bMaxY);
    Double_t rms   = hist->GetRMS();
    Double_t low   = hist->GetXaxis()->GetBinCenter(bMaxY-4);
    hist->GetXaxis()->SetRangeUser(0.2, xMaxY+(fELossNRMS+1)*rms);
    rms  = hist->GetRMS();
    hist->GetXaxis()->SetRange(0,-1);
    TF1*          func = makeLandauGaus(name);
    func->SetParameter(1, xMaxY);
    func->SetLineColor(kGreen+4);
    // func->SetLineStyle(2);
    Double_t high = xMax; // xMaxY+fELossNRMS*rms;
    if (fELossNRMS > 0) high = xMaxY+fELossNRMS*rms;
    
    // Check we don't have an empty fit range 
    if (low >= high) return ret;

    // Check that we have enough counts in the fit range 
    Int_t bLow  = hist->FindBin(low);
    Int_t bHigh = hist->FindBin(high);
    if (bLow >= bHigh || hist->Integral(bLow, bHigh) < fELossMinEntries)
      return ret;

    // Set our fit function 
    TString fitOpt("QS");
    TFitResultPtr res   = hist->Fit(func, fitOpt, "", low, high);
    Int_t         color = func->GetLineColor();
    UShort_t      qual  = CheckFit(hist, res, low, high, color);

    // Make sure we save the function in the full range of the histogram
    func = hist->GetFunction("landauGaus");
    if (fELossNRMS <= 0) func->SetRange(xMin, xMax);
    // func->SetParent(hist);
    func->Save(xMin, xMax, 0, 0, 0, 0);
    func->SetLineColor(color);

    fitOpt.Append("+");
    res = hist->Fit("pol2", fitOpt, "", fELossMinSharing, low-0.05);
    func = hist->GetFunction("pol2");
    Double_t   s     = Chi2Scale(hist,fELossMinEntries*100);
    Double_t   chi2  = (!res.Get() ? 0 : res->Chi2());
    Int_t      nu    = (!res.Get() ? 1 : res->Ndf());
    Double_t   red   = s * (nu == 0 ? fELossFkupChi2Nu : chi2 / nu);
    if (AliDebugLevel()) 
      printf("FIT: %s, 2, %d, %f, %f\n", hist->GetName(),
	     Int_t(hist->GetEntries()), red, s * red);
    red *= s;
    if (red > fELossFkupChi2Nu) func->SetLineColor(kRed);
    else                        func->SetLineColor(kGreen+4);

    // Now check if this histogram should be cleared or not 
    if (fELossMaxEntries > 0 && hist->GetEntries() > fELossMaxEntries)
      hist->SetBit(AliFMDQADataMakerRec::kResetBit);
    if (qual > 0) { 
      func->SetLineWidth(3);
      func->SetLineStyle(1);
      if (qual > 1) 
	func->SetLineWidth(4);
    }
    ret += qual;
  }

  return ret;
}
//__________________________________________________________________
UShort_t
AliFMDQAChecker::CheckSim(AliRecoParam::EventSpecie_t /* specie*/, 
			  TH1*                        hist) const
{
  // 
  // Check simulated hits 
  // 
  return BasicCheck(hist);
}
//__________________________________________________________________
UShort_t
AliFMDQAChecker::CheckRec(AliRecoParam::EventSpecie_t /* specie*/, 
			  TH1*                        hist) const
{
  // 
  // Check reconstructed data 
  // 
  return BasicCheck(hist);
}

//__________________________________________________________________
void AliFMDQAChecker::AddStatusPave(TH1* hist, Int_t qual, 
				    Double_t xl, Double_t yl, 
				    Double_t xh, Double_t yh) const
{
  //
  // Add a status pave to a plot
  // 
  if (xh < 0) xh = gStyle->GetStatX();
  if (xl < 0) xl = xh - gStyle->GetStatW(); 
  if (yh < 0) yh = gStyle->GetStatY();
  if (yl < 0) yl = xl - gStyle->GetStatH(); 
  
  TPaveText* text = new TPaveText(xl, yl, xh, yh, "brNDC");
  Int_t   bg  = kGreen-10;
  Int_t   fg  = kBlack;
  TString msg = "OK";
  if      (qual >= kWhatTheFk) { bg = kRed+1; fg = kWhite; msg = "Argh!"; }
  else if (qual >= kBad)       { bg = kRed-3; fg = kWhite; msg = "Bad"; }
  else if (qual >= kProblem)   { bg = kOrange-4; msg = "Warning"; }
  text->AddText(msg);
  text->SetTextFont(62);
  text->SetTextColor(fg);
  text->SetFillColor(bg);

  TList*   ll  = hist->GetListOfFunctions();
  TObject* old = ll->FindObject(text->GetName());
  if (old) { 
    ll->Remove(old);
    delete old;
  }
  ll->Add(text);
}

//__________________________________________________________________
void AliFMDQAChecker::Check(Double_t*                   rv, 
			    AliQAv1::ALITASK_t          what, 
			    TObjArray**                 list, 
			    const AliDetectorRecoParam* /*t*/) 
{
  // 
  // Member function called to do the actual checking
  //
  // Parameters: 
  //    rv   Array of return values. 
  //    what What to check 
  //    list Array of arrays of histograms.  There's one arrat for
  //         each 'specie'
  //    t    Reconstruction parameters - not used. 
  //
  // The bounds defined for RV are 
  // 
  //   FATAL:      [-1,   0.0]
  //   ERROR:      (0.0,  0.002]
  //   WARNING:    (0.002,0.5]
  //   INFO:       (0.5,  1.0]
  //
  // Double_t* rv = new Double_t[AliRecoParam::kNSpecies] ; 

  if (!fDidExternal) ProcessExternalParams();

  for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) {
    // Int_t count   = 0;
    rv[specie]    = 0.; 

    if (!AliQAv1::Instance()->IsEventSpecieSet(specie) ) 
      continue ;
    
    if(!list[specie]) continue;
    
    TH1*     hist  = 0;
    Int_t    nHist = list[specie]->GetEntriesFast();

    // Find the status histogram if any 
    TH2*  status  = 0;
    Int_t istatus = AliFMDQADataMakerRec::GetHalfringIndex(4, 'i', 0, 0);
    if (istatus < nHist) 
      status = dynamic_cast<TH2*>(list[specie]->At(istatus));
      
    UShort_t ret   = 0;
    for(Int_t i= 0; i< nHist; i++) {
      if (!(hist = static_cast<TH1*>(list[specie]->At(i)))) continue;
      if (hist == status) continue;
      
      Int_t qual = CheckOne(what, AliRecoParam::ConvertIndex(specie), hist);
      hist->SetUniqueID(Quality2Bit(qual));
      hist->SetStats(0);
      AddStatusPave(hist, qual);
      ret += qual;

      if (!status) continue;

      // Parse out the detector and ring, calculate the bin, and fill
      // status histogram.
      TString nme(hist->GetName());
      Char_t cD   = nme[nme.Length()-2];
      Char_t cR   = nme[nme.Length()-1];
      Int_t  xbin = 0;
      switch (cD) { 
      case '1': xbin = 1; break;
      case '2': xbin = 2 + ((cR == 'i' || cR == 'I') ? 0 : 1); break;
      case '3': xbin = 4 + ((cR == 'i' || cR == 'I') ? 0 : 1); break;
      }
      if (xbin == 0) continue;
      status->Fill(xbin, qual);
		   
    } // for (int i ...)
    rv[specie] = ret;
    // if      (ret > kWhatTheFk) rv[specie] = fLowTestValue[AliQAv1::kFATAL];
    // else if (ret > kBad)       rv[specie] = fUpTestValue[AliQAv1::kERROR]; 
    // else if (ret > kProblem)   rv[specie] = fUpTestValue[AliQAv1::kWARNING]; 
    // else                       rv[specie] = fUpTestValue[AliQAv1::kINFO]; 
    AliDebugF(3, "Combined sum is %d -> %f", ret, rv[specie]);

    if (status) { 
      Int_t nProblem = 0;
      Int_t nBad     = 0;
      for (Int_t i = 1; i < status->GetXaxis()->GetNbins(); i++) { 
	nProblem += status->GetBinContent(i, 3);
	nBad     += status->GetBinContent(i, 4);
      }
      Int_t qual = 0;
      if (nProblem > fMaxNProblem) qual++;
      if (nBad     > fMaxNBad)     qual += 2;
      status->SetUniqueID(Quality2Bit(qual));
      AddStatusPave(status, qual);
    }
    // if (count != 0) rv[specie] /= count;
  }
  // return rv;
}

namespace {
  Int_t CheckForLog(TAxis*       axis,
		    TVirtualPad* pad, 
		    Int_t        xyz)
  {
    Int_t ret = 0;
    TString t(axis->GetTitle());
    if (!t.Contains("[log]", TString::kIgnoreCase)) return 0;
    t.ReplaceAll("[log]", "");
    switch (xyz) { 
    case 1: pad->SetLogx(); ret |= 0x1; break;
    case 2: pad->SetLogy(); ret |= 0x2; break;
    case 3: pad->SetLogz(); ret |= 0x4; break;
    }
    axis->SetTitle(t);
    return ret;
  }
  void RestoreLog(TAxis* axis, Bool_t log) 
  {
    if (!log) return;
    TString t(axis->GetTitle());
    t.Append("[log]");
    axis->SetTitle(t);
  }
}

namespace {
  void FindMinMax(TH1* h, Double_t& min, Double_t& max)
  {
    Double_t tmin = 1e9;
    Double_t tmax = 0;
    for (Int_t i = 1; i <= h->GetNbinsX(); i++) { 
      Double_t c = h->GetBinContent(i);
      if (c < 1e-8) continue;
      tmin = TMath::Min(tmin, c);
      tmax = TMath::Max(tmax, c);
    }
    min = tmin;
    max = tmax;
  }
}

namespace { 
  Int_t GetHalfringPad(TH1* h) {
    TString nme(h->GetName());
    Char_t cD   = nme[nme.Length()-2];
    Char_t cR   = nme[nme.Length()-1];
    Int_t  xbin = 0;
    switch (cD) { 
    case '1': xbin = 1; break;
    case '2': xbin = ((cR == 'i' || cR == 'I') ? 2 : 5); break;
    case '3': xbin = ((cR == 'i' || cR == 'I') ? 3 : 6); break;
    }
    return xbin;
  }
}

//____________________________________________________________________________ 
void 
AliFMDQAChecker::MakeImage(TObjArray** list, 
			   AliQAv1::TASKINDEX_t task, 
			   AliQAv1::MODE_t mode) 
{
  // makes the QA image for sim and rec
  // 
  // Parameters: 
  //    task What to check 
  //    list Array of arrays of histograms.  There's one array for
  //         each 'specie'
  //    t    Reconstruction parameters - not used. 
  // 
  Int_t    nImages = 0 ;
  Double_t max     = 0;
  Double_t min     = 10000;

  // Loop over all species 
  for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) {
    AliRecoParam::EventSpecie_t spe = AliRecoParam::ConvertIndex(specie);
    if (!AliQAv1::Instance(AliQAv1::GetDetIndex(GetName()))
	->IsEventSpecieSet(spe)) 
      continue;
									
    // If nothing is defined for this specie, go on. 
    if(!list[specie] || list[specie]->GetEntriesFast() == 0) continue;

    // Loop over the histograms and figure out how many histograms we
    // have and the min/max 
    TH1* hist  = 0;
    Int_t nHist = list[specie]->GetEntriesFast();
    for(Int_t i= 0; i< nHist; i++) {
      hist = static_cast<TH1F*>(list[specie]->At(i));
      if (hist && hist->TestBit(AliQAv1::GetImageBit())) {
        nImages++; 
	TString name(hist->GetName());
	if (name.Contains("readouterrors", TString::kIgnoreCase) || 
	    name.Contains("status", TString::kIgnoreCase)) continue;

	// Double_t hMax = hist->GetMaximum(); 
	// hist->GetBinContent(hist->GetMaximumBin());
	// Double_t hMin = hist->GetMinimum();
	// hist->GetBinContent(hist->GetMinimumBin());
	Double_t hMax, hMin;
	FindMinMax(hist, hMin, hMax);
	max = TMath::Max(max, hMax);
	min = TMath::Min(min, hMin);
	// AliInfoF("Min/max of %40s: %f/%f, global -> %f/%f", 
	//          hist->GetName(), hMin, hMax, min, max);
      }
    }
    break ; 
  }
  // AliInfoF("Global min/max=%f/%f", min, max);
  min = TMath::Max(1e-1, min);
  max = TMath::Max(1e5,  max);

  // IF no images, go on. 
  if (nImages == 0) {
    AliDebug(AliQAv1::GetQADebugLevel(), 
	     Form("No histogram will be plotted for %s %s\n", GetName(), 
		  AliQAv1::GetTaskName(task).Data()));
    return;
  }

  AliDebug(AliQAv1::GetQADebugLevel(), 
	   Form("%d histograms will be plotted for %s %s\n", 
		nImages, GetName(), AliQAv1::GetTaskName(task).Data()));  
  gStyle->SetOptStat(0);
  
  // Again loop over species and draw a canvas 
  for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) {
    if (!AliQAv1::Instance(AliQAv1::GetDetIndex(GetName()))
	->IsEventSpecieSet(specie)) continue;

    // if Nothing here, go on
    if(!list[specie] || list[specie]->GetEntries() <= 0 || 
       nImages <= 0) continue;

   // Form the title 
    const Char_t * title = Form("QA_%s_%s_%s", GetName(), 
				AliQAv1::GetTaskName(task).Data(), 
				AliRecoParam::GetEventSpecieName(specie)); 
    if (!fImage[specie]) fImage[specie] = new TCanvas(title, title) ;
    fImage[specie]->Clear() ; 
    fImage[specie]->SetTitle(title) ; 
    fImage[specie]->cd() ; 

    // Put something in the canvas - even if empty 
    TPaveText someText(0.015, 0.015, 0.98, 0.98) ;
    someText.AddText(title) ;
    someText.SetFillColor(0);
    someText.SetFillStyle(0);
    someText.SetBorderSize(0);
    someText.SetTextColor(kRed+1);
    someText.Draw() ; 
    TString outName(Form("%s%s%d.%s", AliQAv1::GetImageFileName(), 
			 AliQAv1::GetModeName(mode), 
			 AliQAChecker::Instance()->GetRunNumber(), 
			 AliQAv1::GetImageFileFormat()));
    fImage[specie]->Print(outName, "ps") ; 

    // Now set some parameters on the canvas 
    fImage[specie]->Clear(); 
    fImage[specie]->SetTopMargin(0.10);
    fImage[specie]->SetBottomMargin(0.15);
    fImage[specie]->SetLeftMargin(0.15);
    fImage[specie]->SetRightMargin(0.05);
    
    // Put title on top 
    const char* topT = Form("Mode: %s, Task: %s, Specie: %s, Run: %d",
			    AliQAv1::GetModeName(mode), 
			    AliQAv1::GetTaskName(task).Data(), 
			    AliRecoParam::GetEventSpecieName(specie),
			    AliQAChecker::Instance()->GetRunNumber());
    TLatex* topText = new TLatex(.5, .99, topT);
    topText->SetTextAlign(23);
    topText->SetTextSize(.038);
    topText->SetTextFont(42);
    topText->SetTextColor(kBlue+3);
    topText->SetNDC();
    topText->Draw();

    // Find the status histogram if any 
    TH2*  status  = 0;
    Int_t istatus = AliFMDQADataMakerRec::GetHalfringIndex(4, 'i', 0, 0);
    if (istatus < list[specie]->GetEntriesFast()) 
      status = dynamic_cast<TH2*>(list[specie]->At(istatus));

    // Divide canvas 
    // if (fDoScale) 
    TVirtualPad* plots = fImage[specie];
    TVirtualPad* stat  = 0;
    if (status) {
      // AliWarning("Drawing plots sub-pad");
      TPad* pM = new TPad("plots", "Plots Pad", 0, .2, 1., .9, 0, 0);
      fImage[specie]->cd();
      pM->Draw();
      plots = pM;
      // AliWarning("Drawing status sub-pad");
      TPad* pS = new TPad("status", "Status Pad", 0, 0, 1., .2, 0, 0);
      fImage[specie]->cd();
      pS->Draw();
      pS->SetLogz();
      stat = pS;
      // status->DrawCopy("colz");
    }
    // AliWarningF("fImage[specie]=%p, plots=%p", fImage[specie], plots);
    // plots->cd();
    Int_t nx = 3;
    Int_t ny = (nImages + .5) / nx;
    plots->Divide(nx, ny, 0, 0);
    // else fImage[specie]->Divide(nx, ny);
    
    
    // Loop over histograms 
    TH1*  hist  = 0;
    Int_t nHist = list[specie]->GetEntriesFast();
    for (Int_t i = 0; i < nHist; i++) { 
      hist = static_cast<TH1*>(list[specie]->At(i));
      if (!hist || !hist->TestBit(AliQAv1::GetImageBit())) continue;
      if (hist == status) continue;
      TString name(hist->GetName());
      Bool_t isROE = name.Contains("readouterrors", TString::kIgnoreCase);

      // Go to sub-pad 
      TVirtualPad* pad = 0;
      if      (isROE) pad = plots->cd(4);
      else            pad = plots->cd(GetHalfringPad(hist));
      
      pad->SetRightMargin(0.01);
      if (!fDoScale) { 
	pad->SetLeftMargin(0.10);
	pad->SetBottomMargin(0.10);
      }

      // Check for log scale 
      Int_t logOpts = 0;
      logOpts |= CheckForLog(hist->GetXaxis(), pad, 1);
      logOpts |= CheckForLog(hist->GetYaxis(), pad, 2);
      logOpts |= CheckForLog(hist->GetZaxis(), pad, 3);

      // Figure out special cases 
      TString opt("");
      if (isROE) {
	pad->SetRightMargin(0.15);
	pad->SetBottomMargin(0.10);
	// pad->SetTopMargin(0.02);
	opt="COLZ";
      }
      else {
	pad->SetGridx();
	pad->SetGridy();
	if (fDoScale) { 
	  hist->SetMinimum(min);
	  hist->SetMaximum(max);
	}
	else { 
	  hist->SetMinimum();
	  hist->SetMaximum();
	}
      }
      // Draw (As a copy)
      hist->DrawCopy(opt);
      
      // Special cases 
      if (!name.Contains("readouterrors", TString::kIgnoreCase)) {
	gStyle->SetOptTitle(0);
	TPad* insert = new TPad("insert", "Zoom", 
				.5,.5, .99, .95, 0, 0, 0);
	insert->SetTopMargin(0.01);
	insert->SetRightMargin(0.01);
	insert->SetFillColor(0);
	insert->SetBorderSize(1);
	insert->SetBorderMode(0);
	insert->Draw();
	insert->cd();
	if (logOpts & 0x1) insert->SetLogx();
	if (logOpts & 0x2) insert->SetLogy();
	if (logOpts & 0x4) insert->SetLogz();
	hist->GetXaxis()->SetRange(1, hist->GetNbinsX()/8);
	TH1* copy = hist->DrawCopy(opt);
	copy->GetXaxis()->SetNdivisions(408, false);
	// Restore full range 
	hist->GetXaxis()->SetRange(0, 0);
	gStyle->SetOptTitle(1);
      }
      pad->cd();
      // Possibly restore the log options 
      RestoreLog(hist->GetXaxis(), logOpts & 0x1);
      RestoreLog(hist->GetYaxis(), logOpts & 0x2);
      RestoreLog(hist->GetZaxis(), logOpts & 0x4);
    }
    if (status && stat) {
      stat->cd();
      status->DrawCopy("BOX TEXT");
    }
    // Print to a post-script file 
    fImage[specie]->Print(outName, "ps");
    if (AliDebugLevel() > 0) 
      fImage[specie]->Print(Form("%s_%d.png", 
				 AliRecoParam::GetEventSpecieName(specie), 
				 AliQAChecker::Instance()->GetRunNumber()));
  }
}

//__________________________________________________________________
//
// EOF
//
