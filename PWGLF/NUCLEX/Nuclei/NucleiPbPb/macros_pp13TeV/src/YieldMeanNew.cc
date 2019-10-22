#if !defined (__CINT__) || (defined(__MAKECINT__))
#include <iostream>
#include "TH1.h"
#include "TVirtualFitter.h"
#include "TMath.h"
#include "TFile.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TRandom.h"
#include "TDirectory.h"

#include <memory>

#endif
using namespace std;
/* definition of the fields in the histogram returned */
enum EValue_t {
  kYield = 1,
  kYieldStat,
  kYieldSysUncorr,
  kYieldSysHiCorr,
  kYieldSysLoCorr,
  kMean,
  kMeanStat,
  kMeanSysUncorr,
  kMeanSysHiCorr,
  kMeanSysLoCorr,
  kFitRes
};

void YieldMean_IntegralMean(TH1 *hdata, TH1 *hlo, TH1 *hhi, Double_t &integral, Double_t &mean,Bool_t printinfo=kFALSE);
TH1 * YieldMean_LowExtrapolationHisto(TH1 *h, TF1 *f, Double_t min, Double_t binwidth = 0.01);
TH1 * YieldMean_HighExtrapolationHisto(TH1 *h, TF1 *f, Double_t max, Double_t binwidth = 0.1);
TH1 * YieldMean_ReturnRandom(TH1 *hin);
TH1 * YieldMean_ReturnCoherentRandom(TH1 *hin);
TH1 * YieldMean_ReturnExtremeHisto(TH1 *hin, const char* hout_title, Float_t sign = 1.);
TH1 * YieldMean_ReturnExtremeHardHisto(TH1 *hin, const char* hout_title);
TH1 * YieldMean_ReturnExtremeSoftHisto(TH1 *hin, const char* hout_title);
TH1 * YieldMean_ReturnExtremeLowHisto(TH1 *hin, const char* hout_title);
TH1 * YieldMean_ReturnExtremeHighHisto(TH1 *hin, const char* hout_title);
int Fitter(TH1* histo, TF1* func, Option_t *opt);
void RandomShifter(TH1* hin, TH1* hhi, TH1* hlo, TH1*& hIntegral, TH1*& hMean, float integral_limit, float mean_limit, int nRepetitions = 1000);
void SaveToFile(const char* file_name, const char* dir_name, TCanvas& canvas, bool debug);
//TH1 * MergeHistograms(TH1 *hdata, TH1 *hlo, TH1 *hhi, const char* hout_name);

TH1 *
YieldMeanNew(TH1 *hstat, TH1 *hsys, TH1 *hsys_mult_corr, TH1 *hsys_mult_uncorr, TF1* &fout, TF1 *f = NULL, Double_t min = 0., Double_t max = 10., Double_t loprecision = 0.01, Double_t hiprecision = 0.1, bool store_log = false, const char* logfilename="log.root", const char* path = "", Option_t *opt = "0q", bool debug = false)
{
  /* set many iterations when fitting the data so we don't
     stop minimization with MAX_CALLS */
  TVirtualFitter::SetMaxIterations(1000000);

  /* create output histo */
  double integral, mean, integral_ref, mean_ref;
  TH1 *hout = new TH1D(Form("hout_%s",hstat->GetName()), "", kFitRes, 0, kFitRes);
  TH1 *hlo = nullptr, *hhi = nullptr;

  /* create histo with stat+sys errors */
  TH1 *htot = (TH1 *)hstat->Clone(Form("%sfittedwith%s",hstat->GetName(),f->GetName()));
  for (Int_t ibin = 0; ibin < htot->GetNbinsX(); ibin++) {
    htot->SetBinError(ibin + 1, TMath::Sqrt(hsys->GetBinError(ibin + 1) * hsys->GetBinError(ibin + 1) + hstat->GetBinError(ibin + 1) * hstat->GetBinError(ibin + 1)));
  }

  /*
   *   measure the central value
   */
  int fitres = Fitter(htot,f,opt);
  hout->SetBinContent(kFitRes,fitres);

  fout = new TF1(*f);
  fout->SetLineColor(kRed);
  cout<<"Fit sys+stat for " <<f->GetName()<<endl;
  cout<<"NDF="<<f->GetNDF()<<" Chi^2="<<f->GetChisquare()<<" Chi^2/NDF="<<f->GetChisquare()/f->GetNDF()<<endl;

  hlo = YieldMean_LowExtrapolationHisto(htot, f, min, loprecision);
  hhi = YieldMean_HighExtrapolationHisto(htot, f, max, hiprecision);

  YieldMean_IntegralMean(htot, hlo, hhi, integral, mean,kTRUE);
  hout->SetBinContent(kYield, integral);
  hout->SetBinContent(kMean, mean);
  integral_ref = integral;
  mean_ref = mean;

  TCanvas cCanvasDefault(Form("cCanvasDefault_%s",hstat->GetName()));
  cCanvasDefault.DrawFrame(min,0,max,1.3 * hstat->GetMaximum());
  htot->Draw("pesame");
  //hlo->Draw("pesame");
  //hhi->Draw("pesame");
  fout->Draw("same");


  if(store_log){
    SaveToFile(logfilename,path,cCanvasDefault,debug);
  }  

  /*
   * STATISTICS
   */

  TCanvas cCanvasStat(Form("cCanvasStat_%s",hstat->GetName()));
  cCanvasStat.Divide(2, 1);

  /*
   * measure statistical error
   */

  /* fit with stat error */
  Fitter(hstat,f,opt);
  hlo = YieldMean_LowExtrapolationHisto(hstat, f, min, loprecision);
  hhi = YieldMean_HighExtrapolationHisto(hstat, f, max, hiprecision);

  TH1* hIntegral = nullptr;
  TH1* hMean = nullptr;

  RandomShifter(hstat, hhi, hlo, hIntegral, hMean, integral_ref, mean_ref);

  TF1 *gaus = (TF1 *)gROOT->GetFunction("gaus");
  
  cCanvasStat.cd(1);
  hIntegral->Fit(gaus, "q");
  integral = hout->GetBinContent(kYield) * gaus->GetParameter(2) / gaus->GetParameter(1);
  hout->SetBinContent(kYieldStat, integral);

  cCanvasStat.cd(2);
  hMean->Fit(gaus, "q");
  mean = hout->GetBinContent(kMean) * gaus->GetParameter(2) / gaus->GetParameter(1);
  hout->SetBinContent(kMeanStat, mean);
  
  if(store_log){
    SaveToFile(logfilename,path,cCanvasStat,debug);
  }

  TCanvas cCanvasStatExtra(Form("cCanvasStatExtra_%s",hstat->GetName()));
  cCanvasStatExtra.DrawFrame(min,0,max,1.3 * hstat->GetMaximum());
  hstat->Draw("pesame");
  hlo->Draw("pesame");
  hhi->Draw("pesame");

  if(store_log){
    SaveToFile(logfilename,path,cCanvasStatExtra,debug);
  }
  delete hMean;
  delete hIntegral;

  /*
   * MULTIPLICITY UNCORRELATED SYSTEMATICS
   */

  /* fit with uncorrelated syst error */
  TCanvas cCanvasSystUncorr(Form("cCanvasSystUncorr_%s",hstat->GetName()));
  cCanvasSystUncorr.Divide(2, 1);

  Fitter(hsys_mult_uncorr,f,opt);
  delete hlo;
  delete hhi;
  hlo = YieldMean_LowExtrapolationHisto(hsys_mult_uncorr, f, min, loprecision);
  hhi = YieldMean_HighExtrapolationHisto(hsys_mult_uncorr, f, max, hiprecision);

  TH1* hIntegralSyst = nullptr;
  TH1* hMeanSyst = nullptr;

  RandomShifter(hsys_mult_uncorr, hhi, hlo, hIntegralSyst, hMeanSyst, integral_ref, mean_ref);
  
  cCanvasSystUncorr.cd(1);
  hIntegralSyst->Fit(gaus, "q");
  integral = hout->GetBinContent(kYield) * gaus->GetParameter(2) / gaus->GetParameter(1);
  hout->SetBinContent(kYieldSysUncorr, integral);

  cCanvasSystUncorr.cd(2);
  hMeanSyst->Fit(gaus, "q");
  mean = hout->GetBinContent(kMean) * gaus->GetParameter(2) / gaus->GetParameter(1);
  hout->SetBinContent(kMeanSysUncorr, mean);
  
  if(store_log){
    SaveToFile(logfilename,path,cCanvasSystUncorr,debug);
  }

  TCanvas cCanvasSystUncorrExtra(Form("cCanvasStatExtra_%s",hstat->GetName()));
  cCanvasSystUncorrExtra.DrawFrame(min,0,max,1.3 * hsys_mult_corr->GetMaximum());
  hsys_mult_corr->Draw("pesame");
  hlo->Draw("pesame");
  hhi->Draw("pesame");

  if(store_log){
    SaveToFile(logfilename,path,cCanvasSystUncorrExtra,debug);
  }

  delete hMeanSyst;
  delete hIntegralSyst;

  /*
   * MULTIPLICITY CORRELATED SYSTEMATICS
   */

  TCanvas cCanvasSysCorr(Form("cCanvasYieldSysCorr_%s",hstat->GetName()));
  cCanvasSysCorr.Divide(2, 1);
  cCanvasSysCorr.cd(1)->DrawFrame(min, 0.7 * hsys_mult_corr->GetMinimum(), max, 1.3 * hsys_mult_corr->GetMaximum());
  hsys_mult_corr->SetMarkerStyle(20);
  hsys_mult_corr->SetMarkerColor(1);
  hsys_mult_corr->SetMarkerSize(1);
  hsys_mult_corr->Draw("same");
  cCanvasSysCorr.cd(2)->DrawFrame(min, 0.7 * hsys_mult_corr->GetMinimum() , max, 1.3 * hsys_mult_corr->GetMaximum());
  hsys_mult_corr->Draw("same");

  /*
   * systematic error high
   */

  TH1 *hhigh_mult_corr = YieldMean_ReturnExtremeHighHisto(hsys_mult_corr,Form("%s_extreme_high",hsys_mult_corr->GetName()));
  Fitter(hhigh_mult_corr,f,opt);

  delete hlo;
  delete hhi;
  hlo = YieldMean_LowExtrapolationHisto(hhigh_mult_corr, f, min, loprecision);
  hhi = YieldMean_HighExtrapolationHisto(hhigh_mult_corr, f, max, hiprecision);
  YieldMean_IntegralMean(hhigh_mult_corr, hlo, hhi, integral, mean);
  integral = TMath::Abs(integral - hout->GetBinContent(kYield));
  hout->SetBinContent(kYieldSysHiCorr, integral);

  TCanvas cCanvasHighCorr(Form("cCanvasHighCorr_%s",hstat->GetName()));
  cCanvasHighCorr.DrawFrame(min,0,max,1.3 * hhigh_mult_corr->GetMaximum());
  hhigh_mult_corr->Draw("pesame");
  hlo->Draw("pesame");
  hhi->Draw("pesame");

  if(store_log){
    SaveToFile(logfilename,path,cCanvasHighCorr,debug);
  }

  cCanvasSysCorr.cd(1);
  f->SetLineColor(2);
  f->DrawCopy("same");

  /*
   * systematic error hard
   */

  TH1 *hhard_mult_corr = YieldMean_ReturnExtremeHardHisto(hsys_mult_corr,Form("%s_extreme_hard",hsys_mult_corr->GetName()));
  Fitter(hhard_mult_corr,f,opt);

  delete hlo;
  delete hhi;
  hlo = YieldMean_LowExtrapolationHisto(hhard_mult_corr, f, min, loprecision);
  hhi = YieldMean_HighExtrapolationHisto(hhard_mult_corr, f, max, hiprecision);
  YieldMean_IntegralMean(hhard_mult_corr, hlo, hhi, integral, mean);
  mean = TMath::Abs(mean - hout->GetBinContent(kMean));
  hout->SetBinContent(kMeanSysHiCorr, mean);

  TCanvas cCanvasHardCorr(Form("cCanvasHardCorr_%s",hstat->GetName()));
  cCanvasHardCorr.DrawFrame(min,0,max,1.3 * hhigh_mult_corr->GetMaximum());
  hhigh_mult_corr->Draw("pesame");
  hlo->Draw("pesame");
  hhi->Draw("pesame");

  if(store_log){
    SaveToFile(logfilename,path,cCanvasHardCorr,debug);
  }

  cCanvasSysCorr.cd(2);
  f->SetLineColor(2);
  f->DrawCopy("same");

  /*
   * systematic error low
   */

  TH1 *hlow_mult_corr = YieldMean_ReturnExtremeLowHisto(hsys_mult_corr, Form("%s_extreme_low",hsys_mult_corr->GetName()));
  Fitter(hlow_mult_corr,f,opt);

  delete hlo;
  delete hhi;
  hlo = YieldMean_LowExtrapolationHisto(hlow_mult_corr, f, min, loprecision);
  hhi = YieldMean_HighExtrapolationHisto(hlow_mult_corr, f, max, hiprecision);
  YieldMean_IntegralMean(hlow_mult_corr, hlo, hhi, integral, mean);
  integral = TMath::Abs(integral - hout->GetBinContent(kYield));
  hout->SetBinContent(kYieldSysLoCorr, integral);

  TCanvas cCanvasLowCorr(Form("cCanvasLowCorr_%s",hstat->GetName()));
  cCanvasLowCorr.DrawFrame(min,0,max,1.3 * hlow_mult_corr->GetMaximum());
  hlow_mult_corr->Draw("pesame");
  hlo->Draw("pesame");
  hhi->Draw("pesame");

  if(store_log){
    SaveToFile(logfilename,path,cCanvasLowCorr,debug);
  }

  cCanvasSysCorr.cd(1);
  f->SetLineColor(4);
  f->DrawCopy("same");

  /*
   * systematic error soft
   */

  TH1 *hsoft_mult_corr = YieldMean_ReturnExtremeSoftHisto(hsys_mult_corr, Form("%s_extreme_soft",hsys_mult_corr->GetName()));
  Fitter(hsoft_mult_corr,f,opt);

  delete hlo;
  delete hhi;
  hlo = YieldMean_LowExtrapolationHisto(hsoft_mult_corr, f, min, loprecision);
  hhi = YieldMean_HighExtrapolationHisto(hsoft_mult_corr, f, max, hiprecision);
  YieldMean_IntegralMean(hsoft_mult_corr, hlo, hhi, integral, mean);
  mean = TMath::Abs(mean - hout->GetBinContent(kMean));
  hout->SetBinContent(kMeanSysLoCorr, mean);

  TCanvas cCanvasSoftCorr(Form("cCanvasSoftCorr_%s",hstat->GetName()));
  cCanvasSoftCorr.DrawFrame(min,0,max,1.3 * hlow_mult_corr->GetMaximum());
  hlow_mult_corr->Draw("pesame");
  hlo->Draw("pesame");
  hhi->Draw("pesame");

  if(store_log){
    SaveToFile(logfilename,path,cCanvasSoftCorr,debug);
  }

  cCanvasSysCorr.cd(2);
  f->SetLineColor(4);
  f->DrawCopy("same");

  if(store_log){
    SaveToFile(logfilename,path,cCanvasSysCorr,debug);
  }

  cout << "dN / dy = " << hout->GetBinContent(kYield) << " +- " << hout->GetBinContent(kYieldStat);
  cout << " +- " << hout->GetBinContent(kYieldSysUncorr);
  cout << " + " << hout->GetBinContent(kYieldSysHiCorr) << " - " << hout->GetBinContent(kYieldSysLoCorr);
  cout << "\n<pT> = " << hout->GetBinContent(kMean) << " +- " << hout->GetBinContent(kMeanStat);
  cout << " + " << hout->GetBinContent(kMeanSysUncorr);
  cout << " + " << hout->GetBinContent(kMeanSysHiCorr) << " - " << hout->GetBinContent(kMeanSysLoCorr);
  cout << endl;

  delete hlo;
  delete hhi;
  return hout;
}

TH1 *
YieldMean_LowExtrapolationHisto(TH1 *h, TF1 *f, Double_t min, Double_t binwidth)
{
  /* find lowest edge in histo */
  Int_t binlo = -1;
  Double_t lo = 0.;
  for (Int_t ibin = 1; ibin < h->GetNbinsX() + 1; ibin++) {
    if (h->GetBinContent(ibin) != 0.) {
      binlo = ibin;
      lo = h->GetBinLowEdge(ibin);
      break;
    }
  }
  if (binlo == -1) ::Fatal("YieldMean_LowExtrapolationHisto","Lower bin is undefined !?");

  Int_t nbins = (lo - min) / binwidth;
  TH1 *hlo = new TH1F(Form("%s_lo",h->GetName()), "", nbins, min, lo);

  /* integrate function in histogram bins */
  Double_t cont, err, width;
  for (Int_t ibin = 0; ibin < hlo->GetNbinsX(); ibin++) {
    width = hlo->GetBinWidth(ibin + 1);
    cont = f->Integral(hlo->GetBinLowEdge(ibin + 1), hlo->GetBinLowEdge(ibin + 2), 1.e-6);//(Double_t *)0, 1.e-6);
    err = f->IntegralError(hlo->GetBinLowEdge(ibin + 1), hlo->GetBinLowEdge(ibin + 2), (Double_t *)0, (Double_t *)0, 1.e-6);
    hlo->SetBinContent(ibin + 1, cont / width);
    hlo->SetBinError(ibin + 1, err / width);
  }

  return hlo;
}

TH1 *
YieldMean_HighExtrapolationHisto(TH1 *h, TF1 *f, Double_t max, Double_t binwidth)
{
  /* find highest edge in histo */
  Int_t binhi = -1;
  Double_t hi = 0.;
  for (Int_t ibin = h->GetNbinsX(); ibin > 0; ibin--) {
    if (h->GetBinContent(ibin) != 0.) {
      binhi = ibin + 1;
      hi = h->GetBinLowEdge(ibin + 1);
      break;
    }
  }
  if (binhi == -1) ::Fatal("YieldMean_HighExtrapolationHisto","Higher bin is undefined !?");
  if(max<hi) {
    Printf("Warning! You should probably set a higher max value (Max = %f, hi = %f)", max, hi);
  }
  Int_t nbins = (max - hi) / binwidth;
  TH1 *hhi = new TH1F(Form("%s_hi",h->GetName()), "", nbins, hi, max);

  /* integrate function in histogram bins */
  Double_t cont, err, width;
  for (Int_t ibin = 0; ibin < hhi->GetNbinsX(); ibin++) {
    width = hhi->GetBinWidth(ibin + 1);
    cont = f->Integral(hhi->GetBinLowEdge(ibin + 1), hhi->GetBinLowEdge(ibin + 2), 1.e-6); //(Double_t *)0, 1.e-6);
    err = f->IntegralError(hhi->GetBinLowEdge(ibin + 1), hhi->GetBinLowEdge(ibin + 2), (Double_t *)0, (Double_t *)0, 1.e-6);
    hhi->SetBinContent(ibin + 1, cont / width);
    hhi->SetBinError(ibin + 1, err / width);
  }

  return hhi;
}

TH1 *
YieldMean_ReturnRandom(TH1 *hin)
{
  TH1 *hout = (TH1 *)hin->Clone("hout");
  hout->Reset();
  Double_t cont, err;
  for (Int_t ibin = 0; ibin < hin->GetNbinsX(); ibin++) {
    if (hin->GetBinError(ibin + 1) <= 0.) continue;
    cont = hin->GetBinContent(ibin + 1);
    err = hin->GetBinError(ibin + 1);
    hout->SetBinContent(ibin + 1, gRandom->Gaus(cont, err));
    hout->SetBinError(ibin + 1, err);
  }
  return hout;
}

TH1 *
YieldMean_ReturnCoherentRandom(TH1 *hin)
{
  TH1 *hout = (TH1 *)hin->Clone("hout");
  hout->Reset();
  Double_t cont, err, cohe;
  cohe = gRandom->Gaus(0., 1.);
  for (Int_t ibin = 0; ibin < hin->GetNbinsX(); ibin++) {
    if (hin->GetBinError(ibin + 1) <= 0.) continue;
    cont = hin->GetBinContent(ibin + 1);
    err = hin->GetBinError(ibin + 1);
    hout->SetBinContent(ibin + 1, cont + cohe * err);
    hout->SetBinError(ibin + 1, err);
  }
  return hout;
}

TH1 *
YieldMean_ReturnExtremeHighHisto(TH1 *hin, const char* hout_title)
{
  TH1 *hout = (TH1 *)hin->Clone(hout_title);
  for (Int_t ibin = 0; ibin < hin->GetNbinsX(); ibin++) {
    if (hin->GetBinError(ibin + 1) <= 0.) continue;
    Double_t val = hin->GetBinContent(ibin + 1);
    Double_t err = hin->GetBinError(ibin + 1);
    hout->SetBinContent(ibin + 1, val + err);
  }
  return hout;
}

TH1 *
YieldMean_ReturnExtremeLowHisto(TH1 *hin, const char* hout_title)
{
  TH1 *hout = (TH1 *)hin->Clone(hout_title);
  for (Int_t ibin = 0; ibin < hin->GetNbinsX(); ibin++) {
    if (hin->GetBinError(ibin + 1) <= 0.) continue;
    Double_t val = hin->GetBinContent(ibin + 1);
    Double_t err = hin->GetBinError(ibin + 1);
    hout->SetBinContent(ibin + 1, val - err);
  }
  return hout;
}

TH1 *
YieldMean_ReturnExtremeSoftHisto(TH1 *hin, const char* hout_title)
{
  return YieldMean_ReturnExtremeHisto(hin, hout_title, -1.);
}

TH1 *
YieldMean_ReturnExtremeHardHisto(TH1 *hin, const char* hout_title)
{
  return YieldMean_ReturnExtremeHisto(hin, hout_title, 1.);
}

TH1 *
YieldMean_ReturnExtremeHisto(TH1 *hin, const char* hout_title, Float_t sign)
{
  Double_t ptlow = -1., pthigh = -1.;
  for (Int_t ibin = 0; ibin < hin->GetNbinsX(); ibin++) {
    if (hin->GetBinError(ibin + 1) <= 0.) continue;
    ptlow = hin->GetBinLowEdge(ibin + 1);
    break;
  }
  for (Int_t ibin = hin->GetNbinsX(); ibin >= 0; ibin--) {
    if (hin->GetBinError(ibin + 1) <= 0.) continue;
    pthigh = hin->GetBinLowEdge(ibin + 2);
    break;
  }
  if (ptlow < 0. || pthigh < 0.) ::Fatal("YieldMean_ReturnExtremeHisto","Problem in the determination of the pt region");

  Double_t mean = hin->GetMean();
  Double_t maxdiff = 0.;
  TH1 *hmax = NULL;
  for (Int_t inode = 0; inode < hin->GetNbinsX(); inode++) {

    Double_t ptnode = hin->GetBinCenter(inode + 1);
    TH1 *hout = (TH1 *)hin->Clone(hout_title);

    for (Int_t ibin = 0; ibin < hin->GetNbinsX(); ibin++) {
      if (hin->GetBinError(ibin + 1) <= 0.) continue;
      Double_t val = hin->GetBinContent(ibin + 1);
      Double_t err = hin->GetBinError(ibin + 1);
      Double_t cen = hin->GetBinCenter(ibin + 1);
      if (cen < ptnode)
        err *= -1. + (cen - ptlow) / (ptnode - ptlow);
      else
        err *= (cen - ptnode) / (pthigh - ptnode);

      hout->SetBinContent(ibin + 1, val + sign * err);
    }

    Double_t diff = TMath::Abs(mean - hout->GetMean());
    if (diff > maxdiff) {
      //      printf("found max at %f\n", ptnode);
      if (hmax) delete hmax;
      hmax = (TH1 *)hout->Clone(Form("%s_hmax",hout_title));
      maxdiff = diff;
    }
    delete hout;
  }
  return hmax;
}

void YieldMean_IntegralMean(TH1 *hdata, TH1 *hlo, TH1 *hhi, Double_t &integral, Double_t &mean,Bool_t printinfo)
{

  /*
   * compute integrals
   */

  Double_t cont, err, width, cent;
  Double_t I = 0., IX = 0., Ierr = 0., IXerr = 0., Ilerr = 0., IXlerr = 0.;
  Double_t M = 0., Merr = 0., Mlerr = 0., C;
  Double_t dataonly=0.0;

  /* integrate the data */
  for (Int_t ibin = 0; ibin < hdata->GetNbinsX(); ibin++) {
    cent = hdata->GetBinCenter(ibin + 1);
    width = hdata->GetBinWidth(ibin + 1);
    cont = width * hdata->GetBinContent(ibin + 1);
    err = width * hdata->GetBinError(ibin + 1);
    if (err <= 0.) continue;
    I += cont;
    IX += cont * cent;
  }

  dataonly=I;
  /* integrate low */
  for (Int_t ibin = 0; ibin < hlo->GetNbinsX(); ibin++) {
    cent = hlo->GetBinCenter(ibin + 1);
    width = hlo->GetBinWidth(ibin + 1);
    cont = width * hlo->GetBinContent(ibin + 1);
    err = width * hlo->GetBinError(ibin + 1);
    if (err <= 0.) continue;
    I += cont;
    IX += cont * cent;
  }
  /* integrate high */
  for (Int_t ibin = 0; ibin < hhi->GetNbinsX(); ibin++) {
    cent = hhi->GetBinCenter(ibin + 1);
    width = hhi->GetBinWidth(ibin + 1);
    cont = width * hhi->GetBinContent(ibin + 1);
    err = width * hhi->GetBinError(ibin + 1);
    if (err <= 0.) continue;
    I += cont;
    IX += cont * cent;
  }

  /* set values */
  integral = I;
  mean = IX / I;
  if(printinfo)
  	cout<<"data only = "<<dataonly<<" total = "<<I<<" ratio= "<<dataonly/I<<endl;
}

void RandomShifter(TH1* hin, TH1* hhi, TH1* hlo, TH1*& hIntegral, TH1*& hMean, float integral_limit, float mean_limit, int nRepetitions){

  TH1F hIntegral_tmp("hIntegral_tmp", "", 1000, 0.75 * integral_limit, 1.25 * integral_limit); //NOTE: only used for the range of the next histograms
  TH1F hMean_tmp("hMean_tmp", "", 1000, 0.75 * mean_limit, 1.25 * mean_limit);
  
  double integral, mean;
  
  for (int irnd = 0; irnd < 100; irnd++) {
    /* get random histogram */
    TH1 *hrnd = YieldMean_ReturnRandom(hin); //QUESTION: why not coherent?
    /* fit */
    TH1 *hrndlo = YieldMean_ReturnCoherentRandom(hlo);
    TH1 *hrndhi = YieldMean_ReturnCoherentRandom(hhi);
    /* integrate */
    YieldMean_IntegralMean(hrnd, hrndlo, hrndhi, integral, mean);
    hIntegral_tmp.Fill(integral);
    hMean_tmp.Fill(mean);
    delete hrnd;
    delete hrndlo;
    delete hrndhi;
  }
  if(hIntegral){
    delete hIntegral;
  }
  hIntegral  = new TH1F(Form("hIntegral_%s",hin->GetName()), "", 100, hIntegral_tmp.GetMean() - 10. * hIntegral_tmp.GetRMS(), hIntegral_tmp.GetMean() + 10. * hIntegral_tmp.GetRMS());
  if(hMean){
    delete hMean;
  }
  hMean = new TH1F(Form("hMean_%s",hin->GetName()), "", 100, hMean_tmp.GetMean() - 10. * hMean_tmp.GetRMS(), hMean_tmp.GetMean() + 10. * hMean_tmp.GetRMS());
  for (Int_t irnd = 0; irnd < nRepetitions; irnd++) {
    /* get random histogram */
    TH1 *hrnd = YieldMean_ReturnRandom(hin); //QUESTION: why not coherent?
    /* fit */
    TH1 *hrndlo = YieldMean_ReturnCoherentRandom(hlo);
    TH1 *hrndhi = YieldMean_ReturnCoherentRandom(hhi);
    /* integrate */
    YieldMean_IntegralMean(hrnd, hrndlo, hrndhi, integral, mean);
    hIntegral->Fill(integral);
    hMean->Fill(mean);
    delete hrnd;
    delete hrndlo;
    delete hrndhi;
  }
}

int Fitter(TH1* histo, TF1* func, Option_t *opt){
  int fitres;
  int trials = 0;
  do {
    fitres = histo->Fit(func, opt);
    //Printf("Trial: %d", trials++);
    trials++;
    if(trials > 10) {
      Printf("FIT DOES NOT CONVERGE IN LINE %d",__LINE__);
      break;
    }
  }
  while (fitres != 0);
  return fitres;
}

void SaveToFile(const char* file_name, const char* dir_name, TCanvas& canvas, bool debug){
  TFile file(file_name,"UPDATE");
  bool isDir = file.cd(dir_name);
  if(!isDir){
    file.mkdir(dir_name);
    file.cd(dir_name);
  }
  if(debug) printf("directoy name: %s\n", dir_name);
  canvas.Write();
  file.Close();
}