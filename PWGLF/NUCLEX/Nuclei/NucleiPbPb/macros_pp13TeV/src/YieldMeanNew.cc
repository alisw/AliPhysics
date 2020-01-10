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
#include "src/Utils.h"

#include <memory>

#endif

using namespace utils;
using namespace std;
/* definition of the fields in the histogram returned */

namespace yieldmeannew{

  enum EValue_t {
    kYield = 1,
    kYieldStat,
    kYieldSysTot,
    kYieldSysHiCorr,
    kYieldSysLoCorr,
    kMean,
    kMeanStat,
    kMeanSysTot,
    kMeanSysHardCorr,
    kMeanSysSoftCorr,
    kFitRes
  };

  void YieldMean_IntegralMean(TH1 *hdata, TH1 *hlo, TH1 *hhi, double &integral, double &mean,Bool_t printinfo=kFALSE);
  TH1 * YieldMean_LowExtrapolationHisto(TH1 *h, TF1 *f, double min, double binwidth = 0.01);
  TH1 * YieldMean_HighExtrapolationHisto(TH1 *h, TF1 *f, double max, double binwidth = 0.1);
  TH1 * YieldMean_ReturnRandom(TH1 *hin);
  TH1 * YieldMean_ReturnCoherentRandom(TH1 *hin);
  TH1 * YieldMean_ReturnExtremeHisto(TH1 *hin, const char* hout_title, float sign = 1.);
  TH1 * YieldMean_ReturnExtremeHardHisto(TH1 *hin, const char* hout_title);
  TH1 * YieldMean_ReturnExtremeSoftHisto(TH1 *hin, const char* hout_title);
  TH1 * YieldMean_ReturnExtremeLowHisto(TH1 *hin, const char* hout_title);
  TH1 * YieldMean_ReturnExtremeHighHisto(TH1 *hin, const char* hout_title);
  int Fitter(TH1* histo, TF1* func, Option_t *opt);
  void RandomShifter(TH1* hin, TH1* hhi, TH1* hlo, TH1*& hIntegral, TH1*& hMean, float integral_limit, float mean_limit, int nRepetitions = 1000);
  void SaveToFile(const char* file_name, const char* dir_name, TCanvas& canvas);
  //TH1 * MergeHistograms(TH1 *hdata, TH1 *hlo, TH1 *hhi, const char* hout_name);

  TH1 *
  YieldMeanNew(TH1 *hstat, TH1 *hsys, TH1 *hsys_pt_uncorr, TH1 *hsys_pt_corr, TH1 *hsys_mult_corr, TF1* &fout, TF1 *f = nullptr, double min = 0., double max = 10., double loprecision = 0.01, double hiprecision = 0.1, bool store_log = false, const char* logfilename="log.root", const char* path = "", Option_t *opt = "0q")
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
    for (int ibin = 0; ibin < htot->GetNbinsX(); ibin++) {
      htot->SetBinError(ibin + 1, TMath::Sqrt(hsys->GetBinError(ibin + 1) * hsys->GetBinError(ibin + 1) + hstat->GetBinError(ibin + 1) * hstat->GetBinError(ibin + 1)));
    }

    /*
    *   measure the central value
    */
    int fitres = Fitter(htot,f,opt);
    hout->SetBinContent(kFitRes,fitres);

    fout = new TF1(*f);
    fout->SetLineColor(kRed);
    std::cout<<"Fit sys+stat for " <<f->GetName()<<std::endl;
    std::cout<<"NDF="<<f->GetNDF()<<" Chi^2="<<f->GetChisquare()<<" Chi^2/NDF="<<f->GetChisquare()/f->GetNDF()<<std::endl;

    hlo = YieldMean_LowExtrapolationHisto(htot, f, min, loprecision);
    hhi = YieldMean_HighExtrapolationHisto(htot, f, max, hiprecision);

    YieldMean_IntegralMean(htot, hlo, hhi, integral, mean, true);
    hout->SetBinContent(kYield, integral);
    hout->SetBinContent(kMean, mean);
    integral_ref = integral;
    mean_ref = mean;

    TCanvas cCanvasDefault("cCanvasDefault");
    cCanvasDefault.DrawFrame(min,0,max,1.3 * hstat->GetMaximum(),";#it{p}_{T} (GeV/#it{c});");
    htot->Draw("pesame");
    //hlo->Draw("pesame");
    //hhi->Draw("pesame");
    fout->Draw("same");


    if(store_log){
      SaveToFile(logfilename,path,cCanvasDefault);
    }

    /*
    * STATISTICS
    */

    TCanvas cCanvasStat("cCanvasStat");
    cCanvasStat.Divide(2, 1);

    /*
    * measure statistical error
    */

    /* fit with stat error */
    Fitter(hstat,f,opt);
    delete hlo;
    delete hhi;
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
      SaveToFile(logfilename,path,cCanvasStat);
    }

    TCanvas cCanvasStatExtra("cCanvasStatExtra");
    cCanvasStatExtra.DrawFrame(min,0,max,1.3 * hstat->GetMaximum(),";#it{p}_{T} (GeV/#it{c});");
    hstat->Draw("pesame");
    hlo->Draw("pesame");
    hhi->Draw("pesame");

    if(store_log){
      SaveToFile(logfilename,path,cCanvasStatExtra);
    }
    delete hMean;
    delete hIntegral;

    /*
    * PT-UNCORRELATED SYSTEMATICS
    */

    /* fit with pt-uncorrelated syst error */
    TCanvas cCanvasSystPtUncorr("cCanvasSystPtUncorr");
    cCanvasSystPtUncorr.Divide(2, 1);

    Fitter(hsys_pt_uncorr,f,opt);
    delete hlo;
    delete hhi;
    hlo = YieldMean_LowExtrapolationHisto(hsys_pt_uncorr, f, min, loprecision);
    hhi = YieldMean_HighExtrapolationHisto(hsys_pt_uncorr, f, max, hiprecision);

    TH1* hIntegralSyst = nullptr;
    TH1* hMeanSyst = nullptr;

    RandomShifter(hsys_pt_uncorr, hhi, hlo, hIntegralSyst, hMeanSyst, integral_ref, mean_ref);
    
    cCanvasSystPtUncorr.cd(1);
    hIntegralSyst->Fit(gaus, "q");
    integral = hout->GetBinContent(kYield) * gaus->GetParameter(2) / gaus->GetParameter(1);
    hout->SetBinContent(kYieldSysTot, integral);

    cCanvasSystPtUncorr.cd(2);
    hMeanSyst->Fit(gaus, "q");
    mean = hout->GetBinContent(kMean) * gaus->GetParameter(2) / gaus->GetParameter(1);
    hout->SetBinContent(kMeanSysTot, mean);
    
    if(store_log){
      SaveToFile(logfilename,path,cCanvasSystPtUncorr);
    }

    TCanvas cCanvasSystPtUncorrExtra("cCanvasSystPtUncorrExtra");
    cCanvasSystPtUncorrExtra.DrawFrame(min,0,max,1.3 * hsys_pt_uncorr->GetMaximum(),";#it{p}_{T} (GeV/#it{c});");
    hsys_pt_uncorr->Draw("pesame");
    hlo->Draw("pesame");
    hhi->Draw("pesame");

    if(store_log){
      SaveToFile(logfilename,path,cCanvasSystPtUncorrExtra);
    }

    delete hMeanSyst;
    delete hIntegralSyst;

    /*
    * PT-CORRELATED SYSTEMATICS
    */

    TCanvas cCanvasSysPtCorr("cCanvasYieldSysPtCorr");
    cCanvasSysPtCorr.Divide(2, 1);
    cCanvasSysPtCorr.cd(1)->DrawFrame(min, 0.7 * hsys_pt_corr->GetMinimum(), max, 1.3 * hsys_pt_corr->GetMaximum(),";#it{p}_{T} (GeV/#it{c});");
    hsys_pt_corr->SetMarkerStyle(20);
    hsys_pt_corr->SetMarkerColor(1);
    hsys_pt_corr->SetMarkerSize(1);
    hsys_pt_corr->Draw("same");
    cCanvasSysPtCorr.cd(2)->DrawFrame(min, 0.7 * hsys_pt_corr->GetMinimum() , max, 1.3 * hsys_pt_corr->GetMaximum(),";#it{p}_{T} (GeV/#it{c});");
    hsys_pt_corr->Draw("same");

    /*
    * systematic error high
    */

    TH1 *hhigh_pt_corr = YieldMean_ReturnExtremeHighHisto(hsys_pt_corr,"pt_corr_extreme_high");
    Fitter(hhigh_pt_corr,f,opt);

    delete hlo;
    delete hhi;
    hlo = YieldMean_LowExtrapolationHisto(hhigh_pt_corr, f, min, loprecision);
    hhi = YieldMean_HighExtrapolationHisto(hhigh_pt_corr, f, max, hiprecision);
    YieldMean_IntegralMean(hhigh_pt_corr, hlo, hhi, integral, mean);
    integral = TMath::Abs(integral - hout->GetBinContent(kYield));
    integral = sqrt(Sq(integral)+Sq(hout->GetBinContent(kYieldSysTot)));
    hout->SetBinContent(kYieldSysTot, integral);

    TCanvas cCanvasHighPtCorr("cCanvasHighPtCorr");
    cCanvasHighPtCorr.DrawFrame(min,0,max,1.3 * hhigh_pt_corr->GetMaximum(),";#it{p}_{T} (GeV/#it{c});");
    hhigh_pt_corr->Draw("pesame");
    hlo->Draw("pesame");
    hhi->Draw("pesame");

    if(store_log){
      SaveToFile(logfilename,path,cCanvasHighPtCorr);
    }

    cCanvasSysPtCorr.cd(1);
    f->SetLineColor(2);
    f->DrawCopy("same");

    /*
    * systematic error hard
    */

    TH1 *hhard_pt_corr = YieldMean_ReturnExtremeHardHisto(hsys_pt_corr,"pt_corr_extreme_hard");
    Fitter(hhard_pt_corr,f,opt);

    delete hlo;
    delete hhi;
    hlo = YieldMean_LowExtrapolationHisto(hhard_pt_corr, f, min, loprecision);
    hhi = YieldMean_HighExtrapolationHisto(hhard_pt_corr, f, max, hiprecision);
    YieldMean_IntegralMean(hhard_pt_corr, hlo, hhi, integral, mean);
    mean = TMath::Abs(mean - hout->GetBinContent(kMean));
    mean = sqrt(Sq(mean) + Sq(hout->GetBinContent(kMeanSysTot)));
    hout->SetBinContent(kMeanSysTot, mean);

    TCanvas cCanvasHardPtCorr("cCanvasHardPtCorr");
    cCanvasHardPtCorr.DrawFrame(min,0,max,1.3 * hhard_pt_corr->GetMaximum(),";#it{p}_{T} (GeV/#it{c});");
    hhard_pt_corr->Draw("pesame");
    hlo->Draw("pesame");
    hhi->Draw("pesame");

    if(store_log){
      SaveToFile(logfilename,path,cCanvasHardPtCorr);
    }

    cCanvasSysPtCorr.cd(2);
    f->SetLineColor(2);
    f->DrawCopy("same");

    /*
    * systematic error low
    */

    TH1 *hlow_pt_corr = YieldMean_ReturnExtremeLowHisto(hsys_pt_corr, "pt_corr_extreme_low");
    Fitter(hlow_pt_corr,f,opt);

    delete hlo;
    delete hhi;
    hlo = YieldMean_LowExtrapolationHisto(hlow_pt_corr, f, min, loprecision);
    hhi = YieldMean_HighExtrapolationHisto(hlow_pt_corr, f, max, hiprecision);
    YieldMean_IntegralMean(hlow_pt_corr, hlo, hhi, integral, mean);
    integral = TMath::Abs(integral - hout->GetBinContent(kYield));
    integral = sqrt(Sq(integral) + Sq(hout->GetBinContent(kYieldSysTot)));
    hout->SetBinContent(kYieldSysTot, integral);

    TCanvas cCanvasLowPtCorr("cCanvasLowPtCorr");
    cCanvasLowPtCorr.DrawFrame(min,0,max,1.3 * hlow_pt_corr->GetMaximum(),";#it{p}_{T} (GeV/#it{c});");
    hlow_pt_corr->Draw("pesame");
    hlo->Draw("pesame");
    hhi->Draw("pesame");

    if(store_log){
      SaveToFile(logfilename,path,cCanvasLowPtCorr);
    }

    cCanvasSysPtCorr.cd(1);
    f->SetLineColor(4);
    f->DrawCopy("same");

    /*
    * systematic error soft
    */

    TH1 *hsoft_pt_corr = YieldMean_ReturnExtremeSoftHisto(hsys_pt_corr, "pt_corr_extreme_soft");
    Fitter(hsoft_pt_corr,f,opt);

    delete hlo;
    delete hhi;
    hlo = YieldMean_LowExtrapolationHisto(hsoft_pt_corr, f, min, loprecision);
    hhi = YieldMean_HighExtrapolationHisto(hsoft_pt_corr, f, max, hiprecision);
    YieldMean_IntegralMean(hsoft_pt_corr, hlo, hhi, integral, mean);
    mean = TMath::Abs(mean - hout->GetBinContent(kMean));
    mean = sqrt(Sq(mean) + Sq(hout->GetBinContent(kMeanSysTot)));
    hout->SetBinContent(kMeanSysTot, mean);

    TCanvas cCanvasSoftPtCorr("cCanvasSoftPtCorr");
    cCanvasSoftPtCorr.DrawFrame(min,0,max,1.3 * hsoft_pt_corr->GetMaximum(),";#it{p}_{T} (GeV/#it{c});");
    hsoft_pt_corr->Draw("pesame");
    hlo->Draw("pesame");
    hhi->Draw("pesame");

    if(store_log){
      SaveToFile(logfilename,path,cCanvasSoftPtCorr);
    }

    cCanvasSysPtCorr.cd(2);
    f->SetLineColor(4);
    f->DrawCopy("same");

    if(store_log){
      SaveToFile(logfilename,path,cCanvasSysPtCorr);
    }

    /*
    * MULTIPLICITY CORRELATED SYSTEMATICS
    */

    TCanvas cCanvasSysMultCorr("cCanvasYieldSysMultCorr");
    cCanvasSysMultCorr.Divide(2, 1);
    cCanvasSysMultCorr.cd(1)->DrawFrame(min, 0.7 * hsys_mult_corr->GetMinimum(), max, 1.3 * hsys_mult_corr->GetMaximum(),";#it{p}_{T} (GeV/#it{c});");
    hsys_mult_corr->SetMarkerStyle(20);
    hsys_mult_corr->SetMarkerColor(1);
    hsys_mult_corr->SetMarkerSize(1);
    hsys_mult_corr->Draw("same");
    cCanvasSysMultCorr.cd(2)->DrawFrame(min, 0.7 * hsys_mult_corr->GetMinimum() , max, 1.3 * hsys_mult_corr->GetMaximum(),";#it{p}_{T} (GeV/#it{c});");
    hsys_mult_corr->Draw("same");

    /*
    * systematic error high
    */

    TH1 *hhigh_mult_corr = YieldMean_ReturnExtremeHighHisto(hsys_mult_corr,"mult_corr_extreme_high");
    Fitter(hhigh_mult_corr,f,opt);

    delete hlo;
    delete hhi;
    hlo = YieldMean_LowExtrapolationHisto(hhigh_mult_corr, f, min, loprecision);
    hhi = YieldMean_HighExtrapolationHisto(hhigh_mult_corr, f, max, hiprecision);
    YieldMean_IntegralMean(hhigh_mult_corr, hlo, hhi, integral, mean);
    integral = TMath::Abs(integral - hout->GetBinContent(kYield));
    hout->SetBinContent(kYieldSysHiCorr, integral);

    TCanvas cCanvasHighMultCorr("cCanvasHighMultCorr");
    cCanvasHighMultCorr.DrawFrame(min,0,max,1.3 * hhigh_mult_corr->GetMaximum(),";#it{p}_{T} (GeV/#it{c});");
    hhigh_mult_corr->Draw("pesame");
    hlo->Draw("pesame");
    hhi->Draw("pesame");

    if(store_log){
      SaveToFile(logfilename,path,cCanvasHighMultCorr);
    }

    cCanvasSysMultCorr.cd(1);
    f->SetLineColor(2);
    f->DrawCopy("same");

    /*
    * systematic error hard
    */

    TH1 *hhard_mult_corr = YieldMean_ReturnExtremeHardHisto(hsys_mult_corr,"mult_corr_extreme_hard");
    Fitter(hhard_mult_corr,f,opt);

    delete hlo;
    delete hhi;
    hlo = YieldMean_LowExtrapolationHisto(hhard_mult_corr, f, min, loprecision);
    hhi = YieldMean_HighExtrapolationHisto(hhard_mult_corr, f, max, hiprecision);
    YieldMean_IntegralMean(hhard_mult_corr, hlo, hhi, integral, mean);
    mean = TMath::Abs(mean - hout->GetBinContent(kMean));
    hout->SetBinContent(kMeanSysHardCorr, mean);

    TCanvas cCanvasHardMultCorr("cCanvasHardMultCorr");
    cCanvasHardMultCorr.DrawFrame(min,0,max,1.3 * hhard_mult_corr->GetMaximum(),";#it{p}_{T} (GeV/#it{c});");
    hhard_mult_corr->Draw("pesame");
    hlo->Draw("pesame");
    hhi->Draw("pesame");

    if(store_log){
      SaveToFile(logfilename,path,cCanvasHardMultCorr);
    }

    cCanvasSysMultCorr.cd(2);
    f->SetLineColor(2);
    f->DrawCopy("same");

    /*
    * systematic error low
    */

    TH1 *hlow_mult_corr = YieldMean_ReturnExtremeLowHisto(hsys_mult_corr, "mult_corr_extreme_low");
    Fitter(hlow_mult_corr,f,opt);

    delete hlo;
    delete hhi;
    hlo = YieldMean_LowExtrapolationHisto(hlow_mult_corr, f, min, loprecision);
    hhi = YieldMean_HighExtrapolationHisto(hlow_mult_corr, f, max, hiprecision);
    YieldMean_IntegralMean(hlow_mult_corr, hlo, hhi, integral, mean);
    integral = TMath::Abs(integral - hout->GetBinContent(kYield));
    hout->SetBinContent(kYieldSysLoCorr, integral);

    TCanvas cCanvasLowMultCorr("cCanvasLowMultCorr");
    cCanvasLowMultCorr.DrawFrame(min,0,max,1.3 * hlow_mult_corr->GetMaximum(),";#it{p}_{T} (GeV/#it{c});");
    hlow_mult_corr->Draw("pesame");
    hlo->Draw("pesame");
    hhi->Draw("pesame");

    if(store_log){
      SaveToFile(logfilename,path,cCanvasLowMultCorr);
    }

    cCanvasSysMultCorr.cd(1);
    f->SetLineColor(4);
    f->DrawCopy("same");

    /*
    * systematic error soft
    */

    TH1 *hsoft_mult_corr = YieldMean_ReturnExtremeSoftHisto(hsys_mult_corr, "mult_corr_extreme_soft");
    Fitter(hsoft_mult_corr,f,opt);

    delete hlo;
    delete hhi;

    hlo = YieldMean_LowExtrapolationHisto(hsoft_mult_corr, f, min, loprecision);
    hhi = YieldMean_HighExtrapolationHisto(hsoft_mult_corr, f, max, hiprecision);
 
    YieldMean_IntegralMean(hsoft_mult_corr, hlo, hhi, integral, mean);
    mean = TMath::Abs(mean - hout->GetBinContent(kMean));
    hout->SetBinContent(kMeanSysSoftCorr, mean);

    TCanvas cCanvasSoftMultCorr("cCanvasSoftMultCorr");
    cCanvasSoftMultCorr.DrawFrame(min,0,max,1.3 * hsoft_mult_corr->GetMaximum(),";#it{p}_{T} (GeV/#it{c});");
    hsoft_mult_corr->Draw("pesame");
    hlo->Draw("pesame");
    hhi->Draw("pesame");    

    if(store_log){
      SaveToFile(logfilename,path,cCanvasSoftMultCorr);
    }

    cCanvasSysMultCorr.cd(2);
    f->SetLineColor(4);
    f->DrawCopy("same");

    if(store_log){
      SaveToFile(logfilename,path,cCanvasSysMultCorr);
    }

    std::cout << "dN / dy = " << hout->GetBinContent(kYield) << " +- " << hout->GetBinContent(kYieldStat);
    std::cout << " +- " << hout->GetBinContent(kYieldSysTot);
    std::cout << " ( + " << hout->GetBinContent(kYieldSysHiCorr) << " - " << hout->GetBinContent(kYieldSysLoCorr) << ") ";
    std::cout << "\n<pT> = " << hout->GetBinContent(kMean) << " +- " << hout->GetBinContent(kMeanStat);
    std::cout << " + " << hout->GetBinContent(kMeanSysTot);
    std::cout << " ( + " << hout->GetBinContent(kMeanSysHardCorr) << " - " << hout->GetBinContent(kMeanSysSoftCorr) << ") ";
    std::cout << std::endl;

    delete hlo;
    delete hhi;
    return hout;
  }

  TH1 *
  YieldMean_LowExtrapolationHisto(TH1 *h, TF1 *f, double min, double binwidth)
  {
    /* find lowest edge in histo */
    int binlo = -1;
    double lo = 0.;
    for (int ibin = 1; ibin < h->GetNbinsX() + 1; ibin++) {
      if (h->GetBinContent(ibin) != 0.) {
        binlo = ibin;
        lo = h->GetBinLowEdge(ibin);
        break;
      }
    }
    if (binlo == -1) ::Fatal("YieldMean_LowExtrapolationHisto","Lower bin is undefined !?");

    int nbins = (lo - min) / binwidth;
    TH1 *hlo = new TH1F(Form("%s_lo",h->GetName()), "", nbins, min, lo);

    /* integrate function in histogram bins */
    double cont, err, width;
    for (int ibin = 0; ibin < hlo->GetNbinsX(); ibin++) {
      width = hlo->GetBinWidth(ibin + 1);
      cont = f->Integral(hlo->GetBinLowEdge(ibin + 1), hlo->GetBinLowEdge(ibin + 2), 1.e-6);//(double *)0, 1.e-6);
      err = f->IntegralError(hlo->GetBinLowEdge(ibin + 1), hlo->GetBinLowEdge(ibin + 2), (double *)0, (double *)0, 1.e-6);
      hlo->SetBinContent(ibin + 1, cont / width);
      hlo->SetBinError(ibin + 1, err / width);
    }

    return hlo;
  }

  TH1 *
  YieldMean_HighExtrapolationHisto(TH1 *h, TF1 *f, double max, double binwidth)
  {
    /* find highest edge in histo */
    int binhi = -1;
    double hi = 0.;
    for (int ibin = h->GetNbinsX(); ibin > 0; ibin--) {
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
    int nbins = (max - hi) / binwidth;
    TH1 *hhi = new TH1F(Form("%s_hi",h->GetName()), "", nbins, hi, max);

    /* integrate function in histogram bins */
    double cont, err, width;
    for (int ibin = 0; ibin < hhi->GetNbinsX(); ibin++) {
      width = hhi->GetBinWidth(ibin + 1);
      cont = f->Integral(hhi->GetBinLowEdge(ibin + 1), hhi->GetBinLowEdge(ibin + 2), 1.e-6); //(double *)0, 1.e-6);
      err = f->IntegralError(hhi->GetBinLowEdge(ibin + 1), hhi->GetBinLowEdge(ibin + 2), (double *)0, (double *)0, 1.e-6);
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
    double cont, err;
    for (int ibin = 0; ibin < hin->GetNbinsX(); ibin++) {
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
    double cont, err, cohe;
    cohe = gRandom->Gaus(0., 1.);
    for (int ibin = 0; ibin < hin->GetNbinsX(); ibin++) {
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
    for (int ibin = 0; ibin < hin->GetNbinsX(); ibin++) {
      if (hin->GetBinError(ibin + 1) <= 0.) continue;
      double val = hin->GetBinContent(ibin + 1);
      double err = hin->GetBinError(ibin + 1);
      hout->SetBinContent(ibin + 1, val + err);
      hout->SetBinError(ibin + 1, err);
    }
    return hout;
  }

  TH1 *
  YieldMean_ReturnExtremeLowHisto(TH1 *hin, const char* hout_title)
  {
    TH1 *hout = (TH1 *)hin->Clone(hout_title);
    for (int ibin = 0; ibin < hin->GetNbinsX(); ibin++) {
      if (hin->GetBinError(ibin + 1) <= 0.) continue;
      double val = hin->GetBinContent(ibin + 1);
      double err = hin->GetBinError(ibin + 1);
      hout->SetBinContent(ibin + 1, val - err);
      hout->SetBinError(ibin + 1, err);
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
  YieldMean_ReturnExtremeHisto(TH1 *hin, const char* hout_title, float sign)
  {
    double ptlow = -1., pthigh = -1.;
    for (int ibin = 0; ibin < hin->GetNbinsX(); ibin++) {
      if (hin->GetBinError(ibin + 1) <= 0.) continue;
      ptlow = hin->GetBinLowEdge(ibin + 1);
      break;
    }
    for (int ibin = hin->GetNbinsX(); ibin >= 0; ibin--) {
      if (hin->GetBinError(ibin + 1) <= 0.) continue;
      pthigh = hin->GetBinLowEdge(ibin + 2);
      break;
    }
    if (ptlow < 0. || pthigh < 0.) ::Fatal("YieldMean_ReturnExtremeHisto","Problem in the determination of the pt region");

    double mean = hin->GetMean();
    double maxdiff = 0.;
    TH1 *hmax = nullptr;
    for (int inode = 0; inode < hin->GetNbinsX(); inode++) {

      double ptnode = hin->GetBinCenter(inode + 1);
      TH1 *hout = (TH1 *)hin->Clone("extreme_tmp");

      for (int ibin = 0; ibin < hin->GetNbinsX(); ibin++) {
        if (hin->GetBinError(ibin + 1) <= 0.) continue;
        double val = hin->GetBinContent(ibin + 1);
        double err = hin->GetBinError(ibin + 1);
        double cen = hin->GetBinCenter(ibin + 1);
        double new_err = 0.;
        if (cen < ptnode)
          new_err = err * (-1. + (cen - ptlow) / (ptnode - ptlow));
        else
          new_err *= err * ((cen - ptnode) / (pthigh - ptnode));

        hout->SetBinContent(ibin + 1, val + sign * new_err);
        hout->SetBinError(ibin + 1, err);
      }

      double diff = TMath::Abs(mean - hout->GetMean());
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

  void YieldMean_IntegralMean(TH1 *hdata, TH1 *hlo, TH1 *hhi, double &integral, double &mean,Bool_t printinfo)
  {

    /*
    * compute integrals
    */

    double cont, err, width, cent;
    double I = 0., IX = 0., Ierr = 0., IXerr = 0., Ilerr = 0., IXlerr = 0.;
    double M = 0., Merr = 0., Mlerr = 0., C;
    double dataonly=0.0;

    /* integrate the data */
    for (int ibin = 0; ibin < hdata->GetNbinsX(); ibin++) {
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
    for (int ibin = 0; ibin < hlo->GetNbinsX(); ibin++) {
      cent = hlo->GetBinCenter(ibin + 1);
      width = hlo->GetBinWidth(ibin + 1);
      cont = width * hlo->GetBinContent(ibin + 1);
      err = width * hlo->GetBinError(ibin + 1);
      if (err <= 0.) continue;
      I += cont;
      IX += cont * cent;
    }
    /* integrate high */
    for (int ibin = 0; ibin < hhi->GetNbinsX(); ibin++) {
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
      std::cout<<"data only = "<<dataonly<<" total = "<<I<<" ratio= "<<dataonly/I<<std::endl;
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
    hIntegral  = new TH1F(Form("hIntegral_%s",hin->GetName()), ";d#it{N}/d#it{y};", 100, hIntegral_tmp.GetMean() - 10. * hIntegral_tmp.GetRMS(), hIntegral_tmp.GetMean() + 10. * hIntegral_tmp.GetRMS());
    if(hMean){
      delete hMean;
    }
    hMean = new TH1F(Form("hMean_%s",hin->GetName()), ";<#it{p}_{T}> (GeV/#it{c});", 100, hMean_tmp.GetMean() - 10. * hMean_tmp.GetRMS(), hMean_tmp.GetMean() + 10. * hMean_tmp.GetRMS());
    for (int irnd = 0; irnd < nRepetitions; irnd++) {
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

  void SaveToFile(const char* file_name, const char* dir_name, TCanvas& canvas){
    TFile file(file_name,"UPDATE");
    bool isDir = file.cd(dir_name);
    if(!isDir){
      file.mkdir(dir_name);
      file.cd(dir_name);
    }
    canvas.Write();
    file.Close();
  }
}