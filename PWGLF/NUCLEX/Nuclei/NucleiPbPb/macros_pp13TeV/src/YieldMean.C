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

#endif
using namespace std;
/* definition of the fields in the histogram returned */
namespace yieldmean{

  enum EValue_t {
    kYield = 1,
    kYieldStat,
    kYieldSysHi,
    kYieldSysLo,
    kMean,
    kMeanStat,
    kMeanSysHi,
    kMeanSysLo,
    kFitRes
  };

  void YieldMean_IntegralMean(TH1 *hdata, TH1 *hlo, TH1 *hhi, Double_t &integral, Double_t &mean,Bool_t printinfo=kFALSE);
  TH1 * YieldMean_LowExtrapolationHisto(TH1 *h, TF1 *f, Double_t min, Double_t binwidth = 0.01);
  TH1 * YieldMean_HighExtrapolationHisto(TH1 *h, TF1 *f, Double_t max, Double_t binwidth = 0.1);
  TH1 * YieldMean_ReturnRandom(TH1 *hin);
  TH1 * YieldMean_ReturnCoherentRandom(TH1 *hin);
  TH1 * YieldMean_ReturnExtremeHisto(TH1 *hin, Float_t sign = 1.);
  TH1 * YieldMean_ReturnExtremeHardHisto(TH1 *hin);
  TH1 * YieldMean_ReturnExtremeSoftHisto(TH1 *hin);
  TH1 * YieldMean_ReturnExtremeLowHisto(TH1 *hin);
  TH1 * YieldMean_ReturnExtremeHighHisto(TH1 *hin);
  //TH1 * MergeHistograms(TH1 *hdata, TH1 *hlo, TH1 *hhi, const char* hout_name);

  TH1 *
  YieldMean(TH1 *hstat, TH1 *hsys, TF1* &fout, TF1 *f = NULL, Double_t min = 0., Double_t max = 10., Double_t loprecision = 0.01, Double_t hiprecision = 0.1, bool store_log = false, const char* logfilename="log.root", const char* path = "", Option_t *opt = "0q")
  {
    /* set many iterations when fitting the data so we don't
      stop minimization with MAX_CALLS */
    TVirtualFitter::SetMaxIterations(1000000);

    /* create output histo */
    Double_t integral, mean;
    TH1 *hout = new TH1D(Form("hout_%s",hstat->GetName()), "", 9, 0, 9);
    TH1 *hlo = 0x0, *hhi = 0x0;

    /* create histo with stat+sys errors */
    TH1 *htot = (TH1 *)hstat->Clone(Form("%sfittedwith%s",hstat->GetName(),f->GetName()));
    for (Int_t ibin = 0; ibin < htot->GetNbinsX(); ibin++) {
      htot->SetBinError(ibin + 1, TMath::Sqrt(hsys->GetBinError(ibin + 1) * hsys->GetBinError(ibin + 1) + hstat->GetBinError(ibin + 1) * hstat->GetBinError(ibin + 1)));
    }

    /*
    *   measure the central value
    */
    int fitres;
    int trials = 0;
    do {
      fitres = htot->Fit(f, opt);
      Printf("Trial: %d", trials++);
      if(trials > 10) {
        Printf("FIT DOES NOT CONVERGE IN LINE %d",__LINE__);
        break;
      }
    }
    while (fitres != 0);
    hout->SetBinContent(kFitRes,fitres);

    fout = new TF1(*f);
    cout<<" Fit sys+stat for " <<f->GetName()<<endl;
    cout<<"NDF="<<f->GetNDF()<<" Chi^2="<<f->GetChisquare()<<" Chi^2/NDF="<<f->GetChisquare()/f->GetNDF()<<endl;

    hlo = YieldMean_LowExtrapolationHisto(htot, f, min, loprecision);
    hhi = YieldMean_HighExtrapolationHisto(htot, f, max, hiprecision);

    YieldMean_IntegralMean(htot, hlo, hhi, integral, mean,kTRUE);
    hout->SetBinContent(kYield, integral);
    hout->SetBinContent(kMean, mean);

    TCanvas *cCanvasDefault = new TCanvas(Form("cCanvasDefault_%s",hstat->GetName()));
    cCanvasDefault->DrawFrame(min,0,max,1.3 * hstat->GetMaximum());
    htot->Draw("pesame");
    hlo->Draw("pesame");
    hhi->Draw("pesame");


    if(store_log){
      TFile* filewithfits=TFile::Open(logfilename,"UPDATE");
      filewithfits->mkdir(path);
      filewithfits->cd(path);
      cCanvasDefault->Write();
      filewithfits->Close();
      delete filewithfits;
    }

    /*
    * STATISTICS
    */

    TCanvas *cCanvasStat = new TCanvas(Form("cCanvasStat_%s",hstat->GetName()));
    cCanvasStat->Divide(2, 1);

    /*
    * measure statistical error
    */

    /* fit with stat error */
    trials = 0;
    do {
      fitres = hstat->Fit(f, opt);
      Printf("Trial: %d", trials++);
      if(trials > 10) {
        Printf("FIT DOES NOT CONVERGE IN LINE %d",__LINE__);
        break;
      }
    }
    while (fitres != 0);
    hlo = YieldMean_LowExtrapolationHisto(hstat, f, min, loprecision);
    hhi = YieldMean_HighExtrapolationHisto(hstat, f, max, hiprecision);

    /* random generation with integration (coarse) */
    TH1 *hIntegral_tmp = new TH1F("hIntegral_tmp", "", 1000, 0.75 * integral, 1.25 * integral); //NOTE: only used for the range of the next histograms
    TH1 *hMean_tmp = new TH1F("hMean_tmp", "", 1000, 0.75 * mean, 1.25 * mean);
    for (Int_t irnd = 0; irnd < 100; irnd++) {
      /* get random histogram */
      TH1 *hrnd = YieldMean_ReturnRandom(hstat); //QUESTION: why not coherent?
      /* fit */
      TH1 *hrndlo = YieldMean_ReturnCoherentRandom(hlo);
      TH1 *hrndhi = YieldMean_ReturnCoherentRandom(hhi);
      /* integrate */
      YieldMean_IntegralMean(hrnd, hrndlo, hrndhi, integral, mean);
      hIntegral_tmp->Fill(integral);
      hMean_tmp->Fill(mean);
      delete hrnd;
      delete hrndlo;
      delete hrndhi;
    }
    /* random generation with integration (fine) */
    TH1 *hIntegral = new TH1F(Form("hIntegral_%s",hstat->GetName()), "", 100,
                              hIntegral_tmp->GetMean() - 10. * hIntegral_tmp->GetRMS(),
                              hIntegral_tmp->GetMean() + 10. * hIntegral_tmp->GetRMS());
    TH1 *hMean = new TH1F(Form("hMean_%s",hstat->GetName()), "", 100,
                          hMean_tmp->GetMean() - 10. * hMean_tmp->GetRMS(),
                          hMean_tmp->GetMean() + 10. * hMean_tmp->GetRMS());
    for (Int_t irnd = 0; irnd < 1000; irnd++) {
      /* get random histogram */
      TH1 *hrnd = YieldMean_ReturnRandom(hstat); //QUESTION: why not coherent?
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
    delete hMean_tmp;
    delete hIntegral_tmp;
    TF1 *gaus = (TF1 *)gROOT->GetFunction("gaus");

    cCanvasStat->cd(1);
    hIntegral->Fit(gaus, "q");
    integral = hout->GetBinContent(kYield) * gaus->GetParameter(2) / gaus->GetParameter(1);
    hout->SetBinContent(kYieldStat, integral);

    cCanvasStat->cd(2);
    hMean->Fit(gaus, "q");
    mean = hout->GetBinContent(kMean) * gaus->GetParameter(2) / gaus->GetParameter(1);
    hout->SetBinContent(kMeanStat, mean);

    TCanvas *cCanvasStatExtra = new TCanvas(Form("cCanvasStatExtra_%s",hstat->GetName()));
    cCanvasStatExtra->DrawFrame(min,0,max,1.3 * hstat->GetMaximum());
    hstat->Draw("pesame");
    hlo->Draw("pesame");
    hhi->Draw("pesame");

    if(store_log){
      TFile* filewithfits=TFile::Open(logfilename,"UPDATE");
      filewithfits->cd(path);
      cCanvasStatExtra->Write();
      filewithfits->Close();
      delete filewithfits;
    }

    /*
    * SYSTEMATICS
    */

    TCanvas *cCanvasSys = new TCanvas(Form("cCanvasYieldSys_%s",hstat->GetName()));
    cCanvasSys->Divide(2, 1);
    cCanvasSys->cd(1)->DrawFrame(min, 0.7 * hsys->GetMinimum(), max, 1.3 * hsys->GetMaximum());
    hsys->SetMarkerStyle(20);
    hsys->SetMarkerColor(1);
    hsys->SetMarkerSize(1);
    hsys->Draw("same");
    cCanvasSys->cd(2)->DrawFrame(min, 0.7 * hsys->GetMinimum() , max, 1.3 * hsys->GetMaximum());
    hsys->Draw("same");

    /*
    * systematic error high
    */

    TH1 *hhigh = YieldMean_ReturnExtremeHighHisto(hsys);
    trials = 0;
    do {
      fitres = hhigh->Fit(f, opt);
      Printf("Trial: %d", trials++);
      if(trials > 10) {
        Printf("FIT DOES NOT CONVERGE IN LINE %d",__LINE__);
        break;
      }
    }
    while (fitres != 0);

    delete hlo;
    delete hhi;
    hlo = YieldMean_LowExtrapolationHisto(hhigh, f, min, loprecision);
    hhi = YieldMean_HighExtrapolationHisto(hhigh, f, max, hiprecision);
    YieldMean_IntegralMean(hhigh, hlo, hhi, integral, mean);
    integral = TMath::Abs(integral - hout->GetBinContent(kYield));
    hout->SetBinContent(kYieldSysHi, integral);

    TCanvas *cCanvasHigh = new TCanvas(Form("cCanvasHigh_%s",hstat->GetName()));
    cCanvasHigh->DrawFrame(min,0,max,1.3 * hhigh->GetMaximum());
    hhigh->Draw("pesame");
    hlo->Draw("pesame");
    hhi->Draw("pesame");

    if(store_log){
      TFile* filewithfits=TFile::Open(logfilename,"UPDATE");
      filewithfits->cd(path);
      cCanvasHigh->Write();
      filewithfits->Close();
      delete filewithfits;
    }

    cCanvasSys->cd(1);
    f->SetLineColor(2);
    f->DrawCopy("same");

    /*
    * systematic error hard
    */

    TH1 *hhard = YieldMean_ReturnExtremeHardHisto(hsys);
    trials = 0;
    do {
      fitres = hhard->Fit(f, opt);
      Printf("Trial: %d", trials++);
      if(trials > 10) {
        Printf("FIT DOES NOT CONVERGE IN LINE %d",__LINE__);
        break;
      }
    }
    while (fitres != 0);

    delete hlo;
    delete hhi;
    hlo = YieldMean_LowExtrapolationHisto(hhard, f, min, loprecision);
    hhi = YieldMean_HighExtrapolationHisto(hhard, f, max, hiprecision);
    YieldMean_IntegralMean(hhard, hlo, hhi, integral, mean);
    mean = TMath::Abs(mean - hout->GetBinContent(kMean));
    hout->SetBinContent(kMeanSysHi, mean);

    TCanvas *cCanvasHard = new TCanvas(Form("cCanvasHard_%s",hstat->GetName()));
    cCanvasHard->DrawFrame(min,0,max,1.3 * hhigh->GetMaximum());
    hhigh->Draw("pesame");
    hlo->Draw("pesame");
    hhi->Draw("pesame");

    if(store_log){
      TFile* filewithfits=TFile::Open(logfilename,"UPDATE");
      filewithfits->cd(path);
      cCanvasHard->Write();
      filewithfits->Close();
      delete filewithfits;
    }

    cCanvasSys->cd(2);
    f->SetLineColor(2);
    f->DrawCopy("same");

    /*
    * systematic error low
    */

    TH1 *hlow = YieldMean_ReturnExtremeLowHisto(hsys);
    trials = 0;
    do {
      fitres = hlow->Fit(f, opt);
      Printf("Trial: %d", trials++);
      if(trials > 10) {
        Printf("FIT DOES NOT CONVERGE IN LINE %d",__LINE__);
        break;
      }
    }
    while (fitres != 0);

    delete hlo;
    delete hhi;
    hlo = YieldMean_LowExtrapolationHisto(hlow, f, min, loprecision);
    hhi = YieldMean_HighExtrapolationHisto(hlow, f, max, hiprecision);
    YieldMean_IntegralMean(hlow, hlo, hhi, integral, mean);
    integral = TMath::Abs(integral - hout->GetBinContent(kYield));
    hout->SetBinContent(kYieldSysLo, integral);

    TCanvas *cCanvasLow = new TCanvas(Form("cCanvasLow_%s",hstat->GetName()));
    cCanvasLow->DrawFrame(min,0,max,1.3 * hlow->GetMaximum());
    hlow->Draw("pesame");
    hlo->Draw("pesame");
    hhi->Draw("pesame");

    if(store_log){
      TFile* filewithfits=TFile::Open(logfilename,"UPDATE");
      filewithfits->cd(path);
      cCanvasLow->Write();
      filewithfits->Close();
      delete filewithfits;
    }

    cCanvasSys->cd(1);
    f->SetLineColor(4);
    f->DrawCopy("same");

    /*
    * systematic error soft
    */

    TH1 *hsoft = YieldMean_ReturnExtremeSoftHisto(hsys);
    trials = 0;
    do {
      fitres = hsoft->Fit(f, opt);
      Printf("Trial: %d", trials++);
      if(trials > 10) {
        Printf("FIT DOES NOT CONVERGE IN LINE %d",__LINE__);
        break;
      }
    }
    while (fitres != 0);

    delete hlo;
    delete hhi;
    hlo = YieldMean_LowExtrapolationHisto(hsoft, f, min, loprecision);
    hhi = YieldMean_HighExtrapolationHisto(hsoft, f, max, hiprecision);
    YieldMean_IntegralMean(hsoft, hlo, hhi, integral, mean);
    mean = TMath::Abs(mean - hout->GetBinContent(kMean));
    hout->SetBinContent(kMeanSysLo, mean);

    TCanvas *cCanvasSoft = new TCanvas(Form("cCanvasSoft_%s",hstat->GetName()));
    cCanvasSoft->DrawFrame(min,0,max,1.3 * hlow->GetMaximum());
    hlow->Draw("pesame");
    hlo->Draw("pesame");
    hhi->Draw("pesame");

    if(store_log){
      TFile* filewithfits=TFile::Open(logfilename,"UPDATE");
      filewithfits->cd(path);
      cCanvasSoft->Write();
      filewithfits->Close();
      delete filewithfits;
    }

    cCanvasSys->cd(2);
    f->SetLineColor(4);
    f->DrawCopy("same");

    cout << "dN / dy = " << hout->GetBinContent(kYield) << " +- " << hout->GetBinContent(kYieldStat);
    cout << " + " << hout->GetBinContent(kYieldSysHi) << " - " << hout->GetBinContent(kYieldSysLo);
    cout << "\n<pT> = " << hout->GetBinContent(kMean) << " +- " << hout->GetBinContent(kMeanStat);
    cout << " + " << hout->GetBinContent(kMeanSysHi) << " - " << hout->GetBinContent(kMeanSysLo);
    cout << endl;

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
    TH1 *hlo = new TH1F(Form("hlo_%s",f->GetName()), "", nbins, min, lo);

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
    TH1 *hhi = new TH1F(Form("hhi_%s",f->GetName()), "", nbins, hi, max);

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
  YieldMean_ReturnExtremeHighHisto(TH1 *hin)
  {
    TH1 *hout = (TH1 *)hin->Clone(Form("%s_extremehigh", hin->GetName()));
    for (Int_t ibin = 0; ibin < hin->GetNbinsX(); ibin++) {
      if (hin->GetBinError(ibin + 1) <= 0.) continue;
      Double_t val = hin->GetBinContent(ibin + 1);
      Double_t err = hin->GetBinError(ibin + 1);
      hout->SetBinContent(ibin + 1, val + err);
    }
    return hout;
  }

  TH1 *
  YieldMean_ReturnExtremeLowHisto(TH1 *hin)
  {
    TH1 *hout = (TH1 *)hin->Clone(Form("%s_extremelow", hin->GetName()));
    for (Int_t ibin = 0; ibin < hin->GetNbinsX(); ibin++) {
      if (hin->GetBinError(ibin + 1) <= 0.) continue;
      Double_t val = hin->GetBinContent(ibin + 1);
      Double_t err = hin->GetBinError(ibin + 1);
      hout->SetBinContent(ibin + 1, val - err);
    }
    return hout;
  }

  TH1 *
  YieldMean_ReturnExtremeSoftHisto(TH1 *hin)
  {
    return YieldMean_ReturnExtremeHisto(hin, -1.);
  }

  TH1 *
  YieldMean_ReturnExtremeHardHisto(TH1 *hin)
  {
    return YieldMean_ReturnExtremeHisto(hin, 1.);
  }

  TH1 *
  YieldMean_ReturnExtremeHisto(TH1 *hin, Float_t sign)
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
      TH1 *hout = (TH1 *)hin->Clone(Form("%s_extremehard", hin->GetName()));

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
        hmax = (TH1 *)hout->Clone("hmax");
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
}

// TH1 * MergeHistograms(TH1 *hdata, TH1 *hlo, TH1 *hhi, const char* hout_name){
//
//   std::vector<double> bin_limits;
//
//   for(int i=1; i <= hlo->GetNbinsX(); i++){
//     bin_limits.push_back(hlo->GetBinLowEdge(i));
//   }
//   double hlo_limit = hlo->GetBinLowEdge(hlo->GetNbinsX()+1);
//   for(int i=1; i <= hdata->GetNbinsX(); i++){
//     bin_limits.push_back(hdata->GetBinLowEdge(i));
//   }
//   double hdata_limit = hlo->GetBinLowEdge(hlo->GetNbinsX()+1);
//   for(int i=1; i <= hhi->GetNbinsX(); i++){
//     bin_limits.push_back(hhi->GetBinLowEdge(i));
//   }
//   bin_limits.push_back(hhi->GetBinLowEdge(hhi->GetNbinsX()+1));
//   int n_bins = (int)bin_limits.size();
//
//   TH1 *hout_tmp = (TH1*)hdata->Clone(Form("%s_tmp",hout_name));
//   hout_tmp->Reset();
//   TH1 * hout = (TH1*)hout->Rebin(n_bins,hout_name,bin_limits.data());
//   for(int i=1; i<=hout->GetNbinsX(); i++){
//     double bin_centre = hout->GetBinCenter(i);
//     if( bin_centre < hlo_limit){
//       hout->SetBinContent(i,hlo->GetBinContent(hlo->FindBin(bin_centre)));
//       hout->SetBinError(i,hlo->GetBinError(hlo->FindBin(bin_centre)));
//     } else if(bin_centre < hdata_limit){
//       hout->SetBinContent(i,hdata->GetBinContent(hdata->FindBin(bin_centre)));
//       hout->SetBinError(i,hdata->GetBinError(hdata->FindBin(bin_centre)));
//     } else{
//       hout->SetBinContent(i,hhi->GetBinContent(hhi->FindBin(bin_centre)));
//       hout->SetBinError(i,hhi->GetBinError(hhi->FindBin(bin_centre)));
//     }
//   }
//
//   return hout;
// }
