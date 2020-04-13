#include "src/Plotting.h"
using namespace plotting;

#include "src/Utils.h"
using namespace utils;

#include <TFile.h>
#include <TH2F.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <Riostream.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TMath.h>
#include "src/Common.h"
#include "src/YieldMean.C"
#include "src/YieldMeanNew.cc"
#include "AliPWGFunc.h"
#include "AdditionalFunctions.h"

#include <cstdio>
#include <fstream>

constexpr double kParticleMass = 1.87561;

//Levy-Tsallis
enum e_param_tsallis {e_mass, e_n, e_C, e_norm, e_chi2};
const char * param_names_tsallis[5] = {"mass","n","C","norm","chi2"};

//Levi-Tsallis parameters
const float normal=2e-4, normMin=1e-6, normMax=1.;
const float n=10, nMin=2, nMax=100;
const float C=0.2, CMin=0.01, CMax=0.4;

void Denormalize(TH1 * h) {
  for (int i = 1; i <= h->GetNbinsX(); ++i) {
    h->SetBinContent(i,h->GetBinContent(i) * TMath::TwoPi() * h->GetBinCenter(i));
    h->SetBinError(i,h->GetBinError(i) * TMath::TwoPi() * h->GetBinCenter(i));
  }
}

void CheckCorrelated(bool antimatter_analysys = false) {
  const char*  particlename = (!antimatter_analysys) ? "deuterons" : "antideuterons";
  gStyle->SetOptFit(1);
  gStyle->SetOptStat(0);
  gStyle->SetTitleXOffset(1.3);
  gStyle->SetTitleYOffset(1.);
  TFile *input_file = TFile::Open(kFinalOutput.data());
  TH1D *stat[kCentLength],*syst[kCentLength],*syst_pt_uncorr[kCentLength],*syst_pt_corr[kCentLength],*syst_mult_corr[kCentLength],*syst_mult_uncorr[kCentLength];
  TH1 *results_old[kCentLength], *results_new[kCentLength];

  AliPWGFunc pwgfunc;
  pwgfunc.SetVarType(AliPWGFunc::kdNdpt);
  TF1* fit_function = LevyTsallis("LevyTsallis", kParticleMass);
  //levi-tsallis parameters'
  double tsallis_param_val[kCentLength][5];
  double tsallis_param_err[kCentLength][5];
  const char* tsallis_param_name[5] = {"$m$ ($Gev/$c^{2}$)","$n$","$C$ (GeV)","\\mathrm{d}N/\\mathrm{d}y","\\chi^2 ndf"};
  //
  TFile outputfile(Form("%s/%s_checkcorr.root",kBaseOutputDir.data(),particlename),"recreate");
  TDirectory* datadir = outputfile.mkdir("data");
  for (int iC = 0; iC < kCentLength; ++iC) {
    stat[iC] = (TH1D*)input_file->Get(Form("%s/%i/stat",particlename,iC));
    Requires(stat[iC],Form("%s/%i/stat",particlename,iC));
    syst[iC] = (TH1D*)input_file->Get(Form("%s/%i/syst",particlename,iC));
    Requires(syst[iC],Form("%s/%i/syst",particlename,iC));
    syst_pt_uncorr[iC] = (TH1D*)input_file->Get(Form("%s/%i/syst_pt_uncorr",particlename,iC));
    Requires(syst_pt_uncorr[iC],Form("%s/%i/syst_pt_uncorr",particlename,iC));
    syst_pt_corr[iC] = (TH1D*)input_file->Get(Form("%s/%i/syst_pt_corr",particlename,iC));
    Requires(syst_pt_corr[iC],Form("%s/%i/syst_pt_corr",particlename,iC));
    syst_mult_corr[iC] = (TH1D*)input_file->Get(Form("%s/%i/syst_mult_corr",particlename,iC));
    Requires(syst_mult_corr[iC],Form("%s/%i/syst_mult_corr",particlename,iC));
    
    syst_mult_uncorr[iC] = (TH1D*)input_file->Get(Form("%s/%i/syst_mult_uncorr",particlename,iC));
    Requires(syst_mult_uncorr[iC],Form("%s/%i/syst_mult_uncorr",particlename,iC));
    TH1D *mStat = new TH1D(Form("mStat%d",iC),";#it{p}_{T} (GeV / #it{c}); 1/#it{N}_{ev} d^{2}#it{N}/d#it{p}_{T}d#it{y} (GeV/#it{c})^{-1}",stat[iC]->GetNbinsX(), stat[iC]->GetXaxis()->GetXbins()->GetArray());
    TH1D *mSystTot = new TH1D(Form("mSystTot%d",iC),";#it{p}_{T} (GeV/#it{c});1/#it{N}_{ev} d^{2}#it{N}/d#it{p}_{T}d#it{y} (GeV/#it{c})^{-1}",stat[iC]->GetNbinsX(), stat[iC]->GetXaxis()->GetXbins()->GetArray());
    TH1D *mSystPtUncorr = new TH1D(Form("mSystPtUncorr%d",iC),";#it{p}_{T} (GeV/#it{c});1/#it{N}_{ev} d^{2}#it{N}/d#it{p}_{T}d#it{y} (GeV/#it{c})^{-1}",stat[iC]->GetNbinsX(), stat[iC]->GetXaxis()->GetXbins()->GetArray());
    TH1D *mSystPtCorr = new TH1D(Form("mSystPtCorr%d",iC),";#it{p}_{T} (GeV/#it{c});1/#it{N}_{ev} d^{2}#it{N}/d#it{p}_{T}d#it{y} (GeV/#it{c})^{-1}",stat[iC]->GetNbinsX(), stat[iC]->GetXaxis()->GetXbins()->GetArray());
    TH1D *mSystMultCorr = new TH1D(Form("mSystMultCorr%d",iC),";#it{p}_{T} (GeV/#it{c});1/#it{N}_{ev} d^{2}#it{N}/d#it{p}_{T}d#it{y} (GeV/#it{c})^{-1}",stat[iC]->GetNbinsX(), stat[iC]->GetXaxis()->GetXbins()->GetArray());
    TH1D *mSystMultUncorr = new TH1D(Form("mSystMultUncorr%d",iC),";#it{p}_{T} (GeV/#it{c});1/#it{N}_{ev} d^{2}#it{N}/d#it{p}_{T}d#it{y} (GeV/#it{c})^{-1}",stat[iC]->GetNbinsX(), stat[iC]->GetXaxis()->GetXbins()->GetArray());

    for (int j = 1; j <= stat[iC]->GetNbinsX(); ++j) {
      double x = stat[iC]->GetBinCenter(j);
      if (x < 0.6 || x > kCentPtLimits[iC]) continue;
      float val = stat[iC]->FindBin(mStat->GetBinCenter(j));
      float stat_err = stat[iC]->GetBinError(val);
      float syst_err = syst[iC]->GetBinError(val);
      float syst_pt_uncorr_err = syst_pt_uncorr[iC]->GetBinError(val);
      float syst_pt_corr_err = syst_pt_corr[iC]->GetBinError(val);
      float syst_mult_corr_err = syst_mult_corr[iC]->GetBinError(val);
      float syst_mult_uncorr_err = syst_mult_uncorr[iC]->GetBinError(val);
      mStat->SetBinContent(j,stat[iC]->GetBinContent(val));
      mStat->SetBinError(j,stat_err);
      mSystTot->SetBinContent(j,stat[iC]->GetBinContent(val));
      mSystTot->SetBinError(j,syst_err);
      mSystPtUncorr->SetBinContent(j,stat[iC]->GetBinContent(val));
      mSystPtUncorr->SetBinError(j,syst_pt_uncorr_err);
      mSystPtCorr->SetBinContent(j,stat[iC]->GetBinContent(val));
      mSystPtCorr->SetBinError(j,syst_pt_corr_err);
      mSystMultCorr->SetBinContent(j,stat[iC]->GetBinContent(val));
      mSystMultCorr->SetBinError(j,syst_mult_corr_err);
      mSystMultUncorr->SetBinContent(j,stat[iC]->GetBinContent(val));
      mSystMultUncorr->SetBinError(j,syst_mult_uncorr_err);
      
    }
  
    fit_function->SetParameters(kParticleMass, n, C, normal);
    fit_function->SetParLimits(1, nMin, nMax);
    fit_function->SetParLimits(2, CMin, CMax);
    fit_function->SetParLimits(3, normMin, normMax);
    TF1* fout;
    results_old[iC] = yieldmean::YieldMean(mStat,mSystTot,fout,fit_function,0,10.,0.01,0.1,true,Form("%s/%soldchecklog.root",kBaseOutputDir.data(),particlename),Form("%d",iC));
    results_new[iC] = yieldmeannew::YieldMeanNew(mStat,mSystTot,mSystPtUncorr,mSystPtCorr,mSystMultCorr,fout,fit_function,0,10.,0.01,0.1,true,Form("%s/%schecklog.root",kBaseOutputDir.data(),particlename),Form("%d", iC));
    cout << "\n*****************************" << endl;
    printf("Levy-Tsallis\n");
    printf("iC: %d\n", iC);
    for (int iP = 0; iP < fout->GetNpar(); ++iP) {
      cout << fout->GetParName(iP) << ": " << fout->GetParameter(iP);
      cout << " +/- " << fout->GetParError(iP) << endl;
    }
    cout << "*****************************" << endl << endl;
    outputfile.mkdir(Form("%d",iC));
    outputfile.cd(Form("%d",iC));
    fout->Write(Form("Levy-Tsallis%d",iC));
    for(int iParam=0; iParam<4; iParam++){
      tsallis_param_val[iC][iParam] = fout->GetParameter(iParam);
      tsallis_param_err[iC][iParam] = fout->GetParError(iParam);
    }
    tsallis_param_val[iC][4] = fout->GetChisquare()/fout->GetNDF();
    tsallis_param_err[iC][4] = 0;
    results_new[iC]->Write(Form("result_new%d",iC));
    results_old[iC]->Write(Form("result_old%d",iC));
    datadir->mkdir(Form("%d",iC));
    datadir->cd(Form("%d",iC));
    mStat->Write(Form("stat%d",iC));
    mSystTot->Write(Form("syst%d",iC));
    mSystPtUncorr->Write(Form("syst_mult_uncorr%d",iC));
    mSystMultCorr->Write(Form("syst_mult_corr%d",iC));
    mSystMultUncorr->Write(Form("syst_mult_uncorr%d",iC));
  }

  // yields vs multiplicity

  float dNdEta_tmp, dNdEtaErr_tmp, proton_yields_tmp, proton_yields_stat_tmp, proton_yields_syst_tmp;

  vector<float> dNdEta_vec, dNdEtaErr_vec, proton_yields_vec, proton_yields_stat_vec, proton_yields_syst_vec;

  std::ifstream protonFile(Form("%s/proton_yields.txt",kBaseOutputDir.data()));
  if(!protonFile.is_open()){
    printf("The file %s could not be opened\n", Form("%s/proton_yields.txt",kBaseOutputDir.data()));
    exit(1);
  }
  else{
    while(protonFile >> dNdEta_tmp >> dNdEtaErr_tmp >> proton_yields_tmp >> proton_yields_stat_tmp >> proton_yields_syst_tmp){
      dNdEta_vec.push_back(dNdEta_tmp);
      dNdEtaErr_vec.push_back(dNdEtaErr_tmp);
      proton_yields_vec.push_back(proton_yields_tmp);
      proton_yields_stat_vec.push_back(proton_yields_stat_tmp);
      proton_yields_syst_vec.push_back(proton_yields_syst_tmp);
    }
  }

  double dNdEta[kCentLength-1], dNdEtaShift[kCentLength-1], dNdEtaErr[kCentLength-1], proton_yields[kCentLength-1], proton_yields_stat[kCentLength-1], proton_yields_syst[kCentLength-1];

  int iCent=0;
  for(int i=0; i< (int) dNdEta_vec.size();i++){
    if(i!=4){
      dNdEta[iCent] = dNdEta_vec[i];
      dNdEtaErr[iCent] = dNdEtaErr_vec[i];
      proton_yields[iCent] = proton_yields_vec[i];
      proton_yields_stat[iCent] = proton_yields_stat_vec[i];
      proton_yields_syst[iCent] = proton_yields_syst_vec[i];
      if (i!=3) iCent++;
    }
    else{
      dNdEta[iCent] = (dNdEta_vec[i]+dNdEta[iCent])/2;
      dNdEtaErr[iCent] = TMath::Sqrt(Sq(dNdEtaErr_vec[i]) + Sq(dNdEtaErr[iCent]))/2;
      proton_yields[iCent] = (proton_yields_vec[i]+proton_yields[iCent])/2;
      proton_yields_stat[iCent] = TMath::Sqrt(Sq(proton_yields_stat_vec[i]) + Sq(proton_yields_stat[iCent]))/2;
      proton_yields_syst[iCent] = TMath::Sqrt(Sq(proton_yields_syst_vec[i]) + Sq(proton_yields_syst[iCent]))/2;
      iCent++;
    }
  }
  for(int iC=0; iC<kCentLength; iC++){
    dNdEtaShift[iCent] = dNdEta[iCent] + 1;
  }
  // new approach
  double nucleus_mean_pt_new[kCentLength-1]{0.};
  double nucleus_mean_pt_stat_new[kCentLength-1]{0.};
  double nucleus_mean_pt_syst_new[kCentLength-1]{0.};
  double nucleus_mean_pt_syst_corr_new[kCentLength-1]{0.};
  double nucleus_mean_pt_syst_uncorr_new[kCentLength-1]{0.};
  double nucleus_yield_new[kCentLength-1]{0.};
  double nucleus_yield_stat_new[kCentLength-1]{0.};
  double nucleus_yield_syst_new[kCentLength-1]{0.};
  double nucleus_yield_syst_corr_new[kCentLength-1]{0.};
  double nucleus_yield_syst_uncorr_new[kCentLength-1]{0.};
  // old approach
  double nucleus_mean_pt_old[kCentLength-1]{0.};
  double nucleus_mean_pt_stat_old[kCentLength-1]{0.};
  double nucleus_mean_pt_syst_old[kCentLength-1]{0.};
  double nucleus_yield_old[kCentLength-1]{0.};
  double nucleus_yield_stat_old[kCentLength-1]{0.};
  double nucleus_yield_syst_old[kCentLength-1]{0.};
  for (int iC = 0; iC < kCentLength-1; ++iC) {
    ///////////////////////////////
    //////     new method  ////////
    ///////////////////////////////
    nucleus_yield_new[iC] = results_new[iC]->GetBinContent(yieldmeannew::kYield);
    nucleus_yield_stat_new[iC] = results_new[iC]->GetBinContent(yieldmeannew::kYieldStat);
    nucleus_yield_syst_new[iC] = results_new[iC]->GetBinContent(yieldmeannew::kYieldSysTot);
    nucleus_yield_syst_corr_new[iC] = std::sqrt(Sq(results_new[iC]->GetBinContent(yieldmeannew::kYieldSysHiCorr)) + Sq(results_new[iC]->GetBinContent(yieldmeannew::kYieldSysLoCorr)));
    nucleus_yield_syst_uncorr_new[iC] = std::sqrt(Sq(nucleus_yield_syst_new[iC]) - Sq(nucleus_yield_syst_corr_new[iC]));
    nucleus_mean_pt_new[iC] = results_new[iC]->GetBinContent(yieldmeannew::kMean);
    nucleus_mean_pt_stat_new[iC] = results_new[iC]->GetBinContent(yieldmeannew::kMeanStat);
    nucleus_mean_pt_syst_new[iC] = results_new[iC]->GetBinContent(yieldmeannew::kMeanSysTot);
    nucleus_mean_pt_syst_corr_new[iC] = std::sqrt(Sq(results_new[iC]->GetBinContent(yieldmeannew::kMeanSysHardCorr)) + Sq(results_new[iC]->GetBinContent(yieldmeannew::kMeanSysSoftCorr)));
    nucleus_mean_pt_syst_uncorr_new[iC] = std::sqrt(Sq(nucleus_mean_pt_syst_new[iC]) - Sq(nucleus_mean_pt_syst_corr_new[iC]));
    std::cout << "New method" << std::endl;
    std::cout << std::setprecision(3);
    std::cout << nucleus_yield_new[iC] << " $\\pm$ " << nucleus_yield_stat_new[iC] << " $\\pm$ " << nucleus_yield_syst_new[iC] << " ($\\pm$ " << nucleus_yield_syst_corr_new[iC] << ")\t& ";
    std::cout << nucleus_mean_pt_new[iC] << " $\\pm$ " << nucleus_mean_pt_stat_new[iC] << " $\\pm$ " << nucleus_mean_pt_syst_new[iC]<< " ($\\pm$ " << nucleus_mean_pt_syst_corr_new[iC] << ")" << std::endl;
    
    ///////////////////////////////
    //////     old method  ////////
    ///////////////////////////////
    nucleus_yield_old[iC] = results_old[iC]->GetBinContent(yieldmean::kYield);;
    nucleus_yield_stat_old[iC] = results_old[iC]->GetBinContent(yieldmean::kYieldStat);
    nucleus_yield_syst_old[iC] = std::sqrt(Sq(results_old[iC]->GetBinContent(yieldmean::kYieldSysHi)) + Sq(results_old[iC]->GetBinContent(yieldmean::kYieldSysLo)));
    nucleus_mean_pt_old[iC] = results_old[iC]->GetBinContent(yieldmean::kMean);
    nucleus_mean_pt_stat_old[iC] = results_old[iC]->GetBinContent(yieldmean::kMeanStat);
    nucleus_mean_pt_syst_old[iC] = std::sqrt(Sq(results_old[iC]->GetBinContent(yieldmean::kMeanSysHi)) + Sq(results_old[iC]->GetBinContent(yieldmean::kMeanSysLo)));
    //
    std::cout << "Old method" << std::endl;
    std::cout << std::setprecision(3);
    std::cout << nucleus_yield_old[iC] << " $\\pm$ " << nucleus_yield_stat_old[iC] << " $\\pm$ " << nucleus_yield_syst_old[iC] << "\t& ";
    std::cout << nucleus_mean_pt_old[iC] << " $\\pm$ " << nucleus_mean_pt_stat_old[iC] << " $\\pm$ " << nucleus_mean_pt_syst_old[iC] << std::endl;
  }

  ///////////////////////////////
  //////     new method  ////////
  ///////////////////////////////

  TCanvas cMeanPtNew("cMeanPtNew","cMeanPtNew");
  cMeanPtNew.SetBottomMargin(0.15);
  //
  TGraphErrors mean_pt_gr_stat_new(kCentLength-1,dNdEta,nucleus_mean_pt_new,dNdEtaErr,nucleus_mean_pt_stat_new);
  mean_pt_gr_stat_new.SetTitle(Form("%s ",particlename));
  mean_pt_gr_stat_new.GetYaxis()->SetTitle("#LT #it{p}_{T} #GT (GeV/c)");
  mean_pt_gr_stat_new.GetXaxis()->SetTitle("#LTd#it{N}_{ch}/d#it{#eta}#GT_{|#it{#eta}|<0.5}");
  mean_pt_gr_stat_new.SetMarkerColor(kOrange-3);
  mean_pt_gr_stat_new.SetLineColor(kOrange-3);
  mean_pt_gr_stat_new.SetMarkerStyle(20);
  mean_pt_gr_stat_new.SetFillStyle(0);
  //
  TGraphErrors mean_pt_gr_syst_new(kCentLength-1,dNdEta,nucleus_mean_pt_new,dNdEtaErr,nucleus_mean_pt_syst_new);
  mean_pt_gr_syst_new.SetTitle(Form("%s ",particlename));
  mean_pt_gr_syst_new.GetYaxis()->SetTitle("#LT #it{p}_{T} #GT (GeV/c)");
  mean_pt_gr_syst_new.GetXaxis()->SetTitle("#LTd#it{N}_{ch}/d#it{#eta}#GT_{|#it{#eta}|<0.5}");
  mean_pt_gr_syst_new.SetLineColor(kOrange-3);
  mean_pt_gr_syst_new.SetMarkerColor(kOrange-3);
  mean_pt_gr_syst_new.SetMarkerStyle(20);
  mean_pt_gr_syst_new.SetFillStyle(0);
  //
  TGraphErrors mean_pt_gr_syst_uncorr_new(kCentLength-1,dNdEta,nucleus_mean_pt_new,dNdEtaErr,nucleus_mean_pt_syst_uncorr_new);
  mean_pt_gr_syst_uncorr_new.SetTitle(Form("%s ",particlename));
  mean_pt_gr_syst_uncorr_new.GetYaxis()->SetTitle("#LT #it{p}_{T} #GT (GeV/c)");
  mean_pt_gr_syst_uncorr_new.GetXaxis()->SetTitle("#LTd#it{N}_{ch}/d#it{#eta}#GT_{|#it{#eta}|<0.5}");
  mean_pt_gr_syst_uncorr_new.SetLineColor(kOrange-3);
  mean_pt_gr_syst_uncorr_new.SetMarkerColor(kOrange-3);
  mean_pt_gr_syst_uncorr_new.SetMarkerStyle(20);
  mean_pt_gr_syst_uncorr_new.SetFillStyle(0);
  //
  TGraphErrors mean_pt_gr_syst_corr_new(kCentLength-1,dNdEta,nucleus_mean_pt_new,dNdEtaErr,nucleus_mean_pt_syst_corr_new);
  mean_pt_gr_syst_corr_new.SetTitle(Form("%s ",particlename));
  mean_pt_gr_syst_corr_new.GetYaxis()->SetTitle("#LT #it{p}_{T} #GT (GeV/c)");
  mean_pt_gr_syst_corr_new.GetXaxis()->SetTitle("#LTd#it{N}_{ch}/d#it{#eta}#GT_{|#it{#eta}|<0.5}");
  mean_pt_gr_syst_corr_new.SetLineColor(kOrange-3);
  mean_pt_gr_syst_corr_new.SetMarkerColor(kOrange-3);
  mean_pt_gr_syst_corr_new.SetMarkerStyle(20);
  mean_pt_gr_syst_corr_new.SetFillStyle(3003);
  mean_pt_gr_syst_corr_new.SetFillColor(kOrange-3);
  
  //
  cMeanPtNew.cd();
  mean_pt_gr_syst_new.Draw("AP2");
  mean_pt_gr_syst_corr_new.Draw("P2SAME");
  mean_pt_gr_stat_new.Draw("PSAME");

  TCanvas cYieldNew("cYieldNew","cYieldNew");
  cYieldNew.SetBottomMargin(0.15);
  TGraphErrors yield_gr_stat_new(kCentLength-1,dNdEta,nucleus_yield_new,dNdEtaErr,nucleus_yield_stat_new);
  yield_gr_stat_new.SetTitle(Form("%s ",particlename));
  yield_gr_stat_new.GetYaxis()->SetTitle("d#it{N}/d#it{y}");
  yield_gr_stat_new.GetXaxis()->SetTitle("#LTd#it{N}_{ch}/d#it{#eta}#GT_{|#it{#eta}|<0.5}");
  yield_gr_stat_new.SetMarkerColor(kOrange-3);
  yield_gr_stat_new.SetLineColor(kOrange-3);
  yield_gr_stat_new.SetMarkerStyle(20);
  yield_gr_stat_new.SetFillStyle(0);
  //
  TGraphErrors yield_gr_syst_new(kCentLength-1,dNdEta,nucleus_yield_new,dNdEtaErr,nucleus_yield_syst_new);
  yield_gr_syst_new.SetTitle(Form("%s ",particlename));
  yield_gr_syst_new.GetYaxis()->SetTitle("d#it{N}/d#it{y}");
  yield_gr_syst_new.GetXaxis()->SetTitle("#LTd#it{N}_{ch}/d#it{#eta}#GT_{|#it{#eta}|<0.5}");
  yield_gr_syst_new.SetLineColor(kOrange-3);
  yield_gr_syst_new.SetMarkerColor(kOrange-3);
  yield_gr_syst_new.SetMarkerStyle(20);
  yield_gr_syst_new.SetFillStyle(0);
  //
  TGraphErrors yield_gr_syst_uncorr_new(kCentLength-1,dNdEta,nucleus_yield_new,dNdEtaErr,nucleus_yield_syst_uncorr_new);
  yield_gr_syst_uncorr_new.SetTitle(Form("%s ",particlename));
  yield_gr_syst_uncorr_new.GetYaxis()->SetTitle("d#it{N}/d#it{y}");
  yield_gr_syst_uncorr_new.GetXaxis()->SetTitle("#LTd#it{N}_{ch}/d#it{#eta}#GT_{|#it{#eta}|<0.5}");
  yield_gr_syst_uncorr_new.SetLineColor(kOrange-3);
  yield_gr_syst_uncorr_new.SetMarkerColor(kOrange-3);
  yield_gr_syst_uncorr_new.SetMarkerStyle(20);
  yield_gr_syst_uncorr_new.SetFillStyle(0);
  //
  TGraphErrors yield_gr_syst_corr_new(kCentLength-1,dNdEta,nucleus_yield_new,dNdEtaErr,nucleus_yield_syst_corr_new);
  yield_gr_syst_corr_new.SetTitle(Form("%s ",particlename));
  yield_gr_syst_corr_new.GetYaxis()->SetTitle("d#it{N}/d#it{y}");
  yield_gr_syst_corr_new.GetXaxis()->SetTitle("#LTd#it{N}_{ch}/d#it{#eta}#GT_{|#it{#eta}|<0.5}");
  yield_gr_syst_corr_new.SetLineColor(kOrange-3);
  yield_gr_syst_corr_new.SetMarkerColor(kOrange-3);
  yield_gr_syst_corr_new.SetMarkerStyle(20);
  yield_gr_syst_corr_new.SetFillStyle(3003);
  yield_gr_syst_corr_new.SetFillColor(kOrange-3);
  //
  cYieldNew.cd();
  yield_gr_syst_new.Draw("AP2");
  yield_gr_syst_corr_new.Draw("P2SAME");
  yield_gr_stat_new.Draw("PSAME");

  TDirectory* new_dir = outputfile.mkdir("YieldMeanNew");
  new_dir->cd();
  cMeanPtNew.Write(Form("can_%s_meanpt_new",particlename));
  mean_pt_gr_stat_new.Write(Form("%s_meanpt_stat_new",particlename));
  mean_pt_gr_syst_new.Write(Form("%s_meanpt_syst_new",particlename));
  mean_pt_gr_syst_corr_new.Write(Form("%s_meanpt_syst_corr_new",particlename));
  mean_pt_gr_syst_uncorr_new.Write(Form("%s_meanpt_syst_uncorr_new",particlename));
  cYieldNew.Write(Form("can_%s_yield_new",particlename));
  yield_gr_stat_new.Write(Form("%s_yield_stat_new",particlename));
  yield_gr_syst_new.Write(Form("%s_yield_syst_new",particlename));
  yield_gr_syst_corr_new.Write(Form("%s_yield_syst_corr_new",particlename));
  yield_gr_syst_uncorr_new.Write(Form("%s_yield_syst_uncorr_new",particlename));

  ///////////////////////////////
  //////     old method  ////////
  ///////////////////////////////

  TCanvas cMeanPtOld("cMeanPtOld","cMeanPtOld");
  cMeanPtOld.SetBottomMargin(0.15);
  //
  TGraphErrors mean_pt_gr_stat_old(kCentLength-1,dNdEta,nucleus_mean_pt_old,dNdEtaErr,nucleus_mean_pt_stat_old);
  mean_pt_gr_stat_old.SetTitle(Form("%s ",particlename));
  mean_pt_gr_stat_old.GetYaxis()->SetTitle("#LT #it{p}_{T} #GT (GeV/c)");
  mean_pt_gr_stat_old.GetXaxis()->SetTitle("#LTd#it{N}_{ch}/d#it{#eta}#GT_{|#it{#eta}|<0.5}");
  mean_pt_gr_stat_old.SetMarkerColor(kBlack);
  mean_pt_gr_stat_old.SetLineColor(kBlack);
  mean_pt_gr_stat_old.SetMarkerStyle(20);
  mean_pt_gr_stat_old.SetFillStyle(0);
  //
  TGraphErrors mean_pt_gr_syst_old(kCentLength-1,dNdEta,nucleus_mean_pt_old,dNdEtaErr,nucleus_mean_pt_syst_old);
  mean_pt_gr_syst_old.SetTitle(Form("%s ",particlename));
  mean_pt_gr_syst_old.GetYaxis()->SetTitle("#LT #it{p}_{T} #GT (GeV/c)");
  mean_pt_gr_syst_old.GetXaxis()->SetTitle("#LTd#it{N}_{ch}/d#it{#eta}#GT_{|#it{#eta}|<0.5}");
  mean_pt_gr_syst_old.SetLineColor(kBlack);
  mean_pt_gr_syst_old.SetMarkerColor(kBlack);
  mean_pt_gr_syst_old.SetMarkerStyle(20);
  mean_pt_gr_syst_old.SetFillStyle(0);
  //
  cMeanPtOld.cd();
  mean_pt_gr_syst_old.Draw("AP2");
  mean_pt_gr_stat_old.Draw("PSAME");

  TCanvas cYieldOld("cYieldOld","cYieldOld");
  cYieldOld.SetBottomMargin(0.15);
  //
  TGraphErrors yield_gr_stat_old(kCentLength-1,dNdEta,nucleus_yield_old,dNdEtaErr,nucleus_yield_stat_old);
  yield_gr_stat_old.SetTitle(Form("%s ",particlename));
  yield_gr_stat_old.GetYaxis()->SetTitle("d#it{N}/d#it{y}");
  yield_gr_stat_old.GetXaxis()->SetTitle("#LTd#it{N}_{ch}/d#it{#eta}#GT_{|#it{#eta}|<0.5}");
  yield_gr_stat_old.SetMarkerColor(kBlack);
  yield_gr_stat_old.SetLineColor(kBlack);
  yield_gr_stat_old.SetMarkerStyle(20);
  yield_gr_stat_old.SetFillStyle(0);
  //
  TGraphErrors yield_gr_syst_old(kCentLength-1,dNdEta,nucleus_yield_old,dNdEtaErr,nucleus_yield_syst_old);
  yield_gr_syst_old.SetTitle(Form("%s ",particlename));
  yield_gr_syst_old.GetYaxis()->SetTitle("d#it{N}/d#it{y}");
  yield_gr_syst_old.GetXaxis()->SetTitle("#LTd#it{N}_{ch}/d#it{#eta}#GT_{|#it{#eta}|<0.5}");
  yield_gr_syst_old.SetLineColor(kBlack);
  yield_gr_syst_old.SetMarkerColor(kBlack);
  yield_gr_syst_old.SetMarkerStyle(20);
  yield_gr_syst_old.SetFillStyle(0);
  //
  cYieldOld.cd();
  yield_gr_syst_old.Draw("AP2");
  yield_gr_stat_old.Draw("PSAME");

  TDirectory* old_dir = outputfile.mkdir("YieldMean");
  old_dir->cd();
  cMeanPtOld.Write(Form("can_%s_meanpt_old",particlename));
  mean_pt_gr_stat_old.Write(Form("%s_meanpt_stat_old",particlename));
  mean_pt_gr_syst_old.Write(Form("%s_meanpt_syst_old",particlename));
  cYieldOld.Write(Form("can_%s_yield_old",particlename));
  yield_gr_stat_old.Write(Form("%s_yield_stat_old",particlename));
  yield_gr_syst_old.Write(Form("%s_yield_syst_old",particlename));

  TDirectory* comp_dir = outputfile.mkdir("Comparison");
  TCanvas cYieldComp("cYieldComp","cYieldComp");
  cYieldComp.SetBottomMargin(0.15);
  yield_gr_syst_old.Draw("AP2");
  yield_gr_stat_old.Draw("PSAME");
  yield_gr_syst_new.Draw("P2Same");
  yield_gr_syst_corr_new.Draw("P2SAME");
  yield_gr_stat_new.Draw("PSAME");
  TLegend leg_yield(0.7,0.2,0.9,0.4,"","brNDC");
  leg_yield.AddEntry(&yield_gr_syst_old, "YieldMean.C", "fp");
  leg_yield.AddEntry(&yield_gr_syst_new, "YieldMeanNew.cc", "fp");
  leg_yield.Draw();
  comp_dir->cd();
  cYieldComp.Write(Form("can_%s_yield_comp",particlename));

  TCanvas cMeanComp("cMeanComp","cMeanComp");
  cMeanComp.SetBottomMargin(0.15);
  mean_pt_gr_syst_old.Draw("AP2");
  mean_pt_gr_stat_old.Draw("PSAME");
  mean_pt_gr_syst_new.Draw("P2SAME");
  mean_pt_gr_syst_corr_new.Draw("P2SAME");
  mean_pt_gr_stat_new.Draw("PSAME");
  TLegend leg_mean(0.7,0.2,0.9,0.4,"","brNDC");
  leg_mean.AddEntry(&mean_pt_gr_syst_old, "YieldMean.C", "fp");
  leg_mean.AddEntry(&mean_pt_gr_syst_new, "YieldMeanNew.cc", "fp");
  leg_mean.Draw();
  comp_dir->cd();
  cMeanComp.Write(Form("can_%s_mean_pt_comp",particlename));

  outputfile.Close();
}
