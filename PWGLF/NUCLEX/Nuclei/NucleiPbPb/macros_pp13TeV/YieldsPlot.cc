#include "src/Common.h"
#include "src/Utils.h"
using namespace utils;
#include "src/YieldMean.C"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
using std::string;
#include <fstream>
#include <vector>
using std::vector;

#include <TFile.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TF1.h>
#include <TStyle.h>
#include <TGraphErrors.h>
#include <TLegend.h>

const int kKnownMult = 6;
//constexpr double dNdEta[kCentLength]    = {25.75, 19.83, 16.12, 12.78, 10.11, 8.07, 6.48, 4.64, 2.52, 6.84};
//constexpr double dNdEtaErr[kCentLength] = { 0.40,  0.30,  0.24,  0.14,  0.15, 0.12, 0.10, 0.07, 0.04, 0.01};
constexpr int kNfitFunctions = 3;
const string kFitFunctionNames[kNfitFunctions] = {"LevyTsallis", "Boltzmann", "Mt-exp"};
const string kCentralities[kCentLength-1] = {"0-1%","1-5%","5-10%","10-20%","20-30%","30-40%","40-50%","50-70%", "70-100%"};//, "0-100%"};

template<typename F> F Max(vector<F> &v) {
  double max = -1.e24;
  for (F& val : v) {
    if (val > max) max = val;
  }
  return max;
}

template<typename F> F Min(vector<F> &v) {
  double min = 1.e24;
  for (F& val : v) {
    if (val < min) min = val;
  }
  return min;
}

void YieldsPlot() {
  // /// Compute dNdEta in the bins used by the analysis
  // std::cout << "*** dN/deta in the centrality bins used" << std::endl;
  // double dNdEta[kCentLength]{0.};
  // double dNdEtaErr[kCentLength]{0.};
  // for (int iC = 0; iC < kCentLength; ++iC) {
  //   const int n = kCentBins[iC][1] - kCentBins[iC][0] + 1;
  //   for (int iI0 = kCentBins[iC][0]; iI0 <= kCentBins[iC][1]; ++iI0) {
  //     dNdEta[iC] += dNdEta[iI0-1];
  //     dNdEtaErr[iC] += dNdEtaErr[iI0-1];
  //   }
  //   dNdEta[iC] /= n;
  //   dNdEtaErr[iC] /= n;
  //   std::cout << kCentLabels[iC] << " " << dNdEta[iC] << " +/- " << dNdEtaErr[iC] << std::endl;
  // }

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

  double dNdEta[kCentLength-1], dNdEtaErr[kCentLength-1], proton_yields[kCentLength-1], proton_yields_stat[kCentLength-1], proton_yields_syst[kCentLength-1];

  int iC=0;
  for(int i=0; i< (int) dNdEta_vec.size();i++){
    if(i!=4){
      dNdEta[iC] = dNdEta_vec[i];
      dNdEtaErr[iC] = dNdEtaErr_vec[i];
      proton_yields[iC] = proton_yields_vec[i];
      proton_yields_stat[iC] = proton_yields_stat_vec[i];
      proton_yields_syst[iC] = proton_yields_syst_vec[i];
      if (i!=3) iC++;
    }
    else{
      dNdEta[iC] = (dNdEta_vec[i]+dNdEta[iC])/2;
      dNdEtaErr[iC] = TMath::Sqrt(Sq(dNdEtaErr_vec[i]) + Sq(dNdEtaErr[iC]))/2;
      proton_yields[iC] = (proton_yields_vec[i]+proton_yields[iC])/2;
      proton_yields_stat[iC] = TMath::Sqrt(Sq(proton_yields_stat_vec[i]) + Sq(proton_yields_stat[iC]))/2;
      proton_yields_syst[iC] = TMath::Sqrt(Sq(proton_yields_syst_vec[i]) + Sq(proton_yields_syst[iC]))/2;
      iC++;
    }
  }
  for(int ib=0; ib<kCentLength-1; ib++){
    printf("%f %f %f %f %f\n", dNdEta[ib],dNdEtaErr[ib],proton_yields[ib],proton_yields_stat[ib],proton_yields_syst[ib]);
  }

  std::cout << "*** Nuclei yield and mean pt" << std::endl;
  double nucleus_mean_pt[kCentLength-1]{0.};
  double nucleus_mean_pt_stat[kCentLength-1]{0.};
  double nucleus_mean_pt_syst[kCentLength-1]{0.};
  double nucleus_yield[kCentLength-1]{0.};
  double nucleus_yield_stat[kCentLength-1]{0.};
  double nucleus_yield_syst[kCentLength-1]{0.};
  for(int iS=0; iS<2; iS++){
    const char* filename = (iS==0) ? Form("%sfits.root",kBaseOutputDir.data()) : Form("%santifits.root",kBaseOutputDir.data());
    const char*  particlename = (iS==0) ? "Deuterons" : "Antideuterons";
    TFile nucleus(filename);
    printf("\n\n***************************************************************\n");
    printf("\t\t\t%s\n",particlename);
    printf("***************************************************************\n\n");
    for (int iC = 0; iC < kCentLength-1; ++iC) {
      vector<double> meanpts;
      vector<double> yields;
      float chi2{0.f};
      for (int iF = 0; iF < kNfitFunctions; ++iF) {
        TH1D* res = (TH1D*)nucleus.Get(Form("%s/%i/result%i",kFitFunctionNames[iF].data(),iC,iC));
        Requires(res, Form("%s/%i/result%i",kFitFunctionNames[iF].data(),iC,iC));
        if (res->GetBinContent(yieldmean::kFitRes) > 1.e-10) continue;
        yields.push_back(res->GetBinContent(yieldmean::kYield));
        meanpts.push_back(res->GetBinContent(yieldmean::kMean));
        if (iF==0) {
          nucleus_yield[iC] = yields[0];
          nucleus_yield_stat[iC] = res->GetBinContent(yieldmean::kYieldStat);
          nucleus_yield_syst[iC] = std::sqrt(Sq(res->GetBinContent(yieldmean::kYieldSysHi)) + Sq(res->GetBinContent(yieldmean::kYieldSysLo)));
          nucleus_mean_pt[iC] = meanpts[0];
          nucleus_mean_pt_stat[iC] = res->GetBinContent(yieldmean::kMeanStat);
          nucleus_mean_pt_syst[iC] = std::sqrt(Sq(res->GetBinContent(yieldmean::kMeanSysHi)) + Sq(res->GetBinContent(yieldmean::kMeanSysLo)));
          TF1* bw = (TF1*)nucleus.Get(Form("%s/%i/%s%i",kFitFunctionNames[iF].data(),iC,kFitFunctionNames[iF].data(),iC));
          Requires(bw,Form("%s/%i/%s%i",kFitFunctionNames[iF].data(),iC,kFitFunctionNames[iF].data(),iC));
          chi2 = bw->GetChisquare() / bw->GetNDF();
        }
      }
      nucleus_mean_pt_syst[iC] = std::sqrt(Sq(nucleus_mean_pt_syst[iC]) + Sq(0.5 * (Max(meanpts) - Min(meanpts))));
      nucleus_yield_syst[iC] = std::sqrt(Sq(nucleus_yield_syst[iC]) + Sq(0.5 * (Max(yields) - Min(yields))));
      std::cout << std::setprecision(3);
      std::cout << nucleus_yield[iC] << " $\\pm$ " << nucleus_yield_stat[iC] << " $\\pm$ " << nucleus_yield_syst[iC] << "\t& ";
      std::cout << nucleus_mean_pt[iC] << " $\\pm$ " << nucleus_mean_pt_stat[iC] << " $\\pm$ " << nucleus_mean_pt_syst[iC] << "\t& ";
      std::cout << chi2 << "\t\\\\" << std::endl;
    }

    /// Nucleus mean pt
    TFile mean_pt_file(Form("%s/%s_mean_pt.root",kBaseOutputDir.data(),particlename),"recreate");

    TCanvas cMeanPt("cMeanPt","cMeanPt");
    cMeanPt.SetBottomMargin(0.15);
    TGraphErrors mean_pt_gr_stat(kCentLength-1,dNdEta,nucleus_mean_pt,dNdEtaErr,nucleus_mean_pt_stat);
    mean_pt_gr_stat.SetTitle(Form("%s ",particlename));
    mean_pt_gr_stat.GetYaxis()->SetTitle("#LT #it{p}_{T} #GT (GeV/c)");
    mean_pt_gr_stat.GetXaxis()->SetTitle("#LTd#it{N}_{ch}/d#it{#eta}#GT_{|#it{#eta}|<0.5}");
    TGraphErrors mean_pt_gr_syst(kCentLength-1,dNdEta,nucleus_mean_pt,dNdEtaErr,nucleus_mean_pt_syst);
    mean_pt_gr_syst.SetTitle(Form("%s ",particlename));
    mean_pt_gr_syst.GetYaxis()->SetTitle("#LT #it{p}_{T} #GT (GeV/c)");
    mean_pt_gr_syst.GetXaxis()->SetTitle("#LTd#it{N}_{ch}/d#it{#eta}#GT_{|#it{#eta}|<0.5}");
    mean_pt_gr_stat.SetMarkerColor(kOrange-3);
    mean_pt_gr_stat.SetLineColor(kOrange-3);
    mean_pt_gr_stat.SetMarkerStyle(20);
    mean_pt_gr_stat.SetFillStyle(0);
    mean_pt_gr_syst.SetLineColor(kOrange-3);
    mean_pt_gr_syst.SetMarkerColor(kOrange-3);
    mean_pt_gr_syst.SetMarkerStyle(20);
    mean_pt_gr_syst.SetFillStyle(0);
    mean_pt_gr_stat.Write(Form("%s_meanpt_stat",particlename));
    mean_pt_gr_syst.Write(Form("%s_meanpt_syst",particlename));
    cMeanPt.cd();
    mean_pt_gr_syst.Draw("AP2");
    mean_pt_gr_stat.Draw("PSAME");

    TCanvas cYield("cYield","cYield");
    cYield.SetBottomMargin(0.15);
    TGraphErrors yield_gr_stat(kCentLength-1,dNdEta,nucleus_yield,dNdEtaErr,nucleus_yield_stat);
    TGraphErrors yield_gr_syst(kCentLength-1,dNdEta,nucleus_yield,dNdEtaErr,nucleus_yield_syst);
    yield_gr_syst.SetTitle(Form("%s ",particlename));
    yield_gr_syst.GetYaxis()->SetTitle("d#it{N}/d#it{y}");
    yield_gr_syst.GetXaxis()->SetTitle("#LTd#it{N}_{ch}/d#it{#eta}#GT_{|#it{#eta}|<0.5}");
    yield_gr_stat.SetMarkerColor(kOrange-3);
    yield_gr_stat.SetLineColor(kOrange-3);
    yield_gr_stat.SetMarkerStyle(20);
    yield_gr_stat.SetFillStyle(0);
    yield_gr_syst.SetLineColor(kOrange-3);
    yield_gr_syst.SetMarkerColor(kOrange-3);
    yield_gr_syst.SetMarkerStyle(20);
    yield_gr_syst.SetFillStyle(0);
    yield_gr_syst.Fit("pol1");
    yield_gr_stat.Write(Form("%s_yield_stat",particlename));
    yield_gr_syst.Write(Form("%s_yirld_syst",particlename));
    cYield.cd();
    yield_gr_syst.Draw("AP2");
    yield_gr_stat.Draw("PSAME");

    cMeanPt.Write(Form("can_%s_meanpt",particlename));
    cYield.Write(Form("can_%s_yield",particlename));
    mean_pt_file.Close();
  }

  /// Proton Yield
  std::cout << "*** Proton yields in the centrality bins used" << std::endl;

  for (int iC = 0; iC < kCentLength-1; ++iC) {
    std::cout << kCentralities[iC] << " " << proton_yields[iC] << " +/- " << proton_yields_stat[iC] << " +/- " << proton_yields_syst[iC] <<  std::endl;
  }

  /// Ratio Nucleus / p 2015
  vector<double> ratio(kCentLength-1,0.);
  vector<double> ratio_stat(kCentLength-1,0.);
  vector<double> ratio_syst(kCentLength-1,0.);
  std::cout << "*** d/p in the centrality bins used" << std::endl;
  for (int iC = 0; iC < kCentLength-1; ++iC) {
    printf("****************************************\n");
    printf("iC : %d\n", iC);
    printf("NucleusYied: %f +/- %f +/ %f\n", 2*nucleus_yield[iC], 2*nucleus_yield_stat[iC], 2*nucleus_yield_syst[iC]);
    printf("ProtonYied: %f +/- %f +/ %f\n", proton_yields[iC], proton_yields_stat[iC], proton_yields_syst[iC]);
    ratio[iC] = 2 * nucleus_yield[iC] / proton_yields[iC];
    ratio_stat[iC] = ratio[iC] * std::sqrt(Sq(nucleus_yield_stat[iC] / nucleus_yield[iC]) + Sq(proton_yields_stat[iC] / proton_yields[iC]));
    ratio_syst[iC] = ratio[iC] * std::sqrt(Sq(nucleus_yield_syst[iC] / nucleus_yield[iC]) + Sq(proton_yields_syst[iC] / proton_yields[iC]));
    std::cout /*<< kCentralities[iC] */<< " " << "ratio: "<< ratio[iC] << " +/- " << ratio_stat[iC] << " +/- " << ratio_syst[iC] <<  std::endl;
    printf("*****************************************\n");
  }

  // TCanvas* ratio_cv = new TCanvas("ratio_cv");
  // ratio_cv->DrawFrame(0.,0.,2000.,1.5 * Max(ratio),";#LTd#it{N}_{ch}/d#it{#eta}#GT_{|#it{#eta}|<0.5};2 d / (p + #bar{p})");
  // ratio_cv->SetLogx();
  // TFile f_prev(Form("%s/doverp.root",kBaseOutputDir.data()));
  // TCanvas* ratio_cv = (TCanvas*)f_prev.Get("c1");
  TFile f_prev(Form("%sdoverpPaperProp.root",kBaseOutputDir.data()));
  TCanvas* ratio_cv_input = (TCanvas*)f_prev.Get("c1_n19");
  TPad* pad = (TPad*)ratio_cv_input->GetPrimitive("c1_n19_1");
  //pad->GetListOfPrimitives()->ls();
  //TLegend* legPrel = (TLegend*)pad->GetPrimitive("legPrel");
  TFile f_doverp(Form("%sfinal_doverp.root",kBaseOutputDir.data()),"recreate");
  TCanvas* ratio_cv = new TCanvas("cDopverp","cDopverp");
  TGraphErrors* ratio_gr_stat = new TGraphErrors(kCentLength-1,dNdEta,ratio.data(),0,ratio_stat.data());
  ratio_gr_stat->SetName("pp13TeVstat");
  TGraphErrors* ratio_gr_syst = new TGraphErrors(kCentLength-1,dNdEta,ratio.data(),dNdEtaErr,ratio_syst.data());
  ratio_gr_syst->SetName("pp13TeVsyst");
  // legPrel->AddEntry(ratio_gr_syst,"ALICE, pp, #sqrt{s} = 13 TeV","PF");
  // legPrel->AddEntry((TObject*)nullptr,"V0M Multiplicity Classes","");
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  ratio_gr_syst->SetFillStyle(0);
  ratio_gr_syst->SetMarkerStyle(20);
  ratio_gr_stat->SetMarkerStyle(20);
  ratio_gr_syst->SetMarkerColor(kOrange-3);
  ratio_gr_stat->SetMarkerColor(kOrange-3);
  ratio_gr_syst->SetLineColor(kOrange-3);
  ratio_gr_stat->SetLineColor(kOrange-3);
  pad->cd();
  cout << "PadName: " << gPad->GetName() << endl;
  ratio_gr_stat->Draw("samepz");
  ratio_gr_syst->Draw("samep2");
  pad->Modified();
  pad->Update();
  printf("************************************************\n");
  //pad->GetListOfPrimitives()->ls();
  ratio_cv->cd();
  pad->Draw();
  ratio_cv->Write("cDoverp");
  ratio_gr_stat->Write("ratio_gr_stat");
  ratio_gr_syst->Write("ratio_gr_syst");
  TCanvas* cMio = new TCanvas("hMio","hMio");
  cMio->cd();
  ratio_gr_stat->Draw("apz");
  ratio_gr_syst->Draw("samep2");
  cMio->SaveAs(Form("%smiodoverp.C",kBaseOutputDir.data()));
  cMio->Write("myCanvas");
  //cMio->GetListOfPrimitives()->ls();
}
