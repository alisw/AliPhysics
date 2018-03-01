#include "src/Common.h"
#include "src/Utils.h"
using namespace utils;
#include "YieldMean.C"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
using std::string;

#include <TFile.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TF1.h>
#include <TStyle.h>

constexpr double kdNdEta[10]    = {1942.5, 1585.5,   1180,    786,    512,    318,    183,   96.3,   44.9, 17.52};
constexpr double kdNdEtaErr[10] = {  53.5,   46.0,    31.,    20.,    15.,    12.,     8.,    5.8,    3.4,  1.84};
constexpr int kNfitFunctions = 4;
const string kFitFunctionNames[kNfitFunctions] = {"BlastWave", "Boltzmann", "LevyTsallis", "Mt-exp"};
//const string kCentLabels[kCentLength] = {"0-1%","1-5%","5-10%","10-20%","20-30%","30-40%","40-60%","60-80%", "80-100%", "0-100%"};

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
  //     dNdEta[iC] += kdNdEta[iI0-1];
  //     dNdEtaErr[iC] += kdNdEtaErr[iI0-1];
  //   }
  //   dNdEta[iC] /= n;
  //   dNdEtaErr[iC] /= n;
  //   std::cout << kCentLabels[iC] << " " << dNdEta[iC] << " +/- " << dNdEtaErr[iC] << std::endl;
  // }

  /// Get He3 results
  std::cout << "*** Nuclei yield and mean pt" << std::endl;
  TFile nucleus("fits.root");
  double nucleus_mean_pt[kCentLength]{0.};
  double nucleus_mean_pt_stat[kCentLength]{0.};
  double nucleus_mean_pt_syst[kCentLength]{0.};
  double nucleus_yield[kCentLength]{0.};
  double nucleus_yield_stat[kCentLength]{0.};
  double nucleus_yield_syst[kCentLength]{0.};
  for (int iC = 0; iC < kCentLength; ++iC) {
    vector<double> meanpts;
    vector<double> yields;
    float chi2{0.f};
    for (int iF = 0; iF < kNfitFunctions; ++iF) {
      TH1D* res = (TH1D*)nucleus.Get(Form("%s/result%i",kFitFunctionNames[iF].data(),iC));
      Requires(res, Form("%s/result%i",kFitFunctionNames[iF].data(),iC));
      if (res->GetBinContent(kFitRes) > 1.e-10) continue;
      yields.push_back(res->GetBinContent(kYield));
      meanpts.push_back(res->GetBinContent(kMean));
      if (iF!=2) {
        nucleus_yield[iC] = yields[2];
        nucleus_yield_stat[iC] = res->GetBinContent(kYieldStat);
        nucleus_yield_syst[iC] = std::sqrt(Sq(res->GetBinContent(kYieldSysHi)) + Sq(res->GetBinContent(kYieldSysLo)));
        nucleus_mean_pt[iC] = meanpts[2];
        nucleus_mean_pt_stat[iC] = res->GetBinContent(kMeanStat);
        nucleus_mean_pt_syst[iC] = std::sqrt(Sq(res->GetBinContent(kMeanSysHi)) + Sq(res->GetBinContent(kMeanSysLo)));
        TF1* bw = (TF1*)nucleus.Get(Form("%s/%s%i",kFitFunctionNames[iF].data(),kFitFunctionNames[iF].data(),iC));
        Requires(bw,Form("%s/%s%i",kFitFunctionNames[iF].data(),kFitFunctionNames[iF].data(),iC));
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
}
//
//   /// Nucleus mean pt
//   TFile mean_pt_file("he3_mean_pt.root","recreate");
//   TGraphErrors mean_pt_gr_stat(kCentLength,dNdEta,nucleus_mean_pt,dNdEtaErr,nucleus_mean_pt_stat);
//   TGraphErrors mean_pt_gr_syst(kCentLength,dNdEta,nucleus_mean_pt,dNdEtaErr,nucleus_mean_pt_syst);
//   mean_pt_gr_stat.SetLineColor(kBlue);
//   mean_pt_gr_stat.SetMarkerColor(kBlue);
//   mean_pt_gr_stat.SetMarkerStyle(20);
//   mean_pt_gr_stat.SetFillStyle(0);
//   mean_pt_gr_syst.SetLineColor(kBlue);
//   mean_pt_gr_syst.SetMarkerColor(kBlue);
//   mean_pt_gr_syst.SetMarkerStyle(20);
//   mean_pt_gr_syst.SetFillStyle(0);
//   mean_pt_gr_stat.Write("deuterons_meanpt_stat");
//   mean_pt_gr_syst.Write("deuterons_meanpt_syst");
//   mean_pt_file.Close();
//
//
//   /// Proton Yield
//   std::cout << "*** Proton yields in the centrality bins used" << std::endl;
//   TFile proton("YieldAndMeanPtSumPwSys.root");
//   TH1F* proton_yield_hist_stat = (TH1F*)proton.Get("hStdYieldSummedProton");
//   TH1F* proton_yield_hist_syst = (TH1F*)proton.Get("hStdYieldSysSummedProton");
//   Requires(proton_yield_hist_stat,"hStdYieldSummedProton");
//   Requires(proton_yield_hist_syst,"hStdYieldSysSummedProton");
//   double proton_yields[kCentLength]{0.};
//   double proton_yields_stat[kCentLength]{0.};
//   double proton_yields_syst[kCentLength]{0.};
//   for (int iC = 0; iC < kCentLength; ++iC) {
//     const int n = kCentBins[iC][1] - kCentBins[iC][0] + 1;
//     for (int iI0 = kCentBins[iC][0]; iI0 <= kCentBins[iC][1]; ++iI0) {
//       proton_yields[iC] += proton_yield_hist_stat->GetBinContent(iI0);
//       proton_yields_stat[iC] += Sq(proton_yield_hist_stat->GetBinError(iI0));
//       proton_yields_syst[iC] += Sq(proton_yield_hist_syst->GetBinError(iI0));
//     }
//     proton_yields[iC] /= n;
//     proton_yields_stat[iC] = std::sqrt(proton_yields_stat[iC]) / n;
//     proton_yields_syst[iC] = std::sqrt(proton_yields_syst[iC]) / n;
//     std::cout << kCentLabels[iC] << " " << proton_yields[iC] << " +/- " << proton_yields_stat[iC] << " +/- " << proton_yields_syst[iC] <<  std::endl;
//   }
//
//   /// Ratio Nucleus / p 2015
//   vector<double> ratio(kCentLength,0.);
//   vector<double> ratio_stat(kCentLength,0.);
//   vector<double> ratio_syst(kCentLength,0.);
//   for (int iC = 0; iC < kCentLength; ++iC) {
//     ratio[iC] = 2 * nucleus_yield[iC] / proton_yields[iC];
//     ratio_stat[iC] = ratio[iC] * std::sqrt(Sq(nucleus_yield_stat[iC] / nucleus_yield[iC]) + Sq(proton_yields_stat[iC] / proton_yields[iC]));
//     ratio_syst[iC] = ratio[iC] * std::sqrt(Sq(nucleus_yield_syst[iC] / nucleus_yield[iC]) + Sq(proton_yields_syst[iC] / proton_yields[iC]));
//     std::cout << kCentLabels[iC] << " " << ratio[iC] << " +/- " << ratio_stat[iC] << " +/- " << ratio_syst[iC] <<  std::endl;
//   }
//
//   TCanvas* ratio_cv = new TCanvas("ratio_cv");
//   ratio_cv->DrawFrame(0.,0.,2000.,1.5 * Max(ratio),";#LTd#it{N}_{ch}/d#it{#eta}#GT_{|#it{#eta}|<0.5};2 ^{3}He / (p + #bar{p})");
//   TGraphErrors* ratio_gr_stat = new TGraphErrors(kCentLength,dNdEta,ratio.data(),0,ratio_stat.data());
//   TGraphErrors* ratio_gr_syst = new TGraphErrors(kCentLength,dNdEta,ratio.data(),dNdEtaErr,ratio_syst.data());
//   ratio_gr_syst->SetFillStyle(0);
//   ratio_gr_syst->SetMarkerStyle(20);
//   ratio_gr_stat->SetMarkerStyle(20);
//   ratio_gr_syst->SetMarkerColor(kBlack);
//   ratio_gr_stat->SetMarkerColor(kBlack);
//   ratio_gr_syst->SetLineColor(kBlack);
//   ratio_gr_stat->SetLineColor(kBlack);
//   ratio_gr_stat->Draw("pz");
//   ratio_gr_syst->Draw("p2");
//
//   /// old results
//   double old_x[] = { 1206.7, 266.0 };
//   double old_xerr[] = { 46.67, 9.83};
//   double old_y[] = { 1.0611e-5, 0.8341e-5 };
//   double old_syst[] = { 0.2594952986086646e-5, 0.23774057289406872e-5 };
//   double old_stat[] = { 0.0359e-5, 0.0403e-5 };
//   TGraphErrors* old_gr_stat = new TGraphErrors(2,old_x,old_y,0,old_stat);
//   TGraphErrors* old_gr_syst = new TGraphErrors(2,old_x,old_y,old_xerr,old_syst);
//   old_gr_syst->SetFillStyle(0);
//   old_gr_syst->SetMarkerStyle(21);
//   old_gr_stat->SetMarkerStyle(21);
//   old_gr_syst->SetMarkerColor(kRed);
//   old_gr_stat->SetMarkerColor(kRed);
//   old_gr_syst->SetLineColor(kRed);
//   old_gr_stat->SetLineColor(kRed);
//   old_gr_stat->Draw("pz");
//   old_gr_syst->Draw("p2");
//   TLegend* antani = new TLegend(0.154135,0.161739,0.845865,0.391304);
//   antani->SetMargin(0.1);
//   antani->SetHeader("ALICE Pb-Pb");
//   antani->SetFillStyle(0);
//   antani->AddEntry(ratio_gr_syst,"#sqrt{#it{s}_{NN}} = 5.02 TeV, Preliminary","p");
//   antani->AddEntry(old_gr_syst,"#sqrt{#it{s}_{NN}} = 2.76 TeV, PRC 93 (2016) 2, 024917","p");
//   antani->Draw();
//
//
//   /// deuterons
//   std::cout << "*** Deuteron yields in the centrality bins used" << std::endl;
//   TFile deuteron("deuteron_yields.root");
//   TH1F* deuteron_yield_hist_stat = (TH1F*)deuteron.Get("deuteron_yield_stat");
//   TH1F* deuteron_yield_hist_syst = (TH1F*)deuteron.Get("deuteron_yield_syst");
//   Requires(deuteron_yield_hist_stat,"deuteron_yield_stat");
//   Requires(deuteron_yield_hist_syst,"deuteron_yield_syst");
//   double deuteron_yields[kCentLength]{0.};
//   double deuteron_yields_stat[kCentLength]{0.};
//   double deuteron_yields_syst[kCentLength]{0.};
//   for (int iC = 0; iC < kCentLength; ++iC) {
//     const int n = kCentBins[iC][1] - kCentBins[iC][0] + 1;
//     for (int iI0 = kCentBins[iC][0]; iI0 <= kCentBins[iC][1]; ++iI0) {
//       deuteron_yields[iC] += deuteron_yield_hist_stat->GetBinContent(iI0);
//       deuteron_yields_stat[iC] += Sq(deuteron_yield_hist_stat->GetBinError(iI0));
//       deuteron_yields_syst[iC] += Sq(deuteron_yield_hist_syst->GetBinError(iI0));
//     }
//     deuteron_yields[iC] /= n;
//     deuteron_yields_stat[iC] = std::sqrt(deuteron_yields_stat[iC]) / n;
//     deuteron_yields_syst[iC] = std::sqrt(deuteron_yields_syst[iC]) / n;
//     std::cout << kCentLabels[iC] << " " << deuteron_yields[iC] << " +/- " << deuteron_yields_stat[iC] << " +/- " << deuteron_yields_syst[iC] <<  std::endl;
//   }
//
//   /// Mass
//   gStyle->SetOptFit(1);
//   TCanvas* yields_cv = new TCanvas("yields_cv");
//   yields_cv->DrawFrame(0.72,0.9e-4,3.12,100,";m (GeV/#it{c}^{2});d#it{N}/d#it{y}");
//   yields_cv->SetLogy();
//   double mass[3]{9.38271999359130859e-01,1.87561297416687012e+00,2.80923008918762207e+00};
//   double yields[3]{proton_yields[0] / 2.,deuteron_yields[0],nucleus_yield[0]};
//   double yields_syst[3]{
//     std::sqrt(Sq(proton_yields_syst[0]) + Sq(proton_yields_stat[0])) / 2. ,
//     std::sqrt(Sq(deuteron_yields_syst[0]) + Sq(deuteron_yields_stat[0])),
//     std::sqrt(Sq(nucleus_yield_syst[0]) + Sq(nucleus_yield_stat[0]))};
//   double yields_stat[3]{proton_yields_stat[0] / 2. ,deuteron_yields_stat[0],nucleus_yield_stat[0]};
//   TGraphErrors* yields_gr_stat = new TGraphErrors(3,mass,yields,0x0,yields_stat);
//   TGraphErrors* yields_gr_syst = new TGraphErrors(3,mass,yields,0x0,yields_syst);
//   yields_gr_syst->SetMarkerStyle(21);
//   yields_gr_stat->SetMarkerStyle(21);
//   yields_gr_syst->SetMarkerColor(kBlue);
//   yields_gr_stat->SetMarkerColor(kBlue);
//   yields_gr_syst->SetLineColor(kBlue);
//   yields_gr_stat->SetLineColor(kBlue);
//   yields_gr_stat->Draw("p");
//   yields_gr_syst->Draw("[]");
//   TF1 *f = new TF1("my_exp","[0]*exp(-x/[1])",0,3.2);
//   f->SetParNames("A","T");
//   f->SetParLimits(1,0.1,0.2);
//   auto ptr = yields_gr_syst->Fit(f,"S");
//   cout << "Chi2/ndf: " << ptr.Get()->Chi2() / ptr.Get()->Ndf() <<endl;
//   cout << "Tchem: " << f->GetParameter(1) << " +/- " << f->GetParError(1) <<endl;
//   f->Draw("same");
//   TLatex l;
//   l.SetTextFont(42);
//   l.DrawLatexNDC(.35,.83,"This work, 0-10% Pb-Pb #sqrt{#it{s}_{NN}}=5.02 TeV");
//   l.DrawLatexNDC(.15,.2,Form("#splitline{#it{T} = %3.1f #pm %1.1f MeV/#it{c}}{#chi^{2}/NDF = %1.3f}",f->GetParameter(1)*1000,f->GetParError(1)*1000,ptr.Get()->Chi2() / ptr.Get()->Ndf()));
//   TLegend *leg = new TLegend(0.6,0.66,0.92,0.79);
//   leg->SetLineWidth(0);
//   leg->SetFillStyle(0);
//   leg->AddEntry(yields_gr_syst, "Data","lp");
//   leg->AddEntry(f,"Fit function Ae^{-m/#it{T}}","l");
//   leg->Draw();
//
//   /// Ratios over p
//   double doverp = nucleus_yield[0] / deuteron_yields[0];
//   double doverp_stat = doverp * std::sqrt(Sq(deuteron_yields_stat[0]/deuteron_yields[0]) + Sq(nucleus_yield_stat[0] / nucleus_yield[0]));
//   double doverp_syst = doverp * std::sqrt(Sq(deuteron_yields_syst[0]/deuteron_yields[0]) + Sq(nucleus_yield_syst[0] / nucleus_yield[0]));
//   std::cout << "he3/d " << doverp << " +/- " << doverp_stat << " +/- " << doverp_syst << std::endl;
//   std::cout << "He3/p " << ratio[0] << " +/- " << ratio_stat[0] << " +/- " << ratio_syst[0] << std::endl;
//
// }
