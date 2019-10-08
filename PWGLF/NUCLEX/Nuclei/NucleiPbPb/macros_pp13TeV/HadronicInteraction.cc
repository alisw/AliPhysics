#include "src/Common.h"
#include "src/Utils.h"
using namespace utils;
#include "src/Plotting.h"
using namespace plotting;
#include "src/FitModules.h"

#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <algorithm>

void HadronicInteraction(){
  TFile signalfile(kSignalOutput.data());
  TFile signalMCfile(kSignalMCOutput.data());

  TFile output(Form("%sHadronicInteraction.root",kBaseOutputDir.data()),"recreate");

  TH1D* hTPC[2];
  TH1D* hTPCMC[2];
  TH1D* hRawCounts[2];
  TH1D* hRawCountsMC[2];

  TH1D* fraction[2];
  TH1D* fractionMC[2];

  std::vector<float> correction[2];
  float mean_correction[2] = {0.};
  float syst_correction[2] = {0.};

  TCanvas* cFraction[2];

  const int nFracBins = 3;
  double frac_bins[nFracBins+1] = {0.9,1.0,1.1,1.2};//,1.4};

  TTList* list = (TTList*)signalfile.Get(kFilterListNames.data());
  TTList* listmc = (TTList*)signalMCfile.Get(kFilterListNamesMCasData.data());
  for(int iS=0; iS<2; iS++){
    cFraction[iS] = new TCanvas(Form("cFraction_%c",kLetter[iS]),Form("cFraction_%c",kLetter[iS]));
    //
    hRawCounts[iS] = (TH1D*)signalfile.Get(Form("%s/%s/Fits/hRawCounts%c%i",list->GetName(),kNames[iS].data(),kLetter[iS],9));
    Requires(hRawCounts[iS],Form("%s/%s/Fits/hRawCounts%c%i",list->GetName(),kNames[iS].data(),kLetter[iS],9));
    hRawCounts[iS]->SetDirectory(0);
    hTPC[iS] = (TH1D*)signalfile.Get(Form("%s/%s/TPConly/hRawCountsTPC%c%i",list->GetName(),kNames[iS].data(),kLetter[iS],9));
    Requires(hTPC[iS],Form("%s/%s/TPConly/hRawCountsTPC%c%i",list->GetName(),kNames[iS].data(),kLetter[iS],9));
    hTPC[iS]->SetDirectory(0);
    hRawCountsMC[iS] = (TH1D*)signalMCfile.Get(Form("%s/%s/Fits/hRawCounts%c%i",listmc->GetName(),kNames[iS].data(),kLetter[iS],9));
    Requires(hRawCountsMC[iS],Form("MC: %s/%s/Fits/hRawCounts%c%i",listmc->GetName(),kNames[iS].data(),kLetter[iS],9));
    hRawCountsMC[iS]->SetDirectory(0);
    hTPCMC[iS] = (TH1D*)signalMCfile.Get(Form("%s/%s/TPConly/hRawCountsTPC%c%i",listmc->GetName(),kNames[iS].data(),kLetter[iS],9));
    Requires(hTPCMC[iS],Form("MC: %s/%s/TPConly/hRawCountsTPC%c%i",listmc->GetName(),kNames[iS].data(),kLetter[iS],9));
    hTPCMC[iS]->SetDirectory(0);
    //
    fraction[iS] = new TH1D(Form("hFraction%c",kLetter[iS]),";#it{p}_{T} (GeV/#it{c}); TOF + TPC counts / TPC counts",nFracBins,frac_bins);
    SetHistStyle(fraction[iS],kBlue);
    fraction[iS]->GetYaxis()->SetRangeUser(0.,1.);
    fractionMC[iS] = new TH1D(Form("hFractionMC%c",kLetter[iS]),";#it{p}_{T} (GeV/#it{c}); TOF + TPC counts / TPC counts",nFracBins,frac_bins);
    fractionMC[iS]->GetYaxis()->SetRangeUser(0.,1.);
    SetHistStyle(fractionMC[iS],kRed);
    //
    TLegend leg(0.72,0.22,0.92,0.35);
    leg.SetBorderSize(0);
    leg.AddEntry(fraction[iS],"Data","PE");
    leg.AddEntry(fractionMC[iS],"MC","PE");
    float ratio=0., ratio_err=0., ratio_mc=0., ratio_mc_err=0.;
    for(int iB=1; iB<=kNPtBins; iB++){
      float bin_center = (kPtBins[iB]+kPtBins[iB-1])/2;
      if (bin_center < kPtHadronicInteractionRange[0] || bin_center > kPtHadronicInteractionRange[1]) continue;
      //
      ratio = hRawCounts[iS]->GetBinContent(iB)/hTPC[iS]->GetBinContent(iB);
      ratio_err = ratio * TMath::Sqrt( Sq(hRawCounts[iS]->GetBinError(iB)/hRawCounts[iS]->GetBinContent(iB)) + Sq(hTPC[iS]->GetBinError(iB)/hTPC[iS]->GetBinContent(iB)) );
      fraction[iS]->SetBinContent(fraction[iS]->FindBin(hRawCounts[iS]->GetBinCenter(iB)),ratio);
      fraction[iS]->SetBinError(fraction[iS]->FindBin(hRawCounts[iS]->GetBinCenter(iB)),ratio_err);
      //
      ratio_mc = hRawCountsMC[iS]->GetBinContent(iB)/hTPCMC[iS]->GetBinContent(iB);
      ratio_mc_err = ratio_mc * TMath::Sqrt( Sq(hRawCountsMC[iS]->GetBinError(iB)/hRawCountsMC[iS]->GetBinContent(iB)) + Sq(hTPCMC[iS]->GetBinError(iB)/hTPCMC[iS]->GetBinContent(iB)) );
      fractionMC[iS]->SetBinContent(fractionMC[iS]->FindBin(hRawCountsMC[iS]->GetBinCenter(iB)),ratio_mc);
      fractionMC[iS]->SetBinError(fractionMC[iS]->FindBin(hRawCountsMC[iS]->GetBinCenter(iB)),ratio_mc_err);
      //
      correction[iS].push_back(ratio/ratio_mc);
      printf("bin_center: %f matter: %c ratio: %f ratio_mc: %f ratio/ratio_mc: %f\n", bin_center, kLetter[iS], ratio, ratio_mc, ratio/ratio_mc);
    }
    for(auto p : correction[iS]){
      printf("corr_%c: %f\n",kLetter[iS],p);
      mean_correction[iS] += p;
    }
    mean_correction[iS]/=nFracBins;
    syst_correction[iS] = (*std::max_element(correction[iS].begin(),correction[iS].end()) - *std::min_element(correction[iS].begin(),correction[iS].end()))/2/mean_correction[iS];
    printf("State: %c correction: %f systematic: %f \n",kLetter[iS], mean_correction[iS], syst_correction[iS]);
    output.cd();
    hRawCountsMC[iS]->Write(Form("TofMC_%c",kLetter[iS]));
    hRawCounts[iS]->Write(Form("Tof_%c",kLetter[iS]));
    hTPC[iS]->Write(Form("Tpc_%c",kLetter[iS]));
    hTPCMC[iS]->Write(Form("TpcMC_%c",kLetter[iS]));
    fraction[iS]->Write();
    fractionMC[iS]->Write();
    cFraction[iS]->cd();
    fraction[iS]->Draw();
    fractionMC[iS]->Draw("same");
    leg.Draw();
    cFraction[iS]->Write();
  }
}
