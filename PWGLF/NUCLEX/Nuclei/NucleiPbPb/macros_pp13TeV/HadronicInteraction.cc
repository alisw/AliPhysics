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
  TFile signalMCfile(Form("%ssignalMC.root",kBaseOutputDir.data()));

  TFile output(Form("%sHadronicInteraction.root",kBaseOutputDir.data()),"recreate");

  TH1D* hTPC3sigma[2];
  TH1D* hTPC3sigmaMC[2];
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
  TTList* listmc = (TTList*)signalMCfile.Get(kFilterListNames.data());
  for(int iS=0; iS<2; iS++){
    cFraction[iS] = new TCanvas(Form("cFraction_%c",kLetter[iS]),Form("cFraction_%c",kLetter[iS]));
    //
    hRawCounts[iS] = (TH1D*)signalfile.Get(Form("%s/%s/Fits/hRawCounts%c%i",list->GetName(),kNames[iS].data(),kLetter[iS],9));
    Requires(hRawCounts[iS],Form("%s/%s/Fits/hRawCounts%c%i",list->GetName(),kNames[iS].data(),kLetter[iS],9));
    hRawCounts[iS]->SetDirectory(0);
    hTPC3sigma[iS] = (TH1D*)signalfile.Get(Form("%s/%s/TPC3sigma/hTPC3sigma%c%i",list->GetName(),kNames[iS].data(),kLetter[iS],9));
    Requires(hTPC3sigma[iS],Form("%s/%s/TPC3sigma/hTPC3sigms%c%i",list->GetName(),kNames[iS].data(),kLetter[iS],9));
    hTPC3sigma[iS]->SetDirectory(0);
    hRawCountsMC[iS] = (TH1D*)signalMCfile.Get(Form("%s/%s/Fits/hRawCounts%c%i",list->GetName(),kNames[iS].data(),kLetter[iS],9));
      Requires(hRawCountsMC[iS],Form("%s/%s/Fits/hRawCountsMC%c%i",list->GetName(),kNames[iS].data(),kLetter[iS],9));
    hRawCountsMC[iS]->SetDirectory(0);
    hTPC3sigmaMC[iS] = (TH1D*)signalMCfile.Get(Form("%s/%s/TPC3sigma/hTPC3sigma%c%i",list->GetName(),kNames[iS].data(),kLetter[iS],9));
    Requires(hTPC3sigmaMC[iS],Form("MC: %s/%s/TPC3sigma/hTPC3sigma%c%i",list->GetName(),kNames[iS].data(),kLetter[iS],9));
    hTPC3sigmaMC[iS]->SetDirectory(0);
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
    for(int iB=5; iB<=7; iB++){
      //
      ratio = hRawCounts[iS]->GetBinContent(iB)/hTPC3sigma[iS]->GetBinContent(iB);
      ratio_err = ratio * TMath::Sqrt( Sq(hRawCounts[iS]->GetBinError(iB)/hRawCounts[iS]->GetBinContent(iB)) + Sq(hTPC3sigma[iS]->GetBinError(iB)/hTPC3sigma[iS]->GetBinContent(iB)) );
      fraction[iS]->SetBinContent(fraction[iS]->FindBin(hRawCounts[iS]->GetBinCenter(iB)),ratio);
      fraction[iS]->SetBinError(fraction[iS]->FindBin(hRawCounts[iS]->GetBinCenter(iB)),ratio_err);
      //
      ratio_mc = hRawCountsMC[iS]->GetBinContent(iB)/hTPC3sigmaMC[iS]->GetBinContent(iB);
      ratio_mc_err = ratio_mc * TMath::Sqrt( Sq(hRawCountsMC[iS]->GetBinError(iB)/hRawCountsMC[iS]->GetBinContent(iB)) + Sq(hTPC3sigmaMC[iS]->GetBinError(iB)/hTPC3sigmaMC[iS]->GetBinContent(iB)) );
      fractionMC[iS]->SetBinContent(fractionMC[iS]->FindBin(hRawCountsMC[iS]->GetBinCenter(iB)),ratio_mc);
      fractionMC[iS]->SetBinError(fractionMC[iS]->FindBin(hRawCountsMC[iS]->GetBinCenter(iB)),ratio_mc_err);
      //
      correction[iS].push_back(ratio/ratio_mc);
      printf("matter: %c ratio: %f ratio_mc: %f ratio/ratio_mc: %f\n", kLetter[iS], ratio, ratio_mc, ratio/ratio_mc);
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
    hTPC3sigma[iS]->Write(Form("Tpc_%c",kLetter[iS]));
    hTPC3sigmaMC[iS]->Write(Form("TpcMC_%c",kLetter[iS]));
    fraction[iS]->Write();
    fractionMC[iS]->Write();
    cFraction[iS]->cd();
    fraction[iS]->Draw();
    fractionMC[iS]->Draw("same");
    leg.Draw();
    cFraction[iS]->Write();
  }
}
