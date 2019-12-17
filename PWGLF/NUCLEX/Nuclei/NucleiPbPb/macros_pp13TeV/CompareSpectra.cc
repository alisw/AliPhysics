#include "src/Common.h"
#include "src/Utils.h"
using namespace utils;
#include "src/Plotting.h"
using namespace plotting;

#include <map>
#include <vector>
#include <array>
using std::array;
#include <memory>
#include <cmath>

#include <Riostream.h>
#include <TFile.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TDirectory.h>
#include <TLine.h>

double zcalc(double m0, double m1, double s0, double s1){
  return (m0-m1)/sqrt(s0*s0+s1*s1);
}

void CompareSpectra(){

  TFile spectra_file_2016((kBaseOutputDir + "spectra.root").data());
  TFile spectra_file_2017((kBaseOutputDir + "spectra_2017.root").data());
  TFile comp_file((kBaseOutputDir + "spectra_comp.root").data(),"recreate");
  TDirectory* tpc_dir = comp_file.mkdir("tpc");
  TDirectory* tof_dir = comp_file.mkdir("tof");
  TDirectory* all_dir = comp_file.mkdir("all");

  const char* kRomanLabels[10] = {"I","II","III","IV + V","VI","VII","VIII","IX","X", "I - X"};
  double y_limits[10] = {1e-3,8e-4,7e-4,6e-4,5e-4,4e-4,3e-4,2e-4,1e-4,3e-4};

  TH1F* stat_tof_2016[2][kCentLength];
  TH1F* stat_tof_2017[2][kCentLength];

  TH1F* stat_tpc_2016[2][kCentLength];
  TH1F* stat_tpc_2017[2][kCentLength];

  TH1F* stat_all_2016[2][kCentLength];
  TH1F* stat_all_2017[2][kCentLength];

  TH1F* ztest[2][kCentLength];

  TCanvas* cv_comp[2][kCentLength][3];
  TCanvas* cZtest[2];

  for(int iS = 0; iS < 2; iS++){
    for(int iC=0; iC<kCentLength; iC++){
      string tof_basepath = kFilterListNames + "/" + kNames[iS] + "/" + std::to_string(iC) + "TOF/TOFspectra" + std::to_string(iC);
      stat_tof_2016[iS][iC] = (TH1F*)spectra_file_2016.Get(tof_basepath.data());
      Requires(stat_tof_2016[iS][iC],("Missing 2016" + tof_basepath).data());
      plotting::SetHistStyle(stat_tof_2016[iS][iC],kRed,20);
      stat_tof_2017[iS][iC] = (TH1F*)spectra_file_2017.Get(tof_basepath.data());
      Requires(stat_tof_2017[iS][iC],("Missing 2017" + tof_basepath).data());
      plotting::SetHistStyle(stat_tof_2017[iS][iC],kBlue,20);
      string tpc_basepath = kFilterListNames + "/" + kNames[iS] + "/" + std::to_string(iC) + "TPC/TPCspectra" + std::to_string(iC);
      stat_tpc_2016[iS][iC] = (TH1F*)spectra_file_2016.Get(tpc_basepath.data());
      Requires(stat_tpc_2016[iS][iC],("Missing 2016" + tpc_basepath).data());
      plotting::SetHistStyle(stat_tpc_2016[iS][iC],kRed,21);
      stat_all_2016[iS][iC] = (TH1F*)stat_tof_2016[iS][iC]->Clone(Form("stat_all_2016_%c_%d",kLetter[iS],iC));
      stat_all_2016[iS][iC]->Reset();
      stat_tpc_2017[iS][iC] = (TH1F*)spectra_file_2017.Get(tpc_basepath.data());
      Requires(stat_tpc_2017[iS][iC],("Missing 2017" + tpc_basepath).data());
      plotting::SetHistStyle(stat_tpc_2017[iS][iC],kBlue,21);
      stat_all_2017[iS][iC] = (TH1F*)stat_tof_2017[iS][iC]->Clone(Form("stat_all_2017_%c_%d",kLetter[iS],iC));
      stat_all_2017[iS][iC]->Reset();
      ztest[iS][iC] = (TH1F*)stat_tof_2017[iS][iC]->Clone(Form("ztest_%c_%d",kLetter[iS],iC));
      plotting::SetHistStyle(ztest[iS][iC],kSpectraColors[iC],20);
      for(int iB=1; iB<=kNPtBins; iB++){
        if(stat_all_2016[iS][iC]->GetBinCenter(iB)<0.6 || stat_all_2016[iS][iC]->GetBinCenter(iB)> kCentPtLimits[iC]) continue;
        if(stat_all_2016[iS][iC]->GetBinCenter(iB)<1){
          stat_all_2016[iS][iC]->SetBinContent(iB,stat_tpc_2016[iS][iC]->GetBinContent(iB));
          stat_all_2016[iS][iC]->SetBinError(iB,stat_tpc_2016[iS][iC]->GetBinError(iB));
          stat_all_2017[iS][iC]->SetBinContent(iB,stat_tpc_2017[iS][iC]->GetBinContent(iB));
          stat_all_2017[iS][iC]->SetBinError(iB,stat_tpc_2017[iS][iC]->GetBinError(iB));
        }
        else{
          stat_all_2016[iS][iC]->SetBinContent(iB,stat_tof_2016[iS][iC]->GetBinContent(iB));
          stat_all_2016[iS][iC]->SetBinError(iB,stat_tof_2016[iS][iC]->GetBinError(iB));
          stat_all_2017[iS][iC]->SetBinContent(iB,stat_tof_2017[iS][iC]->GetBinContent(iB));
          stat_all_2017[iS][iC]->SetBinError(iB,stat_tof_2017[iS][iC]->GetBinError(iB));
        }
        ztest[iS][iC]->SetBinContent(iB, zcalc(stat_all_2016[iS][iC]->GetBinContent(iB),stat_all_2017[iS][iC]->GetBinContent(iB),stat_all_2016[iS][iC]->GetBinError(iB),stat_all_2017[iS][iC]->GetBinError(iB)));
      }

      cv_comp[iS][iC][0] = new TCanvas(Form("cv_%c_%d_TOF",kLetter[iS],iC),Form("cv_%c_%d_TOF",kLetter[iS],iC),800,600);
      cv_comp[iS][iC][0]->SetLeftMargin(0.20);
      TH1* hFrame_0 = cv_comp[iS][iC][0]->DrawFrame(
          0.4,
          1e-6,
          3.8,
          y_limits[iC],
          Form("Mult. %s ;#it{p}_{T} (GeV/#it{c});#frac{1}{#it{N}_{ev}} #frac{d^{2}#it{N}}{d#it{p}_{T}d#it{y}} (GeV/#it{c})^{-1}",kRomanLabels[iC])
          );
      hFrame_0->GetYaxis()->SetTitleOffset(1.7);
      stat_tof_2016[iS][iC]->Draw("same");
      stat_tof_2017[iS][iC]->Draw("same");
      TLegend* leg_tof = new TLegend(0.780702,0.650435,0.87218,0.86087);
      leg_tof->SetBorderSize(0);
      leg_tof->SetTextSize(0.028);
      leg_tof->AddEntry(stat_tof_2016[iS][iC],"2016","pe");
      leg_tof->AddEntry(stat_tof_2017[iS][iC],"2017","pe");
      leg_tof->Draw();
      cv_comp[iS][iC][1] = new TCanvas(Form("cv_%c_%d_TPC",kLetter[iS],iC),Form("cv_%c_%d_TPC",kLetter[iS],iC),800,600);
      cv_comp[iS][iC][1]->SetLeftMargin(0.20);
      TH1* hFrame_1 =cv_comp[iS][iC][1]->DrawFrame(
          0.4,
          1e-6,
          3.8,
          y_limits[iC],
          Form("Mult. %s;#it{p}_{T} (GeV/#it{c});#frac{1}{#it{N}_{ev}} #frac{d^{2}#it{N}}{d#it{p}_{T}d#it{y}} (GeV/#it{c})^{-1}",kRomanLabels[iC])
          );
      hFrame_1->GetYaxis()->SetTitleOffset(1.7);
      stat_tpc_2016[iS][iC]->Draw("same");
      stat_tpc_2017[iS][iC]->Draw("same");
      TLegend* leg_tpc = new TLegend(0.780702,0.650435,0.87218,0.86087);
      leg_tpc->SetBorderSize(0);
      leg_tpc->SetTextSize(0.028);
      leg_tpc->AddEntry(stat_tpc_2016[iS][iC],"2016","pe");
      leg_tpc->AddEntry(stat_tpc_2017[iS][iC],"2017","pe");
      leg_tpc->Draw();

      cv_comp[iS][iC][2] = new TCanvas(Form("cv_%c_%d_all",kLetter[iS],iC),Form("cv_%c_%d_all",kLetter[iS],iC),800,600);
      cv_comp[iS][iC][2]->SetLeftMargin(0.20);
      TH1* hFrame_2 =cv_comp[iS][iC][2]->DrawFrame(
          0.4,
          1e-6,
          3.8,
          y_limits[iC],
          Form("Mult. %s;#it{p}_{T} (GeV/#it{c});#frac{1}{#it{N}_{ev}} #frac{d^{2}#it{N}}{d#it{p}_{T}d#it{y}} (GeV/#it{c})^{-1}",kRomanLabels[iC])
          );
      hFrame_2->GetYaxis()->SetTitleOffset(1.7);
      stat_all_2016[iS][iC]->Draw("same");
      stat_all_2017[iS][iC]->Draw("same");
      TLegend* leg_all = new TLegend(0.780702,0.650435,0.87218,0.86087);
      leg_all->SetBorderSize(0);
      leg_all->SetTextSize(0.028);
      leg_all->AddEntry(stat_all_2016[iS][iC],"2016","pe");
      leg_all->AddEntry(stat_all_2017[iS][iC],"2017","pe");
      leg_all->Draw();

      comp_file.cd();
      tof_dir->cd();
      cv_comp[iS][iC][0]->Write();
      stat_tof_2016[iS][iC]->Write();
      stat_tof_2017[iS][iC]->Write();
      tpc_dir->cd();
      cv_comp[iS][iC][1]->Write();
      stat_tpc_2016[iS][iC]->Write();
      stat_tpc_2017[iS][iC]->Write();
      all_dir->cd();
      stat_all_2016[iS][iC]->Write();
      stat_all_2017[iS][iC]->Write();
      cv_comp[iS][iC][2]->Write();
      cv_comp[iS][iC][2]->SaveAs(Form("%scv_%c_%d_all.png",kFiguresFolder.data(),kLetter[iS],iC));
    }

    cZtest[iS] = new TCanvas(Form("cZtest_%c",kLetter[iS]),Form("cZtest_%c",kLetter[iS]),800,600);
    cZtest[iS]->DrawFrame(0.4, -8., 3.8, 8,";#it{p}_{T} (GeV/#it{c});#frac{N_{2016} - N_{2017}}{#sigma}");
    TLegend* z_leg = new TLegend(0.827538,0.596386,0.901947,0.875904);
    z_leg->SetBorderSize(0);
    for (int iC = 0; iC < kCentLength; ++iC) {
      ztest[iS][iC]->Draw("psame");
      z_leg->AddEntry(ztest[iS][iC],kRomanLabels[iC],"p");
    }
    for(int i=-3; i<=3; i++){
      TLine *line = new TLine(0.4,i,3.8,i);
      int absi = (i<0) ? -1*i : i;
      switch(absi){
        case(3):
          line->SetLineColor(kRed);
          break;
        case(2):
          line->SetLineColor(kOrange);
          break;
        case(1):
          line->SetLineColor(kGreen+3);
          break;
        default:
          line->SetLineColor(kBlack);
          break;
      }
      line->SetLineStyle(2);
      line->Draw();
    }
    z_leg->Draw();
    all_dir->cd();
    cZtest[iS]->Write();
    cZtest[iS]->SaveAs(Form("%szTest_%c.png",kFiguresFolder.data(),kLetter[iS]));
  }



}
