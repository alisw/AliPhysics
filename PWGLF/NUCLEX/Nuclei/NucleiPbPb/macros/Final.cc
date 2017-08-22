#include "src/Common.h"
#include "src/Plotting.h"

#include <map>
#include <vector>
#include <array>
using std::array;

#include <Riostream.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TList.h>
#include <TLegend.h>
#include <TLine.h>


const string labels[3] = {"0-10%","10-20%","20-40%"};
void MyFinal() {
  TFile spectra_file(kSpectraOutput.data());
  TFile TOF_systematics_file(kSystematicsOutput.data());
  TFile TPC_systematics_file(kSystematicsOutputTPC.data());
  TFile final_file(kFinalOutput.data(),"recreate");

  const double pt_bin_limits[16] = {0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.4,1.6,1.8,2.0,2.2,2.6,3.0,3.4,3.8};//,4.4};
  const int n_pt_bins = 15;

  const int n_centralities = 1;

  TH1F* stat[2][n_centralities];
  TH1F* syst[2][n_centralities];
  TH1F* stat_tpc[2][n_centralities];
  TH1F* syst_tpc[2][n_centralities];
  for (int iS = 0; iS < 2; ++iS) {
    TDirectory* s_dir = final_file.mkdir(kNames[iS].data());
    for (int iC = 0; iC < n_centralities; ++iC) {
      TDirectory *c_dir = s_dir->mkdir(to_string(iC).data());
      c_dir->cd();
      TH1F* totsyst_tmp = (TH1F*)TOF_systematics_file.Get((to_string(iC) + "/" + kNames[iS] + "/totsyst").data());
      Requires(totsyst_tmp, "Missing totsyt");
      TH1F* totsyst = (TH1F*)totsyst_tmp->Rebin(n_pt_bins,Form("totsyst_%d_%d",iS,iC),pt_bin_limits);
      TH1F* totsyst_tpc_tmp = (TH1F*)TPC_systematics_file.Get((to_string(iC) + "/" + kNames[iS] + "/totsyst_tpc").data());
      Requires(totsyst_tpc_tmp, "Missing totsyt_tpc");
      TH1F* totsyst_tpc = (TH1F*)totsyst_tpc_tmp->Rebin(n_pt_bins,Form("totsyst_tpc_%d_%d",iS,iC),pt_bin_limits);
      string tof_basepath = kFilterListNames + "/" + kNames[iS] + "/TOFspectra" + to_string(iC);
      TH1F* spectra_tof_tmp  = (TH1F*)spectra_file.Get(tof_basepath.data());
      Requires(spectra_tof_tmp,tof_basepath.data());
      stat[iS][iC] =(TH1F*)spectra_tof_tmp->Rebin(n_pt_bins,Form("stat_%d_%d",iS,iC),pt_bin_limits);
      auto ptAxis = stat[iS][iC]->GetXaxis();
      string tpc_basepath = kFilterListNames + "/" + kNames[iS] + "/TPCspectra" + to_string(iC);
      TH1F* spectra_tpc_tmp  = (TH1F*)spectra_file.Get(tpc_basepath.data());
      Requires(spectra_tpc_tmp,tpc_basepath.data());
      stat_tpc[iS][iC] = (TH1F*)spectra_tpc_tmp->Rebin(n_pt_bins,Form("stat_tpc_%d_%d",iS,iC),pt_bin_limits);
      syst[iS][iC]  = (TH1F*)totsyst->Clone(("syst" + to_string(iC)).data());
      syst[iS][iC]->Reset();
      syst_tpc[iS][iC]  = (TH1F*)totsyst->Clone(("syst_tpc" + to_string(iC)).data());
      syst_tpc[iS][iC]->Reset();

      for (int iB = 1; iB <= n_pt_bins; ++iB) {
        if(ptAxis->GetBinCenter(iB)<1.){
          stat[iS][iC]->SetBinContent(iB,0.);
          stat[iS][iC]->SetBinError(iB,0.);
          syst[iS][iC]->SetBinContent(iB,0.);
          syst[iS][iC]->SetBinError(iB,0.);
        }
        else{
          syst[iS][iC]->SetBinContent(iB,stat[iS][iC]->GetBinContent(iB));
          syst[iS][iC]->SetBinError(iB,totsyst->GetBinContent(iB) * stat[iS][iC]->GetBinContent(iB));
        }
        if(ptAxis->GetBinCenter(iB)<1.4){
          syst_tpc[iS][iC]->SetBinContent(iB,stat_tpc[iS][iC]->GetBinContent(iB));
          syst_tpc[iS][iC]->SetBinError(iB,totsyst_tpc->GetBinContent(iB) * stat_tpc[iS][iC]->GetBinContent(iB));
        }
        else{
          stat_tpc[iS][iC]->SetBinContent(iB,0.);
          stat_tpc[iS][iC]->SetBinError(iB,0.);
          syst_tpc[iS][iC]->SetBinContent(iB,0.);
          syst_tpc[iS][iC]->SetBinError(iB,0.);
        }
      }
      SetHistStyle(stat[iS][iC],0);
      SetHistStyle(syst[iS][iC],0);
      stat[iS][iC]->Write("stat_tof");
      syst[iS][iC]->Write("syst_tof");
      SetHistStyle(stat_tpc[iS][iC],4);
      SetHistStyle(syst_tpc[iS][iC],4);
      stat_tpc[iS][iC]->Write("stat_tpc");
      syst_tpc[iS][iC]->Write("syst_tpc");

    }
    s_dir->cd();
    TCanvas spectra("spectra","spectra");
    spectra.DrawFrame(
        0.5 * kPtRange[0],
        0.5 * syst[iS][n_centralities-1]->GetBinContent(syst[iS][n_centralities-1]->GetXaxis()->GetNbins()),
        1.1 * 4.4,
        1.8 * syst_tpc[iS][0]->GetMaximum(),
        ";#it{p}_{T} (GeV/#it{c});#frac{1}{#it{N}_{ev}} #frac{d^{2}#it{N}}{d#it{p}_{T}d#it{y}}"
        );
    TLegend final_leg(0.61,0.60,0.91,0.83);
    final_leg.SetBorderSize(0);
    final_leg.SetHeader(Form("%s, pp #sqrt{s} = 13 TeV",kNames[iS].data()));
    for (int iC = 0; iC < n_centralities; ++iC) {
      stat[iS][iC]->SetMinimum(1e-6);
      stat[iS][iC]->Draw("esamex0");
      syst[iS][iC]->SetMinimum(1e-6);
      syst[iS][iC]->Draw("e2same");
      stat_tpc[iS][iC]->SetMinimum(1e-6);
      stat_tpc[iS][iC]->Draw("esamex0");
      syst_tpc[iS][iC]->SetMinimum(1e-6);
      syst_tpc[iS][iC]->Draw("e2same");
      final_leg.AddEntry(syst[iS][iC],"TOF spectrum","fp");
      final_leg.AddEntry(syst_tpc[iS][iC],"TPC spectrum","fp");
    }
    final_leg.Draw();
    spectra.SetLogy();
    spectra.Write();
    if (kPrintFigures) {
      spectra.SaveAs((kFiguresFolder + "spectraTOF" + kLetter[iS] + ".eps").data());
      spectra.SaveAs((kMacrosFolder + "spectraTOF" + kLetter[iS] + ".C").data());
    }
  }

  TDirectory* r_dir = final_file.mkdir("ratio");
  for (int iC = 0; iC < n_centralities; ++iC) {
    r_dir->mkdir(to_string(iC).data())->cd();
    stat[1][iC]->Divide(stat[0][iC]);
    syst[1][iC]->Divide(syst[0][iC]);
    for(int iB=1; iB<=n_pt_bins; iB++){
      if(stat[1][iC]->GetBinCenter(iB)<1.){
        stat[1][iC]->SetBinContent(iB, 0.);
        stat[1][iC]->SetBinError(iB, 0.);
        syst[1][iC]->SetBinContent(iB, 0.);
        syst[1][iC]->SetBinError(iB, 0.);
      }
    }
    stat[1][iC]->Write("stat_tof");
    syst[1][iC]->Write("syst_tof");

    stat_tpc[1][iC]->Divide(stat_tpc[0][iC]);
    syst_tpc[1][iC]->Divide(syst_tpc[0][iC]);
    for(int iB=1; iB<=n_pt_bins; iB++){
      if(stat_tpc[1][iC]->GetBinCenter(iB)>1.4){
        stat_tpc[1][iC]->SetBinContent(iB, 0.);
        stat_tpc[1][iC]->SetBinError(iB, 0.);
        syst_tpc[1][iC]->SetBinContent(iB, 0.);
        syst_tpc[1][iC]->SetBinError(iB, 0.);
      }
    }
    TCanvas ratio("ratio","ratio");
    ratio.DrawFrame(
        0.5 * kPtRange[0],
        0.1,
        1.05 * kPtRange[1],
        1.9,
        ";#it{p}_{T} (GeV/#it{c});#bar{d}/d"
        );
    stat[1][iC]->Draw("esamex0");
    syst[1][iC]->Draw("e2same");
    stat_tpc[1][iC]->Draw("esamex0");
    syst_tpc[1][iC]->Draw("e2same");
    TLine *line = new TLine(0.5 * kPtRange[0],1,1.05 * kPtRange[1],1);
    line->SetLineColor(kBlack);
    line->Draw();
    if (kPrintFigures) ratio.SaveAs((kFiguresFolder + "ratio.eps").data());
    ratio.Write();
  }
}
