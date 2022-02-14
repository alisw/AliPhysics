#include "src/Common.h"
#include "src/Utils.h"
using namespace utils;
#include "src/Plotting.h"

#include <TCanvas.h>
#include <TLegend.h>
#include <TEfficiency.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraphAsymmErrors.h>
#include <TH2F.h>
#include <TStyle.h>
#include <TError.h>

int kCentBins[10][2] = {{2,2},{3,3},{4,4},{5,5},{6,6},{7,7},{8,8},{9,10},{11,13},{2,13}};
int kNcentBins = 10;
const char* kLegendLabels[2] = {"Deuterons","Antideuterons"};

TH1* DoEff(TH1* hNum, TH1* hDen, string name, char letter, int iC, TArrayD& cent_labels) {
  float lmin=0., lmax=0.;
  if(iC==0){
    lmin = cent_labels[0];
    lmax = cent_labels[cent_labels.GetSize()-1];
  }
  else{
    int iC0 = kCentBins[iC-1][0];
    int iC1 = kCentBins[iC-1][1];
    lmin = cent_labels[iC0-1];
    lmax = cent_labels[iC1];
  }
  TH1* hEff = (TH1*)hNum->Clone(Form("eff%s%c%i",name.data(),letter,iC));
  hEff->SetTitle(Form("%3.0f - %3.0f%%",lmin,lmax));
  hEff->GetYaxis()->SetRangeUser(0.f,1.1f);
  hEff->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  hEff->GetYaxis()->SetTitle("Efficiency x Acceptance");
  hEff->SetMarkerStyle(24);
  hEff->SetMarkerSize(0.7);
  hEff->SetLineColor(kBlack);
  hEff->SetMarkerColor(kBlack);
  for (int iBin = 1; iBin <= hEff->GetNbinsX(); ++iBin) {
    double num = hNum->GetBinContent(iBin);
    double den = hDen->GetBinContent(iBin);
    if (std::abs(den) < 1.e-24) continue;
    double eff = num/den;
    hEff->SetBinContent(iBin, eff);
    hEff->SetBinError(iBin, std::sqrt(eff * (1. - eff) / den));
  }
  hEff->Write();
  return hEff;
}

void Efficiency(bool MBonly = false) {
  /// Taking all the histograms from the MC file
  TFile input_file(kMCfilename.data());
  const char* file_name_MB = (kUseIntegratedForMB) ? kMCfilename.data() : kMCfilenameMB.data();
  TFile input_file_MB(file_name_MB);
  TFile output_file(kEfficiencyOutput.data(),"recreate");

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gErrorIgnoreLevel=kError;

  int iList = 0;
  int counter=0;
  for (auto list_key : *input_file.GetListOfKeys()) {
    /// Preliminary operation to read the list and create an output dir
    // if (string(list_key->GetName()).find(kFilterListNames.data()) == string::npos) continue;
    if (string(list_key->GetName()).find(kFilterListNames.data()) == string::npos) continue;
    if (string(list_key->GetName())==Form("%schisquare0",kFilterListNames.data())) continue;
    if (string(list_key->GetName()).find("_pid2") != string::npos) continue;
    if (string(list_key->GetName()).find("_pid3") != string::npos) continue;
    if (string(list_key->GetName()).find("_sectemplate") != string::npos) continue;
    TTList* list = (TTList*)input_file.Get(list_key->GetName());
    string out_list = list_key->GetName();
    printf("list name: %s\n",out_list.data());
    string list_name_MB = out_list;
    if(!kUseIntegratedForMB) list_name_MB.insert(kFilterListNames.size()-1,"MB");
    printf("list name MB: %s\n",list_name_MB.data());
    TTList* listMB = (TTList*)input_file_MB.Get(list_name_MB.data());
    //replace(out_list,"mpuccio","nuclei");
    TDirectory* base_dir = output_file.mkdir(out_list.data());
    base_dir->cd();
    /// Getting all the histograms
    TH2F  *fITS_TPC[2],*fTotal[2],*fITS_TPC_TOF[2];
    for (int iS = 0; iS < 2; ++iS) {
      fTotal[iS] = dynamic_cast<TH2F*>(list->Get(Form("f%cTotal",kLetter[iS])));
      fITS_TPC[iS] = dynamic_cast<TH2F*>(list->Get(Form("f%cITS_TPC",kLetter[iS])));
      fITS_TPC_TOF[iS] = dynamic_cast<TH2F*>(list->Get(Form("f%cITS_TPC_TOF",kLetter[iS])));
    }

    /// Getting all the histograms for the MB
    TH2F  *fITS_TPC_MB[2],*fTotal_MB[2],*fITS_TPC_TOF_MB[2];
    for (int iS = 0; iS < 2; ++iS) {
      fTotal_MB[iS] = dynamic_cast<TH2F*>(listMB->Get(Form("f%cTotal",kLetter[iS])));
      fITS_TPC_MB[iS] = dynamic_cast<TH2F*>(listMB->Get(Form("f%cITS_TPC",kLetter[iS])));
      fITS_TPC_TOF_MB[iS] = dynamic_cast<TH2F*>(listMB->Get(Form("f%cITS_TPC_TOF",kLetter[iS])));
    }

    /// Taking information about centrality bins
    auto n_centralities = fTotal[0]->GetNbinsX();
    auto cent_labels = *(fTotal[0]->GetXaxis()->GetXbins());

    TCanvas* cEff_MB = new TCanvas("cEff_MB","c_Eff_MB");
    TLegend leg(0.58,0.21,0.93,0.41);
    leg.SetBorderSize(0);

    TCanvas* cEff[kCentLength];
    TLegend* legCent[kCentLength];
    for(int iC=0; iC<kCentLength; iC++){
      cEff[iC] = new TCanvas(Form("cEff_%d",iC),Form("cEff_%d",iC));
      legCent[iC] = new TLegend(0.58,0.21,0.93,0.41);
      legCent[iC]->SetBorderSize(0);
    }

    TH1* effTofMB[2];
    TH1* effTpcMB[2];

    /// Loop to analyse the centrality bins individually
    for (int iS = 0; iS < 2; ++iS) {
      TDirectory* dir = base_dir->mkdir(kNames[iS].data());
      dir->cd();
      TDirectory* MB_dir = dir->mkdir("MB");
      MB_dir->cd();
      TH1D *tpcMB_tmp = (TH1D*)fITS_TPC_MB[iS]->ProjectionY("tpc_MB_tmp");
      TH1D *tpcMB = (TH1D*)tpcMB_tmp->Rebin(kNPtBins,"tpc_MB",kPtBins);
      TH1D *tofMB_tmp = (TH1D*)fITS_TPC_TOF_MB[iS]->ProjectionY("tof_MB_tmp");
      TH1D *tofMB = (TH1D*)tofMB_tmp->Rebin(kNPtBins,"tof_MB",kPtBins);
      TH1D *totMB_tmp = (TH1D*)fTotal_MB[iS]->ProjectionY("tot_MB_tmp");
      TH1D *totMB = (TH1D*)totMB_tmp->Rebin(kNPtBins,"tot_MB",kPtBins);

      effTofMB[iS] = DoEff(tofMB,totMB,"Tof",kLetter[iS],0,cent_labels);
      effTofMB[iS]->SetTitle("");
      effTpcMB[iS] = DoEff(tpcMB,totMB,"Tpc",kLetter[iS],0,cent_labels);
      effTpcMB[iS]->SetTitle("");
      effTpcMB[iS]->SetMarkerStyle(25);
      if(iS==1){
        effTofMB[iS]->SetMarkerColor(kRed);
        effTofMB[iS]->SetLineColor(kRed);
        effTpcMB[iS]->SetMarkerColor(kRed);
        effTpcMB[iS]->SetLineColor(kRed);
      }
      cEff_MB->cd();
      if(iS==0){
        effTofMB[iS]->Draw("PE");
        effTpcMB[iS]->Draw("PESAME");
      }
      else{
        effTofMB[iS]->Draw("PESAME");
        effTpcMB[iS]->Draw("PESAME");
      }
      leg.AddEntry(effTofMB[iS],Form("%s: TPC +  TOF",kLegendLabels[iS]),"PE");
      leg.AddEntry(effTpcMB[iS],Form("%s: TPC",kLegendLabels[iS]),"PE");
      if(iS==1){
        leg.Draw();
        cEff_MB->Write();
      }

      if(MBonly) continue;
      for (int iC = 1; iC <= kNcentBins; ++iC) {

        dir->cd();
        TDirectory* cent_dir = dir->mkdir(Form("C_%d",iC));
        cent_dir->cd();
        TH1D *tpc_tmp = (iC==10) ? (TH1D*)fITS_TPC_MB[iS]->ProjectionY(Form("tpc%i_tmp",iC)) : (TH1D*)fITS_TPC[iS]->ProjectionY(Form("tpc%i_tmp",iC),kCentBins[iC-1][0],kCentBins[iC-1][1]);
        TH1D *tpc = (TH1D*)tpc_tmp->Rebin(kNPtBins,Form("tpc%i",iC),kPtBins);
        TH1D *tof_tmp = (iC==10) ? (TH1D*)fITS_TPC_TOF_MB[iS]->ProjectionY(Form("tof%i_tmp",iC)) : (TH1D*)fITS_TPC_TOF[iS]->ProjectionY(Form("tof%i_tmp",iC),kCentBins[iC-1][0],kCentBins[iC-1][1]);
        TH1D *tof = (TH1D*)tof_tmp->Rebin(kNPtBins,Form("tof%i",iC),kPtBins);
        TH1D *tot_tmp = (iC==10) ?  (TH1D*)fTotal_MB[iS]->ProjectionY(Form("tot%i_tmp",iC)) : (TH1D*)fTotal[iS]->ProjectionY(Form("tot%i_tmp",iC),kCentBins[iC-1][0],kCentBins[iC-1][1]);
        TH1D *tot = (TH1D*)tot_tmp->Rebin(kNPtBins,Form("tot%i",iC),kPtBins);

        int color = (iS==0) ? kBlack : kRed;
        TH1* effTof = DoEff(tof,tot,"Tof",kLetter[iS],iC,cent_labels);
        plotting::SetHistStyle(effTof,color,20);
        TH1* effTpc = DoEff(tpc,tot,"Tpc",kLetter[iS],iC,cent_labels);
        plotting::SetHistStyle(effTpc,color,21);

        tpc->Sumw2();
        tof->Sumw2();
        int iC0 = kCentBins[iC-1][0];
        int iC1 = kCentBins[iC-1][1];
        plotting::SetHistStyle(tpc,plotting::kSpectraColors[iC-1]);
        tpc->SetTitle(Form("%3.0f - %3.0f%%",cent_labels[iC0-1],cent_labels[iC1]));
        tpc->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
        tpc->GetYaxis()->SetTitle("Ratio to MB");
        plotting::SetHistStyle(tof,plotting::kSpectraColors[iC-1]);
        tof->SetTitle(Form("%3.0f - %3.0f%%",cent_labels[iC0-1],cent_labels[iC1]));
        tof->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
        tof->GetYaxis()->SetTitle("Ratio to MB");

        cEff[iC-1]->cd();
        if(iS==0){
          effTof->Draw("PE");
          effTpc->Draw("PESAME");
        }
        else{
          effTof->Draw("PESAME");
          effTpc->Draw("PESAME");
        }
        legCent[iC-1]->AddEntry(effTof,Form("%s: TPC +  TOF",kLegendLabels[iS]),"PE");
        legCent[iC-1]->AddEntry(effTpc,Form("%s: TPC",kLegendLabels[iS]),"PE");
        if(iS==1){
          legCent[iC-1]->Draw();
          cEff[iC-1]->Write();
        }

        tpc->Multiply(totMB);
        tpc->Divide(tot);
        tpc->Divide(tpcMB);
        tpc->Write(Form("Ratio2MBtpc%c%i",kLetter[iS],iC));
        tof->Multiply(totMB);
        tof->Divide(tot);
        tof->Divide(tofMB);
        tof->Write(Form("Ratio2MBtof%c%i",kLetter[iS],iC));
      }
    }
  }
  output_file.Close();
}
