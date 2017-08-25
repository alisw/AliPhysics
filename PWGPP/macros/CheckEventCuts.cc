#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TLegend.h>
#include <TList.h>
#include <TString.h>
#include <TStyle.h>

#include <array>
using std::array;

/// Now it runs only with ROOT6

namespace {
  array<TString,2> raw_sel{"_raw","_selected"};
  array<TString,3> basic_hists_names{"Vtz","DeltaVtz","Centrality"};
  array<TString,2> multiplicity_hists_names{"EstimCorrelation","MultCentCorrelation"};
  array<TString,5> advanced_hists_names{"fTOFvsFB32","fTPCvsAll","fMultvsV0M","fTPCvsTrkl","fVZEROvsTPCout"};
  array<TString,2> compare_titles{"Old","New"};
  array<Color_t,2> compare_colors{kRed,kBlack};
  array<Style_t,2> compare_fills{3345,3354};
};

template<typename T, int N> array<T*,N> getHistograms(TList* list, array<TString,N> names, TString &suffix) {
  array<T*,N> histos{nullptr};
  for (size_t iH = 0; iH < N; ++iH) {
    histos[iH] = static_cast<T*>(list->FindObject((names[iH]+suffix).Data()));
  }
  return histos;
}

void CheckEventCuts(TString input_file_name = "AnalysisResults.root") {
  gStyle->SetOptStat(0);

  TFile input_file(input_file_name.Data());
  TList* input_list = static_cast<TList*>(input_file.Get("_std"));

  TCanvas* cv_basics[2]{nullptr};
  TCanvas* cv_multiplicity[2]{nullptr};
  TCanvas* cv_advanced[2]{nullptr};

  for (int iC = 0; iC < 2; ++iC) {
    cv_basics[iC] = new TCanvas(Form("cv_basics%s",raw_sel[iC].Data()));
    cv_basics[iC]->Divide(2,2);
    auto hist_basics = getHistograms<TH1,3>(input_list,basic_hists_names,raw_sel[iC]);
    for (size_t iN = 0; iN < basic_hists_names.size(); ++iN) {
      cv_basics[iC]->cd(iN+1);
      hist_basics[iN]->DrawCopy();
    }

    cv_multiplicity[iC] = new TCanvas(Form("cv_multiplicity%s",raw_sel[iC].Data()));
    cv_multiplicity[iC]->Divide(2);
    auto hist_multiplicity = getHistograms<TH2,2>(input_list,multiplicity_hists_names,raw_sel[iC]);
    for (size_t iN = 0; iN < hist_multiplicity.size(); ++iN) {
      cv_multiplicity[iC]->cd(iN+1);
      hist_multiplicity[iN]->DrawCopy("colz");
    }

    cv_advanced[iC] = new TCanvas(Form("cv_advanced%s",raw_sel[iC].Data()));
    cv_advanced[iC]->Divide(3,2);
    auto hist_advanced = getHistograms<TH2,5>(input_list,advanced_hists_names,raw_sel[iC]);
    for (size_t iN = 0; iN < hist_advanced.size(); ++iN) {
      cv_advanced[iC]->cd(iN+1);
      hist_advanced[iC]->DrawCopy("colz");
    }
  }

  TCanvas* cv_summary{new TCanvas("cv_summary")};
  cv_summary->Divide(2);
  cv_summary->cd(1);
  TH1* hist_cut_stats = static_cast<TH1*>(input_list->FindObject("fCutStats"));
  hist_cut_stats->DrawCopy();
  cv_summary->cd(2);
  TH1* hist_normalisation = static_cast<TH1*>(input_list->FindObject("fNormalisationHist"));
  hist_normalisation->DrawCopy();
}

void CompareEventSelections(TString input_file_name0 = "oldAOD.root",
                            TString input_file_name1 = "newAOD.root") {
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  array<TFile*,2> input_files{TFile::Open(input_file_name0.Data()),TFile::Open(input_file_name1.Data())};
  array<TList*,2> input_lists{nullptr};

  array<TH1*,2> hists_cut_stats{nullptr};
  array<TH1*,2> hists_normalisation{nullptr};
  array<double,2> number_of_events{-1.};

  array<array<TH1*,3>,2> basics_hist{{nullptr}};
  array<array<TH2*,2>,2> multiplicity_hist{{nullptr}};

  TCanvas* cv_summary{new TCanvas("cv_summary","Summary",1000,700)};
  cv_summary->Divide(2);

  TCanvas* cv_basics = new TCanvas("cv_basics","Basic selections",1000,700);
  cv_basics->Divide(2,2);

  TCanvas* cv_multiplicity = new TCanvas("cv_multiplicity","Multiplicity selections",1000,700);
  cv_multiplicity->Divide(2);

  for (int iF = 0; iF < 2; ++iF) {
    input_lists[iF] = static_cast<TList*>(input_files[iF]->Get("_std"));

    hists_cut_stats[iF] = static_cast<TH1*>(input_lists[iF]->FindObject("fCutStats"));
    hists_normalisation[iF] = static_cast<TH1*>(input_lists[iF]->FindObject("fNormalisationHist"));
    number_of_events[iF] = hists_cut_stats[iF]->GetBinContent(4);

    basics_hist[iF] = getHistograms<TH1,3>(input_lists[iF],basic_hists_names,raw_sel[1]);
    multiplicity_hist[iF] = getHistograms<TH2,2>(input_lists[iF],multiplicity_hists_names,raw_sel[1]);
    const double norm = 1000. / number_of_events[iF];

    cv_summary->cd(1);
    hists_cut_stats[iF]->SetTitle(compare_titles[iF].Data());
    hists_cut_stats[iF]->UseCurrentStyle();
    hists_cut_stats[iF]->Scale(norm);
    hists_cut_stats[iF]->SetLineColor(compare_colors[iF]);
    hists_cut_stats[iF]->SetFillStyle(compare_fills[iF]);
    hists_cut_stats[iF]->SetFillColor(compare_colors[iF]);
    hists_cut_stats[iF]->DrawCopy(iF ? "hist same" : "hist");

    cv_summary->cd(2);
    hists_normalisation[iF]->SetTitle(compare_titles[iF].Data());
    hists_normalisation[iF]->UseCurrentStyle();
    hists_normalisation[iF]->Scale(norm);
    hists_normalisation[iF]->SetLineColor(compare_colors[iF]);
    hists_normalisation[iF]->SetFillStyle(compare_fills[iF]);
    hists_normalisation[iF]->SetFillColor(compare_colors[iF]);
    hists_normalisation[iF]->DrawCopy(iF ? "hist same" : "hist");

    for (size_t iH = 0; iH < basics_hist[iF].size(); ++iH) {
      cv_basics->cd(iH+1);
      basics_hist[iF][iH]->SetTitle(compare_titles[iF].Data());
      basics_hist[iF][iH]->UseCurrentStyle();
      basics_hist[iF][iH]->Scale(norm);
      basics_hist[iF][iH]->SetLineColor(compare_colors[iF]);
      basics_hist[iF][iH]->SetFillStyle(compare_fills[iF]);
      basics_hist[iF][iH]->SetFillColor(compare_colors[iF]);
      basics_hist[iF][iH]->DrawCopy(iF ? "same hist" : "hist");
    }

    for (size_t iH = 0; iH < multiplicity_hist[iF].size(); ++iH) {
      cv_multiplicity->cd(iH+1)->SetLogz();
      if (iH) cv_multiplicity->cd(iH+1)->SetLogy();
      multiplicity_hist[iF][iH]->SetTitle(compare_titles[iF].Data());
      multiplicity_hist[iF][iH]->UseCurrentStyle();
      multiplicity_hist[iF][iH]->Scale(norm);
      multiplicity_hist[iF][iH]->SetMarkerColor(compare_colors[iF]);
      multiplicity_hist[iF][iH]->SetLineColor(compare_colors[iF]);
      multiplicity_hist[iF][iH]->SetFillColor(compare_colors[iF]);
      multiplicity_hist[iF][iH]->SetFillStyle(compare_fills[iF]);
      multiplicity_hist[iF][iH]->DrawCopy(iF ? "same" : "");
    }
  }

  cv_summary->cd(2)->BuildLegend();
  cv_basics->cd(2)->BuildLegend(0.65);
  cv_multiplicity->cd(2)->BuildLegend(0.16,0.16,0.4,0.3);

  for (auto &file : input_files)
    file->Close();
}