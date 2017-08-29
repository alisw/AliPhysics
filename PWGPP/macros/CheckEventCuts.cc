#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TLegend.h>
#include <TList.h>
#include <TString.h>
#include <TStyle.h>

#include <string>
#include <vector>
using std::string;
using std::vector;
#endif

namespace {
  const string raw_sel[2] = {"_raw","_selected"};
  enum {
    n_basic_hists = 3,
    n_multiplicity_hists = 2,
    n_advanced_hists = 5
  };
  const string basic_hists_names[n_basic_hists] = {"Vtz","DeltaVtz","Centrality"};
  const string multiplicity_hists_names[n_multiplicity_hists] = {"EstimCorrelation","MultCentCorrelation"};
  const string advanced_hists_names[n_advanced_hists] = {"fTOFvsFB32","fTPCvsAll","fMultvsV0M","fTPCvsTrkl","fVZEROvsTPCout"};
  const string compare_titles[2] = {"Old","New"};
  const Color_t compare_colors[2] = {kRed,kBlack};
  const Style_t compare_fills[2] = {3345,3354};
}

template<typename T> vector<T* > getHistograms(TList* list, const string* names, const int n_names, const string &suffix) {
  vector<T*> histos(n_names,nullptr);
  for (int iH = 0; iH < n_names; ++iH) {
    histos[iH] = static_cast<T*>(list->FindObject((names[iH]+suffix).data()));
  }
  return histos;
}

void CheckEventCuts(string input_file_name = "AnalysisResults.root") {
  gStyle->SetOptStat(0);

  TFile input_file(input_file_name.data());
  TList* input_list = static_cast<TList*>(input_file.Get("_std"));

  TCanvas* cv_basics[2]{nullptr};
  TCanvas* cv_multiplicity[2]{nullptr};
  TCanvas* cv_advanced[2]{nullptr};

  for (int iC = 0; iC < 2; ++iC) {
    cv_basics[iC] = new TCanvas(Form("cv_basics%s",raw_sel[iC].data()));
    cv_basics[iC]->Divide(2,2);
    auto hist_basics = getHistograms<TH1>(input_list,basic_hists_names,n_basic_hists,raw_sel[iC]);
    for (size_t iN = 0; iN < n_basic_hists; ++iN) {
      cv_basics[iC]->cd(iN+1);
      hist_basics[iN]->DrawCopy();
    }

    cv_multiplicity[iC] = new TCanvas(Form("cv_multiplicity%s",raw_sel[iC].data()));
    cv_multiplicity[iC]->Divide(2);
    auto hist_multiplicity = getHistograms<TH2>(input_list,multiplicity_hists_names,n_multiplicity_hists,raw_sel[iC]);
    for (size_t iN = 0; iN < n_multiplicity_hists; ++iN) {
      cv_multiplicity[iC]->cd(iN+1);
      hist_multiplicity[iN]->DrawCopy("colz");
    }

    cv_advanced[iC] = new TCanvas(Form("cv_advanced%s",raw_sel[iC].data()));
    cv_advanced[iC]->Divide(3,2);
    auto hist_advanced = getHistograms<TH2>(input_list,advanced_hists_names,n_advanced_hists,raw_sel[iC]);
    for (size_t iN = 0; iN < n_advanced_hists; ++iN) {
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

void CompareEventSelections(string input_file_name0 = "oldAOD.root",
                            string input_file_name1 = "newAOD.root") {
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  TFile* input_files[2]{TFile::Open(input_file_name0.data()),TFile::Open(input_file_name1.data())};
  TList* input_lists[2]{nullptr};

  TH1* hists_cut_stats[2]{nullptr};
  TH1* hists_normalisation[2]{nullptr};
  double number_of_events[2]{-1.};

  vector<TH1*> basics_hist[2] = {vector<TH1*>(3,nullptr)};
  vector<TH2*> multiplicity_hist[2] = {vector<TH2*>(2,nullptr)};

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

    basics_hist[iF] = getHistograms<TH1>(input_lists[iF],basic_hists_names,n_basic_hists,raw_sel[1]);
    multiplicity_hist[iF] = getHistograms<TH2>(input_lists[iF],multiplicity_hists_names,n_multiplicity_hists,raw_sel[1]);
    const double norm = 1000. / number_of_events[iF];

    cv_summary->cd(1);
    hists_cut_stats[iF]->SetTitle(compare_titles[iF].data());
    hists_cut_stats[iF]->UseCurrentStyle();
    hists_cut_stats[iF]->Scale(norm);
    hists_cut_stats[iF]->SetLineColor(compare_colors[iF]);
    hists_cut_stats[iF]->SetFillStyle(compare_fills[iF]);
    hists_cut_stats[iF]->SetFillColor(compare_colors[iF]);
    hists_cut_stats[iF]->DrawCopy(iF ? "hist same" : "hist");

    cv_summary->cd(2);
    hists_normalisation[iF]->SetTitle(compare_titles[iF].data());
    hists_normalisation[iF]->UseCurrentStyle();
    hists_normalisation[iF]->Scale(norm);
    hists_normalisation[iF]->SetLineColor(compare_colors[iF]);
    hists_normalisation[iF]->SetFillStyle(compare_fills[iF]);
    hists_normalisation[iF]->SetFillColor(compare_colors[iF]);
    hists_normalisation[iF]->DrawCopy(iF ? "hist same" : "hist");

    for (size_t iH = 0; iH < n_basic_hists; ++iH) {
      cv_basics->cd(iH+1);
      basics_hist[iF][iH]->SetTitle(compare_titles[iF].data());
      basics_hist[iF][iH]->UseCurrentStyle();
      basics_hist[iF][iH]->Scale(norm);
      basics_hist[iF][iH]->SetLineColor(compare_colors[iF]);
      basics_hist[iF][iH]->SetFillStyle(compare_fills[iF]);
      basics_hist[iF][iH]->SetFillColor(compare_colors[iF]);
      basics_hist[iF][iH]->DrawCopy(iF ? "same hist" : "hist");
    }

    for (size_t iH = 0; iH < n_multiplicity_hists; ++iH) {
      cv_multiplicity->cd(iH+1)->SetLogz();
      if (iH) cv_multiplicity->cd(iH+1)->SetLogy();
      multiplicity_hist[iF][iH]->SetTitle(compare_titles[iF].data());
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