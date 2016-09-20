#include "src/Common.h"
#include "src/FitModules.h"

#include <memory>
#include <functional>
using std::function;

#include <TAxis.h>
#include <TDirectory.h>
#include <TEfficiency.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraphAsymmErrors.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TList.h>
#include <TMath.h>

#include <RooArgList.h>
#include <RooDataHist.h>
#include <RooMsgService.h>
#include <RooPlot.h>
#include <RooRealVar.h>

using namespace RooFit;

void Signal() {
  /// Suppressing the output
  RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR) ;

  /// Taking all the histograms from the MC file
  TFile input_file(kDataFilename.data());
  TFile output_file(kSignalOutput.data(),"recreate");

  /// Setting up the fitting environment
  RooRealVar m("dm2","m^{2} - m^{2}_{PDG}",-2.,2.5,"GeV^{2}/c^{4}");
  m.setBins(10000,"cache");
  m.setRange("Full", -2., 2.5);
  FitExpExpTailGaus fExpExpTailGaus(&m);
  //

  for (auto list_key : *input_file.GetListOfKeys()) {
    /// Preliminary operation to read the list and create an output dir
    if (string(list_key->GetName()).find(kFilterListNames.data()) == string::npos) continue;

    TTList* list = (TTList*)input_file.Get(list_key->GetName());
    TDirectory* base_dir = output_file.mkdir(list_key->GetName());
    output_file.cd(list_key->GetName());

    /// Taking all the necessary histogram to perform the analysis
    TH3F *fATOFsignal = (TH3F*)list->Get("fATOFsignal");
    TH3F *fMTOFsignal = (TH3F*)list->Get("fMTOFsignal");
    TH2F *fATPCcounts = (TH2F*)list->Get("fATPCcounts");
    TH2F *fMTPCcounts = (TH2F*)list->Get("fMTPCcounts");

    /// Taking information about centrality bins
    auto n_centralities = fATOFsignal->GetNbinsX();
    auto cent_labels = *(fATOFsignal->GetXaxis()->GetXbins());
    /// Taking information about \f$p_{\mathrm{T}}\f$ bins
    auto n_pt_bins = fATOFsignal->GetNbinsY();
    auto pt_axis = fATOFsignal->GetYaxis();
    auto pt_labels = *(pt_axis->GetXbins());

    /// Now it comes a bit of a complication. To improve fit quality
    /// and to minimize the manual interventions one should do the fits to
    /// the same pT bin for all the species and all the centrality classes
    /// before moving to another pT bin. This requires to build all the
    /// arrays before starting the actual analysis.

    /// Build arrays to analyse both deuteron and anti-deuterons
    /// with the same code
    TH3F* tof_histo[2] = {fMTOFsignal,fATOFsignal};
    TH2F* tpc_histo[2] = {fMTPCcounts,fATPCcounts};

    /// Build arrays to analyse all the centrality classes for
    /// both deuteron and anti-deuterons. Complicate stuff just for fun.
    TH1D* hRawCounts[2][n_centralities];
    TH1D* hSystFit[2][n_centralities];
    TH1D* hSignificance[2][n_centralities];
    TH1D* hTPConly[2][n_centralities];
    //
    for (int iS = 0; iS < 2; ++iS) {
      for (int iC = 0; iC < n_centralities; ++iC) {
        hTPConly[iS][iC] = new TH1D(Form("hTPConly%c%i",kLetter[iS],iC),";p_{T} GeV/c; TPC raw counts",n_pt_bins,pt_labels.GetArray());
        hSignificance[iS][iC] = new TH1D(Form("hSignificance%c%i",kLetter[iS],iC),
                                         Form("%2.0f-%2.0f%%; p_{T}(GeV/c); #frac{S}{#sqrt{S+B}}",
                                              cent_labels[iC],cent_labels[iC + 1]),
                                              n_pt_bins,pt_labels.GetArray());
        hRawCounts[iS][iC] = new TH1D(Form("hRawCounts%c%i",kLetter[iS],iC),
                                      Form("%2.0f-%2.0f%%; p_{T}(GeV/c); RawCounts",
                                            cent_labels[iC],cent_labels[iC + 1]),n_pt_bins,pt_labels.GetArray());
      }
    }

    /// Creating the directories to be used to store the results
    for (int iS = 0; iS < 2; ++iS) {
      TDirectory* dir = base_dir->mkdir(kNames[iS].data());
      dir->cd();
      dir->mkdir("Fits");
      dir->mkdir("Significance");
      dir->mkdir("TPConly");
    }

    for (int iB = 0; iB < n_pt_bins; ++iB) {
      if (pt_axis->GetBinCenter(iB+1) < kPtRange[0] ||
          pt_axis->GetBinCenter(iB+1) > kPtRange[1])
        continue;
      fExpExpTailGaus.SetHighPt(pt_axis->GetBinCenter(iB+1) > 5.);
      for (int iS = 0; iS < 2; ++iS) {
        for (int iC = 0; iC < n_centralities; ++iC) {

          TString iTitle = Form("%1.1f #leq p_{T} < %1.1f",pt_labels[iB],pt_labels[iB + 1]);
          TString iName = Form("d%i_%i",iC,iB);

          TH1D *dat = tof_histo[iS]->ProjectionZ(Form("data%i_%i",iC,iB),iC + 1,iC + 1, iB + 1, iB + 1);
          RooDataHist data("data","data",RooArgList(m),Import(*dat));

          base_dir->cd(Form("%s/Fits",kNames[iS].data()));

          fExpExpTailGaus.FitData(dat, iName, iTitle, "Full");
          fExpExpTailGaus.mPlot->Write();

          const float sVal = fExpExpTailGaus.mSigCounts->getVal();
          const float sErr = fExpExpTailGaus.mSigCounts->getError();

          hRawCounts[iS][iC]->SetBinContent(iB + 1, sVal);
          hRawCounts[iS][iC]->SetBinError(iB + 1, sErr);
          hSignificance[iS][iC]->SetBinContent(iB + 1, sVal / sqrt(sVal + fExpExpTailGaus.mBkgCounts->getVal()));

          /// TPC only raw signal
          hTPConly[iS][iC]->SetBinContent(iB + 1, tpc_histo[iS]->GetBinContent(iC + 1, iB + 1));
        }
      }
    }

    for (int iS = 0; iS < 2; ++iS) {
      for (int iC = 0; iC < n_centralities; ++iC) {
        base_dir->cd(Form("%s/Fits",kNames[iS].data()));
        hRawCounts[iS][iC]->Write();
        base_dir->cd(Form("%s/Significance",kNames[iS].data()));
        hSignificance[iS][iC]->Write();
        base_dir->cd(Form("%s/TPConly",kNames[iS].data()));
        hTPConly[iS][iC]->Write();
      }
    }
    base_dir->Close();
  }
  output_file.Close();
}

