#include "src/Common.h"
#include "src/FitModules.h"

#include <memory>

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
  m.setRange("Range1",-1.75, 2.0);
  m.setRange("Range2",-1.65, 1.65);
  m.setRange("Full", -2., 2.5);
  //
  FitExpExpCB fExpExpCB(&m);
  FitExpExpGaus fExpExpGaus(&m);
  FitExpExpTailGaus fExpExpTailGaus(&m);

  for (auto&& list_key : *input_file.GetListOfKeys()) {
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
    TH1D* hFit0[2][n_centralities];
    TH1D* hFit1[2][n_centralities];
    TH1D* hFit2[2][n_centralities];
    TH1D* hFit3[2][n_centralities];
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
        hSystFit[iS][iC] = new TH1D(Form("hSystFit%c%i",kLetter[iS],iC),
                                    Form("%2.0f-%2.0f%%; p_{T}(GeV/c); Systematics from fit",
                                          cent_labels[iC],cent_labels[iC + 1]),n_pt_bins,pt_labels.GetArray());
        hFit0[iS][iC] = new TH1D(Form("hFit1%c%i",kLetter[iS],iC),";p_{T} GeV/c;Bias",n_pt_bins,pt_labels.GetArray());
        hFit1[iS][iC] = new TH1D(Form("hFit2%c%i",kLetter[iS],iC),";p_{T} GeV/c;Bias",n_pt_bins,pt_labels.GetArray());
        hFit2[iS][iC] = new TH1D(Form("hFit3%c%i",kLetter[iS],iC),";p_{T} GeV/c;Bias",n_pt_bins,pt_labels.GetArray());
        hFit3[iS][iC] = new TH1D(Form("hFit4%c%i",kLetter[iS],iC),";p_{T} GeV/c;Bias",n_pt_bins,pt_labels.GetArray());
      }
    }

    /// Creating the directories to be used to store the results
    for (int iS = 0; iS < 2; ++iS) {
      TDirectory* dir = base_dir->mkdir(kNames[iS].data());
      dir->cd();
      dir->mkdir("Model0");
      dir->mkdir("Model1");
      dir->mkdir("Model2");
      dir->mkdir("Range1");
      dir->mkdir("Range2");
      dir->mkdir("Significance");
      dir->mkdir("Systematics");
      dir->mkdir("TPConly");
    }

    for (int iB = 0; iB < n_pt_bins; ++iB) {
      if (pt_axis->GetBinCenter(iB+1) < kPtRange[0] ||
          pt_axis->GetBinCenter(iB+1) > kPtRange[1])
        continue;
      for (int iS = 0; iS < 2; ++iS) {
        for (int iC = 0; iC < n_centralities; ++iC) {

          TString iTitle = Form("%1.1f #leq p_{T} < %1.1f",pt_labels[iB],pt_labels[iB + 1]);
          TString iName = Form("d%i_%i",iC,iB);

          TH1D *dat = tof_histo[iS]->ProjectionZ(Form("dataM%i_%i",iC,iB),iC + 1,iC + 1, iB + 1, iB + 1);
          RooDataHist data("data","data",RooArgList(m),Import(*dat));

          base_dir->cd(Form("%s/Model0",kNames[iS].data()));
          FitModule &mainModule = pt_labels[iB + 1] > 1.6 ? static_cast<FitModule&>(fExpExpGaus) : static_cast<FitModule&>(fExpExpTailGaus);
          mainModule.FitData(dat, iName, iTitle, "Full");
          mainModule.mPlot->Write();

          const float sVal = mainModule.mSigCounts->getVal();
          const float sErr = mainModule.mSigCounts->getError();

          hRawCounts[iS][iC]->SetBinContent(iB + 1, sVal);
          hRawCounts[iS][iC]->SetBinError(iB + 1, sErr);
          hSignificance[iS][iC]->SetBinContent(iB + 1, sVal / sqrt(sVal + mainModule.mBkgCounts->getVal()));

          float refv = sVal;
          float refe = sErr;
          float values[5];
          values[4] = refv;

          base_dir->cd(Form("%s/Model1",kNames[iS].data()));
          fExpExpCB.FitData(dat,iName,iTitle,"Full");
          fExpExpCB.GetPlot()->Write();
          values[0] = fExpExpCB.mSigCounts->getVal();
          hFit0[iS][iC]->SetBinContent(iB + 1, values[0]);

          base_dir->cd(Form("%s/Model2",kNames[iS].data()));
          fExpExpGaus.FitData(dat,iName,iTitle,"Full");
          fExpExpGaus.GetPlot()->Write();
          values[1] = fExpExpGaus.mSigCounts->getVal();
          hFit1[iS][iC]->SetBinContent(iB + 1, values[1]);

          base_dir->cd(Form("%s/Range1",kNames[iS].data()));
          mainModule.FitData(dat, iName, iTitle, "Range1");
          mainModule.mPlot->Write();
          values[2] = mainModule.mSigCounts->getVal();
          hFit2[iS][iC]->SetBinContent(iB + 1, values[2]);

          base_dir->cd(Form("%s/Range2",kNames[iS].data()));
          mainModule.FitData(dat, iName, iTitle, "Range2");
          mainModule.mPlot->Write();
          values[3] = mainModule.mSigCounts->getVal();
          hFit3[iS][iC]->SetBinContent(iB + 1, values[3]);
          ///TODO: improve the logic for systematics evaluation
          float syst = (TMath::MaxElement(5, values) - TMath::MinElement(5, values)) / TMath::Sqrt(12);
          hSystFit[iS][iC]->SetBinContent(iB + 1, syst / refv);

          /// TPC only raw signal
          hTPConly[iS][iC]->SetBinContent(iB + 1, tpc_histo[iS]->GetBinContent(iC + 1, iB + 1));
        }
      }
    }

    for (int iS = 0; iS < 2; ++iS) {
      for (int iC = 0; iC < n_centralities; ++iC) {
        base_dir->cd(Form("%s/Model0",kNames[iS].data()));
        hRawCounts[iS][iC]->Write();
        base_dir->cd(Form("%s/Model1",kNames[iS].data()));
        hFit0[iS][iC]->Write();
        base_dir->cd(Form("%s/Model2",kNames[iS].data()));
        hFit1[iS][iC]->Write();
        base_dir->cd(Form("%s/Range1",kNames[iS].data()));
        hFit2[iS][iC]->Write();
        base_dir->cd(Form("%s/Range2",kNames[iS].data()));
        hFit3[iS][iC]->Write();
        base_dir->cd(Form("%s/Significance",kNames[iS].data()));
        hSignificance[iS][iC]->Write();
        base_dir->cd(Form("%s/Systematics",kNames[iS].data()));
        hSystFit[iS][iC]->Write();
        base_dir->cd(Form("%s/TPConly",kNames[iS].data()));
        hTPConly[iS][iC]->Write();
      }
    }
    base_dir->Close();
  }
  output_file.Close();
}

