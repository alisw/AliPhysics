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

void FitSystematics() {
  /// Suppressing the output
  RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR) ;

  /// Taking all the histograms from the MC file
  TFile input_file(kDataFilename.data());
  TFile output_file(kFitSystematicsOutput.data(),"recreate");

  /// Setting up the fitting environment
  RooRealVar m("dm2","m^{2} - m^{2}_{PDG}",-2.,2.5,"GeV^{2}/c^{4}");
  m.setBins(10000,"cache");
  m.setRange("Range1",-1.9, 2.0);
  m.setRange("Range2",-1.8, 1.8);
  m.setRange("Full", -2., 2.5);
  //
  FitExpExpCB fExpExpCB(&m);
  FitExpExpGaus fExpExpGaus(&m);
  FitExpExpTailGaus fExpExpTailGaus(&m);
  FitExpPolTailGaus fExpPolTailGaus(&m);


  TTList* list = (TTList*)input_file.Get(kFilterListNames.data());
  TDirectory* base_dir = output_file.GetDirectory("");

  /// Taking all the necessary histogram to perform the analysis
  TH3F *fATOFsignal = (TH3F*)list->Get("fATOFsignal");
  TH3F *fMTOFsignal = (TH3F*)list->Get("fMTOFsignal");

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
  fMTOFsignal->Add(fATOFsignal);
  TH3F* tof_histo = fMTOFsignal;

  /// Build arrays to analyse all the centrality classes for
  /// both deuteron and anti-deuterons. Complicate stuff just for fun.
  TH1D* hRawCounts = new TH1D("hRawCounts","; p_{T}(GeV/c); RawCounts",n_pt_bins,pt_labels.GetArray());
  TH1D* hSystFit = new TH1D("hSystFit","; p_{T}(GeV/c); Systematics from fit",n_pt_bins,pt_labels.GetArray());
  TH1D* hFit0 = new TH1D("hFit1",";p_{T} GeV/c;Bias",n_pt_bins,pt_labels.GetArray());
  TH1D* hFit1 = new TH1D("hFit2",";p_{T} GeV/c;Bias",n_pt_bins,pt_labels.GetArray());
  TH1D* hFit2 = new TH1D("hFit3",";p_{T} GeV/c;Bias",n_pt_bins,pt_labels.GetArray());
  TH1D* hFit3 = new TH1D("hFit4",";p_{T} GeV/c;Bias",n_pt_bins,pt_labels.GetArray());
  TH1D* hTail = new TH1D("hTail",";p_{T} GeV/c;Tail",n_pt_bins,pt_labels.GetArray());
  TH1D* hSigma = new TH1D("hSigma",";p_{T} GeV/c;Sigma",n_pt_bins,pt_labels.GetArray());

  /// Creating the directories to be used to store the results
  TDirectory* dir = base_dir;
  dir->cd();
  dir->mkdir("Model0");
  dir->mkdir("Model1");
  dir->mkdir("Model2");
  dir->mkdir("Range1");
  dir->mkdir("Range2");
  dir->mkdir("Systematics");

  for (int iB = 0; iB < n_pt_bins; ++iB) {
    if (pt_axis->GetBinCenter(iB+1) < kPtRange[0] ||
        pt_axis->GetBinCenter(iB+1) > kPtRange[1])
      continue;
    fExpExpGaus.SetHighPt(pt_axis->GetBinCenter(iB+1) > 5.);
    fExpExpCB.SetHighPt(pt_axis->GetBinCenter(iB+1) > 5.);
    fExpExpTailGaus.SetHighPt(pt_axis->GetBinCenter(iB+1) > 5.);
    fExpPolTailGaus.SetHighPt(pt_axis->GetBinCenter(iB+1) > 5.);

    TString iTitle = Form("%1.1f #leq p_{T} < %1.1f",pt_labels[iB],pt_labels[iB + 1]);
    TString iName = Form("d%i",iB);

    TH1D *dat = tof_histo->ProjectionZ(Form("dataM%i",iB),1,n_centralities, iB + 1, iB + 1);
    RooDataHist data("data","data",RooArgList(m),Import(*dat));

    base_dir->cd("Model0");
    FitModule &mainModule = static_cast<FitModule&>(fExpExpTailGaus);

    mainModule.FitData(dat, iName, iTitle, "Full");
    mainModule.mPlot->Write();

    const float sVal = mainModule.mSigCounts->getVal();
    const float sErr = mainModule.mSigCounts->getError();
    hRawCounts->SetBinContent(iB + 1, sVal);
    hRawCounts->SetBinError(iB + 1, sErr);

    const float sigmaVal = mainModule.mSigma->getVal();
    const float sigmaErr = mainModule.mSigma->getError();
    hSigma->SetBinContent(iB + 1, sigmaVal);
    hSigma->SetBinError(iB + 1, sigmaErr);

    const float tailVal = fExpExpTailGaus.mAlpha0->getVal();
    const float tailErr = fExpExpTailGaus.mAlpha0->getError();
    hTail->SetBinContent(iB + 1, tailVal);
    hTail->SetBinError(iB + 1, tailErr);

    const float refv = sVal;
    const float refe = sErr;
    const float chi2 = mainModule.mChi2;
    vector<float> values(1,refv);

    auto zTest = [&refv,&refe] (FitModule& mod) {
      return ((refv - mod.mSigCounts->getVal()) * (refv - mod.mSigCounts->getVal()) / (refe * refe) < 1) || !kUseBarlowFit;
    };

    base_dir->cd("Model1");
    fExpPolTailGaus.FitData(dat,iName,iTitle,"Full");
    fExpPolTailGaus.GetPlot()->Write();
    if (fExpPolTailGaus.mChi2 < 2 * chi2 && zTest(fExpPolTailGaus)) values.push_back(fExpPolTailGaus.mSigCounts->getVal());
    hFit0->SetBinContent(iB + 1, values.back());

    base_dir->cd("Model2");
    fExpExpGaus.FitData(dat,iName,iTitle,"Full");
    fExpExpGaus.GetPlot()->Write();
    if (fExpExpGaus.mChi2 < 2 * chi2 && zTest(fExpExpGaus)) values.push_back(fExpExpGaus.mSigCounts->getVal());
    hFit1->SetBinContent(iB + 1, values.back());

    base_dir->cd("Range1");
    mainModule.FitData(dat, iName, iTitle, "Range1");
    mainModule.mPlot->Write();
    if (mainModule.mChi2 < 2 * chi2 && zTest(mainModule)) values.push_back(mainModule.mSigCounts->getVal());
    hFit2->SetBinContent(iB + 1, values.back());

    base_dir->cd("Range2");
    mainModule.FitData(dat, iName, iTitle, "Range2");
    mainModule.mPlot->Write();
    if (mainModule.mChi2 < 2 * chi2 && zTest(mainModule)) values.push_back(mainModule.mSigCounts->getVal());
    hFit3->SetBinContent(iB + 1, values.back());
    ///TODO: improve the logic for systematics evaluation
    float syst = TMath::RMS(values.size(),values.data());//(TMath::MaxElement(values.size(), values.data()) - TMath::MinElement(values.size(), values.data())) / TMath::Sqrt(12);
    hSystFit->SetBinContent(iB + 1, syst / refv);

  }

  base_dir->mkdir("Parameters")->cd();
  hTail->Write();
  hSigma->Write();
  base_dir->cd("Model0");
  hRawCounts->Write();
  base_dir->cd("Model1");
  hFit0->Write();
  base_dir->cd("Model2");
  hFit1->Write();
  base_dir->cd("Range1");
  hFit2->Write();
  base_dir->cd("Range2");
  hFit3->Write();
  base_dir->cd("Systematics");
  hSystFit->Write();

  base_dir->Close();
  output_file.Close();
}

