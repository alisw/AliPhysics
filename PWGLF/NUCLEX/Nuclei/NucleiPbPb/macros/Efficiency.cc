#include "src/Common.h"
#include "src/Utils.h"
using namespace utils;
#include "src/Plotting.h"

#include <TCanvas.h>
#include <TEfficiency.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraphAsymmErrors.h>
#include <TH2F.h>

int kCentBins[11][2] = {{1,1},{2,2},{3,3},{4,4},{5,5},{6,6},{7,7},{8,8},{9,9},{10,10},{11,11}};
int kNcentBins = 11;

TH1* DoEff(TH1* tof, TH1* tot, string name, char letter, int iBx, TArrayD& cent_labels) {
  int iC0 = kCentBins[iBx-1][0];
  int iC1 = kCentBins[iBx-1][1];
  float lmin = iBx ? cent_labels[iC0-1] : cent_labels[0];
  float lmax = iBx ? cent_labels[iC1]   : cent_labels[cent_labels.GetSize()-1];

  TH1* efftof = (TH1*)tof->Clone(Form("eff%s%c%i",name.data(),letter,iBx));
  efftof->SetTitle(Form("%3.0f - %3.0f%%",lmin,lmax));
  efftof->GetYaxis()->SetRangeUser(0.f,1.1f);
  efftof->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  efftof->GetYaxis()->SetTitle("Efficiency x Acceptance");
  efftof->SetMarkerStyle(24);
  efftof->SetMarkerSize(0.7);
  efftof->SetLineColor(kBlack);
  efftof->SetMarkerColor(kBlack);
  for (int iBin = 1; iBin <= efftof->GetNbinsX(); ++iBin) {
    double num = tof->GetBinContent(iBin);
    double den = tot->GetBinContent(iBin);
    if (std::abs(den) < 1.e-24) continue;
    double eff = num/den;
    efftof->SetBinContent(iBin, eff);
    efftof->SetBinError(iBin, std::sqrt(eff * (1. - eff) / den));
  }

  efftof->Write();
  return efftof;

}

void Efficiency() {
  /// Taking all the histograms from the MC file
  TFile input_file(kMCfilename.data());
  TFile output_file(kEfficiencyOutput.data(),"recreate");

  /// Bilding the function used to fit the efficiency distribution
  TF1 effModel("effModel","[0]+[1]*exp([2]*x)+[3]/x+[4]/sq(x)+[5]/(x*sq(x))",0,10.);
  effModel.SetParLimits(0, 0.2, 0.8);
  effModel.SetParLimits(1, -5, 5.);
  effModel.SetParLimits(2, -7., 0);
  effModel.SetParLimits(3, -5.,5.);
  effModel.SetParLimits(4, -5., 0.);
  effModel.SetParameters(0.35, -0.98,-2.2,0.001,-0.01);

  int iList = 0;
  int counter=0;
  for (auto list_key : *input_file.GetListOfKeys()) {
    /// Preliminary operation to read the list and create an output dir
    if (string(list_key->GetName()).find(kFilterListNames.data()) == string::npos) continue;
    TTList* list = (TTList*)input_file.Get(list_key->GetName());
    output_file.mkdir(list_key->GetName());
    output_file.cd(list_key->GetName());
    //if (iList) break;
    iList++;
    /// Getting all the histograms
    TH2F  *fITS_TPC[2],*fTotal[2],*fITS_TPC_TOF[2];
    for (int iS = 0; iS < 2; ++iS) {
      fTotal[iS] = dynamic_cast<TH2F*>(list->Get(Form("f%cTotal",kLetter[iS])));
      fITS_TPC[iS] = dynamic_cast<TH2F*>(list->Get(Form("f%cITS_TPC",kLetter[iS])));
      fITS_TPC_TOF[iS] = dynamic_cast<TH2F*>(list->Get(Form("f%cITS_TPC_TOF",kLetter[iS])));
    }

    /// Taking information about centrality bins
    auto n_centralities = fTotal[0]->GetNbinsX();
    auto cent_labels = *(fTotal[0]->GetXaxis()->GetXbins());

    /// Writing a reference to file, just to recover the configuration
    /// number of centrality bins and pT bins later on
    TH2F hReference("hReference",";Centrality (%);#it{p}_{T} (GeV/#it{c})",
        n_centralities,cent_labels.GetArray(),fTotal[0]->GetNbinsY(),
        fTotal[0]->GetYaxis()->GetXbins()->GetArray());
    hReference.Write();

    /// Loop to analyse the centrality bins individually
    for (int iS = 0; iS < 2; ++iS) {
      TH1D *tpcMB = fITS_TPC[iS]->ProjectionY("tpc0");
      TH1D *tofMB = fITS_TPC_TOF[iS]->ProjectionY("tof0");
      TH1D *totMB = fTotal[iS]->ProjectionY("tot0");

      TH1* effTofMB = DoEff(tofMB,totMB,"Tof",kLetter[iS],0,cent_labels);
      TH1* effTpcMB = DoEff(tpcMB,totMB,"Tpc",kLetter[iS],0,cent_labels);
      for (int iBx = 1; iBx <= kNcentBins; ++iBx) {
        TH1D *tpc = fITS_TPC[iS]->ProjectionY(Form("tpc%i",iBx),kCentBins[iBx-1][0],kCentBins[iBx-1][1]);
        TH1D *tof = fITS_TPC_TOF[iS]->ProjectionY(Form("tof%i",iBx),kCentBins[iBx-1][0],kCentBins[iBx-1][1]);
        TH1D *tot = fTotal[iS]->ProjectionY(Form("tot%i",iBx),kCentBins[iBx-1][0],kCentBins[iBx-1][1]);

        TH1* effTof = DoEff(tof,tot,"Tof",kLetter[iS],iBx,cent_labels);
        TH1* effTpc = DoEff(tpc,tot,"Tpc",kLetter[iS],iBx,cent_labels);

        tpc->Sumw2();
        tof->Sumw2();
        int iC0 = kCentBins[iBx-1][0];
        int iC1 = kCentBins[iBx-1][1];
        plotting::SetHistStyle(tpc,plotting::kSpectraColors[iBx-1]);
        tpc->SetTitle(Form("%3.0f - %3.0f%%",cent_labels[iC0-1],cent_labels[iC1]));
        tpc->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
        tpc->GetYaxis()->SetTitle("Ratio to MB");
        plotting::SetHistStyle(tof,plotting::kSpectraColors[iBx-1]);
        tof->SetTitle(Form("%3.0f - %3.0f%%",cent_labels[iC0-1],cent_labels[iC1]));
        tof->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
        tof->GetYaxis()->SetTitle("Ratio to MB");

        TH1* ratioTPC = (TH1*)effTpc->Clone(Form("Ratio2MBtpc%c%i",kLetter[iS],iBx));
        TH1* ratioTOF = (TH1*)effTof->Clone(Form("Ratio2MBtof%c%i",kLetter[iS],iBx));
        ratioTPC->Divide(effTpcMB);
        ratioTOF->Divide(effTofMB);
        ratioTPC->Write();
        ratioTOF->Write();
        tpc->Multiply(totMB);
        tpc->Divide(tot);
        tpc->Divide(tpcMB);
        tpc->Write(Form("Ratio2MBtpc%c%i",kLetter[iS],iBx));
        tof->Multiply(totMB);
        tof->Divide(tot);
        tof->Divide(tofMB);
        tof->Write(Form("Ratio2MBtof%c%i",kLetter[iS],iBx));
      }
    }
  }
  output_file.Close();
}
