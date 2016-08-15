#include "src/Common.h"

#include <TEfficiency.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraphAsymmErrors.h>
#include <TH2F.h>

void Efficiency() {
  /// Taking all the histograms from the MC file
  TFile input_file(kMCfilename.data());
  TFile output_file(kEfficiencyOutput.data(),"recreate");

  /// Bilding the function used to fit the efficiency distribution
  TF1 effModel("effModel","[0]+[1]*exp([2]*x)+[3]/x+[4]/sq(x)",0,10.);
  effModel.SetParLimits(0, 0.2, 0.8);
  effModel.SetParLimits(1, -5, 5.);
  effModel.SetParLimits(2, -7., 0);
  effModel.SetParLimits(3, -5.,5.);
  effModel.SetParLimits(4, -5., 0.);
  effModel.SetParameters(0.35, -0.98,-2.2,0.001,-0.01);

  for (auto&& list_key : *input_file.GetListOfKeys()) {
    /// Preliminary operation to read the list and create an output dir
    if (string(list_key->GetName()).find(kFilterListNames.data()) == string::npos) continue;
    TTList* list = (TTList*)input_file.Get(list_key->GetName());
    output_file.mkdir(list_key->GetName());
    output_file.cd(list_key->GetName());

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
      for (int iBx = 1; iBx <= fITS_TPC[0]->GetNbinsX(); ++iBx) {
        TH1D *tpc = fITS_TPC[iS]->ProjectionY(Form("tpc%i",iBx),iBx,iBx);
        TH1D *tof = fITS_TPC_TOF[iS]->ProjectionY(Form("tof%i",iBx),iBx,iBx);
        TH1D *tot = fTotal[iS]->ProjectionY(Form("tot%i",iBx),iBx,iBx);

        TEfficiency efftpcE(*tpc,*tot),efftofE(*tof,*tot);

        TGraphAsymmErrors *efftof = efftofE.CreateGraph();
        TGraphAsymmErrors *efftpc = efftpcE.CreateGraph();

        efftof->GetYaxis()->SetRangeUser(0.f,1.1f);
        efftof->GetYaxis()->SetTitle("Efficiency x Acceptance");
        efftof->SetMarkerStyle(24);
        efftof->SetMarkerSize(0.7);
        efftof->SetLineColor(kBlue);
        efftof->SetMarkerColor(kBlue);
        efftof->SetNameTitle(Form("effTof%c%i",kLetter[iS],iBx),
            Form("%4.2f-%4.2f%%;p_{T} (GeV/c);Efficiency x Acceptance",
              cent_labels[iBx - 1], cent_labels[iBx]));
        efftof->Fit(&effModel,"MQ");
        efftof->Write();

        efftpc->GetYaxis()->SetRangeUser(0.f,1.1f);
        efftpc->GetYaxis()->SetTitle("Efficiency x Acceptance");
        efftpc->SetMarkerStyle(24);
        efftpc->SetMarkerSize(0.7);
        efftpc->SetLineColor(kRed);
        efftpc->SetMarkerColor(kRed);
        efftpc->SetNameTitle(Form("effTpc%c%i",kLetter[iS],iBx),
            Form("%4.2f-%4.2f%%;p_{T} (GeV/c);Efficiency x Acceptance",
              cent_labels[iBx - 1], cent_labels[iBx]));
        efftpc->Fit(&effModel,"MQ");
        efftpc->Write();
      }
    }
  }
  output_file.Close();
}
