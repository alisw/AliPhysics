#include "src/Common.h"

#include <TDirectoryFile.h>
#include <TEfficiency.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraphAsymmErrors.h>
#include <TH2F.h>

void Normalise(TH1D* h) {
  for (int i = 1; i <= h->GetNbinsX(); ++i) {
    h->SetBinContent(i, h->GetBinContent(i) / h->GetBinWidth(i));
  }
}

void Divide(TH1D* h, TGraphAsymmErrors* gr) {
  for (int i = 1; i <= h->GetNbinsX(); ++i) {
    if (h->GetBinContent(i) == 0 || fabs(gr->Eval(h->GetBinCenter(i))) < 1.e-24) continue;
    h->SetBinContent(i, h->GetBinContent(i) / gr->Eval(h->GetBinCenter(i)));
  }
}

void Spectra() {
  /// Taking all the histograms from the MC file
  TFile signal_file(kSignalOutput.data());
  TFile data_file(kDataFilename.data());
  TFile secondaries_file(kSecondariesOutput.data());
  TFile efficiency_file(kEfficiencyOutput.data());
  TFile output_file(kSpectraOutput.data(),"recreate");

  TAxis *centAxis,*ptAxis;

  TFile g3g4_file("~/Desktop/Repositories/root-files/deuterons/G3G4.root");
  TGraphAsymmErrors* g3g4tof[2]{
    (TGraphAsymmErrors*)g3g4_file.Get("mCorrTof"),
    (TGraphAsymmErrors*)g3g4_file.Get("aCorrTof")
  };
  Requires(g3g4tof[0],"tofA");
  Requires(g3g4tof[1],"tofM");

  TTList* norm_list = (TTList*)data_file.Get(kFilterListNames.data());
  TH1F* number_of_events = (TH1F*)norm_list->Get("fCentralityClasses");
  if (!number_of_events) {
    cout << "Missing normalisation to the number of events " << endl;
    return;
  }
  for (auto&& list_key : *signal_file.GetListOfKeys()) {
    if (string(list_key->GetName()).find(kFilterListNames.data()) == string::npos) continue;

    /// Getting the correct directories
    TDirectory* base_dir = output_file.mkdir(list_key->GetName());
    output_file.cd(list_key->GetName());

    /// Getting the information about centrality and pT bins
    TH2F* hReference = (TH2F*)efficiency_file.Get(Form("%s/hReference",list_key->GetName()));
    auto n_centralities = hReference->GetNbinsX();
    centAxis = hReference->GetXaxis();
    auto centralities = *(hReference->GetXaxis()->GetXbins());
    auto n_pt_bins = hReference->GetNbinsY();
    auto ptbins = *(hReference->GetYaxis()->GetXbins());
    ptAxis = hReference->GetYaxis();

    TF1 *primary_fraction;

    for (int iS = 0; iS < 2; ++iS) {
      TDirectory* particle_dir = base_dir->mkdir(kNames[iS].data());
      particle_dir->cd();
      for (int iC = 0; iC < n_centralities; ++iC) {
        /// Getting efficiencies
        TGraphAsymmErrors* eff_tpc_graph = (TGraphAsymmErrors*)efficiency_file.Get(Form("%s/effTpc%c%i",list_key->GetName(),kLetter[iS],iC+1));
        TGraphAsymmErrors* eff_tof_graph = (TGraphAsymmErrors*)efficiency_file.Get(Form("%s/effTof%c%i",list_key->GetName(),kLetter[iS],iC+1));
        TF1* eff_tpc_func = eff_tpc_graph->GetFunction("effModel");
        TF1* eff_tof_func = eff_tof_graph->GetFunction("effModel");
        /// I am not considering Secondaries here, the analysis is only at high pT
        if (secondaries_file.IsOpen()) {
          TH1F* hResTFF = (TH1F*)secondaries_file.Get(Form("%s/Results/hResTFF%i",list_key->GetName(),iC+1));
          Requires(hResTFF,Form("%s/Results/hResTFF%i",list_key->GetName(),iC+1));
          primary_fraction = hResTFF->GetFunction("fitFrac");
          Requires(primary_fraction,"Missing primary fraction");
        }

        /// Getting raw signals
        TH1D* rawTOF = (TH1D*)signal_file.Get(Form("%s/%s/Fits/hRawCounts%c%i",list_key->GetName(),kNames[iS].data(),kLetter[iS],iC));
        TH1D* rawTPC = (TH1D*)signal_file.Get(Form("%s/%s/TPConly/hTPConly%c%i",list_key->GetName(),kNames[iS].data(),kLetter[iS],iC));
        Requires(rawTOF,"Missing TOF raw counts");
        Requires(rawTPC,"Missing TPC raw counts");

        /// Creating spectra histograms
        TH1D* spectraTOF = new TH1D(*rawTOF);
        spectraTOF->SetTitle(";#it{p}_{T} (GeV/#it{c});#frac{1}{N_{ev}}#frac{dN}{dp_{T}d{y}}");
        Normalise(spectraTOF);
        Divide(spectraTOF,eff_tof_graph);
        Divide(spectraTOF,g3g4tof[iS]);
        if (primary_fraction&&!iS) spectraTOF->Multiply(primary_fraction);
        TH1D* spectraTPC = new TH1D(*rawTPC);
        spectraTPC->SetTitle(";#it{p}_{T} (GeV/#it{c});#frac{1}{N_{ev}}#frac{dN}{dp_{T}d{y}}");
        Normalise(spectraTPC);
        Divide(spectraTPC,eff_tpc_graph);

        spectraTOF->Scale(1. / number_of_events->GetBinContent(iC + 1));
        spectraTPC->Scale(1. / number_of_events->GetBinContent(iC + 1));

        spectraTOF->Write(Form("TOFspectra%i",iC));
        spectraTPC->Write(Form("TPCspectra%i",iC));
      }
      particle_dir->Close();
    }
  }

  output_file.cd();
  ptAxis->Write("pt");
  centAxis->Write("centrality");
  output_file.Close();
}
