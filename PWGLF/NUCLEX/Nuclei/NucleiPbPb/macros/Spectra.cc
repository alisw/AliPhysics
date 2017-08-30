#include "src/Common.h"
#include "src/Utils.h"
using namespace utils;

#include <TDirectoryFile.h>
#include <TEfficiency.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraphAsymmErrors.h>
#include <TH2F.h>
using namespace std;

void Divide(TH1* h, TGraphAsymmErrors* gr, bool printStuff = false) {
  for (int i = 1; i <= h->GetNbinsX(); ++i) {
    if (h->GetBinContent(i) == 0 || fabs(gr->Eval(h->GetBinCenter(i))) < 1.e-24) continue;
    if(printStuff){
      cout << "*************************************************************" << endl;
      cout << "                     Pt : " << h->GetBinCenter(i) << endl;
      cout << "*************************************************************" << endl;
      cout << "Prima : " << h->GetBinContent(i) << endl;
    }
    h->SetBinContent(i, h->GetBinContent(i) / gr->Eval(h->GetBinCenter(i)));
    if(printStuff){
      cout << "Correzione : " << gr->Eval(h->GetBinCenter(i)) << endl;
      cout << "Dopo : " << h->GetBinContent(i) << endl;
    }
  }
}

void Divide(TH1* h, TH1* gr) {
  for (int i = 1; i <= h->GetNbinsX(); ++i) {
    auto center = h->GetBinCenter(i);
    int den_bin = gr->FindBin(center);
    auto den = gr->GetBinContent(den_bin);
    if (h->GetBinContent(i) == 0 || std::abs(den) < 1.e-24) continue;
    h->SetBinContent(i, h->GetBinContent(i) / den);
    h->SetBinError(i,h->GetBinError(i) / den);
  }
}

void Divide(TH1* h, TF1* func, TH1F* h_out) {
  double x_low, x_high;
  func->GetRange(x_low,x_high);
  for (int i = 1; i <= h->GetNbinsX(); ++i) {
    if(h->GetBinCenter(i) < x_low || h->GetBinCenter(i) > x_high) continue;
    if(h->GetBinContent(i) == 0 || fabs(func->Eval(h->GetBinCenter(i))) < 1.e-24) continue;
    h->SetBinContent(i, h->GetBinContent(i) / func->Eval(h->GetBinCenter(i)));
    h_out->SetBinContent(i, func->Eval(h->GetBinCenter(i)));
  }
}


void Spectra() {
  /// Taking all the histograms from the MC file
  TFile signal_file(kSignalOutput.data());
  TFile data_file(kDataFilename.data());
  TFile secondaries_file(kSecondariesOutput.data());
  TFile secondaries_tpc_file(kSecondariesTPCoutput.data());
  TFile efficiency_file(kEfficiencyOutput.data());
  TFile output_file(kSpectraOutput.data(),"recreate");

  TAxis *centAxis=nullptr, *ptAxis=nullptr;

  TH1F* corr_geant_tof[2];
  TH1F* corr_geant_tpc[2];

  string g3g4_path = kBaseOutputDir + "G3G4.root";
  TFile g3g4_file(g3g4_path.data());
  TGraphAsymmErrors* g3g4tof[2]{
    (TGraphAsymmErrors*)g3g4_file.Get("mCorrTof"),
    (TGraphAsymmErrors*)g3g4_file.Get("aCorrTof")
  };
  Requires(g3g4tof[0],"tofA");
  Requires(g3g4tof[1],"tofM");
  TGraphAsymmErrors* g3g4tpc[2]{
    (TGraphAsymmErrors*)g3g4_file.Get("mCorrTpc"),
    (TGraphAsymmErrors*)g3g4_file.Get("aCorrTpc")
  };
  Requires(g3g4tpc[0],"tofA");
  Requires(g3g4tpc[1],"tofM");
  TTList* norm_list_tot = (TTList*)data_file.Get(kNormalisationList.data());
  Requires(norm_list_tot,"norm_list_tot");
  TH1I* fNormalisationHist = (TH1I*)norm_list_tot->Get("fNormalisationHist");
  if (!fNormalisationHist) {
    cout << "Missing normalisation to histogram " << endl;
    return;
  }
  double n_vtx = fNormalisationHist->GetBinContent(3);
  double n_rec = fNormalisationHist->GetBinContent(2);
  double n_sel = fNormalisationHist->GetBinContent(1);
  double n_norm = n_sel * n_vtx / n_rec;
  TTList* norm_list = (TTList*)data_file.Get(kFilterListNames.data());
  Requires(norm_list,"norm_list");
  TH1F* number_of_events = (TH1F*)norm_list->Get("fCentralityClasses");
  if (!number_of_events) {
    cout << "Missing normalisation to the number of events " << endl;
    return;
  }

  for (auto list_key : *signal_file.GetListOfKeys()) {
    if (string(list_key->GetName()).find(kFilterListNames.data()) == string::npos) continue;
    /// Getting the correct directories
    TDirectory* base_dir = output_file.mkdir(list_key->GetName());
    output_file.cd(list_key->GetName());

    /// Getting the information about centrality and pT bins
    TH2F* hReference = (TH2F*)efficiency_file.Get(Form("%s/hReference",list_key->GetName()));
    Requires(hReference,"Reference plot");
    auto n_centralities = 1;
    centAxis = hReference->GetXaxis();
    auto centralities = *(hReference->GetXaxis()->GetXbins());
    auto n_pt_bins = hReference->GetNbinsY();
    auto ptbins = *(hReference->GetYaxis()->GetXbins());
    ptAxis = hReference->GetYaxis();

    TF1 *primary_fraction=nullptr;
    TF1 *primary_fraction_tpc=nullptr;

    for (int iS = 0; iS < 2; ++iS) {
      TDirectory* particle_dir = base_dir->mkdir(kNames[iS].data());
      particle_dir->cd();
      //if(string(list_key->GetName())=="mpuccio_deuterons_") cout << kNames[iS].data() << endl;
      for (int iC = 0; iC < n_centralities; ++iC) {
        /// Getting efficiencies
        TH1* eff_tpc_graph = (TH1*)efficiency_file.Get(Form("%s/effTpc%c%i",list_key->GetName(),kLetter[iS],iC));
        TH1* eff_tof_graph = (TH1*)efficiency_file.Get(Form("%s/effTof%c%i",list_key->GetName(),kLetter[iS],iC));
        Requires(eff_tpc_graph,"eff_tpc_graph");
        Requires(eff_tof_graph,"eff_tof_graph");
        // TF1* eff_tpc_func = eff_tpc_graph->GetFunction("effModel");
        // TF1* eff_tof_func = eff_tof_graph->GetFunction("effModel");
        /// I am not considering Secondaries here, the analysis is only at high pT
        if (secondaries_file.IsOpen()) {
          TH1* hResTFF = (TH1*)secondaries_file.Get(Form("%s/Results/hResTFF",list_key->GetName()));
          Requires(hResTFF,Form("%s/Results/hResTFF",list_key->GetName()));
          primary_fraction = hResTFF->GetFunction("fitFrac");
          Requires(primary_fraction,"Missing primary fraction");
        }
        if (secondaries_tpc_file.IsOpen()) {
          TH1* hResTFF_TPC = (TH1*)secondaries_tpc_file.Get(Form("%s/Results/hResTFF_TPC",list_key->GetName()));
          Requires(hResTFF_TPC,Form("%s/Results/hResTFF_TPC",list_key->GetName()));
          primary_fraction_tpc = hResTFF_TPC->GetFunction("fitFrac");
          Requires(primary_fraction_tpc,"Missing primary fraction for TPC");
        }

        /// Getting raw signals
        TH1D* rawTOF = (TH1D*)signal_file.Get(Form("%s/%s/Fits/hRawCounts%c%i",list_key->GetName(),kNames[iS].data(),kLetter[iS],iC));
        TH1D* rawTPC = (TH1D*)signal_file.Get(Form("%s/%s/TPConly/hTPConly%c%i",list_key->GetName(),kNames[iS].data(),kLetter[iS],iC));
        Requires(rawTOF,"Missing TOF raw counts");
        Requires(rawTPC,"Missing TPC raw counts");

        corr_geant_tof[iS] = (TH1F*)rawTOF->Clone(Form("corr_geant_tof_%c",kLetter[iS]));
        corr_geant_tof[iS]->GetYaxis()->SetTitle("correction");
        corr_geant_tof[iS]->Reset();

        corr_geant_tpc[iS] = (TH1F*)rawTPC->Clone(Form("corr_geant_tpc_%c",kLetter[iS]));
        corr_geant_tpc[iS]->GetYaxis()->SetTitle("correction");
        corr_geant_tpc[iS]->Reset();

        /// Creating spectra histograms
        TH1D* spectraTOF = new TH1D(*rawTOF);
        spectraTOF->SetTitle(";#it{p}_{T} (GeV/#it{c});#frac{1}{N_{ev}}#frac{dN}{dp_{T}d#it{y}}");
        Divide(spectraTOF,eff_tof_graph);
        TF1* function_g3g4_tof;
        if(iS==1) function_g3g4_tof = new TF1(Form("function_g3g4_tof_%c",kLetter[iS]),"pol0",0.9,6.);
        else {
          function_g3g4_tof = new TF1(Form("function;_g3g4_tof_%c",kLetter[iS]),"[0]+[1]*exp([2]*x+[3]*x*x)",0.9,6.);
          function_g3g4_tof->SetParLimits(3,-10.,0.);
        }
        g3g4tof[iS]->SetLineColor(kBlue);
        g3g4tof[iS]->SetMarkerColor(kBlue);
        g3g4tof[iS]->Fit(function_g3g4_tof,"R");
        g3g4tof[iS]->Write();
        Divide(spectraTOF,function_g3g4_tof,corr_geant_tof[iS]);
        corr_geant_tof[iS]->Write();
        if (primary_fraction&&!iS) spectraTOF->Multiply(primary_fraction);
        TH1D* spectraTPC = new TH1D(*rawTPC);
        spectraTPC->SetTitle(";#it{p}_{T} (GeV/#it{c});#frac{1}{N_{ev}}#frac{dN}{dp_{T}d#it{y}}");
        Divide(spectraTPC,eff_tpc_graph);
        TF1 function_g3g4_tpc(Form("function_g3g4_tpc_%c",kLetter[iS]),"pol0",0.6,6.);
        g3g4tpc[iS]->SetLineColor(kBlue);
        g3g4tpc[iS]->SetMarkerColor(kBlue);
        g3g4tpc[iS]->Fit(&function_g3g4_tpc,"R");
        g3g4tpc[iS]->Write();
        Divide(spectraTPC,&function_g3g4_tpc,corr_geant_tpc[iS]);
        corr_geant_tpc[iS]->Write();
        if (primary_fraction_tpc&&!iS) spectraTPC->Multiply(primary_fraction_tpc);
        //spectraTOF->Scale(0.7448 / n_norm,"width");
        //spectraTPC->Scale(0.7448 / n_norm,"width");
        cout << "****************************************************************" << endl << endl;
        cout << "number_of_events : " << number_of_events->GetEntries() <<endl;
        cout << "n_sel : " << n_sel << endl;
        cout << "n_rec : " << n_rec << endl;
        cout << "n_vtx : " << n_vtx << endl;
        cout << "n_norm : " << n_norm << endl;
        spectraTOF->Scale(1. / number_of_events->GetEntries(),"width");
        spectraTPC->Scale(1. / number_of_events->GetEntries(),"width");
        spectraTOF->GetXaxis()->SetRange(1,15);

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
