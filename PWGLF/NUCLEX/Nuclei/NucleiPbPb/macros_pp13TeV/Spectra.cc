#include "src/Common.h"
#include "src/Utils.h"
using namespace utils;

#include <TDirectoryFile.h>
#include <TDirectory.h>
#include <TEfficiency.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraphAsymmErrors.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>

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
    if(i>gr->GetNbinsX()) return;
    auto center = h->GetBinCenter(i);
    int den_bin = gr->FindBin(center);
    auto den = gr->GetBinContent(den_bin);
    if (h->GetBinContent(i) == 0 || std::abs(den) < 1.e-24) continue;
    h->SetBinContent(i, h->GetBinContent(i) / den);
    h->SetBinError(i,h->GetBinError(i) / den);
  }
}

void Multiply(TH1* h, TH1* gr) {
  for (int i = 1; i <= h->GetNbinsX(); ++i) {
    if(i>gr->GetNbinsX()) return;
    auto center = h->GetBinCenter(i);
    int den_bin = gr->FindBin(center);
    auto den = gr->GetBinContent(den_bin);
    if (h->GetBinContent(i) == 0 || std::abs(den) < 1.e-24) continue;
    h->SetBinContent(i, h->GetBinContent(i) * den);
    h->SetBinError(i,h->GetBinError(i) * den);
  }
}

void Divide(TH1* h, TF1* func, TH1F* h_out) {
  double x_low, x_high;
  func->GetRange(x_low,x_high);
  for (int i = 1; i <= h->GetNbinsX(); ++i) {
    if(h->GetBinCenter(i) < x_low || h->GetBinCenter(i) > x_high) continue;
    if(h->GetBinContent(i) == 0 || fabs(func->Eval(h->GetBinCenter(i))) < 1.e-24) continue;
    h->SetBinContent(i, h->GetBinContent(i) / func->Eval(h->GetBinCenter(i)));
    h->SetBinError(i, h->GetBinError(i) / func->Eval(h->GetBinCenter(i)));
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
  TFile signalloss_file(kSignalLossOutput.data());
  TFile output_file(kSpectraOutput.data(),"recreate");

  TDirectory* dir = (TDirectory*)data_file.Get("NucleiPIDqa_default_default");
  Requires(dir, "Directory");
  TTList* norm_list = (TTList*)dir->Get(kNormalisationList.data());
  Requires(norm_list, "norm_list");
  TH1D* hCentrality = (TH1D*)norm_list->Get("Centrality_selected");
  Requires(hCentrality,"Centrality_selected");

  float kCentlLimitBin[kCentLength][2] = {{0.,1.},{1.,5.},{5.,10.},{10.,20.},{20.,30.},{30.,40.},{40.,50},{50,70.},{70.,100.},{0.,100.}};

  double vNevents[kCentLength];
  for(int iC=0; iC<kCentLength; iC++){
    int left_edge_bin = hCentrality->FindBin(kCentlLimitBin[iC][0]+0.5);
    int right_edge_bin = hCentrality->FindBin(kCentlLimitBin[iC][1]+0.5)-1;
    vNevents[iC]= hCentrality->Integral(left_edge_bin, right_edge_bin);
  }

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

  TCanvas* cG3G4[2] = {new TCanvas("G3G4TOF","G3G4TOF"), new TCanvas("G3G4TPC","G3G4TPC")};

  TH1F* hTotMeanSignalLoss[kCentLength];
  for(int iC=0; iC<kCentLength; iC++){
    hTotMeanSignalLoss[iC] = (TH1F*)signalloss_file.Get(Form("%d/hTotMeanSignalLoss_%d",iC,iC));
    Requires(hTotMeanSignalLoss[iC],Form("Missing %d/hMeanSignalLoss_%d",iC,iC));
    hTotMeanSignalLoss[iC]->SetDirectory(0);
    output_file.cd();
  }

  for (auto list_key : *signal_file.GetListOfKeys()) {
    if (string(list_key->GetName()).find(kFilterListNames.data()) == string::npos) continue;
    if (string(list_key->GetName()).find("_MV") != string::npos) continue;
    if (string(list_key->GetName()).find("_pid2") != string::npos) continue;
    if (string(list_key->GetName()).find("_pid3") != string::npos) continue;
    TTList* list = (TTList*)data_file.Get(list_key->GetName());
    /// Getting the correct directories
    TDirectory* base_dir = output_file.mkdir(list_key->GetName());
    output_file.cd(list_key->GetName());

    /// Getting the information about centrality and pT bins
    TH3F *fATOFsignal = (TH3F*)list->Get("fATOFsignal");
    Requires(fATOFsignal,"Requires fATOFsignal");
    auto n_centralities = 1;
    centAxis = fATOFsignal->GetXaxis();
    auto centralities = *(fATOFsignal->GetXaxis()->GetXbins());
    auto n_pt_bins = fATOFsignal->GetNbinsY();
    auto ptbins = *(fATOFsignal->GetYaxis()->GetXbins());
    ptAxis = fATOFsignal->GetYaxis();

    /// Getting information about normalisation
    TH2F* fNormalisationHist = (TH2F*)list->Get("fNormalisationHist");
    Requires(fNormalisationHist,"Requires fNormalisationHist");

    TF1 *primary_fraction=nullptr;
    TF1 *primary_fraction_tpc=nullptr;
    TH1F* hResTFF = nullptr;
    TH1F* hResTFF_TPC = nullptr;

    for (int iS = 0; iS < 2; ++iS) {
      TDirectory* particle_dir = base_dir->mkdir(kNames[iS].data());
      particle_dir->cd();
      for (int iC = 0; iC < kCentLength; ++iC) {
        /// Getting efficiencies
        TH1* eff_tpc_graph = (TH1*)efficiency_file.Get(Form("%s/effTpc%c%i",list_key->GetName(),kLetter[iS],iC+1));
        TH1* eff_tof_graph = (TH1*)efficiency_file.Get(Form("%s/effTof%c%i",list_key->GetName(),kLetter[iS],iC+1));
        Requires(eff_tpc_graph,Form("%s/effTpc%c%i",list_key->GetName(),kLetter[iS],0));
        Requires(eff_tof_graph,Form("%s/effTof%c%i",list_key->GetName(),kLetter[iS],0));
        /// I am not considering Secondaries here, the analysis is only at high pT
        if (secondaries_file.IsOpen()) {
          hResTFF = (TH1F*)secondaries_file.Get(Form("%s/Results/hResTFF_%i",list_key->GetName(),iC));
          Requires(hResTFF,Form("%s/Results/hResTFF_%i",list_key->GetName(),iC));
          primary_fraction = hResTFF->GetFunction("fitFrac");
          Requires(primary_fraction,"Missing primary fraction");
        }
        if (secondaries_tpc_file.IsOpen()) {
          hResTFF_TPC = (TH1F*)secondaries_tpc_file.Get(Form("%s/Results/hResTFF_%i",list_key->GetName(),iC));
          Requires(hResTFF_TPC,Form("%s/Results/hResTFF_%i",list_key->GetName(),iC));
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
        if(iS==0){
          g3g4tof[iS]->SetLineColor(kBlue);
          g3g4tof[iS]->SetMarkerColor(kBlue);
          g3g4tof[iS]->SetMarkerStyle(20);
        }
        else {
          g3g4tof[iS]->SetLineColor(kBlack);
          g3g4tof[iS]->SetMarkerColor(kBlack);
          g3g4tof[iS]->SetMarkerStyle(24);
        }
        g3g4tof[iS]->GetYaxis()->SetTitle("GEANT4 / GEANT3");
        g3g4tof[iS]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
        g3g4tof[iS]->GetYaxis()->SetRangeUser(0.7,1.1);
        g3g4tof[iS]->GetXaxis()->SetRangeUser(1.,6.);
        g3g4tof[iS]->Fit(function_g3g4_tof,"RQ");
        g3g4tof[iS]->Write();
        // for(int iB=1;iB<= spectraTOF->GetNbinsX(); iB++){
        //   float correction = (iS==0) ? 0.906530 : 0.885085;
        //   spectraTOF->SetBinContent(iB,spectraTOF->GetBinContent(iB)/correction);
        //   spectraTOF->SetBinError(iB,spectraTOF->GetBinError(iB)/correction);
        // }
        Divide(spectraTOF,function_g3g4_tof,corr_geant_tof[iS]);
        corr_geant_tof[iS]->Write();
        cG3G4[0]->cd();
        if(iS==0) g3g4tof[iS]->Draw("AP");
        else{
          g3g4tof[iS]->Draw("SAMEP");
          TLegend leg(0.72,0.17,0.90,0.33);
          leg.SetBorderSize(0.);
          leg.AddEntry(g3g4tof[0],"Matter", "PE");
          leg.AddEntry(g3g4tof[1],"Antimatter", "PE");
          leg.Draw();
          cG3G4[0]->Write();
        }
        //if (primary_fraction&&!iS) spectraTOF->Multiply(primary_fraction);
        if(!iS){
          for(int iB=1; iB<= spectraTOF->GetNbinsX(); iB++){
            float center = spectraTOF->GetBinCenter(iB);
            float factor = 0.;
            if(center>kPtRangeMatCorrection[0] && center<kPtRangeMatCorrection[1]){
              factor = hResTFF->GetBinContent(hResTFF->FindBin(center));
            }
            else{
              factor = primary_fraction->Eval(center);
            }
            float prev_val = spectraTOF->GetBinContent(iB);
            float prev_err = spectraTOF->GetBinError(iB);
            spectraTOF->SetBinContent(iB,prev_val*factor);
            spectraTOF->SetBinError(iB,prev_err*factor);
          }
        }
        TH1D* spectraTPC = new TH1D(*rawTPC);
        spectraTPC->SetTitle(";#it{p}_{T} (GeV/#it{c});#frac{1}{N_{ev}}#frac{dN}{dp_{T}d#it{y}}");
        Divide(spectraTPC,eff_tpc_graph);
        TF1 function_g3g4_tpc(Form("function_g3g4_tpc_%c",kLetter[iS]),"pol0",0.6,6.);
        if(iS==0){
          g3g4tpc[iS]->SetLineColor(kBlue);
          g3g4tpc[iS]->SetMarkerColor(kBlue);
          g3g4tpc[iS]->SetMarkerStyle(20);
        }
        else {
          g3g4tpc[iS]->SetLineColor(kBlack);
          g3g4tpc[iS]->SetMarkerColor(kBlack);
          g3g4tpc[iS]->SetMarkerStyle(24);
        }
        g3g4tpc[iS]->GetYaxis()->SetTitle("GEANT4 / GEANT3");
        g3g4tpc[iS]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
        g3g4tpc[iS]->GetYaxis()->SetRangeUser(0.7,1.1);
        g3g4tpc[iS]->GetXaxis()->SetRangeUser(0.6,6.);
        g3g4tpc[iS]->Fit(&function_g3g4_tpc,"RQ");
        g3g4tpc[iS]->Write();
        cG3G4[1]->cd();
        if(iS==0) g3g4tpc[iS]->Draw("AP");
        else{
          g3g4tpc[iS]->Draw("SAMEP");
          TLegend leg(0.72,0.17,0.90,0.33);
          leg.SetBorderSize(0.);
          leg.AddEntry(g3g4tpc[0],"Matter", "PE");
          leg.AddEntry(g3g4tpc[1],"Antimatter", "PE");
          leg.Draw();
          cG3G4[1]->Write();
        }

        Divide(spectraTPC,&function_g3g4_tpc,corr_geant_tpc[iS]);
        //Divide(spectraTOF,&function_g3g4_tpc,corr_geant_tof[iS]);
        corr_geant_tpc[iS]->Write();
        if(!iS){
          for(int iB=1; iB<= spectraTPC->GetNbinsX(); iB++){
            float center = spectraTPC->GetBinCenter(iB);
            float factor = 0.;
            if(center>kPtRangeMatCorrectionTPC[0] && center<kPtRangeMatCorrectionTPC[1]){
              factor = hResTFF_TPC->GetBinContent(hResTFF_TPC->FindBin(center));
            }
            else{
              factor = primary_fraction_tpc->Eval(center);
            }
            float prev_val = spectraTPC->GetBinContent(iB);
            float prev_err = spectraTPC->GetBinError(iB);
            spectraTPC->SetBinContent(iB,prev_val*factor);
            spectraTPC->SetBinError(iB,prev_err*factor);
          }
        }
        // if(!iS){
        //   for (int iB = spectraTPC->FindBin(kPtRangeMatCorrectionTPC[0]); iB <= spectraTPC->FindBin(kPtRangeMatCorrectionTPC[1]); ++iB){
        //     float center = spectraTPC->GetBinCenter(iB);
        //     float factor = hResTFF_TPC->GetBinContent(hResTFF_TPC->FindBin(center));
        //     float prev_val = spectraTPC->GetBinContent(iB);
        //     float prev_err = spectraTPC->GetBinError(iB);
        //     spectraTPC->SetBinContent(iB,prev_val*factor);
        //     spectraTPC->SetBinError(iB,prev_err*factor);
        //   }
        // }//spectraTPC->Multiply(primary_fraction_tpc);
        // if (primary_fraction_tpc&&!iS) spectraTPC->Multiply(primary_fraction_tpc);

        //SignalLoss correction
        if(iC>=6){
          Multiply(spectraTPC,hTotMeanSignalLoss[iC]);
          Multiply(spectraTOF,hTotMeanSignalLoss[iC]);
        }

        /// normalisation
        TH1F* fNormMultClass = (TH1F*)fNormalisationHist->ProjectionY(Form("fNormMultClass"),kCentBinsArray[iC][0],kCentBinsArray[iC][1]);
        double n_vtx = fNormMultClass->GetBinContent(4);
        double n_rec = fNormMultClass->GetBinContent(3);
        double n_sel = fNormMultClass->GetBinContent(2); // Number of selected events
        double n_norm = n_sel * n_vtx / n_rec;
        //printf("Centrality: %d  hNormalisationHist: %f Integral: %f\n", iC, n_vtx, vNevents[iC] );
        spectraTOF->Scale(1. / n_norm,"width"); //0.7448
        spectraTPC->Scale(1. / n_norm,"width");
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
