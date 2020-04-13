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

const bool bGeant = false;

const bool bAbsCorr = false;

const bool bAbsCorrMax = true;

const bool bNewGeantTPC = true;

const bool bNewGeantTOF = false;

int contatore = 0;

const float kNormalisation[10] = {0.999663,0.999173,0.999429,0.9988205,0.996784,0.992304,0.982828,0.953767,0.785473,0.921028};
const float kTriggerEff = 0.7448;

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

void Divide(TH1* h, TF1* func) {
  double x_low, x_high;
  func->GetRange(x_low,x_high);
  for (int i = 1; i <= h->GetNbinsX(); ++i) {
    if(h->GetBinCenter(i) < x_low || h->GetBinCenter(i) > x_high) continue;
    if(h->GetBinContent(i) == 0 || fabs(func->Eval(h->GetBinCenter(i))) < 1.e-24) continue;
    h->SetBinContent(i, h->GetBinContent(i) / func->Eval(h->GetBinCenter(i)));
    h->SetBinError(i, h->GetBinError(i) / func->Eval(h->GetBinCenter(i)));
  }
}


void Spectra() {
  /// Taking all the histograms from the MC file
  TFile signal_file(kSignalOutput.data());
  TFile data_file(kDataFilename.data());
  TFile secondaries_file((kSecondariesOutput).data());
  TFile secondaries_tpc_file(kSecondariesTPCoutput.data());
  TFile efficiency_file(kEfficiencyOutput.data());
  TFile signalloss_file(kSignalLossOutput.data());
  TFile output_file(kSpectraOutput.data(),"recreate");

  for(int iCut=0; iCut<kNsigmaVar; iCut++){
    auto cut_dir = output_file.mkdir(Form("%swidth%i",kFilterListNames.data(),iCut));
    cut_dir->cd();
    for(int iS=0; iS<2; iS++){
      auto part_dir = cut_dir->mkdir(kNames[iS].data());
      for(int iC=kCentLength; iC--;){
        auto cent_dir = part_dir->mkdir(Form("%i",iC));
        cent_dir->mkdir("TOF");
        cent_dir->mkdir("TPC");
        cent_dir->mkdir("Joined");
      }
    }
  }

  for(int iCut=0; iCut<kNshiftVar; iCut++){
    auto cut_dir = output_file.mkdir(Form("%sshift%i",kFilterListNames.data(),iCut));
    cut_dir->cd();
    for(int iS=0; iS<2; iS++){
      auto part_dir = cut_dir->mkdir(kNames[iS].data());
      for(int iC=kCentLength; iC--;){
        auto cent_dir = part_dir->mkdir(Form("%i",iC));
        cent_dir->mkdir("TOF");
        cent_dir->mkdir("TPC");
        cent_dir->mkdir("Joined");
      }
    }
  }

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

  // Old GEANT3/GEANT4 correction

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
  Requires(g3g4tpc[0],"tpcA");
  Requires(g3g4tpc[1],"tpcM");

  TCanvas* cG3G4[2] = {new TCanvas("G3G4TOF","G3G4TOF"), new TCanvas("G3G4TPC","G3G4TPC")};

  TF1* funcG3G4_tof[2] = {nullptr, nullptr};
  TF1* funcG3G4_tpc[2] = {nullptr, nullptr};
  
  auto g3g4_dir = output_file.mkdir("G3G4_dir");
  auto g3g4_tpc_dir = g3g4_dir->mkdir("TPC");
  auto g3g4_tof_dir = g3g4_dir->mkdir("TOF");

  // GEANT3 / GEANT4 correction
  for(int iS = 0; iS<2; iS++){
    // TPC
    g3g4_tpc_dir->cd();
    funcG3G4_tpc[iS] = new TF1(Form("function_g3g4_tpc_%c",kLetter[iS]),"pol0",0.6,6.);
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
    g3g4tpc[iS]->Fit(funcG3G4_tpc[iS],"RQ");
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
    // TOF
    g3g4_tof_dir->cd();
    if(iS==1) funcG3G4_tof[iS] = new TF1(Form("function_g3g4_tof_%c",kLetter[iS]),"pol0",0.9,6.);
    else {
      funcG3G4_tof[iS] = new TF1(Form("function;_g3g4_tof_%c",kLetter[iS]),"[0]+[1]*exp([2]*x+[3]*x*x)",0.9,6.);
      funcG3G4_tof[iS]->SetParLimits(3,-10.,0.);
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
    g3g4tof[iS]->Fit(funcG3G4_tof[iS],"RQ");
    g3g4tof[iS]->Write();
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
  }

  // New GEANT3/GEANT4 correction, together with TOF matching efficiency

  TFile fileG3G4new(kG3G4output.data());
  TH1F* corrG3G4newTPC[2]{
    (TH1F*)fileG3G4new.Get("corrg3g4tpcM"),
    (TH1F*)fileG3G4new.Get("corrg3g4tpcA")
  };
  Requires(corrG3G4newTPC[0],"new G3G4 correction for matter for the TPC");
  Requires(corrG3G4newTPC[1],"new G3G4 correction for antimatter for the TPC");

  TF1* funcG3G4newTPC[2]{
    (TF1*)fileG3G4new.Get("functionM_tpc"),
    (TF1*)fileG3G4new.Get("functionA_tpc")
  };
  Requires(funcG3G4newTPC[0],"new G3G4 correction function for matter for the TPC");
  Requires(funcG3G4newTPC[1],"new G3G4 correction function for antimatter for the TPC");

  TH1F* corrG3G4newTOF[2]{
    (TH1F*)fileG3G4new.Get("corrg3g4tofM"),
    (TH1F*)fileG3G4new.Get("corrg3g4tofA")
  };
  Requires(corrG3G4newTOF[0],"new G3G4 correction for matter for the TOF");
  Requires(corrG3G4newTOF[1],"new G3G4 correction for antimatter for the TOF");

  TF1* funcG3G4newTOF[2]{
    (TF1*)fileG3G4new.Get("functionM_tof"),
    (TF1*)fileG3G4new.Get("functionA_tof")
  };
  Requires(funcG3G4newTPC[0],"new G3G4 correction function for matter for the TOF");
  Requires(funcG3G4newTPC[1],"new G3G4 correction function for antimatter for the TOF");

  TH1F* hTotMeanSignalLoss[kCentLength];
  for(int iC=kCentLength; iC--;){
    hTotMeanSignalLoss[iC] = (TH1F*)signalloss_file.Get(Form("%d/hTotMeanSignalLoss_%d",iC,iC));
    Requires(hTotMeanSignalLoss[iC],Form("Missing %d/hMeanSignalLoss_%d",iC,iC));
    hTotMeanSignalLoss[iC]->SetDirectory(0);
  }

  for (auto list_key : *signal_file.GetListOfKeys()) {
    if (string(list_key->GetName()).find(kFilterListNames.data()) == string::npos) continue;
    if (string(list_key->GetName()).find("_MV") != string::npos) continue;
    if (string(list_key->GetName()).find("_pid2") != string::npos) continue;
    if (string(list_key->GetName()).find("_pid3") != string::npos) continue;
    TTList* list = (TTList*)data_file.Get(list_key->GetName());
    /// Getting the correct directories
    TDirectory* base_dir = output_file.mkdir(list_key->GetName());
    base_dir->cd();

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
      for (int iC=kCentLength; iC--;) {
        TDirectory* cent_dir = particle_dir->mkdir(Form("%d",iC));
        cent_dir->cd();

        ///Efficiencies MB
        TH1* eff_tpc_graph = (TH1*)efficiency_file.Get(Form("%s/%s/MB/effTpc%c%i",list_key->GetName(),kNames[iS].data(),kLetter[iS],0));
        TH1* eff_tof_graph = (TH1*)efficiency_file.Get(Form("%s/%s/MB/effTof%c%i",list_key->GetName(),kNames[iS].data(),kLetter[iS],0));
        Requires(eff_tpc_graph,Form("%s/effTpc%c%i",list_key->GetName(),kLetter[iS],0));
        Requires(eff_tof_graph,Form("%s/effTof%c%i",list_key->GetName(),kLetter[iS],0));
        /// Getting efficiencies
        // TH1* eff_tpc_graph = (TH1*)efficiency_file.Get(Form("%s/%s/C_%i/effTpc%c%i",list_key->GetName(),kNames[iS].data(),iC+1,kLetter[iS],iC+1));
        // TH1* eff_tof_graph = (TH1*)efficiency_file.Get(Form("%s/%s/C_%i/effTof%c%i",list_key->GetName(),kNames[iS].data(),iC+1,kLetter[iS],iC+1));
        // Requires(eff_tpc_graph,Form("%s/%s/C_%i/effTpc%c%i",list_key->GetName(),kNames[iS].data(),iC+1,kLetter[iS],iC+1));
        // Requires(eff_tof_graph,Form("%s/%s/C_%i/effTof%c%i",list_key->GetName(),kNames[iS].data(),iC+1,kLetter[iS],iC+1));

        if (secondaries_file.IsOpen()) {
          hResTFF = (TH1F*)secondaries_file.Get(Form("%s/Results/%d/hResTFF_%i",list_key->GetName(),iC,iC));
          Requires(hResTFF,Form("%s/Results/%d/hResTFF_%i",list_key->GetName(),iC,iC));
          primary_fraction = hResTFF->GetFunction("fitFrac");
          Requires(primary_fraction,"Missing primary fraction");
        }
        if (secondaries_tpc_file.IsOpen()) {
          hResTFF_TPC = (TH1F*)secondaries_tpc_file.Get(Form("%s/Results/%d/hResTFF_%i",list_key->GetName(),iC,iC));
          Requires(hResTFF_TPC,Form("%s/Results/%d/hResTFF_%i",list_key->GetName(),iC,iC));
          primary_fraction_tpc = hResTFF_TPC->GetFunction("fitFrac");
          Requires(primary_fraction_tpc,"Missing primary fraction for TPC");
        }

        /// Getting raw signals
        TH1D* rawTOF = (TH1D*)signal_file.Get(Form("%s/%s/Fits/hRawCounts%c%i",list_key->GetName(),kNames[iS].data(),kLetter[iS],iC));
        TH1D* rawTPC = (TH1D*)signal_file.Get(Form("%s/%s/TPConly/hRawCountsTPC%c%i",list_key->GetName(),kNames[iS].data(),kLetter[iS],iC));
        Requires(rawTOF,"Missing TOF raw counts");
        Requires(rawTPC,"Missing TPC raw counts");

        /// TOF corrected spectra
        TH1D* spectraTOF = new TH1D(*rawTOF);
        spectraTOF->SetTitle(";#it{p}_{T} (GeV/#it{c});#frac{1}{N_{ev}}#frac{dN}{dp_{T}d#it{y}}");
        Divide(spectraTOF,eff_tof_graph);

        TDirectory* TOF_dir = cent_dir->mkdir("TOF");
        TOF_dir->cd();

        if(bGeant){
          Divide(spectraTOF,funcG3G4_tof[iS]);
        }
        if(bAbsCorr){
          for(int iB=1; iB<=spectraTOF->GetNbinsX(); iB++){
            spectraTOF->SetBinContent(iB,spectraTOF->GetBinContent(iB)*kAbsCorr[iS]);
            spectraTOF->SetBinError(iB,spectraTOF->GetBinError(iB)*kAbsCorr[iS]);
          }
        }
        if(bAbsCorrMax){
          for(int iB=1; iB<=spectraTOF->GetNbinsX(); iB++){
            spectraTOF->SetBinContent(iB,spectraTOF->GetBinContent(iB)*kAbsCorrMax[iS]*kCorrG3G4tpc[iS]);
            spectraTOF->SetBinError(iB,spectraTOF->GetBinError(iB)*kAbsCorrMax[iS]*kCorrG3G4tpc[iS]);
          }
        }
        if(bNewGeantTOF){
          for(int iB=1; iB<=spectraTOF->GetNbinsX(); iB++){
            spectraTOF->SetBinContent(iB,spectraTOF->GetBinContent(iB)*kCorrG3G4tof[iS]);
            spectraTOF->SetBinError(iB,spectraTOF->GetBinError(iB)*kCorrG3G4tof[iS]);
          }
        }
        if(!iS){
          for(int iB=1; iB<= spectraTOF->GetNbinsX(); iB++){
            float center = spectraTOF->GetBinCenter(iB);
            float factor = 0.;
            if(center>kPtRangeMatCorrectionTOF[0] && center<kPtRangeMatCorrectionTOF[1]){
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

        //// TPC corrected spectra
        TH1D* spectraTPC = new TH1D(*rawTPC);
        spectraTPC->SetTitle(";#it{p}_{T} (GeV/#it{c});#frac{1}{N_{ev}}#frac{dN}{dp_{T}d#it{y}}");
        Divide(spectraTPC,eff_tpc_graph);

        TDirectory* TPC_dir = cent_dir->mkdir("TPC");

        if(bGeant){
          Divide(spectraTPC,funcG3G4_tpc[iS]);
        }
        if(bNewGeantTPC) {
          for(int iB=1; iB<=spectraTPC->GetNbinsX(); iB++){
            spectraTPC->SetBinContent(iB,spectraTPC->GetBinContent(iB)*kCorrG3G4tpc[iS]);
            spectraTPC->SetBinError(iB,spectraTPC->GetBinError(iB)*kCorrG3G4tpc[iS]);
          }
        }
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
        double n_norm = (iC==9 && !kUseIntegratedForMB) ? (n_sel * n_vtx / n_rec) : n_vtx;
        //printf("Centrality: %d  hNormalisationHist: %f Integral: %f\n", iC, n_vtx, vNevents[iC] );
        if(iC == 9 && !kUseIntegratedForMB){
          spectraTOF->Scale(kTriggerEff/n_norm,"width");
          spectraTPC->Scale(kTriggerEff/n_norm,"width");
        } else {
          spectraTOF->Scale(kNormalisation[iC]/n_norm,"width");
          spectraTPC->Scale(kNormalisation[iC]/n_norm,"width");
        }

        TDirectory* cut_dir = nullptr;
        TDirectory* particle_width_dir = nullptr;
        
        //// Systematics related to bin counting
        if(string(list_key->GetName())==kFilterListNames.data()){
          /// Range-widening
          TH1D* rawTOF_width = nullptr;
          TH1D* rawTPC_width = nullptr;
          TH1D* spectraTOF_width = nullptr;
          TH1D* spectraTPC_width = nullptr;
          for(int iCut=0; iCut<kNsigmaVar; iCut++){
            /// TOF
            rawTOF_width = (TH1D*)signal_file.Get(Form("%s/%s/Range_width/C_%i/hWidenRangeCounts%c%i_%i",kFilterListNames.data(),kNames[iS].data(),iC,kLetter[iS],iC,iCut));
            Requires(rawTOF_width,Form("%s/%s/Range_width/C_%i/hWidenRangeCounts%c%i_%i",kFilterListNames.data(),kNames[iS].data(),iC,kLetter[iS],iC,iCut));
            spectraTOF_width = new TH1D(*rawTOF_width);
            spectraTOF_width->SetTitle(";#it{p}_{T} (GeV/#it{c});#frac{1}{N_{ev}}#frac{dN}{dp_{T}d#it{y}}");
            Divide(spectraTOF_width,eff_tof_graph);
            /// TPC
            rawTPC_width = (TH1D*)signal_file.Get(Form("%s/%s/Range_width_TPC/C_%i/hWidenRangeCountsTPC%c%i_%i",kFilterListNames.data(),kNames[iS].data(),iC,kLetter[iS],iC,iCut));
            Requires(rawTPC_width,Form("%s/%s/Range_width_TPC/C_%i/hWidenRangeCountsTPC%c%i_%i",kFilterListNames.data(),kNames[iS].data(),iC,kLetter[iS],iC,iCut));
            spectraTPC_width = new TH1D(*rawTPC_width);
            spectraTPC_width->SetTitle(";#it{p}_{T} (GeV/#it{c});#frac{1}{N_{ev}}#frac{dN}{dp_{T}d#it{y}}");
            Divide(spectraTPC_width,eff_tpc_graph);

            ///TOF
            if(bGeant){
              Divide(spectraTOF_width,funcG3G4_tof[iS]);
            }
            if(bAbsCorr){
              for(int iB=1; iB<=spectraTOF_width->GetNbinsX(); iB++){
                spectraTOF_width->SetBinContent(iB,spectraTOF_width->GetBinContent(iB)*kAbsCorr[iS]);
                spectraTOF_width->SetBinError(iB,spectraTOF_width->GetBinError(iB)*kAbsCorr[iS]);
              }
            }
            if(bAbsCorrMax){
              for(int iB=1; iB<=spectraTOF_width->GetNbinsX(); iB++){
                spectraTOF_width->SetBinContent(iB,spectraTOF_width->GetBinContent(iB)*kAbsCorrMax[iS]*kCorrG3G4tpc[iS]);
                spectraTOF_width->SetBinError(iB,spectraTOF_width->GetBinError(iB)*kAbsCorrMax[iS]*kCorrG3G4tpc[iS]);
              }
            }
            if(bNewGeantTOF){
              for(int iB=1; iB<=spectraTOF_width->GetNbinsX(); iB++){
                spectraTOF_width->SetBinContent(iB,spectraTOF_width->GetBinContent(iB)*kCorrG3G4tof[iS]);
                spectraTOF_width->SetBinError(iB,spectraTOF_width->GetBinError(iB)*kCorrG3G4tof[iS]);
              }
            }
            if(!iS){
              for(int iB=1;iB<= spectraTOF_width->GetNbinsX(); iB++){
                float center = spectraTOF_width->GetBinCenter(iB);
                float factor = 0.;
                if(center>kPtRangeMatCorrectionTOF[0] && center<kPtRangeMatCorrectionTOF[1]){
                    factor = hResTFF->GetBinContent(hResTFF->FindBin(center));
                } else {
                    factor = primary_fraction->Eval(center);
                }
                float prev_val = spectraTOF_width->GetBinContent(iB);
                float prev_err = spectraTOF_width->GetBinError(iB);
                spectraTOF_width->SetBinContent(iB,prev_val*factor);
                spectraTOF_width->SetBinError(iB,prev_err*factor);
              }
            }
            ///TPC
            if(bGeant){
              Divide(spectraTPC_width,funcG3G4_tpc[iS]);
            }
            if(bNewGeantTPC){
              for(int iB=1; iB<=spectraTPC_width->GetNbinsX(); iB++){
                spectraTPC_width->SetBinContent(iB,spectraTPC_width->GetBinContent(iB)*kCorrG3G4tpc[iS]);
                spectraTPC_width->SetBinError(iB,spectraTPC_width->GetBinError(iB)*kCorrG3G4tpc[iS]);
              }
            }
            if(!iS){
              for(int iB=1; iB<=spectraTPC_width->GetNbinsX(); iB++){
                float center = spectraTPC_width->GetBinCenter(iB);
                float factor = 0.;
                if(center>kPtRangeMatCorrectionTPC[0] && center<kPtRangeMatCorrectionTPC[1]){
                    factor = hResTFF_TPC->GetBinContent(hResTFF_TPC->FindBin(center));
                } else {
                    factor = primary_fraction_tpc->Eval(center);
                }
                float prev_val = spectraTPC_width->GetBinContent(iB);
                float prev_err = spectraTPC_width->GetBinError(iB);
                spectraTPC_width->SetBinContent(iB,prev_val*factor);
                spectraTPC_width->SetBinError(iB,prev_err*factor);
              }
            }
            if(iC>=6){
              Multiply(spectraTPC_width,hTotMeanSignalLoss[iC]);
              Multiply(spectraTOF_width,hTotMeanSignalLoss[iC]);
            }

            if(iC == 9 && !kUseIntegratedForMB){
              spectraTOF_width->Scale(kTriggerEff/n_norm,"width");
              spectraTPC_width->Scale(kTriggerEff/n_norm,"width");
            } else {
              spectraTOF_width->Scale(kNormalisation[iC]/n_norm,"width");
              spectraTPC_width->Scale(kNormalisation[iC]/n_norm,"width");
            }

            output_file.cd(Form("%swidth%i/%s/%i/TOF",kFilterListNames.data(),iCut,kNames[iS].data(),iC));
            spectraTOF_width->Write(Form("TOFspectra%c%i",kLetter[iS],iC));
            output_file.cd(Form("%swidth%i/%s/%i/TPC",kFilterListNames.data(),iCut,kNames[iS].data(),iC));
            spectraTPC_width->Write(Form("TPCspectra%c%i",kLetter[iS],iC));
            TH1D* JoinedSpectra_width = new TH1D(Form("JoiendSpectraWidth%c%i",kLetter[iS],iC),";#it{p}_{T} (GeV/#it{c});#frac{1}{N_{ev}}#frac{dN}{dp_{T}d#it{y}}",kNPtBins,kPtBins);
            for(int iB=1; iB<=kNPtBins; iB++){
              double bin_center = JoinedSpectra_width->GetXaxis()->GetBinCenter(iB);
              if(bin_center<1.){
                int binNumber = spectraTPC_width->FindBin(bin_center);
                JoinedSpectra_width->SetBinContent(iB,spectraTPC_width->GetBinContent(binNumber));
                JoinedSpectra_width->SetBinError(iB,spectraTPC_width->GetBinError(binNumber));
              }else{
                int binNumber = spectraTOF_width->FindBin(bin_center);
                JoinedSpectra_width->SetBinContent(iB,spectraTOF_width->GetBinContent(binNumber));
                JoinedSpectra_width->SetBinError(iB,spectraTOF_width->GetBinError(binNumber));
              }
            }
            output_file.cd(Form("%swidth%i/%s/%i/Joined",kFilterListNames.data(),iCut,kNames[iS].data(),iC));
            JoinedSpectra_width->Write(Form("JoinedSpectra%c%i",kLetter[iS],iC));   
          }
          /// Range-shift
          TH1D* rawTOF_shift = nullptr;
          TH1D* rawTPC_shift = nullptr;
          TH1D* spectraTOF_shift = nullptr;
          TH1D* spectraTPC_shift = nullptr;
          for(int iCut=0; iCut<kNshiftVar; iCut++){
            /// TOF
            rawTOF_shift = (TH1D*)signal_file.Get(Form("%s/%s/Range_shift/C_%i/hShiftRangeCounts%c%i_%i",kFilterListNames.data(),kNames[iS].data(),iC,kLetter[iS],iC,iCut));
            Requires(rawTOF_shift,Form("%s/%s/Range_shift/C_%i/hShiftRangeCounts%c%i_%i",kFilterListNames.data(),kNames[iS].data(),iC,kLetter[iS],iC,iCut));
            spectraTOF_shift = new TH1D(*rawTOF_shift);
            spectraTOF_shift->SetTitle(";#it{p}_{T} (GeV/#it{c});#frac{1}{N_{ev}}#frac{dN}{dp_{T}d#it{y}}");
            Divide(spectraTOF_shift,eff_tof_graph);
            /// TPC
            rawTPC_shift = (TH1D*)signal_file.Get(Form("%s/%s/Range_shift_TPC/C_%i/hShiftRangeCountsTPC%c%i_%i",kFilterListNames.data(),kNames[iS].data(),iC,kLetter[iS],iC,iCut));
            Requires(rawTPC_shift,Form("%s/%s/Range_shift_TPC/C_%i/hShiftRangeCountsTPC%c%i_%i",kFilterListNames.data(),kNames[iS].data(),iC,kLetter[iS],iC,iCut));
            spectraTPC_shift = new TH1D(*rawTPC_shift);
            spectraTPC_shift->SetTitle(";#it{p}_{T} (GeV/#it{c});#frac{1}{N_{ev}}#frac{dN}{dp_{T}d#it{y}}");
            Divide(spectraTPC_shift,eff_tpc_graph);
            ///TOF 
            if(bGeant){
              Divide(spectraTOF_shift,funcG3G4_tof[iS]);
            }
            if(bAbsCorr){
              for(int iB=1; iB<=spectraTOF_shift->GetNbinsX(); iB++){
                spectraTOF_shift->SetBinContent(iB,spectraTOF_shift->GetBinContent(iB)*kAbsCorr[iS]);
                spectraTOF_shift->SetBinError(iB,spectraTOF_shift->GetBinError(iB)*kAbsCorr[iS]);
              }
            }
            if(bAbsCorrMax){
              for(int iB=1; iB<=spectraTOF_shift->GetNbinsX(); iB++){
                spectraTOF_shift->SetBinContent(iB,spectraTOF_shift->GetBinContent(iB)*kAbsCorrMax[iS]*kCorrG3G4tpc[iS]);
                spectraTOF_shift->SetBinError(iB,spectraTOF_shift->GetBinError(iB)*kAbsCorrMax[iS]*kCorrG3G4tpc[iS]);
              }
            }
            if(bNewGeantTOF){
              for(int iB=1; iB<=spectraTOF_shift->GetNbinsX(); iB++){
                spectraTOF_shift->SetBinContent(iB,spectraTOF_shift->GetBinContent(iB)*kCorrG3G4tof[iS]);
                spectraTOF_shift->SetBinError(iB,spectraTOF_shift->GetBinError(iB)*kCorrG3G4tof[iS]);
              }
            }
            if(!iS){
              for(int iB=1;iB<= spectraTOF_shift->GetNbinsX(); iB++){
                float center = spectraTOF_shift->GetBinCenter(iB);
                float factor = 0.;
                if(center>kPtRangeMatCorrectionTOF[0] && center<kPtRangeMatCorrectionTOF[1]){
                    factor = hResTFF->GetBinContent(hResTFF->FindBin(center));
                } else {
                    factor = primary_fraction->Eval(center);
                }
                float prev_val = spectraTOF_shift->GetBinContent(iB);
                float prev_err = spectraTOF_shift->GetBinError(iB);
                spectraTOF_shift->SetBinContent(iB,prev_val*factor);
                spectraTOF_shift->SetBinError(iB,prev_err*factor);
              }
            }
            //TPC
            if(bGeant){
              Divide(spectraTPC_shift,funcG3G4_tpc[iS]);
            }
            if(bNewGeantTPC){
              for(int iB=1; iB<=spectraTPC_shift->GetNbinsX(); iB++){
                spectraTPC_shift->SetBinContent(iB,spectraTPC_shift->GetBinContent(iB)*kCorrG3G4tpc[iS]);
                spectraTPC_shift->SetBinError(iB,spectraTPC_shift->GetBinError(iB)*kCorrG3G4tpc[iS]);
              }
            }
            if(!iS){
              for(int iB=1; iB<=spectraTPC_shift->GetNbinsX(); iB++){
                float center = spectraTPC_shift->GetBinCenter(iB);
                float factor = 0.;
                if(center>kPtRangeMatCorrectionTPC[0] && center<kPtRangeMatCorrectionTPC[1]){
                    factor = hResTFF_TPC->GetBinContent(hResTFF_TPC->FindBin(center));
                } else {
                    factor = primary_fraction_tpc->Eval(center);
                }

                float prev_val = spectraTPC_shift->GetBinContent(iB);
                float prev_err = spectraTPC_shift->GetBinError(iB);
                spectraTPC_shift->SetBinContent(iB,prev_val*factor);
                spectraTPC_shift->SetBinError(iB,prev_err*factor);
              }
            }
            if(iC>=6){
              Multiply(spectraTPC_shift,hTotMeanSignalLoss[iC]);
              Multiply(spectraTOF_shift,hTotMeanSignalLoss[iC]);
            }

            if(iC == 9 && !kUseIntegratedForMB){
              spectraTOF_shift->Scale(kTriggerEff/n_norm,"width");
              spectraTPC_width->Scale(kTriggerEff/n_norm,"width");
            } else {
              spectraTOF_shift->Scale(kNormalisation[iC]/n_norm,"width");
              spectraTPC_shift->Scale(kNormalisation[iC]/n_norm,"width");
            }

            output_file.cd(Form("%sshift%i/%s/%i/TOF",kFilterListNames.data(),iCut,kNames[iS].data(),iC));
            spectraTOF_shift->Write(Form("TOFspectra%c%i",kLetter[iS],iC));
            output_file.cd(Form("%sshift%i/%s/%i/TPC",kFilterListNames.data(),iCut,kNames[iS].data(),iC));
            spectraTPC_shift->Write(Form("TPCspectra%c%i",kLetter[iS],iC));
            TH1D* JoinedSpectra_shift = new TH1D(Form("JoiendSpectraShift%c%i",kLetter[iS],iC),";#it{p}_{T} (GeV/#it{c});#frac{1}{N_{ev}}#frac{dN}{dp_{T}d#it{y}}",kNPtBins,kPtBins);
            for(int iB=1; iB<=kNPtBins; iB++){
              double bin_center = JoinedSpectra_shift->GetXaxis()->GetBinCenter(iB);
              if(bin_center<1.){
                int binNumber = spectraTPC_shift->FindBin(bin_center);
                JoinedSpectra_shift->SetBinContent(iB,spectraTPC_shift->GetBinContent(binNumber));
                JoinedSpectra_shift->SetBinError(iB,spectraTPC_shift->GetBinError(binNumber));
              }else{
                int binNumber = spectraTOF_shift->FindBin(bin_center);
                JoinedSpectra_shift->SetBinContent(iB,spectraTOF_shift->GetBinContent(binNumber));
                JoinedSpectra_shift->SetBinError(iB,spectraTOF_shift->GetBinError(binNumber));
              }
            }
            output_file.cd(Form("%sshift%i/%s/%i/Joined",kFilterListNames.data(),iCut,kNames[iS].data(),iC));
            JoinedSpectra_shift->Write(Form("JoinedSpectra%c%i",kLetter[iS],iC));
          }
          
        }

        TOF_dir->cd();
        spectraTOF->Write(Form("TOFspectra%c%i",kLetter[iS],iC));
        TPC_dir->cd();
        spectraTPC->Write(Form("TPCspectra%c%i",kLetter[iS],iC));

        TDirectory* Joined_dir = cent_dir->mkdir("Joined");
        Joined_dir->cd();
        
        TH1D* JoinedSpectra = new TH1D(Form("JoiendSpectra%c%i",kLetter[iS],iC),";#it{p}_{T} (GeV/#it{c});#frac{1}{N_{ev}}#frac{dN}{dp_{T}d#it{y}}",kNPtBins,kPtBins);
        for(int iB=1; iB<=kNPtBins; iB++){
          double bin_center = JoinedSpectra->GetXaxis()->GetBinCenter(iB);
          if(bin_center<1.){
            int binNumber = spectraTPC->FindBin(bin_center);
            JoinedSpectra->SetBinContent(iB,spectraTPC->GetBinContent(binNumber));
            JoinedSpectra->SetBinError(iB,spectraTPC->GetBinError(binNumber));
          }else{
            int binNumber = spectraTOF->FindBin(bin_center);
            JoinedSpectra->SetBinContent(iB,spectraTOF->GetBinContent(binNumber));
            JoinedSpectra->SetBinError(iB,spectraTOF->GetBinError(binNumber));
          }
        }
        JoinedSpectra->Write(Form("JoinedSpectra%c%i",kLetter[iS],iC));     
      }
      particle_dir->Close();
    }
  }
  output_file.Close();
}
