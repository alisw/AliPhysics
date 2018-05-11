#include "src/Common.h"
#include "src/Plotting.h"
#include "src/Utils.h"
using namespace utils;

#include <map>
#include <vector>
#include <array>
using std::array;
#include <memory>

#include <Riostream.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TList.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TLine.h>
#include <TF1.h>

#include <AliPID.h>
#include <AliPWGFunc.h>

#include <RooRealVar.h>
#include <RooArgSet.h>
#include <RooDataSet.h>
#include <RooPlot.h>
#include <RooGaussian.h>
#include <TVirtualFitter.h>
#include <TPaveText.h>

const char* kPrefix[2] = {"","anti"};

const int kScaleFactor[10] = {256,128,64,32,16,8,4,2,1,1};

void Final(const char* fitFunctionName = "LevyTsallis") {
  TFile spectra_file(kSpectraOutput.data());
  TFile TOF_systematics_file(kSystematicsOutput.data());
  TFile TPC_systematics_file(kSystematicsOutputTPC.data());
  TFile final_file(kFinalOutput.data(),"recreate");

  const double pt_bin_limits[16] = {0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.4,1.6,1.8,2.0,2.2,2.6,3.0,3.4,3.8};//,4.4};
  const int n_pt_bins = 15;

  const int n_centralities = kCentLength;

  const char* syst_names[5] = {"cutsyst","countsyst","abssyst","matsyst","shiftsyst"};

  enum syst_enum { cutsyst, countsyst, abssyst, matsyst, shiftsyst};

  //histograms fot TOF analysis
  TH1F* stat_tof[2][n_centralities];
  TH1F* syst_tof[2][n_centralities];
  TH1F* vec_syst_tof[2][n_centralities][5];
  //histograms for TPC analysis
  TH1F* stat_tpc[2][n_centralities];
  TH1F* syst_tpc[2][n_centralities];
  TH1F* vec_syst_tpc[2][n_centralities][5];
  //inclusive histogram with both tpc and TOF
  TH1F* stat_all[2][n_centralities];
  TH1F* syst_all[2][n_centralities];

  TH1F* all[2][n_centralities];

  TH1F* pull_fit_data[2][n_centralities];
  TH1F* ratio_fit_data[2][n_centralities];

  TH1F* comp[2][n_centralities];

  for (int iS = 0; iS < 2; ++iS) {
    TDirectory* s_dir = final_file.mkdir(kNames[iS].data());
    for (int iC = 0; iC < n_centralities; ++iC) {
      TDirectory *c_dir = s_dir->mkdir(to_string(iC).data());
      c_dir->cd();
      pull_fit_data[iS][iC] = new TH1F(Form("pull_fit_data_%c_%d",kLetter[iS],iC),Form("%2.0f - %2.0f %%;#it{p}_{T} (GeV/#it{c});#frac{fit - data}{#sigma_{data}}",kCentLabels[iC][0],kCentLabels[iC][1]),n_pt_bins,pt_bin_limits);
      pull_fit_data[iS][iC]->GetYaxis()->SetRangeUser(-3.,3.);
      ratio_fit_data[iS][iC] = new TH1F(Form("ratio_fit_data_%c_%d",kLetter[iS],iC),Form("%2.0f - %2.0f %%;#it{p}_{T} (GeV/#it{c});data/fit",kCentLabels[iC][0],kCentLabels[iC][1]),n_pt_bins,pt_bin_limits);
      ratio_fit_data[iS][iC]->GetYaxis()->SetRangeUser(0.,2.);
      TH1F* totsyst_tmp = (TH1F*)TOF_systematics_file.Get((to_string(iC) + "/" + kNames[iS] + "/totsyst").data());
      Requires(totsyst_tmp, "Missing totsyt");
      TH1F* totsyst = (TH1F*)totsyst_tmp->Rebin(n_pt_bins,Form("totsyst_%d_%d",iS,iC),pt_bin_limits);
      TH1F* totsyst_tpc_tmp = (TH1F*)TPC_systematics_file.Get((to_string(iC) + "/" + kNames[iS] + "/totsyst_tpc").data());
      Requires(totsyst_tpc_tmp, "Missing totsyt_tpc");
      TH1F* totsyst_tpc = (TH1F*)totsyst_tpc_tmp->Rebin(n_pt_bins,Form("totsyst_tpc_%d_%d",iS,iC),pt_bin_limits);
      string tof_basepath = kFilterListNames + "/" + kNames[iS] + "/TOFspectra" + to_string(iC);
      TH1F* spectra_tof_tmp  = (TH1F*)spectra_file.Get(tof_basepath.data());
      Requires(spectra_tof_tmp,tof_basepath.data());
      spectra_tof_tmp->Write("spectra_tof_tmp");
      stat_tof[iS][iC] =(TH1F*)spectra_tof_tmp->Rebin(n_pt_bins,Form("stat_%d_%d",iS,iC),pt_bin_limits);
      auto ptAxis = stat_tof[iS][iC]->GetXaxis();
      string tpc_basepath = kFilterListNames + "/" + kNames[iS] + "/TPCspectra" + to_string(iC);
      TH1F* spectra_tpc_tmp  = (TH1F*)spectra_file.Get(tpc_basepath.data());
      Requires(spectra_tpc_tmp,tpc_basepath.data());
      spectra_tpc_tmp->Write("spectra_tpc_tmp");
      stat_tpc[iS][iC] = (TH1F*)spectra_tpc_tmp->Rebin(n_pt_bins,Form("stat_tpc_%d_%d",iS,iC),pt_bin_limits);
      syst_tof[iS][iC]  = (TH1F*)stat_tpc[iS][iC]->Clone(("syst_tof" + to_string(iC)).data());
      syst_tof[iS][iC]->Reset();
      syst_tpc[iS][iC]  = (TH1F*)stat_tpc[iS][iC]->Clone(("syst_tpc" + to_string(iC)).data());
      syst_tpc[iS][iC]->Reset();
      stat_all[iS][iC] = (TH1F*)stat_tpc[iS][iC]->Clone(("stat_all" + to_string(iC)).data());
      stat_all[iS][iC]->Reset();
      syst_all[iS][iC] = (TH1F*)stat_tpc[iS][iC]->Clone(("syst_all" + to_string(iC)).data());
      syst_all[iS][iC]->Reset();
      all[iS][iC] = (TH1F*)stat_tpc[iS][iC]->Clone(("all" + to_string(iC)).data());
      all[iS][iC]->Reset();
      //syst_all[iS][iC] = (TH1F*)totsyst->Clone(("syst_all" + to_string(iC)).data());
      //syst_all[iS][iC]->Reset();
      for(int iSyst=0; iSyst<5; iSyst++){
        TH1F* partial_syst_tmp = (TH1F*) TOF_systematics_file.Get((to_string(iC) + "/" + kNames[iS] + "/" + syst_names[iSyst]).data());
        vec_syst_tof[iS][iC][iSyst] = (TH1F*) partial_syst_tmp->Rebin(n_pt_bins,Form("%s_%d_%d",syst_names[iSyst],iS,iC),pt_bin_limits);
        TH1F* partial_syst_tpc_tmp = (TH1F*) TPC_systematics_file.Get((to_string(iC) + "/" + kNames[iS] + "/" + syst_names[iSyst] + "_tpc").data());
        vec_syst_tpc[iS][iC][iSyst] = (TH1F*) partial_syst_tpc_tmp->Rebin(n_pt_bins,Form("%s_tpc_%d_%d",syst_names[iSyst],iS,iC),pt_bin_limits);
      }
      // Compute the minimum of the systematic uncertainty related to the material budget
      float min_material_tof = 1.,min_material_tpc = 1.;
      float tof_tmp = 0., tpc_tmp = 0.;
      for (int iB = 1; iB <= n_pt_bins; ++iB) {
        tof_tmp = vec_syst_tof[iS][iC][matsyst]->GetBinContent(iB);
        tpc_tmp = vec_syst_tpc[iS][iC][matsyst]->GetBinContent(iB);
        if(tof_tmp>1e-8){
          if(tof_tmp<min_material_tof) min_material_tof = tof_tmp;
        }
        if(tpc_tmp>1e-8){
          if(tpc_tmp<min_material_tpc) min_material_tpc = tpc_tmp;
        }
      }
      comp[iS][iC] = new TH1F(Form("comp_%d_%d",iS,iC),";#it{p}_{T} (GeV/#it{c}); #frac{N_{TPC} - N_{TOF}}{#sigma}",3,0.8,1.1);
      plotting::SetHistStyle(comp[iS][iC],plotting::kSpectraColors[iC]);
      for (int iB = 1; iB <= n_pt_bins; ++iB) {
        if(ptAxis->GetBinCenter(iB)<0.8){
          stat_tof[iS][iC]->SetBinContent(iB,0.);
          stat_tof[iS][iC]->SetBinError(iB,0.);
          syst_tof[iS][iC]->SetBinContent(iB,0.);
          syst_tof[iS][iC]->SetBinError(iB,0.);
          //
          syst_tpc[iS][iC]->SetBinContent(iB,stat_tpc[iS][iC]->GetBinContent(iB));
          syst_tpc[iS][iC]->SetBinError(iB,totsyst_tpc->GetBinContent(iB) * stat_tpc[iS][iC]->GetBinContent(iB));
          //
          stat_all[iS][iC]->SetBinContent(iB,stat_tpc[iS][iC]->GetBinContent(iB));
          stat_all[iS][iC]->SetBinError(iB,stat_tpc[iS][iC]->GetBinError(iB));
          syst_all[iS][iC]->SetBinContent(iB,syst_tpc[iS][iC]->GetBinContent(iB));
          syst_all[iS][iC]->SetBinError(iB,syst_tpc[iS][iC]->GetBinError(iB));
          //
          all[iS][iC]->SetBinContent(iB,stat_tpc[iS][iC]->GetBinContent(iB));
          all[iS][iC]->SetBinError(iB,TMath::Sqrt(stat_tpc[iS][iC]->GetBinError(iB)*stat_tpc[iS][iC]->GetBinError(iB) +                      syst_tpc[iS][iC]->GetBinError(iB)*syst_tpc[iS][iC]->GetBinError(iB)));//-kAbsSyst[iS]*kAbsSyst[iS]*stat_tpc[iS][iC]->GetBinContent(iB)*stat_tpc[iS][iC]->GetBinContent(iB)-min_material_tpc*min_material_tpc*stat_tpc[iS][iC]->GetBinContent(iB)*stat_tpc[iS][iC]->GetBinContent(iB)));
        }
        else if(ptAxis->GetBinCenter(iB)>0.8 && ptAxis->GetBinCenter(iB)<1.1){
          syst_tpc[iS][iC]->SetBinContent(iB,stat_tpc[iS][iC]->GetBinContent(iB));
          syst_tpc[iS][iC]->SetBinError(iB,totsyst_tpc->GetBinContent(iB) * stat_tpc[iS][iC]->GetBinContent(iB));
          syst_tof[iS][iC]->SetBinContent(iB,stat_tof[iS][iC]->GetBinContent(iB));
          syst_tof[iS][iC]->SetBinError(iB,totsyst->GetBinContent(iB) * stat_tof[iS][iC]->GetBinContent(iB));
          float z = stat_tpc[iS][iC]->GetBinContent(iB) - stat_tof[iS][iC]->GetBinContent(iB);
          float err2 = (stat_tof[iS][iC]->GetBinError(iB)-stat_tpc[iS][iC]->GetBinError(iB))*
                       (stat_tof[iS][iC]->GetBinError(iB)-stat_tpc[iS][iC]->GetBinError(iB)) +

                      (vec_syst_tof[iS][iC][cutsyst]->GetBinContent(iB)*stat_tof[iS][iC]->GetBinContent(iB) - vec_syst_tpc[iS][iC][cutsyst]->GetBinContent(iB)*stat_tpc[iS][iC]->GetBinContent(iB))*
                      (vec_syst_tof[iS][iC][cutsyst]->GetBinContent(iB)*stat_tof[iS][iC]->GetBinContent(iB) - vec_syst_tpc[iS][iC][cutsyst]->GetBinContent(iB)*stat_tpc[iS][iC]->GetBinContent(iB)) +

                      (vec_syst_tof[iS][iC][abssyst]->GetBinContent(iB)*stat_tof[iS][iC]->GetBinContent(iB) - vec_syst_tpc[iS][iC][abssyst]->GetBinContent(iB)*stat_tpc[iS][iC]->GetBinContent(iB))*
                      (vec_syst_tof[iS][iC][abssyst]->GetBinContent(iB)*stat_tof[iS][iC]->GetBinContent(iB) - vec_syst_tpc[iS][iC][abssyst]->GetBinContent(iB)*stat_tpc[iS][iC]->GetBinContent(iB)) +

                      (vec_syst_tof[iS][iC][matsyst]->GetBinContent(iB)*stat_tof[iS][iC]->GetBinContent(iB) - vec_syst_tpc[iS][iC][matsyst]->GetBinContent(iB)*stat_tpc[iS][iC]->GetBinContent(iB))*
                      (vec_syst_tof[iS][iC][matsyst]->GetBinContent(iB)*stat_tof[iS][iC]->GetBinContent(iB) - vec_syst_tpc[iS][iC][matsyst]->GetBinContent(iB)*stat_tpc[iS][iC]->GetBinContent(iB)) +

                      vec_syst_tof[iS][iC][countsyst]->GetBinContent(iB)*stat_tof[iS][iC]->GetBinContent(iB)*
                      vec_syst_tof[iS][iC][countsyst]->GetBinContent(iB)*stat_tof[iS][iC]->GetBinContent(iB) +

                      vec_syst_tpc[iS][iC][countsyst]->GetBinContent(iB)*stat_tpc[iS][iC]->GetBinContent(iB)*
                      vec_syst_tpc[iS][iC][countsyst]->GetBinContent(iB)*stat_tpc[iS][iC]->GetBinContent(iB) +

                      vec_syst_tof[iS][iC][shiftsyst]->GetBinContent(iB)*stat_tof[iS][iC]->GetBinContent(iB)*
                      vec_syst_tof[iS][iC][shiftsyst]->GetBinContent(iB)*stat_tof[iS][iC]->GetBinContent(iB) +

                      vec_syst_tpc[iS][iC][shiftsyst]->GetBinContent(iB)*stat_tpc[iS][iC]->GetBinContent(iB)*
                      vec_syst_tpc[iS][iC][shiftsyst]->GetBinContent(iB)*stat_tpc[iS][iC]->GetBinContent(iB);
          z /= TMath::Sqrt(err2);
          comp[iS][iC]->SetBinError(comp[iS][iC]->FindBin(ptAxis->GetBinCenter(iB)),0.);
          comp[iS][iC]->SetBinContent(comp[iS][iC]->FindBin(ptAxis->GetBinCenter(iB)),z);
          if(ptAxis->GetBinCenter(iB)>0.8 && ptAxis->GetBinCenter(iB)<1.){
            stat_all[iS][iC]->SetBinContent(iB,stat_tpc[iS][iC]->GetBinContent(iB));
            stat_all[iS][iC]->SetBinError(iB,stat_tpc[iS][iC]->GetBinError(iB));
            syst_all[iS][iC]->SetBinContent(iB,syst_tpc[iS][iC]->GetBinContent(iB));
            syst_all[iS][iC]->SetBinError(iB,syst_tpc[iS][iC]->GetBinError(iB));
            all[iS][iC]->SetBinContent(iB,stat_tpc[iS][iC]->GetBinContent(iB));
            all[iS][iC]->SetBinError(iB,TMath::Sqrt(stat_tpc[iS][iC]->GetBinError(iB)*stat_tpc[iS][iC]->GetBinError(iB) +                      syst_tpc[iS][iC]->GetBinError(iB)*syst_tpc[iS][iC]->GetBinError(iB)));
          }
          else{
            stat_all[iS][iC]->SetBinContent(iB,stat_tof[iS][iC]->GetBinContent(iB));
            stat_all[iS][iC]->SetBinError(iB,stat_tof[iS][iC]->GetBinError(iB));
            syst_all[iS][iC]->SetBinContent(iB,syst_tof[iS][iC]->GetBinContent(iB));
            syst_all[iS][iC]->SetBinError(iB,syst_tof[iS][iC]->GetBinError(iB));
            //
            all[iS][iC]->SetBinContent(iB,stat_tof[iS][iC]->GetBinContent(iB));
            all[iS][iC]->SetBinError(iB,TMath::Sqrt(stat_tof[iS][iC]->GetBinError(iB)*stat_tof[iS][iC]->GetBinError(iB) +            syst_tof[iS][iC]->GetBinError(iB)*syst_tof[iS][iC]->GetBinError(iB)));
          }
        }
        else{
          stat_tpc[iS][iC]->SetBinContent(iB,0.);
          stat_tpc[iS][iC]->SetBinError(iB,0.);
          syst_tpc[iS][iC]->SetBinContent(iB,0.);
          syst_tpc[iS][iC]->SetBinError(iB,0.);
          //
          syst_tof[iS][iC]->SetBinContent(iB,stat_tof[iS][iC]->GetBinContent(iB));
          syst_tof[iS][iC]->SetBinError(iB,totsyst->GetBinContent(iB) * stat_tof[iS][iC]->GetBinContent(iB));
          //
          stat_all[iS][iC]->SetBinContent(iB,stat_tof[iS][iC]->GetBinContent(iB));
          stat_all[iS][iC]->SetBinError(iB,stat_tof[iS][iC]->GetBinError(iB));
          syst_all[iS][iC]->SetBinContent(iB,syst_tof[iS][iC]->GetBinContent(iB));
          syst_all[iS][iC]->SetBinError(iB,syst_tof[iS][iC]->GetBinError(iB));
          //
          all[iS][iC]->SetBinContent(iB,stat_tof[iS][iC]->GetBinContent(iB));
          all[iS][iC]->SetBinError(iB,TMath::Sqrt(stat_tof[iS][iC]->GetBinError(iB)*stat_tof[iS][iC]->GetBinError(iB) +            syst_tof[iS][iC]->GetBinError(iB)*syst_tof[iS][iC]->GetBinError(iB)));//-kAbsSyst[iS]*kAbsSyst[iS]*stat_tof[iS][iC]->GetBinContent(iB)*stat_tof[iS][iC]->GetBinContent(iB)-min_material_tof*min_material_tof*stat_tof[iS][iC]->GetBinContent(iB)*stat_tof[iS][iC]->GetBinContent(iB)));
        }
      }
      plotting::SetHistStyle(stat_tof[iS][iC],plotting::kSpectraColors[iC]);
      plotting::SetHistStyle(syst_tof[iS][iC],plotting::kSpectraColors[iC]);
      stat_tof[iS][iC]->Write("stat_tof");
      stat_tof[iS][iC]->Scale(kScaleFactor[iC]);
      syst_tof[iS][iC]->Write("syst_tof");
      syst_tof[iS][iC]->Scale(kScaleFactor[iC]);
      plotting::SetHistStyle(stat_tpc[iS][iC],plotting::kSpectraColors[iC],25);
      plotting::SetHistStyle(syst_tpc[iS][iC],plotting::kSpectraColors[iC],25);
      stat_tpc[iS][iC]->Write("stat_tpc");
      stat_tpc[iS][iC]->Scale(kScaleFactor[iC]);
      syst_tpc[iS][iC]->Write("syst_tpc");
      syst_tpc[iS][iC]->Scale(kScaleFactor[iC]);
      stat_all[iS][iC]->Write("stat_all");
      syst_all[iS][iC]->Write("syst_all");
      stat_all[iS][iC]->Scale(kScaleFactor[iC]);
      syst_all[iS][iC]->Scale(kScaleFactor[iC]);
      all[iS][iC]->Scale(kScaleFactor[iC]);
      //all[iS][iC]->Fit(Form("Levi-Tsallis_%d_%d",iS,iC),"IQ");
      plotting::SetHistStyle(stat_all[iS][iC],plotting::kSpectraColors[iC]);
      plotting::SetHistStyle(syst_all[iS][iC],plotting::kSpectraColors[iC]);
    }


    TCanvas spectraComp("spectraComp","spectraComp",3200,2400);
    spectraComp.DrawFrame(
        0.4,
        1e-6,
        3.8,
        0.5,
        ";#it{p}_{T} (GeV/#it{c});#frac{1}{#it{N}_{ev}} #frac{d^{2}#it{N}}{d#it{p}_{T}d#it{y}} (GeV/#it{c})^{-1}"
        );
    spectraComp.SetLeftMargin(0.1425591);
    spectraComp.SetRightMargin(0.027121);
    spectraComp.SetTopMargin(0.06053269);
    spectraComp.SetBottomMargin(0.1598063);
    TLegend final_leg_comp(0.70,0.18,0.94,0.53);
    final_leg_comp.SetBorderSize(0);
    final_leg_comp.SetHeader(Form("%s, pp #sqrt{s} = 13 TeV",kNames[iS].data()));
    //TLegendEntry *header = (TLegendEntry*)final_leg.GetListOfPrimitives()->First();
    //header->SetTextAlign(22);
    final_leg_comp.SetNColumns(2);
    final_leg_comp.AddEntry((TObject*)nullptr,"TPC","");
    final_leg_comp.AddEntry((TObject*)nullptr,"TPC + TOF","");
    for (int iC = 0; iC < n_centralities -1; ++iC) {
      stat_tof[iS][iC]->Draw("esamex0");
      syst_tof[iS][iC]->Draw("e2same");
      stat_tpc[iS][iC]->Draw("esamex0");
      syst_tpc[iS][iC]->Draw("e2same");
      // stat_all[iS][iC]->Draw("esamex0");
      // syst_all[iS][iC]->Draw("e2same");
      //all[iS][iC]->Draw("esamex0");
      final_leg_comp.AddEntry(syst_tpc[iS][iC],Form("%4.0f - %2.0f %% (#times %d)",kCentLabels[iC][0],kCentLabels[iC][1],kScaleFactor[iC]),"fp");
      final_leg_comp.AddEntry(syst_tof[iS][iC],Form("%4.0f - %2.0f %% (#times %d)",kCentLabels[iC][0],kCentLabels[iC][1],kScaleFactor[iC]),"fp");
    }
    final_leg_comp.Draw();
    spectraComp.SetLogy();
    s_dir->cd();
    spectraComp.Write();

    TCanvas spectra("spectra","spectra",3200,2400);
    spectra.DrawFrame(
        0.4,
        1e-6,
        3.8,
        0.5,
        ";#it{p}_{T} (GeV/#it{c});#frac{1}{#it{N}_{ev}} #frac{d^{2}#it{N}}{d#it{p}_{T}d#it{y}} (GeV/#it{c})^{-1}"
        );
    spectra.SetLeftMargin(0.1425591);
    spectra.SetRightMargin(0.027121);
    spectra.SetTopMargin(0.06053269);
    spectra.SetBottomMargin(0.1598063);
    TLegend final_leg(0.70,0.18,0.94,0.53);
    final_leg.SetBorderSize(0);
    final_leg.SetHeader(Form("%s, pp #sqrt{s} = 13 TeV",kNames[iS].data()));
    //TLegendEntry *header = (TLegendEntry*)final_leg.GetListOfPrimitives()->First();
    //header->SetTextAlign(22);
    final_leg.SetNColumns(2);
    //final_leg.AddEntry((TObject*)nullptr,"TPC","");
    //final_leg.AddEntry((TObject*)nullptr,"TPC + TOF","");
    for (int iC = 0; iC < n_centralities -1; ++iC) {
      // stat_tof[iS][iC]->Draw("esamex0");
      // syst_tof[iS][iC]->Draw("e2same");
      // stat_tpc[iS][iC]->Draw("esamex0");
      // syst_tpc[iS][iC]->Draw("e2same");
      stat_all[iS][iC]->Draw("esamex0");
      syst_all[iS][iC]->Draw("e2same");
      //all[iS][iC]->Draw("esamex0");
      final_leg.AddEntry(syst_all[iS][iC],Form("%4.0f - %2.0f %% (#times %d)",kCentLabels[iC][0],kCentLabels[iC][1],kScaleFactor[iC]),"fp");
      //final_leg.AddEntry(syst_tof[iS][iC],Form("%4.0f - %2.0f %% (#times %d)",kCentLabels[iC][0],kCentLabels[iC][1],1<<(kCentLength-iC-1)),"fp");
    }
    TFile bwfile(Form("%s%sfits.root",kBaseOutputDir.data(),kPrefix[iS]),"read");
    if (bwfile.IsOpen()) {
      TF1* bw = nullptr;
      TH1F* scaledbw = nullptr;
      for (int iC = 0; iC < n_centralities -1 ; ++iC) {
        bw = (TF1*)bwfile.Get(Form("%s/%s%i",fitFunctionName,fitFunctionName,iC));
        if (!bw) continue;
        scaledbw = new TH1F(Form("scaledbw%i",iC),"",1000,0.5,1.05*kCentPtLimits[iC]);
        scaledbw->Add(bw);
        scaledbw->Scale(kScaleFactor[iC]);
        scaledbw->SetLineStyle(kDashed);
        scaledbw->SetLineColor(kBlack);
        spectra.cd();
        scaledbw->Draw("lsame");
        for(int i=1; i<=n_pt_bins; i++){
          if(all[iS][iC]->GetBinCenter(i)>kCentPtLimits[iC]) continue;
          float data_val = all[iS][iC]->GetBinContent(i);
          float fit_val = scaledbw->GetBinContent(scaledbw->FindBin(all[iS][iC]->GetBinCenter(i)));
          float data_err = all[iS][iC]->GetBinError(i);
          float pull = (data_val-fit_val)/data_err;
          printf("iC: %d, iB: %d, data_val: %f, fit_val: %f, data_err: %f, pull: %f\n", iC, i, data_val, fit_val, data_err, pull);
          if( all[iS][iC]->GetBinCenter(i)<kCentPtLimits[iC]){
            pull_fit_data[iS][iC]->SetBinContent(i,pull);
            ratio_fit_data[iS][iC]->SetBinContent(i,data_val/fit_val);
            ratio_fit_data[iS][iC]->SetBinError(i,data_err/fit_val);
          }
        }
      }
      final_leg.AddEntry(scaledbw,"Individual fit","l");
    }
    final_leg.Draw();
    spectra.SetLogy();
    s_dir->cd();
    spectra.Write();
    if (kPrintFigures) {
      spectra.SaveAs((kFiguresFolder + "spectraTOF" + kLetter[iS] + ".eps").data());
      spectra.SaveAs((kMacrosFolder + "spectraTOF" + kLetter[iS] + ".C").data());
    }
    TCanvas cComp("cComp","cComp",3200,2400);
    cComp.DrawFrame(0.75,-8.,1.15,8.,";#it{p}_{T} (GeV/#it{c});#frac{N_{TPC} - N_{TOF}}{#sigma}");
    TLegend comp_leg(0.80,0.17,0.94,0.45);
    comp_leg.SetBorderSize(0);
    comp_leg.SetHeader(Form("%s, pp #sqrt{s} = 13 TeV",kNames[iS].data()));
    for (int iC = 0; iC < n_centralities -1; ++iC) {
      comp[iS][iC]->Draw("psame");
      comp_leg.AddEntry(comp[iS][iC],Form("%4.0f - %2.0f %%",kCentLabels[iC][0],kCentLabels[iC][1]),"p");
    }
    for(int i=-3; i<=3; i++){
      TLine *line = new TLine(0.75,i,1.15,i);
      int absi = (i<0) ? -1*i : i;
      switch(absi){
        case(3):
          line->SetLineColor(kRed);
          break;
        case(2):
          line->SetLineColor(kOrange);
          break;
        case(1):
          line->SetLineColor(kGreen+3);
          break;
        default:
          line->SetLineColor(kBlack);
          break;
      }
      line->SetLineStyle(2);
      line->Draw();
    }
    comp_leg.Draw();
    cComp.Write();

    TCanvas* cPullFit = new TCanvas(Form("cPullFit_%c",kLetter[iS]),Form("cPullFit_%c",kLetter[iS]));
    cPullFit->Divide(3,3);
    TCanvas* cRatioFit = new TCanvas(Form("cRatioFit_%c",kLetter[iS]),Form("cRatioFit_%c",kLetter[iS]));
    cRatioFit->Divide(3,3);
    for(int iC=0; iC<kCentLength-1; iC++){
      cPullFit->cd(iC+1);
      plotting::SetHistStyle(pull_fit_data[iS][iC],plotting::kSpectraColors[iC]);
      pull_fit_data[iS][iC]->GetXaxis()->SetRangeUser(0.6,kCentPtLimits[iC]);
      pull_fit_data[iS][iC]->Draw("p");
      for(int i=-2; i<=2; i++){
        TLine *line = new TLine(0.6,i,kCentPtLimits[iC],i);
        int absi = (i<0) ? -1*i : i;
        switch(absi){
          case(2):
            line->SetLineColor(kRed);
            break;
          case(1):
            line->SetLineColor(kGreen+3);
            break;
          default:
            line->SetLineColor(kBlack);
            break;
        }
        line->SetLineStyle(2);
        line->Draw();
      }
      // TF1 *g = new TF1(Form("g_%c_%d",kLetter[iS],iC),"gaus",-5,5);
      // g->SetParName(0,"Norm");
      // g->SetParName(1,"#mu");
      // g->SetParameter(1,0);
      // g->SetParName(2,"#sigma");
      // g->SetParameter(2,1);
      // pull_fit_data[iS][iC]->Fit(Form("g_%c_%d",kLetter[iS],iC));
      cRatioFit->cd(iC+1);
      plotting::SetHistStyle(ratio_fit_data[iS][iC],plotting::kSpectraColors[iC]);
      ratio_fit_data[iS][iC]->GetXaxis()->SetRangeUser(0.6,kCentPtLimits[iC]);
      ratio_fit_data[iS][iC]->Draw("p");
      TLine *line_one = new TLine(0.6,1.,kCentPtLimits[iC],1.);
      line_one->SetLineColor(kBlack);
      line_one->SetLineStyle(2);
      line_one->Draw();
    }
    cPullFit->Write();
    cRatioFit->Write();
  }

  TDirectory* r_dir = final_file.mkdir("ratio");
  TF1* funcpol[n_centralities];

  r_dir->cd();
  TPad* pad[9] = {nullptr};
  TCanvas ratio("ratio","ratio",3200,3200);
  plotting::CanvasPartition(&ratio,pad,3,3);

  float mean_ratio[n_centralities] = {0.};
  float mean_ratio_stat[n_centralities] = {0.};
  float mean_ratio_syst[n_centralities] = {0.};
  float mean_ratio_tot[n_centralities] = {0.};

  for (int iC = 0; iC < n_centralities -1; ++iC) {
    r_dir->mkdir(to_string(iC).data())->cd();
    stat_tof[1][iC]->Divide(stat_tof[0][iC]);
    syst_tof[1][iC]->Divide(syst_tof[0][iC]);
    all[1][iC]->Divide(all[0][iC]);
    stat_all[1][iC]->Divide(stat_all[0][iC]);
    syst_all[1][iC]->Divide(syst_all[0][iC]);
    funcpol[iC] = new TF1(Form("ratiopol_%d",iC),"pol0",0.6,kCentPtLimits[iC]);
    int nx = iC/3;
    int ny = iC%3;
    pad[iC]->cd();
    //stat_all[1][iC]->Fit(Form("ratiopol_%d",iC),"Q");
    for(int iB=1; iB<=n_pt_bins; iB++){
      if(stat_tof[1][iC]->GetBinCenter(iB)<kTOFminPt){
        stat_tof[1][iC]->SetBinContent(iB, 0.);
        stat_tof[1][iC]->SetBinError(iB, 0.);
        syst_tof[1][iC]->SetBinContent(iB, 0.);
        syst_tof[1][iC]->SetBinError(iB, 0.);
      }
    }
    stat_tof[1][iC]->Write("ratio_tof");
    syst_tof[1][iC]->Write("ratio_tof");

    stat_tpc[1][iC]->Divide(stat_tpc[0][iC]);
    syst_tpc[1][iC]->Divide(syst_tpc[0][iC]);
    for(int iB=1; iB<=n_pt_bins; iB++){
      if(stat_tpc[1][iC]->GetBinCenter(iB)>kTPCmaxPt){
        stat_tpc[1][iC]->SetBinContent(iB, 0.);
        stat_tpc[1][iC]->SetBinError(iB, 0.);
        syst_tpc[1][iC]->SetBinContent(iB, 0.);
        syst_tpc[1][iC]->SetBinError(iB, 0.);
      }
    }

    // mean value of the ratio
    for(int iB=1; iB<n_pt_bins; iB++){
      if(stat_tof[1][iC]->GetBinCenter(iB)>kCentPtLimits[iC]) continue;
      float tmp_val=0., tmp_stat=0., tmp_syst=0.;
      if(stat_tof[1][iC]->GetBinCenter(iB)<1.0){
        tmp_val = stat_tpc[1][iC]->GetBinContent(iB);
        tmp_stat = stat_tpc[1][iC]->GetBinError(iB);
        tmp_syst = syst_tpc[1][iC]->GetBinError(iB);
      }
      else{
        tmp_val = stat_tof[1][iC]->GetBinContent(iB);
        tmp_stat = stat_tof[1][iC]->GetBinError(iB);
        tmp_syst = syst_tof[1][iC]->GetBinError(iB);
      }
      //printf("iC %d ; tmp_val %f ; tmp_stat %f ; tmp_syst %f ; \n", iC, tmp_val, tmp_stat, tmp_syst);
      mean_ratio_stat[iC] += 1/Sq(tmp_stat);
      mean_ratio_syst[iC] += 1/Sq(tmp_syst);
      mean_ratio_tot[iC] += 1/(Sq(tmp_stat)+Sq(tmp_syst));
      mean_ratio[iC] += tmp_val/(Sq(tmp_stat)+Sq(tmp_syst));
    }
    mean_ratio[iC] /= mean_ratio_tot[iC];
    mean_ratio_stat[iC] = TMath::Sqrt(1/mean_ratio_stat[iC]);
    mean_ratio_syst[iC] = TMath::Sqrt(1/mean_ratio_syst[iC]);
    std::cout << mean_ratio[iC] << " $\\pm$ " << mean_ratio_stat[iC] << " $\\pm$ " << mean_ratio_syst[iC] << "\\\\" << std::endl;

    //Fill the unique canvas with all the manti-matter/matter ratios
    pad[iC]->cd();
    double XaxisEdge = 0.;
    switch (nx) {
      case 0:
        XaxisEdge=3.9;
        break;
      case 1:
        XaxisEdge=3.2;
        break;
      default:
        XaxisEdge=2.4;
        break;
    }
    gPad->DrawFrame(
        0.35,
        0.1,
        XaxisEdge,
        2.1,
    ";#it{p}_{T} (GeV/#it{c});#bar{d}/d"
    );
    // stat_tof[1][iC]->Draw("esamex0");
    // syst_tof[1][iC]->Draw("e2same");
    // stat_tpc[1][iC]->Draw("esamex0");
    // syst_tpc[1][iC]->Draw("e2same");
    stat_all[1][iC]->Draw("esamex0");
    syst_all[1][iC]->Draw("e2same");
    //funcpol[iC]->Draw("same");
    TLine *line_one = new TLine(0.35,1.,XaxisEdge,1.);
    line_one->SetLineColor(kBlack);
    line_one->SetLineStyle(2);
    line_one->Draw();
    // TLegend* ratio_leg_one = new TLegend(0.49,1.51,1.61,2.01,Form("%4.0f - %2.0f %%",kCentLabels[iC][0],kCentLabels[iC][1]),"br");
    // ratio_leg_one->SetFillStyle(0);
    // ratio_leg_one->SetBorderSize(0);
    // ratio_leg_one->AddEntry(syst_tof[0][iC],"TPC + TOF","p");
    // ratio_leg_one->AddEntry(syst_tpc[0][iC],"TPC","p");
    // ratio_leg_one->AddEntry((TObject*)0, Form("p0: %.2f #pm %.2f", funcpol[iC]->GetParameter(0),funcpol[iC]->GetParError(0)), "");
    // ratio_leg_one->AddEntry((TObject*)0, Form("#chi^{2}/NDF: %.2f / %d",funcpol[iC]->GetChisquare(),funcpol[iC]->GetNDF()),"");
    //
    // ratio_leg_one->Draw();

    TPaveText* pav = new TPaveText(0.55,1.51,0.55+(XaxisEdge-0.5)*0.2,2.01,"bl");
    pav->SetBorderSize(0);
    pav->SetFillColor(0);
    pav->AddText(Form("%4.0f - %2.0f %%",kCentLabels[iC][0],kCentLabels[iC][1]));
    pav->Draw();
  }
  r_dir->cd();
  if (kPrintFigures) ratio.SaveAs((kFiguresFolder + "ratio.eps").data());
  ratio.Write();

}
