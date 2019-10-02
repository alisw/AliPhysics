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
#include <TLatex.h>

const char* kPrefix[2] = {"","anti"};

const char* kParticleNames[2] = {"Deuterons", "Antideuterons"};

void Final_check(const char* fitFunctionName = "LevyTsallis") {
  TFile spectra_file(kSpectraOutput.data());
  TFile TOF_systematics_file(kSystematicsOutput.data());
  TFile TPC_systematics_file(kSystematicsOutputTPC.data());
  TFile final_file(Form("%sfinal_check.root",kBaseOutputDir.data()),"recreate");

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
      TDirectory *c_dir = s_dir->mkdir(std::to_string(iC).data());
      c_dir->cd();
      pull_fit_data[iS][iC] = new TH1F(Form("pull_fit_data_%c_%d",kLetter[iS],iC),Form("%2.0f - %2.0f %%;#it{p}_{T} (GeV/#it{c});#frac{fit - data}{#sigma_{data}}",kCentLabels[iC][0],kCentLabels[iC][1]),kNPtBins,kPtBins);
      pull_fit_data[iS][iC]->GetYaxis()->SetRangeUser(-3.,3.);
      ratio_fit_data[iS][iC] = new TH1F(Form("ratio_fit_data_%c_%d",kLetter[iS],iC),Form("%2.0f - %2.0f %%;#it{p}_{T} (GeV/#it{c});data/fit",kCentLabels[iC][0],kCentLabels[iC][1]),kNPtBins,kPtBins);
      ratio_fit_data[iS][iC]->GetYaxis()->SetRangeUser(0.,2.);
      std::string tot_syst_path = kNames[iS] + "/" + std::to_string(iC) + "/totsyst_" + std::to_string(iC);
      TH1F* totsyst = (TH1F*)TOF_systematics_file.Get(tot_syst_path.data());
      Requires(totsyst, tot_syst_path.data());
      std::string tpc_syst_path = kNames[iS] + "/" + std::to_string(iC) + "/totsyst_" + std::to_string(iC);
      TH1F* totsyst_tpc = (TH1F*)TPC_systematics_file.Get(tpc_syst_path.data());
      Requires(totsyst_tpc, tpc_syst_path.data());
      string tof_basepath = kFilterListNames + "/" + kNames[iS] + "/" + std::to_string(iC) + "/TOF/TOFspectra" + kLetter[iS] + std::to_string(iC);
      TH1F* spectra_tof_tmp  = (TH1F*)spectra_file.Get(tof_basepath.data());
      Requires(spectra_tof_tmp,tof_basepath.data());
      spectra_tof_tmp->Write("spectra_tof_tmp");
      stat_tof[iS][iC] = (TH1F*)spectra_tof_tmp->Rebin(kNPtBins,Form("spectra_tof_%d_%c",iC,kLetter[iS]),kPtBins);
      auto ptAxis = stat_tof[iS][iC]->GetXaxis();
      string tpc_basepath = kFilterListNames + "/" + kNames[iS] + "/" + std::to_string(iC) + "/TPC/TPCspectra" + kLetter[iS] + std::to_string(iC);
      TH1F* spectra_tpc_tmp  = (TH1F*)spectra_file.Get(tpc_basepath.data());
      Requires(spectra_tpc_tmp,tpc_basepath.data());
      spectra_tpc_tmp->Write("spectra_tpc_tmp");
      stat_tpc[iS][iC] = (TH1F*)spectra_tpc_tmp->Rebin(kNPtBins,Form("spectra_tpc_%d_%c",iC,kLetter[iS]),kPtBins);
      syst_tof[iS][iC]  = (TH1F*)stat_tpc[iS][iC]->Clone(("syst_tof" + std::to_string(iC)).data());
      syst_tof[iS][iC]->Reset();
      syst_tpc[iS][iC]  = (TH1F*)stat_tpc[iS][iC]->Clone(("syst_tpc" + std::to_string(iC)).data());
      syst_tpc[iS][iC]->Reset();
      stat_all[iS][iC] = (TH1F*)stat_tpc[iS][iC]->Clone(("stat_all" + std::to_string(iC)).data());
      stat_all[iS][iC]->Reset();
      syst_all[iS][iC] = (TH1F*)stat_tpc[iS][iC]->Clone(("syst_all" + std::to_string(iC)).data());
      syst_all[iS][iC]->Reset();
      all[iS][iC] = (TH1F*)stat_tpc[iS][iC]->Clone(("all" + std::to_string(iC)).data());
      all[iS][iC]->Reset();
      for(int iSyst=0; iSyst<5; iSyst++){
        string cut_syst_path = kNames[iS] + "/" + std::to_string(iC) + "/" + syst_names[iSyst] + "_" + std::to_string(iC);
        vec_syst_tof[iS][iC][iSyst] = (TH1F*) TOF_systematics_file.Get(cut_syst_path.data());
        Requires(vec_syst_tof[iS][iC][iSyst], ("TOF: " +std::to_string(iC) + "/" + kNames[iS] + "/" + syst_names[iSyst] + "_" + std::to_string(iC)).data());
        vec_syst_tpc[iS][iC][iSyst] = (TH1F*) TPC_systematics_file.Get(cut_syst_path.data());
        Requires(vec_syst_tpc[iS][iC][iSyst], ("TPC: " + std::to_string(iC) + "/" + kNames[iS] + "/" + syst_names[iSyst] + "_" + std::to_string(iC)).data());
      }
      // Compute the minimum of the systematic uncertainty related to the material budget
      float min_material_tof = 1.,min_material_tpc = 1.;
      float tof_tmp = 0., tpc_tmp = 0.;
      for (int iB = 1; iB <= kNPtBins; ++iB) {
        tof_tmp = vec_syst_tof[iS][iC][matsyst]->GetBinContent(iB);
        tpc_tmp = vec_syst_tpc[iS][iC][matsyst]->GetBinContent(iB);
        if(tof_tmp>1e-8){
          if(tof_tmp<min_material_tof) min_material_tof = tof_tmp;
        }
        if(tpc_tmp>1e-8){
          if(tpc_tmp<min_material_tpc) min_material_tpc = tpc_tmp;
        }
      }
      comp[iS][iC] = new TH1F(Form("comp_%d_%d",iS,iC),";#it{p}_{T} (GeV/#it{c}); #frac{N_{TPC} - N_{TOF}}{#sigma}",3,0.9,1.2);
      plotting::SetHistStyle(comp[iS][iC],plotting::kSpectraColors[iC]);
      for (int iB = 1; iB <= kNPtBins; ++iB) {
        if(ptAxis->GetBinCenter(iB)<0.9){
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
        else if(ptAxis->GetBinCenter(iB)>0.9 && ptAxis->GetBinCenter(iB)<1.2){
          syst_tpc[iS][iC]->SetBinContent(iB,stat_tpc[iS][iC]->GetBinContent(iB));
          syst_tpc[iS][iC]->SetBinError(iB,totsyst_tpc->GetBinContent(iB) * stat_tpc[iS][iC]->GetBinContent(iB));
          syst_tof[iS][iC]->SetBinContent(iB,stat_tof[iS][iC]->GetBinContent(iB));
          syst_tof[iS][iC]->SetBinError(iB,totsyst->GetBinContent(iB) * stat_tof[iS][iC]->GetBinContent(iB));
          float z = stat_tpc[iS][iC]->GetBinContent(iB) - stat_tof[iS][iC]->GetBinContent(iB);
          // float err2 = (stat_tof[iS][iC]->GetBinError(iB)-stat_tpc[iS][iC]->GetBinError(iB))*
          //              (stat_tof[iS][iC]->GetBinError(iB)-stat_tpc[iS][iC]->GetBinError(iB)) +

          //             (vec_syst_tof[iS][iC][cutsyst]->GetBinContent(iB)*stat_tof[iS][iC]->GetBinContent(iB) - vec_syst_tpc[iS][iC][cutsyst]->GetBinContent(iB)*stat_tpc[iS][iC]->GetBinContent(iB))*
          //             (vec_syst_tof[iS][iC][cutsyst]->GetBinContent(iB)*stat_tof[iS][iC]->GetBinContent(iB) - vec_syst_tpc[iS][iC][cutsyst]->GetBinContent(iB)*stat_tpc[iS][iC]->GetBinContent(iB)) +

          //             (vec_syst_tof[iS][iC][abssyst]->GetBinContent(iB)*stat_tof[iS][iC]->GetBinContent(iB)* vec_syst_tof[iS][iC][abssyst]->GetBinContent(iB)*stat_tof[iS][iC]->GetBinContent(iB))+

          //             // (vec_syst_tof[iS][iC][abssyst]->GetBinContent(iB)*stat_tof[iS][iC]->GetBinContent(iB) - vec_syst_tpc[iS][iC][abssyst]->GetBinContent(iB)*stat_tpc[iS][iC]->GetBinContent(iB))*
          //             // (vec_syst_tof[iS][iC][abssyst]->GetBinContent(iB)*stat_tof[iS][iC]->GetBinContent(iB) - vec_syst_tpc[iS][iC][abssyst]->GetBinContent(iB)*stat_tpc[iS][iC]->GetBinContent(iB)) +

          //             (vec_syst_tof[iS][iC][matsyst]->GetBinContent(iB)*stat_tof[iS][iC]->GetBinContent(iB) - vec_syst_tpc[iS][iC][matsyst]->GetBinContent(iB)*stat_tpc[iS][iC]->GetBinContent(iB))*
          //             (vec_syst_tof[iS][iC][matsyst]->GetBinContent(iB)*stat_tof[iS][iC]->GetBinContent(iB) - vec_syst_tpc[iS][iC][matsyst]->GetBinContent(iB)*stat_tpc[iS][iC]->GetBinContent(iB)) +

          //             vec_syst_tof[iS][iC][countsyst]->GetBinContent(iB)*stat_tof[iS][iC]->GetBinContent(iB)*
          //             vec_syst_tof[iS][iC][countsyst]->GetBinContent(iB)*stat_tof[iS][iC]->GetBinContent(iB) +

          //             vec_syst_tpc[iS][iC][countsyst]->GetBinContent(iB)*stat_tpc[iS][iC]->GetBinContent(iB)*
          //             vec_syst_tpc[iS][iC][countsyst]->GetBinContent(iB)*stat_tpc[iS][iC]->GetBinContent(iB) +

          //             vec_syst_tof[iS][iC][shiftsyst]->GetBinContent(iB)*stat_tof[iS][iC]->GetBinContent(iB)*
          //             vec_syst_tof[iS][iC][shiftsyst]->GetBinContent(iB)*stat_tof[iS][iC]->GetBinContent(iB) +

          //             vec_syst_tpc[iS][iC][shiftsyst]->GetBinContent(iB)*stat_tpc[iS][iC]->GetBinContent(iB)*
          //             vec_syst_tpc[iS][iC][shiftsyst]->GetBinContent(iB)*stat_tpc[iS][iC]->GetBinContent(iB);
          float err2 = (stat_tof[iS][iC]->GetBinError(iB)*stat_tof[iS][iC]->GetBinError(iB)) +
                       (stat_tpc[iS][iC]->GetBinError(iB)*stat_tpc[iS][iC]->GetBinError(iB)) +
                       (syst_tof[iS][iC]->GetBinError(iB)*syst_tof[iS][iC]->GetBinError(iB)) +
                       (syst_tpc[iS][iC]->GetBinError(iB)*syst_tpc[iS][iC]->GetBinError(iB));

                    
          z /= TMath::Sqrt(err2);
          comp[iS][iC]->SetBinError(comp[iS][iC]->FindBin(ptAxis->GetBinCenter(iB)),0.);
          comp[iS][iC]->SetBinContent(comp[iS][iC]->FindBin(ptAxis->GetBinCenter(iB)),z);
          if(ptAxis->GetBinCenter(iB)>0.9 && ptAxis->GetBinCenter(iB)<2.){
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
        2e-7,
        5.,
        6,
        ";#it{p}_{T} (GeV/#it{c});#frac{1}{#it{N}_{ev}} #frac{d^{2}#it{N}}{d#it{p}_{T}d#it{y}} (GeV/#it{c})^{-1}"
        );
    spectraComp.SetLeftMargin(0.1425591);
    spectraComp.SetRightMargin(0.027121);
    spectraComp.SetTopMargin(0.06053269);
    spectraComp.SetBottomMargin(0.1598063);
    TLegend final_leg_comp(0.73,0.17,0.94,0.38);
    final_leg_comp.SetBorderSize(0);
    final_leg_comp.SetHeader(Form("%s, pp #sqrt{s} = 13 TeV",kParticleNames[iS]));
    //TLegendEntry *header = (TLegendEntry*)final_leg.GetListOfPrimitives()->First();
    //header->SetTextAlign(22);
    final_leg_comp.SetNColumns(2);
    final_leg_comp.AddEntry((TObject*)nullptr,"TPC","");
    final_leg_comp.AddEntry((TObject*)nullptr,"TPC + TOF","");
    for (int iC = 0; iC < n_centralities; ++iC) {
      stat_tof[iS][iC]->Draw("esamex0");
      syst_tof[iS][iC]->Draw("e2same");
      stat_tpc[iS][iC]->Draw("esamex0");
      syst_tpc[iS][iC]->Draw("e2same");
      // stat_all[iS][iC]->Draw("esamex0");
      // syst_all[iS][iC]->Draw("e2same");
      //all[iS][iC]->Draw("esamex0");
      final_leg_comp.AddEntry(syst_tpc[iS][iC],Form("%4.0f - %2.0f %% (#times 2^{ %d})",kCentLabels[iC][0],kCentLabels[iC][1],kExponent[iC]),"fp");
      final_leg_comp.AddEntry(syst_tof[iS][iC],Form("%4.0f - %2.0f %% (#times 2^{%d})",kCentLabels[iC][0],kCentLabels[iC][1],kExponent[iC]),"fp");
    }
    final_leg_comp.Draw();
    spectraComp.SetLogy();
    s_dir->cd();
    spectraComp.Write();

    TCanvas spectra("spectra","spectra",800,600);
    TH1* hFrame = spectra.DrawFrame(
        0.4,
        2e-7,
        5.,
        6,
        ";#it{p}_{T} (GeV/#it{c});#frac{1}{#it{N}_{ev}} #frac{d^{2}#it{N}}{d#it{p}_{T}d#it{y}} (GeV/#it{c})^{-1}"
        );
    spectra.SetLeftMargin(0.15);
    spectra.SetRightMargin(0.03);
    spectra.SetTopMargin(0.1);
    spectra.SetBottomMargin(0.14);
    hFrame->GetYaxis()->SetTitleOffset(1.3);
    TLatex text;
    text.SetTextFont(63);
    text.SetTextSize(22);
    //text.DrawText(0.5,0.75,"ALICE Preliminary");
    float name_position = (iS==0) ? 3.3 : 3.0;
    text.DrawLatex(name_position,7.5,Form("#bf{%s, pp, #sqrt{#it{s}} = 13 TeV}",kNames[iS].data()));
    text.DrawLatex(3.35461,0.41997,Form("#bf{#color[%d]{#LTd#it{N}_{ch} / d#it{#eta}_{lab}#GT = 26.22}}",plotting::kSpectraColors[0]));
    text.DrawLatex(0.641832,1.02914e-06,Form("#bf{#color[%d]{#LTd#it{N}_{ch} / d#it{#eta}_{lab}#GT = 2.42}}",plotting::kSpectraColors[8]));
    TLegend final_leg(0.710526,0.257391,0.922306,0.58087);
    final_leg.SetBorderSize(0);
    final_leg.SetTextSize(0.027);
    final_leg.SetHeader("V0M Multiplicity\nClasses");
    TLegendEntry *header = (TLegendEntry*)final_leg.GetListOfPrimitives()->First();
    header->SetTextAlign(22);
    final_leg.SetNColumns(2);
    //final_leg.AddEntry((TObject*)nullptr,"TPC","");
    //final_leg.AddEntry((TObject*)nullptr,"TPC + TOF","");
    for (int iC = 0; iC < n_centralities; ++iC) {
      // stat_tof[iS][iC]->Draw("esamex0");
      // syst_tof[iS][iC]->Draw("e2same");
      // stat_tpc[iS][iC]->Draw("esamex0");
      // syst_tpc[iS][iC]->Draw("e2same");
      stat_all[iS][iC]->Draw("esamex0");
      syst_all[iS][iC]->Draw("e2same");
      //all[iS][iC]->Draw("esamex0");
      final_leg.AddEntry(syst_all[iS][iC],Form("%s (#times 2^{ %d})",kRomanLabels[iC],kExponent[iC]),"fp");
      //final_leg.AddEntry(syst_tof[iS][iC],Form("%4.0f - %2.0f %% (#times %d)",kCentLabels[iC][0],kCentLabels[iC][1],1<<(kCentLength-iC-1)),"fp");
    }
    TFile bwfile(Form("%s%sfits.root",kBaseOutputDir.data(),kPrefix[iS]),"read");
    if (bwfile.IsOpen()) {
      TF1* bw = nullptr;
      TH1F* scaledbw = nullptr;
      for (int iC = 0; iC < n_centralities; ++iC) {
        bw = (TF1*)bwfile.Get(Form("%s/%i/%s%i",fitFunctionName,iC,fitFunctionName,iC));
        Requires(bw, Form("%s/%i/%s%i",fitFunctionName,iC,fitFunctionName,iC));
        if (!bw) continue;
        scaledbw = new TH1F(Form("scaledbw%i",iC),"",1000,0.5,1.05*kCentPtLimits[iC]);
        scaledbw->Add(bw);
        scaledbw->Scale(kScaleFactor[iC]);
        scaledbw->SetLineStyle(kDashed);
        scaledbw->SetLineColor(kBlack);
        spectra.cd();
        scaledbw->Draw("lsame");
        for(int i=1; i<=kNPtBins; i++){
          if(all[iS][iC]->GetBinCenter(i)>kCentPtLimits[iC]) continue;
          float data_val = all[iS][iC]->GetBinContent(i);
          float fit_val = scaledbw->GetBinContent(scaledbw->FindBin(all[iS][iC]->GetBinCenter(i)));
          float data_err = stat_all[iS][iC]->GetBinError(i);
          float pull = (data_val-fit_val)/data_err;
          //printf("iC: %d, iB: %d, data_val: %f, fit_val: %f, data_err: %f, pull: %f\n", iC, i, data_val, fit_val, data_err, pull);
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
    // TPaveText* pav_fin = new TPaveText(0.69,0.403478,0.951128,0.530435,"blNDC");
    // pav_fin->SetBorderSize(0);
    // pav_fin->SetFillColor(10);
    // pav_fin->SetTextSize(0.03);
    // pav_fin->SetTextAlign(31);
    // pav_fin->AddText("ALICE Preliminary");
    // pav_fin->AddText(Form("#bf{%s, pp, #sqrt{#it{s}} = 13 TeV}",kNames[iS].data()));
    // pav_fin->AddText("#bf{V0M Multiplicity Classes}");
    // pav_fin->Draw();
    spectra.SetLogy();
    s_dir->cd();
    spectra.Write();
    //spectra.SaveAs(Form("~/Desktop/Figure_QM_dopo_norm/spectra_%c.pdf",kLetter[iS]));
    TCanvas cComp("cComp","cComp",3200,2400);
    cComp.DrawFrame(0.85,-8.,1.25,8.,";#it{p}_{T} (GeV/#it{c});#frac{N_{TPC} - N_{TOF}}{#sigma}");
    TLegend comp_leg(0.80,0.17,0.94,0.45);
    comp_leg.SetBorderSize(0);
    comp_leg.SetHeader(Form("%s, pp #sqrt{s} = 13 TeV",kParticleNames[iS]));
    for (int iC = 0; iC < n_centralities -1; ++iC) {
      comp[iS][iC]->Draw("psame");
      comp_leg.AddEntry(comp[iS][iC],Form("%4.0f - %2.0f %%",kCentLabels[iC][0],kCentLabels[iC][1]),"p");
    }
    for(int i=-3; i<=3; i++){
      TLine *line = new TLine(0.85,i,1.25,i);
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
}
