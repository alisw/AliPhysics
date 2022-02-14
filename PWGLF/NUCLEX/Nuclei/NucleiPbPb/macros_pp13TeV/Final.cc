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

const char* kPrefix[2] = {"deuterons","antideuterons"};

const char* kParticleNames[2] = {"Deuterons", "Antideuterons"};

void Final(const char* fitFunctionName = "LevyTsallis", bool no_function = false) {
  TFile spectra_file(kSpectraOutput.data());
  TFile systematics_file(kJoinSystematicsOutput.data());
  TFile final_file(kFinalOutput.data(),"recreate");

  //histograms with corrected spectra
  TH1F* stat[2][kCentLength];
  TH1F* syst[2][kCentLength];
  TH1F* syst_pt_uncorr[2][kCentLength];
  TH1F* syst_pt_corr[2][kCentLength];
  TH1F* syst_mult_corr[2][kCentLength];
  TH1F* syst_mult_uncorr[2][kCentLength];

  for (int iS = 0; iS < 2; ++iS) {
    TDirectory* s_dir = final_file.mkdir(kNames[iS].data());
    for (int iC = 0; iC < kCentLength; ++iC) {
      TDirectory *c_dir = s_dir->mkdir(std::to_string(iC).data());
      c_dir->cd();
      string basepath = kFilterListNames + "/" + kNames[iS] + "/" + std::to_string(iC) + "/Joined/JoinedSpectra" + kLetter[iS] + std::to_string(iC);
      stat[iS][iC] = (TH1F*)spectra_file.Get(basepath.data());
      Requires(stat[iS][iC],basepath.data());
      auto ptAxis = stat[iS][iC]->GetXaxis();
      //
      string systpath = kNames[iS] + "/" + std::to_string(iC) + "/joined/totsyst_" + kLetter[iS] + std::to_string(iC);
      TH1F* syst_rel = (TH1F*)systematics_file.Get(systpath.data());
      Requires(syst_rel,systpath.data());
      syst[iS][iC]  = (TH1F*)stat[iS][iC]->Clone(("syst" + std::to_string(iC)).data());
      //
      string pt_uncorr_sytpath = kNames[iS] + "/" + std::to_string(iC) + "/joined/pt_uncorr_" + kLetter[iS] + std::to_string(iC);
      TH1F* pt_uncorr_syst_rel = (TH1F*)systematics_file.Get(pt_uncorr_sytpath.data());
      Requires(pt_uncorr_syst_rel, pt_uncorr_sytpath.data());
      syst_pt_uncorr[iS][iC]  = (TH1F*)stat[iS][iC]->Clone(("syst_pt_uncorr" + std::to_string(iC)).data());
      //
      string pt_corr_sytpath = kNames[iS] + "/" + std::to_string(iC) + "/joined/pt_corr_" + kLetter[iS] + std::to_string(iC);
      TH1F* pt_corr_syst_rel = (TH1F*)systematics_file.Get(pt_corr_sytpath.data());
      Requires(pt_corr_syst_rel, pt_corr_sytpath.data());
      syst_pt_corr[iS][iC]  = (TH1F*)stat[iS][iC]->Clone(("syst_pt_corr" + std::to_string(iC)).data());
      //
      string mult_corr_sytpath = kNames[iS] + "/" + std::to_string(iC) + "/joined/mult_corr_" + kLetter[iS] + std::to_string(iC);
      TH1F* mult_corr_syst_rel = (TH1F*)systematics_file.Get(mult_corr_sytpath.data());
      Requires(mult_corr_syst_rel,mult_corr_sytpath.data());
      syst_mult_corr[iS][iC]  = (TH1F*)stat[iS][iC]->Clone(("syst_mult_corr" + std::to_string(iC)).data());
      //
      string mult_uncorr_sytpath = kNames[iS] + "/" + std::to_string(iC) + "/joined/mult_uncorr_" + kLetter[iS] + std::to_string(iC);
      TH1F* mult_uncorr_syst_rel = (TH1F*)systematics_file.Get(mult_uncorr_sytpath.data());
      Requires(mult_uncorr_syst_rel,mult_uncorr_sytpath.data());
      syst_mult_uncorr[iS][iC]  = (TH1F*)stat[iS][iC]->Clone(("syst_mult_uncorr" + std::to_string(iC)).data());
      //
      for (int iB = 1; iB <= kNPtBins; ++iB) {
        if (ptAxis->GetBinCenter(iB) < kPtRange[0]|| ptAxis->GetBinCenter(iB) > kCentPtLimits[iC]){
          continue;
        }
        syst[iS][iC]->SetBinError(iB,syst_rel->GetBinContent(iB) * stat[iS][iC]->GetBinContent(iB));
        syst_pt_uncorr[iS][iC]->SetBinError(iB,pt_uncorr_syst_rel->GetBinContent(iB) * stat[iS][iC]->GetBinContent(iB));
        syst_pt_corr[iS][iC]->SetBinError(iB,pt_corr_syst_rel->GetBinContent(iB) * stat[iS][iC]->GetBinContent(iB));
        syst_mult_corr[iS][iC]->SetBinError(iB,mult_corr_syst_rel->GetBinContent(iB) * stat[iS][iC]->GetBinContent(iB));
        syst_mult_uncorr[iS][iC]->SetBinError(iB,mult_uncorr_syst_rel->GetBinContent(iB) * stat[iS][iC]->GetBinContent(iB));
      }
      plotting::SetHistStyle(stat[iS][iC],plotting::kSpectraColors[iC]);
      plotting::SetHistStyle(syst[iS][iC],plotting::kSpectraColors[iC]);
      plotting::SetHistStyle(syst_pt_uncorr[iS][iC],plotting::kSpectraColors[iC]);
      plotting::SetHistStyle(syst_pt_corr[iS][iC],plotting::kSpectraColors[iC]);
      plotting::SetHistStyle(syst_mult_corr[iS][iC],plotting::kSpectraColors[iC]);
      plotting::SetHistStyle(syst_mult_uncorr[iS][iC],plotting::kSpectraColors[iC]);
      stat[iS][iC]->Write("stat");
      stat[iS][iC]->Scale(kScaleFactor[iC]);
      syst[iS][iC]->Write("syst");
      syst[iS][iC]->Scale(kScaleFactor[iC]);
      syst_pt_uncorr[iS][iC]->Write("syst_pt_uncorr");
      syst_pt_uncorr[iS][iC]->Scale(kScaleFactor[iC]);
      syst_pt_corr[iS][iC]->Write("syst_pt_corr");
      syst_pt_corr[iS][iC]->Scale(kScaleFactor[iC]);
      syst_mult_corr[iS][iC]->Write("syst_mult_corr");
      syst_mult_corr[iS][iC]->Scale(kScaleFactor[iC]);
      syst_mult_uncorr[iS][iC]->Write("syst_mult_uncorr");
      syst_mult_uncorr[iS][iC]->Scale(kScaleFactor[iC]);
    }

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
    //text.DrawText(0.5,7.5,"This work");
    float name_position = (iS==0) ? 3.3 : 3.0;
    text.DrawLatex(name_position,7.5,Form("#bf{%s, pp, #sqrt{#it{s}} = 13 TeV}",kNamePlot[iS].data()));
    text.DrawLatex(3.35461,0.41997,Form("#bf{#color[%d]{#LTd#it{N}_{ch} / d#it{#eta}#GT = 26.02}}",plotting::kSpectraColors[0]));
    text.DrawLatex(0.641832,1.02914e-06,Form("#bf{#color[%d]{#LTd#it{N}_{ch} / d#it{#eta}#GT = 2.55}}",plotting::kSpectraColors[8]));
    TLegend final_leg(0.70,0.23,0.94,0.56);
    final_leg.SetBorderSize(0);
    final_leg.SetTextSize(0.027);
    final_leg.SetHeader("V0M Multiplicity Classes");
    TLegendEntry *header = (TLegendEntry*)final_leg.GetListOfPrimitives()->First();
    header->SetTextAlign(22);
    final_leg.SetNColumns(2);
    TFile bwfile(Form("%s%s_fits.root",kBaseOutputDir.data(),kPrefix[iS]),"read");
    TF1* bw = nullptr;
    TH1F* scaledbw = nullptr;
    if(!no_function && bwfile.IsOpen()) {  
      for (int iC = 0; iC < kCentLength; ++iC) {
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
      }
    }
    for (int iC = 0; iC < kCentLength; ++iC) {
      stat[iS][iC]->Draw("esamex0");
      syst[iS][iC]->Draw("e2same");
      final_leg.AddEntry(syst[iS][iC],Form("%s (#times 2^{ %d})",kRomanLabels[iC],kExponent[iC]),"fp");
    }
    if(!no_function && bwfile.IsOpen()) {
      final_leg.AddEntry(scaledbw,"Individual fit","l");
    }
    final_leg.Draw();
    spectra.SetLogy();
    s_dir->cd();
    spectra.Write();
  }
}
