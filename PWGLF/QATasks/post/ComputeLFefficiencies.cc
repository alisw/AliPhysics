#include <array>
#include <cmath>
#include <iostream>
#include <string>
#include <utility>

#include <TCanvas.h>
#include <TEfficiency.h>
#include <TFile.h>
#include <TList.h>
#include <TH1D.h>
#include <TH3D.h>

#include <AliPID.h>

#include "AliAnalysisTaskLFefficiencies.h"

const char* TOF_cut_names[7] = {"","","","FB5 + TOF matching", "FB5 + TOF pid", "FB5 + TOF matching - TOF mismatch","FB5 + TOF matching - TOF mismatch + TOF pid"};
const int kNoCheckTOFmismatch = 5;

void DivideBinomial(TH1* num, const TH1* den) {
  for (int iBin = 1; iBin <= num->GetNbinsX(); ++iBin) {
    const double n = num->GetBinContent(iBin);
    const double d = den->GetBinContent(iBin);
    const double eff = (d > 1.e-24) ? n / d : 0.;
    const double err_eff = (d > 1.e-24 && eff < 1) ? std::sqrt(eff * (1. - eff) / d) : 0.;
    num->SetBinContent(iBin, eff);
    num->SetBinError(iBin, err_eff);
  }
}

bool endsWith(const std::string &mainStr, const std::string &&toMatch)
{
	if(mainStr.size() >= toMatch.size() &&
			mainStr.compare(mainStr.size() - toMatch.size(), toMatch.size(), toMatch) == 0)
			return true;
		else
			return false;
}

std::pair<std::array<std::vector<TH1D*>,AliPID::kSPECIESC * 2 - 4>,std::array<std::vector<TH1D*>,AliPID::kSPECIESC * 2 - 4>>
ComputeLFefficiencies(std::string filename = "AnalysisResults") {
  if (endsWith(filename,".root")) {
    for (int i = 0; i < 5; ++i) filename.pop_back();
  }

  TFile input_file((filename + ".root").data());
  TList* input_list = static_cast<TList*>(input_file.Get("PWGLF_QA/efficiencies"));
  if (!input_list) {
    ::Fatal("ComputeLFefficiencies","Missing input list, check the input file.");
  }

  std::pair<std::array<std::vector<TH1D*>,AliPID::kSPECIESC * 2 - 4>,std::array<std::vector<TH1D*>,AliPID::kSPECIESC * 2 - 4>> vec;
  TCanvas cv("canvas");
  TCanvas cv_eta("canvas_eta");
  TCanvas cv_rec("reconstructed");
  TCanvas cv_gen("generated");
  cv.Print((filename + "_y.pdf[").data(),"pdf");
  cv_eta.Print((filename + "_eta.pdf[").data(),"pdf");

  TCanvas cv_tof("canvas_tof");
  cv_tof.Print((filename + "_tof_y.pdf[").data(),"pdf");
  TCanvas cv_tof_eta("canvas_tof_eta");
  cv_tof_eta.Print((filename + "_tof_eta.pdf[").data(),"pdf");

  TFile output_file((filename + "_out.root").data(),"recreate");

  for (int iSpecies = 2; iSpecies < AliPID::kSPECIESC; ++iSpecies) {
    for (int iCharge = 0; iCharge < 2; ++iCharge) {
      cv.cd();
      cv.DrawFrame(0.,0.,6.,1.1,Form("%s %s;#it{p}_{T} (GeV/#it{c}); Efficiency",AliPID::ParticleLatexName(iSpecies),AliAnalysisTaskLFefficiencies::fPosNeg[iCharge].data()));
      TH3D* gen3D = static_cast<TH3D*>(input_list->FindObject(Form("Gen_%s_%s",AliPID::ParticleShortName(iSpecies),AliAnalysisTaskLFefficiencies::fPosNeg[iCharge].data())));
      TH1D* gen = gen3D->ProjectionZ(Form("Gen_%s_%s_pz",AliPID::ParticleShortName(iSpecies),AliAnalysisTaskLFefficiencies::fPosNeg[iCharge].data()),3,7);

      cv_eta.cd();
      cv_eta.DrawFrame(0.,0.,6.,1.1,Form("%s %s |#eta|<0.8;#it{p}_{T} (GeV/#it{c}); Efficiency",AliPID::ParticleLatexName(iSpecies),AliAnalysisTaskLFefficiencies::fPosNeg[iCharge].data()));
      TH3D* gen3Deta = static_cast<TH3D*>(input_list->FindObject(Form("GenEta_%s_%s",AliPID::ParticleShortName(iSpecies),AliAnalysisTaskLFefficiencies::fPosNeg[iCharge].data())));
      TH1D* genEta = gen3Deta->ProjectionZ(Form("GenEta_%s_%s_pz",AliPID::ParticleShortName(iSpecies),AliAnalysisTaskLFefficiencies::fPosNeg[iCharge].data()),2,9);

      TH1D *rec[AliAnalysisTaskLFefficiencies::fNcuts];
      TH1D *recEta[AliAnalysisTaskLFefficiencies::fNcuts];

      for (int iCut = 0; iCut < AliAnalysisTaskLFefficiencies::fNcuts; ++iCut) {
        cv.cd();
        TH3D* rec3D = static_cast<TH3D*>(input_list->FindObject(Form("Rec_%s_%s_%i",AliPID::ParticleShortName(iSpecies),AliAnalysisTaskLFefficiencies::fPosNeg[iCharge].data(),iCut)));
        rec[iCut] = rec3D->ProjectionZ(Form("Rec_%s_%s_%i",AliPID::ParticleShortName(iSpecies),AliAnalysisTaskLFefficiencies::fPosNeg[iCharge].data(),iCut),3,7);
        TH1D* eff = (TH1D*)rec[iCut]->Clone(Form("Eff_%s_%s_%i",AliPID::ParticleShortName(iSpecies),AliAnalysisTaskLFefficiencies::fPosNeg[iCharge].data(),iCut));
        if(iCut==5){
          eff->SetTitle("FB5 + TOF matching - mismatch");
        }
        DivideBinomial(eff,gen);
        if(iCut < kNoCheckTOFmismatch) eff->Draw("same PLC PMC");
        vec.first[(iSpecies - 2) * 2 + iCharge].push_back(static_cast<TH1D*>(eff->Clone()));
        vec.first[(iSpecies - 2) * 2 + iCharge].back()->SetDirectory(0);

        cv_eta.cd();
        TH3D* rec3Deta = static_cast<TH3D*>(input_list->FindObject(Form("RecEta_%s_%s_%i",AliPID::ParticleShortName(iSpecies),AliAnalysisTaskLFefficiencies::fPosNeg[iCharge].data(),iCut)));
        recEta[iCut] = rec3Deta->ProjectionZ(Form("RecEta_%s_%s_%i",AliPID::ParticleShortName(iSpecies),AliAnalysisTaskLFefficiencies::fPosNeg[iCharge].data(),iCut),2,9);
        TH1D* effEta = (TH1D*)recEta[iCut]->Clone(Form("EffEta_%s_%s_%i",AliPID::ParticleShortName(iSpecies),AliAnalysisTaskLFefficiencies::fPosNeg[iCharge].data(),iCut));
        DivideBinomial(effEta,genEta);
        if(iCut < kNoCheckTOFmismatch) effEta->Draw("same PLC PMC");

        output_file.cd();
        eff->Write();
        effEta->Write();
        vec.second[(iSpecies - 2) * 2 + iCharge].push_back(static_cast<TH1D*>(eff->Clone()));
        vec.second[(iSpecies - 2) * 2 + iCharge].back()->SetDirectory(0);
      }

      cv.cd();
      cv.BuildLegend();
      cv.Print((filename + "_y.pdf").data(),"pdf");

      cv_eta.cd();
      cv_eta.BuildLegend();
      cv_eta.Print((filename + "_eta.pdf").data(),"pdf");

      // TOF matching efficiency

      cv_tof.cd();
      cv_tof.DrawFrame(0.,0.,6.,1.1,Form("%s %s;#it{p}_{T} (GeV/#it{c}); TOF Efficiency",AliPID::ParticleLatexName(iSpecies),AliAnalysisTaskLFefficiencies::fPosNeg[iCharge].data()));

      cv_tof_eta.cd();
      cv_tof_eta.DrawFrame(0.,0.,6.,1.1,Form("%s %s |#eta|<0.8;#it{p}_{T} (GeV/#it{c}); TOF Efficiency",AliPID::ParticleLatexName(iSpecies),AliAnalysisTaskLFefficiencies::fPosNeg[iCharge].data()));

      for(int iCut = 3; iCut < AliAnalysisTaskLFefficiencies::fNcuts; ++iCut) {
        cv_tof.cd();
        TH1D* hTOFmatchEff = (TH1D*)rec[iCut]->Clone(Form("TOFmatchEff_%s_%s_%i",AliPID::ParticleShortName(iSpecies),AliAnalysisTaskLFefficiencies::fPosNeg[iCharge].data(),iCut));
        hTOFmatchEff->SetTitle(Form("(%s) / FB5",TOF_cut_names[iCut]));
        DivideBinomial(hTOFmatchEff,rec[0]);
        hTOFmatchEff->Draw("same PLC PMC");

        cv_tof_eta.cd();
        TH1D* hTOFmatchEffEta = (TH1D*)recEta[iCut]->Clone(Form("TOFmatchEffEta_%s_%s_%i",AliPID::ParticleShortName(iSpecies),AliAnalysisTaskLFefficiencies::fPosNeg[iCharge].data(),iCut));
        if(iCut==5){
          hTOFmatchEffEta->Add(recEta[3],-1);
          hTOFmatchEffEta->Scale(-1);
        }
        hTOFmatchEffEta->SetTitle(Form("(%s) / FB5",TOF_cut_names[iCut]));
        DivideBinomial(hTOFmatchEffEta,recEta[0]);
        hTOFmatchEffEta->Draw("same PLC PMC");

      }

      cv_tof.cd();
      cv_tof.BuildLegend();
      cv_tof.Print((filename + "_tof_y.pdf").data(),"pdf");

      cv_tof_eta.cd();
      cv_tof_eta.BuildLegend();
      cv_tof_eta.Print((filename + "_tof_eta.pdf").data(),"pdf");
    }
  }
  cv.Print((filename + "_y.pdf]").data(),"pdf");
  cv_eta.Print((filename + "_eta.pdf]").data(),"pdf");

  cv_tof.Print((filename + "_tof_y.pdf]").data(),"pdf");
  cv_tof_eta.Print((filename + "_tof_eta.pdf]").data(),"pdf");

  output_file.Close();
  return vec;
}

void ComputeLFefficiencies(std::string filename0, std::string filename1) {
  std::array<std::string, 2> filenames { filename0, filename1 };
  std::pair<std::array<std::vector<TH1D*>,AliPID::kSPECIESC * 2 - 4>,std::array<std::vector<TH1D*>,AliPID::kSPECIESC * 2 - 4>> efficiencies[2];

  for (size_t iFile = 0; iFile < filenames.size(); ++iFile) {
    std::string& fname = filenames[iFile];
    efficiencies[iFile] = ComputeLFefficiencies(fname);
  }

  TCanvas cv_y("canvas_y");
  TCanvas cv_eta("canvas_eta");
  cv_y.Print("Comparison_y.pdf[","pdf");
  cv_eta.Print("Comparison_eta.pdf[","pdf");

  for (int iSpecies = 2; iSpecies < AliPID::kSPECIESC; ++iSpecies) {
    for (int iCharge = 0; iCharge < 2; ++iCharge) {
      cv_y.cd();
      cv_y.DrawFrame(0.,0.,6.,2.1,Form("%s %s;#it{p}_{T} (GeV/#it{c}); Ratio",AliPID::ParticleLatexName(iSpecies),AliAnalysisTaskLFefficiencies::fPosNeg[iCharge].data()));
      cv_eta.cd();
      cv_eta.DrawFrame(0.,0.,6.,2.1,Form("%s %s |#eta|<0.8;#it{p}_{T} (GeV/#it{c}); Ratio",AliPID::ParticleLatexName(iSpecies),AliAnalysisTaskLFefficiencies::fPosNeg[iCharge].data()));

      for (size_t iHist = 0; iHist < efficiencies[0].first[(iSpecies-2)*2+iCharge].size(); ++iHist) {
        cv_y.cd();
        efficiencies[0].first[(iSpecies-2)*2+iCharge][iHist]->Divide(efficiencies[1].first[(iSpecies-2)*2+iCharge][iHist]);
        efficiencies[0].first[(iSpecies-2)*2+iCharge][iHist]->Draw("PMC PLC same");

        cv_eta.cd();
        efficiencies[0].second[(iSpecies-2)*2+iCharge][iHist]->Divide(efficiencies[1].second[(iSpecies-2)*2+iCharge][iHist]);
        efficiencies[0].second[(iSpecies-2)*2+iCharge][iHist]->Draw("PMC PLC same");

      }
      cv_y.BuildLegend();
      cv_eta.BuildLegend();
      cv_y.Print("Comparison_y.pdf","pdf");
      cv_eta.Print("Comparison_eta.pdf","pdf");

    }
  }

  cv_y.Print("Comparison_y.pdf]","pdf");
  cv_eta.Print("Comparison_eta.pdf]","pdf");
}
