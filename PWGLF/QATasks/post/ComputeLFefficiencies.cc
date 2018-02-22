#include <cmath>
#include <iostream>
#include <string>

#include <TCanvas.h>
#include <TEfficiency.h>
#include <TFile.h>
#include <TList.h>
#include <TH1D.h>
#include <TH3D.h>

#include <AliPID.h>

#include "AliAnalysisTaskLFefficiencies.h"

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

void ComputeLFefficiencies(std::string filename = "AnalysisResults") {
  if (endsWith(filename,".root")) {
    for (int i = 0; i < 5; ++i) filename.pop_back();
  }

  TFile input_file((filename + ".root").data());
  TList* input_list = static_cast<TList*>(input_file.Get("PWGLF_QA/efficiencies"));
  if (!input_list) {
    ::Fatal("ComputeLFefficiencies","Missing input list, check the input file.");
  }


  TCanvas cv("canvas");
  TCanvas cv_eta("canvas_eta");
  TCanvas cv_rec("reconstructed");
  TCanvas cv_gen("generated");
  cv.Print((filename + "_y.pdf[").data(),"pdf");
  cv_eta.Print((filename + "_eta.pdf[").data(),"pdf");

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

      for (int iCut = 0; iCut < AliAnalysisTaskLFefficiencies::fNcuts; ++iCut) {
        cv.cd();
        TH3D* rec3D = static_cast<TH3D*>(input_list->FindObject(Form("Rec_%s_%s_%i",AliPID::ParticleShortName(iSpecies),AliAnalysisTaskLFefficiencies::fPosNeg[iCharge].data(),iCut)));
        TH1D* eff = rec3D->ProjectionZ(Form("Eff_%s_%s_%i",AliPID::ParticleShortName(iSpecies),AliAnalysisTaskLFefficiencies::fPosNeg[iCharge].data(),iCut),3,7);
        DivideBinomial(eff,gen);
        eff->Draw("same PLC PMC");

        cv_eta.cd();
        TH3D* rec3Deta = static_cast<TH3D*>(input_list->FindObject(Form("RecEta_%s_%s_%i",AliPID::ParticleShortName(iSpecies),AliAnalysisTaskLFefficiencies::fPosNeg[iCharge].data(),iCut)));
        TH1D* effEta = rec3Deta->ProjectionZ(Form("EffEta_%s_%s_%i",AliPID::ParticleShortName(iSpecies),AliAnalysisTaskLFefficiencies::fPosNeg[iCharge].data(),iCut),2,9);
        DivideBinomial(effEta,genEta);
        effEta->Draw("same PLC PMC");

        output_file.cd();
        eff->Write();
        effEta->Write();
      }
      cv.cd();
      cv.BuildLegend();
      cv.Print((filename + "_y.pdf").data(),"pdf");

      cv_eta.cd();
      cv_eta.BuildLegend();
      cv_eta.Print((filename + "_eta.pdf").data(),"pdf");
    }
  }
  cv.Print((filename + "_y.pdf]").data(),"pdf");
  cv_rec.Print((filename + "_eta.pdf]").data(),"pdf");

  output_file.Close();
}
