#include <cmath>
#include <iostream>

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

void ComputeLFefficiencies() {
  TFile input_file("AnalysisResults.root");
  TList* input_list = static_cast<TList*>(input_file.Get("PWGLF_QA/efficiencies"));
  if (!input_list) {
    ::Fatal("ComputeLFefficiencies","Missing input list, check the input file.");
  }

  TCanvas cv("canvas");
  cv.Print("output.pdf[","pdf");
  for (int iSpecies = 2; iSpecies < AliPID::kSPECIESC; ++iSpecies) {
    for (int iCharge = 0; iCharge < 2; ++iCharge) {
      cv.DrawFrame(0.,0.,6.,1.1,Form("%s %s;#it{p}_{T} (GeV/#it{c}); Efficiency",AliPID::ParticleLatexName(iSpecies),AliAnalysisTaskLFefficiencies::fPosNeg[iCharge].data()));
      TH3D* gen3D = static_cast<TH3D*>(input_list->FindObject(Form("Gen_%s_%s",AliPID::ParticleShortName(iSpecies),AliAnalysisTaskLFefficiencies::fPosNeg[iCharge].data())));
      TH1D* gen = gen3D->ProjectionZ(Form("Gen_%s_%s_pz",AliPID::ParticleShortName(iSpecies),AliAnalysisTaskLFefficiencies::fPosNeg[iCharge].data()),3,7);
      for (int iCut = 0; iCut < AliAnalysisTaskLFefficiencies::fNcuts; ++iCut) {
        TH3D* rec3D = static_cast<TH3D*>(input_list->FindObject(Form("Rec_%s_%s_%i",AliPID::ParticleShortName(iSpecies),AliAnalysisTaskLFefficiencies::fPosNeg[iCharge].data(),iCut)));
        TH1D* eff = rec3D->ProjectionZ(Form("Eff_%s_%s_%i",AliPID::ParticleShortName(iSpecies),AliAnalysisTaskLFefficiencies::fPosNeg[iCharge].data(),iCut),3,7);
        DivideBinomial(eff,gen);
        eff->Draw("C same PLC PMC");
      }
      cv.BuildLegend();
      cv.Print("output.pdf","pdf");
    }
  }
  cv.Print("output.pdf]","pdf");
}
