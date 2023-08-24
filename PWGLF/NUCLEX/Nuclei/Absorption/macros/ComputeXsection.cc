#include "TCanvas.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TF1.h"
#include "TLegend.h"
#include "TGraph.h"
#include "TList.h"
#include "TStyle.h"

#include <iostream>
#include <string>

#include "AliPID.h"
#include "AliAnalysisTaskDeuteronAbsorption.h"

constexpr int kRebin = 4;
constexpr double kRhoDx = 0.0006 * 0.794;

constexpr double Sq(double x)
{
  return x * x;
}

double FromRatioToCrossSection(double x)
{
  return -TMath::Log(x) / kRhoDx;
}

double FromRatioToCrossSectionUnc(double x, double ex)
{
  return ex / (x * kRhoDx);
}

void ComputeXsection(std::string fileName = "AnalysisResults.root")
{
  TFile paramFile("param.root");
  TGraph *paramLine{nullptr};
  if (paramFile.IsOpen())
  {
    paramLine = static_cast<TGraph *>(paramFile.Get("xsection_diff"));
  }

  TFile inputFile(fileName.data());
  TList *inputList = (TList *)inputFile.Get("DeuteronAbsorption/100");

  constexpr double massInterval = 0.2; /// 20% around the expected mass

  TH2F *fHist2Matching[kNabsSpecies][2][2];
  TH1D *fHist1MomSpectrum[kNabsSpecies][2][2];
  TH1D *fHist1MomSpectrumTPC[kNabsSpecies][2][2];
  TH1D *fHist1RatioTRDnoTRD[kNabsSpecies][2];
  TH1D *fHist1DoubleRatio[kNabsSpecies];
  TH1D *fHist1CrossSection[kNabsSpecies];
  TH1F *fHist1AcceptanceAll[2][2][2];
  TH1F *fHist1AcceptanceTRDnoTRD[2][2];

  std::string wTRD[2] = {"woTRD", "wTRD"};
  std::string wTOF[2] = {"woTOF", "wTOF"};
  std::string pos_neg[2] = {"neg", "pos"};
  std::string labels[kNabsSpecies][2] = {
      {"K^{-}", "K^{+}"},
      {"#bar{p}", "p"},
      {"^{2}#bar{H}", "^{2}H"},
      {"^{3}#bar{H}", "^{3}H"},
      {"^{3}#bar{He}", "^{3}He"}
  };
  double maxp[kNabsSpecies]{1.,1.,2.,2.5,7.0};

  TFile outputFile(Form("xsection_%s", fileName.data()), "recreate");
  TDirectory *spectraDir = outputFile.mkdir("Spectra");
  TDirectory *matchingDir = outputFile.mkdir("Matching");
  TDirectory *trdNoTrdDir = outputFile.mkdir("TRDnoTRD");
  TDirectory *acceptanceDir = outputFile.mkdir("Acceptance");
  TDirectory *doubleRatioDir = outputFile.mkdir("DoubleRatios");
  TDirectory *xSectionDir = outputFile.mkdir("xSections");

  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  TCanvas cvMatching("cvMatching", "Matching efficiency", 1600, 800);
  cvMatching.Divide(2, 2);
  TCanvas cvDoubleRatio("cvDoubleRatio", "Double ratios", 1600, 800);
  cvDoubleRatio.Divide(2, 2);
  for (int iCharge = 0; iCharge < 2; ++iCharge)
  {
    for (int iTRD = 0; iTRD < 2; ++iTRD)
    {
      for (int iTOF = 0; iTOF < 2; ++iTOF)
      {
        fHist1AcceptanceAll[iCharge][iTRD][iTOF] = static_cast<TH1F *>(inputList->FindObject(Form("fHist1AcceptanceAll%s_%s_%s", pos_neg[iCharge].data(), wTRD[iTRD].data(), wTOF[iTOF].data())));
        fHist1AcceptanceAll[iCharge][iTRD][iTOF]->Rebin(kRebin);
        fHist1AcceptanceAll[iCharge][iTRD][iTOF]->Sumw2();
      }
      fHist1AcceptanceAll[iCharge][iTRD][0]->Add(fHist1AcceptanceAll[iCharge][iTRD][1]);
    }
    acceptanceDir->cd();
    for (int iTOF = 0; iTOF < 2; ++iTOF)
    {
      fHist1AcceptanceTRDnoTRD[iCharge][iTOF] = static_cast<TH1F *>(fHist1AcceptanceAll[iCharge][1][iTOF]->Clone(Form("fHist1AcceptanceTRDnoTRD%s_%s", pos_neg[iCharge].data(), wTOF[iTOF].data())));
      fHist1AcceptanceTRDnoTRD[iCharge][iTOF]->Divide(fHist1AcceptanceAll[iCharge][0][iTOF]);
      fHist1AcceptanceTRDnoTRD[iCharge][iTOF]->Write();
    }
  }

  for (int iSpecies = 0; iSpecies < kNabsSpecies; ++iSpecies)
  {
    for (int iCharge = 0; iCharge < 2; ++iCharge)
    {
      for (int iTRD = 0; iTRD < 2; ++iTRD)
        fHist2Matching[iSpecies][iCharge][iTRD] = static_cast<TH2F *>(inputList->FindObject(Form("fHist2Matching%s_%s_%s", AliAnalysisTaskDeuteronAbsorption::fgkParticleNames[iSpecies].data(), pos_neg[iCharge].data(), wTRD[iTRD].data())));
      const double minMass = Sq(AliPID::ParticleMassZ(AliAnalysisTaskDeuteronAbsorption::fgkSpecies[iSpecies])) * (1. - massInterval);
      const int minBin = fHist2Matching[iSpecies][iCharge][0]->GetYaxis()->FindBin(minMass);
      const double maxMass = Sq(AliPID::ParticleMassZ(AliAnalysisTaskDeuteronAbsorption::fgkSpecies[iSpecies])) * (1. + massInterval);
      const int maxBin = fHist2Matching[iSpecies][iCharge][0]->GetYaxis()->FindBin(maxMass);
      std::cout << "Mass^2 limits used for " << AliAnalysisTaskDeuteronAbsorption::fgkParticleNames[iSpecies].data() << ": " << minMass << ", " << maxMass << std::endl;
      for (int iTRD = 0; iTRD < 2; ++iTRD)
      {
        spectraDir->cd();
        fHist1MomSpectrumTPC[iSpecies][iCharge][iTRD] = fHist2Matching[iSpecies][iCharge][iTRD]->ProjectionX(Form("fHist1MomSpectrumTPC%s_%s_%s", AliAnalysisTaskDeuteronAbsorption::fgkParticleNames[iSpecies].data(), pos_neg[iCharge].data(), wTRD[iTRD].data()));
        fHist1MomSpectrumTPC[iSpecies][iCharge][iTRD]->Rebin(kRebin);
        fHist1MomSpectrumTPC[iSpecies][iCharge][iTRD]->Sumw2();
        fHist1MomSpectrumTPC[iSpecies][iCharge][iTRD]->Write();
        fHist1MomSpectrum[iSpecies][iCharge][iTRD] = fHist2Matching[iSpecies][iCharge][iTRD]->ProjectionX(Form("fHist1MomSpectrum%s_%s_%s", AliAnalysisTaskDeuteronAbsorption::fgkParticleNames[iSpecies].data(), pos_neg[iCharge].data(), wTRD[iTRD].data()), minBin, maxBin);
        fHist1MomSpectrum[iSpecies][iCharge][iTRD]->Rebin(kRebin);
        fHist1MomSpectrum[iSpecies][iCharge][iTRD]->Sumw2();
        fHist1MomSpectrum[iSpecies][iCharge][iTRD]->Write();
        matchingDir->cd();
        fHist1MomSpectrumTPC[iSpecies][iCharge][iTRD]->Divide(fHist1MomSpectrum[iSpecies][iCharge][iTRD], fHist1MomSpectrumTPC[iSpecies][iCharge][iTRD], 1., 1., "B");
        fHist1MomSpectrumTPC[iSpecies][iCharge][iTRD]->Write(Form("fHist1TOFmatching%s_%s_%s", AliAnalysisTaskDeuteronAbsorption::fgkParticleNames[iSpecies].data(), pos_neg[iCharge].data(), wTRD[iTRD].data()));
        cvMatching.cd(iSpecies + 1);
        fHist1MomSpectrumTPC[iSpecies][iCharge][iTRD]->GetXaxis()->SetRangeUser(0,maxp[iSpecies]);
        fHist1MomSpectrumTPC[iSpecies][iCharge][iTRD]->Draw("same PLC PMC");
      }
      trdNoTrdDir->cd();
      fHist1RatioTRDnoTRD[iSpecies][iCharge] = static_cast<TH1D *>(fHist1MomSpectrum[iSpecies][iCharge][1]->Clone(Form("fHist1RatioTRDnoTRD%s_%s", AliAnalysisTaskDeuteronAbsorption::fgkParticleNames[iSpecies].data(), pos_neg[iCharge].data())));
      fHist1RatioTRDnoTRD[iSpecies][iCharge]->Divide(fHist1MomSpectrum[iSpecies][iCharge][0]);
      fHist1RatioTRDnoTRD[iSpecies][iCharge]->SetTitle(Form("%s; #it{p} (GeV/#it{c}); TRD / no TRD", labels[iSpecies][iCharge].data()));
      fHist1RatioTRDnoTRD[iSpecies][iCharge]->Write();
      fHist1RatioTRDnoTRD[iSpecies][iCharge]->Divide(fHist1AcceptanceTRDnoTRD[iCharge][1]);
    }
    doubleRatioDir->cd();
    fHist1DoubleRatio[iSpecies] = static_cast<TH1D *>(fHist1RatioTRDnoTRD[iSpecies][0]->Clone(Form("fHist1DoubleRatio%s", AliAnalysisTaskDeuteronAbsorption::fgkParticleNames[iSpecies].data())));
    fHist1DoubleRatio[iSpecies]->Divide(fHist1RatioTRDnoTRD[iSpecies][1]);
    fHist1DoubleRatio[iSpecies]->SetTitle(Form("; #it{p} (GeV/#it{c}); (TRD / no TRD)_{%s} / (TRD / no TRD)_{%s}", labels[iSpecies][0].data(), labels[iSpecies][1].data()));
    fHist1DoubleRatio[iSpecies]->Write();
    cvDoubleRatio.cd(iSpecies + 1);
    fHist1DoubleRatio[iSpecies]->Draw("same PLC PMC");
    fHist1DoubleRatio[iSpecies]->GetXaxis()->SetRangeUser(0,maxp[iSpecies]);
    fHist1CrossSection[iSpecies] = static_cast<TH1D *>(fHist1DoubleRatio[iSpecies]->Clone(Form("fHist1CrossSection%s", AliAnalysisTaskDeuteronAbsorption::fgkParticleNames[iSpecies].data())));
    for (int iBin = 1; iBin <= fHist1CrossSection[iSpecies]->GetNbinsX(); ++iBin)
    {
      double ratio = fHist1CrossSection[iSpecies]->GetBinContent(iBin);
      if (ratio < 1.e-8)
        continue;
      fHist1CrossSection[iSpecies]->SetBinContent(iBin, FromRatioToCrossSection(ratio));
      fHist1CrossSection[iSpecies]->SetBinError(iBin, FromRatioToCrossSectionUnc(ratio, fHist1CrossSection[iSpecies]->GetBinError(iBin)));
    }
    xSectionDir->cd();
    fHist1CrossSection[iSpecies]->SetTitle(Form("; #it{p} (GeV/#it{c}); #sigma_{%s} - #sigma_{%s} (mb)", labels[iSpecies][0].data(), labels[iSpecies][1].data()));
    fHist1CrossSection[iSpecies]->Write();
  }
  for (int iPad = 0; iPad < kNabsSpecies; ++iPad) {
    cvMatching.cd(iPad + 1)->BuildLegend(0.15,0.65,0.38,0.88);
    cvDoubleRatio.cd(iPad + 1)->BuildLegend(0.15,0.65,0.38,0.88);
  }
  cvMatching.SaveAs("QA/matching.pdf");
  cvDoubleRatio.SaveAs("QA/doubleRatio.pdf");
  outputFile.Close();
}
