#include "src/Common.h"
#include "src/Utils.h"
using namespace utils;
#include "src/Plotting.h"
using namespace plotting;

#include "TFile.h"
#include "TDirectory.h"
#include "TGraphAsymmErrors.h"
#include "TH1F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"

#include <string>

void HadronicSyst()
{
  TH1F *corr_geant_tpc[2];

  string g3g4_path = kBaseOutputDir + "G3G4.root";
  TFile g3g4_file(g3g4_path.data());
  TGraphAsymmErrors *g3g4tpc[2]{
      (TGraphAsymmErrors *)g3g4_file.Get("mCorrTpc"),
      (TGraphAsymmErrors *)g3g4_file.Get("aCorrTpc")};
  Requires(g3g4tpc[0], "tpcA");
  Requires(g3g4tpc[1], "tpcM");

  TFile fOut(Form("%sHadSyst.root", kBaseOutputDir.data()), "recreate");

  TCanvas *cG3G4 = new TCanvas("G3G4TPC", "G3G4TPC");
  TF1* function_g3g4_tpc[2];

  for (int iS = 0; iS < 2; iS++)
  {
    function_g3g4_tpc[iS]= new TF1(Form("function_g3g4_tpc_%c", kLetter[iS]), "pol0", 0.6, 1.4);
    if (iS == 0)
    {
      g3g4tpc[iS]->SetLineColor(kBlue);
      g3g4tpc[iS]->SetMarkerColor(kBlue);
      g3g4tpc[iS]->SetMarkerStyle(20);
      function_g3g4_tpc[iS]->SetLineColor(kBlue);
      function_g3g4_tpc[iS]->SetLineStyle(kDashed);
    }
    else
    {
      g3g4tpc[iS]->SetLineColor(kBlack);
      g3g4tpc[iS]->SetMarkerColor(kBlack);
      g3g4tpc[iS]->SetMarkerStyle(24);
      function_g3g4_tpc[iS]->SetLineColor(kBlack);
      function_g3g4_tpc[iS]->SetLineStyle(kDashed);
    }
    g3g4tpc[iS]->GetYaxis()->SetTitle("GEANT4 / GEANT3");
    g3g4tpc[iS]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    g3g4tpc[iS]->GetYaxis()->SetRangeUser(0.8, 1.05);
    g3g4tpc[iS]->GetXaxis()->SetRangeUser(0.6, 1.4);
    g3g4tpc[iS]->Fit(function_g3g4_tpc[iS], "RQ");
    g3g4tpc[iS]->Write();
    cG3G4->cd();
    if (iS == 0)
      g3g4tpc[iS]->Draw("AP");
    else
    {
      g3g4tpc[iS]->Draw("SAMEP");
      TLegend* leg = new TLegend(0.13,0.22,0.63,0.35 ,"","brNDC");
      leg->SetBorderSize(0.);
      leg->AddEntry(g3g4tpc[0], Form("%s : %.3f\n", kNames[0].data(), function_g3g4_tpc[0]->GetParameter(0)), "PE");
      leg->AddEntry(g3g4tpc[1], Form("%s : %.3f\n", kNames[iS].data(), function_g3g4_tpc[1]->GetParameter(0)), "PE");
      leg->Draw();
      cG3G4->Write();
    }
  }
}