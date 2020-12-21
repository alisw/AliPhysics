#include "TArrayD.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH2.h"
#include "TH2D.h"
#include "TH3.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TStyle.h"

TH1D* GetEfficiency(const Char_t* particleType = "Pion");
TH1* calculateEffRatio(const Char_t* particleType = "Pion");

TH1D* GetEfficiency(const Char_t* particleType)
{
  //
  // get efficiency for a certain particle type
  //
  TFile* inFileMC = TFile::Open("./Output_2018_03_17_MC/MergeOutput_MC.root"); // open analysis results folder
  TList* listMC = (TList*)inFileMC->Get("chist");
  //
  TH1F* histPartGen = (TH1F*)listMC->FindObject(Form("fHist%sGen", particleType));
  TH1F* histPartReco = (TH1F*)listMC->FindObject(Form("fHist%sReco", particleType));
  histPartGen->Sumw2();
  histPartReco->Sumw2();
  //
  TH1D* effPart = (TH1D*)histPartReco->Clone(Form("eff%s", particleType));
  effPart->GetXaxis()->SetTitle("p_{T}(GeV/c)");
  effPart->GetYaxis()->SetTitle("All Reco. tracks/All Gen.tracks");
  effPart->GetXaxis()->SetRangeUser(0, 10);
  effPart->Divide(histPartReco, histPartGen, 1.0, 1.0, "B");
  //TCanvas* canvPartEff = new TCanvas(Form("canv%sEff", particleType), Form("canv%sEff", particleType));
  //effPart->DrawCopy("EP");
  //canvPartEff->SaveAs(Form("%s.eps", effPart->GetName()));
  effPart->SetDirectory(0);
 // delete listMC;
 // inFileMC->Close();
  return effPart;
}
TH1 * calculateEffRatio(const Char_t* particleType)
{
  //
  // calculate efficiencies Ratios
  //
  TH1* hnum = GetEfficiency(particleType);
  TH1* hden = GetEfficiency(Form("Anti%s", particleType));
  TCanvas* EffPart = new TCanvas(Form("eff%s",particleType),Form("eff%s",particleType));
  EffPart->DrawFrame(0,0,1,1,";pt;eff");
  TLegend *leg = new TLegend(.1, .7, .5, .9);
  hnum->SetLineColor(kRed);
  hnum->DrawCopy("EPsame");
  leg->AddEntry(hnum, Form("%s Eff.", particleType));
  hden->SetLineColor(kBlue);
  hden->DrawCopy("EPsame");
  leg->AddEntry(hden, Form("Anti%s Eff.", particleType));
  leg->Draw();

  hnum->Divide(hden);
  hnum->SetMaximum(1.3);
  hnum->SetMinimum(0.7);
  hnum->SetDirectory(0);
  hnum->SetName(Form("EffRatio%s", particleType));
  hnum->SetTitle(Form(";#it{p}_{T} (GeV/#it{c});RatioEff %s/Anti%s", particleType, particleType));
  TCanvas* canvRatio = new TCanvas(Form("cRatioEff%s", particleType), particleType);
  canvRatio->DrawFrame(0,0.5,1,1.5,";pT;pos/neg effRatios");
  hnum->Draw("EPsame");
  return hnum;
}
void DrawRatioEfficiency()
{
  gStyle->SetOptTitle(0);
  //calculateEffRatio("Pion");
  //calculateEffRatio("Kaon");
  calculateEffRatio("Prot");
}
