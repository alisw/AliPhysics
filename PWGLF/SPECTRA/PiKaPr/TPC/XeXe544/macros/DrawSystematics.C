//
//calculates the systematic uncertainties by taking ratio of multi-gaussian fit to fixed sigma cut fit
//
#include "TCanvas.h"
#include "TColor.h"
#include "TFile.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH2.h"
#include "TH2D.h"
#include "TH3.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TMath.h"
#include "TStyle.h"
#include "TGraph.h"
#include <vector>


TH1 *GetGaus(const Char_t * particleType ="Pion", Int_t cent = 0);
TH1 *GetFixCut(const Char_t * particleType ="Pion", Int_t cent = 0);
void SetColor(TH1*h, Int_t cent = 0);
const Char_t *GetCent(Int_t cent = 0);
TH1* DoRatio(const Char_t * particleType ="Pion", Int_t cent = 0);
TCanvas* DrawRatioSpectra(const Char_t * particleType ="Pion");

TH1 *GetGaus(const Char_t * particleType, Int_t cent){
  TVirtualPad *currpad1 = gPad;
  TFile fGaus(Form("FinalSpectra_forQM/canvPtDist%s.root",particleType), "READ");
  Printf("%s",fGaus.GetName());
  TCanvas *cGaus = (TCanvas*) fGaus.Get(Form("canvPtDist%s",particleType));
  Printf("%s",cGaus->GetName());
  TList * lGaus = cGaus->GetListOfPrimitives();
  lGaus->ls();
  TH1 * hGaus = (TH1*)lGaus->At(cent)->Clone(Form("%s%i", particleType, cent));
  hGaus->SetDirectory(0);
  fGaus.Close();
  currpad1->cd();
  return hGaus;
}

TH1 *GetFixCut(const Char_t * particleType, Int_t cent){
  TVirtualPad *currpad = gPad;
  TFile fFixCut(Form("Fix_Cut_Spectra/canvPtDist%s.root",particleType), "READ");
  Printf("%s",fFixCut.GetName());
  TCanvas *cFixCut = (TCanvas*) fFixCut.Get(Form("canvPtDist%s",particleType));
  Printf("%s",cFixCut->GetName());
  TList * lFixCut = cFixCut->GetListOfPrimitives();
  lFixCut->ls();
  TH1 * hFixCut = (TH1*)lFixCut->At(cent)->Clone(Form("%s%i", particleType, cent));
  hFixCut->SetDirectory(0);
  fFixCut.Close();
  currpad->cd();
  return hFixCut;
}

void SetColor(TH1*h, Int_t cent){
  const TString col[10] = {"#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a"};
  //http://colorbrewer2.org/#type=qualitative&scheme=Paired&n=10
  h->SetLineColor(TColor::GetColor(col[cent]));
  h->SetMarkerColor(TColor::GetColor(col[cent]));
}

const Char_t *GetCent(Int_t cent){
  Float_t min = 0,width = 0;
  for (Int_t i = 0; i <= cent; i++) {
    if (i < 2)
      width = 5;
    if (i >= 2 && i <= 7)
      width = 10;
    else if (i > 7)
      width = 20;
    if (cent == i)
      break;
    min += width;
  }

  return Form("%.0f-%.0f",min, min +width);
}

TH1* DoRatio(const Char_t * particleType, Int_t cent){
  TH1*hnum = GetGaus(particleType, cent);
  Printf("%s", hnum->GetName());
  TH1*hden = GetFixCut(particleType, cent);
  Printf("%s", hden->GetName());
  hnum->Divide(hden);
  hnum->SetMaximum(1.3);
  hnum->SetMinimum(0.7);
  hnum->SetDirectory(0);
  hnum->SetName(Form("Ratio%s%i",particleType,cent));
  hnum->SetTitle(Form("%s%%;#it{p}_{T} (GeV/#it{c});Relative Error",GetCent(cent) ));
  SetColor(hnum, cent);
  for (Int_t bin = 1; bin <= hnum->GetNbinsX(); bin++)
    hnum->SetBinContent(bin, TMath::Abs(hnum->GetBinContent(bin)-1));
    //hnum->Print("all");
  return hnum;
}

TCanvas* DrawRatioSpectra(const Char_t * particleType) {
  TCanvas *canvRatio = new TCanvas(Form("c%s", particleType), particleType);
  canvRatio->DrawFrame(0.1, 0, 1, .1, "Centrality;pT;error");
  for(Int_t i = 0;i<9; i++)
  DoRatio(particleType, i)->Draw("LPsame");
  canvRatio->BuildLegend(0.75,0.6,0.9,0.9);
  TLatex *label = new TLatex(.8, .91, particleType);
  label->SetNDC();
  label->Draw();
  canvRatio->SetGridy();
  canvRatio->SaveAs(Form("FinalSpectra_forQM/RelativeError_%s.png", particleType));
  canvRatio->SaveAs(Form("FinalSpectra_forQM/RelativeError_%s.root", particleType));
  return canvRatio;
}

void DrawSystematics() {
  gStyle->SetOptTitle(0);
  DrawRatioSpectra("Pion");
  DrawRatioSpectra("Kaon");
  DrawRatioSpectra("Proton");
  DrawRatioSpectra("AntiPion");
  DrawRatioSpectra("AntiKaon");
  DrawRatioSpectra("AntiProton");

}