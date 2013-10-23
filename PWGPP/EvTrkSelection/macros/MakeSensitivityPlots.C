#include "TH1D.h"
#include "TCanvas.h"
#include "TFile.h"
#include "THn.h"
#include "TList.h"
#include "TStyle.h"
#include "TLegend.h"

void MakeSensitivityPlots();
TH1D * GetAcceptedFractionNclCut(Int_t nclCut = 80, const Char_t * inFileName = "output/Data_LHC13b.root");

void MakeSensitivityPlots() {
  //
  //  -> THIS MACRO SHOULD BE COMPILABLE.
  //  -> ALL PLOTS SHOULD BE LABELED (ESPECIALLY THE AXES).
  //  -> DATA RED AND MC BLUE.
  //

  Int_t nclCut = 120;

  TH1D * nclAcceptedData80 = GetAcceptedFractionNclCut(nclCut, "output/Data_LHC13b.root");
  nclAcceptedData80->SetNameTitle(Form("nr clusters cut %d",nclCut),Form("nr clusters cut %d",nclCut));
  nclAcceptedData80->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  nclAcceptedData80->GetYaxis()->SetTitle("accepted fraction");
  nclAcceptedData80->GetXaxis()->SetTitleOffset(1.2);
  nclAcceptedData80->GetYaxis()->SetTitleOffset(1.2);
  nclAcceptedData80->GetYaxis()->SetTitleSize(0.045);
  nclAcceptedData80->SetMaximum(1.07);
  nclAcceptedData80->SetMinimum(0.57);
  nclAcceptedData80->SetLineColor(kRed -3);
  //
  TH1D * nclAcceptedMc80   = GetAcceptedFractionNclCut(nclCut, "output/MC_LHC13b.root");
  nclAcceptedMc80->SetLineColor(kBlue -3);
  //
  TCanvas * canvNclCut = new TCanvas("canvNclCut","sensitivity to ncl cut",600,800);
  gStyle->SetOptStat(0);
  canvNclCut->Divide(1,2,0,0);
  canvNclCut->cd(1)->SetLogx();
  gPad->SetTicky(2);
  gPad->SetFillStyle(0);
  //
  nclAcceptedData80->GetYaxis()->SetLabelFont(62);
  nclAcceptedData80->DrawCopy();
  nclAcceptedMc80->DrawCopy("SAME");
    
  TLegend *leg = new TLegend(0.75,0.2,.9,0.3);
  leg->AddEntry(nclAcceptedData80,"Data","f");
  leg->AddEntry(nclAcceptedMc80,"MC","f");
  leg->SetBorderSize(0);
  gStyle->SetFillColor(0);
  leg->Draw();
  //
  //
  //
    
  TH1D * nclAcceptedMcDataRatio = (TH1D*)nclAcceptedData80->Clone();
  //
  canvNclCut->cd(2)->SetLogx();
  gPad->SetTicky(2);
  gPad->SetFillStyle(0);
  //
  nclAcceptedMcDataRatio->Divide(nclAcceptedMc80);
  nclAcceptedMcDataRatio->GetYaxis()->SetTitle("ratio");
  nclAcceptedMcDataRatio->GetYaxis()->SetTitleSize(0.045);
  nclAcceptedMcDataRatio->GetYaxis()->SetTitleOffset(1.1);
  nclAcceptedMcDataRatio->SetMaximum(1.18);
  nclAcceptedMcDataRatio->SetMinimum(0.8);
  nclAcceptedMcDataRatio->GetYaxis()->SetLabelFont(62);
  nclAcceptedMcDataRatio->DrawCopy();

}


TH1D * GetAcceptedFractionNclCut(Int_t nclCut, const Char_t * inFileName) {
  //
  // accepted fraction of tracks for ncl cut vs. pT
  //
  TFile * inFileData = TFile::Open(inFileName);
  TList * l = (TList * ) inFileData->Get("akalweit_TrackingUncert");
  THnF * histNcl = (THnF *) l->FindObject("histNcl");
  //  histNcl->GetListOfAxes()->Print();
  //
  // determine sensitivities
  //
  TH1D * hAll = histNcl->Projection(1);
  hAll->SetNameTitle("hAll","hAll");  
  //
  const Int_t kVeryBig = 10000;
  histNcl->GetAxis(0)->SetRangeUser(nclCut, kVeryBig);
  TH1D * hAccepted = histNcl->Projection(1);
  hAccepted->SetNameTitle(Form("hAccepted%d",nclCut),Form("hAccepted%d",nclCut));
  //
  histNcl->GetAxis(0)->SetRangeUser(0,nclCut);
  TH1D * hRejected = histNcl->Projection(1);
  hRejected->SetNameTitle("hRejected","hRejected");
  //
  //
  hAccepted->Divide(hAll);
  hRejected->Divide(hAll);
  //
  // some cosmetics
  //
  hAccepted->SetLineWidth(2);
  //
  return hAccepted;

  
}
