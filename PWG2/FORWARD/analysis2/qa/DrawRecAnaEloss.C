/**
 * @file   DrawRecAnaEloss.C
 * @author Christian Holm Christensen <cholm@dalsgaard.hehi.nbi.dk>
 * @date   Thu Jul  7 10:58:50 2011
 * 
 * @brief  Draw energ-loss before/after merging and used in the
 * density calculations 
 *
 * @ingroup pwg2_forward_scripts_qa
 */

#ifndef __CINT__
#include <TList.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TLine.h>
#include <TFile.h>
#include <TError.h>
#include <TParameter.h>
#include <TStyle.h>
#else
class TLatex;
#endif

/** 
 * Draw some text
 * 
 * @param l 
 * @param x 
 * @param y 
 * @param c1 
 * @param c2 
 *
 * @ingroup pwg2_forward_scripts_qa
 */
void 
DrawText(TLatex* l, Double_t x, Double_t& y, const char* c1, const char* c2)
{
  y -= 1.2*l->GetTextSize();
  l->DrawLatex(x,    y, c1);
  l->DrawLatex(x+.4, y, c2);
}
/** 
 * Draw the energy loss before/after mergin for a single ring
 * 
 * @param p      List 1 
 * @param p2     List 2
 * @param lowCut Low cut
 * @param d      Detector
 * @param r      Ring 
 *
 * @ingroup pwg2_forward_scripts_qa
 */
void
DrawRingRecAnaEloss(TList* p, TList* p2, Double_t lowCut, UShort_t d, Char_t r)
{
  if (!p) return;

  TList* ring = static_cast<TList*>(p->FindObject(Form("FMD%d%c",d,r)));
  if (!ring) { 
    Error("DrawRecAnaEloss", "List FMD%d%c not found in %s",d,r,p->GetName());
    return;
  }
  TList* ring2 = static_cast<TList*>(p2->FindObject(Form("FMD%d%c",d,r)));
  if (!ring2){
    Error("DrawRecAnaEloss","List FMD%d%c not found in %s",d,r,p2->GetName());
    return;
  }

  TH1* before = static_cast<TH1D*>(ring->FindObject("esdEloss"));
  if (!before) { 
    Error("DrawRingRecAnaEloss", "Histogram esdEloss not found in FMD%d%c",
	  d, r);
    return;
  }
  TH1* after = static_cast<TH1D*>(ring->FindObject("anaEloss"));
  if (!after) { 
    Error("DrawRingRecAnaEloss", "Histogram anaEloss not found in FMD%d%c",
	  d, r);
    return;
  }
  TH1* presented = static_cast<TH1D*>(ring2->FindObject("eloss"));
  if (!presented) {
    Error("DrawRingRecAnaEloss", "Histogram eloss not found in FMD%d%c",
	  d, r);
    return;
  }
  TH1* used = static_cast<TH1D*>(ring2->FindObject("elossUsed"));
  if (!used) {
    Error("DrawRingRecAnaEloss", "Histogram elossUsed not found in FMD%d%c",
	  d, r);
    return;
  }


  Int_t low = before->GetXaxis()->FindBin(lowCut);
  Int_t ib  = Int_t(before->Integral(low,before->GetNbinsX()));
  Int_t ia  = Int_t(after->Integral(low,after->GetNbinsX()));
  Int_t ip  = Int_t(presented->Integral(low,presented->GetNbinsX()));
  Int_t iu  = Int_t(used->Integral(low,used->GetNbinsX()));
  // before->Scale(1. / ib);
  // after->Scale(1. / ib);
  // presented->Scale(1. / ib);
  // used->Scale(1. / ib);

  gPad->SetLogy();
  gPad->SetFillColor(0);
  before->SetTitle(Form("FMD%d%c",d,r));
  before->Draw("");
  after->Draw("same");
  presented->Draw("same");
  used->Draw("same");
  
  // ib           = before->Integral(low,before->GetNbinsX());
  // ia           = after->Integral(low,after->GetNbinsX());
  // ip           = presented->Integral(low,presented->GetNbinsX());
  // iu           = used->Integral(low,used->GetNbinsX());
  Double_t ts  = 0.03;
  Double_t x   = gPad->GetLeftMargin() + .25;
  Double_t y   = 1-gPad->GetTopMargin()-gStyle->GetTitleH()+ts;
  TLatex*  ltx = new TLatex(x, y, Form("FMD%d%c", d, r));
  ltx->SetNDC();
  ltx->SetTextAlign(13);
  ltx->SetTextSize(ts);
  // ltx->Draw();
  // ltx->SetTextSize(.05);
  TString inte(Form("Integral [%4.2f,#infty]", lowCut));
  DrawText(ltx, x, y, Form("%s before:", inte.Data()), Form("%9d", ib));
  DrawText(ltx, x, y, Form("%s after:",  inte.Data()), Form("%9d", ia));
  DrawText(ltx, x, y, Form("%s input:",  inte.Data()), Form("%9d", ip));
  DrawText(ltx, x, y, Form("%s user:",   inte.Data()), Form("%9d", iu));
  TLine* l = new TLine;
  l->SetLineWidth(1);
  l->DrawLineNDC(x, y-0.9*ts, 1-gPad->GetRightMargin()-0.01, y-0.9*ts);
  if (ib != 0 && ia != 0) {
    DrawText(ltx, x, y, "Change (merging)", Form("%5.1f%%", (100.*(ia-ib))/ib));
    DrawText(ltx, x, y, "Change (input)",   Form("%5.1f%% (%5.1f%%)", 
						 (100.*(ip-ia))/ia,
						 (100.*(ip-ib))/ib));
    DrawText(ltx, x, y, "Change (use)",     Form("%5.1f%% (%5.1f%%)", 
						 (100.*(iu-ip))/ip,
						 (100.*(iu-ib))/ib));
  }
  before->GetXaxis()->SetRangeUser(0, 4);
  gPad->cd();
}

/** 
 * Draw energy loss before/after merging 
 * 
 * @param filename 
 *
 * @ingroup pwg2_forward_scripts_qa
 */
void
DrawRecAnaEloss(const char* filename="forward.root",
		const char* folder="ForwardResults")
{
  gStyle->SetPalette(1);
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(1);
  gStyle->SetTitleW(.4);
  gStyle->SetTitleH(.1);
  gStyle->SetTitleColor(0);
  gStyle->SetTitleStyle(0);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleX(.6);
  
  TFile* file = TFile::Open(filename, "READ");
  if (!file) { 
    Error("DrawRecAnaEloss", "failed to open %s", filename);
    return;
  }

  TList* forward = static_cast<TList*>(file->Get(folder));
  if (!forward) { 
    Error("DrawRecAnaEloss", "List %s not found in %s", folder, filename);
    return;
  }

  TList* sf = static_cast<TList*>(forward->FindObject("fmdSharingFilter"));
  if (!sf) { 
    Error("DrawRecAnaEloss", "List fmdSharingFilter not found in Forward");
    return;
  }
  TList* dc = static_cast<TList*>(forward->FindObject("fmdDensityCalculator"));
  if (!dc) {
    Error("DrawRecAnaEloss","List fmdDensityCalculator not found in Forward");
    return;
  }

  TParameter<double>* lowCut = 
    static_cast<TParameter<double>*>(sf->FindObject("lowCut"));
  Double_t low = (lowCut ? lowCut->GetVal() : 0.15);
  if (!lowCut)
    Warning("DrawRecAnaEloss", "Low cut not found in %s, assuming %f",
	    sf->GetName(), low);
  TCanvas* c = new TCanvas("recAnaELoss", 
			   "Reconstructed and Analysed energy loss", 900, 700);
  c->SetFillColor(0);
  c->SetBorderSize(0);
  c->SetLeftMargin(0.15);
  c->SetRightMargin(0.02);
  c->SetTopMargin(0.02);
  c->Divide(3, 2, 0, 0);
  
  c->cd(1); DrawRingRecAnaEloss(sf, dc, low, 1, 'I');
  c->cd(2); DrawRingRecAnaEloss(sf, dc, low, 2, 'I');
  c->cd(5); DrawRingRecAnaEloss(sf, dc, low, 2, 'O');
  c->cd(3); DrawRingRecAnaEloss(sf, dc, low, 3, 'I');
  c->cd(6); DrawRingRecAnaEloss(sf, dc, low, 3, 'O');
  TVirtualPad* p = c->cd(4);
  // p->SetTopMargin(0.05);
  p->SetRightMargin(0.15);
  p->SetFillColor(0);
  TH2D* highCuts = static_cast<TH2D*>(sf->FindObject("highCuts"));
  if (highCuts) highCuts->Draw("colz");
  c->cd();
  c->SaveAs("recAnaELoss.png");
}

  
  
 
//
// EOF
//
