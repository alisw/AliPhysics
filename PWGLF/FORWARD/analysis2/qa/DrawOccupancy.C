/**
 * @file   DrawOccupancy.C
 * @author Christian Holm Christensen <cholm@dalsgaard.hehi.nbi.dk>
 * @date   Thu Jul  7 10:24:58 2011
 * 
 * @brief  A script to draw the occupancy as given by Poisson method
 *
 * @deprecated Use QATrender instead
 * @ingroup pwg2_forward_scripts_qa
 * 
 */
#ifndef __CINT__
# include <TH1.h>
# include <TH2.h>
# include <TList.h>
# include <TFile.h>
# include <TString.h>
# include <TError.h>
# include <TPad.h>
# include <TCanvas.h>
# include <TLine.h>
# include <TLatex.h>
# include <TStyle.h>
#else
class TList;
#endif

/** 
 * Draw the Poisson estimate of the occupancy in a given ring.
 * 
 * @param p            List 
 * @param d            Detector
 * @param r            Ring
 * 
 * @return The occupancy (in percent)
 *
 * @deprecated Use QATrender instead
 * @ingroup pwg2_forward_scripts_qa
 */
Double_t
DrawRingOccupancy(TList* p, UShort_t d, Char_t r)
{
  if (!p) return 0;

  TList* ring = static_cast<TList*>(p->FindObject(Form("FMD%d%c",d,r)));
  if (!ring) { 
    Error("DrawOccupancy", "List FMD%d%c not found in %s",d,r,p->GetName());
    return 0;
  }
  
  TH1* corr = static_cast<TH1*>(ring->FindObject("occupancy"));
  if (!corr) { 
    Error("DrawRingOccupancy", "Histogram occupancy not found in FMD%d%c",
	  d, r);
    return 0;
  }
  corr->Rebin(4);

  TPad* pad = static_cast<TPad*>(gPad);
  pad->SetGridy();
  pad->SetGridx();
  pad->SetLogy();
  pad->SetFillColor(0);
    pad->SetRightMargin(0.01);
#if 0
  if (d == 3) { 
    pad->SetPad(pad->GetXlowNDC(), pad->GetYlowNDC(), .99, 
		 pad->GetYlowNDC()+pad->GetHNDC());
    pad->SetRightMargin(0.15);
  }
#endif

  corr->Draw("hist");

  TLatex* ltx = new TLatex(.95, .95, Form("FMD%d%c", d, r));
  ltx->SetNDC();
  ltx->SetTextAlign(33);
  ltx->SetTextSize(.08);
  ltx->Draw();

  return corr->GetMean();
}

/** 
 * Draw the Poisson estimate of the occupancy
 * 
 * @param filename Input file name 
 * @param folder   Input folder name in file
 * 
 * @deprecated Use QATrender instead
 * @ingroup pwg2_forward_scripts_qa
 */
void
DrawOccupancy(const char* filename="forward.root", 
	      const char* folder="ForwardResults")
{
  gStyle->SetPalette(1);
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetTitleW(.4);
  gStyle->SetTitleH(.1);
  gStyle->SetTitleX(.4);
  // gStyle->SetTitleY(.1);
  gStyle->SetTitleColor(0);
  gStyle->SetTitleStyle(0);
  gStyle->SetTitleBorderSize(0);
  
  TFile* file = TFile::Open(filename, "READ");
  if (!file) { 
    Error("DrawOccupancy", "failed to open %s", filename);
    return;
  }

  TList* forward = static_cast<TList*>(file->Get(folder));
  if (!forward) { 
    Error("DrawOccupancy", "List %s not found in %s", folder, filename);
    return;
  }

  TList* dc = static_cast<TList*>(forward->FindObject("fmdDensityCalculator"));
  if (!dc) { 
    Error("DrawOccupancy", "List fmdDensityCalculator not found in Forward");
    return;
  }
  
  TCanvas* c = new TCanvas("occupancy", 
			   "Mean Occupancy", 900, 700);
  c->SetFillColor(0);
  c->SetBorderSize(0);
  c->SetBorderMode(0);
  c->SetHighLightColor(0);
  c->SetBottomMargin(.15);
  c->SetTopMargin(.02);
  c->SetRightMargin(.02);
  c->SetLeftMargin(.15);
  c->Divide(3, 2, 0, 0);
  

  Double_t corrs[5];
  c->cd(1); corrs[0] = DrawRingOccupancy(dc, 1, 'I');
  c->cd(2); corrs[1] = DrawRingOccupancy(dc, 2, 'I');
  c->cd(5); corrs[2] = DrawRingOccupancy(dc, 2, 'O');
  c->cd(3); corrs[3] = DrawRingOccupancy(dc, 3, 'I');
  c->cd(6); corrs[4] = DrawRingOccupancy(dc, 3, 'O');

  TVirtualPad* p = c->cd(4);
  p->SetTopMargin(0.05);
  p->SetRightMargin(0.10);
  p->SetLeftMargin(0.15);
  p->SetBottomMargin(0.15);
  p->SetFillColor(0);

  TH1D* hc = new TH1D("occ", "Mean occupancy", 5, .5, 5.5);
  hc->SetFillColor(kRed+1);
  hc->SetFillStyle(3001);
  hc->SetMinimum(0.0);
  hc->GetXaxis()->SetBinLabel(1,"FMD1i"); hc->SetBinContent(1,corrs[0]);
  hc->GetXaxis()->SetBinLabel(2,"FMD2i"); hc->SetBinContent(2,corrs[1]);
  hc->GetXaxis()->SetBinLabel(3,"FMD2o"); hc->SetBinContent(3,corrs[2]);
  hc->GetXaxis()->SetBinLabel(4,"FMD3i"); hc->SetBinContent(4,corrs[3]);
  hc->GetXaxis()->SetBinLabel(5,"FMD3o"); hc->SetBinContent(5,corrs[4]);
  hc->GetXaxis()->SetLabelSize(0.08);
  hc->GetYaxis()->SetTitle("#bar{occupancy}");
  hc->SetMarkerSize(1.5);
  hc->Draw("text hist");
  hc->SetMaximum(hc->GetMaximum()*1.5);

  // TH2D* highCuts = static_cast<TH2D*>(dc->FindObject("highCuts"));
  // if (highCuts) highCuts->Draw("colz");
  c->cd();
  c->SaveAs("occupancy.png");
}
//
// EOF
//
