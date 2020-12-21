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

void PlotFraction(const Char_t* particleType = "Pion");
TH1D* GetFraction(const Char_t* particleType, Float_t binlow = 0, Float_t binUp = 5);
void SetColor(TH1* h, Int_t cent = 0);

void SetColor(TH1* h, Int_t cent)
{
  const TString col[10] = { "#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a" };
  //http://colorbrewer2.org/#type=qualitative&scheme=Paired&n=10
  h->SetLineColor(TColor::GetColor(col[cent]));
  h->SetMarkerColor(TColor::GetColor(col[cent]));
}

TH1D* GetFraction(const Char_t* particleType, Float_t binlow, Float_t binUp)
{
  TVirtualPad* currpad = gPad;
  Printf("%f %f", binlow, binUp);
  TFile newfile(Form("FeedDownPlots/%s_%.0f_%.0f_fraction.root", particleType, binlow, binUp), "READ");
  Printf("%s", newfile.GetName());
  newfile.ls();
  TCanvas* c = (TCanvas*)newfile.Get(Form("%s_%.0f_%.0f_fraction", particleType, binlow, binUp));
  Printf("%s", c->GetName());
  TList* l = c->GetListOfPrimitives();
  l->ls();
 // TH1D* hFrac = (TH1D*)l->At(2)->Clone(Form("%s_%s", l->At(2)->GetName(), particleType));
  TH1D* hFrac = (TH1D*)l->FindObject(Form("DCASpectrum_hist0to0_pTvsDCA_%sPrim", particleType));
  hFrac->SetMarkerStyle(4);
  hFrac->SetDirectory(0);
  for (Int_t i = 1; i < hFrac->GetNbinsX(); i++)
    hFrac->SetBinError(i, 0);
  newfile.Close();
  currpad->cd();
  return hFrac;
}

void PlotFraction(const Char_t* particleType)
{
  gStyle->SetOptTitle(0);
  const Int_t nCents = 9;
  Double_t centMin[nCents] = { 0.5, 5.5, 10.5, 20.5, 30.5, 40.5, 50.5, 60.5, 70.5 };
  Double_t centMax[nCents] = { 4.5, 9.5, 19.5, 29.5, 39.5, 49.5, 59.5, 69.5, 89.5 };
  TH1D* frac[nCents] = { 0x0 };
  //Bool_t skip[nCents] = { kTRUE };
  //skip[0] = kFALSE;
  TString strType(particleType);
  TCanvas* primFracAllBins = new TCanvas(Form("primFracAllBins_%s", particleType));
  primFracAllBins->DrawFrame(0, 0, 1, 1, Form("centrality;pT;Prim Frac %s", particleType));
  for (Int_t icent = 0; icent < nCents; icent++) {
    //if (skip[icent])
    // continue;
    frac[icent] = GetFraction(strType, centMin[icent] - .5, centMax[icent] + .5);
    primFracAllBins->cd();
    frac[icent]->Draw("EPsame");
    SetColor(frac[icent], icent);
    frac[icent]->SetTitle(Form("%.0f-%.0f", centMin[icent] - .5, centMax[icent] + .5));
  }
  primFracAllBins->BuildLegend(0.7, 0.5, 0.9, 0.9, "");
}
void DrawPrimFrac()
{
 // PlotFraction("Proton");
 // PlotFraction("AntiProton");
    PlotFraction("Pion");
  PlotFraction("AntiPion"); 
}
