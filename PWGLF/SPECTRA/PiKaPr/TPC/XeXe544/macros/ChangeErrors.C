//
//to include statistical and systematic error in the pT spectra
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

TH1* GetGaus(const Char_t* particleType = "Pion", Int_t cent = 0);
TH1* GetError(const Char_t* particleType = "Pion", Int_t cent = 0);
TH1* CalcSigma(const Char_t* particleType = "Pion", Int_t cent = 0);
TCanvas* DrawError(const Char_t* particleType = "Pion");
void SetColor(TH1* h, Int_t cent = 0);
const Char_t* GetCent(Int_t cent = 0);

void ChangeErrors()
{
  gStyle->SetOptTitle(0);
  //DrawError("Pion");
  //DrawError("AntiPion");
  DrawError("Kaon");
  DrawError("AntiKaon");
  //DrawError("Proton");
  //DrawError("AntiProton");
}
TH1* GetGaus(const Char_t* particleType, Int_t cent)
{
  TVirtualPad* currpad1 = gPad;
  TFile fGaus(Form("FinalSpectra_forQM/canvPtDist%s.root", particleType), "READ");
  Printf("%s",fGaus.GetName());
  TCanvas* cGaus = (TCanvas*)fGaus.Get(Form("canvPtDist%s", particleType));
  Printf("%s",cGaus->GetName());
  TList* lGaus = cGaus->GetListOfPrimitives();
  lGaus->ls();
  TH1* hGaus = (TH1*)lGaus->At(cent)->Clone(Form("%s%i", particleType, cent));  
  hGaus->SetDirectory(0);
  fGaus.Close();
  currpad1->cd();
  return hGaus;
}
TH1* GetError(const Char_t* particleType, Int_t cent)
{
  TVirtualPad* currpad1 = gPad;
  Printf("Getting error for %s in bin %i", particleType, cent);
  TFile fError(Form("FinalSpectra_forQM/RelativeError_%s.root", particleType), "READ");
  Printf("%s", fError.GetName());
  fError.ls();
  TCanvas* cError = (TCanvas*)fError.Get(Form("c%s", particleType));
  Printf("%s", cError->GetName());
  TList* lError = cError->GetListOfPrimitives();
  lError->ls();
  TH1* hError = (TH1*)(lError->FindObject(Form("Ratio%s%i",particleType,cent)))->Clone(Form("%s%i", particleType, cent));
  hError->SetDirectory(0);
  fError.Close();
  currpad1->cd();
  return hError;
}
TH1* CalcSigma(const Char_t* particleType, Int_t cent)
{
  TH1* y = GetGaus(particleType, cent);
  Printf("Spectra %s", y->GetName());
  y->Print("all");
  TH1* relErr = GetError(particleType, cent); 
  Printf("RelativeError %s", relErr->GetName());
  relErr->Print("all");
  for (Int_t bin = 1; bin <= y->GetNbinsX(); bin++){
      if(relErr->GetBinContent(bin)<=0) continue;
  Printf("binErr   %f and bincontent %f",relErr->GetBinContent(bin),y->GetBinContent(bin));
      if(y->GetBinContent(bin)<=0) continue;
      Double_t relError = relErr->GetBinContent(bin); // from ratio of Gauss to fixed cut
      relError = TMath::Sqrt(relError*relError + 0.02*0.02); // add 2% in quadrature as baseline uncertainty (muon contamination etc.)
    y->SetBinError(bin, relError * y->GetBinContent(bin));
  }
  y->SetDirectory(0);
  y->SetName(Form("SpectraWsys_%s%i", particleType, cent));
  y->SetTitle(Form("%s%%;#it{p}_{T} (GeV/#it{c});Relative Error", GetCent(cent)));
  SetColor(y, cent);
  return y;
}
TCanvas* DrawError(const Char_t* particleType)
{
  TCanvas* canvErr = new TCanvas(Form("canvErr%s", particleType), particleType);
  canvErr->DrawFrame(0.1, 0, 1, 140, ";pT;1/N_{ev} #frac{d^{2}N}{dp_{T}d#it{y}}(GeV/#it{c})"); // K
  //canvErr->DrawFrame(0.1, 0, 1, 15, ";pT;1/N_{ev} #frac{d^{2}N}{dp_{T}d#it{y}}(GeV/#it{c})"); // p
  //canvErr->DrawFrame(0.1, 0, 1, 900, "Centrality;pT;1/N_{ev} #frac{d^{2}N}{dp_{T}d#it{y}}(GeV/#it{c})"); // Pi
  
  for (Int_t i = 0; i < 9; i++) {
    TH1 * specSys = CalcSigma(particleType, i);
    TH1 * specStat = GetGaus(particleType, i);
    specSys->SetFillStyle(0);
    specSys->SetLineColor(specSys->GetMarkerColor());
    specSys->SetFillColor(specSys->GetMarkerColor());
    specSys->Draw("E2,same");
    specStat->SetLineColor(specSys->GetMarkerColor());
    specStat->SetMarkerColor(specSys->GetMarkerColor());
    specStat->Draw("same");
  }
    //CalcSigma(particleType, 0)->Draw("E2,same");
  //canvErr->BuildLegend(0.1,0.5,0.3,0.9,"Centrality");
  //canvErr->BuildLegend(0.7,0.5,0.9,0.9,"Centrality");
  TLatex* label = new TLatex(.8, .91, particleType);
  label->SetNDC();
  label->Draw();
  canvErr->SetGridy();
  canvErr->SaveAs(Form("FinalSpectra_forQM/ChangeErr_%s.png", particleType));
  canvErr->SaveAs(Form("FinalSpectra_forQM/ChangeErr_%s.root", particleType));
  return canvErr;
}
void SetColor(TH1* h, Int_t cent)
{
  const TString col[10] = { "#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a" };
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
