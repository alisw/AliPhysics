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

#include "AliTOFTemplateFitter.h"

void SetStyle(TH1* h);
void Draw(TObjArray* arr, TString opt);
void DrawFract(TObjArray* arr, TString opt);
TH2D* GetFeedDownInCentralityBin(TH3D* histPartDCAcent, const Char_t* particleType = "Pion", const Char_t* particleSrc = "", Float_t binlow = 0.5, Float_t binUp = 4.5);
TObjArray* GetCorrection(TH2D* pTvsDCA, TH2D* pTvsDCAPrim = 0x0, TH2D* pTvsDCASecSt = 0x0, TH2D* pTvsDCASecMat = 0x0, const Char_t* particleType = "", Double_t centLow = 0.0, Double_t centHigh = 0.0);
void GetFeedDown(const Char_t* particleType = "Pion", Double_t centLow = 0.5, Double_t centHigh = 4.5);

void PlotFeedDown(const Char_t* particleType)
{
  GetFeedDown(particleType, 0., 5);
  GetFeedDown(particleType, 5, 10);
  GetFeedDown(particleType, 10, 20);
  GetFeedDown(particleType, 20, 30);
  GetFeedDown(particleType, 30, 40);
  GetFeedDown(particleType, 40, 50);
  GetFeedDown(particleType, 50, 60);
  GetFeedDown(particleType, 60, 70);
  GetFeedDown(particleType, 70, 90);
}

void PlotFeedDown()
{
  //PlotFeedDown("Pion");
  //PlotFeedDown("AntiPion");
  PlotFeedDown("Proton");
  //PlotFeedDown("AntiProton");

  //PlotFeedDown("Kaon");
  //PlotFeedDown("AntiKaon");
}

//void GetFeedDown(const Char_t* particleType){
void GetFeedDown(const Char_t* particleType, Double_t centLow, Double_t centHigh)
{

  //
  TFile* inFile = TFile::Open("./Output_2018_03_17data/MergeData.root"); // open analysis results folder
  TList* list = (TList*)inFile->Get("chist");                            // make my output container to store results
  //
  TH1F* fHistCent = (TH1F*)list->FindObject("fHistCent");
  //
  //TCanvas* canvDCA3D = new TCanvas(Form("canvDCA%s3D", particleType), Form("canvDCA%s3D", particleType));
  TH3D* fDCAPart = (TH3D*)list->FindObject(Form("fDCA%s", particleType));
  //fDCAPart->DrawCopy();
  fDCAPart->Sumw2();
  //

  //TH2D* histDCAPart = GetFeedDownInCentralityBin(fDCAPart, particleType, "", 0.5, 4.5);
  TH2D* histDCAPart = GetFeedDownInCentralityBin(fDCAPart, particleType, "", centLow + 0.5, centHigh - 0.5);
  //
  TFile* inFileMC = TFile::Open("./Output_2018_03_17_MC/MergeOutput_MC.root"); // open analysis results folder
  TList* listMC = (TList*)inFileMC->Get("chist");
  //
  TH3D* fDCAPartMC = (TH3D*)listMC->FindObject(Form("fDCA%sMC", particleType));
  TH2D* histDCAMCPrimPart = GetFeedDownInCentralityBin(fDCAPartMC, particleType, "Prim", 0., 0.);
  TH2D* histDCAMCSecStPart = GetFeedDownInCentralityBin(fDCAPartMC, particleType, "SecSt", 1., 1.);
  TH2D* histDCAMCSecMatPart = GetFeedDownInCentralityBin(fDCAPartMC, particleType, "SecMat", 2., 2.);

  if (strcmp(particleType, "AntiProton") == 0)
    histDCAMCSecMatPart = 0x0; //commented for protons
  TObjArray* fit = GetCorrection(histDCAPart, histDCAMCPrimPart, histDCAMCSecStPart, histDCAMCSecMatPart, particleType, centLow, centHigh);
  //TCanvas* canvFeedDown = new TCanvas(Form("canvFeedDown%s", particleType), Form("canvFeedDown%s", particleType));
  TCanvas* canvFeedDown = new TCanvas(Form("canvFeedDown%s_%.0fto%.0f", particleType, centLow, centHigh), Form("canvFeedDown%s_%.0fto%.0f", particleType, centLow, centHigh));
  canvFeedDown->DrawFrame(0, 1, 1, 1e6, ";#it{p}_{T} (GeV/#it{c});Counts");
  gPad->SetLogy();
  Draw(fit, "Epsame");
  TLatex* part = new TLatex(0.78, 0.86, Form("%s", particleType));
  part->SetNDC();
  part->SetTextSize(0.04);
  part->SetTextFont(42);
  part->SetLineWidth(2);
  part->Draw();
  //canvFeedDown->SaveAs(Form("FeedDownPlots/%s_yield.png", particleType));
  canvFeedDown->SaveAs(Form("FeedDownPlots/%s_%.0f_%.0f_yield.png", particleType, centLow, centHigh));

  //TCanvas* canvPrimFract = new TCanvas(Form("canvPrimFract%s", particleType), Form("canvPrimFract%s", particleType));
  TCanvas* canvPrimFract = new TCanvas(Form("%s_%.0f_%.0f_fraction", particleType, centLow, centHigh), Form("%s_%.0f_%.0f_fraction", particleType, centLow, centHigh));
  canvPrimFract->cd();
  canvPrimFract->DrawFrame(0, 0, 1, 1, ";#it{p}_{T} (GeV/#it{c});Fraction");
  DrawFract(fit, "Epsame");
  part->Draw();
  canvPrimFract->SaveAs(Form("FeedDownPlots/%s_%.0f_%.0f_fraction.png", particleType, centLow, centHigh));
  canvPrimFract->SaveAs(Form("FeedDownPlots/%s_%.0f_%.0f_fraction.root", particleType, centLow, centHigh));

  TCanvas* canvprimFracAllBin = new TCanvas("canvprimFracAllBin", "canvprimFracAllBin");
  canvprimFracAllBin->cd();
  //for (Int_t i=centLow; i<=centHigh; i++ )
  //{
  DrawFract(fit, "Epsame");
  // }
}

TH2D* GetFeedDownInCentralityBin(TH3D* histPartDCAcent, const Char_t* particleType, const Char_t* particleSrc, Float_t binlow, Float_t binUp)
{
  histPartDCAcent->GetZaxis()->SetRangeUser(binlow, binUp);
  TH2D* pTvsDCA = static_cast<TH2D*>(histPartDCAcent->Project3D("yx"));
  TString name = Form("hist%ito%i_pTvsDCA_%s%s", TMath::Nint(binlow - 0.5), TMath::Nint(binUp + 0.5), particleType, particleSrc);

  TString title = Form("%s;#it{p}_{T} (GeV/#it{c});DCA_{xy}", name.Data());
  pTvsDCA->SetNameTitle(name, title);
  return pTvsDCA;
}

Double_t FindDCAxyCut(Double_t pt, const Double_t factor = -1.0)
{
  const Double_t dca = 0.0105 + 0.0350 / (TMath::Power(pt, 1.1));
  Printf("DCAxy cut at %.1f GeV/c is %.3f", pt, dca);
  if (factor < 0)
    return dca;
  return dca * factor;
}

TObjArray* GetCorrection(TH2D* pTvsDCA, TH2D* pTvsDCAPrim, TH2D* pTvsDCASecSt, TH2D* pTvsDCASecMat, const Char_t* particleType, Double_t centLow, Double_t centHigh)
{
  //
  // get the primary fraction from the 2D pt vs DCAxy histogram
  //

  //pTvsDCA->RebinY(2); //6 for 40to50, 2 for 50to70,2 for 80 Pions
  //pTvsDCAPrim->RebinY(2);
  //pTvsDCASecSt->RebinY(2);
  //pTvsDCASecMat->RebinY(4);
  TCanvas* canvQADCAFit = new TCanvas(Form("canvQADCAFit_%s", pTvsDCA->GetName()), Form("canvQADCAFit_%s", pTvsDCA->GetName()));
  pTvsDCA->Draw("colz");
  Int_t nsrc = 0;
  TH1D* dcaSpectrumPrim = 0x0;
  if (pTvsDCAPrim) {
    dcaSpectrumPrim = pTvsDCA->ProjectionX(Form("DCASpectrum_%s", pTvsDCAPrim->GetName()));
    SetStyle(dcaSpectrumPrim);
    dcaSpectrumPrim->Reset();
    nsrc++;
  }
  TH1D* dcaSpectrumSecSt = 0x0;
  if (pTvsDCASecSt) {
    dcaSpectrumSecSt = pTvsDCA->ProjectionX(Form("DCASpectrum_%s", pTvsDCASecSt->GetName()));
    SetStyle(dcaSpectrumSecSt);
    dcaSpectrumSecSt->Reset();
    nsrc++;
  }
  TH1D* dcaSpectrumSecMat = 0x0;
  if (pTvsDCASecMat) {
    dcaSpectrumSecMat = pTvsDCA->ProjectionX(Form("DCASpectrum_%s", pTvsDCASecMat->GetName()));
    SetStyle(dcaSpectrumSecMat);
    dcaSpectrumSecMat->Reset();
    nsrc++;
  }
  TObjArray* array = new TObjArray(nsrc);
  if (dcaSpectrumPrim)
    array->Add(dcaSpectrumPrim);
  if (dcaSpectrumSecSt)
    array->Add(dcaSpectrumSecSt);
  if (dcaSpectrumSecMat)
    array->Add(dcaSpectrumSecMat);

  //
  TCanvas* canvQADCAFitSingleBin = new TCanvas("DCAfit", "DCAfit");
  pTvsDCA->Draw("colz");
  //pTvsDCA->SetBit(TH2::kNoTitle);
  pTvsDCA->SetBit(TH2::kNoStats);
  canvQADCAFitSingleBin->SaveAs(Form("FeedDownPlots/%s_%ito%i_fit.pdf(", particleType, TMath::Nint(centLow), TMath::Nint(centHigh)), "pdf");
  gPad->SetLogy();
  for (Int_t iBin = 0; iBin < pTvsDCA->GetXaxis()->GetNbins(); iBin++) {
    //Int_t  iBin = 17;
    //{
    Double_t pT = pTvsDCA->GetXaxis()->GetBinCenter(iBin);
    if (pT < 0.18 || pT > 0.8)
      continue;
    Double_t DCAcut = FindDCAxyCut(pT);
    //
    TH1D* DCAHist = pTvsDCA->ProjectionY(Form("DCA_%s_%f", pTvsDCA->GetName(), pT), iBin, iBin);
    SetStyle(DCAHist);
    TH1D* DCAPrim = 0x0;
    Int_t nsrc = 0;
    if (pTvsDCAPrim) {
      DCAPrim = pTvsDCAPrim->ProjectionY(Form("DCA_%s_%f", pTvsDCAPrim->GetName(), pT), iBin, iBin);
      SetStyle(DCAPrim);
      nsrc++;
    }
    TH1D* DCASecSt = 0x0;
    if (pTvsDCASecSt) {
      DCASecSt = pTvsDCASecSt->ProjectionY(Form("DCA_%s_%f", pTvsDCASecSt->GetName(), pT), iBin, iBin);
      SetStyle(DCASecSt);
      nsrc++;
    }
    TH1D* DCASecMat = 0x0;
    if (pTvsDCASecMat) {
      DCASecMat = pTvsDCASecMat->ProjectionY(Form("DCA_%s_%f", pTvsDCASecMat->GetName(), pT), iBin, iBin);
      SetStyle(DCASecMat);
      nsrc++;
    }
    if (nsrc > 0) {
      TObjArray* mc = new TObjArray(nsrc);
      if (DCAPrim)
        mc->Add(DCAPrim);
      if (DCASecSt)
        mc->Add(DCASecSt);
      if (DCASecMat)
        mc->Add(DCASecMat);
      Double_t range[2] = { -DCAcut, DCAcut };
      Double_t fitrange[2] = { -3, 3 };
      TArrayD fraction(nsrc);
      TArrayD fractionErr(nsrc);
      TObjArray* prediction = new TObjArray(nsrc);
      Bool_t status = PerformFitWithTFF(DCAHist, mc, range, fitrange, fraction, fractionErr, prediction);
      canvQADCAFitSingleBin->Clear();
      canvQADCAFitSingleBin->cd();
      DCAHist->DrawCopy();//->GetXaxis()->SetRangeUser(range[0], range[1]);
      for (Int_t i = 0; i < nsrc; i++) {
        ((TH1D*)array->At(i))->SetBinContent(iBin, fraction[i]);
        ((TH1D*)array->At(i))->SetBinError(iBin, fractionErr[i]);
      }
      //
      for (Int_t i = 0; i < nsrc + 1 && status; i++)
        ((TH1D*)prediction->At(i))->Draw("SAME");
      //
      Draw(prediction, "Epsame");

      //
      TLatex* ptlabel = new TLatex(0.65, 0.7, Form("#it{p}_{T} = %.3f GeV/#it{c}", pT));
      ptlabel->SetTextSize(0.04);
      ptlabel->SetTextFont(42);
      ptlabel->SetLineWidth(2);
      ptlabel->SetNDC();
      ptlabel->Draw();
      /*TLatex* part = new TLatex(0.7, 0.8, Form("for %s", particleType));
      part->SetTextSize(0.04);
      part->SetTextFont(42);
      part->SetLineWidth(2);
      part->SetNDC();
      part->Draw();*/
      TLatex* cent = new TLatex(0.65, 0.8, Form("centrality bin= %.0f-%.0f%%", centLow,centHigh));
      cent->SetTextSize(0.04);
      cent->SetTextFont(42);
      cent->SetLineWidth(2);
      cent->SetNDC();
      cent->Draw();
      canvQADCAFitSingleBin->SaveAs(Form("FeedDownPlots/%s_%ito%i_fit.pdf", particleType, TMath::Nint(centLow), TMath::Nint(centHigh)), "pdf");
    }
  }
  canvQADCAFitSingleBin->DrawFrame(0, 1, 0, 1, "EMPTY;EMPTY;EMPTY");
  canvQADCAFitSingleBin->SaveAs(Form("FeedDownPlots/%s_%ito%i_fit.pdf)", particleType, TMath::Nint(centLow), TMath::Nint(centHigh)), "pdf");
  delete canvQADCAFitSingleBin;
  //
  return array;
}
void SetStyle(TH1* h)
{
  TString c[3] = { "#e41a1c", "#377eb8", "#4daf4a" };
  TString t[3] = { "Prim.", "Sec.St.", "Sec.Mat" };
  TString name = h->GetName();
  Int_t coli = -1;
  if (name.Contains("Prim"))
    coli = 0;
  else if (name.Contains("SecSt"))
    coli = 1;
  else if (name.Contains("SecMat"))
    coli = 2;
  Int_t col = kBlack;
  if (coli >= 0) {
    h->SetTitle(t[coli]);
    col = TColor::GetColor(c[coli]);
  }
  h->SetBit(TH1::kNoTitle);
  h->SetBit(TH1::kNoStats);
  h->SetLineColor(col);
  h->SetLineWidth(2);
  h->SetMarkerColor(col);
}
void Draw(TObjArray* arr, TString opt)
{
  TLegend* leg = new TLegend(0.1, 0.7, 0.35, 0.85);

  for (Int_t i = 0; i < arr->GetEntries(); i++) {

    if (i == 1)
      opt.Append("same");
    arr->At(i)->Draw(opt);
    leg->AddEntry(arr->At(i));
  }
  
  leg->Draw();
}
void DrawFract(TObjArray* arr, TString opt)
{
  TLegend* leg = new TLegend(0.1, 0.7, 0.35, 0.85);
  TH1* hsum = 0x0;
  TH1* harr[10] = { 0x0 };
  Int_t counter = 0;
  for (Int_t i = 0; i < arr->GetEntries(); i++) {
    if (i == 1)
      opt.Append("same");
    harr[counter++] = ((TH1*)arr->At(i))->DrawCopy(opt);
    harr[counter - 1]->SetDirectory(0);
    leg->AddEntry(harr[counter - 1]);
    //
    if (!hsum)
      hsum = (TH1*)harr[counter - 1]->Clone("hsum");
    else
      hsum->Add(harr[counter - 1]);
  }
  for (Int_t i = 0; i < arr->GetEntries(); i++) {
    harr[i]->Divide(hsum);
  }
  leg->Draw();
  delete hsum;
}
