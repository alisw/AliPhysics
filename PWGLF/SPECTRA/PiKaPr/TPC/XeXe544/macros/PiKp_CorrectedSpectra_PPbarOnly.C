//
// raw spectranormalized with (1 / Nev)*  dN/(dy*dpt)
// applied efficiency correction
// positive and negative particles both
// pT range cut applied
// multigaussian fit applied
// feed down correction applied
//
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
#include "DrawPrimFrac.C"


void Cosmetics(TH1D* hist);
void CutPtRange(TH1D* inputHist, const Char_t* particleType = "Pion");
TH1D* GetRawSpectrum(TH2D* pTvsNsigma, const Char_t* particleType = "Pion",Float_t binlow = 0, Float_t binUp = 5);
TH1D* GetSpectrumInCentralityBinMultGaus(TH3D* histPartCent, const Char_t* particleType = "Pion", Float_t binlow = 0.5, Float_t binUp = 4.5);
TH1D* GetSpectrumInCentralityBinFixCut(TH3D* histPartCent, const Char_t* particleType = "Pion", Float_t binlow = 0.5, Float_t binUp = 4.5);
void NormalizeSpectrum(TH1D* spectrum, Float_t dy, Float_t numberOfEvents);
//TH1D* GetFraction(const Char_t* particleType = "Pion", Float_t binlow = 0, Float_t binUp = 5);
TH1D* GetEfficiency(const Char_t* particleType = "Muon");
void CalcProjPerCent(const Char_t* particleType, TH1D* efficiency);

void PiKp_CorrectedSpectra_PPbarOnly()
{

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  //
  TFile* inFile = TFile::Open("./Output_2018_03_17data/MergeData.root"); // open analysis results folder
  TList* list = (TList*)inFile->Get("chist");                            // make my output container to store results
  //
  // QA plot for centrality
  //
  TH1F* fHistCent = (TH1F*)list->FindObject("fHistCent");
  TCanvas* canvCentrality = new TCanvas("canvCentrality", "Centrality");
  fHistCent->DrawCopy("EP");
  inFile->Close();
  //
  // smooth primary fraction with a fit
  //
  /*
  TF1* fitFraction = new TF1("fitFraction", "1-[0]*TMath::Exp(-[1]*x)", 0.45, 0.75);
  fitFraction->SetParameters(1, 1);
  fracAntiProton->Fit(fitFraction, "QNR");
  for (Int_t iBin = 1; iBin < fracAntiProton->GetXaxis()->GetNbins(); iBin++) {
    if (fracAntiProton->GetBinContent(iBin) < 0.4) {
      fracAntiProton->SetBinContent(iBin, fitFraction->Eval(fracAntiProton->GetBinCenter(iBin)));
    }
  } */
  //
  // calculate efficiencies
  //
  //TH1D * effPartPion = GetEfficiency("Pion");
  //TH1D * effPartAntiPion = GetEfficiency("AntiPion");
  //TH1D * effPartKaon = GetEfficiency("Kaon");
  //TH1D * effPartAntiKaon = GetEfficiency("AntiKaon");
  TH1D* effPartProt = GetEfficiency("Prot");
  TH1D* effPartAntiProt = GetEfficiency("AntiProt");
   TF1* antiPrCorr = new TF1("antiPrCorr", "1+[0]*TMath::Exp(-[1]*x)", 0.4,0.85);
  antiPrCorr->SetParameters(0.25993,3.27138);
  for (Int_t iBin = 1; iBin < effPartAntiProt->GetXaxis()->GetNbins(); iBin++) {
    Double_t binContent = effPartAntiProt->GetBinContent(iBin);
    Double_t pbarCorr = antiPrCorr->Eval(effPartAntiProt->GetBinCenter(iBin));
    Double_t newEfficiency  = binContent*pbarCorr ;
    effPartAntiProt->SetBinContent(iBin,newEfficiency );
  }
  TCanvas* newpbarEff = new TCanvas("newpbarEff", "newpbarEff");
  effPartAntiProt->Draw("EP");
  effPartAntiProt->SaveAs("EffwithCorrection.root"); 
  //
  // Generate and correct Raw Spectra
  //
  //CalcProjPerCent("Pion", effPartPion);
  //CalcProjPerCent("AntiPion", effPartAntiPion);
  //CalcProjPerCent("Kaon", effPartKaon);
  //CalcProjPerCent("AntiKaon", effPartAntiKaon);
  CalcProjPerCent("Proton", effPartProt);
  CalcProjPerCent("AntiProton", effPartAntiProt);
}
//******************************************************************************************************************

void CalcProjPerCent(const Char_t * particleType, TH1D * efficiency)
{

  TFile* inFile = TFile::Open("./Output_2018_03_17data/MergeData.root"); // open analysis results folder
  TList* list = (TList*)inFile->Get("chist");                            // make my output container to store results
  //
  // Plot centrality
  //
  TH1F* fHistCent = (TH1F*)list->FindObject("fHistCent");
  //
  // Plot 3D raw spectra histograms
  //
  TH3D* fPartCent = (TH3D*)list->FindObject(Form("f%sDeDxCent", particleType));
  fPartCent->Sumw2();
    /*  
      TF1* fitFraction = new TF1("fitFraction", "1-[0]*TMath::Exp(-[1]*x)", 0.45, 0.75);
  fitFraction->SetParameters(1, 1);
    if (strType.EqualTo("AntiProton")) {
      frac_40to50->Fit(fitFraction, "QNR");
  for (Int_t iBin = 1; iBin < frac_40to50->GetXaxis()->GetNbins(); iBin++) {
    if (frac_40to50->GetBinContent(iBin) < 0.7) {
      frac_40to50->SetBinContent(iBin, fitFraction->Eval(frac_40to50->GetBinCenter(iBin)));
    }
  }
    } */
    const Int_t nCents = 9;
    Double_t centMin[nCents] = { 0.5, 5.5, 10.5, 20.5, 30.5, 40.5, 50.5, 60.5, 70.5 };
    Double_t centMax[nCents] = { 4.5, 9.5, 19.5, 29.5, 39.5, 49.5, 59.5, 69.5, 89.5 };
    const TString col[10] = { "#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a" }; //http://colorbrewer2.org/#type=qualitative&scheme=Paired&n=10
    TH1D* histPart[nCents] = { 0x0 };
    TH1D* frac[nCents] = { 0x0 };
    TString strType(particleType);
    if (strType.EqualTo("AntiPion"))
      strType = "Pion";

    for (Int_t icent = 0; icent < nCents && !strType.Contains("Kaon"); icent++) {
      frac[icent] = GetFraction(strType, centMin[icent] - .5, centMax[icent] + .5);
    }
    //for (Int_t icent = 0; icent < nCents && strType.EqualTo("Proton"); icent++) {
    //frac[icent] = GetFraction(strType, 70, 90);
    //}

    for (Int_t icent = 0; icent < nCents; icent++) {
      // histPart[icent] = GetSpectrumInCentralityBinFixCut(fPartCent,
      //     particleType,
      //     centMin[icent],
      //     centMax[icent]);
      histPart[icent] = GetSpectrumInCentralityBinMultGaus(fPartCent,
        particleType,
          centMin[icent],
          centMax[icent]);
      NormalizeSpectrum(histPart[icent],
          1.0,
          fHistCent->Integral(
              fHistCent->GetXaxis()->FindBin(centMin[icent]),
              fHistCent->GetXaxis()->FindBin(centMax[icent])));
      histPart[icent]->Divide(efficiency);
      if (frac[icent])
       histPart[icent]->Multiply(frac[icent]);
      CutPtRange(histPart[icent], particleType);
      //histPart[icent]->Scale(1e9);
      histPart[icent]->SetLineColor(TColor::GetColor(col[icent]));
      histPart[icent]->SetMarkerColor(TColor::GetColor(col[icent]));
      Cosmetics(histPart[icent]);
      //histPart[icent]->Print("all");
    }
    //
    //Plot the raw spectra
    //
    TCanvas* canvPtDist = new TCanvas(Form("canvPtDist%s", particleType), Form("canvPtDist%s", particleType));
    canvPtDist->SetTitle(Form("%sCentralityClasses", particleType));
    gPad->SetLogy();

    TLegend* legend = new TLegend(0.75, 0.6, 0.9, 0.9, "Xe-Xe,5.44TeV");
    legend->SetHeader(Form("centrality classes"));
    legend->SetTextSize(0.02);
    for (Int_t icent = 0; icent < nCents; icent++) {
      histPart[icent]->GetYaxis()->SetTitle("1/N_{ev} #frac{d^{2}N}{dp_{T}d#it{y}}(GeV^{-1}/c)");
      histPart[icent]->GetXaxis()->SetTitle("p_{T}(GeV/c)");
      histPart[icent]->GetYaxis()->SetRangeUser(1e-1, 100);
      histPart[icent]->GetXaxis()->SetRangeUser(0.2, 1.1);
      histPart[icent]->DrawCopy(icent == 0 ? "Ep" : "Ep,SAME");
      legend->AddEntry(histPart[icent], Form("%.0f-%.0f%%", centMin[icent] - .5, centMax[icent] + .5), "lp");
    }

    legend->Draw();
    canvPtDist->SaveAs(Form("FinalSpectra_forQM/%s.root", canvPtDist->GetName()));
    canvPtDist->SaveAs(Form("FinalSpectra_forQM/%s.png", canvPtDist->GetName()));
    //canvPtDist->SaveAs(Form("Fix_Cut_Spectra/%s.root", canvPtDist->GetName()));
    //canvPtDist->SaveAs(Form("Fix_Cut_Spectra/%s.png", canvPtDist->GetName()));
}

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
  effPart->GetXaxis()->SetRangeUser(0, 1);
  effPart->Divide(histPartReco, histPartGen, 1.0, 1.0, "B");
  TCanvas* canvPartEff = new TCanvas(Form("canv%sEff", particleType), Form("canv%sEff", particleType));
  effPart->DrawCopy("EP");
  canvPartEff->SaveAs(Form("Efficiency_Plots/Eff%s.png", effPart->GetName()));
  canvPartEff->SaveAs(Form("Efficiency_Plots/%s.root", effPart->GetName()));

  effPart->SetDirectory(0);
  delete listMC;
  inFileMC->Close();
  return effPart;
}
/*
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
  TH1D* hFrac = (TH1D*)l->FindObject(Form("DCASpectrum_hist0to0_pTvsDCA_%sPrim", particleType));
  hFrac->SetDirectory(0);
  newfile.Close();
  currpad->cd();
  return hFrac;
}*/

void NormalizeSpectrum(TH1D* spectrum, Float_t dy, Float_t numberOfEvents)
{
  //
  // make (1 / Nev)*  dN/(dy*dpt)
  //
  spectrum->Scale(1. / dy);
  spectrum->Scale(1. / numberOfEvents);
  //
  Int_t nBins = spectrum->GetNbinsX();
  for (Int_t i = 0; i < nBins + 1; i++) {
    Double_t Content = spectrum->GetBinContent(i);
    Double_t error = spectrum->GetBinError(i);
    if (spectrum->GetBinWidth(i) == 0)
      continue;
    spectrum->SetBinContent(i, Content / spectrum->GetBinWidth(i));
    spectrum->SetBinError(i, error / spectrum->GetBinWidth(i));
  }
}

TH1D* GetSpectrumInCentralityBinFixCut(TH3D* histPartCent, const Char_t* particleType, Float_t binlow, Float_t binUp)
{

  TH1D* hist = histPartCent->ProjectionX(Form("hist%ito%i_%s", TMath::Nint(binlow - 0.5), TMath::Nint(binUp + 0.5), particleType), histPartCent->GetYaxis()->FindBin(-3.5), histPartCent->GetYaxis()->FindBin(3.5), histPartCent->GetZaxis()->FindBin(binlow), histPartCent->GetZaxis()->FindBin(binUp));
  return hist;
}

TH1D* GetSpectrumInCentralityBinMultGaus(TH3D* histPartCent, const Char_t* particleType, Float_t binlow, Float_t binUp)
{

  histPartCent->GetZaxis()->SetRangeUser(binlow, binUp);
  TH2D* pTvsNsigma = (TH2D*)histPartCent->Project3D("yx");
  pTvsNsigma->SetNameTitle(Form("hist%ito%i_pTvsNsigma_%s", TMath::Nint(binlow - 0.5), TMath::Nint(binUp + 0.5), particleType), Form("hist%ito%i_pTvsNsigma_%s", TMath::Nint(binlow - 0.5), TMath::Nint(binUp + 0.5), particleType));
  TH1D* hist = GetRawSpectrum(pTvsNsigma,particleType, binlow, binUp);

  return hist;
}

TH1D* GetRawSpectrum(TH2D* pTvsNsigma, const Char_t* particleType, Float_t binlow, Float_t binUp)
{
  //
  // get the raw yield from pT vs Nisgma for a given centrality class
  //
  //TCanvas * canvQAmultGausFit = new TCanvas(Form("canvQAmultGausFit_%s",pTvsNsigma->GetName()),Form("canvQAmultGausFit_%s",pTvsNsigma->GetName()));
  //pTvsNsigma->Draw("colz");
  //canvQAmultGausFit->SaveAs(Form("canvQAmultGausFit_%s.root",pTvsNsigma->GetName()));
  TH1D* rawSpectrum = pTvsNsigma->ProjectionX(Form("rawSpectrum_%s", pTvsNsigma->GetName()));
  rawSpectrum->Reset();
  //
  TCanvas* canvQAmultGausFitSingleBin = new TCanvas(Form("canvQAmultGausFitSingleBin_%s", pTvsNsigma->GetName()), Form("canvQAmultGausFitSingleBin_%s", pTvsNsigma->GetName()));
  canvQAmultGausFitSingleBin->Print(Form("Gaus_Fit/canvQAmultGausFitSingleBin_%s.pdf(", pTvsNsigma->GetName()), "pdf");
  for (Int_t iBin = 0; iBin < pTvsNsigma->GetXaxis()->GetNbins(); iBin++) {
    Double_t pT = pTvsNsigma->GetXaxis()->GetBinCenter(iBin);
    if (pT < 0.1 || pT > 0.8)
      continue;

    TH1D* multGausHist = pTvsNsigma->ProjectionY(Form("multGaus_%s_%f", pTvsNsigma->GetName(), pT), iBin, iBin);
    multGausHist->RebinX(2);
    multGausHist->SetMarkerStyle(20);
    canvQAmultGausFitSingleBin->cd();
    multGausHist->DrawCopy("Ep");

    TF1* funcGausFit = new TF1("funcGausFit", "(x <= ([3] + [1])) * TMath::Abs([0]) * TMath::Gaus(x, [1], [2]) + (x > ([3] + [1])) * [0] * TMath::Gaus([3] + [1], [1], [2]) * TMath::Exp(-([3]) * (x - [3] - [1]) / ([2] * [2]))   +  [4]*TMath::Exp(-x/[5]) ", -5.0, 8.2); //0:norm, 1:mean, 2:sigma, 3:tail
/*
    TString name = multGausHist->GetName();
    if (name.Contains("Pion")) {
      funcGausFit->SetParameters(68231.705, -0.197, 1.199, 1.337, -0.042, 0.8);
      if ((binlow >= 40) && (binUp < 50)) funcGausFit->SetParameters(30165.682, -0.220, 1.073, 1.436, -0.022, 0.800);
      if ((binlow >= 60) && (binUp < 70)) funcGausFit->SetParameters(8958.026, -0.109, 1.064, 1.362, -0.007, 0.800);
      else if ((binlow >= 70) && (binUp < 90)) funcGausFit->SetParameters(10053.502, -0.101, 1.067, 1.492, -0.008, 0.800);

      if (pT < 0.7)
        funcGausFit->SetParLimits(3, 1.4, 3.0);
      funcGausFit->SetParLimits(5, 0, 0.8);

      //if (pT < 0.6) funcGausFit->SetRange(-8.0,8.0);
      //if (pT < 0.2) funcGausFit->SetParLimits(5, -10.0, -2.0);
    }*/

    //if (name.Contains("Prot")) {
      funcGausFit->SetParameters(1480.234, -0.121, 1.095, 1.773, 0.0007, 0.723);
      if (pT < 0.6)
        funcGausFit->SetRange(-8.0, 8.0);

      funcGausFit->SetParLimits(1, -0.3, 0.3);
      funcGausFit->SetParLimits(2, 0.7, 1.3);
      funcGausFit->SetParLimits(3, 0, 3.0);
      funcGausFit->SetParLimits(5, 0, 0.8);
    //}

    //funcGausFit->SetParLimits(7, 0, 5.0);
    multGausHist->Fit(funcGausFit, "QNRL");
    funcGausFit->Draw("SAME");
    Printf("Parameters: pT = %4.3f; par = %4.3f,%4.3f, %4.3f,%4.3f, %4.3f,%4.3f", pT, funcGausFit->GetParameter(0), funcGausFit->GetParameter(1), funcGausFit->GetParameter(2), funcGausFit->GetParameter(3), funcGausFit->GetParameter(4), funcGausFit->GetParameter(5));
    //
    TF1* funcSignal = new TF1("funcSignal", "(x <= ([3] + [1])) * [0] * TMath::Gaus(x, [1], [2]) + (x > ([3] + [1])) * [0] * TMath::Gaus([3] + [1], [1], [2]) * TMath::Exp(-([3]) * (x - [3] - [1]) / ([2] * [2]))", -10.0, 10.0);
    funcSignal->SetParameters(funcGausFit->GetParameter(0), funcGausFit->GetParameter(1), funcGausFit->GetParameter(2), funcGausFit->GetParameter(3));
    for (int i = 0; i < 4; i++)
      funcSignal->SetParError(i, funcGausFit->GetParError(i));
    funcSignal->SetLineColor(kBlue);
    funcSignal->Draw("SAME");
    Double_t yieldSignal = funcSignal->Integral(-10.0, 10.0) / multGausHist->GetBinWidth(10); // random
   //Double_t yieldSignalError = funcSignal->IntegralError(-10.0, 10.0) / multGausHist->GetBinWidth(10);

    TF1* funcBackGround = new TF1("funcBackGround", " [0]*TMath::Exp(-x/[1])", -10.0, 10.0);
    funcBackGround->SetParameters(funcGausFit->GetParameter(4), funcGausFit->GetParameter(5));
    for (int i = 0; i < 2; i++)
      funcBackGround->SetParError(i + 4, funcGausFit->GetParError(i + 4));
    funcBackGround->SetLineColor(kOrange);
    funcBackGround->Draw("SAME");

    TLatex* ptlabel = new TLatex(0.7, 0.8, Form("pT = %f", pT));
    ptlabel->SetNDC();
    ptlabel->SetTextSize(0.04);
    ptlabel->SetTextFont(42);
    ptlabel->SetLineWidth(2);
    ptlabel->Draw();

    TLatex* params = new TLatex(0.01, 0.7, Form("%4.3f,%4.3f,%4.3f,%4.3f,%4.3f,%4.3f", funcGausFit->GetParameter(0), funcGausFit->GetParameter(1), funcGausFit->GetParameter(2), funcGausFit->GetParameter(3), funcGausFit->GetParameter(4), funcGausFit->GetParameter(5)));
    params->SetNDC();
    params->SetTextSize(0.04);
    params->SetTextFont(42);
    params->SetLineWidth(2);
    params->Draw();

    //canvQAmultGausFitSingleBin->SaveAs(Form("canvQAmultGausFitSingleBin_%s_%f.png",pTvsNsigma->GetName(),pT));
    Double_t yieldBck = funcBackGround->Integral(-4.0, 4.0) / multGausHist->GetBinWidth(4); // random
    //Double_t yieldBckError = funcBackGround->IntegralError(-4.0, 4.0) / multGausHist->GetBinWidth(4);
    //
    Double_t yieldError=0;
    Double_t yield = multGausHist->IntegralAndError(multGausHist->GetXaxis()->FindBin(-4.0), multGausHist->GetXaxis()->FindBin(+4.0),yieldError) ;
    yield -= yieldBck;
    //yieldError= TMath::Sqrt(yieldError*yieldError+yieldBckError*yieldBckError);
    yieldError= TMath::Sqrt(yieldError + 2*yieldBck);
    //
    //rawSpectrum->SetBinContent(iBin, yield);
    rawSpectrum->SetBinContent(iBin, yieldSignal);
    rawSpectrum->SetBinError(iBin, yieldError);
    canvQAmultGausFitSingleBin->SetLogy();
    canvQAmultGausFitSingleBin->Print(Form("Gaus_Fit/canvQAmultGausFitSingleBin_%s.pdf", pTvsNsigma->GetName()), "pdf");

    //canvQAmultGausFitSingleBin->SaveAs(Form("Gaus_Fit/canvQAmultGausFitSingleBin_%s_%f.pdf", pTvsNsigma->GetName(), pT), "pdf");
  }
  canvQAmultGausFitSingleBin->DrawFrame(0, 0, 1, 1, "EMPTY;EMPTY;EMPTY");
  canvQAmultGausFitSingleBin->Print(Form("Gaus_Fit/canvQAmultGausFitSingleBin_%s.pdf)", pTvsNsigma->GetName()), "pdf");
  //delete canvQAmultGausFitSingleBin;
  //
  return rawSpectrum;
  //return hist;
}
void CutPtRange(TH1D* inputHist, const Char_t* particleType)
{
  //
  // remove parts were PID does not work
  //
  TH1D* oldHist = (TH1D*)inputHist->Clone("oldHist");
  inputHist->Reset();
  TString particleString(particleType);
  Double_t pTmin = 0.0;
  Double_t pTmax = 10.0;
  if (particleString.Contains("Pion")) {
    pTmin = 0.1;
    pTmax = 0.75;
  }
  if (particleString.Contains("Kaon")) {
    pTmin = 0.2;
    pTmax = 0.55;
  }
  if (particleString.Contains("Proton")) {
    pTmin = 0.4;
    pTmax = 0.85;
  }
  //
  for (Int_t iBin = 0; iBin < inputHist->GetXaxis()->GetNbins(); iBin++) {
    Double_t pT = inputHist->GetXaxis()->GetBinCenter(iBin);
    if (pT > pTmax || pT < pTmin)
      continue;
    inputHist->SetBinContent(iBin, oldHist->GetBinContent(iBin));
    inputHist->SetBinError(iBin, oldHist->GetBinError(iBin));
  }
  delete oldHist;
}

void Cosmetics(TH1D* hist)
{

  hist->SetMarkerStyle(20);
  hist->SetMarkerSize(0.4);
  // hist->GetXaxis()->SetRangeUser(0,1.5);
}
