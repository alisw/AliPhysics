#include "TF2.h"
#include "TH2F.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include "TFile.h"
#include "TLatex.h"
#include "TROOT.h"
#include "TLegend.h"

void AddPoint(TGraphErrors* graph, Float_t x, Float_t y, Float_t xe, Float_t ye)
{
	graph->SetPoint(graph->GetN(), x, y);
	graph->SetPointError(graph->GetN() - 1, xe, ye);
}

void DrawLatex(Float_t x, Float_t y, Int_t color, const char* text)
{
  TLatex* latex = new TLatex(x, y, text);
  latex->SetNDC();
  latex->SetTextSize(0.06);
  latex->SetTextColor(color);
  latex->Draw();
}

/* 
// old one with 1 Gaussian
Double_t DeltaPhiWidth2DFitFunction(Double_t *x, Double_t *par)
{
  // params: 0: gaussian amplitude, 1: phi width, 2: eta width
  //         3..bins+2 constants as fct of phi
  
  Float_t phi = x[0];
  if (phi < 0)
    phi = -phi;
  if (phi > TMath::Pi())
    phi = TMath::TwoPi() - phi;
  Int_t phiBin = (Int_t) (phi / TMath::Pi() * 36);
//   phiBin = 0;
  
  return par[3+phiBin]+par[0]*TMath::Exp(-0.5*((x[0]/par[1])*(x[0]/par[1])+(x[1]/par[2])*(x[1]/par[2])));
}*/

Double_t DeltaPhiWidth2DFitFunction(Double_t *x, Double_t *par)
{
  // params: 0: gaussian amplitude, 1: phi width I, 2: eta width I
  //         3: fraction for first Gaussian, 4: phi width II, 5: eta width II
  //         6..bins+5 constants as fct of phi
  
  Float_t phi = x[0];
  if (phi < 0)
    phi = -phi;
  if (phi > TMath::Pi())
    phi = TMath::TwoPi() - phi;
  Int_t phiBin = (Int_t) (phi / TMath::Pi() * 36);
//   phiBin = 0;
  
  return par[6+phiBin]+par[0]*(par[3]*TMath::Exp(-0.5*((x[0]/par[1])*(x[0]/par[1])+(x[1]/par[2])*(x[1]/par[2]))) + (1-par[3]) * TMath::Exp(-0.5*((x[0]/par[4])*(x[0]/par[4])+(x[1]/par[5])*(x[1]/par[5]))));
}

void SubtractEtaGap(TH2* hist, Float_t etaLimit, Float_t outerLimit, Bool_t scale)
{
  TString histName(hist->GetName());

  TH1D* etaGap = hist->ProjectionX(histName + "_1", TMath::Max(1, hist->GetYaxis()->FindBin(-outerLimit + 0.01)), hist->GetYaxis()->FindBin(-etaLimit - 0.01));
  Int_t etaBins = hist->GetYaxis()->FindBin(-etaLimit - 0.01) - TMath::Max(1, hist->GetYaxis()->FindBin(-outerLimit + 0.01)) + 1;

  TH1D* tracksTmp = hist->ProjectionX(histName + "_2", hist->GetYaxis()->FindBin(etaLimit + 0.01), TMath::Min(hist->GetYaxis()->GetNbins(), hist->GetYaxis()->FindBin(outerLimit - 0.01)));
  etaBins += TMath::Min(hist->GetYaxis()->GetNbins(), hist->GetYaxis()->FindBin(outerLimit - 0.01)) - hist->GetYaxis()->FindBin(etaLimit + 0.01) + 1;
  
  etaGap->Add(tracksTmp);

  // get per bin result
  etaGap->Scale(1.0 / etaBins);
  
//   new TCanvas; etaGap->DrawCopy();
  
  TH2* histTmp2D = (TH2*) hist->Clone("histTmp2D");
  histTmp2D->Reset();
  
  for (Int_t xbin=1; xbin<=histTmp2D->GetNbinsX(); xbin++)
    for (Int_t y=1; y<=histTmp2D->GetNbinsY(); y++)
      histTmp2D->SetBinContent(xbin, y, etaGap->GetBinContent(xbin));
    
  if (scale)
  {
    // mixed event does not reproduce away-side perfectly
    // --> extract scaling factor on the away-side from ratios of eta gap and central region
    TH1D* centralRegion = hist->ProjectionX(histName + "_3", hist->GetYaxis()->FindBin(-etaLimit + 0.01), hist->GetYaxis()->FindBin(etaLimit - 0.01));
    etaBins = hist->GetYaxis()->FindBin(etaLimit - 0.01) - hist->GetYaxis()->FindBin(-etaLimit + 0.01) + 1;
    centralRegion->Scale(1.0 / etaBins);
    
//     new TCanvas; centralRegion->DrawCopy(); etaGap->SetLineColor(2); etaGap->DrawCopy("SAME");
    centralRegion->Divide(etaGap);
//     new TCanvas; centralRegion->Draw();
    centralRegion->Fit("pol0", "0", "", TMath::Pi() - 1, TMath::Pi() + 1);
    Float_t scalingFactor = centralRegion->GetFunction("pol0")->GetParameter(0);
    Printf("  scalingFactor = %f", scalingFactor);
    histTmp2D->Scale(scalingFactor);
  }
    
//   new TCanvas; hist->DrawCopy("SURF1");

  hist->Add(histTmp2D, -1);  
}

void FitDeltaPhiEtaGap2D(TH2* hist, Bool_t scale,  TCanvas* canvas, Int_t canvasPos, TGraphErrors** width, Float_t x, Float_t yPosChi2, TGraphErrors* chi2_1, TGraphErrors* chi2_2)
{
  Float_t etaLimit = 1.0;
  Float_t outerLimit = 1.5;
  
  SubtractEtaGap(hist, etaLimit, outerLimit, scale);

//   new TCanvas; hist->DrawCopy("SURF1");

  hist->GetYaxis()->SetRangeUser(-outerLimit+0.01, outerLimit-0.01);
  hist->GetXaxis()->SetRangeUser(-TMath::Pi() / 2 + 0.01, TMath::Pi() * 0.5 - 0.01);
  
  canvas->cd(canvasPos);
  hist->SetStats(0);
  hist->DrawCopy("SURF1");
  
  Float_t min = hist->GetMinimum();
  Float_t max = hist->GetMaximum();
  
  // ranges are to exclude eta gap region from fit
  TF2* func = new TF2("func", "[0]+[1]*exp(-0.5*((x/[2])**2+(y/[3])**2))", -0.5 * TMath::Pi(), 0.5 * TMath::Pi(), -outerLimit, outerLimit);
  func->SetParameters(0, 1, 0.3, 0.3);
  func->SetParLimits(1, 0, 10);
  func->SetParLimits(2, 0.05, 1);
  func->SetParLimits(3, 0.05, 1);
  func->FixParameter(0, 0);

  /*
  TF2* func = new TF2("func", DeltaPhiWidth2DFitFunction, -TMath::Pi() / 2, TMath::Pi() * 1.5, -1, 1, 4);
  func->SetParameters(1, 0.3, 0.3, 0);
  func->SetParLimits(0, 0, 10);
  func->SetParLimits(1, 0.1, 10);
  func->SetParLimits(2, 0.1, 10);
  */
  
  hist->Fit(func, "0R", "");
  hist->Fit(func, "I0R", "");

//   func->SetRange(-0.5 * TMath::Pi(), 0.5 * TMath::Pi(), -outerLimit, outerLimit);
  
  canvas->cd(canvasPos + 1);
  TH2* funcHist = (TH2*) hist->Clone("funcHist");
  funcHist->Reset();
  funcHist->Add(func);
  funcHist->SetMinimum(min);
  funcHist->SetMaximum(max);
  funcHist->Draw("SURF1");
  
  canvas->cd(canvasPos + 2);
  hist->Add(func, -1);
  hist->SetMinimum(min);
  hist->SetMaximum(max);
  hist->DrawCopy("SURF1");
  
  width[0]->SetPoint(width[0]->GetN(), x, TMath::Abs(func->GetParameter(2)));
  width[0]->SetPointError(width[0]->GetN()-1, 0, func->GetParError(2));
    
  width[1]->SetPoint(width[1]->GetN(), x, TMath::Abs(func->GetParameter(3)));
  width[1]->SetPointError(width[1]->GetN()-1, 0, func->GetParError(3));

  Float_t chi2 = 0;
  Int_t ndf = 0;
  for (Int_t i=hist->GetXaxis()->FindBin(-0.8); i<=hist->GetXaxis()->FindBin(0.8); i++)
    for (Int_t j=hist->GetYaxis()->FindBin(-0.8); j<=hist->GetYaxis()->FindBin(0.8); j++)
    {
      if (hist->GetBinError(i, j) > 0)
      {
	chi2 += TMath::Power(hist->GetBinContent(i, j) / hist->GetBinError(i, j), 2);
	ndf++;
      }
    }
  ndf -= func->GetNumberFreeParameters();
  
  printf("#chi^{2}/ndf = %.1f/%d = %.1f  ", func->GetChisquare(), func->GetNDF(), func->GetChisquare() / func->GetNDF());
  Printf("#chi^{2}/ndf = %.1f/%d = %.1f", chi2, ndf, chi2 / ndf);

  DrawLatex(0.5, yPosChi2, 1, Form("#chi^{2}/ndf = %.1f/%d = %.1f", func->GetChisquare(), func->GetNDF(), func->GetChisquare() / func->GetNDF()));
  DrawLatex(0.5, yPosChi2 - 0.05, 1, Form("#chi^{2}/ndf = %.1f/%d = %.1f", chi2, ndf, chi2 / ndf));

  chi2_1->SetPoint(chi2_1->GetN(), x, func->GetChisquare() / func->GetNDF());
  chi2_2->SetPoint(chi2_2->GetN(), x, chi2 / ndf);
}

void FitDeltaPhi2DOneFunction(TH2* hist, TCanvas* canvas, Int_t canvasPos, TGraphErrors** width, Float_t x, Float_t yPosChi2, TGraphErrors* chi2_1, TGraphErrors* chi2_2, Bool_t quick = kFALSE, Int_t histId = 0, Int_t limits = 0)
{
  Float_t etaLimit = 1.0;
  Float_t outerLimit = 1.79;

  Float_t mean = hist->Integral() / hist->GetNbinsX() / hist->GetNbinsY();
  hist->GetYaxis()->SetRangeUser(-outerLimit+0.01, outerLimit-0.01);
  hist->GetXaxis()->SetRangeUser(-TMath::Pi() / 2 + 0.01, TMath::Pi() * 0.5 - 0.01);
  
  // sums
  TGraphErrors* sumSummary = new TGraphErrors;
  TGraphErrors* phiWidthSummary = new TGraphErrors;
  TGraphErrors* etaWidthSummary = new TGraphErrors;
  
  canvas->cd(canvasPos++);
  hist->SetStats(0);
  hist->DrawCopy("SURF1");
  sumSummary->SetPoint(sumSummary->GetN(), sumSummary->GetN(), hist->Integral());
  
  Float_t min = hist->GetMinimum();
  Float_t max = hist->GetMaximum();

  TH2* subtractFlow = 0;
  if (!quick)
  {
    Int_t bins = hist->GetNbinsX() / 2 / 2;
    
    TF2* func = new TF2("func", DeltaPhiWidth2DFitFunction, -TMath::Pi() / 2, TMath::Pi() * 0.5, -outerLimit, outerLimit, bins+6);
    func->SetParameters(1, 0.3, 0.3, 0.25, 0.6, 0.6);
    for (Int_t i=6; i<bins+6; i++)
      func->SetParameter(i, mean);

    func->SetParLimits(0, 0, 10);
    func->SetParLimits(1, 0.05, 1);
    func->SetParLimits(2, 0.05, 1);
    func->SetParLimits(3, 0.1, 0.9);
    func->SetParLimits(4, 0.05, 1);
    func->SetParLimits(5, 0.05, 1);

    hist->Fit(func, "0R", "");
    //   hist->Fit(func, "0RI", "");

    //   width[0]->SetPoint(width[0]->GetN(), x, TMath::Abs(func->GetParameter(1)));
    //   width[0]->SetPointError(width[0]->GetN()-1, 0, func->GetParError(1));
    //     
    //   width[1]->SetPoint(width[1]->GetN(), x, TMath::Abs(func->GetParameter(2)));
    //   width[1]->SetPointError(width[1]->GetN()-1, 0, func->GetParError(2));

    AddPoint(phiWidthSummary, phiWidthSummary->GetN(), func->GetParameter(1), 0, func->GetParError(1));
    AddPoint(phiWidthSummary, phiWidthSummary->GetN(), func->GetParameter(4), 0, func->GetParError(4));
    
    AddPoint(etaWidthSummary, etaWidthSummary->GetN(), func->GetParameter(2), 0, func->GetParError(2));
    AddPoint(etaWidthSummary, etaWidthSummary->GetN(), func->GetParameter(5), 0, func->GetParError(5));

    canvas->cd(canvasPos++);
    TH2* funcHist = (TH2*) hist->Clone("funcHist");
    funcHist->Reset();
    funcHist->Add(func);
    funcHist->SetMinimum(min);
    funcHist->SetMaximum(max);
    funcHist->Draw("SURF1");
    sumSummary->SetPoint(sumSummary->GetN(), sumSummary->GetN(), funcHist->Integral());
    
    canvas->cd(canvasPos++);
    TH2* residuals = (TH2*) hist->Clone("residuals");
    residuals->Add(func, -1);
    residuals->SetMinimum(-(max - min) / 2);
    residuals->SetMaximum((max - min) / 2);
    residuals->Draw("SURF1");

    //   Float_t chi2 = 0;
    //   Int_t ndf = 0;
    //   for (Int_t i=hist->GetXaxis()->FindBin(-0.8); i<=hist->GetXaxis()->FindBin(0.8); i++)
    //     for (Int_t j=hist->GetYaxis()->FindBin(-0.8); j<=hist->GetYaxis()->FindBin(0.8); j++)
    //     {
    //       if (residuals->GetBinError(i, j) > 0)
    //       {
    // 	chi2 += TMath::Power(residuals->GetBinContent(i, j) / residuals->GetBinError(i, j), 2);
    // 	ndf++;
    //       }
    //     }
    //   ndf -= func->GetNumberFreeParameters();
    //   
    //   printf("#chi^{2}/ndf = %.1f/%d = %.1f  ", func->GetChisquare(), func->GetNDF(), func->GetChisquare() / func->GetNDF());
    //   Printf("#chi^{2}/ndf = %.1f/%d = %.1f", chi2, ndf, chi2 / ndf);
    // 
    //   DrawLatex(0.5, yPosChi2, 1, Form("#chi^{2}/ndf = %.1f/%d = %.1f", func->GetChisquare(), func->GetNDF(), func->GetChisquare() / func->GetNDF()));
    //   DrawLatex(0.5, yPosChi2 - 0.05, 1, Form("#chi^{2}/ndf = %.1f/%d = %.1f", chi2, ndf, chi2 / ndf));
    // 
    //   chi2_1->SetPoint(chi2_1->GetN(), x, func->GetChisquare() / func->GetNDF());
    //   chi2_2->SetPoint(chi2_2->GetN(), x, chi2 / ndf);

    // eta gap subtraction
    canvas->cd(canvasPos++);
    func->SetParameter(0, 0);
    subtractFlow = (TH2*) hist->Clone("subtractFlow");
    subtractFlow->Add(func, -1);
    subtractFlow->SetMinimum(-(max - min) / 2);
    subtractFlow->SetMaximum((max - min) / 2);
    subtractFlow->DrawCopy("SURF1");
    sumSummary->SetPoint(sumSummary->GetN(), sumSummary->GetN(), subtractFlow->Integral());
  }

  canvas->cd(canvasPos++);
  SubtractEtaGap(hist, etaLimit, outerLimit, kFALSE);
  hist->SetMinimum(-(max - min) / 2);
  hist->SetMaximum((max - min) / 2);
  hist->DrawCopy("SURF1");
  sumSummary->SetPoint(sumSummary->GetN(), sumSummary->GetN(), hist->Integral());
    
  if (!quick)
  {
    canvas->cd(canvasPos++);
    TH2* difference = (TH2*) hist->Clone("difference");
    difference->Add(subtractFlow, -1);
    difference->SetMinimum(-(max - min) / 2);
    difference->SetMaximum((max - min) / 2);
    difference->DrawCopy("SURF1");

    canvas->cd(canvasPos++);
    TF2* func2 = new TF2("func2", "[0]+[1]*exp(-0.5*((x/[2])**2+(y/[3])**2))", -0.5 * TMath::Pi(), 0.5 * TMath::Pi(), -outerLimit, outerLimit);
    func2->SetParameters(0, 1, 0.3, 0.3);
    func2->SetParLimits(1, 0, 10);
    func2->SetParLimits(2, 0.05, 1);
    func2->SetParLimits(3, 0.05, 1);
    func2->FixParameter(0, 0);
    
    hist->Fit(func2, "0R", "");
    //   hist->Fit(func2, "I0R", "");
    
    AddPoint(phiWidthSummary, phiWidthSummary->GetN(), func2->GetParameter(2), 0, func2->GetParError(2));
    AddPoint(etaWidthSummary, etaWidthSummary->GetN(), func2->GetParameter(3), 0, func2->GetParError(3));
  
    TH2* funcHist2 = (TH2*) hist->Clone("funcHist2");
    funcHist2->Reset();
    funcHist2->Add(func2);
    funcHist2->SetMinimum(-(max - min) / 2);
    funcHist2->SetMaximum((max - min) / 2);
    funcHist2->Draw("SURF1");
    sumSummary->SetPoint(sumSummary->GetN(), sumSummary->GetN(), funcHist2->Integral());
    
    canvas->cd(canvasPos++);
    TH2* residuals2 = (TH2*) hist->Clone("residuals");
    residuals2->Add(funcHist2, -1);
    residuals2->SetMinimum(-(max - min) / 2);
    residuals2->SetMaximum((max - min) / 2);
    residuals2->Draw("SURF1");
  }
  
  Float_t momentLimit = 0.8;
  TH1* projx2 = hist->ProjectionX(Form("%s_projx2", hist->GetName()), hist->GetYaxis()->FindBin(-momentLimit+0.01), hist->GetYaxis()->FindBin(momentLimit-0.01));
  projx2->GetXaxis()->SetRangeUser(-1, 1);
  
  TH1* projy2 = hist->ProjectionY(Form("%s_projy2", hist->GetName()), hist->GetXaxis()->FindBin(-momentLimit+0.01), hist->GetXaxis()->FindBin(momentLimit-0.01));
  projy2->GetXaxis()->SetRangeUser(-1, 1);
  
  // calculate moments
  for (Int_t axis=0; axis<2; axis++)
  {
    TH1* proj = (axis == 0) ? projx2 : projy2;
    
    for (Int_t n=2; n <= 4; n++)
    {
      Float_t moment = 0;
      Float_t sum = 0;
      for (Int_t bin = proj->GetXaxis()->FindBin(-momentLimit+0.01); bin <= proj->GetXaxis()->FindBin(momentLimit-0.01); bin++)
      {
	moment += proj->GetBinContent(bin) * TMath::Power(proj->GetMean() - proj->GetXaxis()->GetBinCenter(bin), n);
	sum += proj->GetBinContent(bin);
      }
      moment /= sum;
      
      Float_t error = 0;
      for (Int_t bin = proj->GetXaxis()->FindBin(-momentLimit+0.01); bin <= proj->GetXaxis()->FindBin(momentLimit-0.01); bin++)
      {
	error += proj->GetBinError(bin) * proj->GetBinError(bin) * 
	  TMath::Power(TMath::Power(proj->GetMean() - proj->GetXaxis()->GetBinCenter(bin), n) / sum 
	    - moment / sum, 2);
      }
      
      AddPoint(width[6+(n-2)*2 + axis], x, moment, 0, TMath::Sqrt(error));
  //     Printf("%d %f +- %f <-> %f +- %f", n, moment, error, projx2->GetRMS() * projx2->GetRMS(), 2 * projx2->GetRMSError() * projx2->GetRMSError());
    }
  }
  
  projx2->GetXaxis()->SetRangeUser(-momentLimit, momentLimit);
  AddPoint(width[12], x, projx2->GetKurtosis(1), 0, projx2->GetKurtosis(11));

  projy2->GetXaxis()->SetRangeUser(-momentLimit, momentLimit);
  AddPoint(width[13], x, projy2->GetKurtosis(1), 0, projy2->GetKurtosis(11));

//   return;
  
  if (!quick)
  {
    canvas->cd(canvasPos++);
    TH1* projx1 = subtractFlow->ProjectionX(Form("%s_projx1", hist->GetName()), hist->GetYaxis()->FindBin(-0.8), hist->GetYaxis()->FindBin(0.8));
    projx1->Draw();
    projx1->Fit("gaus", "I");
    
  //   TF1* twoGauss = new TF1("twoGauss", "gaus(0)+gaus(3)", -2, 2);
  //   twoGauss->SetParameters(1, 0, 0.3, 1, 0, 0.6);
  //   twoGauss->FixParameter(1, 0);
  //   twoGauss->FixParameter(4, 0);
  //   twoGauss->SetLineColor(4);
  //   projx1->Fit("twoGauss", "I+", "SAME");
    
    projx2->SetLineColor(2);
    projx2->Draw("SAME");
    projx2->Fit("gaus", "I+", "SAME");
    projx2->GetFunction("gaus")->SetLineColor(2);

    canvas->cd(canvasPos++);
    TH1* projy1 = subtractFlow->ProjectionY(Form("%s_projy1", hist->GetName()), hist->GetXaxis()->FindBin(-0.8), hist->GetXaxis()->FindBin(0.8));
    projy1->Draw();
    projy1->Fit("gaus", "I");
    
    projy2->SetLineColor(2);
    projy2->Draw("SAME");
    projy2->Fit("gaus", "I+", "SAME");
    projy2->GetFunction("gaus")->SetLineColor(2);
    
    // 1d fit (lots of params)
    AddPoint(phiWidthSummary, phiWidthSummary->GetN(), projx2->GetFunction("gaus")->GetParameter(2), 0, projx2->GetFunction("gaus")->GetParError(2));
    
    // 1d fit (eta gap subtraction)
    AddPoint(phiWidthSummary, phiWidthSummary->GetN(), projx1->GetFunction("gaus")->GetParameter(2), 0, projx1->GetFunction("gaus")->GetParError(2));

    // 1d fit (lots of params)
    AddPoint(etaWidthSummary, etaWidthSummary->GetN(), projy2->GetFunction("gaus")->GetParameter(2), 0, projy2->GetFunction("gaus")->GetParError(2));
    
    // 1d fit (eta gap subtraction)
    AddPoint(etaWidthSummary, etaWidthSummary->GetN(), projy1->GetFunction("gaus")->GetParameter(2), 0, projy1->GetFunction("gaus")->GetParError(2));
  }

  // 2d fit with two gaussians
  canvas->cd(canvasPos++);
//   TF2* func3 = new TF2("func3", "[0]+[1]*exp(-0.5*((x/[2])**2+(y/[3])**2))+[4]*exp(-0.5*((x/[5])**2+(y/[6])**2))", -0.5 * TMath::Pi(), 0.5 * TMath::Pi(), -outerLimit, outerLimit);
//   func3->SetParameters(0, 1, 0.3, 0.3, 1, 0.6, 0.6);
//   func3->SetParLimits(4, 0, 10);
  Float_t etaFitLimit = outerLimit;
//   Float_t etaFitLimit = 0.5;
  TF2* func3 = new TF2("func3", "[0]+[1]*([4]*exp(-0.5*((x/[2])**2+(y/[3])**2))+(1-[4])*exp(-0.5*((x/[5])**2+(y/[6])**2)))", -0.5 * TMath::Pi(), 0.5 * TMath::Pi(), -etaFitLimit, etaFitLimit);
  func3->SetParameters(0, 1, 0.3, 0.3, 0.25, 0.6, 0.6);
  func3->SetParLimits(4, 0.1, 0.9);

  func3->SetParLimits(1, 0, 10);
  func3->SetParLimits(2, 0.05, 1);
  func3->SetParLimits(3, 0.05, 1);
  func3->SetParLimits(5, 0.05, 1);
  func3->SetParLimits(6, 0.05, 1);
  func3->FixParameter(0, 0);

  if (0 && histId == 0)
  {
    // central --> go to 1 Gaussian only
    func3->SetParLimits(4, 1, 1);
    func3->FixParameter(4, 1);
    func3->FixParameter(5, 0.05);
    func3->FixParameter(6, 0.05);
  }
  
  // set errors very large for bins potentially affected by two-track effects
  hist->SetBinError(hist->GetXaxis()->FindBin(0.0001), hist->GetYaxis()->FindBin(0.0001), 1e5);
  hist->SetBinError(hist->GetXaxis()->FindBin(0.0001), hist->GetYaxis()->FindBin(-0.0001), 1e5);
  hist->SetBinError(hist->GetXaxis()->FindBin(-0.0001), hist->GetYaxis()->FindBin(0.0001), 1e5);
  hist->SetBinError(hist->GetXaxis()->FindBin(-0.0001), hist->GetYaxis()->FindBin(-0.0001), 1e5);

  hist->Fit(func3, "0R", "");
//   hist->Fit(func3, "I0R", "");

  Int_t first = 2;
  Int_t second = 5;
  if (func3->GetParameter(2) < func3->GetParameter(5))
  {
    first = 5;
    second = 2;
  }
  //dphi
  AddPoint(width[1], x, TMath::Abs(func3->GetParameter(first)), 0, func3->GetParError(first));
  AddPoint(width[4], x, TMath::Abs(func3->GetParameter(second)), 0, func3->GetParError(second));

  AddPoint(phiWidthSummary, phiWidthSummary->GetN(), TMath::Abs(func3->GetParameter(first)), 0, func3->GetParError(first));
  AddPoint(phiWidthSummary, phiWidthSummary->GetN(), TMath::Abs(func3->GetParameter(second)), 0, func3->GetParError(second));
    
  //deta
  AddPoint(width[2], x, TMath::Abs(func3->GetParameter(first+1)), 0, func3->GetParError(first+1));
  AddPoint(width[5], x, TMath::Abs(func3->GetParameter(second+1)), 0, func3->GetParError(second+1));
  
  AddPoint(etaWidthSummary, etaWidthSummary->GetN(), TMath::Abs(func3->GetParameter(first+1)), 0, func3->GetParError(first+1));
  AddPoint(etaWidthSummary, etaWidthSummary->GetN(), TMath::Abs(func3->GetParameter(second+1)), 0, func3->GetParError(second+1));
  
  // norm
/*  AddPoint(width[0], x, func3->GetParameter(first-1), 0, func3->GetParError(first-1));
  AddPoint(width[3], x, func3->GetParameter(second-1), 0, func3->GetParError(second-1));*/
  AddPoint(width[0], x, func3->GetParameter(1), 0, func3->GetParError(1));
  if (first < second)
    AddPoint(width[3], x, func3->GetParameter(4), 0, func3->GetParError(4));
  else
    AddPoint(width[3], x, 1.0 - func3->GetParameter(4), 0, func3->GetParError(4));

  TH2* funcHist3 = (TH2*) hist->Clone("funcHist3");
  funcHist3->Reset();
  funcHist3->Add(func3);
  funcHist3->SetMinimum(-(max - min) / 2);
  funcHist3->SetMaximum((max - min) / 2);
  funcHist3->Draw("SURF1");
  sumSummary->SetPoint(sumSummary->GetN(), sumSummary->GetN(), funcHist3->Integral());
  
  canvas->cd(canvasPos++);
  TH2* residuals3 = (TH2*) hist->Clone("residuals");
  residuals3->Add(funcHist3, -1);
  residuals3->SetMinimum(-(max - min) / 2);
  residuals3->SetMaximum((max - min) / 2);
  residuals3->Draw("SURF1");

  Float_t chi2 = 0;
  Int_t ndf = 0;
  for (Int_t i=hist->GetXaxis()->FindBin(-0.8); i<=hist->GetXaxis()->FindBin(0.8); i++)
    for (Int_t j=hist->GetYaxis()->FindBin(-etaFitLimit); j<=hist->GetYaxis()->FindBin(etaFitLimit); j++)
    {
      if (residuals3->GetBinError(i, j) > 0)
      {
	chi2 += TMath::Power(residuals3->GetBinContent(i, j) / residuals3->GetBinError(i, j), 2);
	ndf++;
      }
    }
  ndf -= func3->GetNumberFreeParameters();
  
  if (func3->GetNDF() > 0)
  {
    printf("#chi^{2}/ndf = %.1f/%d = %.1f  ", func3->GetChisquare(), func3->GetNDF(), func3->GetChisquare() / func3->GetNDF());
    DrawLatex(0.2, yPosChi2, 1, Form("#chi^{2}/ndf = %.1f/%d = %.1f", func3->GetChisquare(), func3->GetNDF(), func3->GetChisquare() / func3->GetNDF()));
    chi2_1->SetPoint(chi2_1->GetN(), x, func3->GetChisquare() / func3->GetNDF());
  }
  if (ndf)
  {
    Printf("#chi^{2}/ndf = %.1f/%d = %.1f", chi2, ndf, chi2 / ndf);
    DrawLatex(0.2, yPosChi2 - 0.05, 1, Form("#chi^{2}/ndf = %.1f/%d = %.1f", chi2, ndf, chi2 / ndf));
    chi2_2->SetPoint(chi2_2->GetN(), x, chi2 / ndf);
  }

  canvas->cd(canvasPos++);
  phiWidthSummary->SetMarkerStyle(20);
  phiWidthSummary->Draw("AP");
  gPad->SetGridy();
    
  etaWidthSummary->SetMarkerStyle(21);
  etaWidthSummary->SetLineColor(2);
  etaWidthSummary->SetMarkerColor(2);
  etaWidthSummary->Draw("PSAME");  
  
  phiWidthSummary->GetYaxis()->SetRangeUser(0.1, 0.9);
  
  canvas->cd(canvasPos++);
  sumSummary->Draw("*A");
  gPad->SetGridy();
}

// norm,dphi,deta,norm2,dphi2,deta2,chi2_1,chi2_2,moment2phi,moment2eta,moment3phi,moment3eta,moment4phi,moment4eta,kurtosisphi,kurtosiseta
const Int_t NGraphs = 16;
const Int_t NHists = 6;
TGraphErrors*** graphs = 0;
TList* labelList = 0;

void CreateGraphStructure()
{
  if (!graphs)
  {
    graphs = new TGraphErrors**[NGraphs];
    for (Int_t i=0; i<NGraphs; i++)
    {
      graphs[i] = new TGraphErrors*[NHists];
      for (Int_t j=0; j<NHists; j++)
	graphs[i][j] = new TGraphErrors;
    }
    
    labelList = new TList;
  }
}

void WriteGraphs()
{
  TFile::Open("graphs.root", "RECREATE");
  for (Int_t i=0; i<NGraphs; i++)
    for (Int_t j=0; j<NHists; j++)
      graphs[i][j]->Write(Form("graph_%d_%d", i, j));

  labelList->Write("labelList", TObject::kSingleKey);
    
  gFile->Close();
}

void ReadGraphs(const char* fileName = "graphs.root")
{
  CreateGraphStructure();
  TFile::Open(fileName);
  for (Int_t i=0; i<NGraphs; i++)
    for (Int_t j=0; j<NHists; j++)
      graphs[i][j] = (TGraphErrors*) gFile->Get(Form("graph_%d_%d", i, j));
    
  labelList = (TList*) gFile->Get("labelList");
}

void AnalyzeDeltaPhiEtaGap2D(const char* fileName, Int_t method)
{
  gROOT->SetBatch(kTRUE);
  if (!gROOT->IsBatch())
  {
    Printf("Not in batch mode. Exiting!");
    return;
  }
  
  CreateGraphStructure();

  TFile::Open(fileName);
  
  Int_t maxLeadingPt = 3;
  Int_t maxAssocPt = 6;

  for (Int_t histId = 0; histId < NHists; histId++)
  {
    for (Int_t i=0; i<maxLeadingPt; i++)
    {
      for (Int_t j=1; j<maxAssocPt; j++)
      {
// 	if (i != 1 || j != 2)
// 	  continue;
    
	TH2* hist1 = (TH2*) gFile->Get(Form("dphi_%d_%d_%d", i, j+1, histId));
	if (!hist1)
	  continue;
	
	TCanvas* canvas = new TCanvas(Form("DeltaPhi_%d_%d_%d_%d", histId, i, j, method), Form("DeltaPhi_%d_%d_%d_%d", histId, i, j, method), 1400, 1100);
	canvas->Divide(5, 3);
	
	for (Int_t k=1; k<=3; k++)
	{
	  canvas->cd(3 * j + k);
	  gPad->SetLeftMargin(0.15);
	  gPad->SetBottomMargin(0.2);
	  gPad->SetTopMargin(0.01);
	  gPad->SetRightMargin(0.01);
	}
	
	if (histId == 0)
	  labelList->Add(new TObjString(hist1->GetTitle()));
	
// 	hist1->Rebin2D(2, 1);
	
	Float_t xPos = j*8+i;
	
	TGraphErrors* width[14] = { graphs[0][histId], graphs[1][histId], graphs[2][histId], graphs[3][histId], graphs[4][histId], graphs[5][histId], graphs[8][histId], graphs[9][histId], graphs[10][histId], graphs[11][histId], graphs[12][histId], graphs[13][histId], graphs[14][histId], graphs[15][histId] };

	if (method == 0)
	  FitDeltaPhiEtaGap2D((TH2*) hist1, kFALSE, canvas, 1, width, xPos, 0.9, graphs[6][histId], graphs[7][histId]);
	else
	  FitDeltaPhi2DOneFunction((TH2*) hist1, canvas, 1, width, xPos, 0.9, graphs[6][histId], graphs[7][histId], kTRUE, histId, (i > 0) ? 1 : 0);
	
// 	canvas->SaveAs(Form("%s.png", canvas->GetName()));
// 	delete canvas;
// 	break;
// return;
      }
      
//       break;
    }
    
//     break;
  }
  
  WriteGraphs();
}

void AnalyzeDeltaPhiEtaGap2DExample(const char* fileName, Int_t i, Int_t j, Int_t histId)
{
  CreateGraphStructure();

  TFile::Open(fileName);
  
  TH2* hist1 = (TH2*) gFile->Get(Form("dphi_%d_%d_%d", i, j+1, histId));
  if (!hist1)
    return;
  
  TCanvas* canvas = new TCanvas(Form("DeltaPhi_%d_%d_%d_%d", histId, i, j, 1), Form("DeltaPhi_%d_%d_%d_%d", histId, i, j, 1), 1400, 1100);
  canvas->Divide(5, 3);
  
  Float_t xPos = j*8+i;
  
  TGraphErrors* width[14] = { graphs[0][histId], graphs[1][histId], graphs[2][histId], graphs[3][histId], graphs[4][histId], graphs[5][histId], graphs[8][histId], graphs[9][histId], graphs[10][histId], graphs[11][histId], graphs[12][histId], graphs[13][histId], graphs[14][histId], graphs[15][histId] };

  FitDeltaPhi2DOneFunction((TH2*) hist1, canvas, 1, width, xPos, 0.9, graphs[6][histId], graphs[7][histId], kFALSE, histId);
}

void SqrtAll(Int_t nHists, TGraphErrors** graph)
{
  for (Int_t histId = 0; histId < nHists; histId++)
  {
    for (Int_t i=0; i<graph[histId]->GetN(); i++)
    {
      if (graph[histId]->GetY()[i] > 0)
      {
	graph[histId]->GetEY()[i] *= 0.5 / TMath::Sqrt(graph[histId]->GetY()[i]);
	graph[histId]->GetY()[i] = TMath::Sqrt(graph[histId]->GetY()[i]);
      }
    }
  }
}

Int_t marker[6] = { 24, 25, 26, 30, 27, 28 };
Int_t marker2[6] = { 20, 21, 22, 29, 33, 34 };
const char* labels[6] = { "0-5%", "60-90%", "pp", "20-30%", "10-20%", "40-60%" };

void Draw(const char* canvasName, Int_t nHists, TGraphErrors** graph, TGraphErrors** graph2 = 0, Bool_t normalize = kFALSE, Float_t min = 0, Float_t max = 0)
{
  Int_t colors[6] = { 1, 2, 4, 6, 7, 8 };
  Bool_t found = kTRUE;
  TCanvas* c1 = (TCanvas*) gROOT->GetListOfCanvases()->FindObject(canvasName);
  if (!c1)
  {
    c1 = new TCanvas(canvasName, canvasName, 800, 600);
//     colors[0] = 1;
    found = kFALSE;
  }
  c1->cd();
  
  if (normalize && graph2)
  {
    // relative norms
    for (Int_t histId = 0; histId < nHists; histId++)
    {
      for (Int_t i=0; i<graph[histId]->GetN(); i++)
      {
	Float_t sum = graph[histId]->GetY()[i] + graph2[histId]->GetY()[i];
	graph[histId]->GetY()[i] /= sum;
	graph2[histId]->GetY()[i] /= sum;
	graph[histId]->GetEY()[i] /= sum;
	graph2[histId]->GetEY()[i] /= sum;
      }
    }
  }
  
  for (Int_t histId = 0; histId < nHists; histId++)
  {
    graph[histId]->SetMarkerStyle((found) ? marker2[histId] : marker[histId]);
    graph[histId]->SetMarkerColor(colors[histId]);
    graph[histId]->SetLineColor(colors[histId]);
    if (normalize)
      graph[histId]->GetYaxis()->SetRangeUser(0, 1);
    if (max > min)
      graph[histId]->GetYaxis()->SetRangeUser(min, max);
    graph[histId]->GetXaxis()->SetRangeUser(7, 38);
    graph[histId]->DrawClone((histId == 0 && !found) ? "AP" : "PSAME");
    DrawLatex(0.7, 0.8 - 0.05 * histId, colors[histId], labels[histId]);
  }
  
  if (graph2)
  {
    for (Int_t histId = 0; histId < nHists; histId++)
    {
      graph2[histId]->SetMarkerStyle(marker2[histId]);
      graph2[histId]->SetMarkerColor(colors[histId]);
      graph2[histId]->SetLineColor(colors[histId]);
      graph2[histId]->DrawClone("PSAME");
    }
  }
  
  gPad->SetGridx();
  gPad->SetGridy();
  c1->SaveAs(Form("%s.png", canvasName));
}

void DrawRMS(const char* canvasName, Int_t nHists, TGraphErrors** graph, TGraphErrors** graph2, TGraphErrors** weight, Float_t min = 0, Float_t max = 0)
{
Int_t colors[6] = { 2, 2, 4, 6, 7, 8 };
  Bool_t found = kTRUE;
  TCanvas* c1 = (TCanvas*) gROOT->GetListOfCanvases()->FindObject(canvasName);
  if (!c1)
  {
    c1 = new TCanvas(canvasName, canvasName, 800, 600);
    colors[0] = 1;
    found = kFALSE;
  }
  c1->cd();
  
  for (Int_t histId = 0; histId < nHists; histId++)
  {
    graph[histId]->SetMarkerStyle(marker[histId]);
    graph[histId]->SetMarkerColor(colors[histId]);
    graph[histId]->SetLineColor(colors[histId]);
    if (max > min)
      graph[histId]->GetYaxis()->SetRangeUser(min, max);
    graph[histId]->DrawClone((histId == 0 && !found) ? "AP" : "PSAME");
    DrawLatex(0.7, 0.8 - 0.05 * histId, colors[histId], labels[histId]);

    graph2[histId]->SetMarkerStyle(marker2[histId]);
    graph2[histId]->SetMarkerColor(colors[histId]);
    graph2[histId]->SetLineColor(colors[histId]);
    graph2[histId]->DrawClone("PSAME");
  }

  // rms
  for (Int_t histId = 0; histId < nHists; histId++)
  {
    for (Int_t i=0; i<graph[histId]->GetN(); i++)
    {
      Float_t rms = graph[histId]->GetY()[i] * weight[histId]->GetY()[i] + graph2[histId]->GetY()[i] * (1.0 - weight[histId]->GetY()[i]);
      graph[histId]->GetY()[i] = rms;
      // error**2 = weight**2 * error_sigma**2 + (weight-1)**2 * error_sigma2**2 + (sigma1 + sigma2)**2 * error_weight**2
      // TODO this neglects some correlations
      graph[histId]->GetEY()[i] = TMath::Sqrt(TMath::Power(weight[histId]->GetY()[i] * graph[histId]->GetEY()[i], 2) + 
	TMath::Power((weight[histId]->GetY()[i] - 1) * graph2[histId]->GetEY()[i], 2) + TMath::Power((graph[histId]->GetY()[i] + graph2[histId]->GetY()[i]) * weight[histId]->GetEY()[i], 2));
    }
    graph[histId]->DrawClone("*SAME");
  }

  gPad->SetGridx();
  gPad->SetGridy();
  c1->SaveAs(Form("%s.png", canvasName));
}

void DrawCentrality(const char* canvasName, Int_t nHists, TGraphErrors** graph, Float_t min = 0, Float_t max = 0)
{
//const char* labels[6] = { "0-5%", "60-90%", "pp", "20-30%", "10-20%", "40-60%" };

  Int_t colors[6] = { 1, 2, 4, 6, 7, 8 };
  Bool_t found = kTRUE;
  TCanvas* c1 = (TCanvas*) gROOT->GetListOfCanvases()->FindObject(canvasName);
  if (!c1)
  {
    c1 = new TCanvas(canvasName, canvasName, 800, 600);
    found = kFALSE;
  }
  c1->cd();
  
  TLegend* legend = new TLegend(0.5, 0.5, 0.8, 0.8);
  
  Int_t beginWidth = 1;
  for (Int_t i=beginWidth; i<graph[0]->GetN(); i++)
  {
    if (i == 9 || i == 11)
      continue;
    
    TGraphErrors* graphcentrality = new TGraphErrors;
    graphcentrality->GetXaxis()->SetTitle("Centrality");
    graphcentrality->SetMarkerStyle(20+i);
    graphcentrality->SetMarkerColor(i+1);
    graphcentrality->SetLineColor(i+1);
    
    Float_t centralityAxisMapping[] = { 2.5, 75, 100, 25, 15, 50 };
    Float_t centralityAxisMappingE[] = { 2.5, 15, 0, 5, 5, 10 };
    
    for (Int_t histId = 0; histId < nHists; histId++)
    {
      AddPoint(graphcentrality, centralityAxisMapping[histId], graph[histId]->GetY()[i], centralityAxisMappingE[histId], graph[histId]->GetEY()[i]);
//       Printf("%f %f", graph[histId]->GetY()[i], graph[histId]->GetEY()[i]);
//       graphcentrality->Print();
    }
    
//     return;
    
    if (max > min)
      graphcentrality->GetYaxis()->SetRangeUser(min, max);
    graphcentrality->DrawClone((i == beginWidth && !found) ? "AP" : "PSAME");
    legend->AddEntry(graphcentrality->Clone(), (labelList) ? labelList->At(i)->GetName() : "", "PL");
  }
  
  legend->Draw();
  
  gPad->SetGridx();
  gPad->SetGridy();
//   c1->SaveAs(Form("%s.png", canvasName));
}

void DrawResults(const char* fileName = "graphs.root")
{
  ReadGraphs(fileName);
  
//   return;

  Int_t nHists = 6;
  
  if (1)
  {
    Draw("width_phi", nHists, graphs[1], graphs[4], kFALSE, 0, 0.8);
    Draw("width_eta", nHists, graphs[2], graphs[5], kFALSE, 0, 0.8);
    Draw("norm", nHists, graphs[0], graphs[3], kFALSE, 0, 1);
    Draw("chi2", nHists, graphs[6], graphs[7], kFALSE, 0.5, 2.5);
    DrawRMS("width_phi_rms", nHists, graphs[1], graphs[4], graphs[3], 0, 0.8);
    DrawRMS("width_eta_rms", nHists, graphs[2], graphs[5], graphs[3], 0, 0.8);
  }
  
  Draw("moment2", nHists, graphs[8], graphs[9], kFALSE, 0, 0.3);
  Draw("moment3", nHists, graphs[10], graphs[11], kFALSE, -0.05, 0.05);
  Draw("moment4", nHists, graphs[12], graphs[13], kFALSE, -0.1, 0.2);
  Draw("kurtosis", nHists, graphs[14], graphs[15], kFALSE, -2, 2);
  Draw("kurtosisphi", nHists, graphs[14], 0, kFALSE, -2, 2);
  Draw("kurtosiseta", nHists, graphs[15], 0, kFALSE, -2, 2);
}

void DrawResultsCentrality(const char* fileName = "graphs.root")
{
  ReadGraphs(fileName);

  Int_t nHists = 6;

  DrawCentrality("width_phi1_centrality", nHists, graphs[1], 0, 0.8);
  DrawCentrality("width_phi2_centrality", nHists, graphs[4], 0, 0.8);
  DrawCentrality("width_eta1_centrality", nHists, graphs[2], 0, 0.8);
  DrawCentrality("width_eta2_centrality", nHists, graphs[5], 0, 0.8);
  
  DrawRMS("width_phi_rms", nHists, graphs[1], graphs[4], graphs[3], 0, 0.8);
  DrawRMS("width_eta_rms", nHists, graphs[2], graphs[5], graphs[3], 0, 0.8);

  DrawCentrality("phi_rms_centrality", nHists, graphs[1], 0, 0.8);
  DrawCentrality("eta_rms_centrality", nHists, graphs[2], 0, 0.8);

  DrawCentrality("chi2_1", nHists, graphs[6], 0.5, 2.5);
  DrawCentrality("chi2_2", nHists, graphs[7], 0.5, 2.5);

  // sqrt(moment2) = sigma
  SqrtAll(nHists, graphs[8]);
  SqrtAll(nHists, graphs[9]);
  DrawCentrality("sigma_phi", nHists, graphs[8], 0, 0.8);
  DrawCentrality("sigma_eta", nHists, graphs[9], 0, 0.8);
  
  DrawCentrality("kurtosisphi_centrality", nHists, graphs[14], -2, 2);
  DrawCentrality("kurtosiseta_centrality", nHists, graphs[15], -2, 2);
}

void CompareRMS(const char* graphFileName1, const char* graphFileName2)
{
  ReadGraphs(graphFileName1);

  Int_t nHists = 3;

  DrawRMS("width_phi_rms1", nHists, graphs[1], graphs[4], graphs[3], 0, 0.8);
  DrawRMS("width_eta_rms1", nHists, graphs[2], graphs[5], graphs[3], 0, 0.8);
  
  Draw("width_phi_rms_compare", nHists, graphs[1], 0, kFALSE, 0, 0.8);
  Draw("width_eta_rms_compare", nHists, graphs[2], 0, kFALSE, 0, 0.8);
  
  ReadGraphs(graphFileName2);
  DrawRMS("width_phi_rms2", nHists, graphs[1], graphs[4], graphs[3], 0, 0.8);
  DrawRMS("width_eta_rms2", nHists, graphs[2], graphs[5], graphs[3], 0, 0.8);

  Draw("width_phi_rms_compare", nHists, graphs[1], 0, kFALSE, 0, 0.8);
  Draw("width_eta_rms_compare", nHists, graphs[2], 0, kFALSE, 0, 0.8);
}

void DrawExamples(const char* histFileName, const char* graphFileName, Bool_t drawFunc = kTRUE)
{
  Float_t etaLimit = 1.0;
  Float_t outerLimit = 1.79;
  Float_t projectLimit = 0.8;

  ReadGraphs(graphFileName);

  TFile::Open(histFileName);
  
  Int_t j=1;

  Int_t exColors[] = { 1, 2, 4, 6 };
  
  for (Int_t i=0; i<6; i++)
  {
    TCanvas* c = new TCanvas(Form("c_%d", i), Form("c_%d", i), 1000, 600);
    c->Divide(2, 1);
    Int_t nHists = 3;
    for (Int_t histId = 0; histId < nHists; histId++)
    {
      TH2* hist = (TH2*) gFile->Get(Form("dphi_%d_%d_%d", i, j+1, histId));
      if (!hist)
	continue;

      Float_t xPos = j*8+i;
      Int_t n = 0;
      while (graphs[0][histId]->GetX()[n] != xPos && n < graphs[0][histId]->GetN())
	n++;
      if (n == graphs[0][histId]->GetN())
      {
	Printf("ERROR: Pos in graph not found: %d %d", i, histId);
	continue;
      }
      
      SubtractEtaGap(hist, etaLimit, outerLimit, kFALSE);
    
      c->cd(1);
      TH1* proj = hist->ProjectionX(Form("%s_proj1", hist->GetName()), hist->GetYaxis()->FindBin(-projectLimit+0.01), hist->GetYaxis()->FindBin(projectLimit-0.01));
      proj->SetLineColor(exColors[histId]);
      proj->GetXaxis()->SetRangeUser(-1.49, 1.49);
      proj->GetYaxis()->SetRangeUser(proj->GetMinimum() * 1.2, proj->GetMaximum() * 1.5);
      proj->Draw((histId == 0) ? "" : "SAME");
    
      if (drawFunc)
      {
	// integral over y
	TF1* func = new TF1("func", "[0]+[1]*([4]*[3]*TMath::Sqrt(TMath::TwoPi())*exp(-0.5*(x/[2])**2)+(1-[4])*[6]*TMath::Sqrt(TMath::TwoPi())*exp(-0.5*(x/[5])**2))", -0.5 * TMath::Pi(), 0.5 * TMath::Pi());
	func->SetParameter(0, 0);
	for (Int_t k=0; k<6; k++)
	  func->SetParameter(k+1, graphs[k][histId]->GetY()[n]);
	// scale by bin width (to compare with projection)
	func->SetParameter(1, func->GetParameter(1) / hist->GetYaxis()->GetBinWidth(1));
	func->SetLineColor(exColors[histId]);
	func->SetLineWidth(1);
	func->DrawCopy("SAME");
	
	// draw contributions
	Float_t scale = func->GetParameter(1);
	Float_t weighting = func->GetParameter(4);
	func->SetParameter(1, scale * weighting);
	func->SetParameter(4, 1);
	func->SetLineStyle(2);
	func->DrawCopy("SAME");

	func->SetParameter(1, scale * (1.0-weighting));
	func->SetParameter(4, 0);
	func->SetLineStyle(3);
	func->DrawCopy("SAME");
      }

      DrawLatex(0.15, 0.8 - 0.05 * histId, exColors[histId], labels[histId]);

      c->cd(2);
      proj = hist->ProjectionY(Form("%s_proj2b", hist->GetName()), hist->GetXaxis()->FindBin(-projectLimit+0.01), hist->GetXaxis()->FindBin(projectLimit-0.01));
      proj->SetLineColor(exColors[histId]);
      proj->GetXaxis()->SetRangeUser(-1.49, 1.49);
      proj->GetYaxis()->SetRangeUser(proj->GetMinimum() * 1.2, proj->GetMaximum() * 2);
      proj->Draw((histId == 0) ? "" : "SAME");
      
      if (0)
      {
	proj = hist->ProjectionY(Form("%s_proj2", hist->GetName()), hist->GetXaxis()->FindBin(-0.5 * TMath::Pi()), hist->GetXaxis()->FindBin(0.5 * TMath::Pi()));
	proj->SetLineColor(exColors[histId]);
	proj->GetXaxis()->SetRangeUser(-1.2, 1.2);
	proj->GetYaxis()->SetRangeUser(proj->GetMinimum() * 1.2, proj->GetMaximum() * 2);
	proj->SetLineStyle(2);
	proj->Draw("SAME");
      }

      if (drawFunc)
      {
	// integral over x
	TF1* func = new TF1("func", "[0]+[1]*([4]*[2]*TMath::Sqrt(TMath::TwoPi())*exp(-0.5*(x/[3])**2)+(1-[4])*[5]*TMath::Sqrt(TMath::TwoPi())*exp(-0.5*(x/[6])**2))", -outerLimit, outerLimit);
	func->SetParameter(0, 0);
	for (Int_t k=0; k<6; k++)
	  func->SetParameter(k+1, graphs[k][histId]->GetY()[n]);
	// scale by bin width (to compare with projection)
	func->SetParameter(1, func->GetParameter(1) / hist->GetXaxis()->GetBinWidth(1));
	func->SetLineColor(exColors[histId]);
	func->SetLineWidth(1);
	func->DrawCopy("SAME");

	// draw contributions
	Float_t scale = func->GetParameter(1);
	Float_t weighting = func->GetParameter(4);
	func->SetParameter(1, scale * weighting);
	func->SetParameter(4, 1);
	func->SetLineStyle(2);
	func->DrawCopy("SAME");

	func->SetParameter(1, scale * (1.0-weighting));
	func->SetParameter(4, 0);
	func->SetLineStyle(3);
	func->DrawCopy("SAME");   
      }

      DrawLatex(0.15, 0.8 - 0.05 * histId, exColors[histId], labels[histId]);
    }
  }
}

void CompareExamples(const char* histFileName1, const char* histFileName2, Int_t i, Int_t j, Int_t histId)
{
  Float_t etaLimit = 1.0;
  Float_t outerLimit = 1.49;
  Float_t projectLimit = 0.8;

  Int_t exColors[] = { 1, 2, 4, 6 };
  
  TCanvas* c = new TCanvas(Form("c_%d", i), Form("c_%d", i), 1000, 600);
  c->Divide(2, 1);

  for (Int_t k=0; k<2; k++)
  {
    if (k == 0)
      TFile::Open(histFileName1);
    else
      TFile::Open(histFileName2);
  
    TH2* hist = (TH2*) gFile->Get(Form("dphi_%d_%d_%d", i, j+1, histId));
    if (!hist)
      continue;

    SubtractEtaGap(hist, etaLimit, outerLimit, kFALSE);
  
    c->cd(1);
    TH1* proj = hist->ProjectionX(Form("%s_proj1%d", hist->GetName(), k), hist->GetYaxis()->FindBin(-projectLimit+0.01), hist->GetYaxis()->FindBin(projectLimit-0.01));
    proj->SetLineColor(exColors[k]);
    proj->GetXaxis()->SetRangeUser(-1.49, 1.49);
    proj->GetYaxis()->SetRangeUser(proj->GetMinimum() * 1.2, proj->GetMaximum() * 1.5);
    proj->Draw((k == 0) ? "" : "SAME");
  
    c->cd(2);
    proj = hist->ProjectionY(Form("%s_proj2_%d", hist->GetName(), k), hist->GetXaxis()->FindBin(-projectLimit+0.01), hist->GetXaxis()->FindBin(projectLimit-0.01));
    proj->SetLineColor(exColors[k]);
    proj->GetXaxis()->SetRangeUser(-1.49, 1.49);
    proj->GetYaxis()->SetRangeUser(proj->GetMinimum() * 1.2, proj->GetMaximum() * 2);
    proj->Draw((k == 0) ? "" : "SAME");
  }
}

void TestMomentCode()
{
  Float_t momentLimit = 1;
  TH1* proj = new TH1F("hist", "", 100, -momentLimit, momentLimit);
  TF1* func = new TF1("func", "gaus(0)", -momentLimit, momentLimit);
  func->SetParameters(1, 0, 0.2);
  proj->FillRandom("func", 100);
  proj->Sumw2();
  proj->Draw();
  
  for (Int_t n=2; n <= 4; n++)
  {
    Float_t moment = 0;
    Float_t sum = 0;
    for (Int_t bin = proj->GetXaxis()->FindBin(-momentLimit); bin <= proj->GetXaxis()->FindBin(momentLimit); bin++)
    {
      moment += proj->GetBinContent(bin) * TMath::Power(proj->GetMean() - proj->GetXaxis()->GetBinCenter(bin), n);
      sum += proj->GetBinContent(bin);
    }
    moment /= sum;
    
    Float_t error = 0;
    for (Int_t bin = proj->GetXaxis()->FindBin(-momentLimit); bin <= proj->GetXaxis()->FindBin(momentLimit); bin++)
    {
      error += proj->GetBinError(bin) * proj->GetBinError(bin) * 
	TMath::Power(TMath::Power(proj->GetMean() - proj->GetXaxis()->GetBinCenter(bin), n) / sum 
	  - moment / sum, 2);
    }
    
    Printf("%d %f %f", n, moment, TMath::Sqrt(error));
  }
}

void DivideGraphs(TGraphErrors* graph1, TGraphErrors* graph2)
{
  graph1->Print();
  graph2->Print();

  for (Int_t bin1 = 0; bin1 < graph1->GetN(); bin1++)
  {
    Float_t x = graph1->GetX()[bin1];

    Int_t bin2 = 0;
    for (bin2 = 0; bin2<graph2->GetN(); bin2++)
      if (graph2->GetX()[bin2] >= x)
        break;

    if (bin2 == graph2->GetN())
            bin2--;

    if (bin2 > 0)
      if (TMath::Abs(graph2->GetX()[bin2-1] - x) < TMath::Abs(graph2->GetX()[bin2] - x))
        bin2--;

    if (graph1->GetY()[bin1] == 0 || graph2->GetY()[bin2] == 0 || bin2 == graph2->GetN())
    {
      Printf("%d %d removed", bin1, bin2);
      graph1->RemovePoint(bin1--);
      continue;
    }

    Float_t graph2Extrapolated = graph2->GetY()[bin2];
    if (TMath::Abs(x - graph2->GetX()[bin2]) > 0.0001)
    {
      Printf("%f %f %d %d not found", x, graph2->GetX()[bin2], bin1, bin2);
      graph1->RemovePoint(bin1--);
      continue;
    }

    Float_t value = graph1->GetY()[bin1] / graph2Extrapolated;
    Float_t error = value * TMath::Sqrt(TMath::Power(graph1->GetEY()[bin1] / graph1->GetY()[bin1], 2) + TMath::Power(graph2->GetEY()[bin2] / graph2->GetY()[bin2], 2));

    graph1->GetY()[bin1] = value;
    graph1->GetEY()[bin1] = error;

    Printf("%d %d %f %f %f %f", bin1, bin2, x, graph2Extrapolated, value, error);
  }
}

void CompareGraph(const char* fileName1, const char* fileName2, Int_t graph, Int_t centrality)
{
  ReadGraphs(fileName1);
  TGraphErrors* graph1 = (TGraphErrors*) graphs[graph][centrality]->Clone("graph1");
  
  ReadGraphs(fileName2);
  TGraphErrors* graph2 = (TGraphErrors*) graphs[graph][centrality]->Clone("graph2");
  
  graph1->Sort();
  graph2->Sort();
  
  TCanvas* c = new TCanvas(Form("%s_%s", fileName1, fileName2), Form("%s_%s", fileName1, fileName2), 800, 800);
  c->Divide(1, 2);
  
  c->cd(1);
  gPad->SetGridx();
  gPad->SetGridy();
  graph1->SetMarkerStyle(24);
  graph1->GetXaxis()->SetRangeUser(7, 38);
  graph1->DrawClone("AP");
  graph2->SetLineColor(2);
  graph2->SetMarkerColor(2);
  graph2->SetMarkerStyle(25);
  graph2->DrawClone("PSAME");
  
  c->cd(2);
  gPad->SetGridx();
  gPad->SetGridy();
  DivideGraphs(graph1, graph2);
  graph1->DrawClone("AP");
}

void Compare2D(const char* fileName1, const char* fileName2, const char* histName)
{
  TFile::Open(fileName1);
  TH1* hist1 = (TH1*) gFile->Get(histName);

  TFile::Open(fileName2);
  TH1* hist2 = (TH1*) gFile->Get(histName);
  
  hist1->Divide(hist2);
  hist1->GetZaxis()->SetRangeUser(0.5, 1.5);
  hist1->Draw("COLZ");
}
