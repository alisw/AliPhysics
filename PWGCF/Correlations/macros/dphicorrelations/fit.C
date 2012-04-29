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
#include "TLegendEntry.h"
#include "TMultiGraph.h"
#include "TPaveText.h"

void AddPoint(TGraphErrors* graph, Float_t x, Float_t y, Float_t xe, Float_t ye)
{
	graph->SetPoint(graph->GetN(), x, y);
	graph->SetPointError(graph->GetN() - 1, xe, ye);
}

void DrawLatex(Float_t x, Float_t y, Int_t color, const char* text, Float_t textSize = 0.06)
{
  TLatex* latex = new TLatex(x, y, text);
  latex->SetNDC();
  latex->SetTextSize(textSize);
  latex->SetTextColor(color);
  latex->Draw();
}

// 0    1    2    3     4     5     6      7      8          9          10         11         12         13         14          15
// norm,dphi,deta,norm2,dphi2,deta2,chi2_1,chi2_2,moment2phi,moment2eta,moment3phi,moment3eta,moment4phi,moment4eta,kurtosisphi,kurtosiseta,
//            16   17   18   19    20    21    22     23     24         25         26         27         28         29         30          31
// second fit:norm,dphi,deta,norm2,dphi2,deta2,chi2_1,chi2_2,moment2phi,moment2eta,moment3phi,moment3eta,moment4phi,moment4eta,kurtosisphi,kurtosiseta,
const Int_t NGraphs = 32;
const Int_t NHists = 6*4; // pt index
TGraphErrors*** graphs = 0;

void CreateGraphStructure()
{
  graphs = new TGraphErrors**[NGraphs];
  for (Int_t i=0; i<NGraphs; i++)
  {
    graphs[i] = new TGraphErrors*[NHists];
    for (Int_t j=0; j<NHists; j++)
      graphs[i][j] = new TGraphErrors;
  }
}

void WriteGraphs()
{
  TFile::Open("graphs.root", "RECREATE");
  for (Int_t i=0; i<NGraphs; i++)
    for (Int_t j=0; j<NHists; j++)
      graphs[i][j]->Write(Form("graph_%d_%d", i, j));

  gFile->Close();
}

void ReadGraphs(const char* fileName = "graphs.root")
{
  CreateGraphStructure();
  TFile* file = TFile::Open(fileName);
  for (Int_t i=0; i<NGraphs; i++)
    for (Int_t j=0; j<NHists; j++)
      graphs[i][j] = (TGraphErrors*) file->Get(Form("graph_%d_%d", i, j));
}

void DivideGraphs(TGraphErrors* graph1, TGraphErrors* graph2)
{
//   graph1->Print();
//   graph2->Print();

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

//     Printf("%d %d %f %f %f %f", bin1, bin2, x, graph2Extrapolated, value, error);
  }
}

void SubtractGraphs(TGraphErrors* graph1, TGraphErrors* graph2)
{
//   graph1->Print();
//   graph2->Print();

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

    Float_t value = graph1->GetY()[bin1] - graph2Extrapolated;
    Float_t error = TMath::Sqrt(TMath::Power(graph1->GetEY()[bin1], 2) + TMath::Power(graph2->GetEY()[bin2], 2));

    graph1->GetY()[bin1] = value;
    graph1->GetEY()[bin1] = error;

//     Printf("%d %d %f %f %f %f", bin1, bin2, x, graph2Extrapolated, value, error);
  }
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

const Double_t k1OverSqrtTwoPi = 1.0 / TMath::Sqrt(TMath::TwoPi());

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
  
  return par[6+phiBin] + par[0]*( 
      par[3]/TMath::TwoPi()/par[1]/par[2] * 
	TMath::Exp(-0.5*((x[0]/par[1])*(x[0]/par[1])+(x[1]/par[2])*(x[1]/par[2]))) 
      + (1-par[3])/TMath::TwoPi()/par[4]/par[5] 
	* TMath::Exp(-0.5*((x[0]/par[4])*(x[0]/par[4])+(x[1]/par[5])*(x[1]/par[5]))) 
    );
}

void SubtractEtaGap(TH2* hist, Float_t etaLimit, Float_t outerLimit, Bool_t scale, Bool_t drawEtaGapDist = kFALSE)
{
  TString histName(hist->GetName());
  Int_t etaBins = 0;

  TH1D* etaGap = hist->ProjectionX(histName + "_1", TMath::Max(1, hist->GetYaxis()->FindBin(-outerLimit + 0.01)), hist->GetYaxis()->FindBin(-etaLimit - 0.01));
  Printf("%f", etaGap->GetEntries());
  if (etaGap->GetEntries() > 0)
    etaBins += hist->GetYaxis()->FindBin(-etaLimit - 0.01) - TMath::Max(1, hist->GetYaxis()->FindBin(-outerLimit + 0.01)) + 1;

  TH1D* tracksTmp = hist->ProjectionX(histName + "_2", hist->GetYaxis()->FindBin(etaLimit + 0.01), TMath::Min(hist->GetYaxis()->GetNbins(), hist->GetYaxis()->FindBin(outerLimit - 0.01)));
  Printf("%f", tracksTmp->GetEntries());
  if (tracksTmp->GetEntries() > 0)
    etaBins += TMath::Min(hist->GetYaxis()->GetNbins(), hist->GetYaxis()->FindBin(outerLimit - 0.01)) - hist->GetYaxis()->FindBin(etaLimit + 0.01) + 1;
  
  etaGap->Add(tracksTmp);

  // get per bin result
  etaGap->Scale(1.0 / etaBins);
  
  if (drawEtaGapDist)
  {
    TH1D* centralRegion = hist->ProjectionX(histName + "_3", hist->GetYaxis()->FindBin(-etaLimit + 0.01), hist->GetYaxis()->FindBin(etaLimit - 0.01));
    
    centralRegion->Scale(1.0 / (hist->GetYaxis()->FindBin(etaLimit - 0.01) - hist->GetYaxis()->FindBin(-etaLimit + 0.01) + 1));

    TCanvas* c = new TCanvas("SubtractEtaGap", "SubtractEtaGap", 800, 800);
    centralRegion->SetStats(0);
    centralRegion->Draw();
    etaGap->DrawCopy("SAME")->SetLineColor(2);
    c->SaveAs("etagap_proj.eps");
  }
  
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

/*void FitDeltaPhiEtaGap2D(TH2* hist, Bool_t scale,  TCanvas* canvas, Int_t canvasPos, TGraphErrors** width, Float_t x, Float_t yPosChi2, TGraphErrors* chi2_1, TGraphErrors* chi2_2)
{
  Float_t etaLimit = 1.0;
  Float_t outerLimit = 1.6;
  
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

//   TF2* func = new TF2("func", DeltaPhiWidth2DFitFunction, -TMath::Pi() / 2, TMath::Pi() * 1.5, -1, 1, 4);
//   func->SetParameters(1, 0.3, 0.3, 0);
//   func->SetParLimits(0, 0, 10);
//   func->SetParLimits(1, 0.1, 10);
//   func->SetParLimits(2, 0.1, 10);
  
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
*/

void CalculateMomentsKurtosis(Float_t momentCalcLimit, TH1* proj, Int_t graphIDStart, Int_t histId, Float_t x, Float_t xE)
{
  //   momentCalcLimit = 0.6;
  for (Int_t n=2; n <= 4; n++)
  {
    Float_t moment = 0;
    Float_t sum = 0;
    for (Int_t bin = proj->GetXaxis()->FindBin(-momentCalcLimit); bin <= proj->GetXaxis()->FindBin(momentCalcLimit); bin++)
    {
      moment += proj->GetBinContent(bin) * TMath::Power(proj->GetMean() - proj->GetXaxis()->GetBinCenter(bin), n);
      sum += proj->GetBinContent(bin);
    }
    Printf("%f %f", moment, sum);
    moment /= sum;
    
    Float_t error = 0;
    for (Int_t bin = proj->GetXaxis()->FindBin(-momentCalcLimit); bin <= proj->GetXaxis()->FindBin(momentCalcLimit); bin++)
    {
      error += proj->GetBinError(bin) * proj->GetBinError(bin) * 
	TMath::Power(TMath::Power(proj->GetMean() - proj->GetXaxis()->GetBinCenter(bin), n) / sum 
	  - moment / sum, 2);
    }
    
    AddPoint(graphs[graphIDStart+(n-2)*2][histId], x, moment, xE, TMath::Sqrt(error));
    Printf("%d %d %f +- %f <-> %f +- %f", n, graphIDStart+(n-2)*2, moment, TMath::Sqrt(error), proj->GetRMS() * proj->GetRMS(), 2 * proj->GetRMSError() / proj->GetRMS() * proj->GetRMS() * proj->GetRMS());
  }
  
  proj->GetXaxis()->SetRangeUser(-momentCalcLimit, momentCalcLimit);
  AddPoint(graphs[graphIDStart+6][histId], x, proj->GetKurtosis(1), xE, proj->GetKurtosis(11));
}

void FitDeltaPhi2DOneFunction(TH2* hist, TCanvas* canvas, Int_t canvasPos, Int_t graphID, Float_t x, Float_t xE, Float_t yPosChi2, Bool_t quick, Int_t histId, Int_t limits, Bool_t twoTrack)
{
  Float_t etaLimit = 1.0;
  Float_t outerLimit = 1.59;
  Float_t sigmaFitLimit = 0.1 - limits * 0.05;
  Float_t etaFitUpperLimit = 0.8;
  Float_t initSigma = 0.6;
  if (histId == 2) // pp
  {
    etaFitUpperLimit = 0.6;
    initSigma = 0.4;
  }

  hist->GetYaxis()->SetRangeUser(-outerLimit+0.01, outerLimit-0.01);
  hist->GetXaxis()->SetRangeUser(-TMath::Pi() / 2 + 0.01, TMath::Pi() * 0.5 - 0.01);
  
  Float_t mean = hist->Integral(hist->GetYaxis()->FindBin(-TMath::Pi() / 2), hist->GetYaxis()->FindBin(TMath::Pi() / 2), hist->GetYaxis()->FindBin(1.0), hist->GetYaxis()->FindBin(outerLimit)) / (hist->GetYaxis()->FindBin(TMath::Pi() / 2) - hist->GetYaxis()->FindBin(-TMath::Pi() / 2)) / (hist->GetYaxis()->FindBin(outerLimit) - hist->GetYaxis()->FindBin(1.0) + 1);
//   Printf("%f", mean);

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

  Int_t bins = hist->GetNbinsX() / 2 / 2;
    
  TF2* func = new TF2("func", DeltaPhiWidth2DFitFunction, -TMath::Pi() / 2, TMath::Pi() * 0.5, -outerLimit, outerLimit, bins+6);
  func->SetParameters(1, 0.3, 0.3, 0.25, initSigma, initSigma);
  for (Int_t i=6; i<bins+6; i++)
    func->SetParameter(i, mean);

  func->SetParLimits(0, 0, 10);
  func->SetParLimits(1, sigmaFitLimit, 0.6);
  func->SetParLimits(2, sigmaFitLimit, etaFitUpperLimit);
  func->SetParLimits(3, 0.1, 0.9);
  func->SetParLimits(4, sigmaFitLimit, 0.6);
  func->SetParLimits(5, sigmaFitLimit, etaFitUpperLimit);

  Int_t fitResult = hist->Fit(func, "0R", "");
  Printf("Fit result: %d", fitResult);

  // if both parameters are within 1%, refit with 1 Gaussian only
  if (TMath::Abs(1.0 - func->GetParameter(1) / func->GetParameter(4)) < 0.01 && TMath::Abs(1.0 - func->GetParameter(2) / func->GetParameter(5)) < 0.01)
  {
    Printf("Parameters within 1%%. Refitting with 1 Gaussian...");
    
    func->SetParLimits(3, 1, 1);
    func->FixParameter(3, 1);
    func->FixParameter(4, sigmaFitLimit);
    func->FixParameter(5, sigmaFitLimit);
    
    fitResult = hist->Fit(func, "0R", "");
    Printf("Fit result: %d", fitResult);
  }
  
  Int_t first = 1;
  Int_t second = 4;
  if (func->GetParameter(1) < func->GetParameter(4))
  {
    first = 4;
    second = 1;
  }
  //dphi
  AddPoint(graphs[1+16][graphID], x, TMath::Abs(func->GetParameter(first)), xE, func->GetParError(first));
  AddPoint(graphs[4+16][graphID], x, TMath::Abs(func->GetParameter(second)), xE, func->GetParError(second));

  //deta
  AddPoint(graphs[2+16][graphID], x, TMath::Abs(func->GetParameter(first+1)), xE, func->GetParError(first+1));
  AddPoint(graphs[5+16][graphID], x, TMath::Abs(func->GetParameter(second+1)), xE, func->GetParError(second+1));
  
  // norm
  AddPoint(graphs[0+16][graphID], x, func->GetParameter(0), xE, func->GetParError(0));
  if (first < second)
    AddPoint(graphs[3+16][graphID], x, func->GetParameter(3), xE, func->GetParError(3));
  else
    AddPoint(graphs[3+16][graphID], x, 1.0 - func->GetParameter(3), xE, func->GetParError(3));    
  
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

  Float_t chi2 = 0;
  Int_t ndf = 0;
  for (Int_t i=hist->GetXaxis()->FindBin(-0.8); i<=hist->GetXaxis()->FindBin(0.8); i++)
    for (Int_t j=hist->GetYaxis()->FindBin(-etaLimit); j<=hist->GetYaxis()->FindBin(etaLimit); j++)
    {
      if (residuals->GetBinError(i, j) > 0)
      {
	chi2 += TMath::Power(residuals->GetBinContent(i, j) / residuals->GetBinError(i, j), 2);
	ndf++;
      }
    }
  ndf -= func->GetNumberFreeParameters();
  
  if (func->GetNDF() > 0)
  {
    printf("#chi^{2}/ndf = %.1f/%d = %.1f  ", func->GetChisquare(), func->GetNDF(), func->GetChisquare() / func->GetNDF());
    DrawLatex(0.2, yPosChi2, 1, Form("#chi^{2}/ndf = %.1f/%d = %.1f", func->GetChisquare(), func->GetNDF(), func->GetChisquare() / func->GetNDF()));
  }
  if (ndf)
  {
    Printf("#chi^{2}/ndf = %.1f/%d = %.1f", chi2, ndf, chi2 / ndf);
    DrawLatex(0.2, yPosChi2 - 0.05, 1, Form("#chi^{2}/ndf = %.1f/%d = %.1f", chi2, ndf, chi2 / ndf));
  }
  
  // draw gaussian only
  TF2* funcClone = new TF2("funcClone", DeltaPhiWidth2DFitFunction, -TMath::Pi() / 2, TMath::Pi() * 0.5, -outerLimit, outerLimit, bins+6);
  for (Int_t i=0; i<6; i++)
    funcClone->SetParameter(i, func->GetParameter(i));
  for (Int_t i=6; i<bins+6; i++)
    funcClone->SetParameter(i, 0);
//   funcClone->Print();
  canvas->cd(canvasPos++);
  funcHist = (TH2*) hist->Clone("funcHistb");
  funcHist->Reset();
  funcHist->Add(funcClone);
  funcHist->SetMinimum(-(max - min) / 2);
  funcHist->SetMaximum((max - min) / 2);
  funcHist->Draw("SURF1");
  
  // eta gap subtraction
  canvas->cd(canvasPos++);
  func->SetParameter(0, 0);
  TH2* subtractFlow = (TH2*) hist->Clone("subtractFlow");
  subtractFlow->Add(func, -1);
  subtractFlow->SetMinimum(-(max - min) / 2);
  subtractFlow->SetMaximum((max - min) / 2);
  subtractFlow->DrawCopy("SURF1");
  sumSummary->SetPoint(sumSummary->GetN(), sumSummary->GetN(), subtractFlow->Integral());

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
    func2->SetParLimits(2, sigmaFitLimit, 1);
    func2->SetParLimits(3, sigmaFitLimit, 1);
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
    
    if (func2->GetNDF() > 0)
    {
      Printf("#chi^{2}/ndf = %.1f/%d = %.1f", func2->GetChisquare(), func2->GetNDF(), func2->GetChisquare() / func2->GetNDF());
      DrawLatex(0.2, yPosChi2, 1, Form("#chi^{2}/ndf = %.1f/%d = %.1f", func2->GetChisquare(), func2->GetNDF(), func2->GetChisquare() / func2->GetNDF()));
    }

    canvas->cd(canvasPos++);
    TH2* residuals2 = (TH2*) hist->Clone("residuals");
    residuals2->Add(funcHist2, -1);
    residuals2->SetMinimum(-(max - min) / 2);
    residuals2->SetMaximum((max - min) / 2);
    residuals2->Draw("SURF1");
  }
  
  Float_t momentFitLimit = 0.8 - 1e-4;
  TH1* projx2 = hist->ProjectionX(Form("%s_projx2", hist->GetName()), hist->GetYaxis()->FindBin(-momentFitLimit+0.01), hist->GetYaxis()->FindBin(momentFitLimit-0.01));
  projx2->GetXaxis()->SetRangeUser(-momentFitLimit, momentFitLimit);
  
  TH1* projy2 = hist->ProjectionY(Form("%s_projy2", hist->GetName()), hist->GetXaxis()->FindBin(-momentFitLimit+0.01), hist->GetXaxis()->FindBin(momentFitLimit-0.01));
  projy2->GetXaxis()->SetRangeUser(-momentFitLimit, momentFitLimit);
  
  CalculateMomentsKurtosis(momentFitLimit, projx2, 8, graphID, x, xE);
  CalculateMomentsKurtosis(momentFitLimit, projy2, 9, graphID, x, xE);

  //   return;
  
  TH1* projx1 = subtractFlow->ProjectionX(Form("%s_projx1", hist->GetName()), hist->GetYaxis()->FindBin(-0.8), hist->GetYaxis()->FindBin(0.8));
  projx1->GetXaxis()->SetRangeUser(-momentFitLimit, momentFitLimit);

  TH1* projy1 = subtractFlow->ProjectionY(Form("%s_projy1", hist->GetName()), hist->GetXaxis()->FindBin(-0.8), hist->GetXaxis()->FindBin(0.8));
  projy1->GetXaxis()->SetRangeUser(-momentFitLimit, momentFitLimit);

  CalculateMomentsKurtosis(momentFitLimit, projx1, 8+16, graphID, x, xE);
  CalculateMomentsKurtosis(momentFitLimit, projy1, 9+16, graphID, x, xE);

//   TF1* twoGauss = new TF1("twoGauss", "gaus(0)+gaus(3)", -2, 2);
//   twoGauss->SetParameters(1, 0, 0.3, 1, 0, 0.6);
//   twoGauss->FixParameter(1, 0);
//   twoGauss->FixParameter(4, 0);
//   twoGauss->SetLineColor(4);
//   projx1->Fit("twoGauss", "I+", "SAME");
  
  canvas->cd(canvasPos++);
  projx1->Draw();
  projx1->Fit("gaus", "I");

  projx2->SetLineColor(2);
  projx2->Draw("SAME");
  projx2->Fit("gaus", "I+", "SAME");
  projx2->GetFunction("gaus")->SetLineColor(2);

  canvas->cd(canvasPos++);
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

  // set errors large for bins potentially affected by two-track effects
  if (0 && twoTrack)
  {
    Printf("NOTE : Skipping bins at (0, 0)");
//     for (Int_t binx = hist->GetXaxis()->FindBin(-0.25); binx <= hist->GetXaxis()->FindBin(0.25); binx++)
    for (Int_t binx = hist->GetXaxis()->FindBin(-0.01); binx <= hist->GetXaxis()->FindBin(0.01); binx++)
      for (Int_t biny = hist->GetYaxis()->FindBin(-0.01); biny <= hist->GetYaxis()->FindBin(0.01); biny++)
      {
// 	hist->SetBinContent(binx, biny,  0);
	hist->SetBinError(binx, biny,  1e5);
      }
  }

  // 2d fit with two gaussians
  canvas->cd(canvasPos++);
//   TF2* func3 = new TF2("func3", "[0]+[1]*exp(-0.5*((x/[2])**2+(y/[3])**2))+[4]*exp(-0.5*((x/[5])**2+(y/[6])**2))", -0.5 * TMath::Pi(), 0.5 * TMath::Pi(), -outerLimit, outerLimit);
//   func3->SetParameters(0, 1, 0.3, 0.3, 1, 0.6, 0.6);
//   func3->SetParLimits(4, 0, 10);
  Float_t etaFitLimit = outerLimit;
//   Float_t etaFitLimit = 0.5;
  TF2* func3 = new TF2("func3", "[0]+[1]*([4]/TMath::TwoPi()/[2]/[3]*exp(-0.5*((x/[2])**2+(y/[3])**2))+(1-[4])/TMath::TwoPi()/[5]/[6]*exp(-0.5*((x/[5])**2+(y/[6])**2)))", -0.5 * TMath::Pi(), 0.5 * TMath::Pi(), -etaFitLimit, etaFitLimit);
  func3->SetParameters(0, 1, 0.3, 0.3, 0.25, initSigma, initSigma);
  func3->SetParLimits(4, 0.1, 0.9);

  func3->SetParLimits(1, 0, 10);
  func3->SetParLimits(2, sigmaFitLimit, 0.6);
  func3->SetParLimits(3, sigmaFitLimit, etaFitUpperLimit);
  func3->SetParLimits(5, sigmaFitLimit, 0.6);
  func3->SetParLimits(6, sigmaFitLimit, etaFitUpperLimit);
  func3->FixParameter(0, 0);
  
  for (Int_t i=0; i<6; i++)
    func3->SetParameter(i+1, func->GetParameter(i));

  if (0 && histId == 0)
  {
    // central --> go to 1 Gaussian only
    func3->SetParLimits(4, 1, 1);
    func3->FixParameter(4, 1);
    func3->FixParameter(5, sigmaFitLimit);
    func3->FixParameter(6, sigmaFitLimit);
  }
  
  // set errors 20% of the value for bins potentially affected by two-track effects
//   hist->SetBinError(hist->GetXaxis()->FindBin(0.0001),  hist->GetYaxis()->FindBin(0.0001),  0.1 * hist->GetBinContent(hist->GetXaxis()->FindBin(0.0001),  hist->GetYaxis()->FindBin(0.0001)));
//   hist->SetBinError(hist->GetXaxis()->FindBin(0.0001),  hist->GetYaxis()->FindBin(-0.0001), 0.1 * hist->GetBinContent(hist->GetXaxis()->FindBin(0.0001),  hist->GetYaxis()->FindBin(-0.0001)));
//   hist->SetBinError(hist->GetXaxis()->FindBin(-0.0001), hist->GetYaxis()->FindBin(0.0001),  0.1 * hist->GetBinContent(hist->GetXaxis()->FindBin(-0.0001), hist->GetYaxis()->FindBin(0.0001)));
//   hist->SetBinError(hist->GetXaxis()->FindBin(-0.0001), hist->GetYaxis()->FindBin(-0.0001), 0.1 * hist->GetBinContent(hist->GetXaxis()->FindBin(-0.0001), hist->GetYaxis()->FindBin(-0.0001)));

  fitResult = hist->Fit(func3, "0R", "");
  Printf("Fit result: %d", fitResult);
  
  // if both parameters are within 1%, refit with 1 Gaussian only
  if (TMath::Abs(1.0 - func3->GetParameter(2) / func3->GetParameter(5)) < 0.01 && TMath::Abs(1.0 - func3->GetParameter(3) / func3->GetParameter(6)) < 0.01)
  {
    Printf("Parameters within 1%%. Refitting with 1 Gaussian...");
    
    func3->SetParLimits(4, 1, 1);
    func3->FixParameter(4, 1);
    func3->FixParameter(5, sigmaFitLimit);
    func3->FixParameter(6, sigmaFitLimit);
    
    fitResult = hist->Fit(func3, "0R", "");
    Printf("Fit result: %d", fitResult);
  }
  
//   hist->Fit(func3, "I0R", "");

  first = 2;
  second = 5;
  if (func3->GetParameter(2) < func3->GetParameter(5))
  {
    first = 5;
    second = 2;
  }
  //dphi
  AddPoint(graphs[1][graphID], x, TMath::Abs(func3->GetParameter(first)), xE, func3->GetParError(first));
  AddPoint(graphs[4][graphID], x, TMath::Abs(func3->GetParameter(second)), xE, func3->GetParError(second));

  AddPoint(phiWidthSummary, phiWidthSummary->GetN(), TMath::Abs(func3->GetParameter(first)), 0, func3->GetParError(first));
  AddPoint(phiWidthSummary, phiWidthSummary->GetN(), TMath::Abs(func3->GetParameter(second)), 0, func3->GetParError(second));
    
  //deta
  AddPoint(graphs[2][graphID], x, TMath::Abs(func3->GetParameter(first+1)), xE, func3->GetParError(first+1));
  AddPoint(graphs[5][graphID], x, TMath::Abs(func3->GetParameter(second+1)), xE, func3->GetParError(second+1));
  
  AddPoint(etaWidthSummary, etaWidthSummary->GetN(), TMath::Abs(func3->GetParameter(first+1)), 0, func3->GetParError(first+1));
  AddPoint(etaWidthSummary, etaWidthSummary->GetN(), TMath::Abs(func3->GetParameter(second+1)), 0, func3->GetParError(second+1));
  
  // norm
  AddPoint(graphs[0][graphID], x, func3->GetParameter(1), xE, func3->GetParError(1));
  if (first < second)
    AddPoint(graphs[3][graphID], x, func3->GetParameter(4), xE, func3->GetParError(4));
  else
    AddPoint(graphs[3][graphID], x, 1.0 - func3->GetParameter(4), xE, func3->GetParError(4));

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

  chi2 = 0;
  ndf = 0;
  for (Int_t i=hist->GetXaxis()->FindBin(-0.8); i<=hist->GetXaxis()->FindBin(0.8); i++)
    for (Int_t j=hist->GetYaxis()->FindBin(-etaLimit); j<=hist->GetYaxis()->FindBin(etaLimit); j++)
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
    AddPoint(graphs[6][graphID], x, func3->GetChisquare() / func3->GetNDF(), xE, 0);
  }
  if (ndf)
  {
    Printf("#chi^{2}/ndf = %.1f/%d = %.1f", chi2, ndf, chi2 / ndf);
    DrawLatex(0.2, yPosChi2 - 0.05, 1, Form("#chi^{2}/ndf = %.1f/%d = %.1f", chi2, ndf, chi2 / ndf));
    AddPoint(graphs[7][graphID], x, chi2 / ndf, xE, 0);
  }

  canvas->cd(canvasPos++);
  phiWidthSummary->SetMarkerStyle(20);
  phiWidthSummary->Draw("AP");
  gPad->SetGridy();
    
  etaWidthSummary->SetMarkerStyle(21);
  etaWidthSummary->SetLineColor(2);
  etaWidthSummary->SetMarkerColor(2);
  etaWidthSummary->Draw("PSAME");  
  
  phiWidthSummary->GetYaxis()->SetRangeUser(0, 0.9);
  
  canvas->cd(canvasPos++);
  sumSummary->Draw("*A");
  gPad->SetGridy();
}

void AnalyzeDeltaPhiEtaGap2D(const char* fileName, Int_t method = 1)
{
  gROOT->SetBatch(kTRUE);
  if (!gROOT->IsBatch())
  {
    Printf("Not in batch mode. Exiting!");
    return;
  }
  
  CreateGraphStructure();

  TFile::Open(fileName);
  
  Int_t maxLeadingPt = 4;
  Int_t maxAssocPt = 6;

  for (Int_t i=0; i<maxLeadingPt; i++)
  {
    for (Int_t j=1; j<maxAssocPt; j++)
    {
      // only process when first is filled
      TH2* hist1 = (TH2*) gFile->Get(Form("dphi_%d_%d_%d", i, j+1, 0));
      if (!hist1)
	continue;
      if (hist1->GetEntries() < 1e4)
      {
	Printf("Only %f entries. Skipping...", hist1->GetEntries());
	continue;
      }
      
      for (Int_t histId = 0; histId < NHists; histId++)
      {
// 	if (i != 1 || j != 2)
// 	  continue;
    
	hist1 = (TH2*) gFile->Get(Form("dphi_%d_%d_%d", i, j+1, histId));
	if (!hist1)
	  continue;
	
	if (hist1->GetEntries() < 1e4)
	{
	  Printf("Only %f entries. Skipping...", hist1->GetEntries());
	  continue;
	}
	
	TCanvas* canvas = new TCanvas(Form("DeltaPhi_%d_%d_%d_%d", histId, i, j, method), Form("DeltaPhi_%d_%d_%d_%d", histId, i, j, method), 1400, 1100);
	canvas->Divide(5, 3);
	
	Printf("\n\n>>> %d %d %d", i, j, histId);
	
	for (Int_t k=1; k<=3; k++)
	{
	  canvas->cd(3 * j + k);
	  gPad->SetLeftMargin(0.15);
	  gPad->SetBottomMargin(0.2);
	  gPad->SetTopMargin(0.01);
	  gPad->SetRightMargin(0.01);
	}
	
	Int_t graphID = i * (maxAssocPt - 1) + j - 1;

	if (histId == 0)
	  for (Int_t k=0; k<NGraphs; k++)
	    graphs[k][graphID]->SetTitle(hist1->GetTitle());
	
	Float_t centralityAxisMapping[] = { 5, 75, 100, 30, 15, 50 };
	Float_t centralityAxisMappingE[] = { 5, 15, 0, 10, 5, 10 };

	FitDeltaPhi2DOneFunction((TH2*) hist1, canvas, 1, graphID, centralityAxisMapping[histId], centralityAxisMappingE[histId], 0.9, kTRUE, histId, (i > 0) ? 1 : 0, (j <= 2 && histId != 2));
	
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

void AnalyzeDeltaPhiEtaGap2DExample(const char* fileName, Int_t i, Int_t j, Int_t histId, Bool_t drawDetails = kFALSE)
{
  CreateGraphStructure();

  TFile::Open(fileName);
  
  TH2* hist1 = (TH2*) gFile->Get(Form("dphi_%d_%d_%d", i, j+1, histId));
  if (!hist1)
    return;
  
  Printf("Entries: %f %s", hist1->GetEntries(), hist1->GetTitle());

  hist1->SetTitle("");
  hist1->GetYaxis()->SetRangeUser(-1.79, 1.79);
  hist1->GetXaxis()->SetTitleOffset(1.5);
  hist1->GetYaxis()->SetTitleOffset(2);
  hist1->SetStats(kFALSE);

  if (hist1->GetEntries() < 1e4)
    return;
  
  TCanvas* canvas = new TCanvas(Form("DeltaPhi_%d_%d_%d_%d", histId, i, j, 1), Form("DeltaPhi_%d_%d_%d_%d", histId, i, j, 1), 1400, 1100);
  canvas->Divide(5, 3);
  
  FitDeltaPhi2DOneFunction((TH2*) hist1, canvas, 1, 0, 0, 0, 0.9, kFALSE, histId, (i > 0) ? 1 : 0, (j <= 2 && histId != 2));
  
  if (!drawDetails)
    return;
  
  TVirtualPad* pad = canvas->cd(6);
  TCanvas* c = new TCanvas("c", "c", 800, 800);
  pad->SetPad(0.05, 0, 1, 1);
  pad->Draw();
  c->SaveAs("fit_subtracted.eps");

  pad = canvas->cd(12);
  c = new TCanvas("c2", "c2", 800, 800);
  pad->SetPad(0.05, 0, 1, 1);
  pad->Draw();
  c->SaveAs("fit_fit1.eps");

  pad = canvas->cd(13);
  c = new TCanvas("c3", "c3", 800, 800);
  gPad->SetLeftMargin(0.15);
  pad->SetPad(0.05, 0, 1, 1);
  pad->Draw();
  c->SaveAs("fit_residual1.eps");

  pad = canvas->cd(4);
  c = new TCanvas("c4", "c4", 800, 800);
  pad->SetPad(0.05, 0, 1, 1);
  pad->Draw();
  c->SaveAs("fit_fit2.eps");

  pad = canvas->cd(3);
  c = new TCanvas("c5", "c5", 800, 800);
  pad->SetPad(0.05, 0, 1, 1);
  pad->Draw();
  c->SaveAs("fit_residual2.eps");

  pad = canvas->cd(7);
  c = new TCanvas("c6", "c6", 800, 800);
  pad->SetPad(0.05, 0, 1, 1);
  pad->Draw();
  c->SaveAs("fit_differences.eps");
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
      else
      {
	Printf("WARNING negative value %d %d, removing from histogram", histId, i);
	graph[histId]->RemovePoint(i);
	i--;
      }
    }
  }
}

void SqrtAll2(Int_t nHists, TGraphErrors** graph, TGraphErrors** graph2)
{
  for (Int_t histId = 0; histId < nHists; histId++)
  {
    if (!graph[histId])
      continue;
 
    for (Int_t i=0; i<graph[histId]->GetN(); i++)
    {
      if (graph[histId]->GetY()[i] > 0 && graph2[histId]->GetY()[i])
      {
	graph[histId]->GetEY()[i] *= 0.5 / TMath::Sqrt(graph[histId]->GetY()[i]);
	graph[histId]->GetY()[i] = TMath::Sqrt(graph[histId]->GetY()[i]);

	graph2[histId]->GetEY()[i] *= 0.5 / TMath::Sqrt(graph2[histId]->GetY()[i]);
	graph2[histId]->GetY()[i] = TMath::Sqrt(graph2[histId]->GetY()[i]);
      }
      else
      {
	Printf("WARNING negative value %d %d, removing from histogram", histId, i);
	graph[histId]->RemovePoint(i);
	graph2[histId]->RemovePoint(i);
	i--;
      }
    }
  }
}

Int_t marker[6] = { 24, 25, 26, 30, 27, 28 };
Int_t marker2[6] = { 20, 21, 22, 29, 33, 34 };
const char* labels[6] = { "0-10%", "60-90%", "pp", "20-40%", "10-20%", "40-60%" };

void CalculateRMS(Int_t nHists, TGraphErrors** graph, TGraphErrors** graph2, TGraphErrors** weight)
{
  // rms
  for (Int_t histId = 0; histId < nHists; histId++)
  {
    if (!graph[histId])
      continue;
    
    if (graph[histId]->GetN() != graph2[histId]->GetN() || graph[histId]->GetN() != weight[histId]->GetN())
    {
      Printf("E-CalculateRMS: Different number of points for hist %d: %d %d %d", histId, graph[histId]->GetN(), graph2[histId]->GetN(), weight[histId]->GetN());
      continue;
    }
    
    for (Int_t i=0; i<graph[histId]->GetN(); i++)
    {
      Float_t rms = graph[histId]->GetY()[i] * weight[histId]->GetY()[i] + graph2[histId]->GetY()[i] * (1.0 - weight[histId]->GetY()[i]);
      graph[histId]->GetY()[i] = rms;
      // error**2 = weight**2 * error_sigma**2 + (weight-1)**2 * error_sigma2**2 + (sigma1 + sigma2)**2 * error_weight**2
      // TODO this neglects some correlations
      graph[histId]->GetEY()[i] = TMath::Sqrt(TMath::Power(weight[histId]->GetY()[i] * graph[histId]->GetEY()[i], 2) + 
	TMath::Power((weight[histId]->GetY()[i] - 1) * graph2[histId]->GetEY()[i], 2) + TMath::Power((graph[histId]->GetY()[i] + graph2[histId]->GetY()[i]) * weight[histId]->GetEY()[i], 2));
    }
  }
}

Bool_t SkipGraph(Int_t i)
{
  return (i == 1 || i == 7 || i == 9 || i == 13 || i == 17 || i == 18 || i == 19 || i == 20);
}

void PrepareGraphs(Int_t nHists, TGraphErrors** graph, TGraphErrors** systematic, TMultiGraph** multiGraph, TMultiGraph** multiGraphSyst, Int_t uncertaintyID)
{
  Int_t colors[11] =  { 1, 2, 3, 4, 6, 7, 8, 9, 10, 11 };
  Int_t markers[11] = { 20, 21, 22, 23, 24, 26, 27, 28, 30, 31 };
  Int_t fillStyle[11] = { 3001, 3002, 3003, 3004, 3005, 3006, 3007, 3008, 3009, 3010, 3011 };
  Int_t count = 0;
  
  if (*multiGraph == 0)
    *multiGraph = new TMultiGraph;
  if (*multiGraphSyst == 0)
    *multiGraphSyst = new TMultiGraph;
  
  for (Int_t i=0; i<nHists; i++)
  {
    if (SkipGraph(i))
      continue;
    
    TGraphErrors* graphcentrality = (TGraphErrors*) graph[i]->Clone();
    graphcentrality->Sort();
    if (graphcentrality->GetN() <= 0)
      continue;

    if (systematic)
    {
      TGraphErrors* graphsystematics = (TGraphErrors*) systematic[i]->Clone();
      graphsystematics->Sort();
      
      if (graphcentrality->GetN() != graphsystematics->GetN())
      {
	Printf("Different number of points %d %d", graphcentrality->GetN(), graphsystematics->GetN());
	return;
      }
      
      for (Int_t j=0; j<graphsystematics->GetN(); j++)
      {
	// uncertaintyID
	Double_t yMin = graphcentrality->GetY()[j];
	Double_t yMax = graphcentrality->GetY()[j];
	
	if (uncertaintyID == 1)
	{
	  // 10%
	  yMin *= 0.9;
	  yMax *= 1.1;
	}
	else if (uncertaintyID == 2)
	{
	  yMin -= 0.15;
	  yMax += 0.15;
	}
	
	yMin = TMath::Min(yMin, graphsystematics->GetY()[j]);
	yMax = TMath::Max(yMax, graphsystematics->GetY()[j]);
	
	graphsystematics->GetEY()[j] = TMath::Abs(yMin - yMax) / 2;
	graphsystematics->GetY()[j]  = (yMin + yMax) / 2;
	graphsystematics->GetEX()[j] = 1;
      }
      
//       graphsystematics->SetFillColor(kGray);
      graphsystematics->SetFillColor(colors[count]);
//       graphsystematics->SetFillStyle(1001);
      graphsystematics->SetFillStyle(fillStyle[count]);
      graphsystematics->SetMarkerStyle(0);
      graphsystematics->SetLineColor(0);
      (*multiGraphSyst)->Add(graphsystematics, "2");
    }
    
    graphcentrality->SetMarkerStyle(markers[count]);
    graphcentrality->SetMarkerColor(colors[count]);
    graphcentrality->SetLineColor(colors[count]);
    
    (*multiGraph)->Add(graphcentrality);

    TString label = graphcentrality->GetTitle();
    if (label.Length() > 0)
    {
      TObjArray* tokens = label.Tokenize("-");
      label.Form("%s-%s", tokens->At(0)->GetName(),tokens->At(1)->GetName());
      label.ReplaceAll(".00", "");
      label.ReplaceAll(".0", "");
    }
    graphcentrality->SetTitle(label);
    Printf("%d %s %d", i, label.Data(), colors[count]);
    
    count++;
  }
}

void DrawCentrality(const char* canvasName, Int_t nHists, TGraphErrors** graph, Float_t min = 0, Float_t max = 0, const char* yLabel = "", TGraphErrors** systematic = 0, TGraphErrors** graph2 = 0, TGraphErrors** systematic2 = 0, Int_t uncertaintyID = -1)
{
  Bool_t found = kTRUE;
  TCanvas* c1 = (TCanvas*) gROOT->GetListOfCanvases()->FindObject(canvasName);
  if (!c1)
  {
    c1 = new TCanvas(canvasName, canvasName, 800, 600);
    found = kFALSE;
  }
  c1->cd();
  
  TMultiGraph* multiGraph = 0;
  TMultiGraph* multiGraphSyst = 0;
  
  PrepareGraphs(nHists, graph, systematic, &multiGraph, &multiGraphSyst, uncertaintyID);
  TLegend* legend = new TLegend(0.65, 0.65, 0.95, 0.95);
  legend->SetFillColor(0);
  for (Int_t i=0; i<multiGraph->GetListOfGraphs()->GetEntries(); i++)
    legend->AddEntry(multiGraph->GetListOfGraphs()->At(i), 0, "PL");

  if (graph2)
  {
    TMultiGraph* multiGraph2 = 0;
    PrepareGraphs(nHists, graph2, systematic2, &multiGraph2, &multiGraphSyst, uncertaintyID);
    for (Int_t i=0; i<multiGraph2->GetListOfGraphs()->GetEntries(); i++)
      ((TGraph*)multiGraph2->GetListOfGraphs()->At(i))->SetLineWidth(2);
    
    while (multiGraph2->GetListOfGraphs()->GetEntries() > multiGraph->GetListOfGraphs()->GetEntries())
      multiGraph2->GetListOfGraphs()->RemoveAt(multiGraph->GetListOfGraphs()->GetEntries());
      
    multiGraphSyst->Add(multiGraph2, "L");
  }

  multiGraphSyst->Add(multiGraph, (found) ? "LX" : "P");
  
  TString drawString("A");
  multiGraphSyst->Draw(drawString);
  multiGraphSyst->GetXaxis()->SetTitle("Centrality");
  multiGraphSyst->GetYaxis()->SetTitle(yLabel);
  if (max > min)
    multiGraphSyst->GetYaxis()->SetRangeUser(min, max);
    
  legend->Draw();
  
  gPad->SetGridx();
  gPad->SetGridy();
//   c1->SaveAs(Form("results/%s.png", canvasName));
//   c1->SaveAs(Form("%s.eps", canvasName));
}

void CalculateRMSSigma()
{
  CalculateRMS(NHists, graphs[1], graphs[4], graphs[3]);
  CalculateRMS(NHists, graphs[2], graphs[5], graphs[3]);

  CalculateRMS(NHists, graphs[1+16], graphs[4+16], graphs[3+16]);
  CalculateRMS(NHists, graphs[2+16], graphs[5+16], graphs[3+16]);

  // sqrt(moment2) = sigma
  SqrtAll2(NHists, graphs[8], graphs[8+16]);
  SqrtAll2(NHists, graphs[9], graphs[9+16]);
}

void DrawResultsCentrality(const char* fileName = "graphs.root")
{
  ReadGraphs(fileName);
  
  Int_t nHists = NHists;

  if (1)
  {
    DrawCentrality("norm", nHists, graphs[0], 0, 0.2, "N (a.u.)", graphs[0+16], 0, 0, 0);
    
//     return;
    
    DrawCentrality("width_phi1_centrality", nHists, graphs[1], 0, 0.8, "#sigma_{#phi, 1} (rad.)", graphs[1+16]);
//     return;
    DrawCentrality("width_phi2_centrality", nHists, graphs[4], 0, 0.8, "#sigma_{#phi, 2} (rad.)", graphs[4+16]);
    DrawCentrality("width_eta1_centrality", nHists, graphs[2], 0, 0.8, "#sigma_{#eta, 1} (rad.)", graphs[2+16]);
    DrawCentrality("width_eta2_centrality", nHists, graphs[5], 0, 0.8, "#sigma_{#eta, 2} (rad.)", graphs[5+16]);
    
    CalculateRMSSigma();

    DrawCentrality("phi_rms_centrality", nHists, graphs[1], 0, 0.8, "#sigma_{#phi} (fit) (rad.)", graphs[1+16], 0, 0, 1);
    DrawCentrality("eta_rms_centrality", nHists, graphs[2], 0, 0.8, "#sigma_{#eta} (fit)", graphs[2+16], 0, 0, 1);

    DrawCentrality("chi2_1", nHists, graphs[6], 0.5, 2.5, "#chi^{2}/ndf (full region)");
    DrawCentrality("chi2_2", nHists, graphs[7], 0.5, 2.5, "#chi^{2}/ndf (peak region)");
    
    DrawCentrality("sigma_phi", nHists, graphs[8], 0, 0.8, "#sigma_{#phi} (rad.)", graphs[8+16], 0, 0, 1);
    DrawCentrality("sigma_eta", nHists, graphs[9], 0, 0.8, "#sigma_{#eta}", graphs[9+16], 0, 0, 1);
  }
  
  DrawCentrality("kurtosisphi_centrality", 12, graphs[14], -2, 4, "Kurtosis #phi", graphs[14+16], 0, 0, 2);
  DrawCentrality("kurtosiseta_centrality", 12, graphs[15], -2, 4, "Kurtosis #eta", graphs[15+16], 0, 0, 2);
}

void MCComparison(const char* fileNameData, const char* fileNameHijing)
{
  Int_t nHists = 12; //NHists;

  ReadGraphs(fileNameHijing);
  if (0)
  {
    Int_t n=0;
    while (n++ < graphs[0][0]->GetN())
    {
  //     Printf("%d %f", n, graphs[0][0]->GetX()[n]); continue;
      if (graphs[0][0]->GetX()[n] == 42)
      {
	for (Int_t i=0; i<NGraphs; i++)
	  graphs[i][0]->RemovePoint(n);
	break;
      }
    }
  }
//   return;
  
  CalculateRMSSigma();
  
//   DrawCentrality("phi_rms_centrality", nHists, graphs[1], 0, 0.8, "#sigma_{#phi} (fit) (rad.)");
//   return;

  TGraphErrors*** graphs1 = graphs;

  ReadGraphs(fileNameData);

  //Remove high trigger pT:18.26,34,42
  if (0)
  {
    Int_t n=0;
    while (n++ < graphs[0][0]->GetN())
    {
  //     Printf("%d %f", n, graphs[0][0]->GetX()[n]); continue;
      if (graphs[0][0]->GetX()[n] == 18 || graphs[0][0]->GetX()[n] == 10) //  || graphs[0][0]->GetX()[n] == 34 || graphs[0][0]->GetX()[n] == 42
      {
	for (Int_t i=0; i<NGraphs; i++)
	{
	  graphs[i][0]->RemovePoint(n);
	}
	n--;
      }
    }
  }

  CalculateRMSSigma();
  
  DrawCentrality("phi_rms_centrality_mc", nHists, graphs[1], 0, 0.8, "#sigma_{#phi} (fit) (rad.)", graphs[1+16], graphs1[1], 0, 1);
  DrawCentrality("eta_rms_centrality_mc", nHists, graphs[2], 0, 0.8, "#sigma_{#eta} (fit)", graphs[2+16], graphs1[2], 0, 1);

  DrawCentrality("sigma_phi", nHists, graphs[8], 0, 0.8, "#sigma_{#phi} (rad.)", graphs[8+16], graphs1[8], 0, 1);
  DrawCentrality("sigma_eta", nHists, graphs[9], 0, 0.8, "#sigma_{#eta}", graphs[9+16], graphs1[9], 0, 1);

  DrawCentrality("kurtosisphi_centrality_mc", nHists, graphs[14], -2, 4, "Kurtosis #phi", graphs[14+16], graphs1[14], 0, 2);
  DrawCentrality("kurtosiseta_centrality_mc", nHists, graphs[15], -2, 4, "Kurtosis #eta", graphs[15+16], graphs1[15], 0, 2);
}

Float_t** ExtractSystematics(const char* baseFile, const char* systFile)
{
  ReadGraphs(baseFile);
  CalculateRMSSigma();
 
  TGraphErrors*** graphsBase = graphs;
  if (systFile)
  {
    ReadGraphs(systFile);
    CalculateRMSSigma();
  }
  
  // calculate syst unc for these graphs
  const Int_t NGraphList = 6;
  Int_t graphList[] = { 1, 2, 8, 9, 14, 15 };

  Float_t** results = new Float_t*[NGraphList];

  for (Int_t i=0; i<NGraphList; i++)
  {
    results[i] = new Float_t[2];

    Bool_t kurtosis = kFALSE;
    if (graphList[i] == 14 || graphList[i] == 15)
      kurtosis = kTRUE;
    
    TH1* hist = 0;
    if (!kurtosis)
      hist = new TH1F(Form("hist_%d", i), "", 50, 0.5, 1.5);  
    else
      hist = new TH1F(Form("hist_%d", i), "", 50, -1, 1);  
  
    TCanvas* c = new TCanvas(Form("%d_%d", i, 0), Form("%d_%d", i, 0), 1000, 1000);
    c->Divide(4, 4);

    Int_t count = 1;
    for (Int_t j=0; j<NHists; j++)
    {
      if (SkipGraph(j))
	continue;
      
      // for kurtosis we show only up to 12
      if (kurtosis)
	if (j >= 12)
	  continue;
      
      TGraphErrors* graph1 = graphsBase[graphList[i]][j];
      TGraphErrors* graph2 = (systFile) ? graphs[graphList[i]][j] : graphsBase[graphList[i]+16][j];
      
      if (graph1->GetN() == 0)
	continue;
      
      graph1->Sort();
      graph2->Sort();
      
      c->cd(count++);
      gPad->SetGridx();
      gPad->SetGridy();
      graph1->SetMarkerStyle(24);
      graph1->DrawClone("AP");
      graph2->SetLineColor(2);
      graph2->SetMarkerColor(2);
      graph2->SetMarkerStyle(25);
      graph2->DrawClone("LSAME");
      
      c->cd(count++);
      gPad->SetGridx();
      gPad->SetGridy();
      
      // reset errors on graph2 
/*      for (Int_t k=0; k<graph2->GetN(); k++)
	graph2->GetEY()[k] = 0;*/
      
      if (!kurtosis)
      {
	DivideGraphs(graph1, graph2);
	((TGraphErrors*) graph1->DrawClone("AP"))->GetYaxis()->SetRangeUser(0.8, 1.2);
      }
      else
      {
	SubtractGraphs(graph1, graph2);
	((TGraphErrors*) graph1->DrawClone("AP"));
      }	
      
      for (Int_t k=0; k<graph1->GetN(); k++)
// 	hist->Fill(graph1->GetY()[k], 1.0 / (graph1->GetEY()[k] / graph1->GetY()[k]));
	hist->Fill(graph1->GetY()[k]);
      
      if (count == 37)
	break;
    }
    
    new TCanvas;
    hist->Draw();
    hist->Sumw2();
    
    hist->Fit("gaus", "");
    Float_t mean = hist->GetFunction("gaus")->GetParameter(1);
    Float_t sigma = hist->GetFunction("gaus")->GetParameter(2);
    
    Printf("%d: %.2f %.2f", i, mean, sigma);
    if (!kurtosis)
    {
      mean -= 1;
      mean *= 100;
      sigma *= 100;
    }
    results[i][0] = mean;
    results[i][1] = sigma;
  }
  
  return results;
}

void ExtractSystematicsAll()
{
  gROOT->SetBatch(kTRUE);
  
  const Int_t NEffects = 6;
  
  const char* defaultFile = "graphs_120425.root";
  const char* systFiles[] = { "graphs_120425_eta08.root", "graphs_120425_eta12.root", "graphs_120425_outer14.root", "graphs_120425_wingremoved.root", "graphs_hybrid_120427.root", "graphs_120425_vertex.root" };

  Float_t** results[NEffects];
  
  for (Int_t i=0; i<NEffects; i++)
  {
    results[i] = ExtractSystematics(defaultFile, systFiles[i]);
  }
  
  const char* names[] = { "$\\sigma_{\\Dphi}$ (fit)", "$\\sigma_{\\Deta}$ (fit)", "$\\sigma_{\\Dphi}$", "$\\sigma_{\\Deta}$", "Kurtosis $\\Dphi$", "Kurtosis $\\Deta$" };
  
  // put together 0-2 (into 2)
  for (Int_t j=0; j<6; j++)
  {
    Float_t mean = 0;
    Float_t sigma = 0;
    
    printf("%s \t:", names[j]);
    for (Int_t i=0; i<=2; i++)
    {
      mean = TMath::Max(mean, TMath::Abs(results[i][j][0]));
      printf("%.1f%% ", results[i][j][0]);
    }
    printf("--> %.1f%% \t", mean);
    results[2][j][0] = mean;
      
    for (Int_t i=0; i<=2; i++)
    {
      sigma = TMath::Max(sigma, TMath::Abs(results[i][j][1]));
      printf("%.1f%% ", results[i][j][1]);
    }
    Printf("--> %.1f%% \t", sigma);
    results[2][j][1] = sigma;
  }
  
  for (Int_t j=0; j<6; j++)
  {
    printf("%s \t & ", names[j]);
    for (Int_t i=2; i<NEffects; i++)
    {
      if (j < 4)
	printf("$%.1f\\%% \\pm %.1f\\%%$ \t & ", TMath::Abs(results[i][j][0]), results[i][j][1]);
      else
	printf("$%.2f \\pm %.2f$ \t & ", TMath::Abs(results[i][j][0]), results[i][j][1]);
    }
    Printf("\\\\");
  }
  
  // todo fit method
}

void CompareEtaPhi(const char* graphFileName1)
{
  ReadGraphs(graphFileName1);

  Int_t nHists = 6;

  CalculateRMSSigma();

  DrawCentrality("rms_centrality", nHists, graphs[1], 0, 0.8, "#sigma_{#phi} (fit) (rad.) / #sigma_{#eta}", graphs[1+16], graphs[2], graphs[2+16], 1);
  DrawCentrality("rms_centrality_nosyst", nHists, graphs[1], 0, 0.8, "#sigma_{#phi} (fit) (rad.) / #sigma_{#eta}", 0, graphs[2], 0);
}

void DrawExamples(const char* histFileName, const char* graphFileName, Bool_t drawFunc = kTRUE)
{
  Float_t etaLimit = 1.0;
  Float_t outerLimit = 1.79;
  Float_t projectLimit = 0.8;

  ReadGraphs(graphFileName);

  TFile::Open(histFileName);
  
  Int_t j=1;

  Int_t exColors[] = { 1, 2, 3, 4, 5, 6 };
  
//   nHists = 2;
  
  for (Int_t i=0; i<1; i++)
  {
    Int_t graphID = i * (6 - 1) + j - 1;

    TCanvas* c = new TCanvas(Form("c_%d", i), Form("c_%d", i), 1000, 600);
    c->Divide(2, 1);
    Int_t nHists = 6;
    for (Int_t histId = 0; histId < nHists; histId++)
    {
      if (histId == 2)
	histId = 5;
      TH2* hist = (TH2*) gFile->Get(Form("dphi_%d_%d_%d", i, j+1, histId));
      if (!hist)
	continue;

      if (graphs[0][graphID]->GetN() < histId)
      {
	Printf("ERROR: Pos in graph not found: %d %d", i, histId);
	continue;
      }
      
      SubtractEtaGap(hist, etaLimit, outerLimit, kFALSE, kFALSE);
    
      c->cd(1);
      TH1* proj = hist->ProjectionX(Form("%s_proj1", hist->GetName()), hist->GetYaxis()->FindBin(-projectLimit+0.01), hist->GetYaxis()->FindBin(projectLimit-0.01));
      proj->SetLineColor(exColors[histId]);
      proj->GetXaxis()->SetRangeUser(-1.49, 1.49);
      proj->GetYaxis()->SetRangeUser(proj->GetMinimum() * 1.2, proj->GetMaximum() * 1.5);
      proj->Draw((histId == 0) ? "" : "SAME");
    
      if (drawFunc)
      {
	// integral over y
	TF1* func = new TF1("func", "[0]+[1]*([4]/[2]/TMath::Sqrt(TMath::TwoPi())*exp(-0.5*(x/[2])**2)+(1-[4])/[5]/TMath::Sqrt(TMath::TwoPi())*exp(-0.5*(x/[5])**2))", -0.5 * TMath::Pi(), 0.5 * TMath::Pi());
	func->SetParameter(0, 0);
	for (Int_t k=0; k<6; k++)
	  func->SetParameter(k+1, graphs[k][graphID]->GetY()[histId]);
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
	TF1* func = new TF1("func", "[0]+[1]*([4]/[3]/TMath::Sqrt(TMath::TwoPi())*exp(-0.5*(x/[3])**2)+(1-[4])/[6]/TMath::Sqrt(TMath::TwoPi())*exp(-0.5*(x/[6])**2))", -outerLimit, outerLimit);
	func->SetParameter(0, 0);
	for (Int_t k=0; k<6; k++)
	  func->SetParameter(k+1, graphs[k][graphID]->GetY()[histId]);
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

void DrawExample(const char* histFileName, Int_t i, Int_t j, Int_t histId, TH1** projPhi, TH1** projEta)
{
  Float_t etaLimit = 1.0;
  Float_t outerLimit = 1.6;
  Float_t projectLimit = 0.8;

  TFile::Open(histFileName);
  
  TCanvas* c = new TCanvas(Form("ex_%d_%d_%d", i, j, histId), Form("ex_%d_%d_%d", i, j, histId), 1200, 400);
  c->Divide(3, 1);

  TH2* hist = (TH2*) gFile->Get(Form("dphi_%d_%d_%d", i, j+1, histId));
  if (!hist)
    return;
  
  TString label(hist->GetTitle());
  label.ReplaceAll(".00", ".0");
  TObjArray* objArray = label.Tokenize("-");
  TPaveText* paveText = new TPaveText(0.52, 0.72, 0.95, 0.9, "BRNDC");
  paveText->SetTextSize(0.04);
  paveText->SetFillColor(0);
  paveText->SetShadowColor(0);
  paveText->AddText(objArray->At(0)->GetName());
  paveText->AddText(objArray->At(1)->GetName());
  if (objArray->GetEntries() == 4)
    paveText->AddText(Form("Pb-Pb %s-%s", objArray->At(2)->GetName(), objArray->At(3)->GetName()));
  else
    paveText->AddText(objArray->At(2)->GetName());
  
  c->cd(1);
  hist->GetXaxis()->SetRangeUser(-TMath::Pi() / 2, TMath::Pi() / 2);
  hist->GetYaxis()->SetRangeUser(-outerLimit, outerLimit);
  hist->GetXaxis()->SetTitleOffset(1.5);
  hist->GetYaxis()->SetTitleOffset(1.5);
  hist->SetStats(0);
  hist->SetTitle("a) Correlation");
  hist->DrawCopy("SURF1");
  paveText->Draw();
  
  SubtractEtaGap(hist, etaLimit, outerLimit, kFALSE, kFALSE);
  c->cd(2);
  hist->SetTitle("b) #eta-gap subtracted");
  hist->DrawCopy("SURF1");
    
  c->cd(3);
  TH1* proj = hist->ProjectionX(Form("%s_proj1", hist->GetName()), hist->GetYaxis()->FindBin(-projectLimit+0.01), hist->GetYaxis()->FindBin(projectLimit-0.01));
  TH1* proj2 = hist->ProjectionY(Form("%s_proj2b", hist->GetName()), hist->GetXaxis()->FindBin(-projectLimit+0.01), hist->GetXaxis()->FindBin(projectLimit-0.01));

  proj->SetStats(0);
  proj->SetTitle("c) Projections");
  proj->GetXaxis()->SetRangeUser(-1.49, 1.49);
  proj->GetXaxis()->SetTitleOffset(1);
  proj->GetYaxis()->SetRangeUser(proj->GetMinimum(), proj->GetMaximum() * 1.4);
  TH1* copy = proj->DrawCopy();
  copy->GetXaxis()->SetTitle(Form("%s / %s", proj->GetXaxis()->GetTitle(), proj2->GetXaxis()->GetTitle()));
  DrawLatex(0.3, 0.85, 1, Form("#Delta#phi projection in |#Delta#eta| < %.1f", projectLimit), 0.04);
    
  proj2->SetLineColor(2);
  proj2->SetStats(0);
  proj2->GetXaxis()->SetTitleOffset(1);
  proj2->GetYaxis()->SetRangeUser(proj2->GetMinimum(), proj2->GetMaximum() * 1.6);
  proj2->DrawCopy("SAME");
  DrawLatex(0.3, 0.80, 2, Form("#Delta#eta projection in |#Delta#phi| < %.1f", projectLimit), 0.04);
  
  proj->SetTitle(label);
  proj2->SetTitle(label);
  
  *projPhi = proj;
  *projEta = proj2;

  c->SaveAs(Form("ex/%s.png", c->GetName()));
  c->SaveAs(Form("ex/%s.eps", c->GetName()));

  if (0)
  {
    c = new TCanvas(Form("ex_%d_%d_%d_a", i, j, histId), Form("ex_%d_%d_%d_a", i, j, histId), 400, 400);
    hist->SetTitle("");
    hist->DrawCopy("SURF1");
    paveText->Draw();
    c->SaveAs(Form("ex/%s.png", c->GetName()));
    c->SaveAs(Form("ex/%s.eps", c->GetName()));
  }
}

void DrawExampleAll(const char* histFileName)
{
  Int_t colors[] = { 1, 2, 4 };
  
  Float_t exampleI[] = { 0, 1, 2 };
  Float_t exampleJ[] = { 1, 2, 3 };
  
  TH1* projectionsPhi[9];
  TH1* projectionsEta[9];
  
  Int_t count = 0;
  for (Int_t i=0; i<3; i++)
  {
    for (Int_t histId = 0; histId<3; histId++)
    {
      DrawExample(histFileName, exampleI[i], exampleJ[i], histId, &projectionsPhi[count], &projectionsEta[count]);
      count++;
    }
    
    TCanvas* c = new TCanvas(Form("centralities_%d", i), Form("centralities_%d", i), 800, 400);
    c->Divide(2, 1);
    
    TLegend* legend = new TLegend(0.62, 0.7, 0.88, 0.88);
    legend->SetFillColor(0);
    
    for (Int_t histId = 0; histId<3; histId++)
    {
      c->cd(1);
      TH1* clone = projectionsPhi[count-3+histId]->DrawCopy((histId > 0) ? "SAME" : "");
      clone->SetLineColor(colors[histId]);
      
      TString label(clone->GetTitle());
      TObjArray* objArray = label.Tokenize("-");
      clone->SetTitle(Form("%s-%s", objArray->At(0)->GetName(), objArray->At(1)->GetName()));
      
      legend->AddEntry(clone, (objArray->GetEntries() == 4) ? Form("%s-%s", objArray->At(2)->GetName(), objArray->At(3)->GetName()) : objArray->At(2)->GetName(), "L");

      c->cd(2);
      clone = projectionsEta[count-3+histId]->DrawCopy((histId > 0) ? "SAME" : "");
      clone->SetTitle(Form("%s-%s", objArray->At(0)->GetName(), objArray->At(1)->GetName()));
      clone->SetLineColor(colors[histId]);
    }
    
    c->cd(1);
    legend->Draw();
    c->SaveAs(Form("ex/%s.png", c->GetName()));
    c->SaveAs(Form("ex/%s.eps", c->GetName()));
  }
  
  for (Int_t histId = 0; histId<3; histId++)
  {
    TCanvas* c = new TCanvas(Form("pt_%d", histId), Form("pt_%d", histId), 800, 400);
    c->Divide(2, 1);
    
    TLegend* legend = new TLegend(0.15, 0.7, 0.88, 0.88);
    legend->SetFillColor(0);
    
    for (Int_t i=2; i>=0; i--)
    {
      c->cd(1);
      TH1* clone = projectionsPhi[i*3+histId]->DrawCopy((i < 2) ? "SAME" : "");
      clone->SetLineColor(colors[i]);
      clone->GetYaxis()->SetRangeUser(clone->GetMinimum(), clone->GetMaximum() * 1.2);
      
      TString label(clone->GetTitle());
      TObjArray* objArray = label.Tokenize("-");
      clone->SetTitle((objArray->GetEntries() == 4) ? Form("%s-%s", objArray->At(2)->GetName(), objArray->At(3)->GetName()) : objArray->At(2)->GetName());
      
      legend->GetListOfPrimitives()->AddFirst(new TLegendEntry(clone, Form("%s-%s", objArray->At(0)->GetName(), objArray->At(1)->GetName()), "L"));
      
      c->cd(2);
      clone = projectionsEta[i*3+histId]->DrawCopy((i < 2) ? "SAME" : "");
      clone->SetTitle((objArray->GetEntries() == 4) ? Form("%s-%s", objArray->At(2)->GetName(), objArray->At(3)->GetName()) : objArray->At(2)->GetName());
      clone->SetLineColor(colors[i]);
//       clone->GetYaxis()->SetRangeUser(clone->GetMinimum(), clone->GetMaximum() * 1.2);
    }
    
    c->cd(1);
    legend->Draw();

    c->SaveAs(Form("ex/%s.png", c->GetName()));
    c->SaveAs(Form("ex/%s.eps", c->GetName()));
  }
}

void DrawDoubleHump(const char* histFileName)
{
  Float_t exampleI[] = { 0, 0, 0, 1, 1, 1};
  Float_t exampleJ[] = { 0, 1, 2, 0, 1, 2};
  
  TH1* projectionsPhi[9];
  TH1* projectionsEta[9];
  
  TCanvas* c = new TCanvas("DrawDoubleHump", "DrawDoubleHump", 1200, 800);
  c->Divide(3, 2);
  
  TLegend* legend = new TLegend(0.15, 0.7, 0.88, 0.88);
  legend->SetFillColor(0);
  
  for (Int_t i=0; i<6; i++)
  {
    DrawExample(histFileName, exampleI[i], exampleJ[i], 0, &projectionsPhi[i], &projectionsEta[i]);

    c->cd(i+1);
    TH1* clone = projectionsPhi[i]->DrawCopy("");
    clone->SetLineColor(1);
    clone->GetYaxis()->SetRangeUser(clone->GetMinimum(), clone->GetMaximum() * 1.2);
//     clone->GetXaxis()->SetTitle(Form("%s / %s", clone->GetXaxis()->GetTitle(), "#Delta#eta"));
    clone->GetXaxis()->SetTitle(Form("%s / %s", "#Delta#phi (rad.)", "#Delta#eta"));
    
    clone = projectionsEta[i]->DrawCopy("SAME");
    clone->SetLineColor(2);

    DrawLatex(0.3, 0.85, 1, Form("#Delta#phi projection in |#Delta#eta| < %.1f", 0.8), 0.04);
    DrawLatex(0.3, 0.80, 2, Form("#Delta#eta projection in |#Delta#phi| < %.1f", 0.8), 0.04);
  }
}

void DrawFullCentralityDependence(const char* histFileName)
{
  Int_t colors[] = { 1, 2, 4, 3, 5, 6 };
  
  TH1* projectionsPhi[9];
  TH1* projectionsEta[9];
  
  Int_t count = 0;
  for (Int_t histId = 0; histId<6; histId++)
  {
    DrawExample(histFileName, 0, 1, histId, &projectionsPhi[count], &projectionsEta[count]);
    count++;
  }
  
  TCanvas* c = new TCanvas(Form("centralities_%d", 0), Form("centralities_%d", 0), 800, 400);
  c->Divide(2, 1);
  
  TLegend* legend = new TLegend(0.62, 0.7, 0.88, 0.88);
  legend->SetFillColor(0);
  
  for (Int_t histId = 0; histId<6; histId++)
  {
    c->cd(1);
    TH1* clone = projectionsPhi[count-6+histId]->DrawCopy((histId > 0) ? "SAME" : "");
    clone->SetLineColor(colors[histId]);
    
    TString label(clone->GetTitle());
    TObjArray* objArray = label.Tokenize("-");
    clone->SetTitle(Form("%s-%s", objArray->At(0)->GetName(), objArray->At(1)->GetName()));
    
    legend->AddEntry(clone, (objArray->GetEntries() == 4) ? Form("%s-%s", objArray->At(2)->GetName(), objArray->At(3)->GetName()) : objArray->At(2)->GetName(), "L");

    c->cd(2);
    clone = projectionsEta[count-6+histId]->DrawCopy((histId > 0) ? "SAME" : "");
    clone->SetTitle(Form("%s-%s", objArray->At(0)->GetName(), objArray->At(1)->GetName()));
    clone->SetLineColor(colors[histId]);
  }
    
  c->cd(1);
  legend->Draw();
  c->SaveAs(Form("ex/%s.png", c->GetName()));
  c->SaveAs(Form("ex/%s.eps", c->GetName()));
}

void CompareExamples(const char* histFileName1, const char* histFileName2, Int_t i, Int_t j, Int_t histId, Bool_t twoTrack = kFALSE)
{
  Float_t etaLimit = 1.0;
  Float_t outerLimit = 1.79;
  Float_t projectLimit = 0.8;

  Int_t exColors[] = { 1, 2, 4, 6 };
  
  TCanvas* c = new TCanvas(Form("c_%d", i), Form("c_%d", i), 1000, 600);
  c->Divide(2, 1);
  
  TCanvas* c2 = new TCanvas(Form("c_%d_ratio", i), Form("c_%d_ratio", i), 1000, 600);
  c2->Divide(2, 1);
  
  TH1* projFirst[2];
  
  for (Int_t k=0; k<2; k++)
  {
    if (k == 0)
      TFile::Open(histFileName1);
    else
      TFile::Open(histFileName2);
  
    TH2* hist = (TH2*) gFile->Get(Form("dphi_%d_%d_%d", i, j+1, histId));
    if (!hist)
      continue;
//     if (k == 1)
//       hist->Scale(2);
//     hist->Rebin2D(2, 2);

//     new TCanvas; hist->DrawCopy("COLZ"); return;

    SubtractEtaGap(hist, etaLimit, outerLimit, kFALSE);
    
  //     new TCanvas; hist->DrawCopy("COLZ");

    // remove bins potentially affected by two-track effects
    if (twoTrack)
    {
      Printf("NOTE : Skipping bins at (0, 0)");
      for (Int_t binx = hist->GetXaxis()->FindBin(-0.25); binx <= hist->GetXaxis()->FindBin(0.25); binx++)
	for (Int_t biny = hist->GetYaxis()->FindBin(-0.01); biny <= hist->GetYaxis()->FindBin(0.01); biny++)
	{
	  hist->SetBinContent(binx, biny,  0);
	  hist->SetBinError(binx, biny, 0);
	}
    }
  
    c->cd(1);
    TH1* proj = hist->ProjectionX(Form("%s_proj1%d", hist->GetName(), k), hist->GetYaxis()->FindBin(-projectLimit+0.01), hist->GetYaxis()->FindBin(projectLimit-0.01));
    proj->SetLineColor(exColors[k]);
    proj->GetXaxis()->SetRangeUser(-1.49, 1.49);
    proj->GetYaxis()->SetRangeUser(proj->GetMinimum() * 1.2, proj->GetMaximum() * 1.5);
    proj->DrawCopy((k == 0) ? "" : "SAME");
  
    if (k == 1)
    {
      c2->cd(1);
      proj->Divide(projFirst[0]);
      proj->Draw();
      proj->GetYaxis()->SetRangeUser(0.5, 1.5);
    }
    else
      projFirst[0] = proj;

    c->cd(2);
    proj = hist->ProjectionY(Form("%s_proj2_%d", hist->GetName(), k), hist->GetXaxis()->FindBin(-projectLimit+0.01), hist->GetXaxis()->FindBin(projectLimit-0.01));
    proj->SetLineColor(exColors[k]);
    proj->GetXaxis()->SetRangeUser(-1.49, 1.49);
    proj->GetYaxis()->SetRangeUser(proj->GetMinimum() * 1.2, proj->GetMaximum() * 2);
    proj->DrawCopy((k == 0) ? "" : "SAME");

    if (k == 1)
    {
      c2->cd(2);
      proj->Divide(projFirst[1]);
      proj->Draw();
      proj->GetYaxis()->SetRangeUser(0.5, 1.5);
    }
    else
      projFirst[1] = proj;
  }
}

void TestMomentCode()
{
  Float_t momentLimit = 1;
  TH1* proj = new TH1F("hist", "", 100, -momentLimit, momentLimit);
  TF1* func = new TF1("func", "gaus(0)", -momentLimit, momentLimit);
  func->SetParameters(1, 0, 0.2);
  proj->FillRandom("func", 100000);
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
    
    Printf("%d %f +- %f <-> %f +- %f", n, moment, TMath::Sqrt(error), proj->GetRMS() * proj->GetRMS(), 2 * proj->GetRMSError() / proj->GetRMS() * proj->GetRMS() * proj->GetRMS());
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

void TestTwoGaussian()
{
//   TF1* func = new TF1("func", "[0]/TMath::Sqrt(TMath::TwoPi())/[2]*exp(-0.5*((x-[1])/[2])**2) + [3]/TMath::Sqrt(TMath::TwoPi())/[5]*exp(-0.5*((x-[4])/[5])**2)", -2, 2);
  TF1* func = new TF1("func", "[0]*exp(-0.5*((x-[1])/[2])**2) + [3]*exp(-0.5*((x-[4])/[5])**2)", -2, 2);
  func->SetParameters(0.25, 0, 0.5, 0.75, 0, 0.2);
  
  TH1* hist = new TH1F("hist", "", 100, -2, 2);
  hist->FillRandom("func", 10000000);
  
  hist->Draw();
  
  TH1* proj = hist;
  
  Float_t momentLimit = 1.99;
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

void DrawEtaGapExample(const char* histFileName, Int_t i = 1, Int_t j = 2, Int_t histId = 0)
{
  Float_t etaLimit = 1.0;
  Float_t outerLimit = 1.59;
  Float_t projectLimit = 0.8;

  TFile::Open(histFileName);

  TH2* hist = (TH2*) gFile->Get(Form("dphi_%d_%d_%d", i, j, histId));
  if (!hist)
    return;
  
  TCanvas* c = new TCanvas("c", "c", 800, 800);
  gPad->SetLeftMargin(0.15);
  TH2* clone = (TH2*) hist->DrawCopy("SURF1");
  clone->SetTitle("");
  clone->GetYaxis()->SetRangeUser(-1.79, 1.79);
  clone->GetXaxis()->SetTitleOffset(1.5);
  clone->GetYaxis()->SetTitleOffset(2);
  clone->SetStats(kFALSE);
  c->SaveAs("etagap_raw.eps");
  
  SubtractEtaGap(hist, etaLimit, outerLimit, kFALSE, kTRUE);
  
  c = new TCanvas("c3", "c3", 800, 800);
  TH1* proj = hist->ProjectionX(Form("%s_proj1%d", hist->GetName(), 0), hist->GetYaxis()->FindBin(-projectLimit+0.01), hist->GetYaxis()->FindBin(projectLimit-0.01));
  proj->GetXaxis()->SetRangeUser(-1.49, 1.49);
  proj->SetStats(0);
  proj->GetXaxis()->SetTitle(Form("%s / %s", proj->GetXaxis()->GetTitle(), "#Delta#eta"));
  proj->GetYaxis()->SetRangeUser(proj->GetMinimum() * 1.2, proj->GetMaximum() * 1.5);
  proj->Draw("");

  proj = hist->ProjectionY(Form("%s_proj2_%d", hist->GetName(), 0), hist->GetXaxis()->FindBin(-projectLimit+0.01), hist->GetXaxis()->FindBin(projectLimit-0.01));
  proj->GetXaxis()->SetRangeUser(-1.49, 1.49);
  proj->GetYaxis()->SetRangeUser(proj->GetMinimum() * 1.2, proj->GetMaximum() * 2);
  proj->SetLineColor(2);
  proj->Draw("SAME");
  c->SaveAs("etagap_subtracted_proj.eps");

  c = new TCanvas("c2", "c2", 800, 800);
  gPad->SetLeftMargin(0.15);
  hist->SetTitle("");
  hist->GetYaxis()->SetRangeUser(-1.79, 1.79);
  hist->GetXaxis()->SetTitleOffset(1.5);
  hist->GetYaxis()->SetTitleOffset(2);
  hist->SetStats(kFALSE);
  hist->Draw("SURF1");
  c->SaveAs("etagap_subtracted.eps");
}
