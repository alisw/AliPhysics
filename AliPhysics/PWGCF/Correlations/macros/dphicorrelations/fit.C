#include "TF2.h"
#include "TH2F.h"
#include "TH3.h"
#include "THn.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TFile.h"
#include "TLatex.h"
#include "TROOT.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TPaveText.h"
#include "TRandom.h"
#include "TSystem.h"
#include "TASImage.h"
#include "TObjString.h"
#include "TFitResult.h"
#include "TMatrixTBase.h"
#include "TMultiGraph.h"
#include "TVirtualFitter.h"
#include "Math/DistFunc.h"
#include "TPaletteAxis.h"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include "TStyle.h"
#include "TColor.h"
#include "TSystem.h"
#include "TGaxis.h"

bool gSTARBinning = false;
int energy = 1; // 0: 2.76 TeV, 1: 5.02 TeV, -1: Mixed

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

//            1     2                    3                                  4        5  6  7  8 11                   12     13     14               16   19       20         21         22     23     24     25     26      27
// second fit:yield,dip divided by yield,dip from all bins divided by yield,beta eta,v1,v2,v3,v4,dip divided by norm,phi CP,eta CP,dip for all bins,norm,beta phi,dphi sigma,deta sigma,chi2_1,chi2_2,chi2_3,chi2_4,phi rms,eta rms
const Int_t NGraphs = 51;
const Int_t NHists = 6*4; // pt index
TGraphErrors*** graphs = 0;
TF2*** fitFunctions = 0;
const char* kCorrFuncTitle = "1/N_{trig} dN_{assoc}/d#Delta#etad#Delta#varphi (1/rad.)";
const char* kProjYieldTitlePhi = "1/N_{trig} dN_{assoc}/d#Delta#varphi (1/rad.)";
const char* kProjYieldTitleEta = "1/N_{trig} dN_{assoc}/d#Delta#eta";
const char* kProjYieldTitlePhiOrEta = "1/N_{trig} dN_{assoc}/d#Delta#varphi (1/rad.) , dN_{assoc}/d#Delta#eta";
TString fgFolder = "tmpresults";
const char* fitLabel = "fit";
Int_t vnmax = 4; //Set to 5 if energy = 1 in Analyze and in AnalyzeExample

void CreateGraphStructure()
{
  graphs = new TGraphErrors**[NGraphs];
  for (Int_t i=0; i<NGraphs; i++)
  {
    graphs[i] = new TGraphErrors*[NHists];
    for (Int_t j=0; j<NHists; j++)
      graphs[i][j] = new TGraphErrors;
  }
  fitFunctions = new TF2**[NHists];
  for (Int_t j=0; j<NHists; j++)
  {
    fitFunctions[j] = new TF2*[11];
    for (Int_t i=0; i<11; i++)
      fitFunctions[j][i] = new TF2;
  }
}

void WriteGraphs(const char* outputFileName = "graphs.root")
{
  TFile::Open(outputFileName, "RECREATE");
  for (Int_t i=0; i<NGraphs; i++)
    for (Int_t j=0; j<NHists; j++)
      graphs[i][j]->Write(Form("graph_%d_%d", i, j));
  for (Int_t j=0; j<NHists; j++)
    for (Int_t i=0; i<9; i++)
      fitFunctions[j][i]->Write(Form("fitFunction_%d_%d", j, i));
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

void DivideGraphs(TGraphErrors* graph1, TGraphErrors* graph2, float corr = 0)
{

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
    Float_t error = value * TMath::Sqrt(TMath::Power(graph1->GetEY()[bin1] / graph1->GetY()[bin1], 2) + TMath::Power(graph2->GetEY()[bin2] / graph2->GetY()[bin2], 2) - corr * 2.*graph1->GetEY()[bin1]*graph2->GetEY()[bin1]/graph1->GetY()[bin1]/graph2->GetY()[bin2]);

    graph1->GetY()[bin1] = value;
    graph1->GetEY()[bin1] = error;

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

const Double_t k1OverSqrtTwoPi = 1.0 / TMath::Sqrt(TMath::TwoPi());

Double_t DeltaPhiWidth2DFitFunction(Double_t *x, Double_t *par)
{
  Float_t phi = x[0];
  return par[5] + par[6] * TMath::Cos(2. * phi) + (vnmax > 2)*par[7] * TMath::Cos(3. * phi) + (vnmax > 3)*par[8] * TMath::Cos(4.*phi) + (vnmax > 4)*par[9] * TMath::Cos(5.*phi) + (vnmax > 5)*par[10] * TMath::Cos(6.*phi) 
    + par[0]*( 
	      par[3]*par[4]/4./par[1]/par[2]/TMath::Gamma(1./par[3])/TMath::Gamma(1./par[4]) * 
	      TMath::Exp(-1.*(TMath::Power(TMath::Abs(x[0]/par[1]),par[3])+TMath::Power(TMath::Abs(x[1]/par[2]),par[4]))) 
	     );
}

TH2* SubtractEtaGap(TH2* hist, Float_t etaLimit, Float_t outerLimit, Bool_t scale, Bool_t drawEtaGapDist = kFALSE)
{
  TString histName(hist->GetName());
  Int_t etaBins = 0;

  TH1D* etaGap = hist->ProjectionX(histName + "_1", TMath::Max(1, hist->GetYaxis()->FindBin(-outerLimit + 0.01)), hist->GetYaxis()->FindBin(-etaLimit - 0.01));
//   Printf("%f", etaGap->GetEntries());
  if (etaGap->GetEntries() > 0)
    etaBins += hist->GetYaxis()->FindBin(-etaLimit - 0.01) - TMath::Max(1, hist->GetYaxis()->FindBin(-outerLimit + 0.01)) + 1;

  TH1D* tracksTmp = hist->ProjectionX(histName + "_2", hist->GetYaxis()->FindBin(etaLimit + 0.01), TMath::Min(hist->GetYaxis()->GetNbins(), hist->GetYaxis()->FindBin(outerLimit - 0.01)));
//   Printf("%f", tracksTmp->GetEntries());
  if (tracksTmp->GetEntries() > 0)
    etaBins += TMath::Min(hist->GetYaxis()->GetNbins(), hist->GetYaxis()->FindBin(outerLimit - 0.01)) - hist->GetYaxis()->FindBin(etaLimit + 0.01) + 1;
  
  etaGap->Add(tracksTmp);

  // get per bin result
  if (etaBins > 0)
    etaGap->Scale(1.0 / etaBins);
 
  if (drawEtaGapDist)
  {
    TH1D* centralRegion = hist->ProjectionX(histName + "_3", hist->GetYaxis()->FindBin(-etaLimit + 0.01), hist->GetYaxis()->FindBin(etaLimit - 0.01));
    
//    centralRegion->Scale(1.0 / (hist->GetYaxis()->FindBin(etaLimit - 0.01) - hist->GetYaxis()->FindBin(-etaLimit + 0.01) + 1));
    centralRegion->Scale(hist->GetXaxis()->GetBinWidth(1));

    TCanvas* c = new TCanvas("SubtractEtaGap", "SubtractEtaGap", 800, 800);
    gPad->SetLeftMargin(0.13);
    centralRegion->SetStats(0);
    TString label(centralRegion->GetTitle());
    label.ReplaceAll(".00", " GeV/c");
    label.ReplaceAll(".0", " GeV/c");
    centralRegion->SetTitle(label);
    centralRegion->SetLineColor(3);
    centralRegion->Draw();
    centralRegion->GetYaxis()->SetTitle(kProjYieldTitlePhi);
    centralRegion->GetYaxis()->SetTitleOffset(1.6);
    TH1* copy = etaGap->DrawCopy("SAME");
    copy->Scale(hist->GetXaxis()->GetBinWidth(1));
    copy->Scale((hist->GetYaxis()->FindBin(etaLimit - 0.01) - hist->GetYaxis()->FindBin(-etaLimit + 0.01) + 1));
    copy->SetLineColor(2);
    TLegend* legend = new TLegend(0.41, 0.73, 0.69, 0.85);
    legend->SetFillColor(0);
    legend->AddEntry(centralRegion, Form("|#Delta#eta| < %.1f", etaLimit), "L");
    legend->AddEntry(copy, Form("%.1f < |#Delta#eta| < %.1f (scaled)", etaLimit, outerLimit), "L");
    legend->Draw();
    
    DrawLatex(0.705, 0.62, 1, "Pb-Pb 2.76 TeV", 0.025);
    DrawLatex(0.705, 0.58, 1, "Stat. unc. only", 0.025);
    
    c->SaveAs("note/etagap_proj.eps");
    c->SaveAs("note/etagap_proj.png");
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
  return histTmp2D;
}

void SubtractEtaGap1D(TH1* projPhi, TH1* projPhiSubtractPositive, TH1* projPhiSubtractNegative, TH1* projEta, Float_t etaLimit, Float_t outerLimit, Bool_t draw = kFALSE )
{
  if (draw)
  {
    TCanvas* c1 = new TCanvas("c1", "eta", 800, 600);
    c1->cd();
    projEta->GetYaxis()->SetRangeUser(0,1.1);
    projEta->DrawCopy(); 
  }
  Float_t background = 0;
  Float_t backgroundError = 0;
  Int_t nBins = 0;

  for (Int_t i=projEta->FindBin(-outerLimit + 0.01); i<=projEta->FindBin(-etaLimit - 0.01); i++)
  {
    nBins++;
    background += projEta->GetBinContent(i);
    backgroundError += TMath::Power(projEta->GetBinError(i),2);
  }
  for (Int_t i=projEta->FindBin(etaLimit + 0.01); i<=projEta->FindBin(outerLimit - 0.01); i++)
  { 
    nBins++;
    background += projEta->GetBinContent(i);
    backgroundError += TMath::Power(projEta->GetBinError(i),2);
  }
  background = background/nBins;
  backgroundError = TMath::Sqrt(backgroundError)/nBins;
  for (Int_t i=1; i<=projEta->GetNbinsX(); i++)
  {
    Float_t content = projEta->GetBinContent(i);
    Float_t contentError = projEta->GetBinError(i);
    projEta->SetBinContent(i,content-background);
    projEta->SetBinError(i,TMath::Sqrt(TMath::Power(contentError,2)+TMath::Power(backgroundError,2)));
  }
  if (draw) projEta->DrawCopy("SAME");
  background = 0;
  backgroundError = 0;
  if (draw)
  {
    TCanvas* c2 = new TCanvas("c2", "phi", 800, 600);
    c2->cd();
    projPhi->GetYaxis()->SetRangeUser(0,1.5);
    projPhi->DrawCopy();
    projPhiSubtractPositive->DrawCopy("SAME");
    projPhiSubtractNegative->DrawCopy("SAME");
  }
  projPhiSubtractNegative->Add(projPhiSubtractPositive);
  projPhiSubtractNegative->Scale(1.0/2.0);
  projPhi->Add(projPhiSubtractNegative,-1);
  if (draw) projPhi->DrawCopy("SAME");
}


void AddRMSGeneralized(TGraphErrors* graph, Float_t x, Float_t xE, TF1* func, TMatrixDSym& cov, Int_t sigmaIndex, Int_t betaIndex)
{
  double sigma = func->GetParameter(sigmaIndex);
  double beta = func->GetParameter(betaIndex);
  double rms = TMath::Sqrt(sigma*sigma*TMath::Gamma(3./beta)/TMath::Gamma(1./beta));
  double sigmaDer = TMath::Sqrt(TMath::Gamma(3./beta)/TMath::Gamma(1./beta));
  TF1* tmp = new TF1("tmp","TMath::Sqrt(TMath::Gamma(3./x)/TMath::Gamma(1./x))",1,2);
  double betaDer = sigma*tmp->Derivative(beta);
  double rmsError = 
    TMath::Power(sigmaDer * func->GetParError(sigmaIndex), 2) +
    TMath::Power(betaDer * func->GetParError(betaIndex), 2) +
    2. * sigmaDer * betaDer * cov(sigmaIndex, betaIndex);
//  cerr << "Error: " << TMath::Power(sigmaDer * func->GetParError(sigmaIndex), 2) << "\t" << TMath::Power(betaDer * func->GetParError(betaIndex), 2) << "\t" << 2. * sigmaDer * betaDer * cov(sigmaIndex, betaIndex) << "\t" << beta << endl;
//  cerr << sigmaDer << "\t" << func->GetParError(sigmaIndex) << "\t" << betaDer << "\t" << func->GetParError(betaIndex) << "\t" << cov(sigmaIndex, betaIndex) << endl;
  rmsError = TMath::Sqrt(rmsError);
  AddPoint(graph, x, rms, xE, rmsError);
  cerr << "RMS: " << rms << endl;
}

Bool_t gExcludeRegion = kTRUE; // exclude region around 0, 0 for fit
Int_t gNBin = 0;
Int_t gNBinFit = 0;
double gPhiPeripheralX, gPhiPeripheralY, gPhiPeripheralE, gEtaPeripheralX, gEtaPeripheralY, gEtaPeripheralE;

Bool_t FitDeltaPhi2DOneFunction(TH2* hist, TCanvas* canvas, Int_t graphID, Float_t x, Float_t xE, Int_t histId)
{
  Int_t canvasPos = 1;
  Float_t yPosChi2 = 0.9;
  Bool_t success = kTRUE;

  Float_t etaLimit = 1.4;
  Float_t outerLimit = 1.59;
//  if (energy == 1 ) outerLimit = 1.79;
  // for pt,T < 4 and pt,a < 3 OR 4 < pT,t < 8, 1 < pT,a < 2
  if (graphID <= 11 || graphID == 15)
    etaLimit = 1.39;
  else
    etaLimit = 1.2;
  
  // for pT,a > 4
  if (graphID == 18 || graphID == 23 || graphID == 24)
  {
    etaLimit = 0.5;
    outerLimit = 0.99;
  }
  Float_t sigmaFitLimit = 0.1;
  if (graphID >= 10)
    sigmaFitLimit -= 0.05;

  Float_t etaFitUpperLimit = 0.8;
  Float_t initSigma = 0.6;
  if (histId == 2) // pp
  {
    etaFitUpperLimit = 0.6;
    initSigma = 0.4;
  }
  
  // fit babysitting
  if ((graphID == 10 && histId == 4) || (graphID == 5 && histId == 0) || (graphID == 5 && histId == 3))
    etaFitUpperLimit += 0.1;
  
  if (graphID == 0)
    etaFitUpperLimit += 0.2;
  cerr << "graphID = " << graphID << endl;
  Printf("Limits (%d %d): etaLimit = %f outerLimit = %f sigmaFitLimit = %f etaFitUpperLimit = %f initSigma = %f", graphID, histId, etaLimit, outerLimit, sigmaFitLimit, etaFitUpperLimit, initSigma);

  TH2* dipHist = (TH2*) hist->Clone("dipHist");
  TH2* peakHist = (TH2*) hist->Clone("peakHist");
  // set errors large for bins around (0,0)
  Float_t exclusionRegion = 0.019;
  if (!gSTARBinning && gExcludeRegion && graphID <= 10 && histId != 1 && histId != 2)
  {
    exclusionRegion = 0.15;
    if (graphID == 0)
      exclusionRegion = 0.29;
    if (graphID == 10 || graphID == 6)
      exclusionRegion = 0.075;
    if (graphID == 0 && histId >= 7)
      exclusionRegion = 0.039;
    if (graphID != 0 && histId >= 7) 
      exclusionRegion = 0.019;
    for (Int_t binx = hist->GetXaxis()->FindBin(-exclusionRegion)-gNBinFit; binx <= hist->GetXaxis()->FindBin(exclusionRegion)+gNBinFit; binx++)
      for (Int_t biny = hist->GetYaxis()->FindBin(-exclusionRegion)-gNBinFit; biny <= hist->GetYaxis()->FindBin(exclusionRegion)+gNBinFit; biny++)
      {
// 	hist->SetBinContent(binx, biny,  0);
	hist->SetBinError(binx, biny,  1e5);
      }
  }
  if (exclusionRegion != 0)  Printf("NOTE : Skipping bins at (0, 0) with exclusion region %0.2f",exclusionRegion);
  if (gNBinFit > 0 && exclusionRegion < 0.02)
  {
    for (Int_t binx = hist->GetXaxis()->FindBin(-0.04)-(gNBinFit-1); binx <= hist->GetXaxis()->FindBin(0.04)+(gNBinFit-1); binx++)
      for (Int_t biny = hist->GetYaxis()->FindBin(-0.04)-(gNBinFit-1); biny <= hist->GetYaxis()->FindBin(0.04)+(gNBinFit-1); biny++)
      {
// 	hist->SetBinContent(binx, biny,  0);
	hist->SetBinError(binx, biny,  1e5);
      }
  }
  

  hist->GetYaxis()->SetRangeUser(-outerLimit+0.01, outerLimit-0.01);
  hist->GetXaxis()->SetRangeUser(-TMath::Pi() / 2 + 0.01, TMath::Pi() * 0.5 - 0.01);
  
  Float_t mean = hist->Integral(hist->GetXaxis()->FindBin(-TMath::Pi() / 2), hist->GetXaxis()->FindBin(TMath::Pi() / 2), hist->GetYaxis()->FindBin(etaLimit), hist->GetYaxis()->FindBin(outerLimit)) / (hist->GetXaxis()->FindBin(TMath::Pi() / 2) - hist->GetXaxis()->FindBin(-TMath::Pi() / 2)) / (hist->GetYaxis()->FindBin(outerLimit) - hist->GetYaxis()->FindBin(etaLimit) + 1);
//   Printf("%f", mean);

  canvas->cd(canvasPos++);
  hist->SetStats(0);
  hist->DrawCopy("SURF1");
  
  Float_t min = hist->GetMinimum();
  Float_t max = hist->GetMaximum();

  Int_t bins = hist->GetNbinsX() / 2 / 2;
   
  float norm = hist->GetBinContent(hist->GetXaxis()->FindBin(0.0), hist->GetYaxis()->FindBin(0.0)) - mean;
  //cerr << "norm " << norm  << ", bin content: " << hist->GetBinContent(hist->GetXaxis()->FindBin(0.0)) << endl;
  if (norm < 0) norm = hist->GetBinContent(hist->GetMaximumBin());
  TF2* func = new TF2("func", DeltaPhiWidth2DFitFunction, -TMath::Pi() / 2, TMath::Pi() * 0.5, -outerLimit, outerLimit, 5+vnmax);
  func->SetParameters(norm, initSigma, initSigma, 1.9, 1.9);
  double betaUpperLimit = 10.;
  double betaLowerLimit = 0.1;
  func->SetParameter(5, mean);
  func->SetParameter(6, 1.e-2);
  func->SetParameter(7, 1.e-2);
  if (vnmax > 3) func->SetParameter(8, 1.e-2);
  if (vnmax > 4) func->SetParameter(9, 1.e-2);
  if (vnmax > 5) func->SetParameter(10, 1.e-2);
  func->SetParLimits(5,mean*0.8,100);
  if (1)
  {
    TVirtualFitter::SetMaxIterations(15000); 
    printf("STEP 1: fit only flow using one delta eta side");
    for (Int_t i=0; i<5; i++)
      func->FixParameter(i, func->GetParameter(i));
    hist->GetYaxis()->SetRangeUser(etaLimit+0.01, outerLimit-0.01);
    for (Int_t i=6; i<5+vnmax; i++)
      func->SetParLimits(i, 1.e-15, 100.);
    Int_t fitResult = hist->Fit(func, "0", "");
    Printf("Fit result: %d; Chi2/ndf: %f/%d", fitResult, func->GetChisquare(), func->GetNDF());
    bool isFlow = true;
    if (func->GetParameter(6) < 1e-6 && func->GetParameter(7) < 1e-6 && func->GetParameter(8) < 1e-6)
     isFlow = false;
//    if (fitResult != 0)
//      success = kFALSE;

    printf("STEP2 : fit only Gaussian in central region");
    for (Int_t i=0; i<5; i++)
      func->ReleaseParameter(i);
    for (Int_t i=5; i<5+vnmax; i++)
      func->FixParameter(i, func->GetParameter(i));
    func->SetParLimits(1, 0.15, 2.);
    func->SetParLimits(2, 0.15, 2.);
    func->SetParLimits(3, betaLowerLimit, betaUpperLimit);
    func->SetParLimits(4, betaLowerLimit, betaUpperLimit);
    hist->GetYaxis()->SetRangeUser(-etaLimit+0.01, etaLimit-0.01);
    fitResult = hist->Fit(func, "0", "");
    Printf("Fit result: %d; Chi2/ndf: %f/%d", fitResult, func->GetChisquare(), func->GetNDF());
//    if (fitResult != 0)
//      success = kFALSE;
//    for (Int_t i=0; i<6+vnmax; i++) func->SetParLimits(i, 0, 0);

    printf("STEP3: fit everything, with limits");
    for (Int_t i=0; i<5+vnmax; i++)
    {
      func->ReleaseParameter(i);
      if (func->GetParameter(i) > 0)
	func->SetParLimits(i, func->GetParameter(i) * 0.8, func->GetParameter(i) * 1.2);
      else
	func->SetParLimits(i, func->GetParameter(i) * 1.2, func->GetParameter(i) * 0.8);
    }
    if (func->GetParameter(3) * 0.8 < betaLowerLimit && func->GetParameter(3) * 1.2 > betaUpperLimit) func->SetParLimits(3,betaLowerLimit,betaUpperLimit);
    else if (func->GetParameter(3) * 0.8 < betaLowerLimit) func->SetParLimits(3,betaLowerLimit,func->GetParameter(3) * 1.2);
    else if (func->GetParameter(3) * 1.2 > betaUpperLimit) func->SetParLimits(3,func->GetParameter(3) * 0.8,betaUpperLimit);
    if (func->GetParameter(4) * 0.8 < betaLowerLimit && func->GetParameter(4) * 1.2 > betaUpperLimit) func->SetParLimits(4,betaLowerLimit,betaUpperLimit);
    else if (func->GetParameter(4) * 0.8 < betaLowerLimit) func->SetParLimits(4,betaLowerLimit,func->GetParameter(4) * 1.2);
    else if (func->GetParameter(4) * 1.2 > betaUpperLimit) func->SetParLimits(4,func->GetParameter(4) * 0.8,betaUpperLimit);
    if (!isFlow)
//    if (func->GetParameter(5) < 1e-3)
      for (Int_t i=6; i<5+vnmax; i++)
      {
        func->SetParLimits(i, 0, 0);
        func->FixParameter(i, 0);
      }
    if ((func->GetParameter(6) < 1e-6 && func->GetParameter(7) < 1e-6) || (func->GetParameter(6) < 1e-6 && func->GetParameter(8) < 1e-6) || (func->GetParameter(7) < 1e-6 && func->GetParameter(8) < 1e-6))
      for (Int_t i=6; i<5+vnmax; i++)
        func->SetParLimits(i,0,0);

    hist->GetYaxis()->SetRangeUser(-outerLimit+0.01, outerLimit-0.01);
    fitResult = hist->Fit(func, "0", "");
    Printf("Fit result: %d; Chi2/ndf: %f/%d", fitResult, func->GetChisquare(), func->GetNDF());
    if (func->GetParameter(6) < 1e-6 && func->GetParameter(7) < 1e-6 && func->GetParameter(8) < 1e-6)
     isFlow = false;
//    if (fitResult != 0)
//      success = kFALSE;

    printf("STEP4: fit everything, without limits");
    for (Int_t i=0; i<5+vnmax; i++) func->SetParLimits(i, 0, 0);
    func->SetParLimits(3, betaLowerLimit, betaUpperLimit);
    func->SetParLimits(4, betaLowerLimit, betaUpperLimit);

    for (Int_t i=6; i<5+vnmax; i++)
    {
      if (func->GetParameter(i) < 0) func->SetParameter(i,0);
      func->SetParLimits(i,0,100);
    }
    if (func->GetParameter(6) < 1e-6 && func->GetParameter(7) < 1e-6 && func->GetParameter(8) < 1e-6 && (vnmax > 4 && func->GetParameter(9) < 1e-6) && (vnmax > 5 && func->GetParameter(10) < 1e-6))
      isFlow = false;
    isFlow = true;
    if (!isFlow)
      for (Int_t i=6; i<5+vnmax; i++)
      {
        func->SetParLimits(i, 0, 0);
        func->FixParameter(i, 0);
      }
//    for (Int_t i=7; i<5+vnmax; i++)
//     func->SetParLimits(i, 0., 2.);
  } // if (1)
  TFitResultPtr fitResultPtr = hist->Fit(func, "S0", "");
  int fitResult = fitResultPtr;
  if (fitResult != 0)
  {
    printf("STEP5: fit everything, without limits again, but with no flow\n");
    for (Int_t i=0; i<5+vnmax; i++) func->SetParLimits(i, 0, 0);
    func->SetParLimits(3, betaLowerLimit, betaUpperLimit);
    func->SetParLimits(4, betaLowerLimit, betaUpperLimit);
 //   if (func->GetParameter(6) < 1e-6 || func->GetParameter(7) < 1e-6 || func->GetParameter(8) < 1e-6)
    for (Int_t i=6; i<5+vnmax; i++)
      if (func->GetParameter(i) < 1e-8)
        func->FixParameter(i, 0);
    fitResultPtr = hist->Fit(func, "S0", "");
    fitResult = fitResultPtr;
  }

  if (fitResult != 0)
  {
    success = kFALSE;
    Printf("Finished with %d", success);
    return success;
  }
  else success = kTRUE;
  TMatrixDSym cov = fitResultPtr->GetCovarianceMatrix();
  cov.Print();
  
  func->SetName(Form("fitFunction_%d_%d", graphID, histId));
  func->Copy(*fitFunctions[graphID][histId]);
  //dphi
  AddRMSGeneralized(graphs[10+16][graphID], x, xE, func, cov,1,3);
  AddPoint(graphs[4+16][graphID], x, TMath::Abs(func->GetParameter(1)), xE, func->GetParError(1));
  //deta
  AddRMSGeneralized(graphs[11+16][graphID], x, xE, func, cov,2,4);
  AddPoint(graphs[5+16][graphID], x, TMath::Abs(func->GetParameter(2)), xE, func->GetParError(2));
  if (func->GetParameter(6) >= 0)
    AddPoint(graphs[5][graphID], x, TMath::Sqrt(func->GetParameter(6)/2./func->GetParameter(5)), xE, func->GetParError(6)); //v2
  else 
    AddPoint(graphs[5][graphID], x, -1.*TMath::Sqrt(-1.*func->GetParameter(6)/2./func->GetParameter(5)), xE, func->GetParError(6)); 
  if (func->GetParameter(7) >= 0)
    AddPoint(graphs[6][graphID], x, TMath::Sqrt(func->GetParameter(7)/2./func->GetParameter(5)), xE, func->GetParError(7)); //v3
  else
    AddPoint(graphs[6][graphID], x, -1.*TMath::Sqrt(-1.*func->GetParameter(7)/2./func->GetParameter(5)), xE, func->GetParError(7));
  if (func->GetParameter(8) >= 0)
    AddPoint(graphs[7][graphID], x, TMath::Sqrt(func->GetParameter(8)/2./func->GetParameter(5)), xE, func->GetParError(8)); //v4
  else 
    AddPoint(graphs[7][graphID], x, -1.*TMath::Sqrt(-1.*func->GetParameter(8)/2./func->GetParameter(5)), xE, func->GetParError(8));
  if (func->GetParameter(9) >= 0)
    AddPoint(graphs[8][graphID], x, TMath::Sqrt(func->GetParameter(9)/2./func->GetParameter(5)), xE, func->GetParError(9)); //v4
  else 
    AddPoint(graphs[8][graphID], x, -1.*TMath::Sqrt(-1.*func->GetParameter(9)/2./func->GetParameter(5)), xE, func->GetParError(9));
  AddPoint(graphs[0][graphID], x, func->GetParameter(0), xE, func->GetParError(0));
  AddPoint(graphs[28][graphID], x, func->GetParameter(6)/2., xE, func->GetParError(6)/2.); //v2
  AddPoint(graphs[29][graphID], x, func->GetParameter(7)/2., xE, func->GetParError(7)/2.); //v3
  AddPoint(graphs[30][graphID], x, func->GetParameter(8)/2., xE, func->GetParError(8)/2.); //v4
  AddPoint(graphs[32][graphID], x, func->GetParameter(9)/2., xE, func->GetParError(9)/2.); //v5
  AddPoint(graphs[31][graphID], x, func->GetParameter(5), xE, func->GetParError(5)); //norm

//  cerr << "x: " << x << endl;
  if (x == 65)
  {
   graphs[10+16][graphID]->GetPoint(graphs[10+16][graphID]->GetN()-1,gPhiPeripheralX,gPhiPeripheralY);
   gPhiPeripheralE = graphs[10+16][graphID]->GetErrorY(graphs[10+16][graphID]->GetN()-1);
   graphs[11+16][graphID]->GetPoint(graphs[11+16][graphID]->GetN()-1,gEtaPeripheralX,gEtaPeripheralY);
   gEtaPeripheralE = graphs[11+16][graphID]->GetErrorY(graphs[11+16][graphID]->GetN()-1);
//   cerr << gPhiPeripheralX << "\t" << gPhiPeripheralY << "\t" << gPhiPeripheralE  << "\t" << gEtaPeripheralX << "\t" << gEtaPeripheralY << "\t" << gEtaPeripheralE << endl;
  }
  else if (x == 5)
  {
    double phiCentralX, phiCentralY, phiCentralE, etaCentralX, etaCentralY, etaCentralE;;
    graphs[10+16][graphID]->GetPoint(graphs[10+16][graphID]->GetN()-1,phiCentralX,phiCentralY);
    phiCentralE = graphs[10+16][graphID]->GetErrorY(graphs[10+16][graphID]->GetN()-1);
    AddPoint(graphs[12][graphID], x, phiCentralY/gPhiPeripheralY, 0, TMath::Sqrt(phiCentralY*phiCentralY/gPhiPeripheralY/gPhiPeripheralY*(phiCentralE*phiCentralE/phiCentralY/phiCentralY+gPhiPeripheralE*gPhiPeripheralE/gPhiPeripheralY/gPhiPeripheralY)));
    graphs[11+16][graphID]->GetPoint(graphs[11+16][graphID]->GetN()-1,etaCentralX,etaCentralY);
    etaCentralE = graphs[11+16][graphID]->GetErrorY(graphs[11+16][graphID]->GetN()-1);
    AddPoint(graphs[13][graphID], x, etaCentralY/gEtaPeripheralY, 0, TMath::Sqrt(etaCentralY*etaCentralY/gEtaPeripheralY/gEtaPeripheralY*(etaCentralE*etaCentralE/etaCentralY/etaCentralY+gEtaPeripheralE*gEtaPeripheralE/gEtaPeripheralY/gEtaPeripheralY)));
  }
  TF2* BG = new TF2("BG", "[0] + [1] * TMath::Cos(2. * x) + [2] * TMath::Cos(3. * x) + [3] * TMath::Cos(4.*x) + 0*y", -TMath::Pi() / 2, TMath::Pi() /2, -outerLimit, outerLimit);
  for (int i=0; i<=3; i++)
    BG->SetParameter(i,func->GetParameter(i+5));
  int scale = 10;
  TH2* funcHighResolution = new TH2F("funcHighResolution","funcHighResolution",hist->GetXaxis()->GetNbins()*scale,-TMath::Pi() / 2,3*TMath::Pi() / 2,hist->GetYaxis()->GetNbins()*scale,hist->GetYaxis()->GetXmin(),hist->GetYaxis()->GetXmax());
  funcHighResolution->Eval(func,"A");
  funcHighResolution->Rebin2D(scale,scale);
  funcHighResolution->Scale(1./scale/scale);
  peakHist->GetYaxis()->SetRangeUser(-outerLimit+0.01, outerLimit-0.01);
  peakHist->GetXaxis()->SetRangeUser(-TMath::Pi() / 2 + 0.01, TMath::Pi() * 0.5 - 0.01);
  dipHist->GetYaxis()->SetRangeUser(-outerLimit+0.01, outerLimit-0.01);
  dipHist->GetXaxis()->SetRangeUser(-TMath::Pi() / 2 + 0.01, TMath::Pi() * 0.5 - 0.01);
  peakHist->Sumw2();
  peakHist->Add(BG, -1);
  dipHist->Sumw2();
  funcHighResolution->Sumw2();
//  dipHist->Add(funcHighResolution, -1);
  dipHist->Add(func, -1);
  bool binsExcluded = true;
  if (exclusionRegion < 0.02 ) binsExcluded = false;
  int nBin = gNBin; 
  if (exclusionRegion < 0.02 || hist->FindBin(-exclusionRegion)-gNBin > hist->FindBin(-exclusionRegion)+gNBin)
  {
    exclusionRegion = 0.08;
    if (nBin > 0) nBin--;
    else if (nBin < 0) nBin = 0;
  }
//  cerr << "exclusionRegion: " << exclusionRegion << endl;
  double yieldIntegral = peakHist->Integral("width");
  AddPoint(graphs[1][graphID], x, yieldIntegral, xE, 0);
  double integral = 0, integralError = 0;
  double integralAll = 0, integralErrorAll = 0;
  double centralExclusionForDip = 0.04; //0.08 for normal binning
  if (exclusionRegion > centralExclusionForDip)
  {
    TH2* dipHistClone = (TH2*) dipHist->Clone("dipHistClone");
    for (int xBin=dipHistClone->GetXaxis()->FindBin(-centralExclusionForDip); xBin<=dipHistClone->GetXaxis()->FindBin(centralExclusionForDip); xBin++)
      for (int y=dipHistClone->GetYaxis()->FindBin(-centralExclusionForDip); y<=dipHistClone->GetYaxis()->FindBin(centralExclusionForDip); y++) 
      {
	dipHistClone->SetBinContent(xBin, y, 0);
	dipHistClone->SetBinError(xBin, y, 0);
      }
      
    integral = dipHistClone->IntegralAndError(dipHistClone->GetXaxis()->FindBin(-exclusionRegion)-nBin,dipHistClone->GetXaxis()->FindBin(exclusionRegion)+nBin,dipHist->GetYaxis()->FindBin(-exclusionRegion)-nBin,dipHist->GetYaxis()->FindBin(exclusionRegion)+nBin,integralError,"width");
  
    if (binsExcluded)
    {    
      AddPoint(graphs[2][graphID], x, -1.*integral/yieldIntegral, xE, integralError/yieldIntegral);
      AddPoint(graphs[11][graphID], x, -1.*integral/func->GetParameter(0), xE, TMath::Sqrt(integral*integral/func->GetParameter(0)/func->GetParameter(0)*(integralError*integralError/integral/integral+func->GetParError(0)*func->GetParError(0)/func->GetParameter(0)/func->GetParameter(0))));
      Printf("Dip yield: %f, (dip/yield=%f/%f)", -1.*integral/func->GetParameter(0),-1.*integral,func->GetParameter(0));
    }
    else if (integral != 0) AddPoint(graphs[14][graphID], x, -1.*integral/func->GetParameter(0), xE, TMath::Sqrt(integral*integral/func->GetParameter(0)/func->GetParameter(0)*(integralError*integralError/integral/integral+func->GetParError(0)*func->GetParError(0)/func->GetParameter(0)/func->GetParameter(0))));
    delete dipHistClone;
  }
  integralAll = dipHist->IntegralAndError(dipHist->GetXaxis()->FindBin(-TMath::Pi() / 2 + 0.01),dipHist->GetXaxis()->FindBin(TMath::Pi() * 0.5 - 0.01),dipHist->GetYaxis()->FindBin(-outerLimit+0.01),dipHist->GetYaxis()->FindBin(outerLimit-0.01),integralErrorAll,"width");
  AddPoint(graphs[3][graphID], x, -1.*integralAll/yieldIntegral, xE, integralErrorAll/yieldIntegral);
  // norm
  // average of baseline
  Float_t avg = 0;
  for (Int_t i=6; i<bins+6; i++)
    avg += func->GetParameter(i);
  avg /= bins;  
  AddPoint(graphs[0+16][graphID], x, avg, xE, 0);
  // beta
  AddPoint(graphs[3+16][graphID], x, func->GetParameter(3), xE, func->GetParError(3));
  AddPoint(graphs[4][graphID], x, func->GetParameter(4), xE, func->GetParError(4));
  
  canvas->cd(canvasPos++);
  TH2* funcHist = (TH2*) hist->Clone("funcHist");
  funcHist->Reset();
  funcHist->Add(func);
  funcHist->SetMinimum(min);
  funcHist->SetMaximum(max);
  funcHist->Draw("SURF1");
  
  canvas->cd(canvasPos++);
  TH2* residuals = (TH2*) hist->Clone("residuals");
  residuals->Add(func, -1);
  residuals->Draw("COLZ");

  canvas->cd(canvasPos++);
  TH2* residuals_norm = (TH2*) hist->Clone("residuals_norm");
  residuals_norm->Reset();
  for (Int_t i=1; i<=residuals_norm->GetNbinsX(); i++)
    for (Int_t j=1; j<=residuals_norm->GetNbinsY(); j++)
      if (residuals->GetBinError(i, j) > 0)
	residuals_norm->SetBinContent(i,j,residuals->GetBinContent(i, j) / residuals->GetBinError(i, j));
  residuals_norm->Draw("COLZ");

  Double_t chi2 = 0, chi2WithoutMiddle = 0, chi2FullRegionWithoutMiddle = 0;
  Int_t ndf = 0, ndfWithoutMiddle = 0, ndfFullRegionWithoutMiddle = 0;
  for (Int_t i=hist->GetXaxis()->FindBin(-0.8); i<=hist->GetXaxis()->FindBin(0.8); i++)
    for (Int_t j=hist->GetYaxis()->FindBin(-etaLimit); j<=hist->GetYaxis()->FindBin(etaLimit); j++)
//  for (Int_t i=hist->GetXaxis()->FindBin(-1.2); i<=hist->GetXaxis()->FindBin(1.2); i++)
//    for (Int_t j=hist->GetYaxis()->FindBin(-etaFitUpperLimit); j<=hist->GetYaxis()->FindBin(etaFitUpperLimit); j++)
    {
      if (residuals->GetBinError(i, j) > 0)
      {
	chi2 += TMath::Power(residuals->GetBinContent(i, j) / residuals->GetBinError(i, j), 2);
	ndf++;
      }
    }
  ndf -= func->GetNumberFreeParameters();
  for (Int_t i=hist->GetXaxis()->FindBin(-0.8); i<=hist->GetXaxis()->FindBin(0.8); i++)
    for (Int_t j=hist->GetYaxis()->FindBin(-etaLimit); j<=hist->GetYaxis()->FindBin(etaLimit); j++)
//  for (Int_t i=hist->GetXaxis()->FindBin(-1.2); i<=hist->GetXaxis()->FindBin(1.2); i++)
//    for (Int_t j=hist->GetYaxis()->FindBin(-etaFitUpperLimit); j<=hist->GetYaxis()->FindBin(etaFitUpperLimit); j++)
    {
      if (hist->GetXaxis()->FindBin(-exclusionRegion)-nBin < i && hist->GetXaxis()->FindBin(exclusionRegion)+nBin > i && hist->GetYaxis()->FindBin(-exclusionRegion)-nBin < j && hist->GetYaxis()->FindBin(exclusionRegion)+nBin > j) continue;
      if (residuals->GetBinError(i, j) > 0)
      {
        chi2WithoutMiddle += TMath::Power(residuals->GetBinContent(i, j) / residuals->GetBinError(i, j), 2);
        ndfWithoutMiddle++;
      }
    }
  ndfWithoutMiddle -= func->GetNumberFreeParameters();

  for (Int_t i=1; i<=hist->GetNbinsX(); i++)
    for (Int_t j=1; j<=hist->GetNbinsY(); j++)
    {
      if (hist->GetXaxis()->FindBin(-exclusionRegion)-nBin < i && hist->GetXaxis()->FindBin(exclusionRegion)+nBin > i && hist->GetYaxis()->FindBin(-exclusionRegion)-nBin < j && hist->GetYaxis()->FindBin(exclusionRegion)+nBin > j) 
      {
        if (residuals->GetBinError(i, j) > 0)
        {
          chi2FullRegionWithoutMiddle += TMath::Power(residuals->GetBinContent(i, j) / residuals->GetBinError(i, j), 2);
          ndfFullRegionWithoutMiddle++;
        }
      }
    }
  ndfFullRegionWithoutMiddle -= func->GetNumberFreeParameters();
  if (func->GetNDF() > 0)

  {
    printf("#chi^{2}/ndf = %.1f/%d = %.1f  ", func->GetChisquare(), func->GetNDF(), func->GetChisquare() / func->GetNDF());
    DrawLatex(0.2, yPosChi2, 1, Form("#chi^{2}/ndf = %.1f/%d = %.1f", func->GetChisquare(), func->GetNDF(), func->GetChisquare() / func->GetNDF()));
    AddPoint(graphs[6+16][graphID], x, func->GetChisquare() / func->GetNDF(), xE, 0);
  }
  if (ndf)
  {
    Printf("#chi^{2}/ndf = %.1f/%d = %.1f", chi2, ndf, chi2 / ndf);
    DrawLatex(0.2, yPosChi2 - 0.05, 1, Form("#chi^{2}/ndf = %.1f/%d = %.1f", chi2, ndf, chi2 / ndf));
    AddPoint(graphs[7+16][graphID], x, chi2 / ndf, xE, 0);
  }
  if (ndfWithoutMiddle)
  {
    Printf("#chi^{2}/ndf = %.1f/%d = %.1f", chi2WithoutMiddle, ndfWithoutMiddle, chi2WithoutMiddle / ndfWithoutMiddle);
//    DrawLatex(0.2, yPosChi2 - 0.05, 1, Form("#chi^{2}/ndf = %.1f/%d = %.1f", chi2WithoutMiddle, ndfWithoutMiddle, chi2WithoutMiddle / ndfWithoutMiddle));
    AddPoint(graphs[8+16][graphID], x, chi2WithoutMiddle / ndfWithoutMiddle, xE, 0);
  }
  if (ndfFullRegionWithoutMiddle)
  {
    Printf("#chi^{2}/ndf = %.1f/%d = %.1f", func->GetChisquare()-chi2FullRegionWithoutMiddle, func->GetNDF()-ndfFullRegionWithoutMiddle, (func->GetChisquare()-chi2FullRegionWithoutMiddle) / (func->GetNDF()-ndfFullRegionWithoutMiddle));
//    DrawLatex(0.2, yPosChi2 - 0.05, 1, Form("#chi^{2}/ndf = %.1f/%d = %.1f", func->GetChisquare()-chi2FullRegionWithoutMiddle, func->GetNDF()-ndfFullRegionWithoutMiddle, (func->GetChisquare()-chi2FullRegionWithoutMiddle) / (func->GetNDF()-ndfFullRegionWithoutMiddle)));
    AddPoint(graphs[9+16][graphID], x, (func->GetChisquare()-chi2FullRegionWithoutMiddle) / (func->GetNDF()-ndfFullRegionWithoutMiddle), xE, 0);
  }
  
  // draw gaussian only
  TF2* funcClone = new TF2("funcClone", DeltaPhiWidth2DFitFunction, -TMath::Pi() / 2, TMath::Pi() * 0.5, -outerLimit, outerLimit, bins+6);
  for (Int_t i=0; i<6; i++)
    funcClone->SetParameter(i, func->GetParameter(i));
  for (Int_t i=6; i<bins+6; i++)
    funcClone->SetParameter(i, 0);
//   funcClone->Print();
/*  canvas->cd(canvasPos++);
  funcHist = (TH2*) hist->Clone("funcHistb");
  funcHist->Reset();
  funcHist->Add(funcClone);
  funcHist->SetMinimum(0);
  funcHist->SetMaximum(max - min);
  funcHist->Draw("SURF1");
  
  canvas->cd(canvasPos++);
  func->SetParameter(0, 0);
  TH2* subtractFlow = (TH2*) hist->Clone("subtractFlow");
  subtractFlow->Add(func, -1);
  subtractFlow->SetMinimum(0);
  subtractFlow->SetMaximum(max - min);
  subtractFlow->DrawCopy("SURF1");

  // eta gap subtraction
  SubtractEtaGap(hist, etaLimit, outerLimit, kFALSE);

  canvas->cd(canvasPos++);
  hist->SetMinimum(-(max - min) / 2);
  hist->SetMaximum((max - min) / 2);
  hist->DrawCopy("SURF1");

  canvas->cd(canvasPos++);
  TH2* difference = (TH2*) hist->Clone("difference");
  difference->Add(subtractFlow, -1);
  difference->SetMinimum(-(max - min) / 2);
  difference->SetMaximum((max - min) / 2);
  difference->DrawCopy("SURF1");
*/
  Printf("Finished with %d", success);

  return success;
}

TH2* RebinTTR(TH2* src, bool ratio)
{
  TH2* target = new TH2F(Form("%s_rebin", src->GetName()), src->GetTitle(), 72, -TMath::Pi()/2, TMath::Pi()*3/2, 40, -2, 2);

  int xBegin = 0.25 * target->GetNbinsX();
  int xEnd   = 0.25 * target->GetNbinsX() + 1;
  int yBegin = target->GetNbinsY() / 2;
  int yEnd   = target->GetNbinsY() / 2 + 1;

  cerr << xBegin << "\t" << xEnd << "\t" << yBegin << "\t" << yEnd << endl;

  for (int i=1; i<=target->GetNbinsX(); i++)
    for (int j=1; j<=target->GetNbinsY(); j++)
    {
      if ((i >= xBegin && i <= xEnd) || (j >= yBegin && j <= yEnd))
        continue;

      int src_i = src->GetXaxis()->FindBin(target->GetXaxis()->GetBinCenter(i));
      int src_j = src->GetYaxis()->FindBin(target->GetYaxis()->GetBinCenter(j));

      target->SetBinContent(i, j, src->GetBinContent(src_i, src_j));
      target->SetBinError(i, j, src->GetBinError(src_i, src_j));
    }

  for (int i=1 ; i<=target->GetNbinsX(); i++)
    for (int j=1; j<=target->GetNbinsY(); j++)
    {
      if ((i<xBegin || i>xEnd) && (j<yBegin || j>yEnd)) continue;
      int lowEdgeX  = src->GetXaxis()->FindBin(target->GetXaxis()->GetBinLowEdge(i) + 1e-4);
      int highEdgeX = src->GetXaxis()->FindBin(target->GetXaxis()->GetBinUpEdge(i) - 1e-4);
      int lowEdgeY  = src->GetYaxis()->FindBin(target->GetYaxis()->GetBinLowEdge(j) + 1e-4);
      int highEdgeY = src->GetYaxis()->FindBin(target->GetYaxis()->GetBinUpEdge(j) - 1e-4);
      double error;
      double value = src->IntegralAndError(lowEdgeX, highEdgeX, lowEdgeY, highEdgeY, error, ratio?"width":"");
      if (ratio)
      {
        value /= target->GetXaxis()->GetBinWidth(i) * target->GetYaxis()->GetBinWidth(j);
        error /= target->GetXaxis()->GetBinWidth(i) * target->GetYaxis()->GetBinWidth(j);
      }
      target->SetBinContent(i, j, value);
      target->SetBinError(i, j, error);
    }
  return target;
}

Int_t Analyze(const char* fileName, const char* outputFileName = "graphs.root")
{
  if (energy == 1)
    vnmax = 5;
  Int_t binsFailed = 0;
  gROOT->SetBatch(kTRUE);
  
  CreateGraphStructure();
  TFile* histFile = new TFile(fileName,"READONLY");
  if (!histFile || histFile->IsZombie())
    return -1;
  TList checkList;
  
  Int_t maxLeadingPt = 4;
  Int_t maxAssocPt = 6;

  for (Int_t i=0; i<maxLeadingPt; i++)
  {
//    if (i < 2) continue;
    int jMin = 1;
    if (gSTARBinning) jMin = 0;
    for (Int_t j=jMin; j<maxAssocPt; j++)
    {
//      if (j < 2) continue;
      // only process when first is filled
      TH2* hist1 = (TH2*) histFile->Get(Form("dphi_%d_%d_%d", i, j+jMin, 0));
      if (!hist1)
      {
        cout << "Hist1 does not exist." << endl;
	continue;
      }

      if (hist1->GetEntries() < 10)
      {
	Printf("%d %d Only %f entries. Skipping...", i, j, hist1->GetEntries());
	continue;
      }
     
      vector<int> orderConvert(11);
      orderConvert[0] = 2;
      orderConvert[1] = 1;
      orderConvert[2] = 5;
      orderConvert[3] = 3;
      orderConvert[4] = 4;
      orderConvert[5] = 0;
      orderConvert[6] = 6;
      orderConvert[7] = 7;
      orderConvert[8] = 8;
      orderConvert[9] = 9;
      orderConvert[10] = 10;
      for (Int_t histId = 0; histId < 11; histId++)
      {
	hist1 = (TH2*) histFile->Get(Form("dphi_%d_%d_%d", i, j+jMin, orderConvert[histId]));
	if (!hist1)
        {
          cout << "Hist1 does not exist." << endl;
          continue;
        }
        hist1->Scale(1.0 / hist1->GetYaxis()->GetBinWidth(1));

        for (int xBin=1; xBin<hist1->GetXaxis()->GetNbins()+1; xBin++)
	
	if (hist1->GetEntries() < 10)
	{
	  Printf("%d %d %d Only %f entries. Skipping...", i, j, orderConvert[histId], hist1->GetEntries());
	  continue;
	}
//        if (i == 4 || (i > 1 && histId > 5)) hist1 = RebinTTR(hist1,true);
	TCanvas* canvas = new TCanvas(Form("DeltaPhi_%d_%d_%d_%d", orderConvert[histId], i, j, 1), Form("DeltaPhi_%d_%d_%d_%d", orderConvert[histId], i, j, 1), 1400, 1100);
	canvas->Divide(5, 3);
	
	for (Int_t k=1; k<=3; k++)
	{
	  canvas->cd(3 * j + k);
	  gPad->SetLeftMargin(0.15);
	  gPad->SetBottomMargin(0.2);
	  gPad->SetTopMargin(0.01);
	  gPad->SetRightMargin(0.01);
	}
	
	Int_t graphID = i * (maxAssocPt - 1) + j - jMin;

	Printf("\n\n>>> %d %d %d %d", i, j, orderConvert[histId], graphID);
	
	if (orderConvert[histId] == 0)
	  for (Int_t k=0; k<NGraphs; k++)
	    graphs[k][graphID]->SetTitle(hist1->GetTitle());
	
//        Float_t centralityAxisMapping[] =  { 5, 15, 25, 35, 45, 55, 65, 75, 85 };
//        Float_t centralityAxisMappingE[] = { 5,  5,   5,  5,  5,  5,  5,  5,  5 };
        Float_t centralityAxisMapping[] =  { 5, 95, 100, 25, 15, 35, 45, 55, 65, 75, 85 };
        Float_t centralityAxisMappingE[] = { 5,  5, 0,    5,  5,  5,  5,  5,  5,  5, 5 };
//        Float_t centralityAxisMapping[] =  { 5, 200, 100, 25, 15, 35, 45, 55, 65, 75 };
//        Float_t centralityAxisMappingE[] = { 5,  5, 0,    5,  5,  5,  5,  5,  5,  5 };
//        Float_t centralityAxisMapping[] =  { 2.5, 7.5, 100, 12.5, 17.5, 22.5, 27.5, 32.5, 37.5, 42.5};
//        Float_t centralityAxisMappingE[] = { 2.5, 2.5,   0,  2.5,  2.5,  2.5,  2.5,  2.5,  2.5,  2.5};
//        Float_t centralityAxisMapping[] =  { 5, 65, 100, 25, 15, 40};
//        Float_t centralityAxisMappingE[] = { 5, 15, 0,    5,  5, 10};

	Bool_t success = FitDeltaPhi2DOneFunction((TH2*) hist1, canvas, graphID, centralityAxisMapping[orderConvert[histId]], centralityAxisMappingE[orderConvert[histId]], orderConvert[histId]);
	if (!success)
        {
          binsFailed++;
	  checkList.Add(new TObjString(Form("AnalyzeDeltaPhiEtaGap2DExample(\"%s\", %d, %d, %d)", fileName, i, j, orderConvert[histId])));
        }
      }
    }
  }
  
  WriteGraphs(outputFileName);

  for (Int_t i=0; i<checkList.GetEntries(); i++)
    Printf("%s", checkList.At(i)->GetName());
  return binsFailed;
}

void AnalyzeExample(const char* fileName, Int_t i, Int_t j, Int_t histId, Bool_t drawDetails = kFALSE)
{
  if (energy == 1)
    vnmax = 5;
  CreateGraphStructure();

  TFile::Open(fileName);
  
  TH2* hist1 = (TH2*) gFile->Get(Form("dphi_%d_%d_%d", i, j+1, histId));
  if (!hist1)
  {
    cerr << "Histogram doesn't exist!" << endl;
    return;
  }

  Printf("Entries: %f %s", hist1->GetEntries(), hist1->GetTitle());

  TString label(hist1->GetTitle());
  if (drawDetails)
    hist1->SetTitle("");
  hist1->GetYaxis()->SetRangeUser(-1.79, 1.79);
  hist1->GetXaxis()->SetTitleOffset(1.5);
  hist1->GetYaxis()->SetTitleOffset(2);
  hist1->GetZaxis()->SetTitle(kCorrFuncTitle);
  hist1->GetZaxis()->SetTitleOffset(1.8);
  hist1->SetStats(kFALSE);
  hist1->Scale(1.0 / hist1->GetYaxis()->GetBinWidth(1));

//  if (hist1->GetEntries() < 1e4)
  if (hist1->GetEntries() < 10)
  {
    cerr << "Number of entries smaller than 10." << endl;
    return;
  }
  // NOTE fix normalization. these 2d correlations come out of AliUEHist normalized by dphi bin width, but not deta
  if (drawDetails)
    hist1->Scale(1.0 / hist1->GetYaxis()->GetBinWidth(1));

  TCanvas* canvas = new TCanvas(Form("DeltaPhi_%d_%d_%d_%d", histId, i, j, 1), Form("DeltaPhi_%d_%d_%d_%d", histId, i, j, 1), 1400, 1100);
  canvas->Divide(2, 2);
//  if (i == 4 || (i > 1 && histId > 5)) hist1 = RebinTTR(hist1,true);
  Int_t graphID = i * (6 - 1) + j - 1;
  FitDeltaPhi2DOneFunction((TH2*) hist1, canvas, graphID, 0, 0, histId);
  
}

void AnalyzeAll(int start=0, int finish=100)
{
  // 2.76 TeV data analysis
  const int nFiles = 19;
  gNBin = 0;
  gNBinFit = 0;
  string names[nFiles] = {"Data/Default/merge/","Data/Default/merge/","Data/Default/merge/","Data/Systematics/merge/def2/Original/","Data/Systematics/merge/def2/TTRLarger/","Data/Systematics/merge/def2/TTRTPCVolume/","Data/Systematics/merge/Eta0.7/","Data/Systematics/merge/Vertex/","Data/Systematics/merge/def1/Original/","Data/Systematics/merge/def3/Original/","Data/Systematics/merge/def3/PairCutHalf/","Data/Systematics/merge/def3/PairCutDouble/","Data/Systematics/merge/Global/","Data/Systematics/merge/def4/MagneticField/","Data/Systematics/merge/def4/Original/","AMPT/Eta0.8/Deta1.59/smON_resON/","AMPT/Eta0.8/Deta1.59/smON_resOFF/","AMPT/Eta0.8/Deta1.59/smOFF_resON/","Hijing/"};
  vector<Int_t> binsFailed(nFiles,0); 
  for (int i=0; i<nFiles; i++)
  {
    if (i < start || i > finish) continue;
    string dphiName = names[i] + "/dphi_corr_norm_wing_removed.root";
    string graphName = names[i] + "/graphs_wing_removed_noBetaLimit.root";
/*    if (i == 14) 
    {
      dphiName = names[i] + "/dphi_corr_norm_wing_removed_1-2pTt_1-2pTa.root";
      graphName = names[i] + "/graphs_wing_removed_1-2pTt_1-2pTa_noBetaLimit.root";
    }
    else if (i == 15) 
    {
      dphiName = names[i] + "/dphi_corr_norm_wing_removed_3-4pTt_2-3pTa.root";
      graphName = names[i] + "/graphs_wing_removed_3-4pTt_2-3pTa_noBetaLimit.root";
    }
*/    if (i == 0)
    {
      gNBin = 1;
      graphName = "Data/Systematics/merge/Exclusion1BinLarger/graphs_wing_removed_onlyForDip_noBetaLimit.root";
    }
    else if (i == 1)
    {
      gNBinFit = 1;
      graphName = "Data/Systematics/merge/Exclusion1BinLarger/graphs_wing_removed_onlyForFit_noBetaLimit.root";
    }
    binsFailed[i] = (Analyze(dphiName.c_str(),graphName.c_str()));
    gNBin = 0;
    gNBinFit = 0;
  }
  for (unsigned int i=0; i<binsFailed.size(); i++)
    if (binsFailed[i] != 0)
    {
      if (binsFailed[i] == -1)
        cerr << "No dphi file in " << names[i] << " with the correct name (i=" << i << ")" << endl;
      else
        cerr << "In " << names[i] << " (i=" << i << ") " << binsFailed[i] << " bins failed" << endl;
    }
}

Bool_t SkipGraph(Int_t i)
{
  return (i == 14);
}

TGraphAsymmErrors* FixGraph(TGraphErrors* graph, Float_t shift)
{
  graph->Sort();
  
  TGraphAsymmErrors* res = new TGraphAsymmErrors;
  res->SetTitle(graph->GetTitle());
  res->GetXaxis()->SetTitle(graph->GetXaxis()->GetTitle());
  res->GetYaxis()->SetTitle(graph->GetYaxis()->GetTitle());

  for (Int_t i=0; i<graph->GetN(); i++)
  {
    res->SetPoint(res->GetN(), graph->GetX()[i] + shift, graph->GetY()[i]);
    res->SetPointError(res->GetN()-1, graph->GetEX()[i], graph->GetEX()[i], graph->GetEY()[i], graph->GetEY()[i]);
    
    // asymmetric if space
    if (graph->GetEX()[i] > TMath::Abs(shift))
    {
      res->GetEXlow()[i] += shift;
      res->GetEXhigh()[i] -= shift;
    }
  }
  
  return res;
}

void GraphRemovePoint(TGraphAsymmErrors* graph1, TGraphAsymmErrors* graph2)
{
  if (graph1->GetN() != graph2->GetN())
  {
    graph1->Sort();
    graph2->Sort();
    Int_t Nbin = 0;
    if (graph1->GetN() < graph2->GetN()) Nbin = graph1->GetN();
    else Nbin = graph2->GetN();
    for (Int_t bin1 = 0; bin1 < Nbin; bin1++)
    {
      Float_t x1 = graph1->GetX()[bin1]; 
      Float_t x2 = graph2->GetX()[bin1];
      if (x1 != x2)
      {
        if (x1<x2) {graph1->RemovePoint(bin1--);cout << x1 << "\t" << x2 << "\t" << bin1 << "\t" << graph1->GetN() << "\t" << graph2->GetN()<< endl;}
        else {graph2->RemovePoint(bin1--);cout << x1 << "\t" << x2 << "\t" << bin1 << "\t" << graph1->GetN() << "\t" << graph2->GetN()<< endl;}
      }
      if (graph1->GetN() < graph2->GetN()) Nbin = graph1->GetN();
      else Nbin = graph2->GetN();
    }
  }
}

void PrepareGraphs(Int_t nHists, TGraphErrors** graph, TGraphErrors** systematicA, TGraphErrors** systematicB, TMultiGraph** multiGraph, TMultiGraph** multiGraphSyst, Int_t uncertaintyID, Float_t offset = 0)
{
  Int_t colors[17] =  { 1, 3, 2, 6, 4, 7, 8, 9, 11, 12, 28, 30, 36, 40, 46, kOrange-3 };
//  Int_t colors[] =  { 1,1,1, 3,3,3, 2,2,2, 6,6,6, 4,4,4, 7,7,7, 8,8,8, 9,9,9, 11,11,11, 12,12,12, 28,28,28, 30,30,30, 36,36,36, 40,40,40, 46,46,46 };
//  Int_t colors[16] =  { 1, kGreen-6, kRed-7, kMagenta-2, kGreen+3, kCyan+1, kBlue, kViolet-9, kGray+1, kOrange+1, 28, 30, 36, 40, 46 };
//  Int_t markers[] = { 24, 21, 27, 24, 21, 27,24, 21, 27,24, 21, 27,24, 21, 27,24, 21, 27,24, 21, 27,24, 21, 27,24, 21, 27,24, 21, 27,24, 21, 27,24, 21, 27,24, 21, 27,24, 21, 27,24, 21, 27};
  Int_t markers[17] = { 20, 21, 22, 23, 24, 25, 26, 27, 28, 30, 31, 32, 33, 34, 2, 5, 3 };
//  Int_t markers[16] = { 20, 21, 34, 31, 24, 25, 33, 27, 28, 30, 31, 32, 33, 34, 2, 5};
//  Int_t fillStyle[11] = { 3001, 3001, 3001, 3001, 3001, 3001, 3001, 3001, 3001, 3001, 3001 };
  Int_t fillStyle[11] = { 3008, 3008, 3008, 3008, 3008, 3008, 3008, 3008, 3008, 3008, 3008 };
//  Int_t fillStyle[11] = { 3001, 3002, 3003, 3004, 3005, 3006, 3007, 3008, 3009, 3010, 3011 };
  Int_t count = 0;
  
  if (*multiGraph == 0)
    *multiGraph = new TMultiGraph;
  if (*multiGraphSyst == 0)
    *multiGraphSyst = new TMultiGraph;
  
//  Float_t shift = -1 + offset;
  Float_t shift = 0;
//`  Float_t shift = -4 + offset;
  for (Int_t i=0; i<nHists; i++)
  {
    if (SkipGraph(i))
      continue;
    TGraphAsymmErrors* graphcentrality = FixGraph((TGraphErrors*) graph[i]->Clone(), shift);
    if (graphcentrality->GetN() <= 0)
      continue;
    
    if (systematicA)
    {
      TGraphAsymmErrors* graphsystematicsA = FixGraph((TGraphErrors*) systematicA[i]->Clone(), shift);
      
      if (graphcentrality->GetN() != graphsystematicsA->GetN())
      {
	Printf("Different number of points %d %d: %s", graphcentrality->GetN(), graphsystematicsA->GetN(), graphcentrality->GetTitle());
        if (graphcentrality->GetN() == 0 || graphsystematicsA->GetN() == 0) continue;
        GraphRemovePoint(graphcentrality,graphsystematicsA);
      }
      
      TGraphErrors* graphsystematicsB = 0;
      if (systematicB)
      {
	graphsystematicsB = (TGraphErrors*) systematicB[i]->Clone();
	graphsystematicsB->Sort();
      
	if (graphcentrality->GetN() != graphsystematicsB->GetN())
	{
	  Printf("Different number of points %d %d", graphcentrality->GetN(), graphsystematicsB->GetN());
          continue;
	}
      }
      
      for (Int_t j=0; j<graphsystematicsA->GetN(); j++)
      {
	// uncertaintyID
	Double_t yMin = graphcentrality->GetY()[j];
	Double_t yMax = graphcentrality->GetY()[j];
        Double_t y = graphcentrality->GetY()[j];

        if (uncertaintyID == 0 && i%3 == 1) //  sigma phi
        {
          yMin *= 0.981;
          yMax *= 1.019;
        }
        else if (uncertaintyID == 1 && i%3 == 1) // sigma eta
        {
          yMin *= 0.958;
          yMax *= 1.042;
        }
        else if (uncertaintyID == 2 && i%3 == 1) // depletion
        {
          if (j == 0) // 0-10%
          {
            yMin *= 0.759;
            yMax *= 1.241;
          }
          else if (j == 1) // 10-20%
          {
            yMin *= 0.735;
            yMax *= 1.265;
          }
          else if (j == 2) // 20-30%
          {
            yMin *= 0.674;
            yMax *= 1.326;
          }
          else if (j == 3) // 30-50%
          {
            yMin *= 0.550;
            yMax *= 1.450;
          }
          if (y-yMin < 0.003)
          {
            yMin = y-0.003;
            yMax = y+0.003;
          }
        }
	
	
	if (graphsystematicsA->GetY()[j] < graphcentrality->GetY()[j])
	  yMin = graphcentrality->GetY()[j] - TMath::Sqrt(TMath::Power(yMin - graphcentrality->GetY()[j], 2) + TMath::Power(graphsystematicsA->GetY()[j] - graphcentrality->GetY()[j], 2));
	if (graphsystematicsA->GetY()[j] > graphcentrality->GetY()[j])
	  yMax = graphcentrality->GetY()[j] + TMath::Sqrt(TMath::Power(yMax - graphcentrality->GetY()[j], 2) + TMath::Power(graphsystematicsA->GetY()[j] - graphcentrality->GetY()[j], 2));
	
	if (graphsystematicsB)
	{
	  if (graphsystematicsB->GetY()[j] < graphcentrality->GetY()[j])
	    yMin = graphcentrality->GetY()[j] - TMath::Sqrt(TMath::Power(yMin - graphcentrality->GetY()[j], 2) + TMath::Power(graphsystematicsB->GetY()[j] - graphcentrality->GetY()[j], 2));
	  if (graphsystematicsB->GetY()[j] > graphcentrality->GetY()[j])
	    yMax = graphcentrality->GetY()[j] + TMath::Sqrt(TMath::Power(yMax - graphcentrality->GetY()[j], 2) + TMath::Power(graphsystematicsB->GetY()[j] - graphcentrality->GetY()[j], 2));
	}

	graphsystematicsA->GetEYlow()[j] = graphcentrality->GetY()[j] - yMin;
	graphsystematicsA->GetEYhigh()[j] = yMax - graphcentrality->GetY()[j];
	graphsystematicsA->GetY()[j]  = graphcentrality->GetY()[j];
	
	graphsystematicsA->GetEXlow()[j] = 1;
	graphsystematicsA->GetEXhigh()[j] = graphsystematicsA->GetEXlow()[j];
      }
      
//       graphsystematicsA->SetFillColor(kGray);
      graphsystematicsA->SetFillColor(colors[count]);
//       graphsystematicsA->SetFillStyle(1001);
      graphsystematicsA->SetFillStyle(fillStyle[count]);
      graphsystematicsA->SetMarkerStyle(0);
      graphsystematicsA->SetLineColor(0);
      (*multiGraphSyst)->Add(graphsystematicsA, "2");
    }
    
    TString label = graphcentrality->GetTitle();
    if (label.Length() > 0)
    {
      cerr << label << endl;
      TObjArray* tokens = label.Tokenize("-");
      label.Form("%s-%s", tokens->At(0)->GetName(),tokens->At(1)->GetName());
      label.ReplaceAll(".00", "");
      label.ReplaceAll(".0", "");
      label.ReplaceAll("p_{T,trig}", "p_{T,t}");
      label.ReplaceAll("p_{T,assoc}", "p_{T,a}");
    }
    label += " GeV/c";
    graphcentrality->SetTitle(label);
    Printf("%d %s %d", i, label.Data(), colors[count]);
    
    graphcentrality->SetMarkerStyle(markers[count]);
    graphcentrality->SetMarkerColor(colors[count]);
    graphcentrality->SetLineColor(colors[count]);
    
    (*multiGraph)->Add(graphcentrality);
//    shift += 0.6;
//    if (shift > 0.7) shift = -1; 
    count++;
  }
}

//Bool_t drawLogo = kFALSE;
const char* MCLabel = 0;

void DrawCentrality(const char* canvasName, Int_t nHists, TGraphErrors** graph, Float_t min = 0, Float_t max = 0, const char* yLabel = "", TGraphErrors** systematicA = 0, TGraphErrors** systematicB = 0, TGraphErrors** graph2 = 0, TGraphErrors** systematic2 = 0, TGraphErrors** graph3 = 0, Int_t uncertaintyID = -1)
{
//   Bool_t found = kTRUE;
  TCanvas* c1 = 0;//(TCanvas*) gROOT->GetListOfCanvases()->FindObject(canvasName);
  if (!c1)
  {
    c1 = new TCanvas(canvasName, canvasName, 800, 600);
//     found = kFALSE;
  }
  c1->cd();
  c1->SetFillStyle(0); 
  TMultiGraph* multiGraph = 0;
  TMultiGraph* multiGraph2 = 0;
  TMultiGraph* multiGraphSyst = 0;
  
  PrepareGraphs(nHists, graph, systematicA, systematicB, &multiGraph, &multiGraphSyst, uncertaintyID);
  if (!multiGraph->GetListOfGraphs())
    return;

  if (graph2)
  {
    PrepareGraphs(nHists, graph2, systematic2, 0, &multiGraph2, &multiGraphSyst, uncertaintyID, 0.5);
    for (Int_t i=0; i<multiGraph2->GetListOfGraphs()->GetEntries(); i++)
      ((TGraph*)multiGraph2->GetListOfGraphs()->At(i))->SetLineWidth(2);
    
    while (multiGraph2->GetListOfGraphs()->GetEntries() > multiGraph->GetListOfGraphs()->GetEntries())
      multiGraph2->GetListOfGraphs()->RemoveAt(multiGraph->GetListOfGraphs()->GetEntries());
    
    TMultiGraph* multiGraph3 = 0;
    if (graph3)
    {
      PrepareGraphs(nHists, graph3, 0, 0, &multiGraph3, &multiGraphSyst, uncertaintyID);
      if (multiGraph3->GetListOfGraphs())
      {
	for (Int_t i=0; i<multiGraph3->GetListOfGraphs()->GetEntries(); i++)
	{
	  ((TGraph*)multiGraph3->GetListOfGraphs()->At(i))->SetLineWidth(2);
	  ((TGraph*)multiGraph3->GetListOfGraphs()->At(i))->SetLineStyle(2);
	}
	
	while (multiGraph3->GetListOfGraphs()->GetEntries() > multiGraph->GetListOfGraphs()->GetEntries())
	  multiGraph3->GetListOfGraphs()->RemoveAt(multiGraph->GetListOfGraphs()->GetEntries());
      }
    }
      
//     multiGraphSyst->Add(multiGraph2, "LX");
    
//     multiGraph->Add(multiGraph2);
    
    if (graph3)
      multiGraph2->Add(multiGraph3, "LX");
  }

//   multiGraphSyst->Add(multiGraph, (found) ? "LX" : "P");
  
  // draw by hand, multigraph draw spoils the pngs for some reason, same thing happens with this long code, have changed to gif below...
  
  TGraphErrors* first = (TGraphErrors*) multiGraph->GetListOfGraphs()->First();
  TH2* dummy = new TH2F("dummy", ";Centrality", 100, -2, 115, 100, min, max);
  dummy->SetStats(0);
  dummy->GetYaxis()->SetTitle(yLabel);
  dummy->GetYaxis()->SetTitleOffset(1.2);
  dummy->DrawCopy();
  
  if (multiGraphSyst->GetListOfGraphs())
  {
    multiGraphSyst->GetListOfGraphs()->First()->Draw((!first) ? "A2" : "2SAME");
    if (!first)
      first = (TGraphErrors*) multiGraphSyst->GetListOfGraphs()->First();
    
    for (Int_t i=1; i<multiGraphSyst->GetListOfGraphs()->GetEntries(); i++)
      multiGraphSyst->GetListOfGraphs()->At(i)->Draw("2SAME");
  }
  
  if (multiGraph2 && multiGraph2->GetListOfGraphs())
  {
    for (Int_t i=0; i<multiGraph2->GetListOfGraphs()->GetEntries(); i++)
    {
      TGraphAsymmErrors* graphTmp = (TGraphAsymmErrors*) multiGraph2->GetListOfGraphs()->At(i);
      for (Int_t j=0; j<graphTmp->GetN(); j++)
      {
	graphTmp->GetEXlow()[j] = 0;
	graphTmp->GetEXhigh()[j] = 0;
      }
      graphTmp->SetMarkerSize(0);
    }

    multiGraph2->GetListOfGraphs()->First()->Draw((!first) ? "ALZ" : "LXSAME");
    if (!first)
      first = (TGraphErrors*) multiGraph2->GetListOfGraphs()->First();
    
    for (Int_t i=1; i<multiGraph2->GetListOfGraphs()->GetEntries(); i++)
      multiGraph2->GetListOfGraphs()->At(i)->Draw("LZSAME");
  }

  if (multiGraph && multiGraph->GetListOfGraphs())
  {
    multiGraph->GetListOfGraphs()->First()->Draw((!first) ? "AP" : "PSAME");
    if (!first)
      first = (TGraphErrors*) multiGraph->GetListOfGraphs()->First();
    
    for (Int_t i=1; i<multiGraph->GetListOfGraphs()->GetEntries(); i++)
      multiGraph->GetListOfGraphs()->At(i)->Draw("PSAME");
  }

//  TPaveText* paveText = new TPaveText(0.826, 0.059, 0.895, 0.094, "BRNDC");
  TPaveText* paveText = new TPaveText(0.77, 0.059, 0.84, 0.094, "BRNDC");
  paveText->SetTextSize(0.04);
  paveText->SetFillColor(0);
  paveText->SetShadowColor(0);
  paveText->SetLineColor(0);
  paveText->AddText("pp");
  paveText->Draw();
    
  TLegend* legend = new TLegend(0.55, 0.65, 0.95, 0.95);
  legend->SetFillColor(0);
  legend->SetTextSize(0.03);
  
  if (graph2)
  {
    legend->SetNColumns(2);
    legend->SetHeader(MCLabel);
  }
  
//  for (Int_t i=1; i<multiGraph->GetListOfGraphs()->GetEntries(); i+=3)
  for (Int_t i=0; i<multiGraph->GetListOfGraphs()->GetEntries(); i++)
  {
    if (graph2)
    {
      legend->AddEntry(multiGraph->GetListOfGraphs()->At(i), " ", "P");
      legend->AddEntry(multiGraph2->GetListOfGraphs()->At(i), 0, "L");
    }
    else
      legend->AddEntry(multiGraph->GetListOfGraphs()->At(i), 0, "P");
  }
  legend->Draw();
/*  TLegend* legend2 = new TLegend(0.37, 0.74, 0.52, 0.87);
  legend2->SetFillColor(0);
  legend2->SetTextSize(0.03);

  legend2->AddEntry(multiGraph->GetListOfGraphs()->At(0), "2010", "P");
  legend2->AddEntry(multiGraph->GetListOfGraphs()->At(1), "Merged", "P");
  legend2->AddEntry(multiGraph->GetListOfGraphs()->At(2), "2011", "P");
  legend2->Draw();
*/
  gPad->SetGridx();
  gPad->SetGridy();
  
  if (energy == 0)
  {
    DrawLatex(0.13, 0.85,  1, "Pb-Pb #sqrt{s_{NN}} = 2.76 TeV", 0.03);
    DrawLatex(0.13, 0.81,  1, "pp #sqrt{s} = 2.76 TeV", 0.03);
    DrawLatex(0.13, 0.77, 1, "|#eta| < 0.8", 0.03);
  }
  else if (energy == 1)
  {
    DrawLatex(0.13, 0.85,  1, "Pb-Pb #sqrt{s_{NN}} = 5.02 TeV", 0.03);
    DrawLatex(0.13, 0.81,  1, "pp #sqrt{s} = 5.02 TeV", 0.03);
    DrawLatex(0.13, 0.77, 1, "|#eta| < 0.8", 0.03);
  }
  TString text;
  TString text2;
  if (TString(canvasName).BeginsWith("sigma_phi"))
  {
    text = "Projected within |#Delta#eta| < 0.80";
    text2 = "Calculated within |#Delta#varphi| < 0.87";
  }
  if (TString(canvasName).BeginsWith("sigma_eta"))
  {
    text = "Projected within |#Delta#varphi| < 0.87";
    text2 = "Calculated within |#Delta#eta| < 0.80";
  }
  if (text.Length() > 0)
  {
    DrawLatex(0.13, 0.73, 1, text, 0.03);
    DrawLatex(0.13, 0.69, 1, text2, 0.03);
  }
  
  if (0 && MCLabel)
  {
    DrawLatex(0.13, 0.18, 1, "Points: Data", 0.03);
    DrawLatex(0.13, 0.14, 1, Form("Lines: %s", MCLabel), 0.03);
  }
}

Float_t** ExtractSystematics(const char* baseFile, const char* systFile, Int_t offset=0, const char* name=NULL)
{
  ReadGraphs(baseFile);
 
  TGraphErrors*** graphsBase = graphs;
  if (systFile)
  {
    ReadGraphs(systFile);
  }
  
  // calculate syst unc for these graphs
  const Int_t NGraphList = 3;
  Int_t graphList[] = { 26, 27, 11, 12, 13};
  string titles[] = { "#Delta#varphi", "#Delta#eta", "Depletion yield", "#sigma_{CP, #Delta#varphi}", "#sigma_{CP, #Delta#eta}"};

  Float_t** results = new Float_t*[NGraphList];

  for (Int_t i=0; i<NGraphList; i++)
  {
    results[i] = new Float_t[2];
    TGraphErrors* graphsTmp[NHists] = {0};


    
    TH1* hist = 0;
    hist = new TH1F(Form("hist_%d", i), "", 50, 0.5, 1.5);  
  
    TCanvas* c = new TCanvas(Form("%s_%d_%d", systFile, i, 0), Form("%s_%d_%d", systFile, i, 0), 1000, 1000);
    c->Divide(4, 5);

    Int_t count = 1;
    for (Int_t j=0; j<NHists; j++)
    {
      if (SkipGraph(j))
	continue;
      
      TGraphErrors* graph1 = graphsBase[graphList[i]+offset][j];
      TGraphErrors* graph2 = (systFile) ? graphs[graphList[i]+offset][j] : graphsBase[graphList[i]+16][j];
      
      if (graph1->GetN() == 0)
	continue;

//      if (i == 2 && j > 5 ) continue; 
      if (i == 2 && j > 10 ) continue; 
//      if (j%5 >= 2) continue;

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
      
      DivideGraphs(graph1, graph2);
      if (i != 2) graph1->GetYaxis()->SetRangeUser(0.9, 1.1);
      if (i == 2) (TGraphErrors*) graph1->DrawClone("AP");
      else ((TGraphErrors*) graph1->DrawClone("AP"))->GetYaxis()->SetRangeUser(0.8, 1.2);
      graphsTmp[j] = (TGraphErrors*) graph1->DrawClone();
      
      for (Int_t k=0; k<graph1->GetN(); k++)
      {
        if (graph1->GetX()[k] == 100) continue;
        hist->Fill(graph1->GetY()[k]);
      }
      
      if (count == 37)
	break;
    }
    
    if (name != NULL) c->SaveAs(Form("%s_%d.pdf", name, i));
    
    new TCanvas;
    hist->Draw();
    hist->Sumw2();
    
    Float_t mean = -1;
    Float_t sigma = -1;
    if (0)
    {
      hist->Fit("gaus", "Q");
      mean = hist->GetFunction("gaus")->GetParameter(1);
      sigma = hist->GetFunction("gaus")->GetParameter(2);
      
      Printf("%d: %.3f +- %.3f %.3f | %.3f +- %.3f", i, mean, sigma, sigma / TMath::Sqrt(hist->GetEntries()), hist->GetMean(), hist->GetMean(11));
    }
    else
    {
      mean = hist->GetMean(1);
      sigma = hist->GetMean(11);

      Printf("%d: %.3f +- %.3f", i, hist->GetMean(), hist->GetMean(11));
    }
    
    mean -= 1;
    mean *= 100;
    sigma *= 100;
    Printf("==> %.1f +- %.1f", mean, sigma);
    results[i][0] = mean;
    results[i][1] = sigma / TMath::Sqrt(hist->GetEntries());
  }
  
  return results;
}

void ExtractSystematicsAll(Int_t offset = 0)
{
  gROOT->SetBatch(kTRUE);
  
  const Int_t NEffects = 11;
  
  const char* defaultFile[] = {"Data/Systematics/merge/AllFilesMerged/graphs_wing_removed_noBetaLimit.root", "Data/Systematics/merge/def2/Original/graphs_wing_removed_noBetaLimit.root", "Data/Systematics/merge/def2/Original/graphs_wing_removed_noBetaLimit.root", "Data/Systematics/merge/def3/Original/graphs_wing_removed_noBetaLimit.root", "Data/Systematics/merge/def3/Original/graphs_wing_removed_noBetaLimit.root", "Data/Systematics/merge/AllFilesMerged/graphs_wing_removed_noBetaLimit.root", "Data/Systematics/merge/AllFilesMerged/graphs_wing_removed_noBetaLimit.root", "Data/Systematics/merge/AllFilesMerged/graphs_wing_removed_noBetaLimit.root", "Data/Default/merge/graphs_wing_removed_noBetaLimit.root", "Data/Default/merge/graphs_wing_removed_noBetaLimit.root", "Data/Systematics/merge/def4/Original/graphs_wing_removed_noBetaLimit.root"};
  const char* systFiles[] = {"Data/Systematics/merge/Global/graphs_wing_removed_noBetaLimit.root", "Data/Systematics/merge/def2/TTRLarger/graphs_wing_removed_noBetaLimit.root", "Data/Systematics/merge/def2/TTRTPCVolume/graphs_wing_removed_noBetaLimit.root", "Data/Systematics/merge/def3/PairCutHalf/graphs_wing_removed_noBetaLimit.root", "Data/Systematics/merge/def3/PairCutDouble/graphs_wing_removed_noBetaLimit.root", "Data/Systematics/merge/Vertex/graphs_wing_removed_noBetaLimit.root", "Data/Systematics/merge/Eta0.7/graphs_wing_removed_noBetaLimit.root", "Data/Systematics/merge/def1/Original/graphs_wing_removed_noBetaLimit.root", "Data/Systematics/merge/Exclusion1BinLarger/graphs_wing_removed_onlyForFit_noBetaLimit.root", "Data/Systematics/merge/Exclusion1BinLarger/graphs_wing_removed_onlyForDip_noBetaLimit.root", "Data/Systematics/merge/def4/MagneticField/graphs_wing_removed_noBetaLimit.root"};


  const char* effectNames[] = {"Track selection and efficiencies", "Small opening angles cut 1", "Small opening angles cut 2", "Neutral-particle decay cut - half", "Neutral-particle decay cut - double", "Vertex range", "$|\\eta|<0.7$", "$|\\eta|<0.9$", "Exclusion region from fit", "Exclusion region from dip region", "Magnetic field" };
  const char* outputFilenames[] = {"Global", "TTRLarger", "TTRTPCVolume", "PairCutHalf", "PairCutDouble", "Vertex", "Eta07", "Eta09", "Exclusion_fit", "Exclusion_dip", "Magnetic_field"};

  Float_t** results[NEffects+2];
  
  for (Int_t i=0; i<NEffects; i++)
  {
    cout << i << ": " <<systFiles[i] << " (" << defaultFile[i] << ") " << endl;
    results[i] = ExtractSystematics(defaultFile[i], systFiles[i], offset, outputFilenames[i]);
  }
  
  const Int_t NParameters = 5;
  const char* names[] = { "$\\sigma_{\\Delta\\varphi}$", "$\\sigma_{\\Delta\\eta}$", "Depletion yield", "$\\sigma_{CP,\\Delta\\varphi}$", "$\\sigma_{CP,\\Delta\\eta}$"};
  // combine in quadrature (only mean)
  results[NEffects+2] = new Float_t*[NParameters];
  for (Int_t j=0; j<NParameters; j++)
  {
    results[NEffects+2][j] = new Float_t[2];
    Float_t mean = 0;
  
    printf("%s \t:", names[j]);
    for (Int_t i=0; i<NEffects; i++)
    {
      mean += results[i][j][0] * results[i][j][0];
      printf("%.1f%% ", results[i][j][0]);
    }
    mean = TMath::Sqrt(mean);
    Printf("--> %.1f%% \t", mean);
    results[NEffects+2][j][0] = mean;
  }
  cerr << endl << endl;
  cerr << "Source" << " &\t";
  for (Int_t j=0; j<NParameters; j++)
    if (j != NParameters-1) cerr <<  names[j] << " &" << "\t";
    else cerr <<  names[j] << "\\\\ \\toprule" << endl;
  for (Int_t i=0; i<NEffects+1; i++)
  {
    if (i == NEffects) cerr << " \\midrule  Total &" << "\t";
    else cerr << effectNames[i] << " &" << "\t";
    for (Int_t j=0; j<NParameters; j++)
    {
      if (i == NEffects) 
      {
        if (j != NParameters-1) printf("%.1f\\%% & ", results[NEffects+2][j][0]);
        else printf("%.1f\\%% \\\\", results[NEffects+2][j][0]);
        continue;
      }
      if (results[i][j][0] >= 0)
      {
        if (j != NParameters-1) printf("%.1f\\%% & ", results[i][j][0]);
        else printf("%.1f\\%% \\\\", results[i][j][0]);
      }
      else
      {
        if (j != NParameters-1) printf("%.1f\\%% & ", -1*results[i][j][0]);
        else printf("%.1f\\%% \\\\", -1*results[i][j][0]);
      }
    }
    cerr << endl;
  }
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

void WingToy(Float_t centralityBegin, Float_t centralityEnd, Float_t trigBegin, Float_t trigEnd, Float_t assocBegin, Float_t assocEnd)
{
  // calculate same/mixed event in absence of physical correlations
  
  TFile::Open("yields_120501.root");
  TH3* hist = (TH3*) gFile->Get("fYields");
  hist->Scale(1.0 / 8891e3);
  
  TH1* same = new TH1D("same", "", 41, -2.05, 2.05);
  TH1* mixed = new TH1D("mixed", "", 41, -2.05, 2.05);
  
  TH1* etaDistArrTrig[hist->GetNbinsX()+1];
  TH1* etaDistArrAssoc[hist->GetNbinsX()+1];
  for (Int_t i=1; i<=hist->GetNbinsX(); i++)
  {
    hist->GetXaxis()->SetRange(i, i);
    hist->GetYaxis()->SetRangeUser(trigBegin+0.01, trigEnd-0.01);
    etaDistArrTrig[i] = hist->Project3D(Form("z%d", i*2));
    hist->GetYaxis()->SetRangeUser(assocBegin+0.01, assocEnd-0.01);
    etaDistArrAssoc[i] = hist->Project3D(Form("z%d", i*2+1));
  }
  
  new TCanvas;
  etaDistArrTrig[1]->DrawCopy()->Scale(1.0 / etaDistArrTrig[1]->Integral());
  etaDistArrTrig[20]->SetLineColor(2); etaDistArrTrig[20]->DrawCopy("SAME")->Scale(1.0 / etaDistArrTrig[20]->Integral());
  etaDistArrTrig[70]->SetLineColor(4); etaDistArrTrig[70]->DrawCopy("SAME")->Scale(1.0 / etaDistArrTrig[70]->Integral());
  
  TH1* lastEtaDist = 0;
//   Float_t lastFactor = 0;
  Int_t nEvents = 5000;
  for (Int_t i=0; i<nEvents; i++)
  {
    if (i % 1000 == 0)
      Printf("%d", i);
    Float_t centrality = gRandom->Uniform(centralityBegin, centralityEnd);
//     Printf("centrality %f", centrality);
    
    TH1* etaDistTrig = etaDistArrTrig[hist->GetYaxis()->FindBin(centrality)];
    TH1* etaDistAssoc = etaDistArrAssoc[hist->GetYaxis()->FindBin(centrality)];
//     TH1* etaDist = etaDistArr[1];
//     Float_t factor = etaDistArr[hist->GetYaxis()->FindBin(centrality)]->Integral() / etaDistArr[1]->Integral();
//     Printf("%f", factor);
//     new TCanvas; etaDistTrig->Draw(); etaDistAssoc->Draw("SAME"); etaDistAssoc->SetLineColor(2); return;
    
    // convolution, same
    for (Int_t j=1; j<=etaDistTrig->GetNbinsX(); j++)
      for (Int_t k=1; k<=etaDistAssoc->GetNbinsX(); k++)
      {
	Double_t weight = etaDistTrig->GetBinContent(j) * etaDistAssoc->GetBinContent(k);
	Double_t deltaEta = etaDistTrig->GetBinCenter(j) - etaDistAssoc->GetBinCenter(k);
// 	Printf("%f", deltaEta);
	
	same->Fill(deltaEta, weight);
      }
      
    // convolution, mixed
    if (lastEtaDist)
    {
      for (Int_t j=1; j<=etaDistTrig->GetNbinsX(); j++)
	for (Int_t k=1; k<=lastEtaDist->GetNbinsX(); k++)
	{
	  Double_t weight = etaDistTrig->GetBinContent(j) * lastEtaDist->GetBinContent(k);
	  Double_t deltaEta = etaDistTrig->GetBinCenter(j) - lastEtaDist->GetBinCenter(k);
	  
	  mixed->Fill(deltaEta, weight);
	}
    }
    lastEtaDist = etaDistAssoc;
//     lastFactor = factor;
    
//     break;
  }
  
  same->Scale(1.0 / nEvents);
  mixed->Scale(1.0 / (nEvents - 1));
  
  same->GetXaxis()->SetRangeUser(-1.6, 1.6);
  new TCanvas; same->DrawCopy(); mixed->SetLineColor(2); mixed->DrawCopy("SAME");
  new TCanvas; same->Divide(mixed); same->Draw();
}

void WingToyUntriggered(Float_t centralityBegin, Float_t centralityEnd)
{
  // calculate same/mixed event in absence of physical correlations
  
  TFile::Open("EtaCent.root");
  TH2* hist = (TH2*) gFile->Get("hist");
  hist->Scale(1e-7 * 500);
  
/*  TH1* same = new TH1D("same", "", 201, -2.01, 2.01);
  TH1* mixed = new TH1D("mixed", "", 201, -2.01, 2.01);*/
  TH1* same = new TH1D("same", "", 41, -2.05, 2.05);
  TH1* mixed = new TH1D("mixed", "", 41, -2.05, 2.05);
  
  TH1* etaDistArr[hist->GetNbinsY()+1];
  for (Int_t i=1; i<=hist->GetNbinsY(); i++)
  {
    etaDistArr[i] = hist->ProjectionX(Form("etaDist_%d", i), i, i);
    etaDistArr[i]->Rebin(5);
  }
  
  TH1* lastEtaDist = 0;
//   Float_t lastFactor = 0;
  Int_t nEvents = 1000;
  for (Int_t i=0; i<nEvents; i++)
  {
    if (i % 1000 == 0)
      Printf("%d", i);
    Float_t centrality = gRandom->Uniform(centralityBegin, centralityEnd);
//     Printf("centrality %f", centrality);
    
    TH1* etaDist = etaDistArr[hist->GetYaxis()->FindBin(centrality)];
//     TH1* etaDist = etaDistArr[1];
//     Float_t factor = etaDistArr[hist->GetYaxis()->FindBin(centrality)]->Integral() / etaDistArr[1]->Integral();
//     Printf("%f", factor);
//     new TCanvas; etaDist->Draw();
    
    // convolution, same
    for (Int_t j=1; j<=etaDist->GetNbinsX(); j++)
      for (Int_t k=1; k<=etaDist->GetNbinsX(); k++)
      {
	Double_t weight = etaDist->GetBinContent(j) * etaDist->GetBinContent(k);
	Double_t deltaEta = etaDist->GetBinCenter(j) - etaDist->GetBinCenter(k);
// 	Printf("%f", deltaEta);
	
	same->Fill(deltaEta, weight);
      }
      
    // convolution, mixed
    if (lastEtaDist)
    {
      for (Int_t j=1; j<=etaDist->GetNbinsX(); j++)
	for (Int_t k=1; k<=etaDist->GetNbinsX(); k++)
	{
	  Double_t weight = etaDist->GetBinContent(j) * lastEtaDist->GetBinContent(k);
	  Double_t deltaEta = etaDist->GetBinCenter(j) - lastEtaDist->GetBinCenter(k);
	  
	  mixed->Fill(deltaEta, weight);
	}
    }
    lastEtaDist = etaDist;
//     lastFactor = factor;
    
//     break;
  }
  
  same->Scale(1.0 / nEvents);
  mixed->Scale(1.0 / (nEvents - 1));
  
  same->GetXaxis()->SetRangeUser(-1.8, 1.8);
  new TCanvas; same->DrawCopy(); mixed->SetLineColor(2); mixed->DrawCopy("SAME");
  new TCanvas; same->Divide(mixed); same->Draw();
}

void RapidityToy()
{
  // calculate effect on same and mixed event due to pseudorapidity dip at 0 from flat rapidity
  
  TH1* same = new TH1D("same", "", 201, -2.01, 2.01);
  TH1* mixed = new TH1D("mixed", "", 201, -2.01, 2.01);
  
  TH1* etaDist = new TH1D("etaDist", "", 201, -2.01, 2.01);

  TF1* yToEta = new TF1("yToEta", "acosh((exp(-x)*sqrt([0]^4 - 2*exp(2*x)*[0]^4 + exp(4*x)*[0]^4 + 2*[0]^2*[1]^2 + 2*exp(4*x)*[0]^2*[1]^2 + [1]^4 + 2*exp(2*x)*[1]^4 + exp(4*x)*[1]^4))/([1]*sqrt(4*[0]^2 + 4*[1]^2)))");
  Double_t pt = 1.5;
//   Double_t mass = 1;
    Double_t mass = 0.140;
  yToEta->SetParameters(mass, pt);
    
  Int_t nEvents = 1000000;
  Double_t lastEta = 0;
  for (Int_t i=0; i<nEvents; i++)
  {
    if (i % 10000 == 0)
      Printf("%d", i);
    
    Double_t y = gRandom->Uniform(-3, 3);
    
    Double_t eta = yToEta->Eval(y);
    if (y < 0)
      eta *= -1;
    etaDist->Fill(eta);
    
    Double_t y2 = y + gRandom->Gaus(0, 0.5);
    Double_t eta2 = yToEta->Eval(y2);
    if (y2 < 0)
      eta2 *= -1;
    
    if (abs(eta) < 1 && abs(eta2) < 1)
      same->Fill(eta - eta2);
    if (i > 0 && abs(lastEta) < 1 && abs(eta2) < 1)
      mixed->Fill(lastEta - eta2);
    lastEta = eta;
  }

  new TCanvas; etaDist->Draw();

  same->Scale(1.0 / nEvents);
  mixed->Scale(1.0 / (nEvents - 1));
  
  same->GetXaxis()->SetRangeUser(-1.8, 1.8);
  new TCanvas; same->DrawCopy(); mixed->SetLineColor(2); mixed->DrawCopy("SAME");
  new TCanvas; same->Divide(mixed); same->Draw();
}

void v2DependencyToy(const char* outputFileName = "hist.root")
{
  TCanvas* c1 = new TCanvas("c1", "Fit", 800, 600);
  TCanvas* c2 = new TCanvas("c2", "Correlation", 1400, 1600);
  TCanvas* c3 = new TCanvas("c3", "dphi", 1400, 1600);
  TCanvas* c4 = new TCanvas("c4", "deta", 1400, 1600);
  c2->Divide(3,3);
  c3->Divide(3,3);
  c4->Divide(3,3);
  TFile::Open(outputFileName, "RECREATE");

  Double_t eta[] = { -3.25, -2.75, -2.25, -1.75, -1.25, -0.75, -0.25, 0.25, 0.75, 1.25, 1.75, 2.25, 2.75, 3.25, 3.75, 4.25, 4.75 };
  Double_t v20[] = { 0.0178046, 0.0210719, 0.0241588, 0.0228301, 0.0224484, 0.0224698, 0.0228943, 0.0228137, 0.0226188, 0.0220402, 0.0225268, 0.0231524, 0.0209274, 0.0190223, 0.0171396, 0.015942, 0.01385 } ;
  Double_t v2_err0[] = { 0.00188915, 0.00223662, 0.00256323, 0.00153577, 0.00158735, 0.00158886, 0.00161888, 0.00161319, 0.00159943, 0.00155848, 0.00155241, 0.00245647, 0.00222054, 0.00201772, 0.00181847, 0.00169292, 0.00147159 };
  Double_t v21[] = { 0.0289574, 0.0329946, 0.0376679, 0.0359883, 0.0357591, 0.0363211, 0.037169, 0.037344, 0.0365328, 0.0352206, 0.0360468, 0.036063, 0.03326, 0.0310704, 0.0278479, 0.0260765, 0.0233602 };
  Double_t v2_err1[] = { 0.00307179, 0.00350245, 0.0039966, 0.00246431, 0.00252859, 0.00256839, 0.00262837, 0.00264062, 0.00258332, 0.0024905, 0.00242361, 0.00382609, 0.00352801, 0.00329574, 0.00295455, 0.00276704, 0.00247998 };
  Double_t v22[] = { 0.0421916, 0.046932, 0.0529053, 0.0505833, 0.0507496, 0.0521991, 0.0535482, 0.0533636, 0.0525199, 0.0505711, 0.0501678, 0.0510414, 0.0470683, 0.0445314, 0.0399421, 0.0372192, 0.0335714 };
  Double_t v2_err2[] = { 0.00447529, 0.00498069, 0.00561346, 0.00350182, 0.00358878, 0.00369109, 0.00378656, 0.00377339, 0.00371389, 0.00357593, 0.00345671, 0.00541562, 0.00499257, 0.00472346, 0.00423742, 0.00394818, 0.00356289 };
  Double_t v23[] = { 0.0535116, 0.0592228, 0.0667268, 0.0643815, 0.0648137, 0.0671077, 0.0689489, 0.0688137, 0.0674753, 0.0649329, 0.064237, 0.0644013, 0.0597028, 0.0561042, 0.0504783, 0.0469539, 0.0432138 };
  Double_t v2_err3[] = { 0.00567636, 0.00628471, 0.0070801, 0.0044467, 0.00458304, 0.00474525, 0.0048755, 0.00486587, 0.00477138, 0.00459154, 0.00440783, 0.00683316, 0.00633296, 0.00595079, 0.00535544, 0.00498189, 0.00458449 };
  Double_t v24[] = { 0.0601229, 0.0657066, 0.0732884, 0.0714449, 0.0732558, 0.0757839, 0.0775992, 0.0777084, 0.0759983, 0.072938, 0.07129, 0.0714013, 0.0660398, 0.0622241, 0.0556125, 0.0515049, 0.0467412 };
  Double_t v2_err4[] = { 0.0063779, 0.00697181, 0.00777585, 0.00498088, 0.00518037, 0.00535884, 0.00548712, 0.00549488, 0.00537406, 0.00515754, 0.00498184, 0.00757661, 0.00700521, 0.00660028, 0.00590042, 0.00546399, 0.00495893 };
  Double_t v25[] = { 0.0606421, 0.0666128, 0.0747658, 0.0732977, 0.0758675, 0.0785067, 0.0804069, 0.0804495, 0.078766, 0.0756848, 0.0737924, 0.0727082, 0.0669812, 0.0627843, 0.0550477, 0.0503637, 0.0455578 };
  Double_t v2_err5[] = { 0.00643316, 0.00706748, 0.0079338, 0.00510717, 0.00536486, 0.00555131, 0.00568569, 0.00568868, 0.00556979, 0.00535181, 0.00512039, 0.00771524, 0.00710494, 0.00665958, 0.005842, 0.00534359, 0.00483393 };
  Double_t v26[] = { 0.0558427, 0.0622843, 0.069759, 0.0699102, 0.0726565, 0.0752472, 0.0772998, 0.0769496, 0.0747283, 0.0731016, 0.0692551, 0.0687594, 0.0630636, 0.0587596, 0.0476679, 0.0433981, 0.0384033 };
  Double_t v2_err6[] = { 0.00592555, 0.00660726, 0.00740837, 0.00497262, 0.00513763, 0.00532121, 0.00546618, 0.00544116, 0.00528409, 0.00516934, 0.0047711, 0.00730395, 0.00668902, 0.0062329, 0.00506006, 0.00460362, 0.00407592 };
  Double_t v27[] = { 0.0485678, 0.0522194, 0.0614082, 0.0617957, 0.067403, 0.069292, 0.0694255, 0.071401, 0.0675572, 0.0677939, 0.0619343, 0.0545356, 0.0462625, 0.0492458, 0.0383392, 0.0321335, 0.0300904 } ;
  Double_t v2_err7[] = { 0.00549521, 0.00590814, 0.0069551, 0.00481929, 0.00524359, 0.00538968, 0.00540076, 0.00555398, 0.00525709, 0.00527562, 0.00471932, 0.0062522, 0.0052376, 0.00557231, 0.00435881, 0.0036367, 0.00340626 };
  Double_t v28[] = { 0.0376917, 0.0425199, 0.0485303, 0.0491033, 0.0579609, 0.0625405, 0.0625229, 0.0626231, 0.0624147, 0.0586143, 0.0527063, 0.0536474, 0.0399211, 0.0385768, 0.0259235, 0.0213431, 0.019931 };
  Double_t v2_err8[] = { 0.00432739, 0.00481277, 0.0054979, 0.00406909, 0.00451077, 0.00486628, 0.00486314, 0.00487131, 0.00485474, 0.00455912, 0.00385339, 0.00611525, 0.00451861, 0.00436453, 0.0029794, 0.00241472, 0.00229507 };
  Double_t v2[9][17];
  Double_t v2_err[9][17];
  for (Int_t i=0;i<17;i++)
  { 
    v2[0][i]=v20[i];
    v2[1][i]=v21[i];
    v2[2][i]=v22[i];
    v2[3][i]=v23[i];
    v2[4][i]=v24[i];
    v2[5][i]=v25[i];
    v2[6][i]=v26[i];
    v2[7][i]=v27[i];
    v2[8][i]=v28[i];
    v2_err[0][i]=v2_err0[i];
    v2_err[1][i]=v2_err1[i];
    v2_err[2][i]=v2_err2[i];
    v2_err[3][i]=v2_err3[i];
    v2_err[4][i]=v2_err4[i];
    v2_err[5][i]=v2_err5[i];
    v2_err[6][i]=v2_err6[i];
    v2_err[7][i]=v2_err7[i];
    v2_err[8][i]=v2_err8[i];
  }
  TH2D* hist = new TH2D("hist","hist",30,-TMath::Pi(),TMath::Pi(),30,-1.8,1.8);
  TH2D* hist2 = new TH2D("hist2","hist2",30,-TMath::Pi(),TMath::Pi(),30,-1.8,1.8);
  for(Int_t i=0; i<9; i++)
  {
    cerr << i << ". entrality bin" << endl;
    TGraphErrors* gr = new TGraphErrors(17, eta, v2[i], 0, v2_err[i]);
    gr->SetMarkerStyle(28-i);
    gr->SetMarkerColor((i==4) ? 15 :9-i);
    gr->SetLineColor((i==4) ? 15 :9-i);
    gr->GetXaxis()->SetRangeUser(-1.3,1.3);
    gr->GetYaxis()->SetRangeUser(0,0.1);
    c1->cd();
    gr->DrawClone((i==0) ? "AP" : "PSAME");
    for (Int_t j=0;j<17;j++)
    {
      if (TMath::Abs(gr->GetX()[j])>1)
        gr->SetPoint(j,gr->GetX()[j],gr->GetY()[j]-gr->GetErrorY(j));
      else if (TMath::Abs(gr->GetX()[j])<0.5)
        gr->SetPoint(j,gr->GetX()[j],gr->GetY()[j]+gr->GetErrorY(j));
      gr->SetPointError(j,0,0);
    }
    gr->DrawClone("PSAME");
    TF1* func = new TF1("func","[0]+[1]*x+[2]*x**2", -1.3,1.3);
    gr->Fit(func,"R","");  
    func->SetLineColor(gr->GetLineColor());
    func->Draw("SAME");
    hist->Reset();
    hist2->Reset();
    TF1* funcArray[30];
    for (Int_t j=0;j<30;j++)
    {
      TF1* func2 = new TF1("func2","1+[0]*2*TMath::Cos(2*x)",-TMath::Pi(),TMath::Pi());
      func2->FixParameter(0,func->Eval(-0.9+j*0.06));
      funcArray[j] = func2;
    }
    Int_t nPart = 10;
    Double_t eta1[nPart];
    Double_t phi1[nPart];
    Double_t phi1m[nPart];
    for (Int_t j=0; j<nPart;j++)
    {
      eta1[j] = gRandom->Uniform(-0.9, 0.9);
      Int_t funcIndex = (Int_t) ((eta1[j]+0.9)/1.8*30);
      phi1[j] = funcArray[(funcIndex!=30) ? funcIndex: 29]->GetRandom();
      phi1m[j] = gRandom->Uniform(-TMath::Pi(),TMath::Pi());
    }
    for (Int_t j=0; j<nPart;j++)
    {
      if (j%1000==0) cerr << j << "\t";
      for (Int_t k=j+1; k<nPart;k++)
      {
        Double_t dphi = phi1[j]-phi1[k];
        Double_t dphim = phi1m[j]-phi1m[k];
        if (dphi > TMath::Pi()) dphi = dphi - TMath::TwoPi();
        else if (dphi < -TMath::Pi()) dphi = dphi + TMath::TwoPi();
        if (dphim > TMath::Pi()) dphim = dphim - TMath::TwoPi();
        else if (dphim < -TMath::Pi()) dphim = dphim + TMath::TwoPi();
        Double_t deta = eta1[j]-eta1[k];
        hist->Fill(dphi,deta);
        hist2->Fill(dphim,deta);
      }
    }
    hist->Sumw2();
    hist2->Sumw2();
    hist->Divide(hist2);
    hist->GetXaxis()->SetTitle(Form("#Delta#varphi"));
    hist->GetYaxis()->SetTitle(Form("#Delta#eta"));
    cout << endl;
    c2->cd(i+1);
    hist->DrawCopy("SURF1");
    hist->Write(Form("hist_%d", i));
    c3->cd(i+1);
    TH1* projPhi = hist->ProjectionX("dphi projection");
    projPhi->Scale(1.0/projPhi->GetNbinsX());
    projPhi->DrawCopy();
    projPhi->Write(Form("projPhi_%d", i));
    c4->cd(i+1);
    TH1D* projEta = hist->ProjectionY("deta projection",hist->GetXaxis()->FindBin(-1),hist->GetXaxis()->FindBin(1));
    Float_t etaBins = hist->GetXaxis()->FindBin(1) - hist->GetXaxis()->FindBin(-1) + 1;
    projEta->Scale(1.0/etaBins);
    projEta->DrawCopy();
    projEta->Write(Form("projEta_%d", i));
  }
  gFile->Close();
}

void v2DependencyToy2(const char* DataFileName, const char* ToyFileName, const char* outputFileName = "hist1.root")
{
  Int_t ptId = 0;
  Int_t maxLeadingPt = 4;
  Int_t maxAssocPt = 6;

  new TCanvas("c1", "", 800, 600);

  TFile* inputToy = TFile::Open(ToyFileName);
  TFile* inputData = TFile::Open(DataFileName);
  TFile* output = TFile::Open(outputFileName, "RECREATE");

  TH2* histToy = (TH2*) inputToy->Get(Form("projEta_6"));
  histToy->GetYaxis()->SetRangeUser(0.9994,1.006);
  histToy->Draw();
  histToy->Fit("pol2","I","SAME");
  TF1* func = new TF1("func","pol2",-1.8,1.8);
  func->SetParameter(0,1-(histToy->GetFunction("pol2")->GetParameter(0)-1));
  func->SetParameter(1,-histToy->GetFunction("pol2")->GetParameter(1));
  func->SetParameter(2,-histToy->GetFunction("pol2")->GetParameter(2));
  func->Draw("SAME");
  for (Int_t i=0; i<maxLeadingPt; i++)
  {
    for (Int_t j=1; j<maxAssocPt; j++)
    {
      // only process when first is filled
      TH2* histData = (TH2*) inputData->Get(Form("dphi_%d_%d_%d", i, j+1, 0));
      if (!histData)
        continue;
      if (histData->GetEntries() < 1e4)
      {
        Printf("%d %d Only %f entries. Skipping...", i, j, histData->GetEntries());
        continue;
      }
      ptId = i*(maxAssocPt-1)+j-1;
      if (SkipGraph(ptId))
        continue;
      for (Int_t histId = 0; histId < NHists; histId++)
      {
        histData = (TH2*) inputData->Get(Form("dphi_%d_%d_%d", i, j+1, histId));
        if (!histData)
          continue;

        if (histData->GetEntries() < 1e4)
        {
          Printf("%d %d %d Only %f entries. Skipping...", i, j, histId, histData->GetEntries());
          continue;
        }
        for (Int_t x=1; x<=histData->GetNbinsX();x++)
        {
          for (Int_t y=1; y<=histData->GetNbinsY();y++)
          {
            Float_t content = histData->GetBinContent(x,y);
            histData->SetBinContent(x,y,content*func->Eval(histData->GetYaxis()->GetBinCenter(y)));
          }
        }
        histData->Write(Form("dphi_%d_%d_%d", i, j+1, histId));
      }
    }
  }
  output->Close();
}


void CalculateIAA(const char* GraphFileName, const char* HistFileName=0, const char* outputFileName = "graphs.root")
{
  CreateGraphStructure();
  TFile* GraphFile = TFile::Open(GraphFileName);
  for (Int_t i=0; i<NGraphs; i++)
    for (Int_t j=0; j<NHists; j++)
    {
      if (i >= 32) continue;
      graphs[i][j] = (TGraphErrors*) GraphFile->Get(Form("graph_%d_%d", i, j));
    }

  Float_t IAA = 0;
  Float_t normP = 0;
  Float_t YerrorP = 0;
  Int_t ptId = 0;

  bool protonPoint = false;

  for (ptId = 0; ptId < NHists; ptId++)
  {
    if (SkipGraph(ptId))
      continue;
    for (Int_t i=0; i<graphs[0][ptId]->GetN(); i++)
      if (graphs[0][ptId]->GetX()[i] == 100)
      {
        normP = graphs[0][ptId]->GetY()[i];
        YerrorP = graphs[0][ptId]->GetErrorY(i);
        protonPoint = true;
      }
    for (Int_t i=0; i<graphs[0][ptId]->GetN() && protonPoint == true; i++)
    {
      IAA = (graphs[0][ptId]->GetY()[i])/normP;
      AddPoint(graphs[32][ptId], graphs[0][ptId]->GetX()[i], IAA, graphs[0][ptId]->GetErrorX(i), sqrt(graphs[0][ptId]->GetErrorY(i)/graphs[0][ptId]->GetY()[i]*graphs[0][ptId]->GetErrorY(i)/graphs[0][ptId]->GetY()[i]+YerrorP/normP*YerrorP/normP)*graphs[0][ptId]->GetY()[i]/normP);
    }
    protonPoint = false;
    graphs[32][ptId]->SetTitle(graphs[0][ptId]->GetTitle());
  }
  cout << "IAAFit graph created." << endl;


/*
  TFile::Open(HistFileName);
  Int_t maxLeadingPt = 4;
  Int_t maxAssocPt = 6;

  Float_t etaLimit = kEtaLimit;
  Float_t outerLimit = kOuterLimit;

  Double_t integral = 0;
  Double_t error = 0;

  Float_t centralityAxisMapping[] = { 5, 65, 100, 30, 15, 50 };
  Float_t centralityAxisMappingE[] = { 5, 5, 0, 10, 5, 10 };

  for (Int_t i=0; i<maxLeadingPt; i++)
  {
    for (Int_t j=1; j<maxAssocPt; j++)
    {
      // only process when first is filled
      TH2* hist = (TH2*) gFile->Get(Form("dphi_%d_%d_%d", i, j+1, 0));
      if (!hist)
        continue;
//      if (hist->GetEntries() < 10)
      if (hist->GetEntries() < 1e4)
      {
        Printf("%d %d Only %f entries. Skipping...", i, j, hist->GetEntries());
//        continue;
      }
      ptId = i*(maxAssocPt-1)+j-1;
      if (SkipGraph(ptId))
        continue;
      Int_t k=0;
      for (Int_t histId = 0; histId < NHists; histId++)
      {
        hist = (TH2*) gFile->Get(Form("dphi_%d_%d_%d", i, j+1, histId));
        if (!hist)
          continue;

        if (hist->GetEntries() < 1e4)
//        if (hist->GetEntries() < 10)
        {
          Printf("%d %d %d Only %f entries. Skipping...", i, j, histId, hist->GetEntries());
//          continue;
        }

        if (ptId >= 15)
        {
          etaLimit = 0.5;
          outerLimit = 0.99;
        }
        SubtractEtaGap(hist, etaLimit, outerLimit, kFALSE);
        Int_t bin1 = hist->FindBin(-etaLimit+0.01,-TMath::Pi() / 2 + 0.01);
        Int_t x1 = 0;
        Int_t y1 = 0;
        Int_t z1 = 0;
        hist->GetBinXYZ(bin1,x1,y1,z1);
        Int_t bin2 = hist->FindBin(etaLimit-0.01,TMath::Pi() / 2 - 0.01);
        Int_t x2 = 0;
        Int_t y2 = 0;
        Int_t z2 = 0;
        hist->GetBinXYZ(bin2,x2,y2,z2);
        integral = hist->IntegralAndError(x1,x2,y1,x2, error);

        AddPoint(graphs[33][ptId], centralityAxisMapping[histId], integral, centralityAxisMappingE[histId], error);
       k++;
        if (histId == 0)
          graphs[33][ptId]->SetTitle(hist->GetTitle());
      }
      for (k=0; k<graphs[33][ptId]->GetN(); k++)
        if (graphs[33][ptId]->GetX()[k] == 100) 
        {
          normP = graphs[33][ptId]->GetY()[k];
          YerrorP = graphs[33][ptId]->GetErrorY(k);
          protonPoint = true;
        }
      for (k=0; k<graphs[33][ptId]->GetN() && protonPoint == true; k++)
      {
        IAA = (graphs[33][ptId]->GetY()[k])/normP;
        AddPoint(graphs[34][ptId], graphs[33][ptId]->GetX()[k], IAA, graphs[33][ptId]->GetErrorX(k),sqrt(graphs[33][ptId]->GetErrorY(k)/graphs[33][ptId]->GetY()[k]*graphs[33][ptId]->GetErrorY(k)/graphs[33][ptId]->GetY()[k]+YerrorP/normP*YerrorP/normP)*graphs[33][ptId]->GetY()[k]/normP);
      }
      protonPoint = false;
      graphs[34][ptId]->SetTitle(graphs[33][ptId]->GetTitle());
    }
  }
  cout << "IAAHist graph created." << endl;
*/
  WriteGraphs(outputFileName);
}

void AcceptanceEfficiencyToy()
{
  // toy MC to study the intermix of a detector efficiency as fct of eta and the mixed event correction
  
//     TF1* eff = new TF1("eff", "0.9 - 0.3 * abs(x)", -1, 1);

  // pA efficiency
  Double_t xAxis1[51] = {-2.5, -2.4, -2.3, -2.2, -2.1, -2, -1.9, -1.8, -1.7, -1.6, -1.5, -1.4, -1.3, -1.2, -1.1, -1, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.1, 2.2, 2.3, 2.4, 2.5}; 
  TH1D *effHist = new TH1D("effHist","step1: projection on #Delta#eta",50, xAxis1);
  effHist->SetBinContent(13,0.271338);
  effHist->SetBinContent(14,0.334383);
  effHist->SetBinContent(15,0.4399013);
  effHist->SetBinContent(16,0.6723888);
  effHist->SetBinContent(17,0.8022655);
  effHist->SetBinContent(18,0.809605);
  effHist->SetBinContent(19,0.8145845);
  effHist->SetBinContent(20,0.817339);
  effHist->SetBinContent(21,0.8210758);
  effHist->SetBinContent(22,0.8228529);
  effHist->SetBinContent(23,0.8246998);
  effHist->SetBinContent(24,0.8249083);
  effHist->SetBinContent(25,0.8162804);
  effHist->SetBinContent(26,0.8148814);
  effHist->SetBinContent(27,0.8237724);
  effHist->SetBinContent(28,0.8247403);
  effHist->SetBinContent(29,0.822931);
  effHist->SetBinContent(30,0.8220786);
  effHist->SetBinContent(31,0.8195946);
  effHist->SetBinContent(32,0.8159786);
  effHist->SetBinContent(33,0.812606);
  effHist->SetBinContent(34,0.8080089);
  effHist->SetBinContent(35,0.6991395);
  effHist->SetBinContent(36,0.4702392);
  effHist->SetBinContent(37,0.3467974);
  effHist->SetBinContent(38,0.2798031);
  effHist->SetBinError(13,0.000586557);
  effHist->SetBinError(14,0.0006215941);
  effHist->SetBinError(15,0.0006542189);
  effHist->SetBinError(16,0.0006190971);
  effHist->SetBinError(17,0.0005263518);
  effHist->SetBinError(18,0.0005197939);
  effHist->SetBinError(19,0.0005148436);
  effHist->SetBinError(20,0.0005131777);
  effHist->SetBinError(21,0.0005095525);
  effHist->SetBinError(22,0.000508042);
  effHist->SetBinError(23,0.0005058743);
  effHist->SetBinError(24,0.0005058348);
  effHist->SetBinError(25,0.000514767);
  effHist->SetBinError(26,0.0005161097);
  effHist->SetBinError(27,0.0005054953);
  effHist->SetBinError(28,0.0005028134);
  effHist->SetBinError(29,0.0005030626);
  effHist->SetBinError(30,0.0005016751);
  effHist->SetBinError(31,0.0005036153);
  effHist->SetBinError(32,0.0005049182);
  effHist->SetBinError(33,0.0005063787);
  effHist->SetBinError(34,0.0005090265);
  effHist->SetBinError(35,0.0005904377);
  effHist->SetBinError(36,0.0006403924);
  effHist->SetBinError(37,0.0006090749);
  effHist->SetBinError(38,0.0005733178);

  //   eff = new TF1("eff", "0.9 - 0.05 * abs(x) - 0.6 * TMath::Floor(abs(x) + 0.2)", -1, 1);
//   eff = new TF1("eff", "0.9 - 0.6 * TMath::Floor(abs(x) + 0.2)", -1, 1);
//   eff = new TF1("eff", "0.9", -1, 1);

//   for (Int_t i=1; i<=50; i++)
//     effHist->SetBinContent(i, 0.2);

  Float_t etaAcceptance = 1.2;

  new TCanvas; effHist->Draw(); effHist->Fit("pol0", "", "", -etaAcceptance+0.01, etaAcceptance-0.01);
  
  Int_t bins = 200;
  TH1D* same = new TH1D("same", "", bins, -2.5, 2.5);
  TH1D* mixed = new TH1D("mixed", "", bins, -2.5, 2.5);
  TH1D* eta =  new TH1D("eta", "", bins, -2.5, 2.5);
  
  TH1D* sameEff = new TH1D("sameEff", "", bins, -2.5, 2.5);
  TH1D* mixedEff = new TH1D("mixedEff", "", bins, -2.5, 2.5);
  TH1D* etaEff =  new TH1D("etaEff", "", bins, -2.5, 2.5);
  TH1D* allEtaEff =  new TH1D("allEtaEff", "", bins, -2.5, 2.5);

  TH1D* etaSource =  new TH1D("etaSource", "", bins, -2.5, 2.5);
  TH1D* etaSourceEff =  new TH1D("etaSourceEff", "", bins, -2.5, 2.5);

  Float_t sigma = 0.5;
  Float_t assoc = 0;
  Bool_t assocTracked = kFALSE;
  
  same->Sumw2();
  mixed->Sumw2();
  eta->Sumw2();
  sameEff->Sumw2();
  mixedEff->Sumw2();
  etaEff->Sumw2();
  allEtaEff->Sumw2();
  etaSource->Sumw2();
  etaSourceEff->Sumw2();
  
  for (Int_t i=0; i<1000000; i++)
  {
    // randomize mean
    Float_t mean = gRandom->Uniform(-5, 5);
    
    Float_t trig = gRandom->Gaus(mean, sigma);
    Bool_t trigTracked = (gRandom->Uniform() < effHist->GetBinContent(effHist->FindBin(trig))); // eff->Eval(trig)
    if (TMath::Abs(trig) < etaAcceptance)
    {
      eta->Fill(trig);
      if (trigTracked)
	etaEff->Fill(trig);
    }
    
    // mixed event
    if (i > 0 && TMath::Abs(trig) < etaAcceptance && TMath::Abs(assoc) < etaAcceptance)
    {
      mixed->Fill(trig - assoc);
    
      if (trigTracked && assocTracked)
	mixedEff->Fill(trig - assoc);
// 	mixedEff->Fill(trig - assoc, 1.0 / effHist->GetBinContent(effHist->FindBin(assoc)));
    }
    
//     mean = gRandom->Uniform(-5, 5);
    assoc = gRandom->Gaus(mean, sigma);
    
    assocTracked = (gRandom->Uniform() < effHist->GetBinContent(effHist->FindBin(assoc)));
    
    if (TMath::Abs(trig) < etaAcceptance && trigTracked)
      allEtaEff->Fill(trig);
    if (TMath::Abs(assoc) < etaAcceptance && assocTracked)
      allEtaEff->Fill(assoc);
    
    // same event
    if (TMath::Abs(trig) < etaAcceptance && TMath::Abs(assoc) < etaAcceptance)
    {
      same->Fill(trig - assoc);
      
//       if (trigTracked)
	etaSource->Fill(assoc);
      
      if (trigTracked && assocTracked)
      {
	sameEff->Fill(trig - assoc);
	etaSourceEff->Fill(assoc);
/*	sameEff->Fill(trig - assoc, 1.0 / effHist->GetBinContent(effHist->FindBin(assoc)));
	etaSourceEff->Fill(assoc, 1.0 / effHist->GetBinContent(effHist->FindBin(assoc)));*/
      }
    }
  }
  
  Float_t mixedConstant = 0;
  if (1)
  {
    new TCanvas;
    TF1* pol1 = new TF1("pol", "pol1(0)", -10, 10);
    mixed->Fit(pol1, "+", "", -1, -0.0001);
    Float_t mixedConstant1 = pol1->Eval(0);
    pol1 = new TF1("pol", "pol1(0)", -10, 10);
    mixed->Fit(pol1, "+", "", 0.0001, 1);
    Float_t mixedConstant2 = pol1->Eval(0);
    mixedConstant = (mixedConstant1 + mixedConstant2) / 2;
    Printf("%f %f %f", mixedConstant1, mixedConstant2, mixedConstant);
  }
  else
    mixedConstant = mixed->Integral(bins / 2, bins / 2 + 1) / 2;

  mixed = (TH1D*) mixed->Clone();
  mixed->Scale(1.0 / mixedConstant);
  
  new TCanvas; same->DrawCopy(); mixed->SetLineColor(2); mixed->DrawCopy("SAME");
  
  same->Divide(mixed);
  same->Scale(1.0 / eta->Integral());
  new TCanvas; same->DrawCopy();

  new TCanvas; eta->DrawCopy(); etaEff->SetLineColor(2); etaEff->DrawCopy("SAME");
  
  new TCanvas; etaSource->DrawCopy(); etaSourceEff->SetLineColor(2); etaSourceEff->DrawCopy("SAME"); allEtaEff->SetLineColor(4); allEtaEff->DrawCopy("SAME");
  
//   etaSource->Multiply(effHist);
  
  if (0)
  {
    new TCanvas;
    TF1* pol1 = new TF1("pol", "pol1(0)", -10, 10);
    mixedEff->Fit(pol1, "+", "", -1, -0.0001);
    Float_t mixedConstant1 = pol1->Eval(0);
    pol1 = new TF1("pol", "pol1(0)", -10, 10);
    mixedEff->Fit(pol1, "+", "", 0.0001, 1);
    Float_t mixedConstant2 = pol1->Eval(0);
    mixedConstant = (mixedConstant1 + mixedConstant2) / 2;
    Printf("%f %f %f", mixedConstant1, mixedConstant2, mixedConstant);
  }
  else
    mixedConstant = mixedEff->Integral(bins / 2, bins / 2 + 1) / 2;
  
  mixedEff = (TH1D*) mixedEff->Clone();
  mixedEff->Scale(1.0 / mixedConstant);
  
  new TCanvas; sameEff->DrawCopy(); mixedEff->SetLineColor(2); mixedEff->DrawCopy("SAME");

  sameEff->Divide(mixedEff);
  sameEff->Scale(1.0 / etaEff->Integral());
  new TCanvas; sameEff->DrawCopy();
  
  new TCanvas; same->DrawCopy(); sameEff->SetLineColor(2); sameEff->DrawCopy("SAME");
  
  sameEff->Divide(same);
  new TCanvas; sameEff->DrawCopy(); sameEff->Fit("pol0");
}

void Acceptance2DToy(Float_t etaAcceptance = 1.0)
{
  // toy MC to study the effect of acceptance on the correlation function

  Int_t bins = 40;
  TH2D* same = new TH2D("same", "", bins, -2.5, 2.5, bins, -TMath::Pi(), TMath::Pi());
  TH2D* mixed = new TH2D("mixed", "", bins, -2.5, 2.5, bins,  -TMath::Pi(), TMath::Pi());
  TH1D* eta =  new TH1D("eta", "", bins, -2.5, 2.5);
  TH1D* phi =  new TH1D("phi", "", bins, -TMath::Pi(), TMath::Pi());
  
  Float_t sigma = 0.4;
  const Int_t nParticles = 8;
  Float_t etas[nParticles];
  Float_t phis[nParticles];
  Float_t lastetas[nParticles];
  Float_t lastphis[nParticles];
  
  TF2* func = new TF2("func", "[0]*exp(-0.5*((x/[1])**2+(y/[2])**2))", -5, 5, -5, 5);
  func->SetParameters(1, sigma, sigma);

  same->Sumw2();
  mixed->Sumw2();
  eta->Sumw2();
  
  for (Int_t i=0; i<1000000; i++)
  {
    for (Int_t j=0; j<nParticles/2; j++)
    {
      // randomize mean
      Float_t meanEta = gRandom->Uniform(-5, 5);
      Float_t meanPhi = gRandom->Uniform(0, TMath::TwoPi());
      
      Double_t gausEta, gausPhi;
      func->GetRandom2(gausEta, gausPhi);
      
      Float_t trigEta = meanEta + gausEta;
      Float_t trigPhi = meanPhi + gausPhi;
      
      if (trigPhi > TMath::Pi())
	trigPhi -= TMath::TwoPi();
      if (trigPhi < -TMath::Pi())
	trigPhi += TMath::TwoPi();

      func->GetRandom2(gausEta, gausPhi);
      
      Float_t assocEta = meanEta + gausEta;
      Float_t assocPhi = meanPhi + gausPhi;
      
      if (assocPhi > TMath::Pi())
	assocPhi -= TMath::TwoPi();
      if (assocPhi < -TMath::Pi())
	assocPhi += TMath::TwoPi();
      
      etas[j*2] = trigEta;
      etas[j*2+1] = assocEta;
      phis[j*2] = trigPhi;
      phis[j*2+1] = assocPhi;
    }
    
    // same event
    for (Int_t j=0; j<nParticles; j++)
    {
      if (TMath::Abs(etas[j]) > etaAcceptance)
	continue;

      eta->Fill(etas[j]);
      phi->Fill(phis[j]);

      for (Int_t k=j+1; k<nParticles; k++)
      {
	if (TMath::Abs(etas[k]) > etaAcceptance)
	  continue;
	
	Float_t deltaEta = etas[j] - etas[k];
	Float_t deltaPhi = phis[j] - phis[k];
      
	if (deltaPhi > TMath::Pi())
	  deltaPhi -= TMath::TwoPi();
	if (deltaPhi < -TMath::Pi())
	  deltaPhi += TMath::TwoPi();
    
	same->Fill(deltaEta, deltaPhi);
      }
    }
    
    // mixed event
    if (i > 0)
    {
      for (Int_t j=0; j<nParticles; j++)
      {
	if (TMath::Abs(etas[j]) > etaAcceptance)
	  continue;

	for (Int_t k=0; k<nParticles; k++)
	{
	  if (TMath::Abs(lastetas[k]) > etaAcceptance)
	    continue;
	  
	  Float_t deltaEta = etas[j] - lastetas[k];
	  Float_t deltaPhi = phis[j] - lastphis[k];
	
	  if (deltaPhi > TMath::Pi())
	    deltaPhi -= TMath::TwoPi();
	  if (deltaPhi < -TMath::Pi())
	    deltaPhi += TMath::TwoPi();
      
	  mixed->Fill(deltaEta, deltaPhi);
	}
      }
      
      for (Int_t j=0; j<nParticles; j++)
      {
	lastetas[j] = etas[j];
	lastphis[j] = phis[j];
      }
    }
/*    // add an uncorrelated particle
    Float_t assocEta2 = gRandom->Uniform(-5, 5);
    Float_t assocPhi2 = gRandom->Uniform(0, TMath::TwoPi());

    deltaEta = trigEta - assocEta2;
    deltaPhi = trigPhi - assocPhi2;

    if (deltaPhi > TMath::Pi())
      deltaPhi -= TMath::TwoPi();
    if (deltaPhi < -TMath::Pi())
      deltaPhi += TMath::TwoPi();

    // same event
    if (TMath::Abs(trigEta) < etaAcceptance && TMath::Abs(assocEta) < etaAcceptance)
      same->Fill(deltaEta, deltaPhi);*/
  }
  
  Float_t mixedConstant = 0;
  if (1)
  {
    TH1* mixedProj = mixed->ProjectionX();
    mixedProj->Scale(1.0 / mixed->GetNbinsY());
    new TCanvas;
    TF1* pol1 = new TF1("pol", "pol1(0)", -10, 10);
    mixedProj->Fit(pol1, "+", "", -1, -0.0001);
    Float_t mixedConstant1 = pol1->Eval(0);
    pol1 = new TF1("pol", "pol1(0)", -10, 10);
    mixedProj->Fit(pol1, "+", "", 0.0001, 1);
    Float_t mixedConstant2 = pol1->Eval(0);
    mixedConstant = (mixedConstant1 + mixedConstant2) / 2;
    Printf("%f %f %f", mixedConstant1, mixedConstant2, mixedConstant);
  }
  else
    mixedConstant = mixed->Integral(bins / 2, bins / 2 + 1) / 2;

  mixed = (TH2D*) mixed->Clone();
  mixed->Scale(1.0 / mixedConstant);
  
  new TCanvas; same->DrawCopy("SURF1"); 
  new TCanvas; mixed->DrawCopy("SURF1");
  
  same->Divide(mixed);
  same->Scale(1.0 / eta->Integral());
  new TCanvas; same->DrawCopy("SURF1");

  new TCanvas; eta->DrawCopy();
  new TCanvas; phi->DrawCopy();
  
  Printf("%f", same->Integral());
  
  TH1* projSame = (TH1*) same->ProjectionX("projSame");

  TF1* func1 = new TF1("func1", "[0]+gaus(1)", -5, 5);
  func1->SetParameters(0.01, 0.01, 0, sigma*2);
  func1->FixParameter(2, 0);
  
  projSame->Fit(func1);
  func1->SetParameter(0, 0);
  Printf("%f", func1->Integral(-2, 2) / projSame->GetBinWidth(1));
}

void Convolute()
{
  TFile::Open("phi.root");
  TH1* phi = (TH1*) gFile->Get("phi");

  phi->SetFillColor(0);
  //TH1* phi = new TH1F("phi", "", 20, 0, 2); phi->SetBinContent(9, 1);  phi->SetBinContent(10, 1); phi->SetBinContent(11, 1); phi->SetBinContent(12, 1);
  
//   TH1* phi = new TH1F("phi", "", 20, 0, TMath::TwoPi()); 
//   for (Int_t x=0; x<phi->GetNbinsX(); x++)    phi->SetBinContent(x+1, 1);
//   phi->SetBinContent(10, 0);
//   phi->SetBinContent(11, 0);
  
//   
//   TH1* phi = new TH1F("phi", "", 4, -1, 1); phi->SetBinContent(1, 1);  phi->SetBinContent(2, 1); phi->SetBinContent(3, 1); phi->SetBinContent(4, 1);
  
  new TCanvas; phi->Draw();
  
  TH1* conv = (TH1*) phi->Clone("conv");
  conv->Reset();
  
  for (Int_t delta=0; delta<phi->GetNbinsX(); delta++)
  {
    Double_t value = 0;
    for (Int_t x=1; x<=phi->GetNbinsX(); x++)
    {
      Int_t y = x + delta;
      if (y < 1)
	y += phi->GetNbinsX();
      if (y > phi->GetNbinsX())
	y -= phi->GetNbinsX();
      value += phi->GetBinContent(x) * phi->GetBinContent(y);
    }
    conv->SetBinContent(delta+1, value);
  }
  
  conv->Scale(1.0 / conv->Integral());
  new TCanvas; conv->Draw();
  
//   return;
  
  // random approach
  
  TH1* rand = (TH1*) phi->Clone("rand");
  rand->Reset();
  
  for (Int_t i=0; i<1000000000; i++)
  {
    Float_t value1 = phi->GetRandom();
    Float_t value2 = phi->GetRandom();
    
    Int_t x = rand->FindBin(value1);
    Int_t y = rand->FindBin(value2);
    
    Int_t bin = x - y;
    if (bin < 0)
      bin += phi->GetNbinsX();
    if (bin >= phi->GetNbinsX())
      bin -= phi->GetNbinsX();
    
    rand->SetBinContent(bin+1, rand->GetBinContent(bin+1) + 1);
  }
  
  rand->Scale(1.0 / rand->Integral());
  rand->SetLineColor(2);
  rand->Draw("SAME");
}

void TwoPlusOneCorrelations_Draw(TH2* hist)
{
  // get rid of triangle
  for (Int_t i=1; i<=hist->GetNbinsX(); i++)
  {
    Float_t weight = 1.0 / (1.0 - TMath::Abs(hist->GetXaxis()->GetBinCenter(i)) / 2.0);
    for (Int_t j=1; j<=hist->GetNbinsY(); j++)
      hist->SetBinContent(i, j, hist->GetBinContent(i, j) * weight);
  }

  TCanvas* c = new TCanvas(Form("%s_c", hist->GetName()), Form("%s_c", hist->GetName()), 600, 1000);
  c->Divide(1, 3);
  
  c->cd(1);
  hist->Draw("SURF1");
  hist->GetXaxis()->SetRangeUser(-1.59, 1.59);
  
  Int_t eta1 = hist->GetXaxis()->FindBin(-1.59);
  Int_t eta2 = hist->GetXaxis()->FindBin(-1.01);
  Int_t eta3 = eta2+1;
  Int_t eta4 = hist->GetXaxis()->FindBin(0.99);
  Int_t eta5 = eta4+1;
  Int_t eta6 = hist->GetXaxis()->FindBin(1.59);
  
  TH1* nsProj1 = hist->ProjectionY(Form("%s_Proj1", hist->GetName()), eta1, eta2);
  TH1* nsProj2 = hist->ProjectionY(Form("%s_Proj2", hist->GetName()), eta3, eta4);
  TH1* nsProj3 = hist->ProjectionY(Form("%s_Proj3", hist->GetName()), eta5, eta6);
  
  nsProj1->Add(nsProj3);
  nsProj1->Scale(1.0 / (eta2 - eta1 + 1 + eta6 - eta5 + 1));
  nsProj2->Scale(1.0 / (eta4 - eta3 + 1));
  
  c->cd(2);
  nsProj2->DrawCopy(); 
  nsProj1->SetLineColor(2); 
  nsProj1->Draw("SAME");
  
  c->cd(3);
  nsProj2->Add(nsProj1, -1);
  nsProj2->DrawCopy(); 
}

void TwoPlusOneCorrelations()
{
  // Check 2 plus 1 correlations for NS and AS and eta gap flow subtraction

  Int_t bins = 40;
  TH2D* ns =  new TH2D("ns", ";delta eta; delta phi", bins, -2, 2, bins, -0.5 * TMath::Pi(), 1.5 * TMath::Pi());
  TH2D* ns2 = new TH2D("ns2", ";delta eta; delta phi", bins, -2, 2, bins, -0.5 * TMath::Pi(), 1.5 * TMath::Pi());
  TH2D* as2 = new TH2D("as2", ";delta eta; delta phi", bins, -2, 2, bins, -0.5 * TMath::Pi(), 1.5 * TMath::Pi());
  TH2D* ns2Same = new TH2D("ns2Same", ";delta eta; delta phi", bins, -2, 2, bins, -0.5 * TMath::Pi(), 1.5 * TMath::Pi());
  TH2D* as2Same = new TH2D("as2Same", ";delta eta; delta phi", bins, -2, 2, bins, -0.5 * TMath::Pi(), 1.5 * TMath::Pi());
//   TH2D* as2AllSame = new TH2D("as2AllSame", ";delta eta; delta phi", bins, -2, 2, bins, -0.5 * TMath::Pi(), 1.5 * TMath::Pi());
  TH2D* ns2Diff = new TH2D("ns2Diff", ";delta eta; delta phi", bins, -2, 2, bins, -0.5 * TMath::Pi(), 1.5 * TMath::Pi());
  TH2D* as2Diff = new TH2D("as2Diff", ";delta eta; delta phi", bins, -2, 2, bins, -0.5 * TMath::Pi(), 1.5 * TMath::Pi());

  Float_t v2 = 0.2;
  
  const Double_t kPi = TMath::Pi();
  const Double_t kTwoPi = 2 * kPi;
  
  // event loop
  for (Int_t i=0; i<500; i++)
  {
    if (i % 10 == 0)
      Printf("%d", i);
    Float_t rpangle = gRandom->Uniform(-kPi, kPi);
    
    const Int_t nParticles = 500;
    Float_t etas[nParticles];
    Float_t phis[nParticles];
    Float_t weights[nParticles];
    Int_t   ids[nParticles];
    Int_t nP = 0;
    
    // some particles from the continuum
    if (0)
    {
      for (Int_t j=0; j<nParticles/2; j++)
      {
	etas[nP] = gRandom->Uniform(-1, 1);
	phis[nP] = gRandom->Uniform(-kPi, kPi);
	ids[nP] = -1;
	if (TMath::Abs(etas[nP]) < 1.0)
	  nP++;
      }
    }
    
    // jets
    if (1)
    {
      Float_t width = 0.15;
      Float_t width2 = 0.3;
      for (Int_t j=0; j<20; j++)
      {
	Float_t jetEta = gRandom->Uniform(-1.5, 1.5);
	Float_t jetPhi = gRandom->Uniform(-kPi, kPi);
	
	Float_t jetASEta = gRandom->Uniform(-1.5, 1.5);
	Float_t jetASPhi = jetPhi + kPi;
	
	// kt
// 	jetASPhi += gRandom->Uniform(-0.2, 0.2);
	
	for (Int_t k=0; k<5; k++)
	{
	  etas[nP] = gRandom->Gaus(jetEta, width);
	  phis[nP] = gRandom->Gaus(jetPhi, width);
	  ids[nP] = j;
	  if (TMath::Abs(etas[nP]) < 1.0)
	    nP++;
	}
	  
	for (Int_t k=0; k<5; k++)
	{
	  etas[nP] = gRandom->Gaus(jetASEta, width2);
	  phis[nP] = gRandom->Gaus(jetASPhi, width2);
	  ids[nP] = j;
	  if (TMath::Abs(etas[nP]) < 1.0)
	    nP++;
	}
      }
    }
    
    // calculate flow weights
    for (Int_t j=0; j<nP; j++)
      weights[j] = 1.0 + v2 * TMath::Cos(2*(rpangle - phis[j]));
    
    // fill, two particle correlation
    for (Int_t j=0; j<nP; j++)
    {
      for (Int_t k=0; k<nP; k++)
      {
	if (j == k)
	  continue;
	Float_t deltaPhi = phis[j] - phis[k];
	while (deltaPhi > 1.5 * kPi)
	  deltaPhi -= kTwoPi;
	while (deltaPhi < -0.5 * kPi)
	  deltaPhi += kTwoPi;
	
	Float_t deltaEta = etas[j] - etas[k];

	Float_t weight = 1;
	// apply flow
	weight *= weights[j] * weights[k];
	
	ns->Fill(deltaEta, deltaPhi, weight);
      }
    }
    
    // fill, two plus one particle correlation
    Float_t alpha = 0.2;
    for (Int_t j=0; j<nP; j++)
    {
      for (Int_t k=0; k<nP; k++)
      {
	if (j == k)
	  continue;
	Float_t deltaPhi = phis[j] - phis[k];
	while (deltaPhi > 1.5 * kPi)
	  deltaPhi -= kTwoPi;
	while (deltaPhi < -0.5 * kPi)
	  deltaPhi += kTwoPi;

	// back to back trigger
	if (TMath::Abs(deltaPhi - kPi) < alpha)
	{
	  for (Int_t l=0; l<nP; l++)
	  {
	    if (j == l || k == l)
	      continue;
	    
	    // NS
	    deltaPhi = phis[j] - phis[l];
	    while (deltaPhi > 1.5 * kPi)
	      deltaPhi -= kTwoPi;
	    while (deltaPhi < -0.5 * kPi)
	      deltaPhi += kTwoPi;
	    
	    Float_t deltaEta = etas[j] - etas[l];

	    Float_t weight = 1;
	    // apply flow
// 	    weight *= weights[j] * weights[k] * weights[l];
	    
	    ns2->Fill(deltaEta, deltaPhi, weight);
	    if (ids[j] == ids[l] && ids[j] == ids[k])
	      ns2Same->Fill(deltaEta, deltaPhi, weight);
	    else
	      ns2Diff->Fill(deltaEta, deltaPhi, weight);
	    
	    // AS
	    deltaPhi = phis[k] - phis[l];
	    while (deltaPhi > 1.5 * kPi)
	      deltaPhi -= kTwoPi;
	    while (deltaPhi < -0.5 * kPi)
	      deltaPhi += kTwoPi;
	    
	    deltaEta = etas[k] - etas[l];

	    weight = 1;
	    // apply flow
// 	    weight *= weights[j] * weights[k] * weights[l];
	    
	    as2->Fill(deltaEta, deltaPhi, weight);
	    if (ids[j] == ids[l] && ids[j] == ids[k])
	      as2Same->Fill(deltaEta, deltaPhi, weight);
	    else
	      as2Diff->Fill(deltaEta, deltaPhi, weight);

/*	    if (ids[k] == ids[l] && ids[j] == ids[l])
	      as2AllSame->Fill(deltaEta, deltaPhi, weight);*/
	  }
	}
      }
    }
  }
  
  TwoPlusOneCorrelations_Draw(ns);
  TwoPlusOneCorrelations_Draw(ns2);
  TwoPlusOneCorrelations_Draw(as2);
  TwoPlusOneCorrelations_Draw(ns2Same);
  TwoPlusOneCorrelations_Draw(as2Same);
//   TwoPlusOneCorrelations_Draw(as2AllSame);
  TwoPlusOneCorrelations_Draw(ns2Diff);
  TwoPlusOneCorrelations_Draw(as2Diff);
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

void GetProjections(TH2* hist, TH1** projPhi, TH1** projEta, Int_t k)
{
  Float_t projectLimit = 0.8;

  *projPhi = hist->ProjectionX(Form("%s_proj1%d", hist->GetName(), k), hist->GetYaxis()->FindBin(-projectLimit+0.01), hist->GetYaxis()->FindBin(projectLimit-0.01));
    
  *projEta = hist->ProjectionY(Form("%s_proj2_%d", hist->GetName(), k), hist->GetXaxis()->FindBin(-projectLimit+0.01), hist->GetXaxis()->FindBin(projectLimit-0.01));
}



/*Drawing functions*/
void CompareCMSResults(const char * fileName, const char * graphFile)
{
  ifstream inputFile;
  inputFile.open(fileName);
  string line;
  getline(inputFile,line);
  double centrLow, centrHigh, pT, deta, dphi;
  vector<TGraphErrors*> graphEtaCMS; // = new TGraph;
  vector<TGraphErrors*> graphPhiCMS; // = new TGraph;
  vector<double> pTV;
  while (inputFile >> centrLow >> centrHigh >> pT >> deta >> dphi)
  {
    int index = -1;
    for (unsigned int i=0; i<pTV.size(); i++)
      if (pTV[i] == pT) 
      {
        index = i;
        break;
      }
    if (index == -1)
    {
      pTV.push_back(pT);
      graphEtaCMS.push_back(new TGraphErrors);
      graphPhiCMS.push_back(new TGraphErrors);
      graphEtaCMS[graphEtaCMS.size()-1]->SetPoint(graphEtaCMS[graphEtaCMS.size()-1]->GetN(),centrLow+(centrHigh-centrLow)/2.,deta);
      graphEtaCMS[graphEtaCMS.size()-1]->SetPointError(graphEtaCMS[graphEtaCMS.size()-1]->GetN()-1,(centrHigh-centrLow)/2.,0);
      graphPhiCMS[graphPhiCMS.size()-1]->SetPoint(graphPhiCMS[graphPhiCMS.size()-1]->GetN(),centrLow+(centrHigh-centrLow)/2.,dphi);
      graphPhiCMS[graphPhiCMS.size()-1]->SetPointError(graphPhiCMS[graphPhiCMS.size()-1]->GetN()-1,(centrHigh-centrLow)/2.,0);
    }
    else
    {
      graphEtaCMS[index]->SetPoint(graphEtaCMS[index]->GetN(),centrLow+(centrHigh-centrLow)/2.,deta);
      graphEtaCMS[index]->SetPointError(graphEtaCMS[index]->GetN()-1,(centrHigh-centrLow)/2.,0);
      graphPhiCMS[index]->SetPoint(graphPhiCMS[index]->GetN(),centrLow+(centrHigh-centrLow)/2.,dphi);
      graphPhiCMS[index]->SetPointError(graphPhiCMS[index]->GetN()-1,(centrHigh-centrLow)/2.,0);
    }
    cerr << centrLow << "\t" << centrHigh << "\t" << pT << "\t" << deta << "\t" << dphi << endl;
  }
  ReadGraphs(graphFile);
  Int_t nHists = 24; //NHists;

  TLegend * legend = new TLegend(0.5,0.7,0.9,0.9);
  legend->SetFillColor(0);

//  TCanvas* etaC = new TCanvas("etaC","etaC",800,600);
  DrawCentrality("eta_rms_centrality", nHists, graphs[11+16], 0, 1.0, Form("#sigma_{#Delta#eta} (%s)", fitLabel));
  for (unsigned int i=0; i<graphEtaCMS.size(); i++)
  {
    graphEtaCMS[i]->SetMarkerStyle(20+i);
    graphEtaCMS[i]->SetMarkerColor(4);
    graphEtaCMS[i]->SetLineColor(4);
    graphEtaCMS[i]->GetYaxis()->SetRangeUser(0,1);
    graphEtaCMS[i]->DrawClone("SAMEP");
    graphEtaCMS[i]->SetFillColor(0);
    legend->AddEntry(graphEtaCMS[i]->Clone(), Form("CMS, p_{T,a} = %.1f",pTV[i]));
  }
  legend->Draw();
//  TCanvas* phiC = new TCanvas("phiC","phiC",800,600);
  DrawCentrality("phi_rms_centrality", nHists, graphs[10+16], 0, 0.7, Form("#sigma_{#Delta#varphi} (%s) (rad.)", fitLabel));
  for (unsigned int i=0; i<graphPhiCMS.size(); i++)
  {
    graphPhiCMS[i]->SetMarkerStyle(20+i);
    graphPhiCMS[i]->SetMarkerColor(4);
    graphPhiCMS[i]->SetLineColor(4);
    graphPhiCMS[i]->GetYaxis()->SetRangeUser(0,1);
    graphPhiCMS[i]->SetFillColor(0);
    graphPhiCMS[i]->DrawClone("SAMEP");
  }
  legend->Draw();

  TLegend * legend2 = new TLegend(0.5,0.7,0.9,0.9);
  legend2->SetFillColor(0);
  new TCanvas("etaC","etaC",800,600);
  for (unsigned int i=0; i<graphEtaCMS.size(); i++)
  {
    graphEtaCMS[i]->SetMarkerStyle(20+i);
    graphEtaCMS[i]->SetMarkerColor(i+1);
    graphEtaCMS[i]->SetLineColor(i+1);
    graphEtaCMS[i]->GetYaxis()->SetRangeUser(0,1);
    graphEtaCMS[i]->Draw(i==0?"AP":"SAMEP");
    graphEtaCMS[i]->SetFillColor(0);
    legend2->AddEntry(graphEtaCMS[i]->Clone(), Form("CMS, p_{T,a} = %.1f",pTV[i]));
  }
  legend2->Draw();
  new TCanvas("phiC","phiC",800,600);
  for (unsigned int i=0; i<graphPhiCMS.size(); i++)
  {
    graphPhiCMS[i]->SetMarkerStyle(20+i);
    graphPhiCMS[i]->SetMarkerColor(i+1);
    graphPhiCMS[i]->SetLineColor(i+1);
    graphPhiCMS[i]->GetYaxis()->SetRangeUser(0,1);
    graphPhiCMS[i]->SetFillColor(0);
    graphPhiCMS[i]->Draw(i==0?"AP":"SAMEP");
  }
  legend2->Draw();
}

void CompareCP(string dataFile, string files, string labels="", double centrality=5.)
{
  vector<string> filesV;
  filesV.push_back(dataFile);
  std::istringstream filesIs(files);
  string filesStr;
  while( filesIs >> filesStr)
    filesV.push_back(filesStr);

  vector<string> labelsV;
  std::istringstream labelsIs(labels);
  string labelsStr;
  while( labelsIs >> labelsStr)
    labelsV.push_back(labelsStr);

  if (filesV.size() != labelsV.size())
  {
    cerr << "Different number of labels and files" << endl;
    if (labelsV.size() != 0) return;
  }
  Int_t nHists = 24; //NHists;

  Int_t colors[16]  = { 1, 3, 2, 6, 4, 7, 8, 9, 11, 12, 28, 30, 36, 40, 46 };
  Int_t markers[16] = { 20, 21, 22, 23, 24, 25, 26, 27, 28, 30, 31, 32, 33, 34, 2, 5};
  vector<Int_t> graphTypes;
  vector<string> yTitles;
  vector<string> yTitlesRatio;
  graphTypes.push_back(26);
  yTitles.push_back("#sigma_{#Delta#varphi} (0-10%)/#sigma_{#Delta#varphi} (50-80%)");
  yTitlesRatio.push_back("#sigma_{#Delta#varphi} (Data)/#sigma_{#Delta#varphi} (MC)");
  graphTypes.push_back(27);
  yTitles.push_back("#sigma_{#Delta#eta} (0-10%)/#sigma_{#Delta#eta} (50-80%)");
  yTitlesRatio.push_back("#sigma_{#Delta#eta} (Data)/#sigma_{#Delta#eta} (MC)");
  if (graphTypes.size() != yTitles.size())
  {
    cerr << "Need one y title for each graph" << endl;
    return;
  }

  vector<TCanvas*> graphC(graphTypes.size());
  vector<TCanvas*> graphFitC(graphTypes.size());
  vector<TCanvas*> ratioC(graphTypes.size());
  double yMax = 2.9;
  double yMin = 0.85;
  vector<TString> labelV;
//  vector<vector<TCanvas*> > graphCentralityC(graphTypes.size());
  for (unsigned int iGraph=0; iGraph<graphTypes.size(); iGraph++)
  {
    graphC[iGraph] = new TCanvas(Form("graph_%d",graphTypes[iGraph]),"graphC",600,600);
    gPad->SetRightMargin(0.02);
    gPad->SetTopMargin(0.02);
    graphFitC[iGraph] = new TCanvas(Form("graphFit_%d",graphTypes[iGraph]),"graphFitC",600,600);
    gPad->SetRightMargin(0.02);
    gPad->SetTopMargin(0.02);
    ratioC[iGraph] = new TCanvas(Form("ratio_%d",graphTypes[iGraph]),"ratioC",600,600);
    gPad->SetRightMargin(0.02);
    gPad->SetTopMargin(0.02);

    TLegend * legend = new TLegend(0.15,0.7,0.5,0.9);
    legend->SetFillColor(0);
    legend->SetLineColor(0);
    legend->SetTextSize(0.03);
    TLegend * legendRatio = new TLegend(0.15,0.2,0.5,0.4);
    legendRatio->SetFillColor(0);
    legendRatio->SetTextSize(0.03);
    TGraphErrors*** dataGraphs = new TGraphErrors**[NGraphs];
    for (Int_t i=0; i<NGraphs; i++)
    {
      dataGraphs[i] = new TGraphErrors*[NHists];
      for (Int_t j=0; j<NHists; j++)
        dataGraphs[i][j] = new TGraphErrors;
    }
  
    for (unsigned int iFile=0; iFile<filesV.size(); iFile++)
    {
      TGraphErrors* graphpT = new TGraphErrors;
      TGraphErrors* graphpTFit = new TGraphErrors;
      TGraphErrors* dataMCRatio = new TGraphErrors();
//      graphCentralityC[iGraph].push_back(new TCanvas(Form("graphCentrality_%d_%d",iFile,graphTypes[iGraph]),"graphCentralityC",800,600));
      ReadGraphs(filesV[iFile].c_str());
      string canvasTitle = yTitles[iGraph] + "_centrality";
      int tmp = 0;
      for (Int_t i=0; i<nHists; i++)
      {
        if (SkipGraph(i))
          continue;
        TGraphAsymmErrors* graphcentrality = (TGraphAsymmErrors*)graphs[graphTypes[iGraph]][i]->Clone();
        if (graphcentrality->GetN() <= 0)
          continue;
        double yC=0, yP=0, yEC=0, yEP=0;
        double* X = graphcentrality->GetX();
        double* Y = graphcentrality->GetY();
        double* yE = graphcentrality->GetEY();
        for (int x=0; x<graphcentrality->GetN(); x++)
        {
          if (X[x] == 5)
          {
            yC = Y[x];
            yEC = yE[x];
          }
          else if (X[x] == 65)
          {
            yP = Y[x];
            yEP = yE[x];
          }
        }
        if (yC==0 || yP==0) {cerr << "Point not found in graph" << endl; return;}
	AddPoint(graphpT, tmp, yC/yP, 0, TMath::Sqrt(yC*yC/yP/yP*(yEC*yEC/yC/yC+yEP*yEP/yP/yP)));
        TF1* lin = new TF1("lin","pol1",0,85);
        double slope = (yP-yC)/60.;
        cerr << slope  << "\t" << slope*0.5 << "\t" << slope*1.5<< endl;
        lin->SetParameter(0,yC);
        lin->SetParameter(1,slope);
        if (slope > 0) lin->SetParLimits(1,slope*0.1,slope*10.0);
//        lin->SetParLimits(0,0,yC*1.2);
//        else lin->SetParLimits(1,slope*10.0,slope*0.1);
        graphcentrality->Fit(lin,"OREX0");
/*        graphCentralityC[iGraph][iFile]->cd();
        graphcentrality->SetMarkerStyle(markers[tmp]);
        graphcentrality->SetMarkerColor(colors[tmp]);
        graphcentrality->SetLineColor(colors[tmp]);
        graphcentrality->GetYaxis()->SetRangeUser(0,3.7);
        graphcentrality->GetYaxis()->SetTitle(yTitles[iGraph].c_str());
        cerr << iFile << "\t" << iGraph << "\t" << i << endl;
        cerr << labelsV[iFile].c_str() << endl;
        graphcentrality->Draw(tmp==0?"AP":"PSAME");
        lin->SetLineColor(colors[tmp]);
        lin->Draw("SAME");
*/        AddPoint(graphpTFit, tmp, lin->Eval(5)/lin->Eval(65), 0, 0);
        if (iFile == filesV.size()-1)
        {
          TString label = graphcentrality->GetTitle();
          if (label.Length() > 0)
          {
            TObjArray* tokens = label.Tokenize("-");
            label.Form("%s-%s", tokens->At(0)->GetName(),tokens->At(1)->GetName());
            label.ReplaceAll(".00", "");
            label.ReplaceAll(".0", "");
            label.ReplaceAll("p_{T,trig}", "p_{T,t}");
            label.ReplaceAll("p_{T,assoc}", "p_{T,a}");
          }
          label += " GeV/c";
          label = Form("%d : ",tmp) + label;
          if (iGraph == 0) labelV.push_back(label);
        }
        if (labelsV.size() != 0) graphcentrality->SetTitle(labelsV[iFile].c_str());
        if (iFile == 0)
          dataGraphs = graphs;
        else
        {
          double yData=0, yEData=0, yMC=0, yEMC=0;
          double* XData = dataGraphs[graphTypes[iGraph]][i]->GetX();
          double* YData = dataGraphs[graphTypes[iGraph]][i]->GetY();
          double* YEData = dataGraphs[graphTypes[iGraph]][i]->GetEY();
          for (int x=0; x<graphcentrality->GetN(); x++)
          {
            bool found = false;
            for (int iData=0; iData<dataGraphs[graphTypes[iGraph]][i]->GetN(); iData++)
            {
              if (X[x] == centrality && XData[iData] == centrality)
              {
                yMC = Y[x];
                yEMC = yE[x];
                yData = YData[iData];
                yEData = YEData[iData];
                found = true;
                break;
              }
            }
            if (found) break;
          }
          if (yMC == 0)
          {
            cerr << "No centrality bin at exactly " << centrality << "%. Please choose between ";
            for (int x=0; x<graphcentrality->GetN(); x++) cerr << X[x] << ", ";
            cerr << endl;
            return;
          }
          AddPoint(dataMCRatio, tmp, yData/yMC, 0, TMath::Sqrt(yData*yData/yMC/yMC*(yEMC*yEMC/yMC/yMC+yEData*yEData/yData/yData)));
        }
        tmp++;
      }
      graphC[iGraph]->cd();
      graphpT->GetYaxis()->SetTitle(yTitles[iGraph].c_str());
      graphpT->GetYaxis()->SetTitleOffset(1.2);
      graphpT->GetXaxis()->CenterTitle(kTRUE);
      graphpT->GetYaxis()->CenterTitle(kTRUE);

//      graphpT->SetTitle("Central (0-10%) / Peripheral (50-80%)");
      graphpT->GetYaxis()->SetRangeUser(yMin,yMax);
      graphpT->GetXaxis()->SetLimits(-0.5,14.5);
      graphpT->SetMarkerStyle(markers[iFile]);
      graphpT->SetMarkerColor(colors[iFile]);
      graphpT->SetLineColor(colors[iFile]);
      graphpT->SetFillColor(0);
      graphpT->SetMarkerSize(1.3);
      if (labelsV.size() != 0) legend->AddEntry(graphpT->Clone(), labelsV[iFile].c_str(),"P");
      graphpT->Draw(iFile==0?"AP":"PSAME");

      graphFitC[iGraph]->cd();
      graphpTFit->GetYaxis()->SetTitle(yTitles[iGraph].c_str());
      graphpTFit->SetTitle("Central (0-10%) / Peripheral (50-80%) - From linear fit");
      graphpTFit->GetYaxis()->SetRangeUser(0.8,5.0);
      graphpTFit->SetMarkerStyle(markers[iFile]);
      graphpTFit->SetMarkerColor(colors[iFile]);
      graphpTFit->SetLineColor(colors[iFile]);
      graphpTFit->SetFillColor(0);
      graphpTFit->SetMarkerSize(1.3);
      graphpTFit->Draw(iFile==0?"AP":"PSAME");

      if (iFile != 0)
      {
        ratioC[iGraph]->cd();
        dataMCRatio->GetYaxis()->SetTitle(yTitlesRatio[iGraph].c_str());
        dataMCRatio->SetTitle(Form("Data / MC, %.0f%%",centrality));
        dataMCRatio->GetYaxis()->SetRangeUser(0.5,2.3);
        dataMCRatio->GetXaxis()->SetLimits(-0.5,9.5);
        dataMCRatio->SetMarkerStyle(markers[iFile]);
        dataMCRatio->SetMarkerColor(colors[iFile]);
        dataMCRatio->SetLineColor(colors[iFile]);
        dataMCRatio->SetFillColor(0);
        dataMCRatio->SetMarkerSize(1.3);
        if (labelsV.size() != 0 && iFile != 0) legendRatio->AddEntry(dataMCRatio->Clone(), labelsV[iFile].c_str(),"P");
        dataMCRatio->Draw(iFile==1?"AP":"PSAME");     
      }
    }
    if (labelsV.size() != 0) 
    {
      graphC[iGraph]->cd();
      TPaveText* paveText = new TPaveText(0.13, 0.9, 0.47, 0.97, "BRNDC");
      paveText->SetTextSize(0.03);
      paveText->SetFillColor(0);
      paveText->SetShadowColor(0);
      paveText->SetBorderSize(0);
      paveText->SetFillStyle(0);
      paveText->SetTextAlign(12);
      paveText->AddText("ALICE");
      paveText->Draw();
      double xLine[4] = {0.5,2.5,5.5,9.5};
      for (int iLine=0; iLine<4; iLine++)
      {
        TLine* l = new TLine(xLine[iLine],yMin,xLine[iLine],yMin+(yMax-yMin)*0.6);
        l->SetLineStyle(2);
        l->Draw();
      }
      legend->Draw();
      graphFitC[iGraph]->cd();
      legend->Draw();
      ratioC[iGraph]->cd();
      legendRatio->Draw();
      TLine* l = new TLine(-0.5,1,9.5,1);
      l->SetLineStyle(2);
      l->Draw();
    }
    graphC[iGraph]->cd();
    for (unsigned int iLabel=0; iLabel<labelV.size(); iLabel++)
      DrawLatex(0.55,0.95-iLabel*0.035,1,labelV[iLabel],0.03);
  }
}

void CompareHistograms(const char* HistFileName1, const char* HistFileName2)
{
  TFile* input1 = TFile::Open(HistFileName1);
  TFile* input2 = TFile::Open(HistFileName2);
  TH2* hist1 = new TH2F;
  TH2* hist2 = new TH2F;

//  Int_t maxLeadingPt = 4;
//  Int_t maxAssocPt = 6;
  Int_t maxLeadingPt = 4;
  Int_t maxAssocPt = 6;
  TGraphErrors** graphsForMean = 0;
  graphsForMean = new TGraphErrors*[NHists];
  for (Int_t i=0; i<NHists; i++)
    graphsForMean[i] = new TGraphErrors;

  Float_t centralityAxisMapping[] =  { 5, 65, 100, 25, 15, 40 };
  Float_t centralityAxisMappingE[] = { 5, 15, 0,    5,  5, 10 };

  for (Int_t i=0; i<maxLeadingPt; i++)
    for (Int_t j=1; j<maxAssocPt; j++)
    {
      if (j-2 >= i) continue; // pTa cannot be larger than pTt
      for (Int_t histId = 0; histId < 6; histId++)
      {
        hist1 = (TH2*) input1->Get(Form("dphi_%d_%d_%d", i, j+1, histId));
        if (!hist1)
        {
          cout << "Hist1 does not exist: " << i << "\t" << j+1 << "\t" << histId << endl;
          continue;
        }
        hist2 = (TH2*) input2->Get(Form("dphi_%d_%d_%d", i, j+1, histId));
        if (!hist2)
        {
          cout << "Hist2 does not exist: " << i << "\t" << j+1 << "\t" << histId << endl;
          continue;
        }
        // NOTE fix normalization. these 2d correlations come out of AliUEHist normalized by dphi bin width, but not deta
        hist1->Scale(1.0 / hist1->GetYaxis()->GetBinWidth(1));
        hist2->Scale(1.0 / hist2->GetYaxis()->GetBinWidth(1));
        TH2* ratio = (TH2*) hist1->Clone("ratio");
//        ratio->Sumw2();
//        hist2->Sumw2();
        ratio->Divide(hist2);
        ratio->GetYaxis()->SetRangeUser(-1.59,1.59);
//        new TCanvas();
//        ratio->Draw("colz");
        double mean = ratio->Integral(1,ratio->GetNbinsX(),ratio->GetYaxis()->FindBin(-1.59),ratio->GetYaxis()->FindBin(1.59))/ratio->GetNbinsX()/(ratio->GetYaxis()->FindBin(1.59)-ratio->GetYaxis()->FindBin(-1.59)+1);
        Int_t graphID = i * (maxAssocPt - 1) + j - 1;
        AddPoint(graphsForMean[graphID], centralityAxisMapping[histId], mean, centralityAxisMappingE[histId], 0);
        graphsForMean[graphID]->SetTitle(hist1->GetTitle());
      }
    }
  DrawCentrality("mean", NHists, graphsForMean, 0., 2., "mean", 0); 
}

void PlotMoreGraphFiles(const char* fileName1, const char* fileName2, const char* fileName3)
{

  Int_t nHists = 24; //NHists;

  TGraphErrors*** graphs1 = 0;
  ReadGraphs(fileName1);
  graphs1 = graphs;

  TGraphErrors*** graphs2 = 0;
  ReadGraphs(fileName2);
  graphs2 = graphs;

  TGraphErrors*** graphs3 = 0;
  ReadGraphs(fileName3);
  graphs3 = graphs;

  TGraphErrors*** graphsMerged = 0;
  graphsMerged = new TGraphErrors**[NGraphs];
  for (Int_t i=0; i<NGraphs; i++)
  {
    graphsMerged[i] = new TGraphErrors*[NHists*3];
    for (Int_t j=0; j<NHists*3; j++)
      graphsMerged[i][j] = new TGraphErrors;
  }

  for (int j=0; j<NGraphs; j++)
  {
    graphsMerged[26][3*j] = graphs1[26][j];
    graphsMerged[26][1+3*j] = graphs2[26][j];
    graphsMerged[26][2+3*j] = graphs3[26][j];
    graphsMerged[27][3*j] = graphs1[27][j];
    graphsMerged[27][1+3*j] = graphs2[27][j];
    graphsMerged[27][2+3*j] = graphs3[27][j];
    graphsMerged[11][3*j] = graphs1[11][j];
    graphsMerged[11][1+3*j] = graphs2[11][j];
    graphsMerged[11][2+3*j] = graphs3[11][j];
  }
  DrawCentrality("phi_rms_centrality_compare", nHists*3, graphsMerged[26], 0.1, 0.6, "#sigma_{#Delta#varphi} (rad.)",0);
  DrawCentrality("eta_rms_centrality_compare", nHists*3, graphsMerged[27], 0.1, 0.9, "#sigma_{#Delta#eta} (rad.)",0);
  DrawCentrality("dip", nHists*3, graphsMerged[11], 0, 0.04, "Depletion yield",0);
  DrawCentrality("phi_rms_centrality_compare_syst", nHists*3, graphsMerged[26], 0.1, 0.6, "#sigma_{#Delta#varphi} (rad.)",graphsMerged[26],0,0,0,0,0);
  DrawCentrality("eta_rms_centrality_compare_syst", nHists*3, graphsMerged[27], 0.1, 0.9, "#sigma_{#Delta#eta} (rad.)",graphsMerged[27],0,0,0,0,1);
  DrawCentrality("dip_syst", nHists*3, graphsMerged[11], 0, 0.04, "Depletion yield", graphsMerged[11],0,0,0,0,2);

}

void DrawResults(const char* fileName = "graphs.root", const char* fileNameWingRemoved = 0, const char* fileName3 = 0, const char* fileName4 = 0, Int_t offset = 16)
{
  TGraphErrors*** graphsWingRemoved = 0;
  if (fileNameWingRemoved)
  {
    ReadGraphs(fileNameWingRemoved);
    graphsWingRemoved = graphs;
  }

  TGraphErrors*** graphs3 = 0;
  if (fileName3)
  {
    ReadGraphs(fileName3);
    graphs3 = graphs;
  }
  
  TGraphErrors*** graphs4 = 0;
  if (fileName4)
  {
    ReadGraphs(fileName4);
    graphs4 = graphs;
  }

  ReadGraphs(fileName);
  
  Int_t nHists = NHists;
  
  Bool_t uncertainties = kFALSE;
//  DrawCentrality("Combinatorial_bg", nHists, graphs[31], 0, 100, "Combinatoral background", (graphsWingRemoved) ? graphsWingRemoved[31] : 0);
  DrawCentrality("yield", nHists, graphs[0], -0.05, 1.4, "N", (graphsWingRemoved) ? graphsWingRemoved[0] : 0, 0, 0, 0, 0, (uncertainties) ? 1 : -1);
//  DrawCentrality("yieldIntegral", nHists, graphs[1], -0.05, 1.0, "Yield from integral", (graphsWingRemoved) ? graphsWingRemoved[1] : 0, 0, 0, 0, 0, (uncertainties) ? 1 : -1);

  DrawCentrality("phi_rms_centrality_compare", nHists, graphs[26], 0, 0.9, "#sigma_{#Delta#varphi} (rad.)", (graphsWingRemoved) ? graphsWingRemoved[10+offset] : 0, 0, (graphs3) ? graphs3[10+offset] : 0, 0, (graphs4) ? graphs4[10+offset] : 0, (uncertainties) ? 1 : -1);
  DrawCentrality("eta_rms_centrality_compare", nHists, graphs[27], 0, 0.9, "#sigma_{#Delta#eta}", (graphsWingRemoved) ? graphsWingRemoved[11+offset] : 0, 0, (graphs3) ? graphs3[11+offset] : 0, 0, (graphs4) ? graphs4[11+offset] : 0, (uncertainties) ? 1 : -1);
//  DrawCentrality("phi_rms_centrality", nHists, graphs[26], 0, 2.0, Form("#sigma_{#Delta#varphi} (rms, %s) (rad.)", fitLabel), (graphsWingRemoved) ? graphsWingRemoved[10+offset] : 0, 0, 0, 0, 0, (uncertainties) ? 1 : -1);
//  DrawCentrality("eta_rms_centrality", nHists, graphs[27], 0, 2.0, Form("#sigma_{#Delta#eta} (rms, %s)", fitLabel), (graphsWingRemoved) ? graphsWingRemoved[11+offset] : 0, 0, 0, 0, 0, (uncertainties) ? 2 : -1);
//  DrawCentrality("phi_eta_rms_compare", nHists, graphs[26], 0, 1.5, "#sigma", graphs[27]);
  DrawCentrality("phi_sigma_centrality", nHists, graphs[4+offset], 0, 1.3, "w_{#Delta#varphi} (rad.)", (graphsWingRemoved) ? graphsWingRemoved[4+offset] : 0, 0, 0, 0, 0, (uncertainties) ? 1 : -1);
  DrawCentrality("eta_sigma_centrality", nHists, graphs[5+offset], 0, 1.3, "w_{#Delta#eta}", (graphsWingRemoved) ? graphsWingRemoved[5+offset] : 0, 0, 0, 0, 0, (uncertainties) ? 2 : -1);
  DrawCentrality("beta_phi", nHists, graphs[19],1.2,2.7,"#gamma_{#varphi}", (graphsWingRemoved) ? graphsWingRemoved[19] : 0);
  DrawCentrality("beta_eta", nHists, graphs[4],1.2,2.7,"#gamma_{#eta}", (graphsWingRemoved) ? graphsWingRemoved[4] : 0);

//  DrawCentrality("chi2_1", nHists, graphs[6+offset], 0.5, 5, "#chi^{2}/ndf (full region)", (graphsWingRemoved) ? graphsWingRemoved[6+offset] : 0);
  DrawCentrality("chi2_2", nHists, graphs[7+offset], 0.5, 10, "#chi^{2}/ndf (peak region)", (graphsWingRemoved) ? graphsWingRemoved[7+offset] : 0);
//  DrawCentrality("chi2_3", nHists, graphs[8+offset], 0.5, 20, "#chi^{2}/ndf (peak region, without center)", (graphsWingRemoved) ? graphsWingRemoved[8+offset] : 0);
//  DrawCentrality("chi2_4", nHists, graphs[9+offset], 0.5, 5, "#chi^{2}/ndf (full region, without center)", (graphsWingRemoved) ? graphsWingRemoved[9+offset] : 0);
    
//  DrawCentrality("v2", nHists, graphs[5], 0, 0.3, "v_{2}", (graphsWingRemoved) ? graphsWingRemoved[5] : 0);
//  DrawCentrality("v3", nHists, graphs[6], 0, 0.2, "v_{3}", (graphsWingRemoved) ? graphsWingRemoved[6] : 0);
//  DrawCentrality("v4", nHists, graphs[7], 0, 0.2, "v_{4}", (graphsWingRemoved) ? graphsWingRemoved[7] : 0);
//  DrawCentrality("v5", nHists, graphs[8], 0, 0.2, "v_{5}", (graphsWingRemoved) ? graphsWingRemoved[8] : 0);
//  DrawCentrality("Yield", nHists, graphs[0], 0, 0.7, "Yield (comparison)", graphs[1]);
//  DrawCentrality("Yield", nHists, graphs[0], 0, 0.7, "Yield (comparison)", (graphsWingRemoved) ? graphsWingRemoved[0] : 0);
  for (int i=0; i<NHists; i++)
    for (int j=0; j<graphs[11][i]->GetN(); j++)
    {
      double x=0, y=0;
      graphs[11][i]->GetPoint(j,x,y);
      graphs[11][i]->SetPoint(j,x,y*100.);
      graphs[11][i]->SetPointError(j,0,graphs[11][i]->GetErrorY(j)*100.);
    }
  for (int i=0; i<NHists && graphsWingRemoved; i++)
    for (int j=0; j<graphsWingRemoved[11][i]->GetN(); j++)
    {
      double x=0, y=0;
      graphsWingRemoved[11][i]->GetPoint(j,x,y);
      graphsWingRemoved[11][i]->SetPoint(j,x,y*100);
      graphsWingRemoved[11][i]->SetPointError(j,0,graphsWingRemoved[11][i]->GetErrorY(j)*100.);
    }
  DrawCentrality("DepletionYield", nHists, graphs[11], 0, 3, "Depletion yield (%)", (graphsWingRemoved) ? graphsWingRemoved[11] : 0);
  DrawCentrality("v2", nHists, graphs[28], 0, 0.8, "V_{2}", (graphsWingRemoved) ? graphsWingRemoved[28] : 0);
  DrawCentrality("v3", nHists, graphs[29], 0, 0.30, "V_{3}", (graphsWingRemoved) ? graphsWingRemoved[29] : 0);
  DrawCentrality("v4", nHists, graphs[30], 0, 0.12, "V_{4}", (graphsWingRemoved) ? graphsWingRemoved[30] : 0);
  DrawCentrality("v5", nHists, graphs[32], 0, 0.04, "V_{5}", (graphsWingRemoved) ? graphsWingRemoved[32] : 0);
  DrawCentrality("norm", nHists, graphs[31], 0, 70, "A", (graphsWingRemoved) ? graphsWingRemoved[31] : 0);
//  DrawCentrality("Dip_all_bins", nHists, graphs[14], -0.1, 0.1, "Dip/yield", (graphsWingRemoved) ? graphsWingRemoved[14] : 0);
//  DrawCentrality("Dip_compare", nHists, graphs[2], -0.15, 0.03, "Dip/yield (comparison)", graphs[11]);
//  DrawCentrality("Dip_1", nHists, graphs[9], -0.15, 0.03, "Dip/yield with one extra bin", (graphsWingRemoved) ? graphsWingRemoved[9] : 0);
//  DrawCentrality("Dip_-1", nHists, graphs[10], -0.15, 0.03, "Dip/yield with one extra bin", (graphsWingRemoved) ? graphsWingRemoved[10] : 0);
//  DrawCentrality("Dip_all", nHists, graphs[3], -0.1, 0.2, "Dip/yield (all bins)", (graphsWingRemoved) ? graphsWingRemoved[3] : 0);
//  DrawCentrality("Normalization", nHists, graphs[0], 0, 0.1, "Normalization", (graphsWingRemoved) ? graphsWingRemoved[0] : 0);
//  DrawCentrality("Integral", nHists, graphs[1], 0, 0.1, "Integral", (graphsWingRemoved) ? graphsWingRemoved[1] : 0);
//  DrawCentrality("phi_rms_CP", nHists, graphs[12], 0, 2, "#sigma_{#Delta#varphi} (C/P)", (graphsWingRemoved) ? graphsWingRemoved[12] : 0);
//  DrawCentrality("eta_rms_CP", nHists, graphs[13], 0, 2, "#sigma_{#Delta#eta} (C/P)", (graphsWingRemoved) ? graphsWingRemoved[13] : 0);
}

void Draw2DExamples(const char* histFileName)
{
  Int_t is[] = { 0, 0, 1, 1, 2, 2, 3 };
  Int_t js[] = { 1, 2, 2, 3, 4, 5, 6 };
  Int_t cs[] = { 0, 3, 5, 1 };
  Float_t outerLimit = 1.59;

  TFile::Open(histFileName);

  TCanvas* canvas = new TCanvas("2d", "2d", 1800, 600);
  canvas->Divide(7, 4);
  Int_t cID = 1;

  for (Int_t c=0; c<4; c++)
  {
    for (Int_t i=0; i<7; i++)
    {
      canvas->cd(cID++);

      TH2* hist = (TH2*) gFile->Get(Form("dphi_%d_%d_%d", is[i], js[i], cs[c]));
      if (!hist)
	continue;

      // NOTE fix normalization. these 2d correlations come out of AliUEHist normalized by dphi bin width, but not deta
      hist->Scale(1.0 / hist->GetYaxis()->GetBinWidth(1));
      hist->GetZaxis()->SetTitle(kCorrFuncTitle);
      hist->GetZaxis()->SetTitleOffset(1.9);
      
      hist->Rebin2D(2, 2);
      hist->Scale(0.25);
  
      TString label(hist->GetTitle());
      label.ReplaceAll(".00", " GeV/c");
      label.ReplaceAll(".0", " GeV/c");
      TObjArray* objArray = label.Tokenize("-");
      TPaveText* paveText = new TPaveText(0.52, 0.72, 0.95, 0.95, "BRNDC");
      paveText->SetTextSize(0.035);
      paveText->SetFillColor(0);
      paveText->SetShadowColor(0);
      paveText->AddText(objArray->At(0)->GetName());
      paveText->AddText(objArray->At(1)->GetName());
      if (objArray->GetEntries() == 4)
	paveText->AddText(Form("Pb-Pb 2.76 TeV %s-%s", objArray->At(2)->GetName(), objArray->At(3)->GetName()));
      else
	paveText->AddText(Form("%s 2.76 TeV", objArray->At(2)->GetName()));
      paveText->AddText("|#eta| < 0.9");
      
      gPad->SetLeftMargin(0.15);
      hist->GetYaxis()->SetRangeUser(-outerLimit, outerLimit);
      hist->GetXaxis()->SetTitleOffset(1.5);
      hist->GetYaxis()->SetTitleOffset(1.7);
      hist->SetStats(0);
      hist->SetTitle("a) Correlation");
      TH2* clone = (TH2*) hist->Clone(Form("%s_clone", hist->GetName()));
    //   clone->GetXaxis()->SetRangeUser(-TMath::Pi() / 2, TMath::Pi() / 2);
      clone->Draw("SURF1");
      paveText->Draw();
    }
  }
}


void CompareGraph(const char* fileName1, const char* fileName2, Int_t iGraph1, Int_t iGraph2, Int_t centrality)
{
  ReadGraphs(fileName1);
  TGraphErrors* graph1 = (TGraphErrors*) graphs[iGraph1][centrality]->Clone("graph1");
  
  ReadGraphs(fileName2);
  TGraphErrors* graph2 = (TGraphErrors*) graphs[iGraph2][centrality]->Clone("graph2");
  
  graph1->Sort();
  graph2->Sort();
  
  TCanvas* c = new TCanvas(Form("%s_%s", fileName1, fileName2), Form("%s_%s", fileName1, fileName2), 800, 800);
  c->Divide(1, 2);
  
  c->cd(1);
  gPad->SetGridx();
  gPad->SetGridy();
  graph1->SetMarkerStyle(24);
  graph1->GetXaxis()->SetRangeUser(0, 110);
  graph1->GetYaxis()->SetRangeUser(0, 2);
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

void Compare2DGraphs(const char* fileName1, const char* fileName2, Int_t ptt, Int_t pta, Int_t histId)
{
  TFile* file1 = TFile::Open(fileName1);
  TFile* file2 = TFile::Open(fileName2);
  int graphID = ptt * 5 + pta - 1;
  TF2* fit1 = (TF2*) file1->Get(Form("fitFunction_%d_%d", graphID, histId));
  TF2* fit2 = (TF2*) file2->Get(Form("fitFunction_%d_%d", graphID, histId));
  TF2* BG1 = new TF2("BG1", "[0] + [1] * TMath::Cos(2. * x) + [2] * TMath::Cos(3. * x) + [3] * TMath::Cos(4.*x) + 0*y", -TMath::Pi() / 2, TMath::Pi() /2, -1.59, 1.59);
  for (int i=0; i<=3; i++)
    BG1->SetParameter(i,fit1->GetParameter(i+5));
  
  TF2* BG2 = new TF2("BG2", "[0] + [1] * TMath::Cos(2. * x) + [2] * TMath::Cos(3. * x) + [3] * TMath::Cos(4.*x) + 0*y", -TMath::Pi() / 2, TMath::Pi() /2, -1.59, 1.59);
  for (int i=0; i<=3; i++)
    BG2->SetParameter(i,fit2->GetParameter(i+5));
  cerr << "Fit1: " << endl;
  fit1->Print();
  cerr << "Fit2: " << endl;
  fit2->Print();

  fit1 = BG1;
  fit2 = BG2;
  fit1->Draw("COLZ");
  TH1* hist1 = fit1->GetHistogram();
  new TCanvas();
  hist1->DrawClone("COLZ");
  new TCanvas();
  fit2->Draw("COLZ");
  TH1* hist2 = fit2->GetHistogram();
  new TCanvas();
  hist2->DrawClone("COLZ");
  hist1->Divide(hist2);
  new TCanvas();
  hist1->DrawClone("COLZ");
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

void SetCOLZHistStyle(TH2* hist, TCanvas* C)
{
  C->SetBottomMargin(0.12);
  C->SetTopMargin(0.02);
  C->SetLeftMargin(0.07);
  C->SetRightMargin(0.22);

  hist->SetTitle("");
  hist->GetXaxis()->SetTitleFont(43);
  hist->GetXaxis()->SetTitleSize(30);
  hist->GetYaxis()->SetTitleFont(43);
  hist->GetYaxis()->SetTitleSize(30);
  hist->GetZaxis()->SetTitleFont(43);
  hist->GetZaxis()->SetTitleSize(30);

  hist->GetXaxis()->SetLabelFont(43);
  hist->GetXaxis()->SetLabelSize(30);
  hist->GetYaxis()->SetLabelFont(43);
  hist->GetYaxis()->SetLabelSize(30);
  hist->GetZaxis()->SetLabelFont(43);
  hist->GetZaxis()->SetLabelSize(30);

  hist->GetXaxis()->SetNdivisions(505);
  hist->GetYaxis()->SetNdivisions(505);
  hist->GetZaxis()->SetNdivisions(504);

  hist->GetXaxis()->SetTitleOffset(0.90);
  hist->GetXaxis()->SetLabelOffset(0.008);
  hist->GetYaxis()->SetTitleOffset(0.50);
  hist->GetYaxis()->SetLabelOffset(0.008);
  hist->GetZaxis()->SetTitleOffset(1.10);
  hist->GetZaxis()->SetLabelOffset(0.008);

  hist->GetXaxis()->CenterTitle(kTRUE);
  hist->GetYaxis()->CenterTitle(kTRUE);
  hist->GetZaxis()->CenterTitle(kTRUE);

  hist->GetYaxis()->SetRangeUser(-1.59,1.59);
  hist->GetXaxis()->SetTitle("#Delta#varphi (rad)");
  hist->GetYaxis()->SetTitle("#Delta#eta");
  hist->GetZaxis()->SetTitle("#frac{1}{#it{N}_{trig}} #frac{d^{2}#it{N}_{assoc}}{d#Delta#etad#Delta#varphi} (rad^{-1})");

//  hist->Scale(1.0 / hist->GetYaxis()->GetBinWidth(1));
  hist->SetStats(0);

  const Int_t nRGBs = 5;
  const Int_t nCont = 99;
  Double_t stops[nRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[nRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[nRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  Double_t blue[nRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  TColor::CreateGradientColorTable(nRGBs, stops, red, green, blue, nCont);
  gStyle->SetNumberContours(nCont);
}

Int_t canvasCount = 0;
void CompareHistDraw(const char* HistFileName1, const char* HistFileName2, Int_t i, Int_t j, Int_t histId, Float_t scaling = 1, Bool_t swapAxis = kFALSE)
{
  Int_t ptId = 0;
  Int_t maxAssocPt = 6;

  TFile* input1 = TFile::Open(HistFileName1);
  TFile* input2 = TFile::Open(HistFileName2);

  TH2* hist1 = (TH2*) input1->Get(Form("dphi_%d_%d_%d", i, j+1, histId));
  if (!hist1)
  {
    cout << "Hist1 does not exist." << endl;
    return;
  }

  TH2* hist2 = (TH2*) input2->Get(Form("dphi_%d_%d_%d", i, j+1, histId));
  if (!hist2)
  {
    cout << "Hist2 does not exist." << endl;
    return;
  }
  // NOTE fix normalization. these 2d correlations come out of AliUEHist normalized by dphi bin width, but not deta
  hist1->Scale(1.0 / hist1->GetYaxis()->GetBinWidth(1));
  hist2->Scale(1.0 / hist2->GetYaxis()->GetBinWidth(1));
  
//   hist2->Scale(1.0 / hist2->GetXaxis()->GetBinWidth(1)); hist1->Rebin2D(4, 2); hist1->Scale(1.0 / 8);

  hist2->Scale(scaling);
  
//   hist1->Rebin2D(2, 2); hist1->Scale(0.25);
//   hist2->Rebin2D(2, 2); hist2->Scale(0.25);
  
  Printf("Integrals: %f %f", hist1->Integral(1, hist1->GetNbinsX(), hist1->GetYaxis()->FindBin(-1.19), hist1->GetYaxis()->FindBin(1.19)), hist2->Integral(1, hist1->GetNbinsX(), hist1->GetYaxis()->FindBin(-1.19), hist1->GetYaxis()->FindBin(1.19)));

  TCanvas *c = new TCanvas(Form("c%d", canvasCount), Form("c%d", canvasCount), 600, 900);
  c->Divide(1,3);
  canvasCount++;

  TH2* ratio = (TH2*) hist1->Clone("ratio");
  if (hist1->GetEntries() < 1e4 || hist2->GetEntries() < 1e4)
  {
    cout << "Hist1: " << hist1->GetEntries() << " entries" << endl << "Hist2: " << hist2->GetEntries() << " entries" << endl;
  }
  ptId = i*(maxAssocPt-1)+j-1;

  if (SkipGraph(ptId))
  {
    cout << "SkipGraph called." << endl;
  }
  for (Int_t x = 1; x<=hist1->GetNbinsX();x++)
  {
    for (Int_t y = 1; y<=hist1->GetNbinsY();y++)
    {
      Int_t binx = hist2->GetXaxis()->FindBin(hist1->GetXaxis()->GetBinCenter(x));
      Int_t biny = hist2->GetYaxis()->FindBin(hist1->GetYaxis()->GetBinCenter(y));
      if (swapAxis)
      {
	biny = hist2->GetYaxis()->FindBin(hist1->GetXaxis()->GetBinCenter(x));
	binx = hist2->GetXaxis()->FindBin(hist1->GetYaxis()->GetBinCenter(y));
      }
      if (hist2->GetBinContent(binx,biny) > 0)
	ratio->SetBinContent(x,y,hist1->GetBinContent(x,y)/hist2->GetBinContent(binx,biny));
    }
  }
  
  hist1->GetYaxis()->SetRangeUser(-1.59, 1.59);
  hist2->GetYaxis()->SetRangeUser(-1.59, 1.59);
  ratio->GetYaxis()->SetRangeUser(-1.59, 1.59);
  
  c->cd(1);
  hist1->Draw("surf1");
  c->cd(2);
  hist2->Draw("surf1");
  c->cd(3);
  ratio->SetStats(0);
  ratio->DrawClone("colz");
  c->SaveAs(Form("ratio_%d_%d_%d.pdf",i,j,histId));  
  TCanvas* C2 = new TCanvas("ratio","ratio",700,500);
  SetCOLZHistStyle(ratio,C2);
  ratio->GetXaxis()->SetRangeUser(-TMath::Pi()/2+0.01,TMath::Pi()/2-0.01);
  ratio->Draw("COLZ");
  ratio->GetZaxis()->SetTitle("Ratio");

  C2->Update();
  TPaletteAxis *palette = (TPaletteAxis*)ratio->GetListOfFunctions()->FindObject("palette");
  palette->SetX1NDC(0.8);
  palette->SetX2NDC(0.85);
  palette->SetY2NDC(0.97);
  C2->Modified();
  C2->Update();

  return;
}

void CompareHistDraw(const char* HistFileName1, Int_t i, Int_t j, Int_t histId, Int_t i2, Int_t j2, Int_t histId2)
{
  TFile* input1 = TFile::Open(HistFileName1);

  TH2* hist1 = (TH2*) input1->Get(Form("dphi_%d_%d_%d", i, j+1, histId));
  if (!hist1)
  {
    cout << "Hist1 does not exist." << endl;
    return;
  }

  TH2* hist2 = (TH2*) input1->Get(Form("dphi_%d_%d_%d", i2, j2+1, histId2));
  if (!hist2)
  {
    cout << "Hist2 does not exist." << endl;
    return;
  }
  // NOTE fix normalization. these 2d correlations come out of AliUEHist normalized by dphi bin width, but not deta
  hist1->Scale(1.0 / hist1->GetYaxis()->GetBinWidth(1));
  hist2->Scale(1.0 / hist2->GetYaxis()->GetBinWidth(1));
  
  Printf("Integrals: %f %f", hist1->Integral(1, hist1->GetNbinsX(), hist1->GetYaxis()->FindBin(-1.19), hist1->GetYaxis()->FindBin(1.19)), hist2->Integral(1, hist1->GetNbinsX(), hist1->GetYaxis()->FindBin(-1.19), hist1->GetYaxis()->FindBin(1.19)));

  TCanvas *c = new TCanvas(Form("c%d", canvasCount), Form("c%d", canvasCount), 600, 900);
  c->Divide(1,3);
  canvasCount++;

  TH2* ratio = (TH2*) hist1->Clone("ratio");
  ratio->Divide(hist2);
  
  if (hist1->GetEntries() < 1e4 || hist2->GetEntries() < 1e4)
  {
    cout << "Hist1: " << hist1->GetEntries() << " entries" << endl << "Hist2: " << hist2->GetEntries() << " entries" << endl;
  }

  hist1->GetYaxis()->SetRangeUser(-1.79, 1.79);
  hist2->GetYaxis()->SetRangeUser(-1.79, 1.79);
  ratio->GetYaxis()->SetRangeUser(-1.79, 1.79);
  
  c->cd(1);
  hist1->Draw("surf1");
  c->cd(2);
  hist2->Draw("surf1");
  c->cd(3);
  ratio->Draw("colz");
}

void CompareSTARpTa(const char* STARFileName, const char* GraphFileName)
{
  TCanvas *c1 = new TCanvas("c1", "STAR comparison in dphi, pTa", 800, 600);
  TCanvas *c2 = new TCanvas("c2", "STAR comparison in deta, pTa", 800, 600);
  TMultiGraph *mg1 = new TMultiGraph();
  TMultiGraph *mg2 = new TMultiGraph();

  c1->SetGrid();
  c2->SetGrid();

  Float_t pTa = 0;
  Float_t dphi[8] = {0};
  Float_t deta[8] = {0};
  Float_t dphierror[8] = {0};
  Float_t detaerror[8] = {0};
  const char* title[8] = {"PYTHIA 62","Cu+Cu 62 0-60%","Au+Au 62 0-80%","PYTHIA 200","d+Au 200 0-95%","Cu+Cu 200 0-60%","STAR Au-Au #sqrt{s_{NN}} = 200 GeV, 3 < p_{T,t} < 6 GeV/#it{c}, 40-80%","STAR Au-Au #sqrt{s_{NN}} = 200 GeV, 3 < p_{T,t} < 6 GeV/#it{c}, 0-12%"};

  ReadGraphs(GraphFileName);
  TGraphErrors*** STARgraphs = new TGraphErrors**[2];
  for (Int_t i=0; i<4; i++)
  {
    STARgraphs[i] = new TGraphErrors*[8];
    for (Int_t j=0; j<8; j++)
      STARgraphs[i][j] = new TGraphErrors;
  }

  ifstream infile(STARFileName);

  TLegend* legend = new TLegend(0.24, 0.65, 0.99, 0.95);
  legend->SetFillColor(0);
  legend->SetTextSize(0.03);

  for (Int_t i=0; i<3; i++)
  {
    infile >> pTa >> dphi[0] >> dphi[1] >> dphierror[1] >> dphi[2] >> dphierror[2] >> dphi[3] >> dphi[4] >> dphierror[4] >> dphi[5] >> dphierror[5] >> dphi[6] >> dphierror[6] >> dphi[7] >> dphierror[7];
    for (Int_t j=0; j<8;j++)
    {
      AddPoint(STARgraphs[0][j], pTa, dphi[j], 0, dphierror[j]); 
      if (i==0)
      {
        if (j==7) 
        {
          STARgraphs[0][j]->SetMarkerStyle(27);
          STARgraphs[0][j]->SetMarkerColor(2);
          STARgraphs[0][j]->SetLineColor(2);
          STARgraphs[0][j]->SetMarkerSize(1.5);
        }
        else 
        {
          STARgraphs[0][j]->SetMarkerStyle(28);
          STARgraphs[0][j]->SetMarkerColor(1);
          STARgraphs[0][j]->SetLineColor(1);
          STARgraphs[0][j]->SetMarkerSize(1.5);
        }
        if (j==7 || j==6) mg1->Add(STARgraphs[0][j]);
        STARgraphs[0][j]->SetTitle(title[j]);
      }
    }
  }

  for (Int_t i=0; i<3; i++)
  {
    infile >> pTa >> deta[0] >> deta[1] >> detaerror[1] >> deta[2] >> detaerror[2] >> deta[3] >> deta[4] >> detaerror[4] >> deta[5] >> detaerror[5] >> deta[6] >> detaerror[6] >> deta[7] >> detaerror[7];
    for (Int_t j=0; j<8;j++) 
    {
      AddPoint(STARgraphs[2][j], pTa, deta[j], 0, detaerror[j]);
      if (i==0)
      {
        if (j==7) 
        {
          STARgraphs[2][j]->SetMarkerStyle(27);
          STARgraphs[2][j]->SetMarkerColor(2);
          STARgraphs[2][j]->SetLineColor(2);
          STARgraphs[2][j]->SetMarkerSize(1.5);
        }
        else 
        {
          STARgraphs[2][j]->SetMarkerStyle(28);
          STARgraphs[2][j]->SetMarkerColor(1);
          STARgraphs[2][j]->SetLineColor(1);
          STARgraphs[2][j]->SetMarkerSize(1.5);
        }
        if (j==7 || j==6) mg2->Add(STARgraphs[2][j]);
        STARgraphs[2][j]->SetTitle(title[j]);
      }
    }
  }

  Float_t pTaAxisMapping[] = { 1.5, 2.5, 3.5, 6.0 };
  Float_t pTaAxisMappingE[] = { 0.5, 0.5, 0.5, 2.0 };

//  Generalized Gaussian
  Int_t graphPhi = 26;
  Int_t graphEta = 27;

  Int_t ptId = 0;
  for (Int_t i = 2; i < 4; i++)
    for (Int_t j = 1; j < 6; j++)
    {
      if ((i == 2 && j > 3) || (i == 3 && j > 4)) continue;
      ptId = i * (6 - 1) + j - 1;
      AddPoint(STARgraphs[1][i-2],pTaAxisMapping[j-1],graphs[graphPhi][ptId]->GetY()[5],pTaAxisMappingE[j-1],graphs[graphPhi][ptId]->GetErrorY(5));
      STARgraphs[1][i-2]->SetMarkerStyle(18+i);
      STARgraphs[1][i-2]->SetMarkerColor(2);
      STARgraphs[1][i-2]->SetLineColor(2);
      STARgraphs[1][i-2]->SetMarkerSize(1.5);
      AddPoint(STARgraphs[1][i],pTaAxisMapping[j-1],graphs[graphPhi][ptId]->GetY()[1],pTaAxisMappingE[j-1],graphs[graphPhi][ptId]->GetErrorY(1));
      STARgraphs[1][i]->SetMarkerStyle(18+i);
      STARgraphs[1][i]->SetMarkerColor(1);
      STARgraphs[1][i]->SetLineColor(1);
      STARgraphs[1][i]->SetMarkerSize(1.5);

      AddPoint(STARgraphs[3][i-2],pTaAxisMapping[j-1],graphs[graphEta][ptId]->GetY()[5],pTaAxisMappingE[j-1],graphs[graphEta][ptId]->GetErrorY(5));
      STARgraphs[3][i-2]->SetMarkerStyle(18+i);
      STARgraphs[3][i-2]->SetMarkerColor(2);
      STARgraphs[3][i-2]->SetLineColor(2);
      STARgraphs[3][i-2]->SetMarkerSize(1.5);
      AddPoint(STARgraphs[3][i],pTaAxisMapping[j-1],graphs[graphEta][ptId]->GetY()[1],pTaAxisMappingE[j-1],graphs[graphEta][ptId]->GetErrorY(1));
      STARgraphs[3][i]->SetMarkerStyle(18+i);
      STARgraphs[3][i]->SetMarkerColor(1);
      STARgraphs[3][i]->SetLineColor(1);
      STARgraphs[3][i]->SetMarkerSize(1.5);
    }
  for (Int_t i = 2; i < 4; i++)
  {
    mg1->Add(STARgraphs[1][i-2]);
    mg2->Add(STARgraphs[3][i-2]);
    legend->AddEntry(mg1->GetListOfGraphs()->At(i), Form("ALICE Pb-Pb #sqrt{s_{NN}} = 2.76 TeV, %.1f < p_{T,t} < %.1f GeV/#it{c}, 0-10%%",pTaAxisMapping[i]-pTaAxisMappingE[i],pTaAxisMapping[i]+pTaAxisMappingE[i]), "P");
  }
  for (Int_t i = 2; i < 4; i++)
  {
    mg1->Add(STARgraphs[1][i]);
    mg2->Add(STARgraphs[3][i]);
    legend->AddEntry(mg1->GetListOfGraphs()->At(i+2), Form("ALICE Pb-Pb #sqrt{s_{NN}} = 2.76 TeV, %.1f < p_{T,t} < %.1f GeV/#it{c}, 50-80%%",pTaAxisMapping[i]-pTaAxisMappingE[i],pTaAxisMapping[i]+pTaAxisMappingE[i]), "P");
  }

  for (Int_t j=7; j>=0;j--)
  {
    if (j==6) legend->AddEntry(mg1->GetListOfGraphs()->At(0), STARgraphs[0][j]->GetTitle(), "P");
    if (j==7) legend->AddEntry(mg1->GetListOfGraphs()->At(1), STARgraphs[0][j]->GetTitle(), "P");
  }
  c1->cd();
  mg1->Draw("AP");
  mg1->GetXaxis()->SetTitle("p_{T,a} (GeV)");
  mg1->GetYaxis()->SetTitle("#sigma_{#Delta#varphi} (rad.)");
  mg1->GetXaxis()->SetRangeUser(1,3);
  mg1->GetYaxis()->SetRangeUser(0.15,0.75);
  mg1->GetYaxis()->SetTitleOffset(1.25);
  legend->Draw();
  for (Int_t i = 2; i < 4; i++)
  {
    // systematic uncertainties
    TGraphErrors* graphSyst1 = (TGraphErrors*) STARgraphs[1][i-2]->Clone();
    TGraphErrors* graphSyst2 = (TGraphErrors*) STARgraphs[1][i]->Clone();
    for (int j=0; j<graphSyst1->GetN(); j++)
    {
      graphSyst1->GetEY()[j] = graphSyst1->GetY()[j] * 0.019;
      graphSyst1->GetEX()[j] = 0.1;
    }
    for (int j=0; j<graphSyst2->GetN(); j++)
    {
      graphSyst2->GetEY()[j] = graphSyst2->GetY()[j] * 0.019;
      graphSyst2->GetEX()[j] = 0.1;
    }
    graphSyst1->SetFillColor(2);
    graphSyst1->SetFillStyle(3002);
    graphSyst1->SetMarkerStyle(0);
    graphSyst1->SetLineColor(0);
    graphSyst2->SetFillColor(1);
    graphSyst2->SetFillStyle(3002);
    graphSyst2->SetMarkerStyle(0);
    graphSyst2->SetLineColor(0);
    graphSyst1->Draw("2 SAME");
    graphSyst2->Draw("2 SAME");
  }
  c2->cd();
  mg2->Draw("AP");
  mg2->GetXaxis()->SetTitle("p_{T,a} (GeV)");
  mg2->GetYaxis()->SetTitle("#sigma_{#Delta#eta}");
  mg2->GetXaxis()->SetRangeUser(1,3);
  mg2->GetYaxis()->SetRangeUser(0.15,0.75);
  mg2->GetYaxis()->SetTitleOffset(1.25);
  legend->Draw();
  for (Int_t i = 2; i < 4; i++)
  {
    // systematic uncertainties
    TGraphErrors* graphSyst1 = (TGraphErrors*) STARgraphs[3][i-2]->Clone();
    TGraphErrors* graphSyst2 = (TGraphErrors*) STARgraphs[3][i]->Clone();
    for (int j=0; j<graphSyst1->GetN(); j++)
    {
      graphSyst1->GetEY()[j] = graphSyst1->GetY()[j] * 0.042;
      graphSyst1->GetEX()[j] = 0.1;
    }
    for (int j=0; j<graphSyst2->GetN(); j++)
    {
      graphSyst2->GetEY()[j] = graphSyst2->GetY()[j] * 0.042;
      graphSyst2->GetEX()[j] = 0.1;
    }
    graphSyst1->SetFillColor(2);
    graphSyst1->SetFillStyle(3002);
    graphSyst1->SetMarkerStyle(0);
    graphSyst1->SetLineColor(0);
    graphSyst2->SetFillColor(1);
    graphSyst2->SetFillStyle(3002);
    graphSyst2->SetMarkerStyle(0);
    graphSyst2->SetLineColor(0);
    graphSyst1->Draw("2 SAME");
    graphSyst2->Draw("2 SAME");
  }
}

void CompareSTARpTt(const char* STARFileName, const char* GraphFileName)
{
  TCanvas *c1 = new TCanvas("c1", "STAR comparison in dphi, pTt", 800, 600);
  TCanvas *c2 = new TCanvas("c2", "STAR comparison in deta, pTt", 800, 600);
  TMultiGraph *mg1 = new TMultiGraph();
  TMultiGraph *mg2 = new TMultiGraph();

  c1->SetGrid();
  c2->SetGrid();

  Float_t pTt = 0;
  Float_t dphi[8] = {0};
  Float_t deta[8] = {0};
  Float_t dphierror[8] = {0};
  Float_t detaerror[8] = {0};
  const char* title[8] = {"PYTHIA 62","Cu+Cu 62 0-60%","Au+Au 62 0-80%","PYTHIA 200","d+Au 200 0-95%","Cu+Cu 200 0-60%","STAR Au-Au #sqrt{s_{NN}} = 200 GeV, 3 < p_{T,t} < 6 GeV/#it{c}, 40-80%","STAR Au-Au #sqrt{s_{NN}} = 200 GeV, 3 < p_{T,t} < 6 GeV/#it{c}, 0-12%"};

  ReadGraphs(GraphFileName);

  TGraphErrors*** STARgraphs = new TGraphErrors**[2];
  for (Int_t i=0; i<4; i++)
  {
    STARgraphs[i] = new TGraphErrors*[20];
    for (Int_t j=0; j<20; j++)
      STARgraphs[i][j] = new TGraphErrors;
  }
  ifstream infile(STARFileName);

  TLegend* legend = new TLegend(0.24, 0.65, 0.99, 0.95);
  legend->SetFillColor(0);
  legend->SetTextSize(0.03);

  for (Int_t i=0; i<3; i++)
  {
    infile >> pTt >> dphi[0] >> dphi[1] >> dphierror[1] >> dphi[2] >> dphierror[2] >> dphi[3] >> dphi[4] >> dphierror[4] >> dphi[5] >> dphierror[5] >> dphi[6] >> dphierror[6] >> dphi[7] >> dphierror[7];
    for (Int_t j=0; j<8;j++)
    {
      AddPoint(STARgraphs[0][j], pTt, dphi[j], 0, dphierror[j]); 
      if (i==0)
      {
        if (j==7) 
        {
          STARgraphs[0][j]->SetMarkerStyle(24);
          STARgraphs[0][j]->SetMarkerColor(2);
          STARgraphs[0][j]->SetLineColor(2);
          STARgraphs[0][j]->SetMarkerSize(1.5);
        }
        else 
        {
          STARgraphs[0][j]->SetMarkerStyle(20+j);
          STARgraphs[0][j]->SetMarkerColor(1);
          STARgraphs[0][j]->SetLineColor(1);
          STARgraphs[0][j]->SetMarkerSize(1.5);
        }
        if (j==7 || j==6) mg1->Add(STARgraphs[0][j]);
        STARgraphs[0][j]->SetTitle(title[j]);
      }
    }
  }

  for (Int_t i=0; i<4; i++)
  {
    infile >> pTt >> dphi[0] >> dphi[1] >> dphierror[1] >> dphi[3] >> dphi[4] >> dphierror[4] >> dphi[5] >> dphierror[5] >> dphi[6] >> dphierror[6] >> dphi[7] >> dphierror[7];
    for (Int_t j=0; j<8;j++)
    {
      if (j==2) continue;
      AddPoint(STARgraphs[0][j], pTt, dphi[j], 0, dphierror[j]);
    }
  }

  for (Int_t i=0; i<3; i++)
  {
    infile >> pTt >> deta[0] >> deta[1] >> detaerror[1] >> deta[2] >> detaerror[2] >> deta[3] >> deta[4] >> detaerror[4] >> deta[5] >> detaerror[5] >> deta[6] >> detaerror[6] >> deta[7] >> detaerror[7];
    for (Int_t j=0; j<8;j++) 
    {
      AddPoint(STARgraphs[2][j], pTt, deta[j], 0, detaerror[j]);
      if (i==0)
      {
        if (j==7) 
        {
          STARgraphs[2][j]->SetMarkerStyle(24);
          STARgraphs[2][j]->SetMarkerColor(2);
          STARgraphs[2][j]->SetLineColor(2);
          STARgraphs[2][j]->SetMarkerSize(1.5);
        }
        else 
        {
          STARgraphs[2][j]->SetMarkerStyle(20+j);
          STARgraphs[2][j]->SetMarkerColor(1);
          STARgraphs[2][j]->SetLineColor(1);
          STARgraphs[2][j]->SetMarkerSize(1.5);
        }
        if (j==7 || j==6) mg2->Add(STARgraphs[2][j]);
        STARgraphs[2][j]->SetTitle(title[j]);
      }
    }
  }

  for (Int_t i=0; i<4; i++)
  {
    infile >> pTt >> deta[0] >> deta[1] >> detaerror[1] >> deta[3] >> deta[4] >> detaerror[4] >> deta[5] >> detaerror[5] >> deta[6] >> detaerror[6] >> deta[7] >> detaerror[7];
    for (Int_t j=0; j<8;j++)
    {
      if (j==2) continue;
      AddPoint(STARgraphs[2][j], pTt, deta[j], 0, detaerror[j]);
    }
  }
  Float_t pTtAxisMapping[] = { 1.5, 2.5, 3.5, 6.0 };
  Float_t pTtAxisMappingE[] = { 0.5, 0.5, 0.5, 2.0 };

// for Generalized Gaussian
  Int_t graphPhi = 26;
  Int_t graphEta = 27;
  for (Int_t i = 2; i < 4; i++)
    for (Int_t j = 1; j < 6; j++)
    {
      if ((i == 2 && j > 3) || (i == 3 && j > 4)) continue;
      Int_t ptId = i * (6 - 1) + j - 1;
      AddPoint(STARgraphs[1][j-1],pTtAxisMapping[i],graphs[graphPhi][ptId]->GetY()[5],pTtAxisMappingE[i],graphs[graphPhi][ptId]->GetErrorY(5));
      STARgraphs[1][j-1]->SetMarkerStyle(19+j);
      STARgraphs[1][j-1]->SetMarkerColor(2);
      STARgraphs[1][j-1]->SetLineColor(2);
      STARgraphs[1][j-1]->SetMarkerSize(1.5);
      AddPoint(STARgraphs[1][5+j-1],pTtAxisMapping[i],graphs[graphPhi][ptId]->GetY()[1],pTtAxisMappingE[i],graphs[graphPhi][ptId]->GetErrorY(1));
      STARgraphs[1][5+j-1]->SetMarkerStyle(19+j);
      STARgraphs[1][5+j-1]->SetMarkerColor(1);
      STARgraphs[1][5+j-1]->SetLineColor(1);
      STARgraphs[1][5+j-1]->SetMarkerSize(1.5);

      AddPoint(STARgraphs[3][j-1],pTtAxisMapping[i],graphs[graphEta][ptId]->GetY()[5],pTtAxisMappingE[i],graphs[graphEta][ptId]->GetErrorY(5));
      STARgraphs[3][j-1]->SetMarkerStyle(19+j);
      STARgraphs[3][j-1]->SetMarkerColor(2);
      STARgraphs[3][j-1]->SetLineColor(2);
      STARgraphs[3][j-1]->SetMarkerSize(1.5);
      AddPoint(STARgraphs[3][5+j-1],pTtAxisMapping[i],graphs[graphEta][ptId]->GetY()[1],pTtAxisMappingE[i],graphs[graphEta][ptId]->GetErrorY(1));
      STARgraphs[3][5+j-1]->SetMarkerStyle(19+j);
      STARgraphs[3][5+j-1]->SetMarkerColor(1);
      STARgraphs[3][5+j-1]->SetLineColor(1);
      STARgraphs[3][5+j-1]->SetMarkerSize(1.5);
    }
  for (Int_t j = 1; j < 6; j++)
  {
    if (j > 4) continue;
    mg1->Add(STARgraphs[1][j-1]);

    mg2->Add(STARgraphs[3][j-1]);
    legend->AddEntry(mg1->GetListOfGraphs()->At(2+j-1), Form("ALICE Pb-Pb #sqrt{s_{NN}} = 2.76 TeV, %.1f < p_{T,t} < %.1f GeV/#it{c}, 0-10%%",pTtAxisMapping[j-1]-pTtAxisMappingE[j-1],pTtAxisMapping[j-1]+pTtAxisMappingE[j-1]), "P");
  }
  for (Int_t j = 1; j < 6; j++)
  {
    if (j > 4) continue;
    mg1->Add(STARgraphs[1][5+j-1]);
    mg2->Add(STARgraphs[3][5+j-1]);
    legend->AddEntry(mg1->GetListOfGraphs()->At(6+j-1), Form("ALICE Pb-Pb #sqrt{s_{NN}} = 2.76 TeV, %.1f < p_{T,t} < %.1f GeV/#it{c}, 0-10%%",pTtAxisMapping[j-1]-pTtAxisMappingE[j-1], pTtAxisMapping[j-1]+pTtAxisMappingE[j-1]), "P");
  }

  for (Int_t j=7; j>=0;j--)
  {
    if (j==6) legend->AddEntry(mg1->GetListOfGraphs()->At(0), STARgraphs[0][j]->GetTitle(), "P");
    if (j==7) legend->AddEntry(mg1->GetListOfGraphs()->At(1), STARgraphs[0][j]->GetTitle(), "P");
  }
  c1->cd();
  mg1->Draw("AP");
  mg1->GetXaxis()->SetTitle("p_{T,t} (GeV)");
  mg1->GetYaxis()->SetTitle(Form("#sigma_{#Delta#varphi} (%s) (rad.)", fitLabel));
  mg1->GetYaxis()->SetRangeUser(0.15,0.6);
  legend->Draw();
  for (Int_t i = 1; i < 6; i++)
  {
    if (i > 4) continue;
    // systematic uncertainties
    TGraphErrors* graphSyst1 = (TGraphErrors*) STARgraphs[1][i-1]->Clone();
    TGraphErrors* graphSyst2 = (TGraphErrors*) STARgraphs[1][5+i-1]->Clone();
    for (int j=0; j<graphSyst1->GetN(); j++)
    {
      graphSyst1->GetEY()[j] = graphSyst1->GetY()[j] * 0.019;
      graphSyst1->GetEX()[j] = 0.1;
    }
    for (int j=0; j<graphSyst2->GetN(); j++)
    {
      graphSyst2->GetEY()[j] = graphSyst2->GetY()[j] * 0.019;
      graphSyst2->GetEX()[j] = 0.1;
    }
    graphSyst1->SetFillColor(2);
    graphSyst1->SetFillStyle(3002);
    graphSyst1->SetMarkerStyle(0);
    graphSyst1->SetLineColor(0);
    graphSyst2->SetFillColor(1);
    graphSyst2->SetFillStyle(3002);
    graphSyst2->SetMarkerStyle(0);
    graphSyst2->SetLineColor(0);
    graphSyst1->Draw("2 SAME");
    graphSyst2->Draw("2 SAME");
  }

  c2->cd();
  mg2->Draw("AP");
  mg2->GetXaxis()->SetTitle("p_{T,t} (GeV)");
  mg2->GetYaxis()->SetTitle("#sigma_{#Delta#eta}");
  mg2->GetYaxis()->SetRangeUser(0.15,0.8);
  legend->Draw();
  for (Int_t i=1; i<6;i++)
  {
    if (i > 4) continue;
    // systematic uncertainties
    TGraphErrors* graphSyst1 = (TGraphErrors*) STARgraphs[3][i-1]->Clone();
    TGraphErrors* graphSyst2 = (TGraphErrors*) STARgraphs[3][5+i-1]->Clone();
    for (int j=0; j<graphSyst1->GetN(); j++)
    {
      graphSyst1->GetEY()[j] = graphSyst1->GetY()[j] * 0.042;
      graphSyst1->GetEX()[j] = 0.1;
    }
    for (int j=0; j<graphSyst2->GetN(); j++)
    {
      graphSyst2->GetEY()[j] = graphSyst2->GetY()[j] * 0.042;
      graphSyst2->GetEX()[j] = 0.1;
    }
    graphSyst1->SetFillColor(2);
    graphSyst1->SetFillStyle(3002);
    graphSyst1->SetMarkerStyle(0);
    graphSyst1->SetLineColor(0);
    graphSyst2->SetFillColor(1);
    graphSyst2->SetFillStyle(3002);
    graphSyst2->SetMarkerStyle(0);
    graphSyst2->SetLineColor(0);
    graphSyst1->Draw("2 SAME");
    graphSyst2->Draw("2 SAME");
  }
}

const int markers[] = { 20, 21, 34, 31, 33, 25, 24, 27, 28, 30, 31, 32, 33, 34, 2, 5};
const int colors[] = { 1, kGreen+1, kRed, kBlue, kOrange-3, kGray+1, kViolet-9, kCyan+1, kMagenta-2, kGreen+3, kGray+1, kOrange+1, 28, 30, 36, 40, 46 };

void SetGraphStyle(TGraph* graph, int iGraph, TCanvas* C)
{
  C->SetBottomMargin(0.15);
  C->SetTopMargin(0.08);
  C->SetLeftMargin(0.13);
  C->SetRightMargin(0.07);
 
  graph->SetTitle("");
  graph->SetMarkerColor(colors[iGraph]);
  graph->SetLineColor(colors[iGraph]);
  graph->SetMarkerStyle(markers[iGraph]);
  graph->SetMarkerSize(1.3);
  if (markers[iGraph] == 27 || markers[iGraph] == 33) graph->SetMarkerSize(1.7);
  graph->GetXaxis()->SetTitleFont(43);
  graph->GetXaxis()->SetTitleSize(30);
  graph->GetYaxis()->SetTitleFont(43);
  graph->GetYaxis()->SetTitleSize(30);
  graph->GetXaxis()->SetLabelFont(43);
  graph->GetXaxis()->SetLabelSize(30);
  graph->GetYaxis()->SetLabelFont(43);
  graph->GetYaxis()->SetLabelSize(30);
  graph->GetXaxis()->SetNdivisions(505);
  graph->GetYaxis()->SetNdivisions(505);
  graph->GetXaxis()->SetTitleOffset(1.10);
  graph->GetXaxis()->SetLabelOffset(0.008);
  graph->GetYaxis()->SetTitleOffset(1.00);
  graph->GetYaxis()->SetLabelOffset(0.008);
  graph->SetFillStyle(0);
}

int marker1[] = {20, 21, 34, 33, 29};
int marker2[] = {24, 25, 28, 27, 30};

void DrawFitFunctionComparison(string names, string legends, int iGraph, TString yTitle, double yLow, double yHigh, bool ratio=false)
{
  vector<TString> fileNames;
  vector<TString> legendNames;
  std::istringstream filesIs(names);
  string filesStr;
  while( filesIs >> filesStr)
    fileNames.push_back(filesStr + "/graphs_wing_removed_noBetaLimit.root");
//    fileNames.push_back(filesStr + "/dphi_corr_norm_wing_removed.root");
  std::istringstream legendsIs(legends);
  string legendsStr;
  while( legendsIs >> legendsStr)
    legendNames.push_back(legendsStr);
  if (fileNames.size() != legendNames.size())
  {
    cerr << "Equal number of files and legends required" << endl;
    return;
  }
  int nFile = fileNames.size();
  TCanvas* c;
  if (!ratio) c = new TCanvas("c","c",700,525);
  else c = new TCanvas("c","c",700,600);
  TPad* pad0 = new TPad("pad0","pad0",0,0.35,1,1);
  TPad* pad1 = new TPad("pad0","pad0",0,0,1,0.35);
  if  (ratio)
  {
    pad0->SetBottomMargin(0);
    pad1->SetBottomMargin(0.32);
    pad0->SetTopMargin(0.08);
    pad1->SetTopMargin(0);
    pad0->SetLeftMargin(0.13);
    pad1->SetLeftMargin(0.13);
    pad0->SetRightMargin(0.07);
    pad1->SetRightMargin(0.07);
    pad0->Draw();
    pad1->Draw();
    pad0->cd();
  }
  TLegend * legend = new TLegend(0.15,0.6,0.95,0.92);
  if (ratio)
  {
    legend->SetY1(0.5);
    legend->SetY2(0.9);
  }
  legend->SetFillColor(0); 
  legend->SetFillStyle(0); 
  legend->SetLineColor(0); 
  legend->SetBorderSize(0);
  legend->SetTextFont(43); 
  legend->SetTextSize(30); 
  if (iGraph != 11) legend->SetNColumns(2);
  TGraphErrors* tmpGHighPT0 = 0;
  TGraphErrors* tmpGLowPT0 = 0;
  for (int i=0; i<nFile; i++)
  {
    if (ratio) pad0->cd();
    cerr << fileNames[i] << endl;
    if (fileNames[i])
      ReadGraphs(fileNames[i]);
    double* X1 = graphs[iGraph][0]->GetX();
    double* Y1 = graphs[iGraph][0]->GetY();
    for (int j=0; j<graphs[iGraph][0]->GetN(); j++)
    {
      graphs[iGraph][0]->SetPoint(j,X1[j]+(i-1)*2,Y1[j]*((iGraph==11)?100.:1.));
      graphs[iGraph][0]->SetPointError(j,0,graphs[iGraph][0]->GetEY()[j]*((iGraph==11)?100.:1.));
    }
    SetGraphStyle(graphs[iGraph][0],i,c);
    graphs[iGraph][0]->SetMarkerColor(1);
    graphs[iGraph][0]->SetLineColor(1);
    graphs[iGraph][0]->SetMarkerStyle(marker1[i]);
    double* X2 = graphs[iGraph][16]->GetX();
    double* Y2 = graphs[iGraph][16]->GetY();
    for (int j=0; j<graphs[iGraph][16]->GetN(); j++)
    {
      graphs[iGraph][16]->SetPoint(j,X2[j]+(i-1)*2,Y2[j]*((iGraph==11)?100.:1.));
      graphs[iGraph][16]->SetPointError(j,0,graphs[iGraph][16]->GetEY()[j]*((iGraph==11)?100.:1.));
    }
    SetGraphStyle(graphs[iGraph][16],i+3,c);
    graphs[iGraph][16]->SetMarkerColor(2);
    graphs[iGraph][16]->SetLineColor(2);
    graphs[iGraph][16]->SetMarkerStyle(marker2[i]);
    graphs[iGraph][16]->SetMarkerSize(1.3);
    graphs[iGraph][0]->GetXaxis()->SetLimits(-5,105);
    graphs[iGraph][0]->GetYaxis()->SetRangeUser(yLow,yHigh);
    graphs[iGraph][0]->GetXaxis()->SetLabelSize(0);
    graphs[iGraph][0]->GetXaxis()->SetTitle("Centrality (%)");
    graphs[iGraph][0]->GetYaxis()->SetTitle(yTitle);
    if (iGraph != 11) 
    {
      legend->AddEntry(graphs[iGraph][0]->Clone(), "     ","P");
      legend->AddEntry(graphs[iGraph][16]->Clone(), legendNames[i],"P");
    }
    else
      legend->AddEntry(graphs[iGraph][0]->Clone(), legendNames[i],"P");
    graphs[iGraph][0]->Draw(i==0?"AP":"PSAME");
    if (iGraph != 11) graphs[iGraph][16]->Draw("PSAME");
    if (ratio)
    {
      if (i == 0)
      {
        tmpGHighPT0 = (TGraphErrors*)graphs[iGraph][16]->Clone("tmpGHighPT0");
        tmpGLowPT0 = (TGraphErrors*)graphs[iGraph][0]->Clone("tmpGLowPT0");
      }
      else
      {
        TGraphErrors* ratioGHighPT0 = (TGraphErrors*)tmpGHighPT0->Clone("ratioGHighPT0");
        ratioGHighPT0->Clear();
        TGraphErrors* ratioGLowPT0 = (TGraphErrors*)tmpGLowPT0->Clone("ratioGLowPT0");
        ratioGLowPT0->Clear();
        TGraphErrors* ratioGHighPT1 = (TGraphErrors*)graphs[iGraph][16]->Clone(Form("ratioGHighPT1_%d",i));
        TGraphErrors* ratioGLowPT1 = (TGraphErrors*)graphs[iGraph][0]->Clone(Form("ratioGLowPT1_%d",i));
        for (int itmp = 0; itmp < tmpGLowPT0->GetN(); itmp++)
        {
          if (iGraph != 11) ratioGHighPT0->SetPoint(itmp,tmpGHighPT0->GetX()[itmp],tmpGHighPT0->GetY()[itmp]/ratioGHighPT1->GetY()[itmp]);
//          ratioGHighPT0->SetPointError();
          ratioGLowPT0->SetPoint(itmp,tmpGLowPT0->GetX()[itmp],tmpGLowPT0->GetY()[itmp]/ratioGLowPT1->GetY()[itmp]);
//          iratioGLowPT0->SetPointError(itmp,0,ratioGLowPT0->GetY()[itmp]*TMath::Sqrt((tmpGLowPT0->GetErrorY(itmp)*tmpGLowPT0->GetErrorY(itmp)/tmpGLowPT0->GetY()[itmp]/tmpGLowPT0->GetY()[itmp]+ratioGLowPT1->GetErrorY(itmp)*ratioGLowPT1->GetErrorY(itmp)/ratioGLowPT1->GetY()[itmp]/ratioGLowPT1->GetY()[itmp])-1.*2*tmpGLowPT0->GetErrorY(itmp)*ratioGLowPT1->GetErrorY(itmp)/tmpGLowPT0->GetY()[itmp]/ratioGLowPT1->GetY()[itmp]));
          ratioGLowPT0->SetPointError(itmp,0,ratioGLowPT0->GetY()[itmp]*TMath::Sqrt((tmpGLowPT0->GetErrorY(itmp)*tmpGLowPT0->GetErrorY(itmp)/tmpGLowPT0->GetY()[itmp]/tmpGLowPT0->GetY()[itmp]+ratioGLowPT1->GetErrorY(itmp)*ratioGLowPT1->GetErrorY(itmp)/ratioGLowPT1->GetY()[itmp]/ratioGLowPT1->GetY()[itmp])-1.*2*tmpGLowPT0->GetErrorY(itmp)*ratioGLowPT1->GetErrorY(itmp)/tmpGLowPT0->GetY()[itmp]/ratioGLowPT1->GetY()[itmp]));
        }
        pad1->cd();
        SetGraphStyle(ratioGLowPT0,i+3,c);
        ratioGLowPT0->GetXaxis()->SetLimits(-5,105);
        ratioGLowPT0->GetYaxis()->SetRangeUser(0.75,1.25);
        ratioGLowPT0->GetXaxis()->SetLabelSize(0);
        ratioGLowPT0->GetXaxis()->SetTitle("Centrality (%)");
        ratioGLowPT0->GetYaxis()->SetTitle("Ratio");
        ratioGLowPT0->GetYaxis()->SetNdivisions(404);
        ratioGLowPT0->GetXaxis()->SetTitleOffset(3);
        ratioGLowPT0->GetXaxis()->SetTickSize(0.055);
        ratioGLowPT0->GetYaxis()->SetTickSize(0.04);
        ratioGLowPT0->SetMarkerColor(1);
        ratioGLowPT0->SetLineColor(1);
        ratioGLowPT0->SetMarkerStyle(marker1[i]);
        ratioGLowPT0->SetMarkerSize(1.3);
        SetGraphStyle(ratioGHighPT0,i+3,c);
        ratioGHighPT0->SetMarkerSize(1.3);
        ratioGHighPT0->SetMarkerColor(2);
        ratioGHighPT0->SetLineColor(2);
        ratioGHighPT0->SetMarkerStyle(marker2[i]);
        ratioGLowPT0->Draw(i==1?"AP":"PSAME");
        ratioGHighPT0->Draw("PSAME");
      }
    }
  }
  if (ratio) pad0->cd();
  TText *t = new TText();
  t->SetTextAlign(32);
  t->SetTextFont(43);
  t->SetTextSize(30);
  int X[] = {2,53,104};
  TString labels[] = {"0","50","pp"};
  if (ratio) pad1->cd();
  for (Int_t i=0;i<3;i++)
  {
    if (!ratio) t->DrawText(X[i],yLow-0.04*(yHigh-yLow),labels[i]);
    else
    {
      t->DrawText(X[i],0.68,labels[i]);
    }
  }
  if (ratio) pad0->cd();
  if (iGraph != 11) legend->SetHeader("Low #it{p}_{T}     #color[2]{High #it{p}_{T}}");
  else legend->SetHeader("Low #it{p}_{T}");
  legend->Draw();

  if (ratio)
  {
    pad1->cd();
    TLine* line = new TLine(-5,1,105,1);
    line->SetLineStyle(2);
//    line->SetLineWidth(1);
    line->Draw();
  }

  c->SaveAs("output.pdf");
}


void DrawSystematics(const char* baseFile, const char* systFile)
{
  ReadGraphs(baseFile);
  TGraphErrors*** graphsBase = graphs;
  if (systFile)
  {
    ReadGraphs(systFile);
  }
  const Int_t NGraphList = 5;
  Int_t graphList[] = { 26, 27, 11, 12, 13};
  string titles[] = { "#sigma_{#Delta#varphi} (def) / #sigma_{#Delta#varphi} (syst)", "#sigma_{#Delta#eta} (def) / #sigma_{#Delta#eta} (syst)", "Depletion yield (def) / Depletion yield (syst)", "#sigma_{CP, #Delta#varphi} (def) / #sigma_{CP, #Delta#varphi} (syst)", "#sigma_{CP, #Delta#eta} (def) / #sigma_{CP, #Delta#eta} (syst)"};
  for (Int_t i=0; i<NGraphList; i++)
  {
    TCanvas* cAll = new TCanvas(Form("All_%d", i), "All", 700, 525);
    for (Int_t j=0; j<NHists; j++)
    {
      if (SkipGraph(j))
        continue;

      TGraphErrors* graph1 = graphsBase[graphList[i]][j];
      TGraphErrors* graph2 = (systFile) ? graphs[graphList[i]][j] : graphsBase[graphList[i]+16][j];
      if (graph1->GetN() == 0)
        continue;

      graph1->Sort();
      graph2->Sort();
      DivideGraphs(graph1, graph2);
      TMultiGraph* multiGraph = 0;
      TMultiGraph* multiGraph2 = 0;
      PrepareGraphs(NHists, graphsBase[graphList[i]], 0, 0, &multiGraph, &multiGraph2, 0);
      TLegend* legend = new TLegend(0.55, 0.65, 0.95, 0.95);
      legend->SetFillColor(0);
      legend->SetFillStyle(0);
      legend->SetLineColor(0);
      legend->SetBorderSize(0);
      legend->SetTextFont(43);
      legend->SetTextSize(30);

      int forEnd = 100;
      if (multiGraph->GetListOfGraphs()->GetEntries() < forEnd) forEnd = multiGraph->GetListOfGraphs()->GetEntries();
      for (Int_t iGraph=0; iGraph<forEnd; iGraph++)
      {
        TGraph* graphTmp = (TGraphErrors*) multiGraph->GetListOfGraphs()->At(iGraph);
        legend->AddEntry(graphTmp, 0, "P");
        graphTmp->GetYaxis()->SetTitle(titles[i].c_str());
        SetGraphStyle(graphTmp,iGraph,cAll);
        graphTmp->GetYaxis()->SetRangeUser(0.8,1.2);
        graphTmp->GetXaxis()->SetLimits(-5, 105);
        graphTmp->DrawClone(iGraph==0?"AP":"SAMEP");
      }
      cAll->cd();
      TLine* line = new TLine(-5,1,105,1);
      line->SetLineStyle(2);
      line->SetLineWidth(3);
      line->Draw();
      legend->Draw();

    }
  }
}


void Set2DHistStyle(TH2* hist, TCanvas* C)
{
  C->SetTheta(23);
  C->SetPhi(31);
  C->SetBottomMargin(0.06);
  C->SetTopMargin(0.02);
  C->SetLeftMargin(0.19);
  C->SetRightMargin(0.02);
  C->SetFillStyle(0);

  hist->SetTitle("");
  hist->GetXaxis()->SetTitleFont(43);
  hist->GetXaxis()->SetTitleSize(30);
  hist->GetYaxis()->SetTitleFont(43);
  hist->GetYaxis()->SetTitleSize(30);
  hist->GetZaxis()->SetTitleFont(43);
  hist->GetZaxis()->SetTitleSize(30);

  hist->GetXaxis()->SetLabelFont(43);
  hist->GetXaxis()->SetLabelSize(30);
  hist->GetYaxis()->SetLabelFont(43);
  hist->GetYaxis()->SetLabelSize(30);
  hist->GetZaxis()->SetLabelFont(43);
  hist->GetZaxis()->SetLabelSize(30);

  hist->GetXaxis()->SetNdivisions(505);
  hist->GetYaxis()->SetNdivisions(505);
  hist->GetZaxis()->SetNdivisions(504);

  hist->GetXaxis()->SetTitleOffset(1.10);
  hist->GetXaxis()->SetLabelOffset(0.008);
  hist->GetYaxis()->SetTitleOffset(1.10);
  hist->GetYaxis()->SetLabelOffset(0.008);
  hist->GetZaxis()->SetTitleOffset(1.30);
  hist->GetZaxis()->SetLabelOffset(0.008);

  hist->GetXaxis()->CenterTitle(kTRUE);
  hist->GetYaxis()->CenterTitle(kTRUE);
  hist->GetZaxis()->CenterTitle(kTRUE);

  hist->GetYaxis()->SetRangeUser(-1.59,1.59);
  hist->GetXaxis()->SetTitle("#Delta#varphi (rad)");
  hist->GetYaxis()->SetTitle("#Delta#eta");
  hist->GetZaxis()->SetTitle("#frac{1}{#it{N}_{trig}} #frac{d^{2}#it{N}_{assoc}}{d#Delta#etad#Delta#varphi} (rad^{-1})");

//  hist->Scale(1.0 / hist->GetYaxis()->GetBinWidth(1));

  hist->SetStats(0);

  const Int_t nRGBs = 5;
  const Int_t nCont = 99;
  Double_t stops[nRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[nRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[nRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  Double_t blue[nRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  TColor::CreateGradientColorTable(nRGBs, stops, red, green, blue, nCont);
  gStyle->SetNumberContours(nCont);
}

void SetHistStyle(TH1* hist, int iHist, TCanvas* C)
{
  hist->SetTitle("");
  hist->SetStats(kFALSE);
  hist->SetMarkerColor(colors[iHist]);
  hist->SetLineColor(colors[iHist]);
  hist->SetMarkerStyle(markers[iHist]);
  hist->SetMarkerSize(1.3);
  if (markers[iHist] == 27 || markers[iHist] == 33) hist->SetMarkerSize(1.7);
  hist->GetXaxis()->SetTitleFont(43);
  hist->GetXaxis()->SetTitleSize(30);
  hist->GetYaxis()->SetTitleFont(43);
  hist->GetYaxis()->SetTitleSize(30);
  hist->GetXaxis()->SetLabelFont(43);
  hist->GetXaxis()->SetLabelSize(30);
  hist->GetYaxis()->SetLabelFont(43);
  hist->GetYaxis()->SetLabelSize(30);
  hist->GetXaxis()->SetNdivisions(505);
  hist->GetYaxis()->SetNdivisions(505);
  hist->GetXaxis()->SetTitleOffset(1.00);
  hist->GetXaxis()->SetLabelOffset(0.008);
  hist->GetYaxis()->SetTitleOffset(0.80);
  hist->GetYaxis()->SetLabelOffset(0.008);
  switch(iHist)
  {
    case 0:
      hist->SetLineStyle(7);
      break;
    case 1:
      hist->SetLineStyle(4);
      break;
    case 2:
      hist->SetLineStyle(10);
      break;
    case 3:
      hist->SetLineStyle(9);
      break;
    case 4:
      hist->SetLineStyle(1);
      break;
  }
  hist->SetLineWidth(4);
  C->SetBottomMargin(0.15);
  C->SetTopMargin(0.08);
  C->SetLeftMargin(0.13);
  C->SetRightMargin(0.07);
}

void DrawExample2D(int trig, int assoc, int cent, bool onlyNearSide, bool noBG = false, TString name = "Data/Default/merge/")
{
  TFile* file = new TFile(name+"/dphi_corr_norm_wing_removed.root","READONLY");
  TH2* hist = (TH2*)file->Get(Form("dphi_%d_%d_%d",trig, assoc+1, cent));
  TCanvas* C1 = new TCanvas("C1","C1",700,500);
  TH2* histRebinned = RebinTTR(hist, true);
//  TH2* histRebinned = (TH2*)hist->Clone();
  cerr << histRebinned->GetTitle() << endl;
  histRebinned->Scale(1.0 / histRebinned->GetYaxis()->GetBinWidth(1));
  Set2DHistStyle(histRebinned,C1);
  if (onlyNearSide) histRebinned->GetXaxis()->SetRangeUser(-TMath::Pi()/2.+0.01,TMath::Pi()/2.);
  if (noBG)
  {
    TFile* graphFile = new TFile(name+"/graphs_wing_removed_noBetaLimit.root","READONLY");
    int graphID = trig * 5 + assoc - 1;
    TF2* fit = (TF2*) graphFile->Get(Form("fitFunction_%d_%d", graphID, cent));
    TF2* BG = new TF2("BG", "[0] + [1] * TMath::Cos(2. * x) + [2] * TMath::Cos(3. * x) + [3] * TMath::Cos(4.*x) + 0*y", -TMath::Pi() / 2, TMath::Pi() /2, -1.59, 1.59);
    for (int i=0; i<=3; i++)
      BG->SetParameter(i,fit->GetParameter(i+5));
    histRebinned->Add(BG,-1.);
    histRebinned->SetMinimum(0);
    histRebinned->GetZaxis()->SetNdivisions(503);
  }

  histRebinned->Draw("SURF2");
  C1->SaveAs("dphi_deta_corr_example.pdf");
  if (!onlyNearSide)
  {
    TCanvas* C2 = new TCanvas("C2","C2",700,500);
    Set2DHistStyle(histRebinned,C2);
    histRebinned->Draw("SURF");
    C2->SaveAs("dphi_deta_corr_example_noColor.pdf");
  }
}

void DrawSameMixed2D(TString path)
{
  TFile* file = new TFile(path+"histFile.root","READONLY");
//  TFile* file = new TFile("Data/Default/2010/histFile.root","READONLY");
  TH2* histSame = (TH2*)file->Get("same_11_4");
  TCanvas* C1 = new TCanvas("C1","C1",700,500);
  TH2* histSameRebinned = RebinTTR(histSame, false);
//  cerr << histSameRebinned->GetTitle() << endl;
//  histSameRebinned->Scale(1.0 / histSameRebinned->GetYaxis()->GetBinWidth(1));
  Set2DHistStyle(histSameRebinned,C1);
  histSameRebinned->GetZaxis()->SetTitle("S(#Delta#varphi,#Delta#eta) (a.u.)");
  histSameRebinned->GetZaxis()->SetTitleOffset(1.5);
//  histSameRebinned->GetZaxis()->SetRangeUser(0,2.2);
//  TGaxis::SetMaxDigits(2);
//  histSameRebinned->GetXaxis()->SetNoExponent(false);
  histSameRebinned->Draw("SURF2");
  C1->SaveAs("corr_hist_same_event.pdf");

  TH2* histMixed = (TH2*)file->Get("mixed_11_4");
  TCanvas* C2 = new TCanvas("C2","C2",700,500);
  TH2* histMixedRebinned = RebinTTR(histMixed, false);
//  histMixedRebinned->Scale(1.0 / histMixedRebinned->GetYaxis()->GetBinWidth(1));
  Set2DHistStyle(histMixedRebinned,C2);
  histMixedRebinned->GetZaxis()->SetTitle("M(#Delta#varphi,#Delta#eta) (a.u.)");
//  histMixedRebinned->GetZaxis()->SetRangeUser(0,2.2);
  histMixedRebinned->Draw("SURF2");
//  latex->Draw();
  C2->SaveAs("corr_hist_mixed_event.pdf");
  TH1* mixedProj = histMixedRebinned->ProjectionX("mixedProj",0,-1,"E");
//  mixedProj->Scale(1./histMixedRebinned->GetYaxis()->GetNbins());
  TCanvas* C3 = new TCanvas("C3","C3",700,525);
  SetHistStyle(mixedProj,0,C3);
  mixedProj->SetLineStyle(0);
  mixedProj->SetLineWidth(1);
  mixedProj->GetYaxis()->SetTitle("M(#Delta#varphi) (a.u.)");
  mixedProj->GetYaxis()->SetTitleOffset(1.1);
  mixedProj->GetYaxis()->SetRangeUser(16.4,17);
  mixedProj->GetXaxis()->CenterTitle(false);
  mixedProj->Draw("");
  C3->SaveAs("corr_hist_mixed_event_proj.pdf");
}

void Divide2D(TString fileName1, TString fileName2, TString name, bool rebin = false)
{
  TFile* file1 = new TFile(fileName1,"READONLY");
  TFile* file2 = new TFile(fileName2,"READONLY");
//  TFile* file = new TFile("Data/Default/2010/histFile.root","READONLY");
  TH2* hist1 = (TH2*)file1->Get(name);
  TH2* hist2 = (TH2*)file2->Get(name);
  TCanvas* C1 = new TCanvas("C1","C1",700,500);
  TCanvas* C2 = new TCanvas("C2","C2",700,500);
  TCanvas* C3 = new TCanvas("C3","C3",700,500);
  TH2* hist1Rebinned = RebinTTR(hist1, rebin);
  TH2* hist2Rebinned = RebinTTR(hist2, rebin);
//  TH2* hist1Rebinned = hist1; //= RebinTTR(hist1, rebin);
//  TH2* hist2Rebinned = hist2; //= RebinTTR(hist2, rebin);

  hist1Rebinned->Scale(1.0 / hist1Rebinned->GetYaxis()->GetBinWidth(1));
  hist2Rebinned->Scale(1.0 / hist2Rebinned->GetYaxis()->GetBinWidth(1));
  Set2DHistStyle(hist1Rebinned,C1);
  Set2DHistStyle(hist2Rebinned,C2);
  Set2DHistStyle(hist2Rebinned,C3);
  C1->cd();
  hist1Rebinned->DrawClone("SURF2");
  C2->cd();
  hist2Rebinned->Draw("SURF2");
  hist1Rebinned->Divide(hist2Rebinned);
  C3->cd();
  hist1Rebinned->Draw("SURF2");
}


void CompareMixed2D(TString path1, TString path2)
{
  TFile* file1 = new TFile(path1+"/histFile.root","READONLY");
  TFile* file2 = new TFile(path2+"/histFile.root","READONLY");

  TH2* histMixed1 = (TH2*)file1->Get("mixed_3_4_3_4_1_6");
  TH2* histMixed2 = (TH2*)file2->Get("mixed_3_4_3_4_1_6");
  TCanvas* C2 = new TCanvas("C2","C2",700,500);
  TH2* histMixedRebinned1 = RebinTTR(histMixed1, false);
  TH2* histMixedRebinned2 = RebinTTR(histMixed2, false);
  histMixedRebinned1->Scale(1.0 / histMixedRebinned1->GetYaxis()->GetBinWidth(1));
  histMixedRebinned2->Scale(1.0 / histMixedRebinned2->GetYaxis()->GetBinWidth(1));
//  Set2DHistStyle(histMixedRebinned1,C2);
//  Set2DHistStyle(histMixedRebinned2,C2);
  histMixedRebinned1->GetYaxis()->SetRangeUser(-1.59,1.59);
  histMixedRebinned2->GetYaxis()->SetRangeUser(-1.59,1.59);
  histMixedRebinned1->SetStats(0);
  histMixedRebinned2->SetStats(0);
  histMixedRebinned1->Divide(histMixedRebinned2);
  histMixedRebinned1->Draw("COLZ");
}


void DrawWing2D(TString fName = "Data/Default/merge/", int ptt = 0, int pta = 1, int cent = 0)
{
  TFile* file = new TFile(fName+"/dphi_corr_norm.root");
  TH2* hist = (TH2*)file->Get(Form("dphi_%d_%d_%d", ptt, pta+1, cent));
  TCanvas* C1 = new TCanvas("C1","C1",700,500);
  TH2* histRebinned = RebinTTR(hist, true);
  cerr << histRebinned->GetTitle() << endl;
  histRebinned->Scale(1.0 / histRebinned->GetYaxis()->GetBinWidth(1));
  SetCOLZHistStyle(histRebinned,C1);
  histRebinned->GetXaxis()->SetRangeUser(1.5,4.7);
  histRebinned->Draw("COLZ");
  C1->SaveAs("wing.pdf");
  TH1* histProj = histRebinned->ProjectionY("histProj",histRebinned->GetXaxis()->FindBin(TMath::Pi()/2),histRebinned->GetXaxis()->FindBin(TMath::Pi()/2*3),"E");
  TCanvas* C2 = new TCanvas("C2","C2",700,525);
  SetHistStyle(histProj,0,C2);
  histProj->SetLineStyle(0);
  histProj->SetLineWidth(1);
  histProj->GetYaxis()->SetNdivisions(503);
  histProj->GetYaxis()->SetTitle("#frac{1}{#it{N}_{trig}} #frac{d#it{N}_{assoc}}{d#Delta#eta}");
  histProj->GetYaxis()->SetTitle("#frac{1}{#it{N}_{trig}} #frac{d#it{N}_{assoc}}{d#Delta#eta}");
//  histProj->GetYaxis()->SetRangeUser(647.8,648.8);
  histProj->GetYaxis()->SetTitleOffset(1.65);
  histProj->Draw();
  TLatex* latex2 = new TLatex(0.48, 0.44, "#frac{#pi}{2} < #Delta#varphi < #frac{3#pi}{2}");
  latex2->SetTextFont(43);
  latex2->SetNDC();
  latex2->SetTextSize(30);
  latex2->SetTextColor(1);
  latex2->Draw();
  TLegend * legend = new TLegend(0.23,0.77,0.62,0.92);
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->SetLineColor(0);
  legend->SetBorderSize(0);
  legend->SetTextFont(43);
  legend->SetTextSize(30);
  legend->AddEntry(histProj,"Without correction");
  TFile* file2 = new TFile(fName+"/dphi_corr_norm_wing_removed.root","READONLY");
  TH2* hist2 = (TH2*)file2->Get(Form("dphi_%d_%d_%d", ptt, pta+1, cent));
  TH2* histRebinned2 = RebinTTR(hist2, true);
  histRebinned2->Scale(1.0 / histRebinned2->GetYaxis()->GetBinWidth(1));
  TH1* histProj2 = histRebinned2->ProjectionY("histProj2",histRebinned2->GetXaxis()->FindBin(TMath::Pi()/2),histRebinned2->GetXaxis()->FindBin(TMath::Pi()/2*3),"E");
  SetHistStyle(histProj2,1,C2);
  C2->SetLeftMargin(0.22);
  histProj2->SetLineStyle(0);
  histProj2->SetLineWidth(1);
  histProj2->Draw("SAME");
  legend->AddEntry(histProj2,"With correction");
  legend->Draw();
  C2->SaveAs("wing_away_side.pdf");
  TH1* histProj3 = histRebinned->ProjectionY("histProj3",histRebinned->GetXaxis()->FindBin(-TMath::Pi()/2),histRebinned->GetXaxis()->FindBin(TMath::Pi()/2),"E");
  TCanvas* C3 = new TCanvas("C3","C3",700,525);
  SetHistStyle(histProj3,0,C3);
  histProj3->SetLineStyle(0);
  histProj3->SetLineWidth(1);
//  histProj3->GetYaxis()->SetRangeUser(664.4,666.6);
  histProj3->GetYaxis()->SetNdivisions(503);
  histProj3->GetYaxis()->SetTitle("#frac{1}{#it{N}_{trig}} #frac{d#it{N}_{assoc}}{d#Delta#eta}");
  histProj3->GetYaxis()->SetTitleOffset(1.65);
  histProj3->Draw();
//  latex->Draw();
  TLatex* latex3 = new TLatex(0.46, 0.44, "- #frac{#pi}{2} < #Delta#varphi < #frac{#pi}{2}");
  latex3->SetTextFont(43);
  latex3->SetNDC();
  latex3->SetTextSize(30);
  latex3->SetTextColor(1);
  latex3->Draw();
  TH1* histProj4 = histRebinned2->ProjectionY("histProj4",histRebinned2->GetXaxis()->FindBin(-TMath::Pi()/2),histRebinned2->GetXaxis()->FindBin(TMath::Pi()/2),"E");
  SetHistStyle(histProj4,1,C3);
  C3->SetLeftMargin(0.22);
  histProj4->SetLineStyle(0);
  histProj4->SetLineWidth(1);
  histProj4->Draw("SAME");
  legend->Draw();
  C3->SaveAs("wing_near_side.pdf");
}

void DrawExclusionRegion2D()
{
  TFile* file = new TFile("Data/2015/MagneticFieldFinal/merge/Default/dphi_corr_norm.root","READONLY");
  TH2* hist = (TH2*)file->Get("dphi_0_2_0");
  TCanvas* C1 = new TCanvas("C1","C1",700,500);
  cerr << hist->GetTitle() << endl;
  hist->Scale(1.0 / hist->GetYaxis()->GetBinWidth(1));
  SetCOLZHistStyle(hist,C1);
  hist->GetXaxis()->SetRangeUser(-0.5,0.5);
  hist->GetYaxis()->SetRangeUser(-0.6,0.6);
  C1->SetLeftMargin(0.09);
  C1->SetRightMargin(0.23);
  hist->GetZaxis()->SetTitleOffset(1.2);
  hist->GetYaxis()->SetTitleOffset(0.7);
  hist->Draw("COLZ");
  C1->SaveAs("exclusionRegion2D.pdf");
}


void DrawGenGaus()
{
  TF1* genGaus = new TF1("geGaus","[0]/2./[1]/TMath::Gamma(1./[0])*TMath::Exp(-1.*(TMath::Power(TMath::Abs(x/[1]),[0])))", -3,3);
  float gamma[] = {1,1.5,2,3};
  TCanvas* C = new TCanvas("C","C",700,525);

  TLegend * legend = new TLegend(0.15,0.59,0.72,0.91);
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->SetLineColor(0);
  legend->SetBorderSize(0);
  legend->SetTextFont(43);
  legend->SetTextSize(30);
//  TF1* twoGaus = new TF1("twoGaus","[0]/TMath::Sqrt(TMath::TwoPi())/[1]*TMath::Exp(-0.5*(x*x/[1]/[1]))+(1-[0])/TMath::Sqrt(TMath::TwoPi())/[2]*TMath::Exp(-0.5*((x*x/[2]/[2])))",-3,3);
//  twoGaus->SetParameter(0,0.5);
//  twoGaus->SetParameter(1,0.5);
//  twoGaus->SetParameter(2,1.5);
/*  TH1F* tmp = new TH1F("tmp","tmp",2000,-3,3);
  for (int j=0;j<2002;j++) tmp->SetBinContent(j,1);
  tmp->TH1::Multiply(twoGaus);
  SetHistStyle(tmp,4,C);
  tmp->GetYaxis()->SetRangeUser(-0.01,0.75);
//  tmp->DrawClone();
*/  for (int i=0; i<4; i++ )
  {
    genGaus->SetTitle("[0]/2./TMath::Gamma(1./[0])*TMath::Exp(-1.*(TMath::Power(TMath::Abs(x/[1]),[0])))");
    genGaus->SetParameter(0,gamma[i]);
    genGaus->SetParameter(1,1);
    TH1F* tmp = new TH1F(Form("tmp_%d",i),"tmp",2000,-3,3);
    for (int j=0;j<2002;j++) tmp->SetBinContent(j,1);
    tmp->TH1::Multiply(genGaus);
    SetHistStyle(tmp,i,C);
    tmp->GetYaxis()->SetRangeUser(-0.01,1);
    tmp->GetXaxis()->SetTitle("x (a.u.)");
    tmp->GetYaxis()->SetTitle("Probability density function");
    if (gamma[i] == 2) legend->AddEntry(tmp->Clone(),Form("#gamma = %0.1f (Gaus)",gamma[i]), "L");
    else if (gamma[i] == 1) legend->AddEntry(tmp->Clone(),Form("#gamma = %0.1f (Exp)",gamma[i]), "L");
    else legend->AddEntry(tmp->Clone(),Form("#gamma = %0.1f",gamma[i]), "L");
    tmp->DrawClone("SAME");
  }
//  legend->AddEntry(tmp->Clone(),"Two Gaus", "L");
  legend->Draw();
  C->SaveAs("gaus.pdf");
}

void DrawFit2D(TString folder = "Data/Default/merge/", int ptt = 0, int pta = 1, int histID = 1)
{
  TFile* file = new TFile(folder + "/dphi_corr_norm_wing_removed.root","READONLY");
  TH2* hist = (TH2*)file->Get(Form("dphi_%d_%d_%d",ptt,pta+1,histID));
  TCanvas* C1 = new TCanvas("C1","C1",700,500);
  TH2* histRebinned = RebinTTR(hist, true);
  cerr << histRebinned->GetTitle() << endl;
  histRebinned->Scale(1.0 / histRebinned->GetYaxis()->GetBinWidth(1));
  Set2DHistStyle(histRebinned,C1);
  histRebinned->Draw("SURF2");
  histRebinned->GetXaxis()->SetRangeUser(-TMath::Pi()/2.,TMath::Pi()/2.);
  float max = histRebinned->GetMaximum();
  histRebinned->SetMaximum(max);
  float min = histRebinned->GetMinimum();
  histRebinned->SetMinimum(min);
  C1->SaveAs("fit1.pdf");
  TCanvas* CProjX = new TCanvas("CProjX","CProjX",700,525);
  TCanvas* CProjY = new TCanvas("CProjY","CProjY",700,525);
  double limitX = TMath::Pi() /2;
  double limitY = 1.6;
  TH1* projX = histRebinned->ProjectionX("projX",histRebinned->GetYaxis()->FindBin(-1.*limitY), histRebinned->GetYaxis()->FindBin(limitY),"e");
  TH1* projY = histRebinned->ProjectionY("projY",histRebinned->GetXaxis()->FindBin(-1.*limitX), histRebinned->GetXaxis()->FindBin(limitX-0.01),"e");
  SetHistStyle(projX,0,CProjX);
  CProjX->SetLeftMargin(0.19);
  CProjX->cd();
  projX->Scale(histRebinned->GetYaxis()->GetBinWidth(1));
  projX->GetYaxis()->SetTitle("#frac{1}{N_{trig}} #frac{dN_{assoc}}{d#Delta#varphi} (rad^{-1})");
  projX->GetYaxis()->SetTitleOffset(1.4);
  projX->Draw("P");
//  latex->SetX(0.73);
//  latex->SetY(0.86);
//  latex->Draw();
  TLatex* latex2 = new TLatex(0.23, 0.86, "|#Delta#eta| < 1.6");
  latex2->SetTextFont(43);
  latex2->SetNDC();
  latex2->SetTextSize(30);
  latex2->SetTextColor(1);
  latex2->Draw();
  SetHistStyle(projY,0,CProjY);
  CProjY->SetLeftMargin(0.19);
  CProjY->cd();
  projY->Scale(histRebinned->GetXaxis()->GetBinWidth(1));
//  projY->GetYaxis()->SetRangeUser(1.99,2.07);
  projY->GetYaxis()->SetTitle("#frac{1}{N_{trig}} #frac{dN_{assoc}}{d#Delta#eta}");
  projY->GetYaxis()->SetTitleOffset(1.4);
  projY->Draw("P");
//  latex->Draw();
  TLatex* latex3 = new TLatex(0.23, 0.86, "|#Delta#varphi| < #pi/2");
  latex3->SetTextFont(43);
  latex3->SetNDC();
  latex3->SetTextSize(30);
  latex3->SetTextColor(1);
  latex3->Draw();
  TLegend * legend = new TLegend(0.40,0.18,0.79,0.43);
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->SetLineColor(0);
  legend->SetBorderSize(0);
  legend->SetTextFont(43);
  legend->AddEntry(projX->Clone(),"Data", "P");
  legend->SetTextSize(30);

  TFile* graphFile = new TFile(folder + "/graphs_wing_removed_noBetaLimit_vnmax5.root","READONLY");
  int graphID = ptt * 5 + pta - 1;
  TF2* fit = (TF2*) graphFile->Get(Form("fitFunction_%d_%d", graphID, histID));
  TCanvas* C8 = new TCanvas("C8","C8",700,500);
  Set2DHistStyle((TH2*)fit->GetHistogram(),C8);
  fit->SetMaximum(max);
  fit->SetMinimum(min);
  fit->Draw("SURF2");
//  latex->Draw();
  C8->SaveAs("fit4.pdf");

  TH2* resGenGaus = (TH2*)histRebinned->Clone("resGenGaus");
  resGenGaus->Add(fit,-1);
  TCanvas* CResGenGaus = new TCanvas("CResGenGaus","CResGenGaus",700,500);
  SetCOLZHistStyle(resGenGaus,CResGenGaus);
  CResGenGaus->SetRightMargin(0.25);
//  resGenGaus->SetMaximum(0.019);
//  resGenGaus->SetMinimum(-0.019);
  resGenGaus->GetZaxis()->SetTitleOffset(1.4);
  resGenGaus->GetZaxis()->SetTitle("Data - fit");
  resGenGaus->Draw("COLZ");
  CResGenGaus->SaveAs("resGenGaus.pdf");

  TF1* fitProjPeakDphi = new TF1("fitProjPeakDphi","[0]*[3]/2./[1]/TMath::Gamma(1./[3])*TMath::Exp(-1.*TMath::Power(TMath::Abs(x/[1]),[3]))", -TMath::Pi() / 2, TMath::Pi() /2);
  for (int i=0; i<=8; i++)
    fitProjPeakDphi->SetParameter(i,fit->GetParameter(i));
  TF1* fitProjPeakDeta = new TF1("fitProjPeakDeta","[0]*[4]/2./[2]/TMath::Gamma(1./[4])*TMath::Exp(-1.*TMath::Power(TMath::Abs(x/[2]),[4]))", -1.6, 1.6);
  for (int i=0; i<=8; i++)
    fitProjPeakDeta->SetParameter(i,fit->GetParameter(i));
  fitProjPeakDeta->SetParameter(0,1);
  fitProjPeakDphi->SetParameter(0,1);
  float integralDeta = fitProjPeakDeta->Integral(-1.6,1.6);
  float integralDphi = fitProjPeakDphi->Integral(-TMath::Pi() / 2, TMath::Pi() /2);
  cerr << integralDphi << "\t" << integralDeta << endl;
  fitProjPeakDeta->SetParameter(0,fit->GetParameter(0)*integralDphi);
  fitProjPeakDphi->SetParameter(0,fit->GetParameter(0)*integralDeta);

  TF1* fitProjDphi = new TF1("fitProjDphi","(3.2*([5] + [6] * TMath::Cos(2. * x) + [7] * TMath::Cos(3. * x) + [8] * TMath::Cos(4.*x)) + [0]*[3]/2./[1]/TMath::Gamma(1./[3])*TMath::Exp(-1.*TMath::Power(TMath::Abs(x/[1]),[3])))", -TMath::Pi() / 2, TMath::Pi() /2);
  for (int i=0; i<=8; i++)
    fitProjDphi->SetParameter(i,fit->GetParameter(i));
  fitProjDphi->SetParameter(0,fit->GetParameter(0)*integralDeta);
  SetHistStyle(fitProjDphi->GetHistogram(),1,CProjX);
  CProjX->SetLeftMargin(0.19);
  CProjX->cd();
  legend->AddEntry(fitProjDphi->GetHistogram()->Clone(),"Gen. Gauss", "L");
  fitProjDphi->GetHistogram()->Draw("SAME");
  CProjX->SaveAs("dphi_fitProj.pdf");

  TF1* fitProjDeta = new TF1("fitProjDeta","(TMath::Pi()*[5] - 2./3.*[7] + [0]*[4]/2./[2]/TMath::Gamma(1./[4])*TMath::Exp(-1.*TMath::Power(TMath::Abs(x/[2]),[4])))", -1.6, 1.6);
  for (int i=0; i<=8; i++)
    fitProjDeta->SetParameter(i,fit->GetParameter(i));
  fitProjDeta->SetParameter(0,fit->GetParameter(0)*integralDphi);
  SetHistStyle(fitProjDeta->GetHistogram(),1,CProjY);
  CProjY->SetLeftMargin(0.19);
  CProjY->cd();
  fitProjDeta->GetHistogram()->Draw("SAME");
  CProjY->SaveAs("deta_fitProj.pdf");

  TF2* BG = new TF2("BG", "[0] + [1] * TMath::Cos(2. * x) + [2] * TMath::Cos(3. * x) + [3] * TMath::Cos(4.*x) + [4] * TMath::Cos(5.*x) + 0*y", -TMath::Pi() / 2, TMath::Pi() /2, -1.59, 1.59);
  for (int i=0; i<=4; i++)
    BG->SetParameter(i,fit->GetParameter(i+5));
  TCanvas* C2 = new TCanvas("C2","C2",700,500);
  Set2DHistStyle((TH2*)BG->GetHistogram(),C2);
  BG->SetMaximum(max);
  BG->SetMinimum(min);
  BG->Draw("SURF2");
//  latex->Draw();
  C2->SaveAs("fit2.pdf");
  TF2* fitPeak = new TF2("fitPeak", "[0]*([3]*[4]/4./[1]/[2]/TMath::Gamma(1./[3])/TMath::Gamma(1./[4]) * TMath::Exp(-1.*(TMath::Power(TMath::Abs(x/[1]),[3])+TMath::Power(TMath::Abs(y/[2]),[4]))))", -TMath::Pi() / 2, TMath::Pi() /2, -1.59, 1.59);
  for (int i=0; i<=4; i++)
    fitPeak->SetParameter(i,fit->GetParameter(i));
  TCanvas* C3 = new TCanvas("C3","C3",700,500);
  Set2DHistStyle((TH2*)fitPeak->GetHistogram(),C3);
  float maxPeak = 0.1;//fitPeak->GetMaximum();
  float minPeak = fitPeak->GetMinimum();
  fitPeak->GetZaxis()->SetNdivisions(502);
  fitPeak->SetMaximum(maxPeak);
  fitPeak->SetMinimum(minPeak);
  fitPeak->Draw("SURF2");
//  latex->Draw();
  C3->SaveAs("fit3.pdf");
/*  TH2* histGenGaus = (TH2*)histRebinned->Clone();
  histGenGaus->Add(BG,-1);
  new TCanvas();
  histGenGaus->Draw("SURF2");*/
  return;
  TFile* graphFile1Gaus = new TFile("Data/Default/merge/graphs_wing_removed_1Gaus.root","READONLY");
  TF2* fit1Gaus = (TF2*) graphFile1Gaus->Get(Form("fitFunction_%d_%d", graphID, histID));
/*  TF2* BG1Gaus = new TF2("BG", "[0] * (1. + [1] * TMath::Cos(x) + [2] * TMath::Cos(2. * x) + [3] * TMath::Cos(3. * x) + [4] * TMath::Cos(4.*x)) + 0*y", -TMath::Pi() / 2, TMath::Pi() /2, -1.59, 1.59);
  for (int i=0; i<=4; i++)
    BG1Gaus->SetParameter(i,fit1Gaus->GetParameter(i+6));
  TH2* hist1Gaus = (TH2*)histRebinned->Clone();
  hist1Gaus->Add(BG1Gaus,-1);
  new TCanvas();
  hist1Gaus->Draw("SURF2");*/
  TF2* peak1Gaus = new TF2("peak1Gaus", "[0]*([3]/TMath::TwoPi()/[1]/[2] * TMath::Exp(-0.5*((x*x/[1]/[1])+(y*y/[2]/[2]))) + (1-[3])/TMath::TwoPi()/[4]/[5] * TMath::Exp(-0.5*((x*x/[4]/[4])+(y*y/[5]/[5]))))", -TMath::Pi() / 2, TMath::Pi() /2, -1.59, 1.59);
  for (int i=0; i<=5; i++)
    peak1Gaus->SetParameter(i,fit1Gaus->GetParameter(i));
//  Set2DHistStyle((TH2*)peak1Gaus->GetHistogram(),C4);
  peak1Gaus->SetMaximum(maxPeak);
  peak1Gaus->SetMinimum(minPeak);
//  peak1GausHist->Divide(fitPeak);
  TCanvas* C5 = new TCanvas("C5","C5",700,500);
  Set2DHistStyle((TH2*)peak1Gaus->GetHistogram(),C5);
  peak1Gaus->GetZaxis()->SetNdivisions(502);
  peak1Gaus->GetHistogram()->DrawClone("SURF2");
//  latex->Draw();
  C5->SaveAs("1Gaus.pdf");

  TH2* res1Gaus = (TH2*)histRebinned->Clone("res1Gaus");
  res1Gaus->Add(fit1Gaus,-1);
  TCanvas* CRes1Gaus = new TCanvas("CRes1Gaus","CRes1Gaus",700,500);
  SetCOLZHistStyle(res1Gaus,CRes1Gaus);
  CRes1Gaus->SetRightMargin(0.25);
  res1Gaus->SetMaximum(0.019);
  res1Gaus->SetMinimum(-0.019);
  res1Gaus->GetZaxis()->SetTitleOffset(1.4);
  res1Gaus->GetZaxis()->SetTitle("Data - fit");
  res1Gaus->Draw("COLZ");
  CRes1Gaus->SaveAs("res1Gaus.pdf");

  TF1* fit1GausProjPeakDphi = new TF1("fit1GausProjPeakDphi","[0]/TMath::Sqrt(TMath::TwoPi())/[1]*TMath::Exp(-0.5*x*x/[1]/[1])", -TMath::Pi() / 2, TMath::Pi() /2);
  for (int i=0; i<=2; i++)
    fit1GausProjPeakDphi->SetParameter(i,peak1Gaus->GetParameter(i));
  TF1* fit1GausProjPeakDeta = new TF1("fit1GausProjPeakDeta","[0]/TMath::Sqrt(TMath::TwoPi())/[2]*TMath::Exp(-0.5*x*x/[2]/[2])", -1.6, 1.6);
  for (int i=0; i<=2; i++)
    fit1GausProjPeakDeta->SetParameter(i,fit1Gaus->GetParameter(i));
  fit1GausProjPeakDeta->SetParameter(0,1);
  fit1GausProjPeakDphi->SetParameter(0,1);
  integralDeta = fit1GausProjPeakDeta->Integral(-1.6,1.6);
  integralDphi = fit1GausProjPeakDphi->Integral(-TMath::Pi() / 2, TMath::Pi() /2);
  cerr << integralDphi << "\t" << integralDeta << endl;
  fit1GausProjPeakDeta->SetParameter(0,fit1Gaus->GetParameter(0)*integralDphi);
  fit1GausProjPeakDphi->SetParameter(0,fit1Gaus->GetParameter(0)*integralDeta);

  TF1* fit1GausProjDphi = new TF1("fit1GausProjDphi","3.2*([6] + [7] * TMath::Cos(2. * x) + [8] * TMath::Cos(3. * x) + [9] * TMath::Cos(4.*x)) + [0]/TMath::Sqrt(TMath::TwoPi())/[1]*TMath::Exp(-0.5*x*x/[1]/[1])", -TMath::Pi() / 2, TMath::Pi() /2);
  for (int i=0; i<=9; i++)
    fit1GausProjDphi->SetParameter(i,fit1Gaus->GetParameter(i));
  fit1GausProjDphi->SetParameter(0,fit1Gaus->GetParameter(0)*integralDeta);
  SetHistStyle(fit1GausProjDphi->GetHistogram(),2,CProjX);
  CProjX->SetLeftMargin(0.19);
  CProjX->cd();
  legend->AddEntry(fit1GausProjDphi->GetHistogram()->Clone(),"1 Gauss", "L");
  fit1GausProjDphi->GetHistogram()->Draw("SAME");
  TF1* fit1GausProjDeta = new TF1("fit1GausProjDeta","TMath::Pi()*[6] - 2./3.*[8] + [0]/TMath::Sqrt(TMath::TwoPi())/[2]*TMath::Exp(-0.5*x*x/[2]/[2])", -1.6, 1.6);
  for (int i=0; i<=9; i++)
    fit1GausProjDeta->SetParameter(i,fit1Gaus->GetParameter(i));
  fit1GausProjDeta->SetParameter(0,fit1Gaus->GetParameter(0)*integralDphi);
  CProjY->cd();
  SetHistStyle(fit1GausProjDeta->GetHistogram(),2,CProjY);
  CProjY->SetLeftMargin(0.19);
  fit1GausProjDeta->GetHistogram()->Draw("SAME");

  TH2* peak1GausHist = (TH2*)peak1Gaus->GetHistogram();
  TCanvas* C4 = new TCanvas("C4","C4",700,500);
  SetCOLZHistStyle((TH2*)peak1Gaus->GetHistogram(),C4);
  peak1GausHist->GetZaxis()->SetTitle("1 Gaus - Gen. Gaus fit");
  peak1GausHist->GetZaxis()->SetTitleOffset(1.5);
  peak1GausHist->Add(fitPeak,-1);
//  peak1GausHist->SetMaximum(0.001);
//  peak1GausHist->SetMinimum(-0.052);
  peak1GausHist->DrawClone("COLZ");
  C4->SaveAs("1Gaus_genGaus.pdf");

  TFile* graphFile2Gaus = new TFile("Data/Default/merge/graphs_wing_removed_2Gaus.root","READONLY");
  TF2* fit2Gaus = (TF2*) graphFile2Gaus->Get(Form("fitFunction_%d_%d", graphID, histID));
  TF2* peak2Gaus = new TF2("peak2Gaus", "[0]*([3]/TMath::TwoPi()/[1]/[2] * TMath::Exp(-0.5*((x*x/[1]/[1])+(y*y/[2]/[2]))) + (1-[3])/TMath::TwoPi()/[4]/[5] * TMath::Exp(-0.5*((x*x/[4]/[4])+(y*y/[5]/[5]))))", -TMath::Pi() / 2, TMath::Pi() /2, -1.59, 1.59);
  for (int i=0; i<=5; i++)
    peak2Gaus->SetParameter(i,fit2Gaus->GetParameter(i));
  peak2Gaus->SetMaximum(maxPeak);
  peak2Gaus->SetMinimum(minPeak);
  TCanvas* C7 = new TCanvas("C7","C7",700,500);
  Set2DHistStyle((TH2*)peak2Gaus->GetHistogram(),C7);
  peak2Gaus->GetZaxis()->SetNdivisions(502);
  peak2Gaus->GetHistogram()->DrawClone("SURF2");
//  latex->Draw();
  C7->SaveAs("2Gaus.pdf");
  TH2* peak2GausHist = (TH2*)peak2Gaus->GetHistogram();
  TCanvas* C6 = new TCanvas("C6","C6",700,500);
  SetCOLZHistStyle((TH2*)peak2Gaus->GetHistogram(),C6);
//  peak2GausHist->Divide(fitPeak);
  peak2GausHist->Add(fitPeak,-1);
  peak2GausHist->GetZaxis()->SetTitle("2 Gaus - Gen. Gaus fit");
  peak2GausHist->GetZaxis()->SetTitleOffset(1.5);
//  peak2GausHist->SetMaximum(0.101);
//  peak2GausHist->SetMinimum(0);
  peak2GausHist->DrawClone("COLZ");
  C6->SaveAs("2Gaus_genGaus.pdf");

  TH2* res2Gaus = (TH2*)histRebinned->Clone("res2Gaus");
  res2Gaus->Add(fit2Gaus,-1);
  TCanvas* CRes2Gaus = new TCanvas("CRes2Gaus","CRes2Gaus",700,500);
  SetCOLZHistStyle(res2Gaus,CRes2Gaus);
  CRes2Gaus->SetRightMargin(0.25);
  res2Gaus->SetMaximum(0.019);
  res2Gaus->SetMinimum(-0.019);
  res2Gaus->GetZaxis()->SetTitleOffset(1.4);
  res2Gaus->GetZaxis()->SetTitle("Data - fit");
  res2Gaus->Draw("COLZ");
  CRes2Gaus->SaveAs("res2Gaus.pdf");
/*
  TF2* BG2Gaus = new TF2("BG", "[0] * (1. + [1] * TMath::Cos(x) + [2] * TMath::Cos(2. * x) + [3] * TMath::Cos(3. * x) + [4] * TMath::Cos(4.*x)) + 0*y", -TMath::Pi() / 2, TMath::Pi() /2, -1.59, 1.59);
  for (int i=0; i<=4; i++)
    BG2Gaus->SetParameter(i,fit2Gaus->GetParameter(i+6));
  TH2* hist2Gaus = (TH2*)histRebinned->Clone();
  hist2Gaus->Add(BG2Gaus,-1);
  new TCanvas();
  hist2Gaus->Draw("SURF2");*/
  TF1* fit2GausProjPeakDphi = new TF1("fit2GausProjPeakDphi","[0]*([3]/TMath::Sqrt(TMath::TwoPi())/[1]*TMath::Exp(-0.5*x*x/[1]/[1])+(1-[3])/TMath::Sqrt(TMath::TwoPi())/[4]*TMath::Exp(-0.5*x*x/[4]/[4]))", -TMath::Pi() / 2, TMath::Pi() /2);
  for (int i=0; i<=5; i++)
    fit2GausProjPeakDphi->SetParameter(i,peak2Gaus->GetParameter(i));
  TF1* fit2GausProjPeakDeta = new TF1("fit2GausProjPeakDeta","[0]*([3]/TMath::Sqrt(TMath::TwoPi())/[2]*TMath::Exp(-0.5*x*x/[2]/[2])+(1-[3])/TMath::Sqrt(TMath::TwoPi())/[5]*TMath::Exp(-0.5*x*x/[5]/[5]))", -1.6, 1.6);
  for (int i=0; i<=5; i++)
    fit2GausProjPeakDeta->SetParameter(i,fit2Gaus->GetParameter(i));
  fit2GausProjPeakDeta->SetParameter(0,1);
  fit2GausProjPeakDphi->SetParameter(0,1);
  integralDeta = fit2GausProjPeakDeta->Integral(-1.6,1.6);
  integralDphi = fit2GausProjPeakDphi->Integral(-TMath::Pi() / 2, TMath::Pi() /2);
  cerr << integralDphi << "\t" << integralDeta << endl;
  fit2GausProjPeakDeta->SetParameter(0,fit2Gaus->GetParameter(0)*integralDphi);
  fit2GausProjPeakDphi->SetParameter(0,fit2Gaus->GetParameter(0)*integralDeta);

  TF1* fit2GausProjDphi = new TF1("fit2GausProjDphi","3.2*([6] + [7] * TMath::Cos(2. * x) + [8] * TMath::Cos(3. * x) + [9] * TMath::Cos(4.*x)) + [0]*([3]/TMath::Sqrt(TMath::TwoPi())/[1]*TMath::Exp(-0.5*x*x/[1]/[1])+(1-[3])/TMath::Sqrt(TMath::TwoPi())/[4]*TMath::Exp(-0.5*x*x/[4]/[4]))", -TMath::Pi() / 2, TMath::Pi() /2);
  for (int i=0; i<=9; i++)
    fit2GausProjDphi->SetParameter(i,fit2Gaus->GetParameter(i));
  fit2GausProjDphi->SetParameter(0,fit2Gaus->GetParameter(0)*integralDeta);
  CProjX->cd();
  SetHistStyle(fit2GausProjDphi->GetHistogram(),3,CProjX);
  CProjX->SetLeftMargin(0.19);
  legend->AddEntry(fit2GausProjDphi->GetHistogram()->Clone(),"2 Gauss", "L");
  fit2GausProjDphi->GetHistogram()->Draw("SAME");
  legend->Draw();
  TF1* fit2GausProjDeta = new TF1("fit2GausProjDeta","TMath::Pi()*[6] - 2./3.*[8] + [0]*([3]/TMath::Sqrt(TMath::TwoPi())/[2]*TMath::Exp(-0.5*x*x/[2]/[2])+(1-[3])/TMath::Sqrt(TMath::TwoPi())/[5]*TMath::Exp(-0.5*x*x/[5]/[5]))", -1.6, 1.6);
  for (int i=0; i<=9; i++)
    fit2GausProjDeta->SetParameter(i,fit2Gaus->GetParameter(i));
  fit2GausProjDeta->SetParameter(0,fit2Gaus->GetParameter(0)*integralDphi);
  CProjY->cd();
  SetHistStyle(fit2GausProjDeta->GetHistogram(),3,CProjY);
  CProjY->SetLeftMargin(0.19);
  fit2GausProjDeta->GetHistogram()->Draw("SAME");
  legend->Draw();
  CProjX->SaveAs("fit_proj_dphi.pdf");
  CProjY->SaveAs("fit_proj_deta.pdf");


}

void CopyAMPT()
{
  TString name[] = {"smON_resON","smON_resOFF","smOFF_resON"};
  TH2* hist = 0;
  TFile* output = new TFile("AMPT.root","RECREATE");
  output->Close();
  for (int i=0; i<3; i++)
  {
    TFile* file = new TFile("AMPT/Eta0.8/Deta1.59/"+name[i]+"/WithPerugia/dphi_corr_norm.root","READONLY");
    hist = (TH2*)file->Get("dphi_0_2_0");
    output = new TFile("AMPT.root","UPDATE");
    hist->Write(name[i]);
    output->Close();
  }
  output->Close();
}

void DrawVsCentrality(TString fileName, int graphID, float yLow, float yHigh, TString yTitle, TString extraLabel="", TString fileName2 = "", float corr = 0)
{
  float X[] = {1.5,53,103};
  TString labels[] = {"0","50","pp"};
  if (extraLabel== "hijing")
  {
    extraLabel = "";
    labels[2] = "100";
  }
  TGraphErrors*** graphs2 = 0;
  if (fileName2 != "")
  {
    ReadGraphs(fileName2);
    graphs2 = graphs;
//    yLow = 0.75;
//    yHigh = 1.25;
  }
  if (fileName)
    ReadGraphs(fileName);

  TString labelpTTPrev = "";
  TString labelpTAPrev = "";

  TCanvas* c = new TCanvas(Form("c_%d",graphID),Form("c_%d",graphID),700,600);
  vector<TPad*> pad;
  pad.push_back(new TPad("pad0","pad",0,0,1,0.69));
  pad.push_back(new TPad("pad1","pad",0,0.65,1,1));
  pad[0]->SetBottomMargin(0.163);
  pad[1]->SetBottomMargin(0);
  for (int iPad=0; iPad<2; iPad++)
  {
    pad[iPad]->SetRightMargin(0.07);
    pad[iPad]->SetLeftMargin(0.13);
    pad[iPad]->SetTopMargin(0);
  }
  pad[1]->Draw();
  pad[0]->Draw();

  int tmp=0;
  int nPT = 1;
  int IPTloc = -1;
  int size = 0;
  vector<TLegend*> legendV;
  vector<TString> labelpTTV, labelpTAV;
  int nHist = 24;
  if (graphID == 11)
  {
    if (fileName2 != "")
    {
//      yLow = 0.5;
//      yHigh = 1.5;
    }
    nHist = 7;
  }
  for (int iPT=0; iPT<nHist; iPT++)
  {
    if (graphs[graphID][iPT]->GetN() == 0) continue;
    double* X1 = graphs[graphID][iPT]->GetX();
    double* Y1 = graphs[graphID][iPT]->GetY();
    graphs[graphID][iPT]->Sort();
    if (fileName2 != "")
    {
      graphs2[graphID][iPT]->Sort();
      DivideGraphs(graphs[graphID][iPT],graphs2[graphID][iPT],corr);
    }
    for (int j=0; j<graphs[graphID][iPT]->GetN(); j++)
    {
      graphs[graphID][iPT]->SetPoint(j,X1[j]+(tmp-5)*0.5,Y1[j]);
      graphs[graphID][iPT]->SetPointError(j,0,graphs[graphID][iPT]->GetEY()[j]);
    }

    TString label = graphs[graphID][iPT]->GetTitle();
    SetGraphStyle(graphs[graphID][iPT],tmp,c);
    graphs[graphID][iPT]->GetXaxis()->SetTitleOffset(1.53);
    if (graphID == 4 || graphID == 19)
      graphs[graphID][iPT]->GetYaxis()->SetTitleOffset(1.1);
    else if (graphID != 12 && graphID != 13)
    {
      graphs[graphID][iPT]->GetXaxis()->SetLimits(-5,105);
      graphs[graphID][iPT]->GetYaxis()->SetTitleOffset(1.23);
    }
    else
    {
      graphs[graphID][iPT]->GetXaxis()->SetLimits(0,10);
      graphs[graphID][iPT]->GetYaxis()->SetTitleOffset(1.2);
      graphs[graphID][iPT]->GetXaxis()->SetNdivisions(0);
    }
    graphs[graphID][iPT]->GetYaxis()->SetRangeUser(yLow,yHigh);
    graphs[graphID][iPT]->GetYaxis()->SetTitle(yTitle+" ");
//      graphs[graphID][iPT]->GetYaxis()->CenterTitle(kTRUE);
    graphs[graphID][iPT]->GetXaxis()->SetLabelSize(0);
    if (graphID != 12 && graphID != 13) graphs[graphID][iPT]->GetXaxis()->SetTitle("Centrality (%)");
    else 
    {
      graphs[graphID][iPT]->GetXaxis()->SetTitle("#it{p}_{T} as indicated by the colors");
      graphs[graphID][iPT]->GetXaxis()->SetTitleOffset(1.);
      graphs[graphID][iPT]->GetXaxis()->CenterTitle(true);
    }
/*    for (int iTmp=0; iTmp<graphs[graphID][iPT]->GetN(); iTmp++)
    {
      double x = 0, y = 0;
      graphs[graphID][iPT]->GetPoint(iTmp,x,y);
      if (x > 95) graphs[graphID][iPT]->RemovePoint(iTmp);
    }
*/    if (label.Length() > 0)
    {
      TObjArray* tokens = label.Tokenize("-");
      TString labelpTT = tokens->At(0)->GetName();
      TString labelpTA = tokens->At(1)->GetName();
      if (labelpTT.CompareTo(labelpTTPrev) != 0)
        label.Form("#lower[0.25]{%s:%sGeV/#it{c}}", tokens->At(0)->GetName(),tokens->At(1)->GetName());
      else if (labelpTA.CompareTo(labelpTAPrev) != 0)
        label.Form("#lower[0.25]{                      %sGeV/#it{c}}", tokens->At(1)->GetName());
      label.ReplaceAll(".00", "");
      label.ReplaceAll(".0", "");
      label.ReplaceAll("p_{T,trig}","#it{p}_{T,trig}");
      label.ReplaceAll("p_{T,assoc}","#it{p}_{T,assoc}");
      labelpTTPrev = labelpTT;
      labelpTAPrev = labelpTA;
    }
    TString labelpTT, labelpTA, labelTmp;
    if (label.Length() > 0)
    {
      label.ReplaceAll("#lower[0.25]{                       ", "");
      label.ReplaceAll("#lower[0.25]{", "");
      label.ReplaceAll(" GeV/#it{c}}", "");
      TObjArray* tokens = label.Tokenize(":");
      if (tokens->GetEntries() > 1)
      {
        labelpTT = tokens->At(0)->GetName();
        labelpTA = tokens->At(1)->GetName();
        tokens = labelpTT.Tokenize("<");
        labelpTT.Form("%s-%s", tokens->At(0)->GetName(),tokens->At(2)->GetName());
        labelpTT.ReplaceAll(" ","");
      }
      else labelpTA = label;
      tokens = labelpTA.Tokenize("<");
      labelpTA.Form("%s-%s", tokens->At(0)->GetName(),tokens->At(2)->GetName());
      labelpTA.ReplaceAll(" ","");
    }
    TString legend;
    if (labelpTT.Length() > 0 && (labelpTTV.size() == 0 || labelpTT.CompareTo(labelpTTV[labelpTTV.size()-1]) != 0))
    {
      IPTloc++;
      labelpTTV.push_back(labelpTT);
      if (IPTloc == 0)      legendV.push_back(new TLegend(0.150,0.58,0.35,0.72));
      else if (IPTloc == 1) legendV.push_back(new TLegend(0.150,0.20,0.35,0.49));
      else if (IPTloc == 2)  legendV.push_back(new TLegend(0.420,0.28,0.62,0.72));
      else                  legendV.push_back(new TLegend(0.700,0.12,0.90,0.72));
      legendV[legendV.size()-1]->SetFillStyle(0);
      legendV[legendV.size()-1]->SetBorderSize(0);
      nPT = 1;
      legend = labelpTTV[IPTloc] + " : " + labelpTA;
    }
    else legend = labelpTTV[IPTloc] + " : " + labelpTA;
    size++;
    legendV[IPTloc]->AddEntry(graphs[graphID][iPT]->Clone(), legend, "P");
    nPT++;
    pad[0]->cd();
    pad[0]->Update();
    graphs[graphID][iPT]->Draw(iPT==0?"AP":"PSAME");
    tmp++;
  }
  for (unsigned int i=0; i<legendV.size(); i++)
  {
    legendV[i]->SetFillColor(0); 
    legendV[i]->SetFillStyle(0); 
    legendV[i]->SetLineColor(0); 
    legendV[i]->SetBorderSize(0);
    legendV[i]->SetTextFont(43); 
    legendV[i]->SetTextSize(30); 
    pad[1]->cd();
    pad[1]->Update();
    legendV[i]->Draw();
  }
  TLine* l = new TLine(0.15,0.75,0.90,0.75);
  l->Draw();
  TLatex* latex = new TLatex(0.25,0.85,"#it{p}_{T,trig} (GeV/#it{c}) : #it{p}_{T,assoc} (GeV/#it{c})");
  latex->SetTextFont(43);
  latex->SetNDC();
  latex->SetTextSize(30);
  latex->SetTextColor(1);
  latex->Draw();
  TBox* box = new TBox(0.130, 0, 0.930, 0.99);
  box->SetFillStyle(0);
  box->SetFillColor(0);
  box->Draw();

  TLatex* latex2 = new TLatex(0.16,0.92,extraLabel);
  latex2->SetTextFont(43);
  latex2->SetNDC();
  latex2->SetTextSize(30);
  latex2->SetTextColor(1);
  pad[0]->cd();
  pad[0]->Update();
  TText *t = new TText();
  t->SetTextAlign(32);
  t->SetTextFont(43);
  t->SetTextSize(30);
  for (Int_t i=0;i<3;i++)
    if (graphID != 12 && graphID != 13) t->DrawText(X[i],yLow-0.045*(yHigh-yLow),labels[i]);
  latex2->DrawClone();
  if (fileName2 != "")
  {
    TLine* line = new TLine(-5,1,105,1);
    if (graphID == 12 || graphID == 13)
    {
      line->SetX1(0);
      line->SetX2(10);
    }
    line->SetLineStyle(2);
    line->SetLineWidth(3);
    line->Draw();
/*    if (graphID == 26)
    {
      TLine* line1 = 0;
      TLine* line2 = 0;
      line1 = new TLine(20,1.041,105,1.041);
      line2 = new TLine(20,0.959,105,0.959);
      line1->SetLineWidth(3);
      line2->SetLineWidth(3);
      line1->SetLineColor(kGreen+2);
      line2->SetLineColor(kGreen+2);
      line1->DrawClone();
      line2->DrawClone();
      line1->SetX1(-5);
      line2->SetX1(-5);
      line1->SetX2(20);
      line2->SetX2(20);
      line1->SetY1(1.071);
      line2->SetY1(0.929);
      line1->SetY2(1.071);
      line2->SetY2(0.929);
      line1->DrawClone();
      line2->DrawClone();
    }
    else if (graphID == 27)
    {
      TLine* line1 = 0;
      TLine* line2 = 0;
      line1 = new TLine(-5,1.215,10,1.215);
      line2 = new TLine(-5,0.785,10,0.785);
      line1->SetLineWidth(3);
      line2->SetLineWidth(3);
      line1->SetLineColor(kGreen+2);
      line2->SetLineColor(kGreen+2);
      line1->DrawClone();
      line2->DrawClone();
      line1->SetX1(10);
      line2->SetX1(10);
      line1->SetX2(20);
      line2->SetX2(20);
      line1->SetY1(1.203);
      line2->SetY1(0.797);
      line1->SetY2(1.203);
      line2->SetY2(0.797);
      line1->DrawClone();
      line2->DrawClone();
      line1->SetX1(20);
      line2->SetX1(20);
      line1->SetX2(105);
      line2->SetX2(105);
      line1->SetY1(1.094);
      line2->SetY1(0.906);
      line1->SetY2(1.094);
      line2->SetY2(0.906);
      line1->DrawClone();
      line2->DrawClone();
    }
    else if (graphID == 11)
    {
      TLine* line1 = 0;
      TLine* line2 = 0;
      line1 = new TLine(-5,1.551,10,1.551);
      line2 = new TLine(-5,0.449,10,0.449);
      line1->SetLineWidth(3);
      line2->SetLineWidth(3);
      line1->SetLineColor(kGreen+2);
      line2->SetLineColor(kGreen+2);
      line1->DrawClone();
      line2->DrawClone();
      line1->SetX1(10);
      line2->SetX1(10);
      line1->SetX2(20);
      line2->SetX2(20);
      line1->SetY1(1.554);
      line2->SetY1(0.446);
      line1->SetY2(1.554);
      line2->SetY2(0.446);
      line1->DrawClone();
      line2->DrawClone();
      line1->SetX1(20);
      line2->SetX1(20);
      line1->SetX2(30);
      line2->SetX2(30);
      line1->SetY1(1.636);
      line2->SetY1(0.364);
      line1->SetY2(1.636);
      line2->SetY2(0.364);
      line1->DrawClone();
      line2->DrawClone();
      line1->SetX1(30);
      line2->SetX1(30);
      line1->SetX2(105);
      line2->SetX2(105);
      line1->SetY1(1.955);
      line2->SetY1(0.045);
      line1->SetY2(1.955);
      line2->SetY2(0.045);
      line1->DrawClone();
      line2->DrawClone();
    }*/
    pad[0]->Update();
  }
  c->SaveAs("output.pdf");
}

void DrawProjection(int type, string names = "", int ptt = 1, int pta = 2, int histID = 5, bool rebin = true, bool noBG = false)
{
  vector<TString> files, legendTitle, folder;
  double limitX = TMath::Pi() /2;
  double limitY = 1.6;
  TLegend * legend = 0;
  switch(type)
  {
    case 0:
      legend = new TLegend(0.22,0.60,0.50,0.75);
      files.push_back("~/Dropbox/peak-shapes/Data/CrossChecks/2010/def1/PosEta/dphi_corr_norm_wing_removed.root");
      legendTitle.push_back("#eta_{trig} > 0");
      files.push_back("~/Dropbox/peak-shapes/Data/CrossChecks/2010/def1/NegEta/dphi_corr_norm_wing_removed.root");
      legendTitle.push_back("#eta_{trig} < 0");
      break;
    case 1:
      legend = new TLegend(0.22,0.60,0.50,0.75);
      files.push_back("Data/Systematics/merge/Eta0.7/dphi_corr_norm_wing_removed.root");
      legendTitle.push_back("|#eta| < 0.7");
      files.push_back("Data/Systematics/merge/AllFilesMerged/dphi_corr_norm_wing_removed.root");
      legendTitle.push_back("|#eta| < 0.8");
      files.push_back("Data/Systematics/merge/def1/Original/dphi_corr_norm_wing_removed.root");
      legendTitle.push_back("|#eta| < 0.9");
      break;
    case 2:
      legend = new TLegend(0.39,0.18,0.67,0.33);
      files.push_back("Data/Default/merge/dphi_corr_norm_wing_removed.root");
      legendTitle.push_back("Default");
      files.push_back("Data/CrossChecks/merge/NoTTR/dphi_corr_norm_wing_removed.root");
      legendTitle.push_back("No two-track cut");
      limitX = 0.3;
      limitY = TMath::Pi()/10.;
      ptt = 1;
      pta = 1;
      histID = 0;
      break;
    case 3:
      if (names == "")
      {
        cerr << "Give a filename" << endl;
        return;
      }
      else
      {
//        TString legendTmp[5] = {"High int. rate","Low int. rate", "|#eta| < 0.9", 	"-+", "+-" };
//        TString legendTmp[5] = {"Hybrid","Global", "|#eta| < 0.7", 	"-+", "+-" };
//        TString legendTmp[5] = {"|#eta| < 0.9","|#eta| < 0.8 (Default)", "|#eta| < 0.7", 	"-+", "+-" };
//        TString legendTmp[5] = {"7 cm (Default)", "3 cm", "Double", 	"-+", "+-" };
//        TString legendTmp[5] = {"Half", "Default", "Double", 	"-+", "+-" };
//        TString legendTmp[5] = {"Default", "Double", "TPC volume", 	"-+", "+-" };
//        TString legendTmp[5] = {"5.02 TeV", "2.76 TeV", "|#eta| < 0.9", 	"-+", "+-" };
        TString legendTmp[5] = {"++", "--", "+-", "-+" };
//        TString legendTmp[5] = {"0.1 < #eta_{trig} < 0.8", "-0.8 < #eta_{trig} < -0.1", "++", 	"-+", "+-" };
//        TString legendTmp[5] = {"Default", "Cut 1", "Cut 2", "Cut 3", "+-" };
        legend = new TLegend(0.61,0.74,0.89,0.89);
        std::istringstream filesIs(names);
        string filesStr;
        int tmp = 0;
        while( filesIs >> filesStr)
        {
          files.push_back(filesStr + "/dphi_corr_norm_wing_removed.root");
          legendTitle.push_back(legendTmp[tmp]);
//          legendTitle.push_back(filesStr);
          folder.push_back(filesStr);
          tmp++;
        }
      }
      break;
  }
  TCanvas* CProjX = new TCanvas("CProjX","CProjX",700,525);
  TCanvas* CProjY = new TCanvas("CProjY","CProjY",700,525);
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->SetLineColor(0);
  legend->SetBorderSize(0);
  legend->SetTextFont(43);
  legend->SetTextSize(30);
  for (unsigned int iFile=0; iFile<files.size(); iFile++)
  {
    TFile* file = new TFile(files[iFile],"READONLY");
    TH2* hist = (TH2*)file->Get(Form("dphi_%d_%d_%d",ptt,pta+1,histID));
    TH2* histRebinned;
    if (rebin) histRebinned = RebinTTR(hist, true);
    else histRebinned = hist;
    cerr << histRebinned->GetTitle() << endl;
    histRebinned->Scale(1.0 / histRebinned->GetYaxis()->GetBinWidth(1));
    histRebinned->GetXaxis()->SetRangeUser(-TMath::Pi()/2.+0.01,TMath::Pi()/2.-0.01);
    Set2DHistStyle(histRebinned,CProjX);
    if (noBG)
    {
      TFile* graphFile = new TFile(folder[iFile] + "/graphs_wing_removed_noBetaLimit_vnmax5.root","READONLY");
      int graphID = ptt * 5 + pta - 1;
      TF2* fit = (TF2*) graphFile->Get(Form("fitFunction_%d_%d", graphID, histID));
      TF2* func = new TF2("func", DeltaPhiWidth2DFitFunction, -1*limitX, limitX, -1*limitY, limitY, 5+vnmax);
      for (int i=0; i<5+vnmax; i++)
        func->SetParameter(i, fit->GetParameter(i));
      func->SetParameter(0,0);
      histRebinned->Add(func,-1.);
    }
    TH1* projX = histRebinned->ProjectionX(Form("projX_%d",iFile),histRebinned->GetYaxis()->FindBin(-1.*limitY+0.01), histRebinned->GetYaxis()->FindBin(limitY-0.01),"e");
    TH1* projY = histRebinned->ProjectionY(Form("projY_%d",iFile),histRebinned->GetXaxis()->FindBin(-1.*limitX), histRebinned->GetXaxis()->FindBin(limitX-0.01),"e");
    SetHistStyle(projX,iFile,CProjX);
    CProjX->SetLeftMargin(0.20);
    CProjX->cd();
    projX->Scale(histRebinned->GetYaxis()->GetBinWidth(1));
//    if (iFile == 0) projX->Scale(1.0/1.8);
//    if (iFile == 1) projX->Scale(1.0/1.6);
//    if (iFile == 2) projX->Scale(1.0/1.4);
    projX->SetLineStyle(0);
    projX->SetLineWidth(1);
    projX->GetYaxis()->SetTitleOffset(1.5);
    projX->GetYaxis()->SetTitle("#frac{1}{N_{trig}} #frac{dN_{assoc}}{d#Delta#varphi} (rad^{-1})");
    projX->Draw(iFile==0?"P":"PSAME");
    projX->Draw(iFile==0?"P":"PSAME");
    legend->AddEntry(projX->Clone(),legendTitle[iFile], "P");
    SetHistStyle(projY,iFile,CProjY);
    CProjY->SetLeftMargin(0.20);
    CProjY->cd();
    projY->Scale(histRebinned->GetXaxis()->GetBinWidth(1));
//    if (iFile == 0) projY->Scale(1.0/1.8);
//    if (iFile == 1) projY->Scale(1.0/1.6);
//    if (iFile == 2) projY->Scale(1.0/1.4);
    projY->SetLineStyle(0);
    projY->SetLineWidth(1);
    if (type < 2) projY->GetYaxis()->SetRangeUser(1.955,2.04);
    projY->GetXaxis()->SetRangeUser(-1.59,1.59);
    projY->GetYaxis()->SetTitle("#frac{1}{N_{trig}} #frac{dN_{assoc}}{d#Delta#eta}");
    projY->GetYaxis()->SetTitleOffset(1.5);
    projY->DrawClone(iFile==0?"P":"PSAME");

  }
  CProjX->cd();
  TLatex* latex2 = 0;
  if (type != 2) latex2 = new TLatex(0.23, 0.86, "|#Delta#eta| < 1.6");
  else latex2 = new TLatex(0.23, 0.86, "|#Delta#eta| < 0.3");
  latex2->SetTextFont(43);
  latex2->SetNDC();
  latex2->SetTextSize(30);
  latex2->SetTextColor(1);
  latex2->Draw();
  legend->Draw();

  CProjY->cd();
//  latex->Draw();
  TLatex* latex3 = 0;
  if (type != 2) latex3 = new TLatex(0.23, 0.86, "|#Delta#varphi| < #pi/2");
  else latex3 = new TLatex(0.23, 0.86, "|#Delta#varphi| < #pi/10");
  latex3->SetTextFont(43);
  latex3->SetNDC();
  latex3->SetTextSize(30);
  latex3->SetTextColor(1);
  latex3->Draw();
  legend->Draw();
  CProjX->SaveAs("dphi_pos_neg_eta.pdf");
  CProjY->SaveAs("deta_pos_neg_eta.pdf");

}

void PlotEff(TString fileName)
{
  TFile* file = new TFile(fileName,"READONLY");
  THn* hist = (THn*)file->Get("correction");
  hist->GetAxis(1)->SetRangeUser(1.01,1.019);
  hist->GetAxis(2)->SetRangeUser(0,9);
  new TCanvas();
  for (int i=0; i<8; i++)
  {
    cerr << -0.79+i*0.2 << "\t" << -0.61+i*0.2 << endl;
    TH1* proj = hist->Projection(0);
    proj->SetName(Form("proj%d",i));
    proj->SetMarkerStyle(20+i);
    proj->DrawClone(i==0?"":"SAME");
  }
}

