#include "TF2.h"
#include "TH2F.h"
#include "TH3.h"
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
#include "TMultiGraph.h"
#include "TPaveText.h"
#include "TRandom.h"
#include "TSystem.h"
#include "TASImage.h"
#include "TObjString.h"
#include "TFitResult.h"
#include "fstream"
#include "TMatrixTBase.h"

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
// 32     33              34      
// IAAFit,Yield(integral),IAAHist,
//        35      36       37    38    39      40       41    42    43       44
// 1D fit:normphi,norm2phi,dphi1,dphi2,normeta,norm2eta,deta1,deta2,chi2phi,chi2eta
// if fitting two 1D Gaussians phi rms is stored in 37 and eta rms in 41 without calling the CalculateRMS function
const Int_t NGraphs = 45;
const Int_t NHists = 6*4; // pt index
TGraphErrors*** graphs = 0;
const char* kCorrFuncTitle = "1/N_{trig} dN_{assoc}/d#Delta#etad#Delta#varphi (1/rad.)";
const char* kProjYieldTitlePhi = "1/N_{trig} dN_{assoc}/d#Delta#varphi (1/rad.)";
const char* kProjYieldTitleEta = "1/N_{trig} dN_{assoc}/d#Delta#eta";
const char* kProjYieldTitlePhiOrEta = "1/N_{trig} dN_{assoc}/d#Delta#varphi (1/rad.) , dN_{assoc}/d#Delta#eta";
TString fgFolder = "tmpresults";
const char* fitLabel = "fit";

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

void WriteGraphs(const char* outputFileName = "graphs.root")
{
  TFile::Open(outputFileName, "RECREATE");
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

void DrawALICELogo(Bool_t prel, Float_t x1, Float_t y1, Float_t x2, Float_t y2, Bool_t debug = kFALSE)
{
  // correct for aspect ratio of figure plus aspect ratio of pad (coordinates are NDC!)
//   Printf("%d %f %d %f", gPad->GetCanvas()->GetWindowHeight(), gPad->GetHNDC(), gPad->GetCanvas()->GetWindowWidth(), gPad->GetWNDC());
//   x2 = x1 + (y2 - y1) * (620. / 671) * gPad->GetCanvas()->GetWindowHeight() * gPad->GetHNDC() / (gPad->GetWNDC() * gPad->GetCanvas()->GetWindowWidth());
  x2 = x1 + (y2 - y1) * (466. / 523) * gPad->GetWh() * gPad->GetHNDC() / (gPad->GetWNDC() * gPad->GetWw());
//   x2 = x1 + (y2 - y1) * (620. / 671) * gPad->GetWh() / gPad->GetWw();

//   Printf("%f %f %f %f", x1, x2, y1, y2);
  
  TPad *myPadLogo = new TPad("myPadLogo", "Pad for ALICE Logo", x1, y1, x2, y2);
  if (debug)
    myPadLogo->SetFillColor(2); // color to first figure out where is the pad then comment !
  myPadLogo->SetLeftMargin(0);
  myPadLogo->SetTopMargin(0);
  myPadLogo->SetRightMargin(0);
  myPadLogo->SetBottomMargin(0);
  myPadLogo->Draw();
  myPadLogo->cd();
//   TASImage *myAliceLogo = new TASImage("~/alice_logo_transparent.png");
  TASImage *myAliceLogo = new TASImage((prel) ? "~/alice_logo_preliminary.eps" : "~/alice_logo_performance.eps");
  myAliceLogo->Draw();
}

void logotest()
{
  TCanvas* c = new TCanvas("c", "c", 800, 200);
  c->Divide(3, 1);
  c->cd(1);
  DrawALICELogo(1, 0.1, 0.1, 0.9, 0.9);
  c->SaveAs("test.eps");  
}

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
    
    DrawALICELogo(kTRUE, 0.7, 0.65, 0.9, 0.85);
    
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

Float_t kEtaLimit = 1.0;
Float_t kOuterLimit = 1.59;


Bool_t FitDeltaPhi2DOneFunction(TH2* hist, TCanvas* canvas, Int_t canvasPos, Int_t graphID, Float_t x, Float_t xE, Float_t yPosChi2, Bool_t quick, Int_t histId, Int_t limits, Bool_t twoTrack)
{
  Bool_t success = kTRUE;
  Float_t etaLimit = kEtaLimit;
  Float_t outerLimit = kOuterLimit;

  if (graphID >= 15)
  {
    etaLimit = 0.5;
    outerLimit = 0.99;
  }
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
  
  Float_t mean = hist->Integral(hist->GetYaxis()->FindBin(-TMath::Pi() / 2), hist->GetYaxis()->FindBin(TMath::Pi() / 2), hist->GetYaxis()->FindBin(etaLimit), hist->GetYaxis()->FindBin(outerLimit)) / (hist->GetYaxis()->FindBin(TMath::Pi() / 2) - hist->GetYaxis()->FindBin(-TMath::Pi() / 2)) / (hist->GetYaxis()->FindBin(outerLimit) - hist->GetYaxis()->FindBin(etaLimit) + 1);
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
  func->SetParameters(hist->GetBinContent(hist->GetXaxis()->FindBin(0.0), hist->GetYaxis()->FindBin(0.0)) - mean, 0.3, 0.3, 0.25, initSigma, initSigma);
  for (Int_t i=6; i<bins+6; i++)
    func->SetParameter(i, mean);
 
 
  if (1)
  {
    // STEP 1: fit only flow using one delta eta side
    for (Int_t i=0; i<6; i++)
      func->FixParameter(i, func->GetParameter(i));
    hist->GetYaxis()->SetRangeUser(etaLimit+0.01, outerLimit-0.01);
    Int_t fitResult = hist->Fit(func, "0", "");
    Printf("Fit result: %d; Chi2/ndf: %f/%d", fitResult, func->GetChisquare(), func->GetNDF());
    if (fitResult != 0)
      success = kFALSE;

    // STEP2 : fit only Gaussian in central region
    for (Int_t i=0; i<6; i++)
      func->ReleaseParameter(i);
    //   func->SetParameters(1, 0.3, 0.3, 0.25, initSigma, initSigma);
    func->SetParLimits(0, 0, 10);
    func->SetParLimits(1, sigmaFitLimit, 0.6);
    func->SetParLimits(2, sigmaFitLimit, etaFitUpperLimit + (((graphID == 5 && histId == 4) || (graphID == 0 && histId == 0)) ? 0.1 : 0));
    func->SetParLimits(3, 0.1, 0.9);
    func->SetParLimits(4, sigmaFitLimit, 0.6);
    func->SetParLimits(5, sigmaFitLimit, etaFitUpperLimit + (((graphID == 5 && histId == 4) || (graphID == 0 && histId == 0)) ? 0.1 : 0));
    for (Int_t i=6; i<bins+6; i++)
      func->FixParameter(i, func->GetParameter(i));
    hist->GetYaxis()->SetRangeUser(-etaLimit+0.01, etaLimit-0.01);
    fitResult = hist->Fit(func, "0", "");
    Printf("Fit result: %d; Chi2/ndf: %f/%d", fitResult, func->GetChisquare(), func->GetNDF());
    if (fitResult != 0)
      success = kFALSE;

    // STEP3: fit everything, with limits
    for (Int_t i=6; i<bins+6; i++)
    {
      func->ReleaseParameter(i);
      func->SetParLimits(i, func->GetParameter(i) * 0.8, func->GetParameter(i) * 1.2);
    }
    hist->GetYaxis()->SetRangeUser(-outerLimit+0.01, outerLimit-0.01);

    // STEP4: fit everything, without limits
    for (Int_t i=6; i<bins+6; i++)
      func->SetParLimits(i, 0, 0);
    fitResult = hist->Fit(func, "0", "");
    Printf("Fit result: %d; Chi2/ndf: %f/%d", fitResult, func->GetChisquare(), func->GetNDF());
    if (fitResult != 0)
      success = kFALSE;
  }
  Int_t fitResult = hist->Fit(func, "0", "");
  Printf("Fit result: %d; Chi2/ndf: %f/%d", fitResult, func->GetChisquare(), func->GetNDF());
  if (fitResult != 0)
    success = kFALSE;

  Printf("Trying 1 Gaussian...");
  TF2* func_clone = new TF2("func_clone", DeltaPhiWidth2DFitFunction, -TMath::Pi() / 2, TMath::Pi() * 0.5, -outerLimit, outerLimit, bins+6);
  for (Int_t i=0; i<bins+6; i++)
  {
    func_clone->SetParameter(i, func->GetParameter(i));
    Double_t parmin, parmax;
    func->GetParLimits(i, parmin, parmax);
    func_clone->SetParLimits(i, parmin, parmax);
  }
  func_clone->SetParLimits(3, 1, 1);
  func_clone->FixParameter(3, 1);
  func_clone->FixParameter(4, sigmaFitLimit);
  func_clone->FixParameter(5, sigmaFitLimit);
  fitResult = hist->Fit(func_clone, "0R", "");
  Printf("Fit result: %d", fitResult);
  
  // if both parameters are within 1%, refit with 1 Gaussian only
  if (TMath::Abs(1.0 - func->GetParameter(1) / func->GetParameter(4)) < 0.01 && TMath::Abs(1.0 - func->GetParameter(2) / func->GetParameter(5)) < 0.01)
  {
    Printf("Parameters within 1%%. Using result from 1 Gaussian...");
    
    func = func_clone;
    
    if (fitResult != 0)
      success = kFALSE;
  }
  else if (func_clone->GetChisquare() * 0.99 < func->GetChisquare())
  {
    Printf("1 Gaussian fit has a at maximum 1%% (%f %f) larger chi2. Using results from 1 Gaussian...", func_clone->GetChisquare(), func->GetChisquare());
    
    func = func_clone;
    
    if (fitResult != 0)
      success = kFALSE;
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
    AddPoint(graphs[6+16][graphID], x, func->GetChisquare() / func->GetNDF(), xE, 0);
  }
  if (ndf)
  {
    Printf("#chi^{2}/ndf = %.1f/%d = %.1f", chi2, ndf, chi2 / ndf);
    DrawLatex(0.2, yPosChi2 - 0.05, 1, Form("#chi^{2}/ndf = %.1f/%d = %.1f", chi2, ndf, chi2 / ndf));
    AddPoint(graphs[7+16][graphID], x, chi2 / ndf, xE, 0);
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
/*
  Float_t etaLimitSubtract = 0;
  Float_t outerLimitSubtract = 0;
  if (graphID >= 15)
  {
    etaLimitSubtract = 0.5;
    outerLimitSubtract = 0.99;
  }
  else 
  {
    etaLimitSubtract = etaLimit;
    outerLimitSubtract = outerLimit;
  }
*/

  Float_t momentFitLimit = 0.8 - 1e-4;
  if (graphID >= 15)
    momentFitLimit = 0.4 - 1e-4;

  TH1* projx3 = hist->ProjectionX(Form("%s_projx3", hist->GetName()), hist->GetYaxis()->FindBin(-etaLimit+0.01), hist->GetYaxis()->FindBin(etaLimit-0.01));
  Float_t nBins = hist->GetYaxis()->FindBin(etaLimit-0.01) - hist->GetYaxis()->FindBin(-etaLimit+0.01)+1;
  projx3->Scale(1.0/nBins);
  
  TH1* projx3SubtractNegative = hist->ProjectionX(Form("%s_projx3SubtractPositive", hist->GetName()), hist->GetYaxis()->FindBin(-outerLimit+0.01), hist->GetYaxis()->FindBin(-etaLimit-0.01));
  nBins = hist->GetYaxis()->FindBin(-etaLimit-0.01) - hist->GetYaxis()->FindBin(-outerLimit+0.01)+1;
  projx3SubtractNegative->Scale(1.0/nBins);

  TH1* projx3SubtractPositive = hist->ProjectionX(Form("%s_projx3SubtractNegative ", hist->GetName()), hist->GetYaxis()->FindBin(etaLimit+0.01), hist->GetYaxis()->FindBin(outerLimit-0.01));
  nBins = hist->GetYaxis()->FindBin(outerLimit-0.01) - hist->GetYaxis()->FindBin(etaLimit+0.01)+1;
  projx3SubtractPositive->Scale(1.0/nBins);

  TH1* projy3 = hist->ProjectionY(Form("%s_projy3", hist->GetName()), hist->GetXaxis()->FindBin(-momentFitLimit+0.01), hist->GetXaxis()->FindBin(momentFitLimit-0.01));
  nBins = hist->GetXaxis()->FindBin(momentFitLimit-0.01) - hist->GetXaxis()->FindBin(-momentFitLimit+0.01)+1;
  projy3->Scale(1.0/nBins);

  SubtractEtaGap1D(projx3, projx3SubtractPositive, projx3SubtractNegative, projy3, etaLimit, outerLimit);
  projy3->GetXaxis()->SetRangeUser(-momentFitLimit+0.01,momentFitLimit-0.01);
  projx3->GetXaxis()->SetRangeUser(-momentFitLimit+0.01, momentFitLimit-0.01);
  SubtractEtaGap(hist, etaLimit, outerLimit, kFALSE);
//  return 0;
  canvas->cd(canvasPos++);
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


  TH1* projx2 = hist->ProjectionX(Form("%s_projx2", hist->GetName()), hist->GetYaxis()->FindBin(-etaLimit+0.01), hist->GetYaxis()->FindBin(etaLimit-0.01));
  projx2->GetXaxis()->SetRangeUser(-momentFitLimit, momentFitLimit);
  nBins = hist->GetYaxis()->FindBin(etaLimit-0.01) - hist->GetYaxis()->FindBin(-etaLimit+0.01)+1;
  projx2->Scale(1.0/nBins);
  
  TH1* projy2 = hist->ProjectionY(Form("%s_projy2", hist->GetName()), hist->GetXaxis()->FindBin(-momentFitLimit+0.01), hist->GetXaxis()->FindBin(momentFitLimit-0.01));
  projy2->GetXaxis()->SetRangeUser(-momentFitLimit, momentFitLimit);
  nBins = hist->GetXaxis()->FindBin(momentFitLimit-0.01) - hist->GetXaxis()->FindBin(-momentFitLimit+0.01)+1;
  projy2->Scale(1.0/nBins);
  
  CalculateMomentsKurtosis(momentFitLimit, projx2, 8, graphID, x, xE);
  CalculateMomentsKurtosis(momentFitLimit, projy2, 9, graphID, x, xE);

//     return success;
  
  TH1* projx1 = subtractFlow->ProjectionX(Form("%s_projx1", hist->GetName()), hist->GetYaxis()->FindBin(-0.79), hist->GetYaxis()->FindBin(0.79));
  projx1->GetXaxis()->SetRangeUser(-momentFitLimit, momentFitLimit);

  TH1* projy1 = subtractFlow->ProjectionY(Form("%s_projy1", hist->GetName()), hist->GetXaxis()->FindBin(-0.79), hist->GetXaxis()->FindBin(0.79));
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
  projx1->GetFunction("gaus")->SetLineColor(1);
  projx1->GetYaxis()->SetRangeUser(0, projx1->GetMaximum() * 1.05);

  projx2->SetLineColor(2);
  projx2->Draw("SAME");
  projx2->Fit("gaus", "I+0", "SAME");
  projx2->GetFunction("gaus")->SetLineColor(2);

  canvas->cd(canvasPos--);
  projy1->Draw();
  projy1->Fit("gaus", "I");
  projy1->GetFunction("gaus")->SetLineColor(1);
  projy1->GetYaxis()->SetRangeUser(0, projy1->GetMaximum() * 1.05);
  
  projy2->SetLineColor(2);
  projy2->Draw("SAME");
  projy2->Fit("gaus", "I+0", "SAME");
  projy2->GetFunction("gaus")->SetLineColor(2);
  
  // 1d fit (eta gap subtraction)
  AddPoint(phiWidthSummary, phiWidthSummary->GetN(), projx2->GetFunction("gaus")->GetParameter(2), 0, projx2->GetFunction("gaus")->GetParError(2));
  
  // 1d fit (lots of params)
  AddPoint(phiWidthSummary, phiWidthSummary->GetN(), projx1->GetFunction("gaus")->GetParameter(2), 0, projx1->GetFunction("gaus")->GetParError(2));

  // 1d fit (eta gap subtraction)
  AddPoint(etaWidthSummary, etaWidthSummary->GetN(), projy2->GetFunction("gaus")->GetParameter(2), 0, projy2->GetFunction("gaus")->GetParError(2));
  
  // 1d fit (lots of params)
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

  Float_t etaFitLimit = outerLimit;

  Bool_t oneGaussian = kTRUE;

  if (!oneGaussian) cout << endl << "Fitting two 1D Gaussians" << endl;
  else cout << endl << "Fitting one 1D Gaussians" << endl;

  TF1* func4phi = new TF1("func4phi", "[0]+[1]*([4]/TMath::Sqrt(TMath::TwoPi())/[2]*exp(-0.5*((x/[2])**2))+(1-[4])/TMath::Sqrt(TMath::TwoPi())/[3]*exp(-0.5*((x/[3])**2)))", -0.5 * TMath::Pi()+0.01, 0.5 * TMath::Pi()-0.01);
  TF1* func4eta = new TF1("func4eta", "[0]+[1]*([4]/TMath::Sqrt(TMath::TwoPi())/[2]*exp(-0.5*((x/[2])**2))+(1-[4])/TMath::Sqrt(TMath::TwoPi())/[3]*exp(-0.5*((x/[3])**2)))", -etaFitLimit+0.01, etaFitLimit-0.01);

  func4phi->SetParLimits(1, 0, 10);
  func4phi->SetParLimits(3, sigmaFitLimit, 0.7);
  func4phi->FixParameter(0, 0);
  if (!oneGaussian)
  {
    func4phi->SetParLimits(4, 0, 1);
    func4phi->SetParLimits(2, sigmaFitLimit, 0.7);
  }
  else
  {
    func4phi->FixParameter(4, 0);
    func4phi->FixParameter(2, 0);
  }

  func4eta->SetParLimits(1, 0, 10);
  func4eta->SetParLimits(3, sigmaFitLimit, etaFitUpperLimit);
  func4eta->FixParameter(0, 0);
  if (!oneGaussian)
  {
    func4eta->SetParLimits(4, 0, 1);
    func4eta->SetParLimits(2, sigmaFitLimit, etaFitUpperLimit);
  }
  else
  {
    func4eta->FixParameter(4, 0);
    func4eta->FixParameter(2, 0);
  }
  canvas->cd(canvasPos++);
  func4phi->SetLineColor(4);
  TFitResultPtr fitResultPtr = projx3->Fit(func4phi, "SIR", "SAME");
  TMatrixDSym cov = fitResultPtr->GetCovarianceMatrix();

  projx3->Draw("SAME");
  projx3->SetLineColor(4);

  if (func4phi->GetNDF() > 0)
  {
    printf("#chi^{2}/ndf = %.1f/%d = %.1f  ", func4phi->GetChisquare(), func4phi->GetNDF(), func4phi->GetChisquare() / func4phi->GetNDF());
    AddPoint(graphs[43][graphID], x, func4phi->GetChisquare() / func4phi->GetNDF(), xE, 0);
  }
  chi2 = 0;
  ndf = 0;

  for (Int_t i=projx3->FindBin(-momentFitLimit); i<=projx3->FindBin(momentFitLimit);i++)
  {
//    Float_t temp = chi2;
    if (projx3->GetBinError(i) > 0)
    {
      chi2 += TMath::Power((projx3->GetBinContent(i)-func4phi->Eval(projx3->GetBinCenter(i)))/projx3->GetBinError(i),2);
      //    cerr << i << "\t" << projx3->GetBinCenter(i) << "\t" << chi2-temp << endl;
      ndf++;
    }
  }
  if (!oneGaussian) ndf = ndf - 4;
  else ndf = ndf - 2;
  printf("Calculated #chi^{2}/ndf = %.1f/%d = %.1f  ", chi2, ndf, chi2/ndf);

  first = 2;
  second = 3;
  if (func4phi->GetParameter(2) < func4phi->GetParameter(3))
  {
    first = 3;
    second = 2;
  }
  
  if (!oneGaussian)
  {
    Float_t vector[3];
    for (Int_t i=2; i<5;i++)
      vector[i-2] = cov[i][2]*func4phi->GetParameter(4) + cov[i][3]*(1-func4phi->GetParameter(4))+cov[i][4]*(TMath::Abs(func4phi->GetParameter(2))-TMath::Abs(func4phi->GetParameter(3)));

    Float_t sigma = TMath::Sqrt(vector[0]*func4phi->GetParameter(4) + vector[1]*(1-func4phi->GetParameter(4))+ vector[2]*(TMath::Abs(func4phi->GetParameter(2))-TMath::Abs(func4phi->GetParameter(3))));
    Float_t rms = TMath::Abs(func4phi->GetParameter(2))*func4phi->GetParameter(4)+TMath::Abs(func4phi->GetParameter(3))*(1-func4phi->GetParameter(4));
    AddPoint(graphs[37][graphID], x, rms, xE, sigma);
  }
  else AddPoint(graphs[37][graphID], x, TMath::Abs(func4phi->GetParameter(first)), xE, func4phi->GetParError(first));
  AddPoint(graphs[38][graphID], x, 0, xE, 0);
  AddPoint(graphs[36][graphID], x, 1, xE, 0);

  AddPoint(graphs[35][graphID], x, func4phi->GetParameter(1), xE, func4phi->GetParError(1));

  AddPoint(phiWidthSummary, phiWidthSummary->GetN(), TMath::Abs(func4phi->GetParameter(first)), 0, func4phi->GetParError(first));
  AddPoint(phiWidthSummary, phiWidthSummary->GetN(), TMath::Abs(func4phi->GetParameter(second)), 0, func4phi->GetParError(second));

  Float_t scale = func4phi->GetParameter(1);
  Float_t weighting = func4phi->GetParameter(4);
  func4phi->SetParameter(1, scale * weighting);
  func4phi->SetParameter(4, 1);
  func4phi->SetLineStyle(2);
  func4phi->DrawCopy("SAME");

  func4phi->SetParameter(1, scale * (1.0-weighting));
  func4phi->SetParameter(4, 0);
  func4phi->SetLineStyle(3);
  func4phi->DrawCopy("SAME");

  fitResultPtr->Print("V");
  if (fitResultPtr != 0)
    success = kFALSE;

  canvas->cd(canvasPos++);
  func4eta->SetLineColor(4);
  fitResultPtr = projy3->Fit(func4eta, "SIR", "SAME");
  cov = fitResultPtr->GetCovarianceMatrix();

  projy3->Draw("SAME");
  projy3->SetLineColor(4);

  chi2 = 0;
  ndf = 0;
  if (func4eta->GetNDF() > 0)
  {
    printf("#chi^{2}/ndf = %.1f/%d = %.1f  ", func4eta->GetChisquare(), func4eta->GetNDF(), func4eta->GetChisquare() / func4eta->GetNDF());
    AddPoint(graphs[44][graphID], x, func4eta->GetChisquare() / func4eta->GetNDF(), xE, 0);
  }
  for (Int_t i=projy3->FindBin(-momentFitLimit); i<=projy3->FindBin(momentFitLimit);i++)
  {
//    Float_t temp = chi2;
    chi2 += TMath::Power((projy3->GetBinContent(i)-func4eta->Eval(projy3->GetBinCenter(i)))/projy3->GetBinError(i),2);
//    cerr << i << "\t" << projy3->GetBinCenter(i) << "\t" << chi2-temp << endl;
    ndf++;
  }
  if (!oneGaussian) ndf = ndf - 4;
  else ndf = ndf - 2;
  printf("Calculated #chi^{2}/ndf = %.1f/%d = %.1f  ", chi2, ndf, chi2/ndf);

  first = 2;
  second = 3;
  if (func4eta->GetParameter(2) < func4eta->GetParameter(3))
  {
    first = 3;
    second = 2;
  }

  AddPoint(graphs[39][graphID], x, func4eta->GetParameter(1), xE, func4eta->GetParError(1));
  if (!oneGaussian)
  {
    Float_t vector[3];
    for (Int_t i=2; i<5;i++)
      vector[i-2] = cov[i][2]*func4eta->GetParameter(4) + cov[i][3]*(1-func4eta->GetParameter(4))+cov[i][4]*(TMath::Abs(func4eta->GetParameter(2))-TMath::Abs(func4eta->GetParameter(3)));

    Float_t sigma = TMath::Sqrt(vector[0]*func4eta->GetParameter(4) + vector[1]*(1-func4eta->GetParameter(4))+ vector[2]*(TMath::Abs(func4eta->GetParameter(2))-TMath::Abs(func4eta->GetParameter(3))));
    Float_t rms = TMath::Abs(func4eta->GetParameter(2))*func4eta->GetParameter(4)+TMath::Abs(func4eta->GetParameter(3))*(1-func4eta->GetParameter(4));
    AddPoint(graphs[41][graphID], x, rms, xE, sigma);
  }
  else AddPoint(graphs[41][graphID], x, TMath::Abs(func4eta->GetParameter(first)), xE, func4eta->GetParError(first));
  AddPoint(graphs[42][graphID], x, 0, xE, 0);
  AddPoint(graphs[40][graphID], x, 1, xE, 0);

  AddPoint(etaWidthSummary, etaWidthSummary->GetN(), TMath::Abs(func4eta->GetParameter(first)), 0, func4eta->GetParError(first));
  AddPoint(etaWidthSummary, etaWidthSummary->GetN(), TMath::Abs(func4eta->GetParameter(second)), 0, func4eta->GetParError(second));

  scale = func4eta->GetParameter(1);
  weighting = func4eta->GetParameter(4);
  func4eta->SetParameter(1, scale * weighting);
  func4eta->SetParameter(4, 1);
  func4eta->SetLineStyle(2);
  func4eta->DrawCopy("SAME");

  func4eta->SetParameter(1, scale * (1.0-weighting));
  func4eta->SetParameter(4, 0);
  func4eta->SetLineStyle(3);
  func4eta->DrawCopy("SAME");
  
  fitResultPtr->Print("V");
  if (fitResultPtr != 0)
    success = kFALSE;

//  return 0;

  // 2d fit with two gaussians
  canvas->cd(canvasPos++);
//   TF2* func3 = new TF2("func3", "[0]+[1]*exp(-0.5*((x/[2])**2+(y/[3])**2))+[4]*exp(-0.5*((x/[5])**2+(y/[6])**2))", -0.5 * TMath::Pi(), 0.5 * TMath::Pi(), -outerLimit, outerLimit);
//   func3->SetParameters(0, 1, 0.3, 0.3, 1, 0.6, 0.6);
//   func3->SetParLimits(4, 0, 10);
//   Float_t etaFitLimit = 0.5;
  TF2* func3 = new TF2("func3", "[0]+[1]*([4]/TMath::TwoPi()/[2]/[3]*exp(-0.5*((x/[2])**2+(y/[3])**2))+(1-[4])/TMath::TwoPi()/[5]/[6]*exp(-0.5*((x/[5])**2+(y/[6])**2)))", -0.5 * TMath::Pi(), 0.5 * TMath::Pi(), -etaFitLimit, etaFitLimit);

  func3->SetParLimits(4, 0.1, 0.9);
  func3->SetParLimits(1, 0, 10);
  func3->SetParLimits(2, sigmaFitLimit, 0.7);
  func3->SetParLimits(3, sigmaFitLimit, etaFitUpperLimit);
  func3->SetParLimits(5, sigmaFitLimit, 0.7);
  func3->SetParLimits(6, sigmaFitLimit, etaFitUpperLimit);
  
  func3->FixParameter(0, 0);
/*
  sigmaFitLimit = 0.09;
  etaFitUpperLimit = 0.12;
  func3->SetParLimits(4, 0.1, 0.9);
  func3->SetParLimits(1, 0, 10);
  func3->SetParLimits(2, sigmaFitLimit, 0.12);
  func3->SetParLimits(3, sigmaFitLimit, etaFitUpperLimit);
  func3->SetParLimits(5, sigmaFitLimit, 0.12);
  func3->SetParLimits(6, sigmaFitLimit, etaFitUpperLimit);
  func3->FixParameter(0, 0);
*/
  func3->SetParameters(0, hist->GetBinContent(hist->GetXaxis()->FindBin(0.0), hist->GetYaxis()->FindBin(0.0)) - mean, 0.3, 0.3, 0.25, initSigma, initSigma);
//   for (Int_t i=0; i<6; i++)
//     func3->SetParameter(i+1, func->GetParameter(i));

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
  if (fitResult != 0)
    success = kFALSE;
  
  Printf("Testing 1 Gaussian...");

  TF2* func3_clone = (TF2*) func3->Clone("func3_clone");
  
  func3_clone->SetParLimits(4, 1, 1);
  func3_clone->FixParameter(4, 1);
  func3_clone->FixParameter(5, sigmaFitLimit);
  func3_clone->FixParameter(6, sigmaFitLimit);
  
  fitResult = hist->Fit(func3_clone, "0R", "");
  Printf("Fit result: %d", fitResult);
  
  // if both parameters were within 1%, use 1 Gaussian parameters
  if (TMath::Abs(1.0 - func3->GetParameter(2) / func3->GetParameter(5)) < 0.01 && TMath::Abs(1.0 - func3->GetParameter(3) / func3->GetParameter(6)) < 0.01)
  {
    Printf("Parameters within 1%%. Using results from 1 Gaussian...");
    
    if (fitResult != 0)
      success = kFALSE;
    func3 = func3_clone;
  }
  else if (func3_clone->GetChisquare() * 0.99 < func3->GetChisquare())
  {
    Printf("1 Gaussian fit has a at maximum 1%% (%f %f) larger chi2. Using results from 1 Gaussian...", func3_clone->GetChisquare(), func3->GetChisquare());
    
    if (fitResult != 0)
      success = kFALSE;
    func3 = func3_clone;
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
  
  Printf("Finished with %d", success);

  return success;
}

void AnalyzeDeltaPhiEtaGap2D(const char* fileName, const char* outputFileName = "graphs.root")
{
  gROOT->SetBatch(kTRUE);
  if (!gROOT->IsBatch())
  {
    Printf("Not in batch mode. Exiting!");
    return;
  }
  
  CreateGraphStructure();

  TFile::Open(fileName);
  TList checkList;
  
  Int_t maxLeadingPt = 4;
  Int_t maxAssocPt = 6;


  for (Int_t i=0; i<maxLeadingPt; i++)
  {
    for (Int_t j=1; j<maxAssocPt; j++)
    {
      // only process when first is filled
      TH2* hist1 = (TH2*) gFile->Get(Form("dphi_%d_%d_%d", i, j+1, 0));
      if (!hist1)
      {
        cout << "Hist1 does not exist." << endl;
	continue;
      }
//      if (hist1->GetEntries() < 1e4)
      if (hist1->GetEntries() < 10)
      {
	Printf("%d %d Only %f entries. Skipping...", i, j, hist1->GetEntries());
	continue;
      }
      
      for (Int_t histId = 0; histId < NHists; histId++)
      {
// 	if (i != 1 || j != 2)
// 	  continue;
    
	hist1 = (TH2*) gFile->Get(Form("dphi_%d_%d_%d", i, j+1, histId));
	if (!hist1)
        {
          cout << "Hist1 does not exist." << endl;
          continue;
        }
	
//	if (hist1->GetEntries() < 1e4)
	if (hist1->GetEntries() < 10)
	{
	  Printf("%d %d %d Only %f entries. Skipping...", i, j, histId, hist1->GetEntries());
	  continue;
	}
	
	TCanvas* canvas = new TCanvas(Form("DeltaPhi_%d_%d_%d_%d", histId, i, j, 1), Form("DeltaPhi_%d_%d_%d_%d", histId, i, j, 1), 1400, 1100);
	canvas->Divide(5, 3);
	
	for (Int_t k=1; k<=3; k++)
	{
	  canvas->cd(3 * j + k);
	  gPad->SetLeftMargin(0.15);
	  gPad->SetBottomMargin(0.2);
	  gPad->SetTopMargin(0.01);
	  gPad->SetRightMargin(0.01);
	}
	
	Int_t graphID = i * (maxAssocPt - 1) + j - 1;

	Printf("\n\n>>> %d %d %d %d", i, j, histId, graphID);
	
	if (histId == 0)
	  for (Int_t k=0; k<NGraphs; k++)
	    graphs[k][graphID]->SetTitle(hist1->GetTitle());
	
	Float_t centralityAxisMapping[] = { 5, 65, 100, 30, 15, 50 };
	Float_t centralityAxisMappingE[] = { 5, 5, 0, 10, 5, 10 };

	Bool_t success = FitDeltaPhi2DOneFunction((TH2*) hist1, canvas, 1, graphID, centralityAxisMapping[histId], centralityAxisMappingE[histId], 0.9, kTRUE, histId, (i > 0) ? 1 : 0, (j <= 2 && histId != 2));
	if (!success)
	  checkList.Add(new TObjString(Form("AnalyzeDeltaPhiEtaGap2DExample(\"%s\", %d, %d, %d)", fileName, i, j, histId)));
	
// 	canvas->SaveAs(Form("%s.png", canvas->GetName()));
// 	delete canvas;
// 	break;
// return;
      }
      
//       break;
    }
    
//     break;
  }
  
  WriteGraphs(outputFileName);

  for (Int_t i=0; i<checkList.GetEntries(); i++)
    Printf("%s", checkList.At(i)->GetName());
}

void AnalyzeDeltaPhiEtaGap2DExample(const char* fileName, Int_t i, Int_t j, Int_t histId, Bool_t drawDetails = kFALSE)
{
  CreateGraphStructure();

  TFile::Open(fileName);
  
  TH2* hist1 = (TH2*) gFile->Get(Form("dphi_%d_%d_%d", i, j+1, histId));
  if (!hist1)
    return;
  
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

//  if (hist1->GetEntries() < 1e4)
  if (hist1->GetEntries() < 10)
    return;
  
  // NOTE fix normalization. these 2d correlations come out of AliUEHist normalized by dphi bin width, but not deta
  if (drawDetails)
    hist1->Scale(1.0 / hist1->GetYaxis()->GetBinWidth(1));

  TCanvas* canvas = new TCanvas(Form("DeltaPhi_%d_%d_%d_%d", histId, i, j, 1), Form("DeltaPhi_%d_%d_%d_%d", histId, i, j, 1), 1400, 1100);
  canvas->Divide(5, 3);
  
  Int_t graphID = i * (6 - 1) + j - 1;
  FitDeltaPhi2DOneFunction((TH2*) hist1, canvas, 1, graphID, 0, 0, 0.9, kFALSE, histId, (i > 0) ? 1 : 0, (j <= 2 && histId != 2));
  
  if (!drawDetails)
    return;
  
  TVirtualPad* pad = canvas->cd(6);
  TCanvas* c = new TCanvas("c", "c", 1500, 1000);
  c->Divide(3, 2);
  c->cd(1);
  pad->SetPad(0, 0, 1, 1);
  pad->SetLeftMargin(0.15);
  pad->Draw();
  pad->cd();

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
  paveText->Draw();
  
  c->cd(1);
  DrawALICELogo(kTRUE, 0.2, 0.7, 0.4, 0.9);
  
//   return;
    
  //   c->SaveAs("fit_subtracted.eps");

  pad = canvas->cd(12);
  c->cd(2);
  pad->SetPad(0, 0, 1, 1);
  pad->SetLeftMargin(0.15);
  pad->Draw();
//   c->SaveAs("fit_fit1.eps");
  
  pad = canvas->cd(13);
  c->cd(3);
  gPad->SetLeftMargin(0.15);
  pad->SetPad(0, 0, 1, 1);
  pad->SetLeftMargin(0.15);
  pad->Draw();
  c->SaveAs("fit_residual1.eps");
  
  TCanvas* c2 = new TCanvas("c3b", "c3b", 800, 800);
  gPad->SetLeftMargin(0.15);
  gPad->SetGridy();
  hist1 = (TH2*) pad->GetListOfPrimitives()->First();
  for (i=0; i<10; i++)
  {
    TH1* proj = hist1->ProjectionX(Form("p_%d", i), 11+i*2, 12+i*2);
    proj->Scale(0.5);
    
    proj->Add(new TF1("func", "1", -10, 10), -0.5 + 0.1 * i);
    proj->SetStats(0);
    proj->GetXaxis()->SetTitleOffset(1.0);
    proj->GetYaxis()->SetTitleOffset(1.4);
    proj->GetYaxis()->SetNdivisions(512);
    proj->SetYTitle("Residuals (projected / shifted)");
    proj->GetYaxis()->SetRangeUser(-0.59, 0.49);
    proj->Draw((i == 0) ? "" : "SAME");    
  }
  c2->SaveAs("fit_residual1_proj.eps");

  pad = canvas->cd(4);
  c->cd(4);
  pad->SetPad(0, 0, 1, 1);
  pad->SetLeftMargin(0.15);
  pad->Draw();
//   c->SaveAs("fit_fit2.eps");

  pad = canvas->cd(3);
  c->cd(5);
  pad->SetPad(0, 0, 1, 1);
  pad->SetLeftMargin(0.15);
  pad->Draw();
//   c->SaveAs("fit_residual2.eps");

  c2 = new TCanvas("c5b", "c5b", 800, 800);
  gPad->SetLeftMargin(0.15);
  gPad->SetGridy();
  hist1 = (TH2*) pad->GetListOfPrimitives()->First();
  for (i=0; i<10; i++)
  {
    TH1* proj = hist1->ProjectionX(Form("p2_%d", i), 11+i*2, 12+i*2);
    proj->Scale(0.5);
    
    proj->Add(new TF1("func", "1", -10, 10), -0.5 + 0.1 * i);
    proj->SetStats(0);
    proj->GetXaxis()->SetTitleOffset(1.0);
    proj->GetYaxis()->SetTitleOffset(1.4);
    proj->GetYaxis()->SetNdivisions(512);
    proj->SetYTitle("Residuals (projected / shifted)");
    proj->GetYaxis()->SetRangeUser(-0.59, 0.49);
    proj->Draw((i == 0) ? "" : "SAME");    
  }
  c2->SaveAs("fit_residual2_proj.eps");

  pad = canvas->cd(7);
  c->cd(6);
  pad->SetPad(0, 0, 1, 1);
  pad->SetLeftMargin(0.15);
  pad->Draw();
  
  c->SaveAs("fit_example.eps");
  c->SaveAs("fit_example.png");
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
      if (graph[histId]->GetY()[i] > 0 && graph2[histId]->GetY()[i] > 0)
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
const char* labels[6] = { "0-10%", "60-70%", "pp", "20-40%", "10-20%", "40-60%" };

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
//  return (i == 1 || i == 7 || i == 9 || i == 13 || i == 14);
  return (i == 7 || i == 9 || i == 13 || i == 14); //For the STAR comparison
}

Int_t skipGraphList[] =  { -1, -1, -1 };
Bool_t SkipGraphForThisPlot(Int_t graphId)
{
  for (Int_t i=0; i<3; i++)
    if (skipGraphList[i] == graphId)
      return kTRUE;
  return kFALSE;
}

TGraphAsymmErrors* FixGraph(TGraphErrors* graph, Float_t shift)
{
  graph->Sort();
//   graph->Print();
  
  if (graph->GetN() > 4 && graph->GetX()[4] == 75)
  {
    graph->GetX()[4] = 65;
    graph->GetEX()[4] = 5;
  }
  
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
    }
  }
}

void PrepareGraphs(Int_t nHists, TGraphErrors** graph, TGraphErrors** systematicA, TGraphErrors** systematicB, TMultiGraph** multiGraph, TMultiGraph** multiGraphSyst, Int_t uncertaintyID, Float_t offset = 0)
{
  Int_t colors[16] =  { 1, 3, 2, 6, 4, 7, 8, 9, 11, 12, 28, 30, 36, 40, 46 };
  Int_t markers[16] = { 20, 21, 22, 23, 24, 25, 26, 27, 28, 30, 31, 32, 33, 34, 2, 5};
//  Int_t fillStyle[11] = { 3003, 3003, 3003, 3003, 3003, 3003, 3003, 3003, 3003, 3003, 3003 };
  Int_t fillStyle[11] = { 3001, 3002, 3003, 3004, 3005, 3006, 3007, 3008, 3009, 3010, 3011 };
  Int_t count = 0;
  
  if (*multiGraph == 0)
    *multiGraph = new TMultiGraph;
  if (*multiGraphSyst == 0)
    *multiGraphSyst = new TMultiGraph;
  
  Float_t shift = -2 + offset;
  for (Int_t i=0; i<nHists; i++)
  {
    if (SkipGraph(i))
      continue;
    
    if (SkipGraphForThisPlot(i))
    {
      count++;
      continue;
    }
    
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
	
	if (uncertaintyID == 1)
	{
	  // 5%
	  yMin *= 0.95;
	  yMax *= 1.05;
	}
	else if (uncertaintyID == 2)
	{
	  yMin -= 0.20;
	  yMax += 0.20;
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
    shift += 1;
   
    count++;
  }
}


Bool_t drawLogo = kFALSE;
const char* MCLabel = 0;

void DrawCentrality(const char* canvasName, Int_t nHists, TGraphErrors** graph, Float_t min = 0, Float_t max = 0, const char* yLabel = "", TGraphErrors** systematicA = 0, TGraphErrors** systematicB = 0, TGraphErrors** graph2 = 0, TGraphErrors** systematic2 = 0, TGraphErrors** graph3 = 0, Int_t uncertaintyID = -1)
{
  Bool_t found = kTRUE;
  TCanvas* c1 = 0;//(TCanvas*) gROOT->GetListOfCanvases()->FindObject(canvasName);
  if (!c1)
  {
    c1 = new TCanvas(canvasName, canvasName, 800, 600);
    found = kFALSE;
  }
  c1->cd();
  
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
  
  gPad->SetGridx();
  gPad->SetGridy();
  
  DrawLatex(0.13, 0.85,  1, "Pb-Pb #sqrt{s_{NN}} = 2.76 TeV", 0.03);
  DrawLatex(0.13, 0.81,  1, "pp #sqrt{s} = 2.76 TeV", 0.03);
  DrawLatex(0.13, 0.77, 1, "|#eta| < 0.9", 0.03);
  TString text;
  TString text2;
  if (TString(canvasName).BeginsWith("sigma_phi") || TString(canvasName).BeginsWith("kurtosisphi"))
  {
    text = "Projected within |#Delta#eta| < 0.80";
    text2 = "Calculated within |#Delta#varphi| < 0.87";
  }
  if (TString(canvasName).BeginsWith("sigma_eta") || TString(canvasName).BeginsWith("kurtosiseta"))
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
  
  if (drawLogo)
    DrawALICELogo(kTRUE, 0.41, 0.75, 0.54, 0.95);
  
  gSystem->mkdir(fgFolder, kTRUE);
  c1->SaveAs(Form("%s/%s.gif", fgFolder.Data(), canvasName));
  c1->SaveAs(Form("%s/%s.eps", fgFolder.Data(), canvasName));
}

void CalculateRMSSigma(TGraphErrors*** graphsTmp = 0)
{
  if (!graphsTmp)
    graphsTmp = graphs;
  
  CalculateRMS(NHists, graphsTmp[1], graphsTmp[4], graphsTmp[3]);
  CalculateRMS(NHists, graphsTmp[2], graphsTmp[5], graphsTmp[3]);

  CalculateRMS(NHists, graphsTmp[1+16], graphsTmp[4+16], graphsTmp[3+16]);
  CalculateRMS(NHists, graphsTmp[2+16], graphsTmp[5+16], graphsTmp[3+16]);

  // sqrt(moment2) = sigma
  SqrtAll2(NHists, graphsTmp[8], graphsTmp[8+16]);
  SqrtAll2(NHists, graphsTmp[9], graphsTmp[9+16]);
}

// void CombineSyst(const char* systFile)
// {
//   TGraphErrors*** current = graphs;
//   
//   ReadGraphs(systFile);
//   CalculateRMSSigma();
//   TGraphErrors*** graphsSyst = graphs;
//   graphs = current;
//   
//   const Int_t NGraphList = 6;
//   Int_t graphList[] = { 1, 2, 8, 9, 14, 15 };
//   for (Int_t i=0; i<NGraphList; i++)
//   {
//     for (Int_t j=0; j<NHists; j++)
//     {
//       
//     graphs[i+16][j]
//     graphsSyst[i][j]
//     }
//   }
//   
// 
// }

void DrawResultsCentrality(const char* fileName = "graphs.root", const char* fileNameWingRemoved = 0, Int_t offset = 0)
{
  TGraphErrors*** graphsWingRemoved = 0;
  if (fileNameWingRemoved)
  {
    ReadGraphs(fileNameWingRemoved);
    graphsWingRemoved = graphs;
  }
  
  ReadGraphs(fileName);
  
  Int_t nHists = 24; //NHists;
  
  if (1)
  {
/*    DrawCentrality("norm", nHists, graphs[0+offset], 0, 0.2, "N (a.u.)", (graphsWingRemoved) ? graphsWingRemoved[0+offset] : 0);
    
//     return;
    
    DrawCentrality("width_phi1_centrality", nHists, graphs[1+offset], 0, 0.8, "#sigma_{#Delta#varphi, 1} (rad.)", (graphsWingRemoved) ? graphsWingRemoved[1+offset] : 0);
//     return;
    DrawCentrality("width_phi2_centrality", nHists, graphs[4+offset], 0, 0.8, "#sigma_{#Delta#varphi, 2} (rad.)", (graphsWingRemoved) ? graphsWingRemoved[4+offset] : 0);
    DrawCentrality("width_eta1_centrality", nHists, graphs[2+offset], 0, 0.8, "#sigma_{#Delta#eta, 1} (rad.)", (graphsWingRemoved) ? graphsWingRemoved[2+offset] : 0);
    DrawCentrality("width_eta2_centrality", nHists, graphs[5+offset], 0, 0.8, "#sigma_{#Delta#eta, 2} (rad.)", (graphsWingRemoved) ? graphsWingRemoved[5+offset] : 0);
*/    
/*    DrawCentrality("norm_phi(1D)", nHists, graphs[35], 0, 1, "N (a.u.)", (graphsWingRemoved) ? graphsWingRemoved[0+offset] : 0);
    DrawCentrality("norm_eta(1D)", nHists, graphs[39], 0, 1, "N (a.u.)", (graphsWingRemoved) ? graphsWingRemoved[0+offset] : 0);
*/    DrawCentrality("width_phi1_centrality(1D)", nHists, graphs[37], 0, 0.8, "#sigma_{#Delta#varphi, 1} (rad.)", (graphsWingRemoved) ? graphsWingRemoved[1+offset] : 0);
//    DrawCentrality("width_phi2_centrality(1D)", nHists, graphs[38], 0, 0.8, "#sigma_{#Delta#varphi, 2} (rad.)", (graphsWingRemoved) ? graphsWingRemoved[4+offset] : 0);
    DrawCentrality("width_eta1_centrality(1D)", nHists, graphs[41], 0, 0.8, "#sigma_{#Delta#eta, 1} (rad.)", (graphsWingRemoved) ? graphsWingRemoved[2+offset] : 0);
//    DrawCentrality("width_eta2_centrality(1D)", nHists, graphs[42], 0, 0.8, "#sigma_{#Delta#eta, 2} (rad.)", (graphsWingRemoved) ? graphsWingRemoved[5+offset] : 0);
//    DrawCentrality("norm2phi(1D)", nHists, graphs[36], 0, 1.5, "N (a.u.)", (graphsWingRemoved) ? graphsWingRemoved[0+offset] : 0);
//    DrawCentrality("norm2eta(1D)", nHists, graphs[40], 0, 1.5, "N (a.u.)", (graphsWingRemoved) ? graphsWingRemoved[0+offset] : 0);

    CalculateRMSSigma();
    if (graphsWingRemoved)
      CalculateRMSSigma(graphsWingRemoved);
/*
    DrawCentrality("phi_rms_centrality", nHists, graphs[1+offset], 0, 0.8, Form("#sigma_{#Delta#varphi} (%s) (rad.)", fitLabel), (graphsWingRemoved) ? graphsWingRemoved[1+offset] : 0, 0, 0, 0, 0, 1);

    DrawCentrality("eta_rms_centrality", nHists, graphs[2+offset], 0, 0.8 + 0.4 / 16 * offset, Form("#sigma_{#Delta#eta} (%s)", fitLabel), (graphsWingRemoved) ? graphsWingRemoved[2+offset] : 0, 0, 0, 0, 0, 1);
*/
//     return;
  DrawCentrality("phi_rms_centrality(1D)", nHists, graphs[37], 0, 0.8, Form("#sigma_{#Delta#varphi} (%s) (rad.)", fitLabel), (graphsWingRemoved) ? graphsWingRemoved[1+offset] : 0, 0, 0, 0, 0, 1);

  DrawCentrality("eta_rms_centrality(1D)", nHists, graphs[41], 0, 0.8 + 0.4 / 16 * offset, Form("#sigma_{#Delta#eta} (%s)", fitLabel), (graphsWingRemoved) ? graphsWingRemoved[2+offset] : 0, 0, 0, 0, 0, 1);
  DrawCentrality("chi2_phi(1D)", nHists, graphs[43], 0, 150, "#chi^{2}/ndf");
  DrawCentrality("chi2_eta(1D)", nHists, graphs[44], 0, 150, "#chi^{2}/ndf");
/*    
    DrawCentrality("chi2_1", nHists, graphs[6+offset], 0.5, 5, "#chi^{2}/ndf (full region)");
    DrawCentrality("chi2_2", nHists, graphs[7+offset], 0.5, 5, "#chi^{2}/ndf (peak region)");
    
    DrawCentrality("sigma_phi", nHists, graphs[8+offset], 0, 1, "#sigma_{#Delta#varphi} (rad.)", (graphsWingRemoved) ? graphsWingRemoved[8+offset] : 0, 0, 0, 0, 0, 1);
    DrawCentrality("sigma_eta", nHists, graphs[9+offset], 0, 1, "#sigma_{#Delta#eta}", (graphsWingRemoved) ? graphsWingRemoved[9+offset] : 0, 0, 0, 0, 0, 1);

//     DrawCentrality("moment3_phi", nHists, graphs[10+offset], -0.5, 0.5, "moment3 #varphi (rad.)", (graphsWingRemoved) ? graphsWingRemoved[10+offset] : 0, 0, 0, 0, 0, 0);
//     DrawCentrality("moment3_eta", nHists, graphs[11+offset], -0.5, 0.5, "moment3 #eta", (graphsWingRemoved) ? graphsWingRemoved[11+offset] : 0, 0, 0, 0, 0, 0);
*/
/*    DrawCentrality("IAAFit", nHists, graphs[32+offset], 0, 20, "I_AA", (graphsWingRemoved) ? graphsWingRemoved[32+offset] : 0);
    DrawCentrality("Yield", nHists, graphs[33+offset], 0, 20, "Yield", (graphsWingRemoved) ? graphsWingRemoved[33+offset] : 0);
    DrawCentrality("IAAHist", nHists, graphs[34+offset], 0, 20, "I_AA", (graphsWingRemoved) ? graphsWingRemoved[34+offset] : 0);
    DrawCentrality("IAA", nHists, graphs[34+offset], 0, 20, "I_AA", graphs[32+offset]);
*/
  }
  
//  DrawCentrality("kurtosisphi_centrality", 12, graphs[14+offset], -1.5, 2.5, "Kurtosis #Delta#varphi", (graphsWingRemoved) ? graphsWingRemoved[14+offset] : 0, 0, 0, 0, 0, 2);
//  DrawCentrality("kurtosiseta_centrality", 12, graphs[15+offset], -1.5, 2.5, "Kurtosis #Delta#eta", (graphsWingRemoved) ? graphsWingRemoved[15+offset] : 0, 0, 0, 0, 0, 2);
}

/*
void DrawResultsCentrality(const char* fileName = "graphs.root", const char* fileNameWingRemoved = 0, Int_t offset = 0)
{
  TGraphErrors*** graphsWingRemoved = 0;
  if (fileNameWingRemoved)
  {
    ReadGraphs(fileNameWingRemoved);
    graphsWingRemoved = graphs;
  }
  
  ReadGraphs(fileName);
  
  Int_t nHists = 13; //NHists;
  
  if (1)
  {
    DrawCentrality("norm", nHists, graphs[0], 0, 0.2, "N (a.u.)", graphs[0+16], (graphsWingRemoved) ? graphsWingRemoved[0] : 0);
    
//     return;
    
    DrawCentrality("width_phi1_centrality", nHists, graphs[1], 0, 0.8, "#sigma_{#varphi, 1} (rad.)", graphs[1+16], (graphsWingRemoved) ? graphsWingRemoved[1] : 0);
//     return;
    DrawCentrality("width_phi2_centrality", nHists, graphs[4], 0, 0.8, "#sigma_{#varphi, 2} (rad.)", graphs[4+16], (graphsWingRemoved) ? graphsWingRemoved[4] : 0);
    DrawCentrality("width_eta1_centrality", nHists, graphs[2], 0, 0.8, "#sigma_{#eta, 1} (rad.)", graphs[2+16], (graphsWingRemoved) ? graphsWingRemoved[2] : 0);
    DrawCentrality("width_eta2_centrality", nHists, graphs[5], 0, 0.8, "#sigma_{#eta, 2} (rad.)", graphs[5+16], (graphsWingRemoved) ? graphsWingRemoved[5] : 0);
    
    CalculateRMSSigma();
    if (graphsWingRemoved)
      CalculateRMSSigma(graphsWingRemoved);

    DrawCentrality("phi_rms_centrality", nHists, graphs[1], 0, 0.8, "#sigma_{#varphi} (fit) (rad.)", graphs[1+16], (graphsWingRemoved) ? graphsWingRemoved[1] : 0, 0, 0, 0, 1);

    DrawCentrality("eta_rms_centrality", nHists, graphs[2], 0, 0.8, "#sigma_{#eta} (fit)", graphs[2+16], (graphsWingRemoved) ? graphsWingRemoved[2] : 0, 0, 0, 0, 1);

    DrawCentrality("chi2_1", nHists, graphs[6], 0.5, 2.5, "#chi^{2}/ndf (full region)");
    DrawCentrality("chi2_2", nHists, graphs[7], 0.5, 2.5, "#chi^{2}/ndf (peak region)");
    
    DrawCentrality("sigma_phi", nHists, graphs[8], 0.1, 0.5, "#sigma_{#varphi} (rad.)", graphs[8+16], (graphsWingRemoved) ? graphsWingRemoved[8] : 0, 0, 0, 0, 1);
    DrawCentrality("sigma_eta", nHists, graphs[9], 0.1, 0.5, "#sigma_{#eta}", graphs[9+16], (graphsWingRemoved) ? graphsWingRemoved[9] : 0, 0, 0, 0, 1);
  }
  
  DrawCentrality("kurtosisphi_centrality", 12, graphs[14], -2, 2, "Kurtosis #varphi", graphs[14+16], (graphsWingRemoved) ? graphsWingRemoved[14] : 0, 0, 0, 0, 2);
  DrawCentrality("kurtosiseta_centrality", 12, graphs[15], -2, 2, "Kurtosis #eta", graphs[15+16], (graphsWingRemoved) ? graphsWingRemoved[15] : 0, 0, 0, 0, 2);
}
*/

void MCComparison(const char* fileNameData, const char* fileNameWingRemoved, const char* fileNameHijing, const char* fileNameAMPT, Int_t offset = 0)
{
  Int_t nHists = 12; //NHists;

  ReadGraphs(fileNameWingRemoved);
  CalculateRMSSigma();
  TGraphErrors*** graphsWingRemoved = graphs;

  ReadGraphs(fileNameHijing);
  CalculateRMSSigma();
  TGraphErrors*** graphs1 = graphs;

  CreateGraphStructure();
  TGraphErrors*** graphs2 = graphs;
  if (fileNameAMPT)
  {
    ReadGraphs(fileNameAMPT);
    CalculateRMSSigma();
    graphs2 = graphs;
  }

  ReadGraphs(fileNameData);
  CalculateRMSSigma();

  if (0)
  {
    Int_t graphList[] = { 1, 2, 8, 9, 14, 15 };
    for (Int_t i=0; i<6; i++)
    {
      graphs[graphList[i]+offset][5] = new TGraphErrors;
      graphs[graphList[i]+offset][10] = new TGraphErrors;
      graphs1[graphList[i]+offset][5] = new TGraphErrors;
      graphs1[graphList[i]+offset][10] = new TGraphErrors;
    }
  }
  
  DrawCentrality("phi_rms_centrality_mc", nHists, graphs[1+offset], 0.2, 0.8, Form("#sigma_{#Delta#varphi} (%s) (rad.)", fitLabel), graphsWingRemoved[1+offset], 0, graphs1[1+offset], 0, graphs2[1+offset], 1);
  
//   return;
  
  DrawCentrality("eta_rms_centrality_mc", nHists, graphs[2+offset], 0.2, 0.8 + 0.4 / 16 * offset, Form("#sigma_{#Delta#eta} (%s)", fitLabel), graphsWingRemoved[2+offset], 0, graphs1[2+offset], 0, graphs2[2+offset], 1);
  
  DrawCentrality("sigma_phi_mc", nHists, graphs[8+offset], 0.2, 0.6, "#sigma_{#Delta#varphi} (rad.)", graphsWingRemoved[8+offset], 0, graphs1[8+offset], 0, graphs2[8+offset], 1);
  DrawCentrality("sigma_eta_mc", nHists, graphs[9+offset], 0.2, 0.6, "#sigma_{#Delta#eta}", graphsWingRemoved[9+offset], 0, graphs1[9+offset], 0, graphs2[9+offset], 1);

  DrawCentrality("kurtosisphi_centrality_mc", nHists, graphs[14+offset], -1.5, 2.5, "Kurtosis #Delta#varphi", graphsWingRemoved[14+offset], 0, graphs1[14+offset], 0, graphs2[14+offset], 2);
  DrawCentrality("kurtosiseta_centrality_mc", nHists, graphs[15+offset], -1.5, 2.5, "Kurtosis #Delta#eta", graphsWingRemoved[15+offset], 0, graphs1[15+offset], 0, graphs2[15+offset], 2);

/*  DrawCentrality("phi_rms_centrality_mc", nHists, graphs[1], 0, 0.8, "#sigma_{#varphi} (fit) (rad.)", graphs[1+16], graphsWingRemoved[1],  graphs1[1], 0, graphs2[1], 1);
  DrawCentrality("eta_rms_centrality_mc", nHists, graphs[2], 0, 0.8, "#sigma_{#eta} (fit)", graphs[2+16], graphsWingRemoved[2], graphs1[2], 0, graphs2[2], 1);

  DrawCentrality("sigma_phi", nHists, graphs[8], 0.1, 0.5, "#sigma_{#varphi} (rad.)", graphs[8+16], graphsWingRemoved[8], graphs1[8], 0, graphs2[8], 1);
  DrawCentrality("sigma_eta", nHists, graphs[9], 0.1, 0.5, "#sigma_{#eta}", graphs[9+16], graphsWingRemoved[9], graphs1[9], 0, graphs2[9], 1);

  DrawCentrality("kurtosisphi_centrality_mc", nHists, graphs[14], -2, 2, "Kurtosis #varphi", graphs[14+16], graphsWingRemoved[14], graphs1[14], 0, graphs2[14], 2);
  DrawCentrality("kurtosiseta_centrality_mc", nHists, graphs[15], -2, 2, "Kurtosis #eta", graphs[15+16], graphsWingRemoved[15], graphs1[15], 0, graphs2[15], 2);*/
}

void ShowWingEffect(const char* fileNameData, const char* fileNameWingRemoved, Int_t offset = 0)
{
  Int_t nHists = 24; //NHists;

  ReadGraphs(fileNameWingRemoved);
  CalculateRMSSigma();
  TGraphErrors*** graphs1 = graphs;

  ReadGraphs(fileNameData);
  CalculateRMSSigma();
  
  DrawCentrality("phi_rms_centrality_mc", nHists, graphs[1+offset], 0, 0.8, "#sigma_{#varphi} (fit) (rad.)", graphs1[1+offset]);
  DrawCentrality("eta_rms_centrality_mc", nHists, graphs[2+offset], 0, 0.8, "#sigma_{#eta} (fit)", graphs1[2+offset]);

  DrawCentrality("sigma_phi", nHists, graphs[8+offset], 0, 0.8, "#sigma_{#varphi} (rad.)", graphs1[8+offset]);
  DrawCentrality("sigma_eta", nHists, graphs[9+offset], 0, 0.8, "#sigma_{#eta}", graphs1[9+offset]);

  DrawCentrality("kurtosisphi_centrality_mc", nHists, graphs[14+offset], -2, 4, "Kurtosis #varphi", graphs1[14+offset]);
  DrawCentrality("kurtosiseta_centrality_mc", nHists, graphs[15+offset], -2, 4, "Kurtosis #eta", graphs1[15+offset]);
}

void ShowFitEffect(const char* fileNameData)
{
  Int_t nHists = 12; //NHists;

  ReadGraphs(fileNameData);
  CalculateRMSSigma();
  
  DrawCentrality("phi_rms_centrality_mc", nHists, graphs[1], 0, 0.8, "#sigma_{#varphi} (fit) (rad.)", graphs[1+16]);
  DrawCentrality("eta_rms_centrality_mc", nHists, graphs[2], 0, 0.8, "#sigma_{#eta} (fit)", graphs[2+16]);

  DrawCentrality("sigma_phi", nHists, graphs[8], 0, 0.8, "#sigma_{#varphi} (rad.)", graphs[8+16]);
  DrawCentrality("sigma_eta", nHists, graphs[9], 0, 0.8, "#sigma_{#eta}", graphs[9+16]);

  DrawCentrality("kurtosisphi_centrality_mc", nHists, graphs[14], -2, 4, "Kurtosis #varphi", graphs[14+16]);
  DrawCentrality("kurtosiseta_centrality_mc", nHists, graphs[15], -2, 4, "Kurtosis #eta", graphs[15+16]);
}

void CompareSigmas(const char* fileNameData)
{
  // compare sigma from fit with sigma from direct calculation but taking for the first the limited range of 0.8 into account
  
  Int_t nHists = 12; //NHists;

  ReadGraphs(fileNameData);
  CalculateRMSSigma();
  TGraphErrors*** graphsOriginal = graphs;
  
  ReadGraphs(fileNameData);
  
  // rms
  for (Int_t histId = 0; histId < nHists; histId++)
  {
    if (!graphs[0][histId])
      continue;
    
    for (Int_t i=0; i<graphs[0][histId]->GetN(); i++)
    {
      // phi
      TF1* func = new TF1("func", "[0]+[1]*([4]/[2]/TMath::Sqrt(TMath::TwoPi())*exp(-0.5*(x/[2])**2)*TMath::Erf([7]/TMath::Sqrt(2)/[3])+(1-[4])/[5]/TMath::Sqrt(TMath::TwoPi())*exp(-0.5*(x/[5])**2)*TMath::Erf([7]/TMath::Sqrt(2)/[6]))", -0.5 * TMath::Pi(), 0.5 * TMath::Pi());
      func->SetParameter(0, 0);
      for (Int_t k=0; k<6; k++)
	func->SetParameter(k+1, graphs[k][histId]->GetY()[i]);
      func->SetParameter(7, 0.8); // project Limit
      
      graphs[8][histId]->GetY()[i] = TMath::Sqrt(func->CentralMoment(2, -0.87, 0.87));
      
      // eta
      func = new TF1("func", "[0]+[1]*([4]/[3]/TMath::Sqrt(TMath::TwoPi())*exp(-0.5*(x/[3])**2)*TMath::Erf([7]/TMath::Sqrt(2)/[2])+(1-[4])/[6]/TMath::Sqrt(TMath::TwoPi())*exp(-0.5*(x/[6])**2)*TMath::Erf([7]/TMath::Sqrt(2)/[5]))", -kOuterLimit, kOuterLimit);
      func->SetParameter(0, 0);
      for (Int_t k=0; k<6; k++)
	func->SetParameter(k+1, graphs[k][histId]->GetY()[i]);
      func->SetParameter(7, 0.87); // project Limit
      
      graphs[9][histId]->GetY()[i] = TMath::Sqrt(func->CentralMoment(2, -0.8, 0.8));
    }
  }
  
  DrawCentrality("phi_comparison1", nHists, graphsOriginal[1], 0.2, 0.8, "#sigma_{#Delta#varphi} (fit)", 0, 0, graphs[8]);
  DrawCentrality("phi_comparison2", nHists, graphsOriginal[8], 0.2, 0.5, "#sigma_{#Delta#varphi} (fit)", 0, 0, graphs[8]);
  DrawCentrality("eta_comparison1", nHists, graphsOriginal[2], 0.2, 0.8, "#sigma_{#Delta#eta} (fit)", 0, 0, graphs[9]);
  DrawCentrality("eta_comparison2", nHists, graphsOriginal[9], 0.2, 0.5, "#sigma_{#Delta#eta} (fit)", 0, 0, graphs[9]);
}  

Float_t** ExtractSystematics(const char* baseFile, const char* systFile, Int_t offset)
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
  
    TCanvas* c = new TCanvas(Form("%s_%d_%d", systFile, i, 0), Form("%s_%d_%d", systFile, i, 0), 1000, 1000);
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
      
      TGraphErrors* graph1 = graphsBase[graphList[i]+offset][j];
      TGraphErrors* graph2 = (systFile) ? graphs[graphList[i]+offset][j] : graphsBase[graphList[i]+16][j];
      
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
      {
// 	hist->Fill(graph1->GetY()[k], 1.0 / (graph1->GetEY()[k] / graph1->GetY()[k]));
	if (kurtosis)
	  hist->Fill(TMath::Abs(graph1->GetY()[k]));
	else
	  hist->Fill(graph1->GetY()[k]);
      }
      
      if (count == 37)
	break;
    }
    
    c->SaveAs(Form("syst/%s_%d.eps", systFile, i));
    
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
    
    if (!kurtosis)
    {
      mean -= 1;
      mean *= 100;
      sigma *= 100;
    }
    results[i][0] = mean;
    results[i][1] = sigma / TMath::Sqrt(hist->GetEntries());
  }
  
  return results;
}

void GetProjections(TH2* hist, TH1** projPhi, TH1** projEta, Int_t k)
{
  Float_t projectLimit = 0.8;

  *projPhi = hist->ProjectionX(Form("%s_proj1%d", hist->GetName(), k), hist->GetYaxis()->FindBin(-projectLimit+0.01), hist->GetYaxis()->FindBin(projectLimit-0.01));
    
  *projEta = hist->ProjectionY(Form("%s_proj2_%d", hist->GetName(), k), hist->GetXaxis()->FindBin(-projectLimit+0.01), hist->GetXaxis()->FindBin(projectLimit-0.01));
}

// uncorrected ones
const char* systBaseFile = "dphi_corr_2d_120429.root";
const char* systTrackCuts1 = "dphi_corr_120425_hybrid.root";
const char* systTrackCuts2 = "dphi_corr_120430_raa.root";
const char* systVertex = "dphi_corr_2d_120425_vertex.root";
const char* systResonances = "dphi_corr_120502_resonances.root";
const char* systTTR = "dphi_corr_2d_120430_widettr.root";

// corrected ones
const char* systBaseFileCorrected = "dphi_corr_2d_120508.root";
const char* systTrackCuts1Corrected = "dphi_corr_120507_hybrid.root";
const char* systTrackCuts2Corrected = "dphi_corr_120507_raa.root";
const char* systWingRemovedCorrected = "dphi_corr_2d_120508_wingremoved.root";

void ExtractSystematicsProjections(Int_t mode, Float_t rangeBegin = -0.4, Float_t rangeEnd = 0.4, Bool_t draw = kFALSE)
{
  if (!draw)
    gROOT->SetBatch(kTRUE);
  
  Int_t maxLeadingPt = 4;
  Int_t maxAssocPt = 6;

  Int_t NEffects = 3;
  
  TFile* file1 = 0;
  const char** systFiles = 0;
  if (mode == 2)
  {
    NEffects = 5;
    file1 = TFile::Open(systBaseFile);
    static const char* systFilesTmp[5] = { systVertex, systResonances, systTTR, systTrackCuts1Corrected, systTrackCuts2Corrected };
    systFiles = systFilesTmp;
  }
  else if (mode == 0)
  {
    file1 = TFile::Open(systBaseFile);
    static const char* systFilesTmp[3] = { systVertex, systResonances, systTTR };
    systFiles = (const char**) systFilesTmp;
  }
  else if (mode == 1)
  {
    NEffects = 3;
    file1 = TFile::Open(systBaseFileCorrected);
    static const char* systFilesTmp[3] = { systTrackCuts1Corrected, systTrackCuts2Corrected, systWingRemovedCorrected };
    systFiles = systFilesTmp;
  }
  else if (mode == 3)
  {
    NEffects = 3;
    file1 = TFile::Open(systBaseFileCorrected);
    static const char* systFilesTmp[3] = { systBaseFileCorrected, systBaseFileCorrected, systBaseFileCorrected };
    systFiles = systFilesTmp;
  }
  
  TH2* uncertainty[4];
  for (Int_t i=0; i<4; i++)
    uncertainty[i] = new TH2F(Form("uncertainty_%d", i), Form("%s %s;ptbin;centrality", (i % 2 == 0) ? "#varphi" : "#eta", (i < 2) ?"average" : "deviation"), 30, 0, 30, 6, 0, 6);
  
  Int_t count = 0;
  Int_t ptBin = 0;
  for (Int_t i=0; i<maxLeadingPt-1; i++)
  {
    for (Int_t j=1; j<maxAssocPt-1; j++)
    {
      if (i == 0 && j == 2)
	continue;
      if (i == 1 && j == 3)
	continue;
      if (i == 2 && j == 4)
	continue;
      
      ptBin++;
      for (Int_t histId = 0; histId < NHists; histId++)
      {
	TH2* hist1 = (TH2*) file1->Get(Form("dphi_%d_%d_%d", i, j+1, histId));
	if (!hist1)
	  continue;
	
	if (hist1->GetEntries() < 1e4)
	{
	  Printf("Only %f entries. Skipping...", hist1->GetEntries());
	  continue;
	}
	
	Printf("%d %d %d %s", i, j, histId, hist1->GetTitle());
	
	TH1* proj1[2];
	TH1* proj2[10][2];
	
	SubtractEtaGap(hist1, kEtaLimit, kOuterLimit, kFALSE);
	GetProjections(hist1, &proj1[0], &proj1[1], count++);

	TCanvas* c = new TCanvas(Form("c_%d_%d_%d", mode, ptBin, histId), Form("c_%d_%d_%d", mode, ptBin, histId), 1200, 1000);
	c->Divide(3, 2);
	  
	for (Int_t k=0; k<2; k++)
	{
	  c->cd(1+k*3);
	  proj1[k]->SetStats(0);
	  proj1[k]->GetXaxis()->SetRangeUser(-1.5, 1.5);
	  proj1[k]->GetYaxis()->SetRangeUser(proj1[k]->GetMinimum() * 1.1, proj1[k]->GetMaximum() * 1.4);
	  proj1[k]->DrawCopy();
	}
	
	Double_t maxAverage[2] = { 0, 0 };
	Double_t maxDev[2] = { 0, 0 };
  
	for (Int_t n=0; n<NEffects; n++)
	{
	  if (histId == 2 && mode == 0 && n > 0)
	    continue;
	  if (histId == 2 && mode == 1 && n > 0)
	    continue;
// 	  Printf("%d %s", n, systFiles[n]);
	  TFile* file2 = TFile::Open(systFiles[n]);
	  hist1 = (TH2*) file2->Get(Form("dphi_%d_%d_%d", i, j+1, histId));
	  if (!hist1)
	    continue;
	  
	  if (mode == 3 && n == 0) // change eta limits
	    SubtractEtaGap(hist1, kEtaLimit-0.2, kOuterLimit, kFALSE);
	  if (mode == 3 && n == 1) // change eta limits
	    SubtractEtaGap(hist1, kEtaLimit+0.2, kOuterLimit, kFALSE);
	  if (mode == 3 && n == 2) // change eta limits
	    SubtractEtaGap(hist1, kEtaLimit, kOuterLimit-0.2, kFALSE);
	  else
	    SubtractEtaGap(hist1, kEtaLimit, kOuterLimit, kFALSE);
	  GetProjections(hist1, &proj2[n][0], &proj2[n][1], count++);
	  
	  for (Int_t k=0; k<2; k++)
	  {
	    c->cd(1+k*3);
	    proj2[n][k]->SetLineColor(n+2);
	    /* TH1* copy = */ proj2[n][k]->DrawCopy("SAME");
	    
	    c->cd(2+k*3);
	    gPad->SetGridx();
	    gPad->SetGridy();
	    TH1* ratio = (TH1*) proj2[n][k]->Clone(Form("%s_ratio", proj2[n][k]->GetName()));
	    ratio->SetStats(0);
	    ratio->Divide(proj1[k]);
	    ratio->GetXaxis()->SetRangeUser(-1.5, 1.5);
	    ratio->Fit("pol0", "0Q", "", -1, 1);
	    Double_t average = ratio->GetFunction("pol0")->GetParameter(0);
// 	    Printf("Average is %f", average);
	    maxAverage[k] = TMath::Max(maxAverage[k], TMath::Abs(average - 1));
// 	    ratio->Scale(1.0 / average);
// 	    copy->Scale(1.0 / average);
	    ratio->GetYaxis()->SetRangeUser(0.5, 1.5);
	    ratio->DrawCopy((n == 0) ? "" : "SAME");
	    
	    if (0)
	    {
	      Double_t sum = 0;
	      Double_t sumCount = 0;
	      for (Int_t bin = ratio->FindBin(rangeBegin); bin <= ratio->FindBin(rangeEnd); bin++)
	      {
		sum += TMath::Abs(ratio->GetBinContent(bin) - 1);
		sumCount++;
	      }
	      maxDev[k] = TMath::Max(maxDev[k], sum / sumCount);
	    }

	    c->cd(3+k*3);
	    gPad->SetGridx();
	    gPad->SetGridy();
	    TH1* diff = (TH1*) proj2[n][k]->Clone(Form("%s_diff", proj2[n][k]->GetName()));
	    diff->Scale(1.0 / average);
	    diff->Add(proj1[k], -1);
	    diff->SetStats(0);
	    diff->GetXaxis()->SetRangeUser(-1.5, 1.5);
	    diff->GetYaxis()->SetRangeUser(-0.05, 0.05);
	    diff->DrawCopy((n == 0) ? "" : "SAME");
	    
	    if (0)
	    {
	      for (Int_t bin = diff->FindBin(rangeBegin); bin <= diff->FindBin(rangeEnd); bin++)
		maxDev[k] = TMath::Max(maxDev[k], TMath::Abs(diff->GetBinContent(bin)));
	    }
	    else
	    {
	      Double_t sum = 0;
	      Double_t sumCount = 0;
	      for (Int_t bin = diff->FindBin(rangeBegin); bin <= diff->FindBin(rangeEnd); bin++)
	      {
		sum += TMath::Abs(diff->GetBinContent(bin));
		sumCount++;
	      }
	      maxDev[k] = TMath::Max(maxDev[k], sum / sumCount);
	    }
	  }
	  
	  delete file2;
	}
	
	// normalize to peak strength
/*	for (Int_t k=0; k<2; k++)
	  maxDev[k] /= proj1[k]->Integral(proj1[k]->GetXaxis()->FindBin(-0.2), proj1[k]->GetXaxis()->FindBin(0.2)) / (proj1[k]->GetXaxis()->FindBin(0.2) - proj1[k]->GetXaxis()->FindBin(-0.2) + 1);*/
	
	Int_t centralityAxisMapping[] = { 0, 4, 3, 5, 1, 2 };
	uncertainty[0]->SetBinContent(ptBin, centralityAxisMapping[histId]+1, maxAverage[0]);
	uncertainty[1]->SetBinContent(ptBin, centralityAxisMapping[histId]+1, maxAverage[1]);
	uncertainty[2]->SetBinContent(ptBin, centralityAxisMapping[histId]+1, maxDev[0]);
	uncertainty[3]->SetBinContent(ptBin, centralityAxisMapping[histId]+1, maxDev[1]);
	
	c->SaveAs(Form("syst_corr/%s.eps", c->GetName()));
      }
      
      if (draw)
	break;
    }
    if (draw)
      break;
  }
  
  TFile::Open("uncertainty.root", "RECREATE");
  for (Int_t i=0; i<4; i++)
    uncertainty[i]->Write();
  gFile->Close();
}

void DrawSystematicsProjections()
{
  TH2* uncertainty[4];
  
  TFile::Open("uncertainty.root");
  for (Int_t i=0; i<4; i++)
    uncertainty[i] = (TH2*) gFile->Get(Form("uncertainty_%d", i));
  
  TCanvas* cResult = new TCanvas("cResult", "cResult", 1000, 1000);
  cResult->Divide(2, 2);
  for (Int_t i=0; i<4; i++)
  {
    cResult->cd(i+1);
    gPad->SetRightMargin(0.15);
    uncertainty[i]->SetStats(0);
    uncertainty[i]->Draw("colz");
  }

  TCanvas* cResult2 = new TCanvas("cResult2", "cResult2", 1000, 1000);
  cResult2->Divide(2, 2);
  for (Int_t i=0; i<4; i++)
  {
    cResult2->cd(i+1);
    gPad->SetRightMargin(0.15);
    Int_t color = 1;
    for (Int_t j=1; j<=uncertainty[i]->GetNbinsX(); j++)
    {
      TH1* proj = uncertainty[i]->ProjectionY(Form("p_%d_%d", i, j), j, j);
      if (proj->Integral() > 0)
      {
	proj->SetLineColor(color++);
	proj->SetStats(0);
	proj->GetYaxis()->SetRangeUser(0, 0.5);
	proj->Draw((color == 2) ? "" : "SAME");
	
	Printf("%d %d: %.3f %.3f", i, color, proj->GetBinContent(1), proj->Integral(2, 6) / 5);
      }
    }
  }
}

void BuildSystematicFiles()
{
  Printf("Hope you know what you are doing... 5 seconds to abort...");
  gSystem->Sleep(5000);

  AnalyzeDeltaPhiEtaGap2D(systBaseFile, "syst_base.root");

  kEtaLimit = 0.8;
  AnalyzeDeltaPhiEtaGap2D(systBaseFile, "syst_eta08.root");
  
  kEtaLimit = 1.2;
  AnalyzeDeltaPhiEtaGap2D(systBaseFile, "syst_eta12.root");
  kEtaLimit = 1.0;
  
  kOuterLimit = 1.39;
  AnalyzeDeltaPhiEtaGap2D(systBaseFile, "syst_outer14.root");
  kOuterLimit = 1.59;
  
  AnalyzeDeltaPhiEtaGap2D(systTrackCuts1, "syst_hybrid.root");
  AnalyzeDeltaPhiEtaGap2D(systTrackCuts2, "syst_raa.root");
  AnalyzeDeltaPhiEtaGap2D(systVertex, "syst_vertex.root");
  AnalyzeDeltaPhiEtaGap2D(systResonances, "syst_resonances.root");
  AnalyzeDeltaPhiEtaGap2D(systTTR, "syst_widettr.root");
}

void ExtractSystematicsAll(Int_t offset = 0)
{
  gROOT->SetBatch(kTRUE);
  
  const Int_t NEffects = 9;
  
//   const char* defaultFile = "graphs_120425.root";
//   const char* systFiles[] = { "graphs_120425_eta08.root", "graphs_120425_eta12.root", "graphs_120425_outer14.root", "graphs_hybrid_120427.root", "graphs_120430_raa.root", "graphs_120425_vertex.root", "graphs_120502_resonances.root", "graphs_120502_widettr.root", "graphs_120425_wingremoved.root" };
  const char* defaultFile = "syst_base.root";
  const char* systFiles[] = { "syst_eta08.root", "syst_eta12.root", "syst_outer14.root", "syst_hybrid.root", "syst_raa.root", "syst_vertex.root", "syst_resonances.root", "syst_widettr.root", "graphs_120429_wingremoved.root" };
  const char* effectIdString[] = { "1a", "1b", "1c", "2a", "2b", "3", "4", "5", "6" };

  Float_t** results[NEffects+3];
  
  for (Int_t i=0; i<NEffects; i++)
  {
    results[i] = ExtractSystematics(defaultFile, systFiles[i], offset);
  }
  
  const Int_t NParameters = 6;
  const char* names[] = { "$\\sigma_{\\Dphi}$ (fit)", "$\\sigma_{\\Deta}$ (fit)", "$\\sigma_{\\Dphi}$", "$\\sigma_{\\Deta}$", "Kurtosis $\\Dphi$", "Kurtosis $\\Deta$" };
  
  for (Int_t j=0; j<NParameters; j++)
  {
    printf("%s \t & ", names[j]);
    for (Int_t i=0; i<NEffects; i++)
    {
      if (j < 4)
	printf("$%.1f\\%% \\pm %.1f\\%%$ \t & ", results[i][j][0], results[i][j][1]);
      else
	printf("$%.2f \\pm %.2f$ \t & ", results[i][j][0], results[i][j][1]);
    }
    Printf("\\\\");
  }

  // put together 0-2 (into NEffects)
  results[NEffects] = new Float_t*[NParameters];
  for (Int_t j=0; j<NParameters; j++)
  {
    results[NEffects][j] = new Float_t[2];
    Float_t mean = 0;
    Float_t sigma = 0;
    
    printf("%s \t:", names[j]);
    for (Int_t i=0; i<=2; i++)
    {
      mean = TMath::Max(mean, TMath::Abs(results[i][j][0]));
      printf("%.1f%% ", results[i][j][0]);
    }
    printf("--> %.1f%% \t", mean);
    results[NEffects][j][0] = mean;
      
    for (Int_t i=0; i<=2; i++)
    {
      sigma = TMath::Max(sigma, TMath::Abs(results[i][j][1]));
      printf("%.1f%% ", results[i][j][1]);
    }
    Printf("--> %.1f%% \t", sigma);
    results[NEffects][j][1] = sigma;
  }

  // put together 3-4 (into NEffects+1)
  results[NEffects+1] = new Float_t*[NParameters];
  for (Int_t j=0; j<NParameters; j++)
  {
    results[NEffects+1][j] = new Float_t[2];
    
    Float_t mean = 0;
    Float_t sigma = 0;
    
    printf("%s \t:", names[j]);
    for (Int_t i=3; i<=4; i++)
    {
      mean = TMath::Max(mean, TMath::Abs(results[i][j][0]));
      printf("%.1f%% ", results[i][j][0]);
    }
    printf("--> %.1f%% \t", mean);
    results[NEffects+1][j][0] = mean;
      
    for (Int_t i=3; i<=4; i++)
    {
      sigma = TMath::Max(sigma, TMath::Abs(results[i][j][1]));
      printf("%.1f%% ", results[i][j][1]);
    }
    Printf("--> %.1f%% \t", sigma);
    results[NEffects+1][j][1] = sigma;
  }
  
  // combine in quadrature (only mean)
  results[NEffects+2] = new Float_t*[NParameters];
  for (Int_t j=0; j<NParameters; j++)
  {
    results[NEffects+2][j] = new Float_t[2];
    Float_t mean = 0;
  
    printf("%s \t:", names[j]);
    for (Int_t i=5; i<NEffects+2; i++)
    {
      mean += results[i][j][0] * results[i][j][0];
      printf("%.1f%% ", results[i][j][0]);
    }
    mean = TMath::Sqrt(mean);
    Printf("--> %.1f%% \t", mean);
    results[NEffects+2][j][0] = mean;
  }
    
  for (Int_t j=0; j<NParameters; j++)
  {
    printf("%s \t & ", names[j]);
    for (Int_t i=0; i<NEffects+3; i++)
    {
      if (i == NEffects || i == NEffects+1)
	continue;
      if (j < 4)
	printf("$%.1f\\%%$ \t ", results[i][j][0]);
      else
	printf("$%.2f$ \t ", results[i][j][0]);
      if (i < NEffects+2)
	printf("& ");
    }
    Printf("\\\\");
  }
  
  // generate LaTeX code
  for (Int_t i=0; i<NEffects; i++)
  {
    for (Int_t j=0; j<NParameters; j++)
    {
      Printf("\\bfig[ht]");
      Printf("  \\includegraphics[width=\\linewidth]{syst/%s_%d.eps}", systFiles[i], j);
      Printf("  \\caption{Systematic effect (%s) - Parameter %s}", effectIdString[i], names[j]);
      Printf("\\efig");
    }
    Printf("\\clearpage");
  }
}

void CompareEtaPhi(const char* fileName, const char* fileNameWingRemoved)
{
  Int_t offset = 0;
  
  TGraphErrors*** graphsWingRemoved = 0;
  if (fileNameWingRemoved)
  {
    ReadGraphs(fileNameWingRemoved);
    CalculateRMSSigma();
    graphsWingRemoved = graphs;
  }
  
  ReadGraphs(fileName);
  
  Int_t nHists = 12; //NHists;
  skipGraphList[0] = 5;
  skipGraphList[1] = 6;
  skipGraphList[2] = 10;
  
  CalculateRMSSigma();
  
  drawLogo = 1;
  
  MCLabel = "#sigma_{#Delta#varphi}     #sigma_{#Delta#eta}";

  DrawCentrality("rms_centrality", nHists, graphs[1+offset], 0.2, 0.8, "#sigma_{#Delta#varphi} (fit) (rad.) , #sigma_{#Delta#eta} (fit)", graphsWingRemoved[1+offset], 0, graphs[2+offset], graphsWingRemoved[2+offset], 0, 1);
//   DrawCentrality("rms_centrality_nosyst", nHists, graphs[1], 0, 0.8, "#sigma_{#varphi} (fit) (rad.) / #sigma_{#eta}", 0, 0,  graphs[2], 0);
}

void DrawExamples(const char* histFileName, const char* graphFileName, Int_t i = 0, Int_t j = 1, Bool_t drawFunc = kTRUE)
{
  Float_t etaLimit = kEtaLimit;
  Float_t outerLimit = kOuterLimit;
  Float_t projectLimit = 0.8;

  ReadGraphs(graphFileName);

  TFile::Open(histFileName);
  
  Int_t exColors[] = { 1, 2, 4, 3, 5, 6 };
  
  Int_t graphID = i * (6 - 1) + j - 1;

  TCanvas* c = new TCanvas(Form("c_%d", i), Form("c_%d", i), 1200, 600);
  c->Divide(2, 1);
  Int_t nHists = 3;
  for (Int_t histId = 0; histId < nHists; histId++)
  {
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
    proj->SetStats(0);
    proj->GetXaxis()->SetRangeUser(-1.49, 1.49);
    proj->GetYaxis()->SetRangeUser(proj->GetMinimum() * 1.2, proj->GetMaximum() * 1.5);
    proj->GetYaxis()->SetTitle(kProjYieldTitlePhi);
    TString label(proj->GetTitle());
    TObjArray* objArray = label.Tokenize("-");
    proj->SetTitle(Form("%s-%s", objArray->At(0)->GetName(), objArray->At(1)->GetName()));
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
    proj->GetYaxis()->SetTitle(kProjYieldTitleEta);
    proj->SetStats(0);
    proj->SetTitle(Form("%s-%s", objArray->At(0)->GetName(), objArray->At(1)->GetName()));
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

TH1* GetSystUnc(TH1* hist, Int_t i, Int_t j, Int_t histId, Int_t etaPhi)
{
//   Float_t uncertainty = (histId == 0) ? 0.1 : 0.08;

  // track cuts
  //   2 2: 0.048 0.020
  //   2 3: 0.078 0.031
  //   2 4: 0.023 0.008
  //   2 5: 0.159 0.056
  //   2 6: 0.043 0.015
  //   2 7: 0.011 0.006
  //   3 2: 0.050 0.022
  //   3 3: 0.089 0.033
  //   3 4: 0.020 0.008
  //   3 5: 0.183 0.056
  //   3 6: 0.042 0.013
  //   3 7: 0.010 0.005
  
  // others
  //   2 2: 0.023 0.007
  //   2 3: 0.024 0.015
  //   2 4: 0.008 0.005
  //   2 5: 0.075 0.024
  //   2 6: 0.019 0.010
  //   2 7: 0.007 0.004
  //   3 2: 0.029 0.009
  //   3 3: 0.032 0.012
  //   3 4: 0.005 0.003
  //   3 5: 0.086 0.013
  //   3 6: 0.015 0.006
  //   3 7: 0.008 0.003

  Float_t uncertainty = 1;
  if (i == 0 && j == 1)
    uncertainty = 0.021;
  if (i == 1 && j == 1)
    uncertainty = 0.032;
  if (i == 1 && j == 2)
    uncertainty = 0.008;
  if (i == 2 && j == 1)
    uncertainty = 0.056;
  if (i == 2 && j == 2)
    uncertainty = 0.014;
  if (i == 2 && j == 3)
    uncertainty = 0.006;
  
  if (histId == 0)
    uncertainty *= 2;
  
  TH1* systUnc = (TH1*) hist->Clone(Form("%s_syst", hist->GetName()));
  for (Int_t n=1; n<=systUnc->GetNbinsX(); n++)
    systUnc->SetBinError(n, uncertainty);
  
  systUnc->SetFillColor(hist->GetLineColor());
  systUnc->SetFillStyle((etaPhi == 0) ? 3004 : 3005);
  systUnc->SetMarkerStyle(0);
  systUnc->SetLineColor(0);
  
  return systUnc;
}

Bool_t disableUncertainties = kTRUE;

void DrawExample(const char* histFileName, Int_t i, Int_t j, Int_t histId, TH1** projPhi, TH1** projEta, TH1** projPhiSyst, TH1** projEtaSyst, Bool_t Ratio = kFALSE)
{
  Float_t etaLimit = kEtaLimit;
  Float_t outerLimit = kOuterLimit;
  Float_t projectLimit = 0.8;
  Int_t maxAssocPt = 6;
  Int_t graphID = i * (maxAssocPt - 1) + j - 1;

  TFile::Open(histFileName);
  
  TCanvas* c = new TCanvas(Form("ex_%d_%d_%d", i, j, histId), Form("ex_%d_%d_%d", i, j, histId), 1800, 600);
  c->Divide(3, 1);

  TH2* hist = (TH2*) gFile->Get(Form("dphi_%d_%d_%d", i, j+1, histId));
  if (!hist)
    return;
  // NOTE fix normalization. these 2d correlations come out of AliUEHist normalized by dphi bin width, but not deta
  hist->Scale(1.0 / hist->GetYaxis()->GetBinWidth(1));
  hist->GetZaxis()->SetTitle(kCorrFuncTitle);
  hist->GetZaxis()->SetTitleOffset(1.9);
  
  if (graphID >= 15)
  {
    etaLimit = 0.5;
    outerLimit = 0.99;
  }
//   hist->Rebin2D(2, 2); hist->Scale(0.25);
  
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
  
  c->cd(1);
  gPad->SetLeftMargin(0.15);
  hist->GetYaxis()->SetRangeUser(-outerLimit, outerLimit);
  hist->GetXaxis()->SetTitleOffset(1.5);
  hist->GetYaxis()->SetTitleOffset(1.7);
  hist->SetStats(0);
  hist->SetTitle("a) Correlation");
  TH2* clone = (TH2*) hist->Clone(Form("%s_clone", hist->GetName()));
  clone->GetXaxis()->SetRangeUser(-TMath::Pi() / 2, TMath::Pi() / 2);
  clone->Draw("SURF1");
  paveText->Draw();
  
  SubtractEtaGap(hist, etaLimit, outerLimit, kFALSE, kFALSE);
  c->cd(2);
  gPad->SetLeftMargin(0.15);
  hist->SetTitle("b) #eta-gap subtracted");
  hist->GetXaxis()->SetRangeUser(-TMath::Pi() / 2, TMath::Pi() / 2);
  hist->DrawCopy("SURF1");
  if (!disableUncertainties)
    DrawALICELogo(kTRUE, 0.7, 0.7, 0.9, 0.9);
    
  c->cd(3);
  gPad->SetLeftMargin(0.13);
  TH1* proj = hist->ProjectionX(Form("%s_proj1", hist->GetName()), hist->GetYaxis()->FindBin(-projectLimit+0.01), hist->GetYaxis()->FindBin(projectLimit-0.01));
  TH1* proj2 = hist->ProjectionY(Form("%s_proj2b", hist->GetName()), hist->GetXaxis()->FindBin(-projectLimit+0.01), hist->GetXaxis()->FindBin(projectLimit-0.01));
  // normalization
  proj->Scale(hist->GetYaxis()->GetBinWidth(1));
  proj2->Scale(hist->GetXaxis()->GetBinWidth(1));

  proj->SetStats(0);
  proj->SetTitle("c) Projections");
  proj->GetXaxis()->SetRangeUser(-1.49, 1.49);
  proj->GetYaxis()->SetTitle(kProjYieldTitlePhi);
  proj->GetYaxis()->SetTitleOffset(1.3);
  proj->GetXaxis()->SetTitleOffset(1);
  proj->GetYaxis()->SetRangeUser(proj->GetMinimum(), proj->GetMaximum() * 1.6);
  TH1* systUncPhi = GetSystUnc(proj, i, j, histId, 0);
  TH1* copy = systUncPhi->DrawCopy("E2 ][");
  copy->GetXaxis()->SetTitle(Form("%s , %s", proj->GetXaxis()->GetTitle(), proj2->GetXaxis()->GetTitle()));
  copy->GetYaxis()->SetTitle(kProjYieldTitlePhiOrEta);
    
  proj2->SetLineColor(2);
  proj2->SetStats(0);
  proj2->GetYaxis()->SetTitle(kProjYieldTitleEta);
  proj2->GetYaxis()->SetTitleOffset(1.2);
  proj2->GetXaxis()->SetTitleOffset(1);
  proj2->GetYaxis()->SetRangeUser(proj2->GetMinimum(), proj2->GetMaximum() * 1.6);
  TH1* systUncEta = GetSystUnc(proj2, i, j, histId, 1);
  systUncEta->Draw("E2 ][ SAME");
  proj->DrawCopy((disableUncertainties) ? "" : "SAME");
  proj2->DrawCopy("SAME");
  if (!disableUncertainties)
  {
    if (histId == 0)
      DrawLatex(0.3, 0.85, 1, "Scale uncertainty: 20%", 0.04);
    else
      DrawLatex(0.3, 0.85, 1, "Scale uncertainty: 10%", 0.04);
  }
  DrawLatex(0.3, 0.80, 1, Form("#Delta#varphi projection in |#Delta#eta| < %.2f", hist->GetYaxis()->GetBinUpEdge(hist->GetYaxis()->FindBin(projectLimit-0.01))), 0.04);
  DrawLatex(0.3, 0.75, 2, Form("#Delta#eta projection in |#Delta#varphi| < %.2f", hist->GetXaxis()->GetBinUpEdge(hist->GetXaxis()->FindBin(projectLimit-0.01))), 0.04);
  
  proj->SetTitle(label);
  proj2->SetTitle(label);
  systUncPhi->SetTitle(label);
  systUncEta->SetTitle(label);
  
  TH1* etaRatio = (TH1*) proj2->Clone("etaRatio");
  etaRatio->GetXaxis()->SetRangeUser(0,2);
  etaRatio->GetYaxis()->SetRangeUser(etaRatio->GetMinimum(), etaRatio->GetMaximum() * 1.6);
  TH1* phiRatio = (TH1*) proj->Clone("phiRatio");
  phiRatio->GetXaxis()->SetRangeUser(0,2);
  phiRatio->GetYaxis()->SetRangeUser(phiRatio->GetMinimum(), phiRatio->GetMaximum() * 1.6);

  for (Int_t x=1; x<=proj->GetNbinsX(); x++)
    if (proj->GetBinCenter(x)>=0)
    {
      phiRatio->SetBinContent(x,proj->GetBinContent(proj->FindBin(proj->GetBinCenter(x)))+proj->GetBinContent(proj->FindBin(-1*proj->GetBinCenter(x))));
      phiRatio->SetBinError(x,sqrt((proj->GetBinError(proj->FindBin(proj->GetBinCenter(x)))*proj->GetBinError(proj->FindBin(proj->GetBinCenter(x)))+proj->GetBinError(proj->FindBin(-1*proj->GetBinCenter(x)))*proj->GetBinError(proj->FindBin(-1*proj->GetBinCenter(x))))));
    }
  for (Int_t x=1; x<=proj2->GetNbinsX(); x++)
    if (proj2->GetBinCenter(x)>=0)
    {
      etaRatio->SetBinContent(x,proj2->GetBinContent(proj2->FindBin(proj2->GetBinCenter(x)))+proj2->GetBinContent(proj2->FindBin(-1*proj2->GetBinCenter(x))));
      etaRatio->SetBinError(x,sqrt((proj2->GetBinError(proj2->FindBin(proj2->GetBinCenter(x)))*proj2->GetBinError(proj2->FindBin(proj2->GetBinCenter(x)))+proj2->GetBinError(proj2->FindBin(-1*proj2->GetBinCenter(x)))*proj2->GetBinError(proj2->FindBin(-1*proj2->GetBinCenter(x))))));
    }

  if (Ratio)
  {
    *projEta = etaRatio;
    *projPhi = phiRatio;
  }
  else
  {
    *projPhi = proj;
    *projEta = proj2;
  }
  *projPhiSyst = systUncPhi;
  *projEtaSyst = systUncEta;

  c->SaveAs(Form("ex/%s.png", c->GetName()));
  c->SaveAs(Form("ex/%s.eps", c->GetName()));

  if (1)
  {
    c = new TCanvas(Form("ex_%d_%d_%d_a", i, j, histId), Form("ex_%d_%d_%d_a", i, j, histId), 800, 800);
    gPad->SetLeftMargin(0.15);
    hist->SetTitle("");
    hist->DrawCopy("SURF1");
    paveText->Draw();
    DrawALICELogo(kTRUE, 0.2, 0.75, 0.4, 0.95);
    c->SaveAs(Form("ex/%s.png", c->GetName()));
    c->SaveAs(Form("ex/%s.eps", c->GetName()));

    c = new TCanvas(Form("ex_%d_%d_%d_b", i, j, histId), Form("ex_%d_%d_%d_b", i, j, histId), 800, 800);
    gPad->SetLeftMargin(0.15);
    clone->SetTitle("");
    clone->DrawCopy("SURF1");
    paveText->Draw();
    DrawALICELogo(kTRUE, 0.2, 0.75, 0.4, 0.95);
    c->SaveAs(Form("ex/%s.png", c->GetName()));
    c->SaveAs(Form("ex/%s.eps", c->GetName()));
  }
}

void DrawExampleAll(const char* histFileName, const char* graphFileName = 0)
{
  if (graphFileName)
    ReadGraphs(graphFileName);
  
  Int_t colors[] = { 1, 2, 4 };
  
  //Float_t exampleI[] = { 0, 2, 2 };
  //Float_t exampleJ[] = { 1, 1, 2 };
  //Int_t exampleI[] = { 0, 1, 2 };
  //Int_t exampleJ[] = { 1, 1, 1 };
  Float_t exampleI[] = { 3, 3, 3 };
  Float_t exampleJ[] = { 3, 4, 5 };
  
  TH1* projectionsPhi[9];
  TH1* projectionsEta[9];
  TH1* projectionsPhiSyst[9];
  TH1* projectionsEtaSyst[9];

  Bool_t Ratio = kFALSE;
  //if Ratio is true it divides by the 60-70% centrality bin the other two 
  //the label on the Y axis is wrong

  disableUncertainties = 1;
  
  Int_t count = 0;
  for (Int_t i=0; i<3; i++)
  {
    for (Int_t histId = 0; histId<3; histId++)
    {
      //it is using the 0-10%, the 20-40% and the 60-70% centrality bins
      DrawExample(histFileName, exampleI[i], exampleJ[i], (histId==2) ? 3 : histId, &projectionsPhi[count], &projectionsEta[count], &projectionsPhiSyst[count], &projectionsEtaSyst[count], Ratio);
      count++;
    }
    
    TCanvas* c = new TCanvas(Form("centralities_%d", i), Form("centralities_%d", i), 1200, 600);
    c->Divide(2, 1);
    
    TLegend* legend = new TLegend(0.62, 0.7, 0.88, 0.88);
    legend->SetFillColor(0);

    // syst first
    for (Int_t histId = 0; histId<3; histId++)
    {
      if (Ratio)
      {
        projectionsPhi[count-3+histId]->Divide(projectionsPhi[count-3+1]);
        projectionsEta[count-3+histId]->Divide(projectionsEta[count-3+1]);

      }
      if(!Ratio || (Ratio && histId!=1))
      {
      c->cd(1);
      TH1* clone;
      gPad->SetLeftMargin(0.12);
      projectionsPhiSyst[count-3+histId]->SetFillStyle(3003+histId);
      if (!disableUncertainties) clone = projectionsPhiSyst[count-3+histId]->DrawCopy((histId > 0) ? "E2 ][ SAME" : "E2 ][");
      else clone = projectionsPhi[count-3+histId]->DrawCopy((histId > 0) ? "E ][ SAME" : "E ][");
      clone->SetFillColor(colors[histId]);
      TString label(clone->GetTitle());
      TObjArray* objArray = label.Tokenize("-");
      clone->SetTitle(Form("%s-%s", objArray->At(0)->GetName(), objArray->At(1)->GetName()));
    
      c->cd(2);
      gPad->SetLeftMargin(0.12);
      projectionsEtaSyst[count-3+histId]->SetFillStyle(3003+histId);
      if (!disableUncertainties) clone = projectionsEtaSyst[count-3+histId]->DrawCopy((histId > 0) ? "E2 ][ SAME" : "E2 ][");
      else clone = projectionsEta[count-3+histId]->DrawCopy((histId > 0) ? "E ][ SAME" : "E ][");
      clone->SetFillColor(colors[histId]);
      clone->SetTitle(Form("%s-%s", objArray->At(0)->GetName(), objArray->At(1)->GetName()));
      }
    }
    
    for (Int_t histId = 0; histId<3; histId++)
    {
      if(!Ratio || (Ratio && histId!=1))
      {
      c->cd(1);
      TH1* clone = projectionsPhi[count-3+histId]->DrawCopy("SAME");
      clone->SetLineColor(colors[histId]);
      TString label(clone->GetTitle());
      TObjArray* objArray = label.Tokenize("-");
      legend->AddEntry(clone, (objArray->GetEntries() == 4) ? Form("%s-%s", objArray->At(2)->GetName(), objArray->At(3)->GetName()) : objArray->At(2)->GetName(), "L");
      
      c->cd(2);
      clone = projectionsEta[count-3+histId]->DrawCopy("SAME");
      clone->SetLineColor(colors[histId]);
      }
    }
    
    c->cd(1);
    legend->Draw();
    DrawLatex(0.15, 0.86,  1, "Pb-Pb #sqrt{s_{NN}} = 2.76 TeV", 0.03);
    DrawLatex(0.15, 0.82,  1, "pp #sqrt{s} = 2.76 TeV", 0.03);
    DrawLatex(0.15, 0.78, 1, "|#eta| < 0.9", 0.03);
    DrawLatex(0.15, 0.74, 1, "Projected within |#Delta#eta| < 0.80", 0.03);

    c->cd(2);
    DrawLatex(0.15, 0.86, 1, "Scale uncertainty: 20% (for 0-10%) / 10% (other bins)", 0.03);
    DrawLatex(0.15, 0.82, 1, "Projected within |#Delta#varphi| < 0.87", 0.03);
    if (!disableUncertainties) DrawALICELogo(kTRUE, 0.75, 0.62, 0.95, 0.78);
    
    c->SaveAs(Form("ex/%s.png", c->GetName()));
    c->SaveAs(Form("ex/%s.eps", c->GetName()));
 
//     return;

    if (graphFileName)
    {
      // fit curves
      Int_t graphID = exampleI[i] * (6 - 1) + exampleJ[i] - 1;
      for (Int_t histId = 0; histId<3; histId++)
      {
	c->cd(1);
	// integral over y
// 	TF1* func = new TF1("func", "[0]+[1]*([4]/[2]/TMath::Sqrt(TMath::TwoPi())*exp(-0.5*(x/[2])**2)+(1-[4])/[5]/TMath::Sqrt(TMath::TwoPi())*exp(-0.5*(x/[5])**2))", -0.5 * TMath::Pi(), 0.5 * TMath::Pi());
	// see /home/jgrosseo/limited2dgaussianintegration.nb
	TF1* func = new TF1("func", "[0]+[1]*([4]/[2]/TMath::Sqrt(TMath::TwoPi())*exp(-0.5*(x/[2])**2)*TMath::Erf([7]/TMath::Sqrt(2)/[3])+(1-[4])/[5]/TMath::Sqrt(TMath::TwoPi())*exp(-0.5*(x/[5])**2)*TMath::Erf([7]/TMath::Sqrt(2)/[6]))", -0.5 * TMath::Pi(), 0.5 * TMath::Pi());

	func->SetParameter(0, 0);
/*	if (histId == 0) func->SetParameter(0, -0.00129472);
	if (histId == 1) func->SetParameter(0, -0.00025859);
	if (histId == 2) func->SetParameter(0, -0.000833301);
*/
	for (Int_t k=0; k<6; k++)
	  func->SetParameter(k+1, graphs[k][graphID]->GetY()[(histId==2) ? 3 : histId]);
	// scale by bin width (to compare with projection)
	// NOTE this scaling looks inconsistent, but is OK. See normalization in DrawExample
	func->SetParameter(1, func->GetParameter(1) / projectionsEta[count-3+histId]->GetXaxis()->GetBinWidth(1));
	func->SetParameter(7, 0.8); // project Limit
	func->SetLineColor(colors[histId]);
	func->SetLineWidth(1);
	func->DrawCopy("SAME");
	
	// draw contributions
	Float_t scale = func->GetParameter(1);
	Float_t weighting = func->GetParameter(4);
	func->SetParameter(1, scale * weighting);
	func->SetParameter(4, 1);
	func->SetLineStyle(2);
	if (weighting > 0.01)
	  func->DrawCopy("SAME");

	func->SetParameter(1, scale * (1.0-weighting));
	func->SetParameter(4, 0);
	func->SetLineStyle(3);
	if (1.0-weighting > 0.01)
	  func->DrawCopy("SAME");
	
	c->cd(2);
	// integral over x
	func = new TF1("func", "[0]+[1]*([4]/[3]/TMath::Sqrt(TMath::TwoPi())*exp(-0.5*(x/[3])**2)*TMath::Erf([7]/TMath::Sqrt(2)/[2])+(1-[4])/[6]/TMath::Sqrt(TMath::TwoPi())*exp(-0.5*(x/[6])**2)*TMath::Erf([7]/TMath::Sqrt(2)/[5]))", -kOuterLimit, kOuterLimit);
	func->SetParameter(0, 0);
	for (Int_t k=0; k<6; k++)
	  func->SetParameter(k+1, graphs[k][graphID]->GetY()[(histId==2) ? 3 : histId]);
	// scale by bin width (to compare with projection)
	// NOTE this scaling looks inconsistent, but is OK. See normalization in DrawExample
	func->SetParameter(1, func->GetParameter(1) / projectionsEta[count-3+histId]->GetXaxis()->GetBinWidth(1));
	func->SetParameter(7, 0.87); // project Limit
	func->SetLineColor(colors[histId]);
	func->SetLineWidth(1);
	func->DrawCopy("SAME");

	// draw contributions
	scale = func->GetParameter(1);
	weighting = func->GetParameter(4);
	func->SetParameter(1, scale * weighting);
	func->SetParameter(4, 1);
	func->SetLineStyle(2);
	if (weighting > 0.01)
	  func->DrawCopy("SAME");

	func->SetParameter(1, scale * (1.0-weighting));
	func->SetParameter(4, 0);
	func->SetLineStyle(3);
	if (1.0-weighting > 0.01)
	  func->DrawCopy("SAME");
      }
      
      c->SaveAs(Form("ex/%s_fits.png", c->GetName()));
      c->SaveAs(Form("ex/%s_fits.eps", c->GetName()));
      
//       return;
    }
  }

  return;
  // TODO implement syst uncertainty plotting
  
  for (Int_t histId = 0; histId<3; histId++)
  {
    TCanvas* c = new TCanvas(Form("pt_%d", histId), Form("pt_%d", histId), 1200, 600);
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

      if (!disableUncertainties) DrawALICELogo(kTRUE, 0.65, 0.65, 0.85, 0.85);
    }
    
    c->cd(1);
    legend->Draw();

    c->SaveAs(Form("ex/%s.png", c->GetName()));
    c->SaveAs(Form("ex/%s.eps", c->GetName()));
  }
}

void DrawDoubleHump(const char* histFileName)
{
  Int_t exampleI[] = { 0, 1, 2, 0, 1, 2};
  Int_t exampleJ[] = { 1, 1, 1, 2, 2, 2};
  
  TH1* projectionsPhi[9];
  TH1* projectionsEta[9];
  TH1* projectionsPhiSyst[9];
  TH1* projectionsEtaSyst[9];
  
  TCanvas* c = new TCanvas("DrawDoubleHump", "DrawDoubleHump", 1200, 800);
  c->Divide(3, 2);
  
  TLegend* legend = new TLegend(0.15, 0.7, 0.88, 0.88);
  legend->SetFillColor(0);
  
  for (Int_t i=0; i<6; i++)
  {
    if (i == 3)
      continue;
    
    DrawExample(histFileName, exampleI[i], exampleJ[i], 0, &projectionsPhi[i], &projectionsEta[i], &projectionsPhiSyst[i], &projectionsEtaSyst[i]);

    c->cd(i+1);
    gPad->SetLeftMargin(0.15);
    projectionsPhiSyst[i]->GetXaxis()->SetTitle(Form("%s , %s", "#Delta#varphi (rad.)", "#Delta#eta"));
    projectionsPhiSyst[i]->GetYaxis()->SetTitle(kProjYieldTitlePhiOrEta);
    projectionsPhiSyst[i]->GetYaxis()->SetTitleOffset(1.8);
//     projectionsPhiSyst[i]->SetTitleSize(0.035);
    
    TString label(projectionsPhiSyst[i]->GetTitle());
    TObjArray* objArray = label.Tokenize("-");
    projectionsPhiSyst[i]->SetTitle(Form("%s-%s", objArray->At(0)->GetName(), objArray->At(1)->GetName()));
    
    projectionsPhiSyst[i]->Draw("E2 ][");
    projectionsEtaSyst[i]->Draw("E2 ][ SAME");
    TH1* clone = projectionsPhi[i]->DrawCopy((disableUncertainties) ? "" : "SAME");
    clone->SetLineColor(1);
    clone->GetYaxis()->SetRangeUser(clone->GetMinimum(), clone->GetMaximum() * 1.2);
//     clone->GetXaxis()->SetTitle(Form("%s , %s", clone->GetXaxis()->GetTitle(), "#Delta#eta"));
    clone->GetXaxis()->SetTitle(Form("%s , %s", "#Delta#varphi (rad.)", "#Delta#eta"));
    clone->GetYaxis()->SetTitle(kProjYieldTitlePhiOrEta);
    
    clone = projectionsEta[i]->DrawCopy("SAME");
    clone->SetLineColor(2);

    if (!disableUncertainties)
      DrawLatex(0.3, 0.85, 1, "Scale uncertainty: 20%", 0.04);
    
    DrawLatex(0.3, 0.81, 1, Form("#Delta#varphi projection in |#Delta#eta| < %.2f", 0.8), 0.04);
    DrawLatex(0.3, 0.77, 2, Form("#Delta#eta projection in |#Delta#varphi| < %.2f", 0.87), 0.04);
//     DrawLatex(0.3, 0.75, 1, Form("%s-%s centrality", objArray->At(2)->GetName(), objArray->At(3)->GetName()), 0.04);
    DrawLatex(0.3, 0.73, 1, "0-10% centrality", 0.04);
    
//     return;
  }
  
  if (!disableUncertainties)
    DrawALICELogo(kTRUE, 0.65, 0.5, 0.85, 0.7);
  
  c->SaveAs(Form("ex/%s.png", c->GetName()));
  c->SaveAs(Form("ex/%s.eps", c->GetName()));
}

void DrawFullCentralityDependence(const char* histFileName)
{
  Int_t colors[] = { 1, 2, 5, 4, 3, 6 };
  
  TH1* projectionsPhi[9];
  TH1* projectionsEta[9];
  TH1* projectionsPhiSyst[9];
  TH1* projectionsEtaSyst[9];
  
  Int_t count = 0;
  for (Int_t histId = 0; histId<6; histId++)
  {
    projectionsPhi[count] = 0;
    DrawExample(histFileName, 0, 1, histId, &projectionsPhi[count], &projectionsEta[count], &projectionsPhiSyst[count], &projectionsEtaSyst[count]);
    count++;
  }
  
  TCanvas* c = new TCanvas(Form("centralities_%d", 0), Form("centralities_%d", 0), 800, 400);
  c->Divide(2, 1);
  
  TLegend* legend = new TLegend(0.62, 0.7, 0.88, 0.88);
  legend->SetFillColor(0);
  
  for (Int_t histId = 0; histId<6; histId++)
  {
    if (histId == 2)
      continue;
    
    if (!projectionsPhi[count-6+histId])
      continue;
    
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

void CompareExamples(const char* histFileName1, const char* histFileName2, Int_t i, Int_t j, Int_t histId)
{
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

    
    SubtractEtaGap(hist, kEtaLimit, kOuterLimit, kFALSE);
    
    // new TCanvas; hist->DrawCopy("COLZ");

    TH1* proj1 = 0;
    TH1* proj2 = 0;
    GetProjections(hist, &proj1, &proj2, k);
    
    c->cd(1);
    proj1->SetLineColor(exColors[k]);
    proj1->GetXaxis()->SetRangeUser(-1.49, 1.49);
    proj1->GetYaxis()->SetRangeUser(proj1->GetMinimum() * 1.2, proj1->GetMaximum() * 1.5);
    proj1->DrawCopy((k == 0) ? "" : "SAME");
  
    if (k == 1)
    {
      c2->cd(1);
      proj1->Divide(projFirst[0]);
      proj1->Draw();
      proj1->GetYaxis()->SetRangeUser(0.5, 1.5);
    }
    else
      projFirst[0] = proj1;

    c->cd(2);
    proj2->SetLineColor(exColors[k]);
    proj2->GetXaxis()->SetRangeUser(-1.49, 1.49);
    proj2->GetYaxis()->SetRangeUser(proj2->GetMinimum() * 1.2, proj2->GetMaximum() * 2);
    
    if (0 && k == 1)
    {
      Printf("Scaling!");
      proj2->Scale(projFirst[1]->Integral(projFirst[1]->FindBin(-0.5), projFirst[1]->FindBin(0.5)) / proj2->Integral(projFirst[1]->FindBin(-0.5), projFirst[1]->FindBin(0.5)));
    }
    
    proj2->DrawCopy((k == 0) ? "" : "SAME");

    if (k == 1)
    {
      c2->cd(2);
      proj2->Divide(projFirst[1]);
      proj2->Draw();
      proj2->GetYaxis()->SetRangeUser(0.5, 1.5);
    }
    else
      projFirst[1] = proj2;
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

void DrawEtaGapExample(const char* histFileName, Int_t i = 2, Int_t j = 2, Int_t histId = 0)
{
  Float_t etaLimit = kEtaLimit;
  Float_t outerLimit = kOuterLimit;
  Float_t projectLimit = 0.8;

  TFile::Open(histFileName);

  TH2* hist = (TH2*) gFile->Get(Form("dphi_%d_%d_%d", i, j, histId));
  if (!hist)
    return;
  hist->GetZaxis()->SetTitle(kCorrFuncTitle);
  hist->GetZaxis()->SetTitleOffset(1.8);
  // NOTE fix normalization. these 2d correlations come out of AliUEHist normalized by dphi bin width, but not deta
  hist->Scale(1.0 / hist->GetYaxis()->GetBinWidth(1));
  
  TCanvas* c = new TCanvas("c", "c", 800, 800);
  gPad->SetLeftMargin(0.15);
  TH2* clone = (TH2*) hist->DrawCopy("SURF1");
  clone->SetTitle("");
  clone->GetYaxis()->SetRangeUser(-1.59, 1.59);
  clone->GetXaxis()->SetTitleOffset(1.5);
  clone->GetYaxis()->SetTitleOffset(2);
  clone->SetStats(kFALSE);
//   clone->Rebin2D(2, 2);
//   clone->Scale(0.25);
//   DrawALICELogo(kTRUE, 0.7, 0.7, 0.9, 0.9);
  c->SaveAs("note/etagap_raw.eps");
  c->SaveAs("note/etagap_raw.png");
  
  SubtractEtaGap(hist, etaLimit, outerLimit, kFALSE, kTRUE);
  
  c = new TCanvas("c3", "c3", 800, 800);
  gPad->SetLeftMargin(0.13);
  TH1* proj = hist->ProjectionX(Form("%s_proj1%d", hist->GetName(), 0), hist->GetYaxis()->FindBin(-projectLimit+0.01), hist->GetYaxis()->FindBin(projectLimit-0.01));
  proj->GetXaxis()->SetRangeUser(-1.49, 1.49);
  proj->SetStats(0);
  proj->GetXaxis()->SetTitle(Form("%s , %s", proj->GetXaxis()->GetTitle(), "#Delta#eta"));
  proj->GetYaxis()->SetRangeUser(proj->GetMinimum() * 1.2, proj->GetMaximum() * 1.5);
  proj->GetYaxis()->SetTitle(kProjYieldTitlePhiOrEta);
  proj->GetYaxis()->SetTitleOffset(1.6);
  proj->Draw("");

  proj = hist->ProjectionY(Form("%s_proj2_%d", hist->GetName(), 0), hist->GetXaxis()->FindBin(-projectLimit+0.01), hist->GetXaxis()->FindBin(projectLimit-0.01));
  proj->GetXaxis()->SetRangeUser(-1.49, 1.49);
  proj->GetYaxis()->SetRangeUser(proj->GetMinimum() * 1.2, proj->GetMaximum() * 2);
  proj->SetLineColor(2);
  proj->Draw("SAME");
  DrawLatex(0.2, 0.85, 1, Form("#Delta#varphi projection in |#Delta#eta| < %.2f", projectLimit), 0.04);
  DrawLatex(0.2, 0.80, 2, Form("#Delta#eta projection in |#Delta#varphi| < %.2f", projectLimit), 0.04);
//   DrawALICELogo(kTRUE, 0.7, 0.65, 0.9, 0.85);
  c->SaveAs("note/etagap_subtracted_proj.eps");
  c->SaveAs("note/etagap_subtracted_proj.png");

  c = new TCanvas("c2", "c2", 800, 800);
  gPad->SetLeftMargin(0.15);
  hist->SetTitle("");
  hist->GetYaxis()->SetRangeUser(-1.59, 1.59);
  hist->GetXaxis()->SetTitleOffset(1.5);
  hist->GetYaxis()->SetTitleOffset(2);
  hist->SetStats(kFALSE);
//   hist->Rebin2D(2, 2);
//   hist->Scale(0.25);
  hist->Draw("SURF1");
//   DrawALICELogo(kTRUE, 0.7, 0.7, 0.9, 0.9);
  c->SaveAs("note/etagap_subtracted.eps");
  c->SaveAs("note/etagap_subtracted.png");
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

void DrawResults()
{
  const char* mainFile = "graphs_120511.root";
  const char* wingFile = "graphs_120511_wingremoved.root";
  const char* hijing = "graphs_hijing_120515.root";
  const char* ampt = "graph_ampt_120516.root";
  
  drawLogo = 1;
//   DrawResultsCentrality(mainFile, wingFile); return;
//   MCLabel = "Data  HIJING 1.36 / Pythia 6 P-0"; MCComparison(mainFile, wingFile, hijing, 0); return;

  gROOT->SetBatch(kTRUE);

  const char* folder = "";
  for (Int_t i=0; i<2; i++)
  {
    if (i == 1)
    {
      skipGraphList[0] = 5;
      skipGraphList[1] = 10;
      folder = "reduced/";
    }

    fitLabel = "fit";
    fgFolder.Form("results/%s%s", folder, "results");
    DrawResultsCentrality(mainFile, wingFile);
    
    fitLabel = "2^{nd} fit";
    fgFolder.Form("results/%s%s", folder, "results-method2");
    DrawResultsCentrality(mainFile, wingFile, 16);

    fitLabel = "fit";
    MCLabel = "Data  HIJING 1.36 / Pythia 6 P-0";
    fgFolder.Form("results/%s%s", folder, "hijing");
    MCComparison(mainFile, wingFile, hijing, 0);

    fitLabel = "2^{nd} fit";
    fgFolder.Form("results/%s%s", folder, "hijing-method2");
    MCComparison(mainFile, wingFile, hijing, 0, 16);
    
    fitLabel = "fit";
    MCLabel = "Data  AMPT 2.25 / Pythia 6 P-0";
    fgFolder.Form("results/%s%s", folder, "ampt");
    MCComparison(mainFile, wingFile, ampt, 0);
    
    fitLabel = "2^{nd} fit";
    fgFolder.Form("results/%s%s", folder, "ampt-method2");
    MCComparison(mainFile, wingFile, ampt, 0, 16);
  }
}

void DrawAMPTComparison()
{
  const char* mainFile = "graphs_120511.root";
  const char* wingFile = "graphs_120511_wingremoved.root";
  const char* ampt2 = "ampt_stringmelting_norescattering/graphs.root";
  const char* ampt3 = "ampt_default/graphs.root";
  
  drawLogo = 1;
//   DrawResultsCentrality(mainFile, wingFile); return;
//   MCLabel = "Data  HIJING 1.36 / Pythia 6 P-0"; MCComparison(mainFile, wingFile, hijing, 0); return;

  gROOT->SetBatch(kTRUE);

  skipGraphList[0] = 5;
  skipGraphList[1] = 10;

  MCLabel = "Data  AMPT 2.25 no resc / Pythia 6 P-0";
  fgFolder.Form("results/%s", "ampt-norescattering");
  MCComparison(mainFile, wingFile, ampt2, 0);

  MCLabel = "Data  AMPT 1.25 / Pythia 6 P-0";
  fgFolder.Form("results/%s", "ampt-default");
  MCComparison(mainFile, wingFile, ampt3, 0);
}

void CalculateIAA(const char* GraphFileName, const char* HistFileName, const char* outputFileName = "graphs.root")
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
  WriteGraphs(outputFileName);
}


void CalculateSignificance(const char* HistFileName, const char* outputFileNameSignalBackgroundRatio = "hist1.root",const char* outputFileNameSignificance = "hist2.root")
{
  Int_t ptId = 0;
  Float_t etaLimit = kEtaLimit;
  Float_t outerLimit = kOuterLimit;

  TFile* output1 = TFile::Open(outputFileNameSignalBackgroundRatio, "RECREATE");
  TFile* output2 = TFile::Open(outputFileNameSignificance, "RECREATE");
  TFile* input = TFile::Open(HistFileName);
  Int_t maxLeadingPt = 4;
  Int_t maxAssocPt = 6;

//  TCanvas *c = new TCanvas("c", "", 1400, 1100);
//  c->Divide(2,2);

  for (Int_t i=0; i<maxLeadingPt; i++)
  {
    for (Int_t j=1; j<maxAssocPt; j++)
    {
      // only process when first is filled
      TH2* hist = (TH2*) input->Get(Form("dphi_%d_%d_%d", i, j+1, 0));
      if (!hist)
        continue;
      TH2* signal = (TH2*) hist->Clone("signal");
      TH2* background = (TH2*) hist->Clone("background");
      TH2* significance = (TH2*) hist->Clone("significance");
      TH2* signalBackground = (TH2*) hist->Clone("signalBackground");
      if (hist->GetEntries() < 1e4)
      {
        Printf("%d %d Only %f entries. Skipping...", i, j, hist->GetEntries());
        continue;
      }
      ptId = i*(maxAssocPt-1)+j-1;
      if (SkipGraph(ptId))
        continue;
      for (Int_t histId = 0; histId < NHists; histId++)
      {
        hist = (TH2*) input->Get(Form("dphi_%d_%d_%d", i, j+1, histId));
        if (!hist)
          continue;
        hist->Scale(1.0 / hist->GetYaxis()->GetBinWidth(1));
        signal = (TH2*) hist->Clone("signal");
        background = (TH2*) hist->Clone("background");
        significance = (TH2*) hist->Clone("significance");
        signalBackground = (TH2*) hist->Clone("signalBackground");

        if (hist->GetEntries() < 1e4)
        {
          Printf("%d %d %d Only %f entries. Skipping...", i, j, histId, hist->GetEntries());
          continue;
        }

        if (ptId >= 15)
        {
          etaLimit = 0.5;
          outerLimit = 0.99;
        }
        background = SubtractEtaGap(signal, etaLimit, outerLimit, kFALSE);
        for (Int_t x = 0; x<=hist->GetNbinsX();x++)
        {
          for (Int_t y = 0; y<=hist->GetNbinsY();y++)
          {
              signalBackground->SetBinContent(x,y,signal->GetBinContent(x,y)/background->GetBinContent(x,y));
              significance->SetBinContent(x,y,signal->GetBinContent(x,y)/sqrt(signal->GetBinContent(x,y)+background->GetBinContent(x,y)));

          }
        }
/*        if (i==2 && j==1 && histId==0)
        {
          c->cd(1);
          hist->Draw("colz");
          c->cd(2);
          signal->Draw("colz");
          c->cd(3);
          background->Draw("colz");
          c->cd(4);
          significance->Draw("colz");
        }
*/
        output1->cd();
        signalBackground->Write(Form("dphi_signal_backgroung_ratio_%d_%d_%d", i, j+1, histId));
        output2->cd();
        significance->Write(Form("dphi_significance_%d_%d_%d", i, j+1, histId));
      }
    }
  }
  output1->Close();
  output2->Close();
}

void CompareGraphsDraw(const char* GraphFile1, const char* GraphFile2, Int_t offset = 0)
{
  Int_t nHists = 24; //NHists;

  ReadGraphs(GraphFile2);
  CalculateRMSSigma();
  TGraphErrors*** graphs1 = graphs;

  ReadGraphs(GraphFile1);
  CalculateRMSSigma();

  //used for 2011 data, where there is not enough data in the most periferal bin
/*  for (Int_t ptId=0; ptId<NHists; ptId++)
    for (Int_t i=0; i<graphs[0][ptId]->GetN(); i++)
      if (graphs[0][ptId]->GetX()[i] == 65)
      {
        graphs[1][ptId]->RemovePoint(i);
        graphs[2][ptId]->RemovePoint(i);
        graphs1[1][ptId]->RemovePoint(i);
        graphs1[2][ptId]->RemovePoint(i);
        cout << "Removing point" << endl;
      }
*/
  DrawCentrality("phi_rms_centrality_compare", nHists, graphs[1+offset], 0, 0.8, Form("#sigma_{#Delta#varphi} (%s) (rad.)", fitLabel), graphs1[1+offset]);
  DrawCentrality("eta_rms_centrality_compare", nHists, graphs[2+offset], 0, 0.8, Form("#sigma_{#Delta#eta} (%s)", fitLabel), graphs1[2+offset]);

  DrawCentrality("kurtosisphi_centrality_compare", 12, graphs[14+offset], -1.5, 2.5, "Kurtosis #Delta#varphi", graphs1[14+offset]);
  DrawCentrality("kurtosiseta_centrality_compare", 12, graphs[15+offset], -1.5, 2.5, "Kurtosis #Delta#eta", graphs1[15+offset]);

  DrawCentrality("chi2_phi(1D)", nHists, graphs[43], 0, 140, "#chi^{2}/ndf", graphs1[43]);
  DrawCentrality("chi2_eta(1D)", nHists, graphs[44], 0, 140, "#chi^{2}/ndf", graphs1[44]);

//  DrawCentrality("IAAHist", nHists, graphs[34+offset], 0, 20, "I_{AA} (Hist)", graphs1[34+offset]);
//  DrawCentrality("IAAFit", nHists, graphs[32+offset], 0, 20, "I_{AA} (Fit)", graphs1[32+offset]);
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

  TCanvas *c = new TCanvas(Form("c%d", canvasCount), Form("c%d", canvasCount), 600, 900);
  c->Divide(1,3);
  canvasCount++;

  TH2* ratio = (TH2*) hist1->Clone("ratio");
  if (hist1->GetEntries() < 1e4 || hist2->GetEntries() < 1e4)
  {
    cout << "Hist1: " << hist1->GetEntries() << " entries" << endl << "Hist2: " << hist2->GetEntries() << " entries" << endl;
  }
  ptId = i*(maxAssocPt-1)+j-1;

  if (0) //SubtractEtaGap
  {
    Float_t etaLimit = kEtaLimit;
    Float_t outerLimit = kOuterLimit;
    if (ptId >= 15)
    {
      etaLimit = 0.5;
      outerLimit = 0.99;
    }
    SubtractEtaGap(hist1, etaLimit, outerLimit, kFALSE);
    SubtractEtaGap(hist2, etaLimit, outerLimit, kFALSE);
  }

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
  
  hist1->GetYaxis()->SetRangeUser(-1.99, 1.99);
  hist2->GetYaxis()->SetRangeUser(-1.99, 1.99);
  ratio->GetYaxis()->SetRangeUser(-1.99, 1.99);
  
  c->cd(1);
  hist1->Draw("surf1");
  c->cd(2);
  hist2->Draw("surf1");
  c->cd(3);
  ratio->Draw("colz");
  
  return;
  
  new TCanvas;
  hist1->ProjectionX("p1", hist1->GetYaxis()->FindBin(-1.5), hist1->GetYaxis()->FindBin(1.5))->Draw();
  hist2->ProjectionY("p2", hist2->GetXaxis()->FindBin(-1.5), hist2->GetXaxis()->FindBin(1.5))->DrawCopy("SAME")->SetLineColor(2);
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
  const char* title[8] = {"PYTHIA 62","Cu+Cu 62 0-60%","Au+Au 62 0-80%","PYTHIA 200","d+Au 200 0-95%","Cu+Cu 200 0-60%","Au+Au 200 40-80%","Au+Au 200 0-12%"};

  CreateGraphStructure();
  TFile* GraphFile = TFile::Open(GraphFileName);
  for (Int_t i=0; i<NGraphs; i++)
    for (Int_t j=0; j<NHists; j++)
    {
      if (i == 32 || i == 33 || i == 34) continue;
      graphs[i][j] = (TGraphErrors*) GraphFile->Get(Form("graph_%d_%d", i, j));
    }

  TGraphErrors*** STARgraphs = new TGraphErrors**[2];
  for (Int_t i=0; i<4; i++)
  {
    STARgraphs[i] = new TGraphErrors*[8];
    for (Int_t j=0; j<8; j++)
      STARgraphs[i][j] = new TGraphErrors;
  }

  ifstream infile(STARFileName);

  TLegend* legend = new TLegend(0.55, 0.65, 0.95, 0.95);
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
          STARgraphs[0][j]->SetMarkerStyle(24);
          STARgraphs[0][j]->SetMarkerColor(2);
          STARgraphs[0][j]->SetLineColor(2);
          STARgraphs[0][j]->SetMarkerSize(1.5);
        }
        else 
        {
          STARgraphs[0][j]->SetMarkerStyle(20+j);
          STARgraphs[0][j]->SetMarkerColor(j);
          STARgraphs[0][j]->SetLineColor(j);
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
          STARgraphs[2][j]->SetMarkerStyle(24);
          STARgraphs[2][j]->SetMarkerColor(2);
          STARgraphs[2][j]->SetLineColor(2);
          STARgraphs[2][j]->SetMarkerSize(1.5);
        }
        else 
        {
          STARgraphs[2][j]->SetMarkerStyle(20+j);
          STARgraphs[2][j]->SetMarkerColor(j);
          STARgraphs[2][j]->SetLineColor(j);
          STARgraphs[2][j]->SetMarkerSize(1.5);
        }
        if (j==7 || j==6) mg2->Add(STARgraphs[2][j]);
        STARgraphs[2][j]->SetTitle(title[j]);
      }
    }
  }

  Float_t pTaAxisMapping[] = { 1.25, 1.75, 2.25, 2.75 };
  Float_t pTaAxisMappingE[] = { 0.25, 0.25, 0.25, 0.25 };

  CalculateRMSSigma();

//  for 1D fits:
  Int_t graphPhi = 37;
  Int_t graphEta = 41;
//  for 2D fits:
//  Int_t graphPhi = 1;
//  Int_t graphEta = 2;

  Int_t ptId = 0;
  for (ptId = 0; ptId < 4; ptId++)
  {
    AddPoint(STARgraphs[1][0],pTaAxisMapping[ptId],graphs[graphPhi][ptId]->GetY()[0],pTaAxisMappingE[ptId],graphs[graphPhi][ptId]->GetErrorY(0));
    STARgraphs[1][0]->SetMarkerStyle(20);
    STARgraphs[1][0]->SetMarkerColor(3);
    STARgraphs[1][0]->SetLineColor(3);
    STARgraphs[1][0]->SetMarkerSize(1.5);
    STARgraphs[1][0]->SetTitle(graphs[graphPhi][ptId]->GetTitle());
    AddPoint(STARgraphs[1][1],pTaAxisMapping[ptId],graphs[graphPhi][ptId]->GetY()[1],pTaAxisMappingE[ptId],graphs[graphPhi][ptId]->GetErrorY(1));
    STARgraphs[1][1]->SetMarkerStyle(22);
    STARgraphs[1][1]->SetMarkerColor(4);
    STARgraphs[1][1]->SetLineColor(4);
    STARgraphs[1][1]->SetMarkerSize(1.5);
    STARgraphs[1][1]->SetTitle(graphs[graphPhi][ptId]->GetTitle());
    AddPoint(STARgraphs[1][2],pTaAxisMapping[ptId],graphs[graphPhi][ptId]->GetY()[5],pTaAxisMappingE[ptId],graphs[graphPhi][ptId]->GetErrorY(5));
    STARgraphs[1][2]->SetMarkerStyle(23);
    STARgraphs[1][2]->SetMarkerColor(1);
    STARgraphs[1][2]->SetLineColor(1);
    STARgraphs[1][2]->SetMarkerSize(1.5);
    STARgraphs[1][2]->SetTitle(graphs[graphPhi][ptId]->GetTitle());

    AddPoint(STARgraphs[3][0],pTaAxisMapping[ptId],graphs[graphEta][ptId]->GetY()[0],pTaAxisMappingE[ptId],graphs[graphEta][ptId]->GetErrorY(0));
    STARgraphs[3][0]->SetMarkerStyle(20);
    STARgraphs[3][0]->SetMarkerColor(3);
    STARgraphs[3][0]->SetLineColor(3);
    STARgraphs[3][0]->SetMarkerSize(1.5);
    STARgraphs[3][0]->SetTitle(graphs[graphEta][ptId]->GetTitle());
    AddPoint(STARgraphs[3][1],pTaAxisMapping[ptId],graphs[graphEta][ptId]->GetY()[1],pTaAxisMappingE[ptId],graphs[graphEta][ptId]->GetErrorY(1));
    STARgraphs[3][1]->SetMarkerStyle(22);
    STARgraphs[3][1]->SetMarkerColor(4);
    STARgraphs[3][1]->SetLineColor(4);
    STARgraphs[3][1]->SetMarkerSize(1.5);
    STARgraphs[3][1]->SetTitle(graphs[graphEta][ptId]->GetTitle());
    AddPoint(STARgraphs[3][2],pTaAxisMapping[ptId],graphs[graphEta][ptId]->GetY()[5],pTaAxisMappingE[ptId],graphs[graphEta][ptId]->GetErrorY(5));
    STARgraphs[3][2]->SetMarkerStyle(23);
    STARgraphs[3][2]->SetMarkerColor(1);
    STARgraphs[3][2]->SetLineColor(1);
    STARgraphs[3][2]->SetMarkerSize(1.5);
    STARgraphs[3][2]->SetTitle(graphs[graphEta][ptId]->GetTitle());
  }
  mg1->Add(STARgraphs[1][0]);
  mg1->Add(STARgraphs[1][1]);
  mg1->Add(STARgraphs[1][2]);

  mg2->Add(STARgraphs[3][0]);
  mg2->Add(STARgraphs[3][1]);
  mg2->Add(STARgraphs[3][2]);

  legend->AddEntry(mg1->GetListOfGraphs()->At(2), "3.0 < P_Tt < 6.0 GeV, 0-10%", "P");
  legend->AddEntry(mg1->GetListOfGraphs()->At(4), "3.0 < P_Tt < 6.0 GeV, 40-60%", "P");
  legend->AddEntry(mg1->GetListOfGraphs()->At(3), "3.0 < P_Tt < 6.0 GeV, 60-70%", "P");

  for (Int_t j=0; j<8;j++)
  {
    if (j==6) legend->AddEntry(mg1->GetListOfGraphs()->At(0), STARgraphs[0][j]->GetTitle(), "P");
    if (j==7) legend->AddEntry(mg1->GetListOfGraphs()->At(1), STARgraphs[0][j]->GetTitle(), "P");
  }
  c1->cd();
  mg1->Draw("AP");
  mg1->GetXaxis()->SetTitle("p_{T,a} (GeV)");
  mg1->GetYaxis()->SetTitle(Form("#sigma_{#Delta#varphi} (%s) (rad.)", fitLabel));
  legend->Draw();

  c2->cd();
  mg2->Draw("AP");
  mg2->GetXaxis()->SetTitle("p_{T,a} (GeV)");
  mg2->GetYaxis()->SetTitle(Form("#sigma_{#Delta#eta} (%s) (rad.)", fitLabel));
  legend->Draw();
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
  const char* title[8] = {"PYTHIA 62","Cu+Cu 62 0-60%","Au+Au 62 0-80%","PYTHIA 200","d+Au 200 0-95%","Cu+Cu 200 0-60%","Au+Au 200 40-80%","Au+Au 200 0-12%"};

  CreateGraphStructure();
  TFile* GraphFile = TFile::Open(GraphFileName);
  for (Int_t i=0; i<NGraphs; i++)
    for (Int_t j=0; j<NHists; j++)
    {
      if (i == 32 || i == 33 || i == 34) continue;
      graphs[i][j] = (TGraphErrors*) GraphFile->Get(Form("graph_%d_%d", i, j));
    }

  TGraphErrors*** STARgraphs = new TGraphErrors**[2];
  for (Int_t i=0; i<4; i++)
  {
    STARgraphs[i] = new TGraphErrors*[8];
    for (Int_t j=0; j<8; j++)
      STARgraphs[i][j] = new TGraphErrors;
  }

  ifstream infile(STARFileName);

  TLegend* legend = new TLegend(0.55, 0.65, 0.95, 0.95);
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
          STARgraphs[0][j]->SetMarkerColor(j);
          STARgraphs[0][j]->SetLineColor(j);
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
          STARgraphs[2][j]->SetMarkerColor(j);
          STARgraphs[2][j]->SetLineColor(j);
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

  Float_t pTtAxisMapping[] = { 2.5, 3.5, 4.5, 5.5 };
  Float_t pTtAxisMappingE[] = { 0.5, 0.5, 0.5, 0.5 };

  CalculateRMSSigma();

//  for 1D fits:
  Int_t graphPhi = 37;
  Int_t graphEta = 41;
//  for 2D fits:
//  Int_t graphPhi = 1;
//  Int_t graphEta = 2;

  for (Int_t i = 0; i < 4; i++)
  {
    Int_t ptId = i * (6 - 1) + 1 - 1;
    AddPoint(STARgraphs[1][0],pTtAxisMapping[i],graphs[graphPhi][ptId]->GetY()[0],pTtAxisMappingE[i],graphs[graphPhi][ptId]->GetErrorY(0));
    STARgraphs[1][0]->SetMarkerStyle(20);
    STARgraphs[1][0]->SetMarkerColor(3);
    STARgraphs[1][0]->SetLineColor(3);
    STARgraphs[1][0]->SetMarkerSize(1.5);
    STARgraphs[1][0]->SetTitle(graphs[graphPhi][ptId]->GetTitle());
    AddPoint(STARgraphs[1][1],pTtAxisMapping[i],graphs[graphPhi][ptId]->GetY()[1],pTtAxisMappingE[i],graphs[graphPhi][ptId]->GetErrorY(1));
    STARgraphs[1][1]->SetMarkerStyle(22);
    STARgraphs[1][1]->SetMarkerColor(4);
    STARgraphs[1][1]->SetLineColor(4);
    STARgraphs[1][1]->SetMarkerSize(1.5);
    STARgraphs[1][1]->SetTitle(graphs[graphPhi][ptId]->GetTitle());
    AddPoint(STARgraphs[1][2],pTtAxisMapping[i],graphs[graphPhi][ptId]->GetY()[5],pTtAxisMappingE[i],graphs[graphPhi][ptId]->GetErrorY(5));
    STARgraphs[1][2]->SetMarkerStyle(23);
    STARgraphs[1][2]->SetMarkerColor(1);
    STARgraphs[1][2]->SetLineColor(1);
    STARgraphs[1][2]->SetMarkerSize(1.5);
    STARgraphs[1][2]->SetTitle(graphs[graphPhi][ptId]->GetTitle());

    AddPoint(STARgraphs[3][0],pTtAxisMapping[i],graphs[graphEta][ptId]->GetY()[0],pTtAxisMappingE[i],graphs[graphEta][ptId]->GetErrorY(0));
    STARgraphs[3][0]->SetMarkerStyle(20);
    STARgraphs[3][0]->SetMarkerColor(3);
    STARgraphs[3][0]->SetLineColor(3);
    STARgraphs[3][0]->SetMarkerSize(1.5);
    STARgraphs[3][0]->SetTitle(graphs[graphEta][ptId]->GetTitle());
    AddPoint(STARgraphs[3][1],pTtAxisMapping[i],graphs[graphEta][ptId]->GetY()[1],pTtAxisMappingE[i],graphs[graphEta][ptId]->GetErrorY(1));
    STARgraphs[3][1]->SetMarkerStyle(22);
    STARgraphs[3][1]->SetMarkerColor(4);
    STARgraphs[3][1]->SetLineColor(4);
    STARgraphs[3][1]->SetMarkerSize(1.5);
    STARgraphs[3][1]->SetTitle(graphs[graphEta][ptId]->GetTitle());
    AddPoint(STARgraphs[3][2],pTtAxisMapping[i],graphs[graphEta][ptId]->GetY()[5],pTtAxisMappingE[i],graphs[graphEta][ptId]->GetErrorY(5));
    STARgraphs[3][2]->SetMarkerStyle(23);
    STARgraphs[3][2]->SetMarkerColor(1);
    STARgraphs[3][2]->SetLineColor(1);
    STARgraphs[3][2]->SetMarkerSize(1.5);
    STARgraphs[3][2]->SetTitle(graphs[graphEta][ptId]->GetTitle());
   }

  mg1->Add(STARgraphs[1][0]);
  mg1->Add(STARgraphs[1][1]);
  mg1->Add(STARgraphs[1][2]);

  mg2->Add(STARgraphs[3][0]);
  mg2->Add(STARgraphs[3][1]);
  mg2->Add(STARgraphs[3][2]);

  legend->AddEntry(mg1->GetListOfGraphs()->At(2), "1.5 GeV < p_{T,a} < p_{T,t}, 0-10%", "P");
  legend->AddEntry(mg1->GetListOfGraphs()->At(4), "1.5 GeV < p_{T,a} < p_{T,t}, 40-60%", "P");
  legend->AddEntry(mg1->GetListOfGraphs()->At(3), "1.5 GeV < p_{T,a} < p_{T,t}, 60-70%", "P");

  for (Int_t j=0; j<8;j++)
  {
    if (j==6) legend->AddEntry(mg1->GetListOfGraphs()->At(0), STARgraphs[0][j]->GetTitle(), "P");
    if (j==7) legend->AddEntry(mg1->GetListOfGraphs()->At(1), STARgraphs[0][j]->GetTitle(), "P");
  }
  c1->cd();
  mg1->Draw("AP");
  mg1->GetXaxis()->SetTitle("p_{T,t} (GeV)");
  mg1->GetYaxis()->SetTitle(Form("#sigma_{#Delta#varphi} (%s) (rad.)", fitLabel));
  legend->Draw();

  c2->cd();
  mg2->Draw("AP");
  mg2->GetXaxis()->SetTitle("p_{T,t} (GeV)");
  mg2->GetYaxis()->SetTitle(Form("#sigma_{#Delta#eta} (%s) (rad.)", fitLabel));
  legend->Draw();
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
  
  for (Int_t x=1; x<=phi->GetNbinsX(); x++)
  {
    Double_t value = 0;
    for (Int_t y=1; y<=phi->GetNbinsX(); y++)
    {
      Int_t bin = x - y;
      if (bin < 1)
	bin += phi->GetNbinsX();
      if (bin > phi->GetNbinsX())
	bin -= phi->GetNbinsX();
      value += phi->GetBinContent(y) * phi->GetBinContent(bin);
    }
    conv->SetBinContent(x, value);
  }
  
  conv->Scale(1.0 / conv->Integral());
  new TCanvas; conv->Draw();
  
  return;
  
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
