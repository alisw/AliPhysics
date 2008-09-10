/* $Id$ */

//
// plots for the note about multplicity measurements
//

#if !defined(__CINT__) || defined(__MAKECINT__)

#include <TCanvas.h>
#include <TPad.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TLine.h>
#include <TF1.h>
#include <TSystem.h>
#include <TFile.h>
#include <TLegend.h>
#include <TStopwatch.h>
#include <TROOT.h>
#include <TGraph.h>
#include <TMath.h>
#include <TPaveText.h>
#include <TImage.h>
#include <TLatex.h>

#include "AliMultiplicityCorrection.h"
#include "AliCorrection.h"
#include "AliCorrectionMatrix3D.h"

#endif

const char* correctionFile = "multiplicity.root";
const char* measuredFile   = "multiplicityMC_1M_3.root";
Int_t etaRange = 2;
Int_t displayRange = 200; // axis range
Int_t ratioRange = 151;   // range to calculate difference
Int_t longDisplayRange = 200;

const char* correctionFileTPC = "multiplicityMC_TPC_1.4M.root";
const char* measuredFileTPC   = "multiplicityMC_TPC_0.6M.root";
Int_t etaRangeTPC = 1;

void loadlibs()
{
  gSystem->Load("libANALYSIS");
  gSystem->Load("libPWG0base");
}

void SetTPC()
{
  correctionFile = correctionFileTPC;
  measuredFile = measuredFileTPC;
  etaRange = etaRangeTPC;
  displayRange = 100;
  ratioRange = 76;
  longDisplayRange = 100;
}

void Smooth(TH1* hist, Int_t windowWidth = 20)
{
  TH1* clone = (TH1*) hist->Clone("clone");
  for (Int_t bin=2; bin<=clone->GetNbinsX(); ++bin)
  {
    Int_t min = TMath::Max(2, bin-windowWidth);
    Int_t max = TMath::Min(clone->GetNbinsX(), bin+windowWidth);
    Float_t average = clone->Integral(min, max) / (max - min + 1);

    hist->SetBinContent(bin, average);
    hist->SetBinError(bin, 0);
  }

  delete clone;
}

void responseMatrixPlot()
{
  gSystem->Load("libPWG0base");

  AliMultiplicityCorrection* mult = new AliMultiplicityCorrection("Multiplicity", "Multiplicity");

  TFile::Open(correctionFile);
  mult->LoadHistograms("Multiplicity");

  TH1* hist = mult->GetCorrelation(etaRange)->Project3D("zy");
  hist->SetStats(kFALSE);

  hist->SetTitle(";true multiplicity;measured multiplicity;Entries");
  hist->GetXaxis()->SetRangeUser(0, longDisplayRange);
  hist->GetYaxis()->SetRangeUser(0, longDisplayRange);

  TCanvas* canvas = new TCanvas("c1", "c1", 800, 600);
  canvas->SetRightMargin(0.15);
  canvas->SetTopMargin(0.05);

  gPad->SetLogz();
  hist->Draw("COLZ");

  canvas->SaveAs("responsematrix.eps");
}

TCanvas* DrawResultRatio(TH1* mcHist, TH1* result, TString epsName)
{
  // normalize unfolded result to mc hist
  result->Scale(1.0 / result->Integral(2, 200));
  result->Scale(mcHist->Integral(2, 200));

  TCanvas* canvas = new TCanvas(epsName, epsName, 800, 600);
  canvas->Range(0, 0, 1, 1);

  TPad* pad1 = new TPad(Form("%s_pad1", epsName.Data()), "", 0, 0.5, 0.98, 0.98);
  pad1->Draw();

  TPad* pad2 = new TPad(Form("%s_pad2", epsName.Data()), "", 0, 0.02, 0.98, 0.5);
  pad2->Draw();

  pad1->SetRightMargin(0.05);
  pad2->SetRightMargin(0.05);

  // no border between them
  pad1->SetBottomMargin(0);
  pad2->SetTopMargin(0);

  pad1->cd();

  mcHist->GetXaxis()->SetLabelSize(0.06);
  mcHist->GetYaxis()->SetLabelSize(0.06);
  mcHist->GetXaxis()->SetTitleSize(0.06);
  mcHist->GetYaxis()->SetTitleSize(0.06);
  mcHist->GetYaxis()->SetTitleOffset(0.6);

  mcHist->GetXaxis()->SetRangeUser(0, displayRange);

  mcHist->SetTitle(";true multiplicity;Entries");
  mcHist->SetStats(kFALSE);

  mcHist->DrawCopy("HIST E");
  gPad->SetLogy();

  result->SetLineColor(2);
  result->DrawCopy("SAME HISTE");

  TLegend* legend = new TLegend(0.6, 0.65, 0.95, 0.9);
  legend->AddEntry(mcHist, "true distribution");
  legend->AddEntry(result, "unfolded distribution");
  legend->SetFillColor(0);
  legend->Draw();

  pad2->cd();
  pad2->SetBottomMargin(0.15);

  // calculate ratio
  mcHist->Sumw2();
  TH1* ratio = (TH1*) mcHist->Clone("ratio");
  result->Sumw2();
  ratio->Divide(ratio, result, 1, 1, "");
  ratio->GetYaxis()->SetTitle("Ratio (true / unfolded)");
  ratio->GetYaxis()->SetRangeUser(0.55, 1.45);

  ratio->DrawCopy();

  // get average of ratio
  Float_t sum = 0;
  for (Int_t i=2; i<=ratioRange; ++i)
  {
    sum += TMath::Abs(ratio->GetBinContent(i) - 1);
  }
  sum /= ratioRange-1;

  printf("Average (2..%d) of |ratio - 1| is %f\n", ratioRange, sum);

  TLine* line = new TLine(0, 1, displayRange, 1);
  line->SetLineWidth(2);
  line->Draw();

  line = new TLine(0, 1.1, displayRange, 1.1);
  line->SetLineWidth(2);
  line->SetLineStyle(2);
  line->Draw();
  line = new TLine(0, 0.9, displayRange, 0.9);
  line->SetLineWidth(2);
  line->SetLineStyle(2);
  line->Draw();

  canvas->Modified();

  canvas->SaveAs(epsName);

  return canvas;
}

TCanvas* Draw2ResultRatio(TH1* mcHist, TH1* result1, TH1* result2, TString epsName)
{
  // draws the 3 plots in the upper plot
  // draws the ratio between result1 and result2 in the lower plot

  // normalize unfolded result to mc hist
  result1->Scale(1.0 / result1->Integral(2, 200));
  result1->Scale(mcHist->Integral(2, 200));
  result2->Scale(1.0 / result2->Integral(2, 200));
  result2->Scale(mcHist->Integral(2, 200));

  TCanvas* canvas = new TCanvas(epsName, epsName, 800, 600);
  canvas->Range(0, 0, 1, 1);

  TPad* pad1 = new TPad(Form("%s_pad1", epsName.Data()), "", 0, 0.5, 0.98, 0.98);
  pad1->Draw();

  TPad* pad2 = new TPad(Form("%s_pad2", epsName.Data()), "", 0, 0.02, 0.98, 0.5);
  pad2->Draw();

  pad1->SetRightMargin(0.05);
  pad2->SetRightMargin(0.05);

  // no border between them
  pad1->SetBottomMargin(0);
  pad2->SetTopMargin(0);

  pad1->cd();

  mcHist->GetXaxis()->SetLabelSize(0.06);
  mcHist->GetYaxis()->SetLabelSize(0.06);
  mcHist->GetXaxis()->SetTitleSize(0.06);
  mcHist->GetYaxis()->SetTitleSize(0.06);
  mcHist->GetYaxis()->SetTitleOffset(0.6);

  mcHist->GetXaxis()->SetRangeUser(0, displayRange);

  mcHist->SetTitle(";true multiplicity;Entries");
  mcHist->SetStats(kFALSE);

  mcHist->DrawCopy("HIST E");
  gPad->SetLogy();

  result1->SetLineColor(2);
  result1->DrawCopy("SAME HISTE");

  result2->SetLineColor(4);
  result2->DrawCopy("SAME HISTE");

  TLegend* legend = new TLegend(0.55, 0.6, 0.95, 0.9);
  legend->AddEntry(mcHist, "true distribution");
  legend->AddEntry(result1, "unfolded distribution (syst)");
  legend->AddEntry(result2, "unfolded distribution (normal)");
  legend->SetFillColor(0);
  legend->Draw();

  pad2->cd();
  pad2->SetBottomMargin(0.15);

  result1->GetXaxis()->SetLabelSize(0.06);
  result1->GetYaxis()->SetLabelSize(0.06);
  result1->GetXaxis()->SetTitleSize(0.06);
  result1->GetYaxis()->SetTitleSize(0.06);
  result1->GetYaxis()->SetTitleOffset(0.6);

  result1->GetXaxis()->SetRangeUser(0, displayRange);

  result1->SetTitle(";true multiplicity;Entries");
  result1->SetStats(kFALSE);

  // calculate ratio
  result1->Sumw2();
  TH1* ratio = (TH1*) result1->Clone("ratio");
  result2->Sumw2();
  ratio->Divide(ratio, result2, 1, 1, "");
  ratio->GetYaxis()->SetTitle("Ratio (syst / normal)");
  ratio->GetYaxis()->SetRangeUser(0.55, 1.45);

  ratio->DrawCopy();

  // get average of ratio
  Float_t sum = 0;
  for (Int_t i=2; i<=ratioRange; ++i)
  {
    sum += TMath::Abs(ratio->GetBinContent(i) - 1);
  }
  sum /= ratioRange-1;

  printf("Average (2..%d) of |ratio - 1| is %f\n", ratioRange, sum);

  TLine* line = new TLine(0, 1, displayRange, 1);
  line->SetLineWidth(2);
  line->Draw();

  line = new TLine(0, 1.1, displayRange, 1.1);
  line->SetLineWidth(2);
  line->SetLineStyle(2);
  line->Draw();
  line = new TLine(0, 0.9, displayRange, 0.9);
  line->SetLineWidth(2);
  line->SetLineStyle(2);
  line->Draw();

  canvas->Modified();

  canvas->SaveAs(epsName);

  return canvas;
}

TCanvas* DrawRatio(TH1* result, Int_t nResultSyst, TH1** resultSyst, TString epsName, Bool_t firstMarker = kFALSE, const char** legendStrings = 0, Bool_t errors = kFALSE)
{
  // compares n results with first results. E.g. one gained with the default response, another with a changed one to study
  // a systematic effect

  // normalize results
  result->Scale(1.0 / result->Integral(2, 200));

  TCanvas* canvas = new TCanvas(epsName, epsName, 800, 400);
  canvas->SetTopMargin(0.05);
  canvas->SetRightMargin(0.05);

  result->GetXaxis()->SetRangeUser(0, displayRange);
  result->GetYaxis()->SetRangeUser(0.55, 1.45);
  result->SetStats(kFALSE);

  // to get the axis how we want it
  TH1* dummy = (TH1*) result->Clone("dummy");
  dummy->Reset();
  dummy->SetTitle(";true multiplicity;Ratio");
  dummy->DrawCopy();
  delete dummy;

  Int_t colors[] = {1, 2, 4, 6, 7, 8, 9, 10};

  TLegend* legend = new TLegend(0.2, 0.75, 0.35, 0.95);
  legend->SetFillColor(0);

  for (Int_t n=0; n<nResultSyst; ++n)
  {
    resultSyst[n]->Scale(1.0 / resultSyst[n]->Integral(2, 200));

    // calculate ratio
    TH1* ratio = (TH1*) result->Clone("ratio");
    ratio->Divide(ratio, resultSyst[n], 1, 1, "");
    ratio->GetXaxis()->SetRangeUser(1, displayRange);

    if (firstMarker)
      ratio->SetMarkerStyle(5);

    ratio->SetLineColor(colors[n / 2]);
    if ((n % 2))
      ratio->SetLineStyle(2);

    TString drawStr("SAME HIST");
    if (n == 0 && firstMarker)
      drawStr = "SAME P";
    if (errors)
      drawStr += " E";

    ratio->DrawCopy(drawStr);

    if (legendStrings && legendStrings[n])
      legend->AddEntry(ratio, legendStrings[n]);

    // get average of ratio
    Float_t sum = 0;
    for (Int_t i=2; i<=ratioRange; ++i)
      sum += TMath::Abs(ratio->GetBinContent(i) - 1);
    sum /= ratioRange-1;

    printf("%d) Average (2..%d) of |ratio - 1| is %f\n", n, ratioRange, sum);
  }

  if (legendStrings)
    legend->Draw();

  TLine* line = new TLine(0, 1, displayRange, 1);
  line->SetLineWidth(2);
  line->Draw();

  line = new TLine(0, 1.1, displayRange, 1.1);
  line->SetLineWidth(2);
  line->SetLineStyle(2);
  line->Draw();
  line = new TLine(0, 0.9, displayRange, 0.9);
  line->SetLineWidth(2);
  line->SetLineStyle(2);
  line->Draw();

  canvas->SaveAs(epsName);
  canvas->SaveAs(Form("%s.gif", epsName.Data()));

  return canvas;
}

TCanvas* DrawRatio(Int_t nResultSyst, TH1** mc, TH1** result, TString epsName, Bool_t smooth = kFALSE, Bool_t dashed = kFALSE)
{
  // draws the ratios of each mc to the corresponding result

  TCanvas* canvas = new TCanvas(epsName, epsName, 800, 400);
  canvas->SetRightMargin(0.05);
  canvas->SetTopMargin(0.05);

  for (Int_t n=0; n<nResultSyst; ++n)
  {
    // normalize
    result[n]->Scale(1.0 / result[n]->Integral(2, 200));
    mc[n]->Scale(1.0 / mc[n]->Integral(2, 200));

    result[n]->GetXaxis()->SetRangeUser(0, displayRange);
    result[n]->SetStats(kFALSE);

    // calculate ratio
    TH1* ratio = (TH1*) result[n]->Clone("ratio");
    ratio->Divide(mc[n], ratio, 1, 1, "B");

    // SetRangeUser(1, ...) would be the same, but the 0 should be still on the axis...
    ratio->SetBinContent(1, 1); ratio->SetBinError(1, 0);

    if (smooth)
      Smooth(ratio);

    ratio->SetTitle(Form(";true multiplicity;Ratio (true / unfolded)%s", ((smooth) ? " (smoothed)" : "")));
    ratio->GetYaxis()->SetRangeUser(0.55, 1.45);

    if (dashed)
    {
      ratio->SetLineColor((n/2)+1);
      ratio->SetLineStyle((n%2)+1);
    }
    else
      ratio->SetLineColor(n+1);

    ratio->DrawCopy((n == 0) ? "HIST" : "SAME HIST");

    // get average of ratio
    Float_t sum = 0;
    for (Int_t i=2; i<=ratioRange; ++i)
      sum += TMath::Abs(ratio->GetBinContent(i) - 1);
    sum /= ratioRange-1;

    printf("%d) Average (2..%d) of |ratio - 1| is %f\n", n, ratioRange, sum);
  }

  TLine* line = new TLine(0, 1, displayRange, 1);
  line->SetLineWidth(2);
  line->Draw();

  line = new TLine(0, 1.1, displayRange, 1.1);
  line->SetLineWidth(2);
  line->SetLineStyle(2);
  line->Draw();
  line = new TLine(0, 0.9, displayRange, 0.9);
  line->SetLineWidth(2);
  line->SetLineStyle(2);
  line->Draw();

  canvas->Modified();

  canvas->SaveAs(epsName);
  canvas->SaveAs(Form("%s.gif", epsName.Data()));

  return canvas;
}

TCanvas* DrawRatioDeduct(TH1* mcBase, TH1* resultBase, Int_t nResultSyst, TH1** mc, TH1** result, TString epsName)
{
  // draws the ratios of each mc to the corresponding result
  // deducts from each ratio the ratio of mcBase / resultBase

  // normalize
  resultBase->Scale(1.0 / resultBase->Integral(2, 200));
  mcBase->Scale(1.0 / mcBase->Integral(2, 200));

  // calculate ratio
  TH1* ratioBase = (TH1*) resultBase->Clone("ratioBase");
  ratioBase->Divide(mcBase, ratioBase, 1, 1, "B");

  TCanvas* canvas = new TCanvas(epsName, epsName, 800, 400);
  canvas->SetRightMargin(0.05);
  canvas->SetTopMargin(0.05);

  for (Int_t n=0; n<nResultSyst; ++n)
  {
    // normalize
    result[n]->Scale(1.0 / result[n]->Integral(2, 200));
    mc[n]->Scale(1.0 / mc[n]->Integral(2, 200));

    result[n]->GetXaxis()->SetRangeUser(0, displayRange);
    result[n]->SetStats(kFALSE);

    // calculate ratio
    TH1* ratio = (TH1*) result[n]->Clone("ratio");
    ratio->Divide(mc[n], ratio, 1, 1, "B");
    ratio->Add(ratioBase, -1);

    ratio->SetTitle(";true multiplicity;Ratio_{syst} (t/u) - Ratio (t/u)");
    ratio->GetYaxis()->SetRangeUser(-1, 1);
    ratio->SetLineColor(n+1);
    ratio->DrawCopy((n == 0) ? "HIST" : "SAME HIST");

    // get average of ratio
    Float_t sum = 0;
    for (Int_t i=2; i<=ratioRange; ++i)
      sum += TMath::Abs(ratio->GetBinContent(i));
    sum /= ratioRange-1;

    printf("%d) Average (2..%d) of |ratio - ratioBase| is %f\n", n, ratioRange, sum);
  }

  TLine* line = new TLine(0, 0, displayRange, 0);
  line->SetLineWidth(2);
  line->Draw();

  line = new TLine(0, 0.1, displayRange, 0.1);
  line->SetLineWidth(2);
  line->SetLineStyle(2);
  line->Draw();
  line = new TLine(0, -0.1, displayRange, -0.1);
  line->SetLineWidth(2);
  line->SetLineStyle(2);
  line->Draw();

  canvas->Modified();

  canvas->SaveAs(epsName);
  canvas->SaveAs(Form("%s.gif", epsName.Data()));

  return canvas;
}

TCanvas* DrawRatioDeductSmooth(TH1* mcBase, TH1* resultBase, Int_t nResultSyst, TH1** mc, TH1** result, TString epsName)
{
  // draws the ratios of each mc to the corresponding result
  // deducts from each ratio the ratio of mcBase / resultBase
  // smoothens the ratios by a sliding window

  // normalize
  resultBase->Scale(1.0 / resultBase->Integral(2, 200));
  mcBase->Scale(1.0 / mcBase->Integral(2, 200));

  // calculate ratio
  TH1* ratioBase = (TH1*) resultBase->Clone("ratioBase");
  ratioBase->Divide(mcBase, ratioBase, 1, 1, "B");

  TCanvas* canvas = new TCanvas(epsName, epsName, 800, 400);
  canvas->SetRightMargin(0.05);
  canvas->SetTopMargin(0.05);

  for (Int_t n=0; n<nResultSyst; ++n)
  {
    // normalize
    result[n]->Scale(1.0 / result[n]->Integral(2, 200));
    mc[n]->Scale(1.0 / mc[n]->Integral(2, 200));

    result[n]->GetXaxis()->SetRangeUser(0, displayRange);
    result[n]->SetStats(kFALSE);

    // calculate ratio
    TH1* ratio = (TH1*) result[n]->Clone("ratio");
    ratio->Divide(mc[n], ratio, 1, 1, "B");
    ratio->Add(ratioBase, -1);

    //new TCanvas; ratio->DrawCopy();
    // clear 0 bin
    ratio->SetBinContent(1, 0); ratio->SetBinError(1, 0);

    Smooth(ratio);

    //ratio->SetLineColor(1); ratio->DrawCopy("SAME");

    canvas->cd();
    ratio->SetTitle(";true multiplicity;Ratio_{syst} (t/u) - Ratio (t/u) (smoothed)");
    ratio->GetYaxis()->SetRangeUser(-0.3, 0.3);
    ratio->SetLineColor((n / 2)+1);
    ratio->SetLineStyle((n % 2)+1);
    ratio->DrawCopy((n == 0) ? "HIST" : "SAME HIST");

    // get average of ratio
    Float_t sum = 0;
    for (Int_t i=2; i<=150; ++i)
      sum += TMath::Abs(ratio->GetBinContent(i));
    sum /= 149;

    printf("%d) Average (2..150) of |ratio - ratioBase| is %f\n", n, sum);
  }

  TLine* line = new TLine(0, 0, displayRange, 0);
  line->SetLineWidth(2);
  line->Draw();

  line = new TLine(0, 0.1, displayRange, 0.1);
  line->SetLineWidth(2);
  line->SetLineStyle(2);
  line->Draw();
  line = new TLine(0, -0.1, displayRange, -0.1);
  line->SetLineWidth(2);
  line->SetLineStyle(2);
  line->Draw();

  canvas->Modified();

  canvas->SaveAs(epsName);
  canvas->SaveAs(Form("%s.gif", epsName.Data()));

  return canvas;
}

void DrawResiduals(TH1* measured, TH1* unfoldedFolded, const char* epsName)
{
  // normalize
  unfoldedFolded->Scale(1.0 / unfoldedFolded->Integral(2, 200));
  unfoldedFolded->Scale(measured->Integral(2, 200));

  TCanvas* canvas = new TCanvas(epsName, epsName, 800, 600);
  canvas->Range(0, 0, 1, 1);

  TPad* pad1 = new TPad(Form("%s_pad1", epsName), "", 0, 0.5, 0.98, 0.98);
  pad1->Draw();
  pad1->SetGridx();
  pad1->SetGridy();

  TPad* pad2 = new TPad(Form("%s_pad2", epsName), "", 0, 0.02, 0.98, 0.5);
  pad2->Draw();
  pad2->SetGridx();
  pad2->SetGridy();

  TPad* pad3 = new TPad(Form("%s_pad3", epsName), "", 0.15, 0.5, 0.35, 0.75);
  pad3->SetGridx();
  pad3->SetGridy();
  pad3->SetRightMargin(0.05);
  pad3->SetTopMargin(0.05);
  pad3->Draw();

  pad1->SetRightMargin(0.05);
  pad2->SetRightMargin(0.05);

  // no border between them
  pad1->SetBottomMargin(0);
  pad2->SetTopMargin(0);

  pad1->cd();

  measured->GetXaxis()->SetLabelSize(0.06);
  measured->GetYaxis()->SetLabelSize(0.06);
  measured->GetXaxis()->SetTitleSize(0.06);
  measured->GetYaxis()->SetTitleSize(0.06);
  measured->GetYaxis()->SetTitleOffset(0.6);

  measured->GetXaxis()->SetRangeUser(0, 150);

  measured->SetTitle(";measured multiplicity;Entries");
  measured->SetStats(kFALSE);

  measured->DrawCopy("HIST");
  gPad->SetLogy();

  unfoldedFolded->SetMarkerStyle(5);
  unfoldedFolded->SetMarkerColor(2);
  unfoldedFolded->SetLineColor(0);
  unfoldedFolded->DrawCopy("SAME P");

  TLegend* legend = new TLegend(0.6, 0.65, 0.95, 0.9);
  legend->AddEntry(measured, "measured distribution");
  legend->AddEntry(unfoldedFolded, "R #otimes unfolded distribution");
  legend->SetFillColor(0);
  legend->Draw();

  pad2->cd();
  pad2->SetBottomMargin(0.15);

  // calculate ratio
  measured->Sumw2();
  TH1* residual = (TH1*) measured->Clone("residual");
  unfoldedFolded->Sumw2();

  residual->Add(unfoldedFolded, -1);

  // projection
  TH1* residualHist = new TH1F("residualHist", ";", 15, -3, 3);

  for (Int_t i=1; i<=residual->GetNbinsX(); ++i)
  {
    if (measured->GetBinError(i) > 0)
    {
      residual->SetBinContent(i, residual->GetBinContent(i) / measured->GetBinError(i));
      residual->SetBinError(i, 1);

      residualHist->Fill(residual->GetBinContent(i));
    }
    else
    {
      residual->SetBinContent(i, 0);
      residual->SetBinError(i, 0);
    }
  }

  residual->GetYaxis()->SetTitle("Residuals   1/e (M - R #otimes U)");
  residual->GetYaxis()->SetRangeUser(-4.5, 4.5);
  residual->DrawCopy();

  TLine* line = new TLine(-0.5, 0, 150.5, 0);
  line->SetLineWidth(2);
  line->Draw();

  pad3->cd();
  residualHist->SetStats(kFALSE);
  residualHist->GetXaxis()->SetLabelSize(0.08);
  residualHist->Fit("gaus");
  residualHist->Draw();

  canvas->Modified();
  canvas->SaveAs(canvas->GetName());

  //const char* epsName2 = "proj.eps";
  //TCanvas* canvas = new TCanvas(epsName2, epsName2, 800, 600);
  //canvas->SetGridx();
  //canvas->SetGridy();

  //canvas->SaveAs(canvas->GetName());
}

void bayesianExample()
{
  TStopwatch watch;
  watch.Start();

  gSystem->Load("libPWG0base");

  TFile::Open(correctionFile);
  AliMultiplicityCorrection* mult = new AliMultiplicityCorrection("Multiplicity", "Multiplicity");
  mult->LoadHistograms("Multiplicity");

  TFile::Open(measuredFile);
  AliMultiplicityCorrection* mult2 = new AliMultiplicityCorrection("Multiplicity2", "Multiplicity2");
  mult2->LoadHistograms("Multiplicity");

  mult->SetMultiplicityESD(etaRange, mult2->GetMultiplicityESD(etaRange));

  mult->ApplyBayesianMethod(etaRange, kFALSE, AliMultiplicityCorrection::kTrVtx, 1, 100);

  TH1* mcHist = mult2->GetMultiplicityVtx(etaRange)->ProjectionY("mymc");
  TH1* result = mult->GetMultiplicityESDCorrected(etaRange);

  //mult->DrawComparison("bayesianExample", etaRange, kFALSE, kTRUE, mcHist, kTRUE);
  DrawResultRatio(mcHist, result, "bayesianExample.eps");

  //Printf("KolmogorovTest says PROB = %f", mcHist->KolmogorovTest(result, "D"));
  //Printf("Chi2Test says PROB = %f", mcHist->Chi2Test(result));

  // draw residual plot

  // TODO take out efficiency correction if other than AliMultiplicityCorrection::kTrVtx
  TH2* convoluted = mult->CalculateMultiplicityESD(result, etaRange);
  TH1* convolutedProj = convoluted->ProjectionY("convolutedProj", -1, -1, "e");

  TH1* measured = mult2->GetMultiplicityESD(etaRange)->ProjectionY("measured");

  DrawResiduals(measured, convolutedProj, "bayesianResiduals.eps");

  watch.Stop();
  watch.Print();
}

void chi2FluctuationResult()
{
  gSystem->Load("libPWG0base");

  TFile::Open(correctionFile);
  AliMultiplicityCorrection* mult = new AliMultiplicityCorrection("Multiplicity", "Multiplicity");
  mult->LoadHistograms("Multiplicity");

  TFile::Open(measuredFile);
  AliMultiplicityCorrection* mult2 = new AliMultiplicityCorrection("Multiplicity2", "Multiplicity2");
  mult2->LoadHistograms("Multiplicity");

  mult->SetMultiplicityESD(etaRange, mult2->GetMultiplicityESD(etaRange));
  mult->SetRegularizationParameters(AliMultiplicityCorrection::kNone, 0);
  mult->ApplyMinuitFit(etaRange, kFALSE, AliMultiplicityCorrection::kTrVtx);

  TH1* mcHist = mult2->GetMultiplicityVtx(etaRange)->ProjectionY("mymc");
  //TH1* result = mult->GetMultiplicityESDCorrected(etaRange);

  mult->DrawComparison("MinuitChi2", etaRange, kFALSE, kTRUE, mcHist, kTRUE);

  TCanvas* canvas = (TCanvas*) gROOT->FindObject("MinuitChi2_DrawComparison_3");
  canvas->SaveAs("chi2FluctuationResult.eps");
}

void chi2Example()
{
  TStopwatch watch;
  watch.Start();

  gSystem->Load("libPWG0base");

  TFile::Open(correctionFile);
  AliMultiplicityCorrection* mult = new AliMultiplicityCorrection("Multiplicity", "Multiplicity");
  mult->LoadHistograms("Multiplicity");

  TFile::Open(measuredFile);
  AliMultiplicityCorrection* mult2 = new AliMultiplicityCorrection("Multiplicity2", "Multiplicity2");
  mult2->LoadHistograms("Multiplicity");

  mult->SetMultiplicityESD(etaRange, mult2->GetMultiplicityESD(etaRange));
  mult->SetRegularizationParameters(AliMultiplicityCorrection::kPol1, 10000);
  mult->ApplyMinuitFit(etaRange, kFALSE, AliMultiplicityCorrection::kTrVtx);

  TH1* mcHist = mult2->GetMultiplicityVtx(etaRange)->ProjectionY("mymc");
  TH1* result = mult->GetMultiplicityESDCorrected(etaRange);

  DrawResultRatio(mcHist, result, "chi2Example.eps");

  watch.Stop();
  watch.Print();
}

void chi2ExampleTPC()
{
  TStopwatch watch;
  watch.Start();

  gSystem->Load("libPWG0base");

  TFile::Open(correctionFileTPC);
  AliMultiplicityCorrection* mult = new AliMultiplicityCorrection("Multiplicity", "Multiplicity");
  mult->LoadHistograms("Multiplicity");

  TFile::Open(measuredFileTPC);
  AliMultiplicityCorrection* mult2 = new AliMultiplicityCorrection("Multiplicity2", "Multiplicity2");
  mult2->LoadHistograms("Multiplicity");

  mult->SetMultiplicityESD(etaRangeTPC, mult2->GetMultiplicityESD(etaRangeTPC));
  mult->SetRegularizationParameters(AliMultiplicityCorrection::kPol1, 10000);
  mult->ApplyMinuitFit(etaRangeTPC, kFALSE, AliMultiplicityCorrection::kTrVtx);

  TH1* mcHist = mult2->GetMultiplicityVtx(etaRangeTPC)->ProjectionY("mymc");
  TH1* result = mult->GetMultiplicityESDCorrected(etaRangeTPC);

  DrawResultRatio(mcHist, result, "chi2ExampleTPC.eps");

  watch.Stop();
  watch.Print();
}

void bayesianNBD()
{
  gSystem->Load("libPWG0base");
  TFile::Open("multiplicityMC_3M.root");

  AliMultiplicityCorrection* mult = new AliMultiplicityCorrection("Multiplicity", "Multiplicity");
  mult->LoadHistograms("Multiplicity");

  TFile::Open("multiplicityMC_3M_NBD.root");
  AliMultiplicityCorrection* mult2 = new AliMultiplicityCorrection("Multiplicity2", "Multiplicity2");
  mult2->LoadHistograms("Multiplicity");

  mult->SetMultiplicityESD(etaRange, mult2->GetMultiplicityESD(etaRange));
  mult->ApplyBayesianMethod(etaRange, kFALSE, AliMultiplicityCorrection::kTrVtx, 1.0, 100);

  TH1* mcHist = mult2->GetMultiplicityVtx(etaRange)->ProjectionY("mymc");

  mcHist->Sumw2();
  TH1* result = mult->GetMultiplicityESDCorrected(etaRange);

  //mult->DrawComparison("bayesianNBD", etaRange, kFALSE, kTRUE, mcHist);
  DrawResultRatio(mcHist, result, "bayesianNBD.eps");
}

void minimizationNBD()
{
  gSystem->Load("libPWG0base");
  TFile::Open("multiplicityMC_3M.root");

  AliMultiplicityCorrection* mult = new AliMultiplicityCorrection("Multiplicity", "Multiplicity");
  mult->LoadHistograms("Multiplicity");

  TFile::Open("multiplicityMC_3M_NBD.root");
  AliMultiplicityCorrection* mult2 = new AliMultiplicityCorrection("Multiplicity2", "Multiplicity2");
  mult2->LoadHistograms("Multiplicity");

  mult->SetMultiplicityESD(etaRange, mult2->GetMultiplicityESD(etaRange));
  mult->SetRegularizationParameters(AliMultiplicityCorrection::kPol1, 10000);
  mult->ApplyMinuitFit(etaRange, kFALSE, AliMultiplicityCorrection::kTrVtx);

  TH1* mcHist = mult2->GetMultiplicityVtx(etaRange)->ProjectionY("mymc");

  mcHist->Sumw2();
  TH1* result = mult->GetMultiplicityESDCorrected(etaRange);

  //mult->DrawComparison("minimizationNBD", etaRange, kFALSE, kTRUE, mcHist);
  DrawResultRatio(mcHist, result, "minimizationNBD.eps");
}

void minimizationInfluenceAlpha()
{
  gSystem->Load("libPWG0base");

  TFile::Open(measuredFile);
  AliMultiplicityCorrection* mult2 = new AliMultiplicityCorrection("Multiplicity2", "Multiplicity2");
  mult2->LoadHistograms("Multiplicity");

  TH1* mcHist = mult2->GetMultiplicityVtx(etaRange)->ProjectionY("mymc");
  mcHist->Scale(1.0 / mcHist->Integral());
  mcHist->GetXaxis()->SetRangeUser(0, 200);
  mcHist->SetStats(kFALSE);
  mcHist->SetTitle(";true multiplicity;P_{N}");

  TCanvas* canvas = new TCanvas("minimizationInfluenceAlpha", "minimizationInfluenceAlpha", 1000, 300);
  canvas->Divide(3, 1);

  TFile::Open("eval-2M-1M/EvaluateChi2MethodDetail.root");

  TH1* hist1 = (TH1*) gFile->Get("MinuitChi2_00_2_100.000000");
  TH1* hist2 = (TH1*) gFile->Get("MinuitChi2_03_2_100000.000000");
  TH1* hist3 = (TH1*) gFile->Get("MinuitChi2_06_2_100000000.000000");

  mcHist->Rebin(2);  mcHist->Scale(0.5);
  hist1->Rebin(2);   hist1->Scale(0.5);
  hist2->Rebin(2);   hist2->Scale(0.5);
  hist3->Rebin(2);   hist3->Scale(0.5);

  mcHist->GetXaxis()->SetRangeUser(0, 200);

  canvas->cd(1);
  gPad->SetLogy();
  mcHist->Draw();
  hist1->SetMarkerStyle(5);
  hist1->SetMarkerColor(2);
  hist1->Draw("SAME PE");

  canvas->cd(2);
  gPad->SetLogy();
  mcHist->Draw();
  hist2->SetMarkerStyle(5);
  hist2->SetMarkerColor(2);
  hist2->Draw("SAME PE");

  canvas->cd(3);
  gPad->SetLogy();
  mcHist->Draw();
  hist3->SetMarkerStyle(5);
  hist3->SetMarkerColor(2);
  hist3->Draw("SAME PE");

  canvas->SaveAs("minimizationInfluenceAlpha.eps");
}

void NBDFit()
{
  gSystem->Load("libPWG0base");

  TFile::Open(correctionFile);
  AliMultiplicityCorrection* mult = new AliMultiplicityCorrection("Multiplicity", "Multiplicity");
  mult->LoadHistograms("Multiplicity");

  TH1* fCurrentESD = mult->GetMultiplicityVtx(etaRange)->ProjectionY();
  fCurrentESD->Sumw2();
  fCurrentESD->Scale(1.0 / fCurrentESD->Integral());

  TF1* func = new TF1("nbd", "[0] * TMath::Binomial([2]+TMath::Nint(x)-1, [2]-1) * pow([1] / ([1]+[2]), TMath::Nint(x)) * pow(1 + [1]/[2], -[2])");
  func->SetParNames("scaling", "averagen", "k");
  func->SetParLimits(0, 0.001, fCurrentESD->GetMaximum() * 1000);
  func->SetParLimits(1, 0.001, 1000);
  func->SetParLimits(2, 0.001, 1000);
  func->SetParameters(fCurrentESD->GetMaximum() * 100, 10, 2);

  TF1* lognormal = new TF1("lognormal", "[0]*exp(-(log(x)-[1])^2/(2*[2]^2))/(x*[2]*TMath::Sqrt(2*TMath::Pi()))", 0.01, 500);
  lognormal->SetParNames("scaling", "mean", "sigma");
  lognormal->SetParameters(1, 1, 1);
  lognormal->SetParLimits(0, 0, 10);
  lognormal->SetParLimits(1, 0, 100);
  lognormal->SetParLimits(2, 1e-3, 10);

  TCanvas* canvas = new TCanvas("c1", "c1", 700, 400);
  fCurrentESD->SetStats(kFALSE);
  fCurrentESD->GetYaxis()->SetTitleOffset(1.3);
  fCurrentESD->SetTitle(";true multiplicity (N);P_{N}");
  fCurrentESD->Draw("HIST");
  fCurrentESD->GetXaxis()->SetRangeUser(0, 200);
  fCurrentESD->Fit(func, "W0", "", 0, 50);
  func->SetRange(0, 100);
  func->Draw("SAME");
  printf("chi2 = %f\n", func->GetChisquare());

  fCurrentESD->Fit(lognormal, "W0", "", 0.01, 100);
  lognormal->SetLineColor(2);
  lognormal->SetLineStyle(2);
  lognormal->SetRange(0, 100);
  lognormal->Draw("SAME");

  canvas->SaveAs("NBDFit.eps");
}

void DifferentSamples()
{
  // data generated by runMultiplicitySelector.C DifferentSamples

  const char* name = "DifferentSamples";

  TFile* file = TFile::Open(Form("%s.root", name));

  TCanvas* canvas = new TCanvas(name, name, 800, 600);
  canvas->Divide(2, 2);

  for (Int_t i=0; i<4; ++i)
  {
    canvas->cd(i+1);
    gPad->SetTopMargin(0.05);
    gPad->SetRightMargin(0.05);
    TH1* chi2Result = (TH1*) file->Get(Form("chi2Result_%d", i));
    TH1* bayesResult = (TH1*) file->Get(Form("bayesResult_%d", i));
    TH1* mc = (TH1*) file->Get(Form("mc_%d", i));
    mc->Sumw2();

    chi2Result->Divide(chi2Result, mc, 1, 1, "");
    bayesResult->Divide(bayesResult, mc, 1, 1, "");

    chi2Result->SetTitle(";true multiplicity;unfolded measured/MC");
    chi2Result->GetXaxis()->SetRangeUser(0, 150);
    chi2Result->GetYaxis()->SetRangeUser(0.5, 1.5);
    chi2Result->GetYaxis()->SetTitleOffset(1.2);
    chi2Result->SetLineColor(1);
    chi2Result->SetStats(kFALSE);

    bayesResult->SetStats(kFALSE);
    bayesResult->SetLineColor(2);

    chi2Result->DrawCopy("HIST");
    bayesResult->DrawCopy("SAME HIST");

    TLine* line = new TLine(0, 1, 150, 1);
    line->Draw();
  }

  canvas->SaveAs(Form("%s.eps", canvas->GetName()));
}

void StartingConditions()
{
  // data generated by runMultiplicitySelector.C StartingConditions

  const char* name = "StartingConditions";

  TFile* file = TFile::Open(Form("%s.root", name));

  TCanvas* canvas = new TCanvas(name, name, 800, 400);
  canvas->Divide(2, 1);

  TH1* mc = (TH1*) file->Get("mc");
  mc->Sumw2();
  mc->Scale(1.0 / mc->Integral());

  //Int_t marker[] = {24, 25, 26, 27, 28, 2, 3, 4, 5};

  TLegend* legend = new TLegend(0.6, 0.7, 0.95, 0.95);
  legend->SetFillColor(0);

  const char* names[] = { "True", "Measured 1", "Measured 2", "Measured 3", "NBD", "Flat" };

  for (Int_t i=0; i<6; ++i)
  {
    Int_t id = i;
    if (id > 2)
      id += 2;

    TH1* chi2Result = (TH1*) file->Get(Form("chi2Result_%d", id));
    TH1* bayesResult = (TH1*) file->Get(Form("bayesResult_%d", id));

    chi2Result->Divide(chi2Result, mc, 1, 1, "");
    bayesResult->Divide(bayesResult, mc, 1, 1, "");

    chi2Result->SetTitle("a) #chi^{2} minimization;true multiplicity;unfolded / MC");
    chi2Result->GetXaxis()->SetRangeUser(0, 150);
    chi2Result->GetYaxis()->SetRangeUser(0.8, 1.2);
    chi2Result->GetYaxis()->SetTitleOffset(1.5);
    //chi2Result->SetMarkerStyle(marker[i]);
    chi2Result->SetLineColor(i+1);
    chi2Result->SetMarkerColor(i+1);
    chi2Result->SetStats(kFALSE);

    bayesResult->SetTitle("b) Bayesian method;true multiplicity;unfolded / MC");
    bayesResult->GetXaxis()->SetRangeUser(0, 150);
    bayesResult->GetYaxis()->SetRangeUser(0.8, 1.2);
    bayesResult->GetYaxis()->SetTitleOffset(1.5);
    bayesResult->SetStats(kFALSE);
    //bayesResult->SetLineColor(2);
    bayesResult->SetLineColor(i+1);

    canvas->cd(1);
    gPad->SetLeftMargin(0.12);
    chi2Result->DrawCopy((i == 0) ? "HIST" : "HIST SAME");

    canvas->cd(2);
    gPad->SetLeftMargin(0.12);
    bayesResult->DrawCopy((i == 0) ? "HIST" : "HIST SAME");

    //TLine* line = new TLine(0, 1, 150, 1);
    //line->Draw();

    legend->AddEntry(chi2Result, names[i]);
  }

  canvas->cd(1);
  legend->Draw();
  canvas->SaveAs(Form("%s.eps", canvas->GetName()));
}

void StatisticsPlot()
{
  const char* name = "StatisticsPlot";

  TFile* file = TFile::Open(Form("%s.root", name));

  TCanvas* canvas = new TCanvas(name, name, 600, 400);

  TGraph* fitResultsChi2 = (TGraph*) file->Get("fitResultsChi2");
  fitResultsChi2->SetTitle(";number of measured events;P_{1}");
  fitResultsChi2->GetYaxis()->SetRangeUser(0, 2);
  fitResultsChi2->Draw("AP");

  TF1* f = new TF1("f", "[0]/x", 1, 1e4);
  fitResultsChi2->Fit(f);

  canvas->SaveAs(Form("%s.eps", canvas->GetName()));

  TH1* mc[5];
  TH1* result[5];

  const char* plotname = "chi2Result";

  name = "StatisticsPlotRatios";
  canvas = new TCanvas(name, name, 600, 400);

  for (Int_t i=0; i<5; ++i)
  {
    mc[i] = (TH1*) file->Get(Form("mc_%d", i));
    result[i] = (TH1*) file->Get(Form("%s_%d", plotname, i));

    result[i]->SetLineColor(i+1);
    result[i]->Draw(((i == 0) ? "" : "SAME"));
  }
}

void SystematicLowEfficiency()
{
  gSystem->Load("libPWG0base");

  TFile::Open(correctionFile);
  AliMultiplicityCorrection* mult = new AliMultiplicityCorrection("Multiplicity", "Multiplicity");
  mult->LoadHistograms("Multiplicity");

  // calculate result with systematic effect
  TFile::Open("multiplicityMC_100k_1_loweff.root");
  AliMultiplicityCorrection* mult2 = new AliMultiplicityCorrection("Multiplicity2", "Multiplicity2");
  mult2->LoadHistograms("Multiplicity");

  mult->SetMultiplicityESD(etaRange, mult2->GetMultiplicityESD(etaRange));
  mult->SetRegularizationParameters(AliMultiplicityCorrection::kPol1, 10000);
  mult->ApplyMinuitFit(etaRange, kFALSE, AliMultiplicityCorrection::kTrVtx);

  TH1* mcHist = mult2->GetMultiplicityVtx(etaRange)->ProjectionY("mymc");
  TH1* result1 = (TH1*) mult->GetMultiplicityESDCorrected(etaRange)->Clone("result1");

  DrawResultRatio(mcHist, result1, "SystematicLowEfficiencyLow.eps");

  // calculate normal result
  TFile::Open("multiplicityMC_100k_1.root");
  mult2->LoadHistograms("Multiplicity");

  mult->SetMultiplicityESD(etaRange, mult2->GetMultiplicityESD(etaRange));
  mult->ApplyMinuitFit(etaRange, kFALSE, AliMultiplicityCorrection::kTrVtx);

  TH1* result2 = (TH1*) mult->GetMultiplicityESDCorrected(etaRange)->Clone("result2");

  DrawResultRatio(mcHist, result2, "SystematicLowEfficiencyOK.eps");

  Draw2ResultRatio(mcHist, result1, result2, "SystematicLowEfficiency.eps");
}

void SystematicMisalignment()
{
  gSystem->Load("libPWG0base");

  TFile::Open(correctionFile);
  AliMultiplicityCorrection* mult = new AliMultiplicityCorrection("Multiplicity", "Multiplicity");
  mult->LoadHistograms("Multiplicity");

  TFile::Open("multiplicityMC_100k_fullmis.root");
  AliMultiplicityCorrection* mult2 = new AliMultiplicityCorrection("Multiplicity2", "Multiplicity2");
  mult2->LoadHistograms("Multiplicity");

  mult->SetMultiplicityESD(etaRange, mult2->GetMultiplicityESD(etaRange));
  mult->SetRegularizationParameters(AliMultiplicityCorrection::kPol1, 10000);
  mult->ApplyMinuitFit(etaRange, kFALSE, AliMultiplicityCorrection::kTrVtx);

  TH1* mcHist = mult2->GetMultiplicityVtx(etaRange)->ProjectionY("mymc");
  TH1* result = mult->GetMultiplicityESDCorrected(etaRange);

  DrawResultRatio(mcHist, result, "SystematicMisalignment.eps");
}

void SystematicMisalignmentTPC()
{
  gSystem->Load("libPWG0base");

  SetTPC();

  TFile::Open(correctionFile);
  AliMultiplicityCorrection* mult = new AliMultiplicityCorrection("Multiplicity", "Multiplicity");
  mult->LoadHistograms("Multiplicity");

  TFile::Open("multiplicityMC_TPC_100k_fullmis.root");
  AliMultiplicityCorrection* mult2 = new AliMultiplicityCorrection("Multiplicity2", "Multiplicity2");
  mult2->LoadHistograms("Multiplicity");

  mult->SetMultiplicityESD(etaRange, mult2->GetMultiplicityESD(etaRange));
  mult->SetRegularizationParameters(AliMultiplicityCorrection::kPol1, 10000);
  mult->ApplyMinuitFit(etaRange, kFALSE, AliMultiplicityCorrection::kTrVtx);

  TH1* mcHist = mult2->GetMultiplicityVtx(etaRange)->ProjectionY("mymc");
  TH1* result = mult->GetMultiplicityESDCorrected(etaRange);

  DrawResultRatio(mcHist, result, "SystematicMisalignmentTPC.eps");
}

void EfficiencySpecies()
{
  loadlibs();

  Int_t marker[] = {24, 25, 26};
  Int_t color[] = {1, 2, 4};

  // SPD TPC
  //const char* fileName[] = { "multiplicityMC_400k_syst.root", "multiplicityMC_TPC_4kfiles_syst.root" };
  const char* fileName[] = { "spd/multiplicity.root", "tpc/multiplicity.root" };
  Float_t etaRange[] = {0.49, 0.9};
  const char* titles[] = { "SPD Tracklets", "TPC Tracks" };

  TCanvas* canvas = new TCanvas("EfficiencySpecies", "EfficiencySpecies", 1000, 500);
  canvas->Divide(2, 1);

  for (Int_t loop=0; loop<2; ++loop)
  {
    Printf("%s", fileName[loop]);

    AliCorrection* correction[4];

    canvas->cd(loop+1);

    gPad->SetGridx();
    gPad->SetGridy();
    gPad->SetRightMargin(0.05);
    //gPad->SetTopMargin(0.05);

    TLegend* legend = new TLegend(0.7, 0.4, 0.85, 0.6);
    legend->SetFillColor(0);
    legend->SetEntrySeparation(0.2);

    Float_t below = 0;
    Float_t total = 0;

    TFile* file = TFile::Open(fileName[loop]);
    if (!file)
    {
      Printf("Could not open %s", fileName[loop]);
      return;
    }

    Float_t sumGen = 0;
    Float_t sumMeas = 0;

    for (Int_t i=0; i<3; ++i)
    {
      Printf("correction %d", i);

      TString name; name.Form("correction_%d", i);
      correction[i] = new AliCorrection(name, name);
      correction[i]->LoadHistograms();

      TH3* gene = correction[i]->GetTrackCorrection()->GetGeneratedHistogram();
      TH3* meas = correction[i]->GetTrackCorrection()->GetMeasuredHistogram();

      // limit vtx axis
      gene->GetXaxis()->SetRangeUser(-3.9, 3.9);
      meas->GetXaxis()->SetRangeUser(-3.9, 3.9);

      // empty over/underflow bin in eta, setting range to +-2 is not enough because this is the maximum range, Project3D takes them into account then (might be a bug)
      /*for (Int_t x = 1; x <= gene->GetNbinsX(); x++)
        for (Int_t z = 1; z <= gene->GetNbinsZ(); z++)
        {
          gene->SetBinContent(x, 0, z, 0);
          gene->SetBinContent(x, gene->GetNbinsY()+1, z, 0);
          meas->SetBinContent(x, 0, z, 0);
          meas->SetBinContent(x, gene->GetNbinsY()+1, z, 0);
        }*/

      // limit eta axis
      gene->GetYaxis()->SetRangeUser(-etaRange[loop], etaRange[loop]);
      meas->GetYaxis()->SetRangeUser(-etaRange[loop], etaRange[loop]);

      TH1* genePt = gene->Project3D(Form("z_%d", i));
      TH1* measPt = meas->Project3D(Form("z_%d", i));

      genePt->Sumw2();
      measPt->Sumw2();

      sumGen += genePt->Integral();
      sumMeas += measPt->Integral();

      TH1* effPt = (TH1*) genePt->Clone(Form("effPt_%d", i));
      effPt->Reset();
      effPt->Divide(measPt, genePt, 1, 1, "B");

      Int_t bin = 0;
      for (bin=20; bin>=1; bin--)
      {
        if (effPt->GetBinContent(bin) < 0.5)
          break;
      }

      Printf("Eff. below 50%% at bin %d, i.e. %.3f GeV/c", bin, effPt->GetXaxis()->GetBinUpEdge(bin));

      Float_t fraction = genePt->Integral(1, bin) / genePt->Integral();
      Printf("%.4f of the particles are below that momentum", fraction);

      below += genePt->Integral(1, bin);
      total += genePt->Integral();

      effPt->SetLineColor(color[i]);
      effPt->SetMarkerColor(color[i]);
      effPt->SetMarkerStyle(marker[i]);

      effPt->GetXaxis()->SetRangeUser(0.06, 1);
      effPt->GetYaxis()->SetRangeUser(0, 1);

      effPt->GetXaxis()->SetTitleOffset(1.1);
      effPt->GetYaxis()->SetTitleOffset(1.2);

      effPt->SetStats(kFALSE);
      effPt->SetTitle(titles[loop]);
      effPt->GetYaxis()->SetTitle("Efficiency");

      effPt->DrawCopy((i == 0) ? "" : "SAME");

      legend->AddEntry(effPt, ((i == 0) ? "#pi^{#pm}" : ((i == 1) ? "K^{#pm}" : "p,#bar{p}")));
    }

    Printf("In total %.4f of the particles are below their effective pt cut off", (Float_t) below / total);

    Printf("%f measured, %f generated, effiency: %f", sumGen, sumMeas, sumMeas / sumGen);

    legend->Draw();
  }

  canvas->SaveAs(Form("%s.eps", canvas->GetName()));
}

void ParticleSpeciesComparison1(Bool_t chi2 = kTRUE, const char* fileNameMC = "multiplicityMC_400k_syst_species.root", const char* fileNameESD = "multiplicityMC_100k_syst.root")
{
  gSystem->Load("libPWG0base");

  TFile::Open(fileNameESD);
  TH2F* hist = (TH2F*) gFile->Get(Form("Multiplicity/fMultiplicityESD%d", etaRange));
  TH2F* hist2 = (TH2F*) gFile->Get(Form("Multiplicity/fMultiplicityVtx%d", etaRange));

  TH1* results[10];

  // loop over cases (normal, enhanced/reduced ratios)
  Int_t nMax = 7;
  for (Int_t i = 0; i<nMax; ++i)
  {
    TString folder;
    folder.Form("Multiplicity_%d", i);

    AliMultiplicityCorrection* mult = new AliMultiplicityCorrection(folder, folder);

    TFile::Open(fileNameMC);
    mult->LoadHistograms();

    mult->SetMultiplicityESD(etaRange, hist);

    if (chi2)
    {
      mult->SetRegularizationParameters(AliMultiplicityCorrection::kPol1, 1e4);
      mult->ApplyMinuitFit(etaRange, kFALSE, AliMultiplicityCorrection::kTrVtx, kFALSE);
      //mult->DrawComparison(Form("ParticleSpeciesComparison_MinuitChi2_%d", i), etaRange, kFALSE, kTRUE, hist2->ProjectionY("mymchist"));
    }
    else
    {
      mult->ApplyBayesianMethod(etaRange, kFALSE, AliMultiplicityCorrection::kTrVtx, 1, 100);
      //mult->DrawComparison(Form("ParticleSpeciesComparison_Bayesian_%d", i), etaRange, kFALSE, kTRUE, hist2->ProjectionY("mymchist2"));
    }

    //Float_t averageRatio = 0;
    //mult->GetComparisonResults(0, 0, 0, &averageRatio);

    results[i] = (TH1*) mult->GetMultiplicityESDCorrected(etaRange)->Clone(Form("result_%d", i));

    //Printf("Case %d. Average ratio is %f", i, averageRatio);
  }

  DrawResultRatio(hist2->ProjectionY("mymchist", -1, -1, "e"), results[0], "ParticleSpeciesComparison1_1.eps");

  TH1* mc = hist2->ProjectionY("mymchist2", -1, -1, "e");

  for (Int_t i=1; i<=results[0]->GetNbinsX(); i++)
  {
    results[0]->SetBinError(i, 0);
    mc->SetBinError(i, 0);
  }

  const char* legendStrings[] = { "#pi^{#pm}", 0, "K^{#pm}", 0, "p,#bar{p}", 0 };

  DrawRatio(results[0], nMax-1, results+1, "ParticleSpeciesComparison1_2.eps", kFALSE, legendStrings);

  //not valid: draw chi2 uncertainty on top!
  /*TFile::Open("bayesianUncertainty_400k_100k_syst.root");
  TH1* errorHist = (TH1*) gFile->Get("errorBoth");
  errorHist->SetLineColor(1);
  errorHist->SetLineWidth(2);
  TH1* errorHist2 = (TH1*) errorHist->Clone("errorHist2");
  for (Int_t i=1; i<=errorHist->GetNbinsX(); i++)
  {
    errorHist->SetBinContent(i, errorHist->GetBinContent(i) + 1);
    errorHist2->SetBinContent(i, 1 - errorHist2->GetBinContent(i));
  }

  errorHist->DrawCopy("SAME");
  errorHist2->DrawCopy("SAME");*/

  //canvas->SaveAs(canvas->GetName());

  DrawRatio(mc, nMax, results, "ParticleSpeciesComparison1_3.eps", kTRUE, 0);

  //errorHist->DrawCopy("SAME");
  //errorHist2->DrawCopy("SAME");

  //canvas2->SaveAs(canvas2->GetName());
}

/*void ParticleSpeciesComparison2()
{
  gSystem->Load("libPWG0base");

  const char* fileNameMC = "multiplicityMC_400k_syst.root";
  const char* fileNameESD = "out.root"; // based on multiplicityMC_100k_syst.root
  Bool_t chi2 = 0;

  TFile::Open(fileNameMC);
  AliMultiplicityCorrection* mult = new AliMultiplicityCorrection("Multiplicity", "Multiplicity");
  mult->LoadHistograms();

  TH1* mc[10];
  TH1* results[10];

  // loop over cases (normal, enhanced/reduced ratios)
  Int_t nMax = 7;
  for (Int_t i = 0; i<nMax; ++i)
  {
    TString folder;
    folder.Form("Multiplicity_%d", i);

    AliMultiplicityCorrection* mult2 = new AliMultiplicityCorrection(folder, folder);

    TFile::Open(fileNameESD);
    mult2->LoadHistograms();

    mult->SetMultiplicityESD(etaRange, mult2->GetMultiplicityESD(etaRange));

    if (chi2)
    {
      mult->SetRegularizationParameters(AliMultiplicityCorrection::kPol1, 1e4);
      mult->ApplyMinuitFit(etaRange, kFALSE, AliMultiplicityCorrection::kTrVtx, kFALSE);
      //mult->DrawComparison(Form("ParticleSpeciesComparison_MinuitChi2_%d", i), etaRange, kFALSE, kTRUE, hist2->ProjectionY("mymchist"));
    }
    else
    {
      mult->ApplyBayesianMethod(etaRange, kFALSE, AliMultiplicityCorrection::kTrVtx, 1, 100);
      //mult->DrawComparison(Form("ParticleSpeciesComparison_Bayesian_%d", i), etaRange, kFALSE, kTRUE, hist2->ProjectionY("mymchist2"));
    }

    //Float_t averageRatio = 0;
    //mult->GetComparisonResults(0, 0, 0, &averageRatio);

    results[i] = (TH1*) mult->GetMultiplicityESDCorrected(etaRange)->Clone(Form("result_%d", i));

    TH2F* hist2 = mult2->GetMultiplicityVtx(etaRange);
    mc[i] = (TH1*) hist2->ProjectionY(Form("mymchist_%d", i), -1, -1, "e");

    //TString fileName; fileName.Form("ParticleSpeciesComparison2_%d.eps", i);
    //DrawResultRatio(hist2->ProjectionY("mymchist", -1, -1, "e"), results[i], fileName);

    //Printf("Case %d. Average ratio is %f", i, averageRatio);
  }

  DrawRatio(nMax, mc, results, "ParticleSpeciesComparison2.eps");
}*/

TH1* Invert(TH1* eff)
{
  // calculate corr = 1 / eff

  TH1* corr = (TH1*) eff->Clone(Form("%s_invert", eff->GetName()));
  corr->Reset();

  for (Int_t i=1; i<=eff->GetNbinsX(); i++)
  {
    if (eff->GetBinContent(i) > 0)
    {
      corr->SetBinContent(i, 1.0 / eff->GetBinContent(i));
      corr->SetBinError(i, eff->GetBinError(i) / eff->GetBinContent(i) * corr->GetBinContent(i));
    }
  }

  return corr;
}

void TriggerVertexCorrection()
{
  //
  // plots the correction performed on the unfolded spectrum to gain the spectrum for the full inelastic sample
  //

  gSystem->Load("libPWG0base");

  TFile::Open(correctionFile);
  AliMultiplicityCorrection* mult = new AliMultiplicityCorrection("Multiplicity", "Multiplicity");
  mult->LoadHistograms("Multiplicity");

  TH1* corrINEL = Invert(mult->GetEfficiency(etaRange, AliMultiplicityCorrection::kINEL));
  TH1* corrMB   = Invert(mult->GetEfficiency(etaRange, AliMultiplicityCorrection::kMB));

  TCanvas* canvas = new TCanvas("TriggerVertexCorrection", "TriggerVertexCorrection", 800, 600);

  corrINEL->SetStats(kFALSE);
  corrINEL->GetXaxis()->SetRangeUser(0, 20);
  corrINEL->GetYaxis()->SetRangeUser(0.5, 2.5);
  corrINEL->SetTitle(";true multiplicity;correction factor");
  corrINEL->SetMarkerStyle(22);
  corrINEL->Draw("PE");

  corrMB->SetStats(kFALSE);
  corrMB->SetLineColor(2);
  corrMB->SetMarkerStyle(25);
  corrMB->SetMarkerColor(2);
  corrMB->Draw("SAME PE");

  TLegend* legend = new TLegend(0.3, 0.5, 0.85, 0.65);
  legend->SetFillColor(0);
  legend->AddEntry(corrINEL, "correction to inelastic sample");
  legend->AddEntry(corrMB, "correction to minimum bias sample");

  legend->Draw();

  canvas->SaveAs(Form("%s.eps", canvas->GetName()));
}

void StatisticalUncertainty(Int_t methodType, Bool_t mc = kFALSE)
{
  gSystem->Load("libPWG0base");

  TFile::Open(correctionFile);
  AliMultiplicityCorrection* mult = new AliMultiplicityCorrection("Multiplicity", "Multiplicity");
  mult->LoadHistograms("Multiplicity");

  TFile::Open(measuredFile);
  AliMultiplicityCorrection* mult2 = new AliMultiplicityCorrection("Multiplicity2", "Multiplicity2");
  mult2->LoadHistograms("Multiplicity");

  mult->SetMultiplicityESD(etaRange, mult2->GetMultiplicityESD(etaRange));

  TH1* mcHist = mult2->GetMultiplicityVtx(etaRange)->ProjectionY("mymc");

  TH1* errorResponse = (TH1*) mult->StatisticalUncertainty((AliMultiplicityCorrection::MethodType) methodType, etaRange, kFALSE, AliMultiplicityCorrection::kTrVtx, kFALSE, kTRUE, ((mc) ? mcHist : 0))->Clone("errorResponse");

  TH1* errorMeasured = (TH1*) mult->StatisticalUncertainty((AliMultiplicityCorrection::MethodType) methodType, etaRange, kFALSE, AliMultiplicityCorrection::kTrVtx, kTRUE, kFALSE, ((mc) ? mcHist : 0))->Clone("errorMeasured");
  TH1* errorBoth = (TH1*) mult->StatisticalUncertainty((AliMultiplicityCorrection::MethodType) methodType, etaRange, kFALSE, AliMultiplicityCorrection::kTrVtx, kTRUE, kTRUE, ((mc) ? mcHist : 0))->Clone("errorBoth");

  if (!mc)
  {
    TH1* result = mult->GetMultiplicityESDCorrected(etaRange);
    DrawResultRatio(mcHist, result, "StatisticalUncertainty2.eps");
  }

  TCanvas* canvas = new TCanvas("StatisticalUncertainty", "StatisticalUncertainty", 600, 400);
  canvas->SetGridx();
  canvas->SetGridy();
  canvas->SetRightMargin(0.05);
  canvas->SetTopMargin(0.05);

  errorResponse->SetLineColor(1);
  errorResponse->GetXaxis()->SetRangeUser(0, 200);
  errorResponse->GetYaxis()->SetRangeUser(0, 0.3);
  errorResponse->SetStats(kFALSE);
  errorResponse->SetTitle(";true multiplicity;Uncertainty");

  errorResponse->Draw();

  errorMeasured->SetLineColor(2);
  errorMeasured->Draw("SAME");

  errorBoth->SetLineColor(4);
  errorBoth->Draw("SAME");

  Printf("Average errorResponse: %f", errorResponse->Integral(2, 150) / 149);
  Printf("Average errorMeasured: %f", errorMeasured->Integral(2, 150) / 149);
  Printf("Average errorBoth: %f", errorBoth->Integral(2, 150) / 149);

  canvas->SaveAs(Form("%s.eps", canvas->GetName()));

  TFile* file = new TFile(Form("%s.root", canvas->GetName()), "RECREATE");
  errorResponse->Write();
  errorMeasured->Write();
  errorBoth->Write();
  file->Close();
}

void StatisticalUncertaintyCompare(const char* det = "SPD")
{
  TFile* file1 = TFile::Open(Form("StatisticalUncertainty%sBayesian.root", det));
  TH1* errorResponse = (TH1*) file1->Get("errorResponse");
  TH1* errorMeasured = (TH1*) file1->Get("errorMeasured");
  TH1* errorBoth = (TH1*) file1->Get("errorBoth");

  TString str;
  str.Form("StatisticalUncertaintyCompare%s", det);

  TCanvas* canvas = new TCanvas(str, str, 600, 400);
  canvas->SetGridx();
  canvas->SetGridy();
  canvas->SetRightMargin(0.05);
  canvas->SetTopMargin(0.05);

  errorResponse->SetLineColor(1);
  errorResponse->GetXaxis()->SetRangeUser(1, (strcmp(det, "TPC") ? 200 : 100));
  errorResponse->GetYaxis()->SetRangeUser(0, 0.3);
  errorResponse->SetStats(kFALSE);
  errorResponse->GetYaxis()->SetTitleOffset(1.2);
  errorResponse->SetTitle(";true multiplicity;#sigma(U-T)/T");

  errorResponse->Draw();

  errorMeasured->SetLineColor(2);
  errorMeasured->Draw("SAME");

  errorBoth->SetLineColor(4);
  errorBoth->Draw("SAME");

  TFile* file2 = TFile::Open(Form("StatisticalUncertainty%sChi2.root", det));
  TH1* errorBoth2 = (TH1*) file2->Get("errorBoth");

  errorBoth2->SetLineColor(4);
  errorBoth2->SetLineStyle(2);
  errorBoth2->Draw("SAME");

  TLegend* legend = new TLegend(0.2, 0.6, 0.6, 0.9);
  legend->SetFillColor(0);
  legend->AddEntry(errorResponse, "response matrix (Bayesian)");
  legend->AddEntry(errorMeasured, "measured (Bayesian)");
  legend->AddEntry(errorBoth, "both (Bayesian)");
  legend->AddEntry(errorBoth2, "both (#chi^{2} minimization)");
  legend->Draw();

  canvas->SaveAs(Form("%s.eps", canvas->GetName()));
}

void EfficiencyComparison(Int_t eventType = 2, Bool_t uncertainty = kTRUE)
{
 const char* files[] = { "multiplicityMC_400k_syst_nd.root", "multiplicityMC_400k_syst_sd.root", "multiplicityMC_400k_syst_dd.root", "multiplicityMC_400k_syst_xsection.root" };

  gSystem->Load("libPWG0base");

  TCanvas* canvas = new TCanvas("EfficiencyComparison", "EfficiencyComparison", 800, 500);
  canvas->SetGridx();
  canvas->SetGridy();
  canvas->SetRightMargin(0.05);
  canvas->SetTopMargin(0.05);

  AliMultiplicityCorrection* data[4];
  TH1* effArray[4];

  Int_t markers[] = { 24, 25, 26, 5 };
  Int_t colors[] = { 1, 2, 3, 4 };

  TLegend* legend = new TLegend(0.45, 0.45, 0.9, 0.7);
  legend->SetFillColor(0);

  TH1* effError = 0;

  for (Int_t i=0; i<4; ++i)
  {
    TString name;
    name.Form("Multiplicity_%d", i);

    TFile::Open(files[i]);
    data[i] = new AliMultiplicityCorrection(name, name);

    if (i < 3)
    {
      data[i]->LoadHistograms("Multiplicity");
    }
    else
      data[i]->LoadHistograms("Multiplicity_0");

    TH1* eff = (TH1*) data[i]->GetEfficiency(etaRange, (AliMultiplicityCorrection::EventType) eventType)->Clone(Form("eff_%d", i));
    effArray[i] = eff;

    eff->GetXaxis()->SetRangeUser(0, 15);
    eff->GetYaxis()->SetRangeUser(0, 1.1);
    eff->SetStats(kFALSE);
    eff->SetTitle(";true multiplicity;Efficiency");
    eff->SetLineColor(colors[i]);
    eff->SetMarkerColor(colors[i]);
    eff->SetMarkerStyle(markers[i]);

    if (i == 3)
    {
      for (Int_t bin=1; bin<=eff->GetNbinsX(); bin++)
        eff->SetBinError(bin, 0);

      // loop over cross section combinations
      for (Int_t j=1; j<7; ++j)
      {
        AliMultiplicityCorrection* mult = new AliMultiplicityCorrection("Multtmp", "Multtmp");
        mult->LoadHistograms(Form("Multiplicity_%d", j));

        TH1* eff2 = mult->GetEfficiency(etaRange, (AliMultiplicityCorrection::EventType) eventType);

        for (Int_t bin=1; bin<=eff->GetNbinsX(); bin++)
        {
          // TODO we could also do asymmetric errors here
          Float_t deviation = TMath::Abs(eff->GetBinContent(bin) - eff2->GetBinContent(bin));

          eff->SetBinError(bin, TMath::Max(eff->GetBinError(bin), (Double_t) deviation));
        }
      }

      for (Int_t bin=1; bin<=20; bin++)
        if (eff->GetBinContent(bin) > 0)
          Printf("Bin %d: Error: %.2f", bin, 100.0 * eff->GetBinError(bin) / eff->GetBinContent(bin));
      
      if (uncertainty) {
	effError = (TH1*) eff->Clone("effError");
	effError->Reset();

	for (Int_t bin=2; bin<=eff->GetNbinsX(); bin++)
	  if (eff->GetBinContent(bin) > 0)
	    effError->SetBinContent(bin, 10.0 * eff->GetBinError(bin) / eff->GetBinContent(bin));

	effError->SetLineColor(1);
	effError->SetMarkerStyle(1);
	effError->DrawCopy("SAME HIST");
      }
    }

    eff->SetBinContent(1, 0);
    eff->SetBinError(1, 0);

    canvas->cd();
    if (i == 0)
    {
      eff->DrawCopy("P");
    }
    else
      eff->DrawCopy("SAME P");

    legend->AddEntry(eff, (((i == 0) ? "non diffractive" : ((i == 1) ? "single diffractive" : ((i == 2) ? "double diffractive" : "Pythia combined")))));
  }

  if (uncertainty)
    legend->AddEntry(effError, "relative syst. uncertainty #times 10");

  legend->Draw();

  canvas->SaveAs(Form("%s.eps", canvas->GetName()));
}

void ModelDependencyPlot()
{
  gSystem->Load("libPWG0base");

  TFile::Open("multiplicityMC_3M.root");
  AliMultiplicityCorrection* mult = new AliMultiplicityCorrection("Multiplicity", "Multiplicity");
  mult->LoadHistograms("Multiplicity");

  TH2* proj = (TH2*) mult->GetCorrelation(3)->Project3D("zy");

  TCanvas* canvas = new TCanvas("ModelDependencyPlot", "ModelDependencyPlot", 800, 400);
  canvas->SetGridx();
  canvas->SetGridy();
  //canvas->SetRightMargin(0.05);
  //canvas->SetTopMargin(0.05);

  canvas->Divide(2, 1);

  canvas->cd(2);
  gPad->SetLogy();
 
  Int_t selectedMult = 30;
  Int_t yMax = 200000;

  TH1* full = proj->ProjectionX("full");
  TH1* selected = proj->ProjectionY("selected", proj->GetXaxis()->FindBin(selectedMult), proj->GetXaxis()->FindBin(selectedMult)); 

  full->SetStats(kFALSE);
  full->GetXaxis()->SetRangeUser(0, 200);
  full->GetYaxis()->SetRangeUser(5, yMax);
  full->SetTitle(";multiplicity");

  selected->SetLineColor(0);
  selected->SetMarkerColor(2);
  selected->SetMarkerStyle(7);

  full->Draw();
  selected->Draw("SAME P");

  TLegend* legend = new TLegend(0.5, 0.65, 0.85, 0.85);
  legend->SetFillColor(0);
  legend->AddEntry(full, "true");
  legend->AddEntry(selected, "measured");
  legend->Draw();
 
  TLine* line = new TLine(selectedMult, 5, selectedMult, yMax);
  line->SetLineWidth(2);
  line->Draw();

  canvas->cd(1);
  gPad->SetLogy();

  full = proj->ProjectionY("full2");
  selected = proj->ProjectionX("selected2", proj->GetYaxis()->FindBin(selectedMult), proj->GetYaxis()->FindBin(selectedMult));

  full->SetStats(kFALSE);
  full->GetXaxis()->SetRangeUser(0, 200);
  full->GetYaxis()->SetRangeUser(5, yMax);
  full->SetTitle(";multiplicity");

  full->SetLineColor(0);
  full->SetMarkerColor(2);
  full->SetMarkerStyle(7);

  full->Draw("P");
  selected->Draw("SAME");

  legend->Draw();

  line = new TLine(selectedMult, 5, selectedMult, yMax);
  line->SetLineWidth(2);
  line->Draw();

  canvas->SaveAs(Form("%s.eps", canvas->GetName()));
}

void SystematicpTSpectrum()
{
  gSystem->Load("libPWG0base");

  TFile::Open("multiplicityMC_400k_syst_ptspectrum.root");
  AliMultiplicityCorrection* mult = new AliMultiplicityCorrection("Multiplicity", "Multiplicity");
  mult->LoadHistograms("Multiplicity");

  TFile::Open("multiplicityMC_100k_syst.root");
  AliMultiplicityCorrection* mult2 = new AliMultiplicityCorrection("Multiplicity2", "Multiplicity2");
  mult2->LoadHistograms("Multiplicity");

  mult->SetMultiplicityESD(etaRange, mult2->GetMultiplicityESD(etaRange));
  mult->SetRegularizationParameters(AliMultiplicityCorrection::kPol1, 10000);
  mult->ApplyMinuitFit(etaRange, kFALSE, AliMultiplicityCorrection::kTrVtx);

  TH1* mcHist = mult2->GetMultiplicityVtx(etaRange)->ProjectionY("mymc");
  TH1* result = mult->GetMultiplicityESDCorrected(etaRange);

  DrawResultRatio(mcHist, result, "SystematicpTSpectrum.eps");
}

// to be deleted
/*void covMatrix(Bool_t mc = kTRUE)
{
  gSystem->Load("libPWG0base");

  TFile::Open(correctionFile);
  AliMultiplicityCorrection* mult = new AliMultiplicityCorrection("Multiplicity", "Multiplicity");
  mult->LoadHistograms("Multiplicity");

  TFile::Open(measuredFile);
  AliMultiplicityCorrection* mult2 = new AliMultiplicityCorrection("Multiplicity2", "Multiplicity2");
  mult2->LoadHistograms("Multiplicity");

  mult->SetMultiplicityESD(etaRange, mult2->GetMultiplicityESD(etaRange));

  TH1* mcHist = mult2->GetMultiplicityVtx(etaRange)->ProjectionY("mymc");

  mult->BayesianStatisticsEffect(etaRange, kFALSE, AliMultiplicityCorrection::kTrVtx, kTRUE, kTRUE, 1, 100, ((mc) ? mcHist : 0));
}*/

Double_t FitPtFunc(Double_t *x, Double_t *par)
{
  Double_t xx = x[0];

  Float_t val1 = par[1] + par[2] * xx + par[3] * xx * xx;
  Float_t val2 = TMath::Exp(par[4] + par[5] * xx);

  const Float_t kTransitionWidth = 0;

  // power law part
  if (xx < par[0] - kTransitionWidth)
  {
    return val1;
  }
  /*else if (xx < par[0] + kTransitionWidth)
  {
    // smooth transition
    Float_t factor = (xx - par[0] + kTransitionWidth) / kTransitionWidth / 2;
    return (1 - factor) * val1 + factor * val2;
  }*/
  else
  {
    return val2;
  }
}

void FitPtNew(const char* fileName = "TruePt14TeV.root")
{
  gSystem->Load("libANALYSIS");
  gSystem->Load("libPWG0base");

  TFile::Open(fileName);

  TH1* genePt = (TH1*) gFile->Get("fHistPt");
  genePt->Sumw2();

  // normalize by bin width
  for (Int_t x=1; x<genePt->GetNbinsX(); x++)
    genePt->SetBinContent(x, genePt->GetBinContent(x) / genePt->GetBinWidth(x));

  genePt->GetXaxis()->SetRangeUser(0.05, 2.0);

  genePt->Scale(1.0 / genePt->Integral());

  TF1* func = new TF1("func", "[0]*TMath::Exp([1]*x*x)", 0.001, 100);
  //func->SetLineColor(2);
  func->SetParameters(1, -1);

  genePt->SetMarkerStyle(25);
  genePt->SetTitle("");
  genePt->SetStats(kFALSE);
  genePt->GetYaxis()->SetRangeUser(1e-4, genePt->GetMaximum() * 1.2);
  //func->Draw("SAME");

  genePt->Fit(func, "0", "", 0.05, 1);

  new TCanvas;
  genePt->DrawCopy("P");
  func->SetRange(0.02, 8);
  func->DrawCopy("SAME");
  gPad->SetLogy();
}

void FitPt(const char* fileName = "firstplots100k_truept.root")
{
  gSystem->Load("libPWG0base");

  TFile::Open(fileName);

  /*
  // merge corrections
  AliCorrection* correction[4];
  TList list;

  for (Int_t i=0; i<4; ++i)
  {
    Printf("correction %d", i);

    TString name; name.Form("correction_%d", i);
    correction[i] = new AliCorrection(name, name);
    correction[i]->LoadHistograms();

    if (i > 0)
      list.Add(correction[i]);
  }

  correction[0]->Merge(&list);

  TH3* gene = correction[0]->GetTrackCorrection()->GetGeneratedHistogram();

  // limit vtx, eta axis
  gene->GetXaxis()->SetRangeUser(-5.9, 5.9);
  gene->GetYaxis()->SetRangeUser(-1.99, 0.99);

  TH1* genePt = gene->Project3D("z");*/
  TH1* genePt = (TH1*) gFile->Get("fdNdpTTrue");
  if (!genePt)
    genePt = (TH1*) gFile->Get("fHistPt");
 
  genePt->Sumw2();

  //genePt->Scale(1.0 / genePt->Integral());

  // normalize by bin width
  for (Int_t x=1; x<genePt->GetNbinsX(); x++)
    genePt->SetBinContent(x, genePt->GetBinContent(x) / genePt->GetBinWidth(x));

  /// genePt->GetXaxis()->GetBinCenter(x));

  genePt->GetXaxis()->SetRangeUser(0, 7.9);
  //genePt->GetYaxis()->SetTitle("a.u.");

  //TF1* func = new TF1("func", "[0]*TMath::Exp([1]*x*x)", 0.001, 100);
  TF1* func = new TF1("func", "[0]*TMath::Exp([1]*x)+[2]/(1+(x*[4])**[3])", 0.001, 100);
  //func->SetLineColor(2);
  func->SetParameters(1, -1, 1, 1, 1);
  func->SetParLimits(3, 1, 10);
  func->SetParLimits(4, 0, 10);

  //TF1* func = new TF1("func", "[1]*x**[0]", 0.001, 100);

  //TF1* func = new TF1("func", &FitPtFunc, 0, 2, 6);
  //func->SetParameters(0.3, -2.34909e-01, 1.54394e+01, -3.04134e+01, 1.41912e+00, -2.79284e+00);
  //func->FixParameter(0, 0.314);
  //func->SetParLimits(0, 0.1, 0.3);

  genePt->SetMarkerStyle(25);
  genePt->SetTitle("");
  genePt->SetStats(kFALSE);
  genePt->GetYaxis()->SetRangeUser(1e-4, genePt->GetMaximum() * 1.2);
  //func->Draw("SAME");

  // fit only exp. part
  func->SetParameters(1, -1);
  func->FixParameter(2, 0);
  func->FixParameter(3, 1);
  func->FixParameter(4, 1);
  genePt->Fit(func, "0", "", 0.2, 1);

  new TCanvas;
  genePt->DrawCopy("P");
  func->SetRange(0.02, 8);
  func->DrawCopy("SAME");
  gPad->SetLogy();

  // now fix exp. parameters and fit second part
  Double_t param0 = func->GetParameter(0);
  Double_t param1 = func->GetParameter(1);
  func->SetParameters(0, -1, 1, 1, 1);
  func->FixParameter(0, 0);
  func->FixParameter(1, -1);
  func->ReleaseParameter(2);
  func->ReleaseParameter(3);
  func->ReleaseParameter(4);
  func->SetParLimits(3, 1, 10);
  func->SetParLimits(4, 0, 10);

  genePt->Fit(func, "0", "", 1.5, 4);

  new TCanvas;
  genePt->DrawCopy("P");
  func->SetRange(0.02, 8);
  func->DrawCopy("SAME");
  gPad->SetLogy();

  // fit both
  func->SetParameter(0, param0);
  func->SetParameter(1, param1);
  func->ReleaseParameter(0);
  func->ReleaseParameter(1);

  new TCanvas;
  genePt->DrawCopy("P");
  func->SetRange(0.02, 5);
  func->DrawCopy("SAME");
  gPad->SetLogy();

  genePt->Fit(func, "0", "", 0.2, 4);

  TCanvas* canvas = new TCanvas("FitPt", "FitPt", 800, 400);
  canvas->Divide(2, 1);
  canvas->cd(1);

  gPad->SetGridx();
  gPad->SetGridy();
  gPad->SetLeftMargin(0.13);
  gPad->SetRightMargin(0.05);
  gPad->SetTopMargin(0.05);

  genePt->GetXaxis()->SetRangeUser(0, 4.9);
  genePt->GetYaxis()->SetRangeUser(1e-2, 1e4);
  genePt->GetYaxis()->SetTitleOffset(1.4);
  genePt->GetXaxis()->SetTitleOffset(1.1);
  genePt->DrawCopy("P");
  func->SetRange(0.02, 5);
  func->DrawCopy("SAME");
  gPad->SetLogy();

  canvas->cd(2);

  TH1* genePtClone = (TH1*) genePt->Clone("genePtClone");
  genePtClone->Reset();
  genePtClone->DrawCopy("P");

  gPad->SetGridx();
  gPad->SetGridy();
  gPad->SetLeftMargin(0.13);
  gPad->SetRightMargin(0.05);
  gPad->SetTopMargin(0.05);

  func->DrawCopy("SAME");
  gPad->SetLogy();

  canvas->SaveAs(Form("%s.eps", canvas->GetName()));

  TH1* first = (TH1*) func->GetHistogram()->Clone("first");

  TCanvas* canvas2 = new TCanvas("FitPt2", "FitPt2", 600, 400);

  TFile* file = TFile::Open("ptspectrum_fit.root", "RECREATE");

  for (Int_t param=0; param<5; param++)
  {
    for (Int_t sign=0; sign<2; sign++)
    {
      TF1* func2 = (TF1*) func->Clone(Form("func_%d_%d", param, sign));  // new TF1(Form("func_%d_%d", param, sign), &FitPtFunc, 0, 2, 6);
      func2->SetParameters(func->GetParameters());
      //TF1* func2 = (TF1*) func->Clone(); // SetParameter after this does not work

      Float_t factor = ((sign == 0) ? 0.9 : 1.1);
      func2->SetParameter(param, func2->GetParameter(param) * factor);
      //func2->Print();

      canvas->cd(2);
      func2->SetLineWidth(1);
      func2->SetLineColor(2);
      func2->DrawCopy("SAME");

      canvas2->cd();
      TH1* second = func2->GetHistogram();
      second->Divide(first);
      second->SetLineColor(param + 1);
      second->GetYaxis()->SetRangeUser(0, 2);
      second->DrawCopy((param == 0 && sign == 0) ? "" : "SAME");
      second->Clone(Form("ptspectrum_%d_%d", param, sign))->Write();
    }
  }

  canvas->SaveAs(Form("%s.eps", canvas->GetName()));
  canvas2->SaveAs(Form("%s.eps", canvas2->GetName()));

  file->Close();
}

void DrawSystematicpT()
{
  TFile* file = TFile::Open("SystematicpT.root");

  TH1* mcHist2 = (TH1*) file->Get("mymc_unity");
  TH1* result2 = (TH1*) file->Get("result_unity");

  TH1* mcHist[12];
  TH1* result[12];

  Int_t nParams = 5;

  for (Int_t id=0; id<nParams*2; ++id)
  {
    mcHist[id] = (TH1*) file->Get(Form("mymc_%d_%d.root", id / 2, id % 2));
    result[id] = (TH1*) file->Get(Form("result_%d_%d.root", id / 2, id % 2));
  }

  DrawResultRatio(mcHist2, result2, "SystematicpT_OK.eps");

  //DrawRatioDeduct(mcHist2, result2, nParams*2, mcHist, result, "SystematicpT_Summary.eps");

  DrawRatio(nParams*2, mcHist, result, "SystematicpT_Ratios.eps", kTRUE, kTRUE);

  //DrawRatioDeductSmooth(mcHist2, result2, nParams*2, mcHist, result, "SystematicpT_Summary.eps");

  // does not make sense: mc is different
  //Draw2ResultRatio(mcHist, result1, result2, "SystematicpT.eps");
}

void SystematicpT(Bool_t chi2 = 1)
{
  gSystem->Load("libPWG0base");

  TFile::Open("ptspectrum900.root");
  AliMultiplicityCorrection* mult = new AliMultiplicityCorrection("Multiplicity", "Multiplicity");
  mult->LoadHistograms("Multiplicity");

  AliMultiplicityCorrection* mult2 = new AliMultiplicityCorrection("Multiplicity2", "Multiplicity2");

  TH1* mcHist[12];
  TH1* result[12];

  Int_t nParams = 5;

  for (Int_t param=0; param<nParams; param++)
  {
    for (Int_t sign=0; sign<2; sign++)
    {
      // calculate result with systematic effect
      TFile::Open(Form("ptspectrum100_%d_%d.root", param, sign));
      mult2->LoadHistograms("Multiplicity");

      mult->SetMultiplicityESD(etaRange, mult2->GetMultiplicityESD(etaRange));

      if (chi2)
      {
        mult->SetRegularizationParameters(AliMultiplicityCorrection::kPol1, 10000);
        mult->ApplyMinuitFit(etaRange, kFALSE, AliMultiplicityCorrection::kTrVtx);
      }
      else
        mult->ApplyBayesianMethod(etaRange, kFALSE, AliMultiplicityCorrection::kTrVtx, 1, 100, 0);

      Int_t id = param * 2 + sign;

      mcHist[id] = mult2->GetMultiplicityVtx(etaRange)->ProjectionY(Form("mymc_%d_%d.root", param, sign));
      result[id] = (TH1*) mult->GetMultiplicityESDCorrected(etaRange)->Clone(Form("result_%d_%d.root", param, sign));

      TString tmp; tmp.Form("SystematicpT_%d_%d.eps", param, sign);
      DrawResultRatio(mcHist[id], result[id], tmp);
    }
  }

  // calculate normal result
  TFile::Open("ptspectrum100_1.root");
  mult2->LoadHistograms("Multiplicity");
  TH1* mcHist2 = mult2->GetMultiplicityVtx(etaRange)->ProjectionY("mymc_unity");

  mult->SetMultiplicityESD(etaRange, mult2->GetMultiplicityESD(etaRange));

  if (chi2)
  {
    mult->ApplyMinuitFit(etaRange, kFALSE, AliMultiplicityCorrection::kTrVtx);
  }
  else
    mult->ApplyBayesianMethod(etaRange, kFALSE, AliMultiplicityCorrection::kTrVtx, 1, 100);

  TH1* result2 = (TH1*) mult->GetMultiplicityESDCorrected(etaRange)->Clone("result_unity");

  TFile* file = TFile::Open("SystematicpT.root", "RECREATE");
  mcHist2->Write();
  result2->Write();
  for (Int_t id=0; id<nParams*2; ++id)
  {
    mcHist[id]->Write();
    result[id]->Write();
  }
  file->Close();

  DrawSystematicpT();
}

void DrawSystematicpT2()
{
  //displayRange = 200;

  // read from file
  TFile* file = TFile::Open("SystematicpT2.root");
  TH1* mcHist = (TH1*) file->Get("mymc");
  TH1* result[12];
  result[0] = (TH1*) file->Get("result_unity");
  Int_t nParams = 5;
  for (Int_t id=0; id<nParams*2; ++id)
    result[id+1] = (TH1*) file->Get(Form("result_%d_%d", id / 2, id % 2));

  DrawResultRatio((TH1*) mcHist->Clone(), (TH1*) result[0]->Clone(), "SystematicpT_OK.eps");
  DrawRatio(mcHist, nParams*2+1, result, "SystematicpT_Ratios_MC.eps", kTRUE);
  DrawRatio(result[0], nParams*2, result+1, "SystematicpT_Ratios.eps");
}

void SystematicpT2(Bool_t tpc = kTRUE, Bool_t chi2 = kTRUE)
{
  gSystem->Load("libPWG0base");

  if (tpc)
  {
    SetTPC();
    TFile::Open("multiplicityMC_TPC_0.6M_syst_pt_unity.root");
  }
  else
    TFile::Open("ptspectrum100_1.root");

  AliMultiplicityCorrection* measured = new AliMultiplicityCorrection("Multiplicity", "Multiplicity");
  measured->LoadHistograms("Multiplicity");
  TH1* mcHist = measured->GetMultiplicityVtx(etaRange)->ProjectionY("mymc");

  TH1* result[12];

  Int_t nParams = 5;

  // -1 = unity change, 0...4 parameters
  for (Int_t id=-1; id<nParams*2; id++)
  {
    Int_t param = id / 2;
    Int_t sign = id % 2;

    TString idStr;
    if (id == -1)
    {
      idStr = "unity";
    }
    else
      idStr.Form("%d_%d", param, sign);

    // calculate result with systematic effect
    if (tpc)
    {
      TFile::Open(Form("multiplicityMC_TPC_1.3M_syst_pt_%s.root", idStr.Data()));
    }
    else
      TFile::Open(Form("ptspectrum900_%s.root", idStr.Data()));

    AliMultiplicityCorrection* response = new AliMultiplicityCorrection("Multiplicity2", "Multiplicity2");
    response->LoadHistograms("Multiplicity");

    response->SetMultiplicityESD(etaRange, measured->GetMultiplicityESD(etaRange));

    if (chi2)
    {
      response->SetRegularizationParameters(AliMultiplicityCorrection::kPol1, 10000);
      response->ApplyMinuitFit(etaRange, kFALSE, AliMultiplicityCorrection::kTrVtx);
    }
    else
      response->ApplyBayesianMethod(etaRange, kFALSE, AliMultiplicityCorrection::kTrVtx, 1, 100, 0);

    result[id+1] = (TH1*) response->GetMultiplicityESDCorrected(etaRange)->Clone(Form("result_%s", idStr.Data()));

    TString tmp; tmp.Form("SystematicpT_%s.eps", idStr.Data());
    DrawResultRatio(mcHist, result[id+1], tmp);
  }

  TFile* file = TFile::Open("SystematicpT2.root", "RECREATE");
  mcHist->Write();
  for (Int_t id=0; id<nParams*2+1; ++id)
    result[id]->Write();
  file->Close();

  DrawSystematicpT2();
}

void SystematicpTCutOff(Bool_t chi2 = kTRUE)
{
  // only needed for TPC
  SetTPC();

  gSystem->Load("libPWG0base");

  TFile::Open("multiplicityMC_TPC_1.3M_syst_pt_unity.root");
  AliMultiplicityCorrection* mult = new AliMultiplicityCorrection("Multiplicity", "Multiplicity");
  mult->LoadHistograms("Multiplicity");

  TFile::Open("multiplicityMC_TPC_0.6M_syst_pt_unity.root");
  AliMultiplicityCorrection* mult2 = new AliMultiplicityCorrection("Multiplicity2", "Multiplicity2");
  mult2->LoadHistograms("Multiplicity");

  // "normal" result
  mult->SetMultiplicityESD(etaRange, mult2->GetMultiplicityESD(etaRange));

  if (chi2)
  {
    mult->SetRegularizationParameters(AliMultiplicityCorrection::kPol1, 10000);
    mult->ApplyMinuitFit(etaRange, kFALSE, AliMultiplicityCorrection::kTrVtx);
  }
  else
    mult->ApplyBayesianMethod(etaRange, kFALSE, AliMultiplicityCorrection::kTrVtx, 1, 100);

  TH1* mcHist = mult2->GetMultiplicityVtx(etaRange)->ProjectionY("mymc");
  TH1* result1 = (TH1*) mult->GetMultiplicityESDCorrected(etaRange)->Clone("result1");

  TH1* syst[2];

  // change of pt spectrum (down)
  TFile::Open("multiplicityMC_TPC_1.3M_syst_pt_red.root");
  AliMultiplicityCorrection* mult3 = new AliMultiplicityCorrection("Multiplicity3", "Multiplicity3");
  mult3->LoadHistograms("Multiplicity");
  mult3->SetMultiplicityESD(etaRange, mult2->GetMultiplicityESD(etaRange));
  if (chi2)
  {
    mult3->SetRegularizationParameters(AliMultiplicityCorrection::kPol1, 10000);
    mult3->ApplyMinuitFit(etaRange, kFALSE, AliMultiplicityCorrection::kTrVtx);
  }
  else
    mult3->ApplyBayesianMethod(etaRange, kFALSE, AliMultiplicityCorrection::kTrVtx, 1, 100);
  syst[0] = (TH1*) mult3->GetMultiplicityESDCorrected(etaRange)->Clone("result2");

  // change of pt spectrum (up)
  TFile::Open("multiplicityMC_TPC_1.3M_syst_pt_inc.root");
  AliMultiplicityCorrection* mult4 = new AliMultiplicityCorrection("Multiplicity4", "Multiplicity4");
  mult4->LoadHistograms("Multiplicity");
  mult4->SetMultiplicityESD(etaRange, mult2->GetMultiplicityESD(etaRange));
  if (chi2)
  {
    mult4->SetRegularizationParameters(AliMultiplicityCorrection::kPol1, 10000);
    mult4->ApplyMinuitFit(etaRange, kFALSE, AliMultiplicityCorrection::kTrVtx);
  }
  else
    mult4->ApplyBayesianMethod(etaRange, kFALSE, AliMultiplicityCorrection::kTrVtx, 1, 100);
  syst[1] = (TH1*) mult4->GetMultiplicityESDCorrected(etaRange)->Clone("result3");

  DrawRatio(result1, 2, syst, "SystematicpTCutOff.eps", kFALSE, 0, kTRUE);

  Draw2ResultRatio(mcHist, result1, syst[0], "SystematicpTCutOff1.eps");
  Draw2ResultRatio(mcHist, result1, syst[1], "SystematicpTCutOff2.eps");
}

TH1* SystematicsSummary(Bool_t tpc = 1)
{
  Int_t nEffects = 7;

  TH1* effects[10];
  const char** names = 0;
  Int_t colors[] = { 1, 2, 3, 4, 6, 7, 8 };
  Int_t markers[] = { 20, 21, 22, 23, 24, 25, 26 };

  for (Int_t i=0; i<nEffects; ++i)
    effects[i] = new TH1F("SystematicsSummary", ";true multiplicity;Effect", 201, -0.5, 200.5);

  if (tpc)
  {
    SetTPC();

    const char* namesTPC[] = { "Unfolding method (#chi^{2})", "Rel. cross-section", "Particle composition", "p_{t} cut off", "Track selection", "Secondaries", "p_{t} spectrum" };
    names = namesTPC;

    // method
    TFile* file = TFile::Open("StatisticalUncertaintyTPCChi2.root");
    TH1* hist = (TH1*) file->Get("errorBoth");

    // smooth a bit, but skip 0 bin
    effects[0]->SetBinContent(2, hist->GetBinContent(2));
    for (Int_t i=3; i<=200; ++i)
      effects[0]->SetBinContent(i, (hist->GetBinContent(i) + hist->GetBinContent(i+1)) / 2);

    // relative x-section
    effects[1]->SetBinContent(2, 0.005);
    effects[1]->SetBinContent(3, 0.0025);
    effects[1]->SetBinContent(4, 0.0025);

    // particle composition
    for (Int_t i=2; i<=101; ++i)
    {
      if (i < 41)
      {
        effects[2]->SetBinContent(i, 0.01);
      }
      else if (i < 76)
      {
        effects[2]->SetBinContent(i, 0.02);
      }
      else
        effects[2]->SetBinContent(i, 0.02 + 0.08 / 25 * (i - 76));
    }

    // pt cut off (only tpc)
    for (Int_t i=2; i<=101; ++i)
    {
      if (i < 11)
      {
        effects[3]->SetBinContent(i, 0.05);
      }
      else if (i < 51)
      {
        effects[3]->SetBinContent(i, 0.01);
      }
      else
        effects[3]->SetBinContent(i, 0.01 + 0.1 / 30 * (i - 51));
    }

    // track selection (only tpc)
    for (Int_t i=2; i<=101; ++i)
      effects[4]->SetBinContent(i, 0.03);

    // secondaries
    for (Int_t i=2; i<=101; ++i)
      effects[5]->SetBinContent(i, 0.01);

    // pt spectrum
    for (Int_t i=2; i<=101; ++i)
    {
      if (i < 21)
      {
        effects[6]->SetBinContent(i, 0.05);
      }
      else if (i < 51)
      {
        effects[6]->SetBinContent(i, 0.02);
      }
      else
        effects[6]->SetBinContent(i, 0.02 + 0.13 / 25 * (i - 51));
    }

  }
  else
  {
    displayRange = 200;
    nEffects = 5;

    const char* namesSPD[] = { "Unfolding Method (#chi^{2})", "Rel. cross-section", "Particle composition", "Secondaries", "p_{t} spectrum"};
    names = namesSPD;

    // method
    TFile* file = TFile::Open("StatisticalUncertaintySPDChi2.root");
    TH1* hist = (TH1*) file->Get("errorBoth");

    // smooth a bit, but skip 0 bin
    effects[0]->SetBinContent(2, hist->GetBinContent(2));
    for (Int_t i=3; i<=201; ++i)
      effects[0]->SetBinContent(i, (hist->GetBinContent(i) + hist->GetBinContent(i+1)) / 2);

    // relative x-section
    effects[1]->SetBinContent(2, 0.01);
    effects[1]->SetBinContent(3, 0.005);

    // particle composition
    for (Int_t i=2; i<=201; ++i)
    {
      if (i < 6)
      {
        effects[2]->SetBinContent(i, 0.3);
      }
      else if (i < 11)
      {
        effects[2]->SetBinContent(i, 0.05);
      }
      else if (i < 121)
      {
        effects[2]->SetBinContent(i, 0.02);
      }
      else if (i < 151)
      {
        effects[2]->SetBinContent(i, 0.02 + 0.04 / 30 * (i - 121));
      }
      else
        effects[2]->SetBinContent(i, 0.06 + 0.1 / 30 * (i - 151));
    }

    // secondaries
    for (Int_t i=2; i<=201; ++i)
      effects[3]->SetBinContent(i, 0.01);

    // pt spectrum
    for (Int_t i=2; i<=201; ++i)
    {
      if (i < 6)
      {
        effects[4]->SetBinContent(i, 1);
      }
      else if (i < 121)
      {
        effects[4]->SetBinContent(i, 0.03);
      }
      else if (i < 151)
      {
        effects[4]->SetBinContent(i, 0.03 + 0.07 / 30 * (i - 121));
      }
      else
        effects[4]->SetBinContent(i, 0.1);
    }
  }

  TCanvas* canvas = new TCanvas("SystematicsSummary.eps", "SystematicsSummary.eps", 800, 400);
  canvas->SetRightMargin(0.25);
  canvas->SetTopMargin(0.05);
  TLegend* legend = new TLegend(0.2, 0.4, 0.5, 0.4 + 0.5 * nEffects / 7);
  legend->SetFillColor(0);

  for (Int_t i=0; i<nEffects; ++i)
  {
    TH1* current = (TH1*) effects[i]->Clone(Form("current_%d", i));
    /*current->Reset();
    for (Int_t j=0; j<nEffects-i; ++j)
      current->Add(effects[j]);*/

    current->SetLineColor(colors[i]);
    //current->SetFillColor(colors[i]);
    current->SetMarkerColor(colors[i]);
    //current->SetMarkerStyle(markers[i]);

    current->SetStats(kFALSE);
    current->GetYaxis()->SetRangeUser(0, 0.4);
    current->GetXaxis()->SetRangeUser(0, displayRange);
    current->DrawCopy(((i == 0) ? "" : "SAME"));
    legend->AddEntry(current, names[i]);

    TLatex* text = new TLatex(displayRange+5, current->GetBinContent(displayRange+1), names[i]);
    text->SetTextColor(colors[i]);
    text->Draw();
  }

  // add total in square
  TH1* total = (TH1*) effects[0]->Clone("total");
  total->Reset();

  for (Int_t i=0; i<nEffects; ++i)
  {
    //Printf("%d %f", i, effects[i]->GetBinContent(20));
    effects[i]->Multiply(effects[i]);
    total->Add(effects[i]);
  }

  for (Int_t i=1; i<=total->GetNbinsX(); ++i)
    if (total->GetBinContent(i) > 0)
      total->SetBinContent(i, TMath::Min(sqrt(total->GetBinContent(i)), 1.0));

  //Printf("%f", total->GetBinContent(20));

  total->SetMarkerStyle(3);
  total->SetMarkerColor(1);
  legend->AddEntry(total, "total");
  total->DrawCopy("SAME P");

  legend->Draw();

  canvas->SaveAs(canvas->GetName());

  return total;
}

void finalPlot(Bool_t tpc = kTRUE, Bool_t chi2 = kTRUE, Bool_t small = kFALSE)
{
  gSystem->Load("libPWG0base");

  if (tpc)
    SetTPC();

  if (!chi2)
    Printf("WARNING: Bayesian set. This is only for test!");

  // systematic error
  TH1* error = SystematicsSummary(tpc);

  TFile::Open(correctionFile);
  AliMultiplicityCorrection* mult = new AliMultiplicityCorrection("Multiplicity", "Multiplicity");
  mult->LoadHistograms("Multiplicity");

  TFile::Open(measuredFile);
  AliMultiplicityCorrection* mult2 = new AliMultiplicityCorrection("Multiplicity2", "Multiplicity2");
  mult2->LoadHistograms("Multiplicity");

  mult->SetMultiplicityESD(etaRange, mult2->GetMultiplicityESD(etaRange));

  if (chi2)
  {
    mult->SetRegularizationParameters(AliMultiplicityCorrection::kPol1, 10000);
    mult->ApplyMinuitFit(etaRange, kFALSE, AliMultiplicityCorrection::kINEL);
  }
  else
    mult->ApplyBayesianMethod(etaRange, kFALSE, AliMultiplicityCorrection::kINEL, 1, 100, 0, kFALSE);

  TH1* mcHist = mult2->GetMultiplicityINEL(etaRange)->ProjectionY("mymc");
  TH1* result = mult->GetMultiplicityESDCorrected(etaRange);

  DrawResultRatio(mcHist, result, "finalPlotCheck.eps");

  // normalize result
  result->Scale(1.0 / result->Integral(2, 200));

  result->GetXaxis()->SetRangeUser(0, ((tpc) ? displayRange : 200));
  result->SetBinContent(1, 0); result->SetBinError(1, 0);
  result->SetTitle(";true multiplicity;Probability");
  result->SetLineColor(1);
  result->SetStats(kFALSE);

  TH1* systError = (TH1*) result->Clone("systError");
  for (Int_t i=2; i<=systError->GetNbinsX(); ++i)
    systError->SetBinError(i, systError->GetBinContent(i) * error->GetBinContent(i));

  // change error drawing style
  systError->SetFillColor(15);

  TCanvas* canvas = new TCanvas("finalPlot.eps", "finalPlot.eps", (small) ? 600 : 800, 400);
  canvas->SetRightMargin(0.05);
  canvas->SetTopMargin(0.05);

  systError->Draw("E2 ][");
  result->DrawCopy("SAME E ][");
  canvas->SetLogy();

  //TPaveText* text = new TPaveText(10, 1e-3, 50, 1e-4, "B");
  TPaveText* text = new TPaveText(0.15, 0.2, 0.5, 0.4, "B NDC");
  text->SetFillColor(0);
  text->SetTextAlign(12);
  text->AddText("Systematic errors summed quadratically");
  text->AddText("0.6 million minimum bias events");
  text->AddText("corrected to inelastic events");
  text->Draw("B");

  TPaveText* text2 = new TPaveText(0.4, 0.7, 0.6, 0.85, "B NDC");
  text2->SetFillColor(0);
  text2->SetTextAlign(12);
  text2->AddText("#sqrt{s} = 14 TeV");
  if (tpc)
  {
    text2->AddText("|#eta| < 0.9");
  }
  else
    text2->AddText("|#eta| < 2.0");
  text2->AddText("simulated data (PYTHIA)");
  text2->Draw("B");

  if (tpc)
  {
    TText* text3 = new TText(0.75, 0.6, "TPC - full tracking");
    text3->SetNDC();
    text3->Draw();
  }
  else
  {
    TText* text3 = new TText(0.75, 0.6, "SPD - Tracklets");
    text3->SetNDC();
    text3->Draw();
  }

  // alice logo
  TPad* pad = new TPad("pad", "pad", 0.8, 0.7, 0.9, 0.9);
  pad->Draw();
  pad->cd();
  TImage* img = TImage::Open("$HOME/alice.png");
  img->SetImageQuality(TAttImage::kImgBest);
  img->Draw();

  canvas->Modified();

/*  TText* text = new TText(10, 1e-4, "Systematic errors summed quadratically");
  text->SetTextSize(0.04);
  text->DrawText(10, 5e-5, "0.6 #cdot 10^{6} minimum bias events");
  text->DrawText(10, 3e-5, "TPC tracks in |#eta| < 0.9");
  text->DrawText(10, 1e-5, "corrected to ineleastic events in |#eta| < 0.9");
  text->Draw();*/


  canvas->SaveAs(canvas->GetName());
}

void BlobelUnfoldingExample()
{
  const Int_t kSize = 20;

  TMatrixD matrix(kSize, kSize);
  for (Int_t x=0; x<kSize; x++)
  {
    for (Int_t y=0; y<kSize; y++)
    {
      if (x == y)
      {
        if (x == 0 || x == kSize -1)
        {
          matrix(x, y) = 0.75;
        }
        else
          matrix(x, y) = 0.5;
      }
      else if (TMath::Abs(x - y) == 1)
      {
        matrix(x, y) = 0.25;
      }
    }
  }

  //matrix.Print();

  TMatrixD inverted(matrix);
  inverted.Invert();

  //inverted.Print();

  TH1F* inputDist = new TH1F("inputDist", ";t;#tilde{T}(t)", kSize, -0.5, (Float_t) kSize - 0.5);
  TVectorD inputDistVector(kSize);
  TH1F* unfolded = inputDist->Clone("unfolded");
  TH1F* measuredIdealDist = inputDist->Clone("measuredIdealDist");
  measuredIdealDist->SetTitle(";m;#tilde{M}(m)");
  TH1F* measuredDist = measuredIdealDist->Clone("measuredDist");

  TF1* gaus = new TF1("func", "gaus(0)", -0.5, kSize);
  // norm: 1/(sqrt(2pi)sigma)
  gaus->SetParameters(10000 / sqrt(2 * TMath::Pi()) / ((Float_t) kSize / 8), (Float_t) kSize / 2, (Float_t) kSize / 8);
  //gaus->Print();

  for (Int_t x=1; x<=inputDist->GetNbinsX(); x++)
  {
    Float_t value = gaus->Eval(inputDist->GetBinCenter(x));
    inputDist->SetBinContent(x, value);
    inputDistVector(x-1) = value;
  }

  TVectorD measuredDistIdealVector = matrix * inputDistVector;
  
  for (Int_t x=1; x<=measuredIdealDist->GetNbinsX(); x++)
    measuredIdealDist->SetBinContent(x, measuredDistIdealVector(x-1));

  measuredDist->FillRandom(measuredIdealDist, 10000);

  // fill error matrix before scaling
  TMatrixD covarianceMatrixMeasured(kSize, kSize);
  for (Int_t x=1; x<=unfolded->GetNbinsX(); x++)
    covarianceMatrixMeasured(x-1, x-1) = TMath::Sqrt(measuredDist->GetBinContent(x));

  TMatrixD covarianceMatrix = inverted * covarianceMatrixMeasured * inverted;
  //covarianceMatrix.Print();

  TVectorD measuredDistVector(kSize);
  for (Int_t x=1; x<=measuredDist->GetNbinsX(); x++)
    measuredDistVector(x-1) = measuredDist->GetBinContent(x);

  TVectorD unfoldedVector = inverted * measuredDistVector;
  for (Int_t x=1; x<=unfolded->GetNbinsX(); x++)
    unfolded->SetBinContent(x, unfoldedVector(x-1));

  TCanvas* canvas = new TCanvas("BlobelUnfoldingExample", "BlobelUnfoldingExample", 1000, 500);
  canvas->SetTopMargin(0.05);
  canvas->Divide(2, 1);

  canvas->cd(1);
  canvas->cd(1)->SetLeftMargin(0.15);
  canvas->cd(1)->SetRightMargin(0.05);
  measuredDist->GetYaxis()->SetTitleOffset(1.7);
  measuredDist->SetStats(0);
  measuredDist->DrawCopy();
  gaus->Draw("SAME");

  canvas->cd(2);
  canvas->cd(2)->SetLeftMargin(0.15);
  canvas->cd(2)->SetRightMargin(0.05);
  unfolded->GetYaxis()->SetTitleOffset(1.7);
  unfolded->SetStats(0);
  unfolded->DrawCopy();
  gaus->Draw("SAME");

  canvas->SaveAs("BlobelUnfoldingExample.eps");
}

void E735Fit()
{
  TH1* fCurrentESD = new TH1F("mult", "mult", 501, -0.5, 500.5);
  fCurrentESD->Sumw2();

  // Open the input stream
  ifstream in;
  in.open("e735data.txt");

  while(in.good())
  {
    Float_t x, y, ye;
    in >> x >> y >> ye;

    //Printf("%f %f %f", x, y, ye);
    fCurrentESD->SetBinContent(fCurrentESD->FindBin(x), y);
    fCurrentESD->SetBinError(fCurrentESD->FindBin(x), ye);
  }

  in.close();

  //new TCanvas; fCurrentESD->DrawCopy(); gPad->SetLogy();

  fCurrentESD->Scale(1.0 / fCurrentESD->Integral());

  TF1* func = new TF1("nbd", "[0] * TMath::Binomial([2]+TMath::Nint(x)-1, [2]-1) * pow([1] / ([1]+[2]), TMath::Nint(x)) * pow(1 + [1]/[2], -[2])");
  func->SetParNames("scaling", "averagen", "k");
  func->SetParLimits(0, 0.001, fCurrentESD->GetMaximum() * 1000);
  func->SetParLimits(1, 0.001, 1000);
  func->SetParLimits(2, 0.001, 1000);
  func->SetParameters(fCurrentESD->GetMaximum() * 100, 10, 2);

  TF1* lognormal = new TF1("lognormal", "[0]*exp(-(log(x)-[1])^2/(2*[2]^2))/(x*[2]*TMath::Sqrt(2*TMath::Pi()))", 0.01, 500);
  lognormal->SetParNames("scaling", "mean", "sigma");
  lognormal->SetParameters(1, 1, 1);
  lognormal->SetParLimits(0, 0, 10);
  lognormal->SetParLimits(1, 0, 100);
  lognormal->SetParLimits(2, 1e-3, 10);

  TCanvas* canvas = new TCanvas("c1", "c1", 700, 400);
  fCurrentESD->SetStats(kFALSE);
  fCurrentESD->GetYaxis()->SetTitleOffset(1.3);
  fCurrentESD->SetTitle(";true multiplicity (N);P_{N}");
  fCurrentESD->Draw("");
  fCurrentESD->GetXaxis()->SetRangeUser(0, 250);
  fCurrentESD->Fit(func, "0", "", 0, 150);
  func->SetRange(0, 250);
  func->Draw("SAME");
  printf("chi2 = %f\n", func->GetChisquare());

  fCurrentESD->Fit(lognormal, "0", "", 0.01, 150);
  lognormal->SetLineColor(2);
  lognormal->SetLineStyle(2);
  lognormal->SetRange(0, 250);
  lognormal->Draw("SAME");

  gPad->SetLogy();

  canvas->SaveAs("E735Fit.eps");
}
