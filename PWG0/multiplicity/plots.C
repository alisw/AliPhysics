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

const char* correctionFile = "multiplicityMC.root";
const char* measuredFile   = "multiplicityESD.root";
Int_t etaRange = 1;
Int_t displayRange = 80; // axis range
Int_t ratioRange = 151;   // range to calculate difference
Int_t longDisplayRange = 120;

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

const char* GetMultLabel(Int_t etaR = -1, Bool_t trueM = kTRUE)
{
	if (etaR == -1)
		etaR = etaRange;
		
	TString tmpStr((trueM) ? "True " : "Measured ");

        tmpStr += "multiplicity";
        //return Form("%s", tmpStr.Data());
		
	if (etaR == 4)
	{
		tmpStr += " (full phase space)";
	}
	else if (etaR == 2)
	{
		tmpStr += " in |#eta| < 1.3";
	}
	else
		tmpStr += Form(" in |#eta| < %.1f", (etaR+1)* 0.5);
	return Form("%s", tmpStr.Data());
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

TH1* responseMatrixPlot(const char* fileName = 0, const char* folder = "Multiplicity")
{
  loadlibs();

  AliMultiplicityCorrection* mult = new AliMultiplicityCorrection(folder, folder);

  if (fileName == 0)
    fileName = correctionFile;
  
  TFile::Open(fileName);
  mult->LoadHistograms();

  // empty under/overflow bins in x, otherwise Project3D takes them into account
  TH2* hist = (TH2*) mult->GetCorrelation(etaRange);
  for (Int_t y=0; y<=hist->GetYaxis()->GetNbins()+1; ++y)
  {
    for (Int_t z=0; z<=hist->GetZaxis()->GetNbins()+1; ++z)
    {
      hist->SetBinContent(0, y, z, 0);
      hist->SetBinContent(hist->GetXaxis()->GetNbins()+1, y, z, 0);
    }
  }
  
  hist = (TH2*) (((TH3*) mult->GetCorrelation(etaRange))->Project3D("zy"));
  hist->SetStats(kFALSE);
  
  if (0)
  {
    // normalize
    // normalize correction for given nPart
    for (Int_t i=1; i<=hist->GetNbinsX(); ++i)
    {
      Double_t sum = hist->Integral(i, i, 1, hist->GetNbinsY());
      if (sum <= 0)
        continue;

      for (Int_t j=1; j<=hist->GetNbinsY(); ++j)
      {
        // npart sum to 1
        hist->SetBinContent(i, j, hist->GetBinContent(i, j) / sum);// * correlation->GetXaxis()->GetBinWidth(i));
        hist->SetBinError(i, j, hist->GetBinError(i, j) / sum);
      }
    }
  }
  
  hist->SetTitle(Form(";%s;%s;Entries", GetMultLabel(), GetMultLabel(etaRange, kFALSE)));
  hist->GetXaxis()->SetRangeUser(0, longDisplayRange);
  hist->GetYaxis()->SetRangeUser(0, longDisplayRange);
  
  hist->GetYaxis()->SetTitleOffset(1.3);
  hist->GetZaxis()->SetTitleOffset(1.2);

  TCanvas* canvas = new TCanvas("c1", "c1", 600, 600);
  canvas->SetRightMargin(0.15);
  canvas->SetTopMargin(0.05);

  gPad->SetLogz();
  hist->Draw("COLZ");

  canvas->SaveAs("responsematrix.eps");
  
  return hist;
}

void multPythiaPhojet()
{
  loadlibs();
  
  TFile::Open("LHC08c11_10TeV_0.5T/mb1/spd/multiplicity.root");
  AliMultiplicityCorrection* mult = new AliMultiplicityCorrection("Multiplicity", "Multiplicity");
  mult->LoadHistograms("Multiplicity");
  hist1 = mult->GetMultiplicityINEL(1)->ProjectionY();
  hist1->Sumw2();
  
  TFile::Open("LHC08c15_10TeV_0.5T_Phojet/mb1/spd/multiplicity.root");
  AliMultiplicityCorrection* mult = new AliMultiplicityCorrection("Multiplicity", "Multiplicity");
  mult->LoadHistograms("Multiplicity");
  hist2 = mult->GetMultiplicityINEL(1)->ProjectionY();
  hist2->Sumw2();
  
  legend = new TLegend(0.6, 0.7, 0.9, 0.9);
  legend->SetFillColor(0);
  legend->SetTextSize(0.04);
  legend->AddEntry(hist1, "Pythia", "L");
  legend->AddEntry(hist2, "Phojet", "L");

  c1 = new TCanvas("c", "c", 600, 600);
  c1->SetTopMargin(0.05);
  c1->SetRightMargin(0.05);
  c1->SetLeftMargin(0.12);
  c1->SetGridx();
  c1->SetGridy(); 
  c1->SetLogy();
  
  //hist1->SetMarkerStyle(20);
  //hist2->SetMarkerStyle(24);
  //hist2->SetMarkerColor(2);
  hist1->SetLineWidth(2);
  hist2->SetLineWidth(2);
  hist2->SetLineStyle(2);
  hist2->SetLineColor(2);
  
  hist1->Scale(1.0 / hist1->Integral());
  hist2->Scale(1.0 / hist2->Integral());
  
  hist1->SetStats(0);
  hist1->GetYaxis()->SetTitleOffset(1.3);
  hist1->GetXaxis()->SetRangeUser(0, 100);
  hist1->SetTitle(";N_{ch};P(N_{ch})");
  
  hist1->Draw("");
  hist2->Draw("SAME");
  legend->Draw();
  
  c1->SaveAs("mult_pythia_phojet.eps");
}

TCanvas* DrawResultRatio(TH1* mcHist, TH1* result, TString epsName)
{
  // normalize unfolded result to mc hist
  result->Scale(1.0 / result->Integral(2, displayRange));
  result->Scale(mcHist->Integral(2, displayRange));

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
  pad1->SetGridx();
  pad1->SetGridy();

  mcHist->GetXaxis()->SetLabelSize(0.06);
  mcHist->GetYaxis()->SetLabelSize(0.06);
  mcHist->GetXaxis()->SetTitleSize(0.06);
  mcHist->GetYaxis()->SetTitleSize(0.06);
  mcHist->GetYaxis()->SetTitleOffset(0.6);

  mcHist->GetXaxis()->SetRangeUser(0, displayRange);

  mcHist->SetTitle(Form(";%s;Entries", GetMultLabel()));
  mcHist->SetStats(kFALSE);

  mcHist->DrawCopy("HIST E");
  gPad->SetLogy();

  result->SetLineColor(2);
  result->SetMarkerColor(2);
  result->SetMarkerStyle(5);
  result->DrawCopy("SAME PE");

  TLegend* legend = new TLegend(0.6, 0.65, 0.95, 0.9);
  legend->AddEntry(mcHist, "True distribution");
  legend->AddEntry(result, "Unfolded distribution", "P");
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
  result1->Scale(1.0 / result1->Integral(2, displayRange));
  result1->Scale(mcHist->Integral(2, displayRange));
  result2->Scale(1.0 / result2->Integral(2, displayRange));
  result2->Scale(mcHist->Integral(2, displayRange));

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
  gPad->SetGridx();
  gPad->SetGridy();

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
  result1->SetMarkerStyle(24);
  result1->DrawCopy("SAME HISTE");

  result2->SetLineColor(4);
  result2->SetMarkerColor(4);
  result2->DrawCopy("SAME HISTE");

  TLegend* legend = new TLegend(0.5, 0.6, 0.95, 0.9);
  legend->AddEntry(mcHist, "True distribution");
  legend->AddEntry(result1, "Unfolded distribution (syst)", "P");
  legend->AddEntry(result2, "Unfolded distribution (normal)", "P");
  legend->SetFillColor(0);
  legend->SetTextSize(0.06);
  legend->Draw();

  pad2->cd();
  pad2->SetBottomMargin(0.15);
  //gPad->SetGridx();
  //gPad->SetGridy();

  result1->GetXaxis()->SetLabelSize(0.06);
  result1->GetYaxis()->SetLabelSize(0.06);
  result1->GetXaxis()->SetTitleSize(0.06);
  result1->GetYaxis()->SetTitleSize(0.06);
  result1->GetYaxis()->SetTitleOffset(0.6);

  result1->GetXaxis()->SetRangeUser(0, displayRange);

  result1->SetTitle(Form(";%s;Entries", GetMultLabel()));
  result1->SetStats(kFALSE);

  // calculate ratio
  result1->Sumw2();
  TH1* ratio = (TH1*) result1->Clone("ratio");
  result2->Sumw2();
  ratio->Divide(ratio, result2, 1, 1, "");
  ratio->SetLineColor(1);
  ratio->SetMarkerColor(1);
  ratio->SetMarkerStyle(0);
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

TCanvas* Draw2ResultRatios(TH1* mcHist, TH1* result1, TH1* result2, TH1* ratio2, TH1* ratio3, TString epsName)
{
  // draws the 3 plots in the upper plot
  // draws the ratio between result1 and result2 in the lower plot
  // also draws ratio2 and ratio3 in the lower plot, uses their name for the legend

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
  gPad->SetGridx();
  gPad->SetGridy();

  mcHist->GetXaxis()->SetLabelSize(0.06);
  mcHist->GetYaxis()->SetLabelSize(0.06);
  mcHist->GetXaxis()->SetTitleSize(0.06);
  mcHist->GetYaxis()->SetTitleSize(0.06);
  mcHist->GetYaxis()->SetTitleOffset(0.6);

  mcHist->GetXaxis()->SetRangeUser(0, displayRange);

  mcHist->SetTitle(";True multiplicity;Entries");
  mcHist->SetStats(kFALSE);

  mcHist->DrawCopy("HIST E");
  gPad->SetLogy();

  result1->SetLineColor(2);
  result1->SetMarkerColor(2);
  result1->SetMarkerStyle(24);
  result1->DrawCopy("SAME HISTE");

  result2->SetLineColor(4);
  result2->SetMarkerColor(4);
  result2->SetMarkerStyle(5);
  result2->DrawCopy("SAME HISTE");

  TLegend* legend = new TLegend(0.55, 0.6, 0.95, 0.9);
  legend->AddEntry(mcHist, "True distribution");
  legend->AddEntry(result1, "Unfolded distribution (5%)", "P");
  legend->AddEntry(result2, "Unfolded distribution (normal)", "P");
  legend->SetFillColor(0);
  legend->SetTextSize(0.06);
  legend->Draw();

  pad2->cd();
  pad2->SetBottomMargin(0.15);
  //gPad->SetGridx();
  //gPad->SetGridy();

  result1->GetXaxis()->SetLabelSize(0.06);
  result1->GetYaxis()->SetLabelSize(0.06);
  result1->GetXaxis()->SetTitleSize(0.06);
  result1->GetYaxis()->SetTitleSize(0.06);
  result1->GetYaxis()->SetTitleOffset(0.6);

  result1->GetXaxis()->SetRangeUser(0, displayRange);

  result1->SetTitle(Form(";%s;Entries", GetMultLabel()));
  result1->SetStats(kFALSE);

  // calculate ratio
  result1->Sumw2();
  TH1* ratio = (TH1*) result1->Clone("ratio");
  result2->Sumw2();
  ratio->Divide(ratio, result2, 1, 1, "");
  ratio->SetLineColor(1);
  ratio->SetMarkerColor(1);
  ratio->SetMarkerStyle(24);
  ratio->GetYaxis()->SetTitle("Ratio (change / normal)");
  ratio->GetYaxis()->SetRangeUser(0.55, 1.45);

  ratio2->SetLineColor(2);
  ratio2->SetMarkerColor(2);
  ratio2->SetMarkerStyle(25);

  ratio3->SetLineColor(4);
  ratio3->SetMarkerColor(4);
  ratio3->SetMarkerStyle(26);
  
  ratio->DrawCopy();
  ratio2->DrawCopy("SAME");
  ratio3->DrawCopy("SAME");
  
  legend2 = new TLegend(0.3, 0.8, 0.8, 0.95);
  legend2->SetNColumns(3);
  legend2->SetFillColor(0);
  legend2->SetTextSize(0.06);
  legend2->AddEntry(ratio, "5% change", "P");
  legend2->AddEntry(ratio2, "2% change", "P");
  legend2->AddEntry(ratio3, "1% change", "P");
  legend2->Draw();

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

void DrawEfficiencyChange()
{
  loadlibs();

  const char* fileName = "chi2a_inel.root";

  base = AliMultiplicityCorrection::Open(Form("spd/%s", fileName));
  low = AliMultiplicityCorrection::Open(Form("spd-loweff/%s", fileName));
  high = AliMultiplicityCorrection::Open(Form("spd-higheff/%s", fileName));
  
  const char* legendStrings[] = { "low efficiency", "high efficiency" };
  
  file = TFile::Open("systunc_detectorefficiency.root", "RECREATE");
  
  for (Int_t etaR = 1; etaR < 2; etaR++)
  {
    base->GetMultiplicityESDCorrected(etaR)->Scale(1.0 / base->GetMultiplicityESDCorrected(etaR)->Integral(2, base->GetMultiplicityESDCorrected(etaR)->GetNbinsX()));
  
    TH1* hists[2];
    hists[0] = low->GetMultiplicityESDCorrected(etaR);
    hists[0]->Scale(1.0 / hists[0]->Integral(2, hists[0]->GetNbinsX()));

    if (high)
    {
      hists[1] = high->GetMultiplicityESDCorrected(etaR);
      hists[1]->Scale(1.0 / hists[1]->Integral(2, hists[1]->GetNbinsX()));
    }
    
    largestError = (TH1*) hists[0]->Clone(Form("detectorefficiency_%d", etaR));
    largestError->Reset();
  
    DrawRatio(base->GetMultiplicityESDCorrected(etaR), (high) ? 2 : 1, hists, Form("eff_%d.png", etaR), kFALSE, legendStrings, kFALSE, largestError);
    
    largestError->Write();
    
    new TCanvas;
    largestError->DrawCopy();
  }
  
  file->Close();
}

TCanvas* DrawRatio(TH1* result, Int_t nResultSyst, TH1** resultSyst, TString epsName, Bool_t firstMarker = kFALSE, const char** legendStrings = 0, Bool_t errors = kFALSE, TH1* largestErrorLow = 0, TH1* largestErrorHigh = 0)
{
  // compares n results with first results. E.g. one gained with the default response, another with a changed one to study
  // a systematic effect

  // normalize results
  //result->Scale(1.0 / result->Integral(1, 200));

  TCanvas* canvas = new TCanvas(epsName, epsName, 800, 500);
  canvas->SetTopMargin(0.05);
  canvas->SetRightMargin(0.05);
  canvas->SetGridx();
  canvas->SetGridy();

  result->GetXaxis()->SetRangeUser(0, displayRange);
  result->GetYaxis()->SetRangeUser(0.55, 1.45);
  result->SetStats(kFALSE);

  // to get the axis how we want it
  TH1* dummy = (TH1*) result->Clone("dummy");
  dummy->Reset();
  dummy->SetTitle(Form(";%s;Ratio", GetMultLabel()));
  dummy->DrawCopy();
  delete dummy;

  Int_t colors[] = {1, 2, 4, 6, 7, 8, 9, 10};

  TLegend* legend = new TLegend(0.2, 0.7, 0.7, 0.93);
  legend->SetFillColor(0);
  //legend->SetTextSize(0.04);
  if (nResultSyst > 6)
    legend->SetNColumns(2);

  for (Int_t n=0; n<nResultSyst; ++n)
  {
    //resultSyst[n]->Scale(1.0 / resultSyst[n]->Integral(1, 200));

    // calculate ratio
    TH1* ratio = (TH1*) result->Clone("ratio");
    ratio->Divide(ratio, resultSyst[n], 1, 1, "");
    ratio->GetXaxis()->SetRangeUser(0, displayRange);

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
      legend->AddEntry(ratio, legendStrings[n], "L");
      
    /*
    s = new TSpline3(ratio);
    s->SetNpx(5);
    s->SetLineColor(kRed);
    s->Draw("same");
    */
    
    if (largestErrorLow && largestErrorHigh)
      for (Int_t i=1; i<=ratio->GetNbinsX(); ++i)
      {
        if (ratio->GetBinContent(i) - 1 > 0)
          largestErrorHigh->SetBinContent(i, TMath::Max(ratio->GetBinContent(i) - 1, largestErrorHigh->GetBinContent(i)));
        if (ratio->GetBinContent(i) - 1 < 0)
          largestErrorLow->SetBinContent(i, TMath::Min(ratio->GetBinContent(i) - 1, largestErrorLow->GetBinContent(i)));
      }

    if (largestErrorLow && !largestErrorHigh)
      for (Int_t i=1; i<=ratio->GetNbinsX(); ++i)
        largestErrorLow->SetBinContent(i, TMath::Max(TMath::Abs(ratio->GetBinContent(i) - 1), largestErrorLow->GetBinContent(i)));

    // get average of ratio
    Float_t sum = 0;
    for (Int_t i=2; i<=ratioRange; ++i)
      sum += TMath::Abs(ratio->GetBinContent(i) - 1);
    sum /= ratioRange-1;

    printf("%d) Average (2..%d) of |ratio - 1| is %f\n", n, ratioRange, sum);
  }

  if (legendStrings)
    legend->Draw();

  TLine* line = new TLine(-0.5, 1, displayRange, 1);
  line->SetLineWidth(2);
  line->Draw();

  line = new TLine(-0.5, 1.1, displayRange, 1.1);
  line->SetLineWidth(2);
  line->SetLineStyle(2);
  line->Draw();
  line = new TLine(-0.5, 0.9, displayRange, 0.9);
  line->SetLineWidth(2);
  line->SetLineStyle(2);
  line->Draw();

  canvas->SaveAs(epsName);

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

void DrawResiduals(const char* fileName, const char* epsName)
{
  loadlibs();

  AliMultiplicityCorrection* mult = AliMultiplicityCorrection::Open(fileName);

  TH1* measured = mult->GetMultiplicityESD(etaRange)->ProjectionY("myesd", 1, 1);
  
  if (0)
  {
    // test how far we are from a normal exponential in the unfolded
    corrected = mult->GetMultiplicityESDCorrected(etaRange);
    new TCanvas; corrected->DrawCopy(); gPad->SetLogy();
    
    expo = new TF1("exp", "expo(0)", 0, 50);
    //expo->SetParameters(10, -0.18);
    expo->SetParameters(9.96, -0.176);
    expo->Draw("SAME");  
    
    for (Int_t i=21; i<=50; i++)
      corrected->SetBinContent(i, expo->Eval(corrected->GetXaxis()->GetBinCenter(i)));
    
    corrected->SetLineColor(2);
    corrected->Draw("SAME");
  }
  
  TH1* unfoldedFolded = mult->GetConvoluted(etaRange, AliMultiplicityCorrection::kINEL);
  //mult->CalculateMultiplicityESD(mult->GetMultiplicityESDCorrected(etaRange), etaRange)->ProjectionY("myfolded", 1, 1);
  
  // normalize
  //unfoldedFolded->Scale(1.0 / unfoldedFolded->Integral(2, displayRange+1));
  //unfoldedFolded->Scale(measured->Integral(2, displayRange+1));

  TCanvas* canvas = new TCanvas(epsName, epsName, 800, 600);
  canvas->Range(0, 0, 1, 1);

  TPad* pad1 = new TPad(Form("%s_pad1", epsName), "", 0, 0.5, 1, 1);
  pad1->Draw();
  pad1->SetGridx();
  pad1->SetGridy();

  TPad* pad2 = new TPad(Form("%s_pad2", epsName), "", 0, 0.02, 1, 0.5);
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

  measured->GetXaxis()->SetRangeUser(0, displayRange);

  measured->SetTitle(Form(";%s;Entries", GetMultLabel(etaRange, kFALSE)));
  measured->SetStats(kFALSE);

  measured->DrawCopy("HIST");
  gPad->SetLogy();

  unfoldedFolded->SetMarkerStyle(5);
  unfoldedFolded->SetMarkerColor(2);
  unfoldedFolded->SetLineColor(2);
  unfoldedFolded->DrawCopy("SAME PHIST");

  TLegend* legend = new TLegend(0.6, 0.65, 0.95, 0.9);
  legend->AddEntry(measured, "Measured distribution", "L");
  legend->AddEntry(unfoldedFolded, "R #otimes unfolded distribution", "P");
  legend->SetFillColor(0);
  legend->SetTextSize(0.06);
  legend->Draw();

  pad2->cd();
  pad2->SetBottomMargin(0.15);

  // calculate ratio
  measured->Sumw2();
  TH1* residual = (TH1*) measured->Clone("residual");
  unfoldedFolded->Sumw2();

  residual->Add(unfoldedFolded, -1);

  // projection
  TH1* residualHist = new TH1F("residualHist", ";", 16, -3, 3);
  residualHist->Sumw2();

  Float_t chi2 = 0;
  Float_t chi2_hump = 0;
  
  for (Int_t i=1; i<=displayRange+1; ++i)
  {
    if (measured->GetBinError(i) > 0)
    {
      residual->SetBinContent(i, residual->GetBinContent(i) / measured->GetBinError(i));
      residual->SetBinError(i, 1);

      residualHist->Fill(residual->GetBinContent(i));
      if (i >= 15 && i <= 23)
        chi2_hump += residual->GetBinContent(i) * residual->GetBinContent(i);

      chi2 += residual->GetBinContent(i) * residual->GetBinContent(i);
    }
    else
    {
      residual->SetBinContent(i, 0);
      residual->SetBinError(i, 0);
    }
  }
  
  Printf("chi2 / ndf = %f / %d = %f", chi2, displayRange+1, chi2 / (displayRange+1));
  Printf("chi2 (hump) / ndf = %f / %d = %f", chi2_hump, 23-15+1, chi2_hump / (23-15+1));

  residual->GetYaxis()->SetTitle("Residuals:   (1/e) (M - R  #otimes U)");
  residual->GetYaxis()->SetRangeUser(-4.5, 4.5);
  residual->DrawCopy();

  TLine* line = new TLine(-0.5, 0, displayRange + 0.5, 0);
  line->SetLineWidth(2);
  line->Draw();

  pad3->cd();
  residualHist->SetStats(kFALSE);
  residualHist->SetLabelSize(0.08, "xy");
  residualHist->Fit("gaus");
  residualHist->Draw("HIST");
  residualHist->FindObject("gaus")->Draw("SAME");

  canvas->Modified();
  canvas->SaveAs(canvas->GetName());

  //const char* epsName2 = "proj.eps";
  //TCanvas* canvas = new TCanvas(epsName2, epsName2, 800, 600);
  //canvas->SetGridx();
  //canvas->SetGridy();

  //canvas->SaveAs(canvas->GetName());
}

void chi2FluctuationResult()
{
  loadlibs();

  TFile::Open("chi2_noregularization.root");
  AliMultiplicityCorrection* mult = new AliMultiplicityCorrection("Multiplicity", "Multiplicity");
  mult->LoadHistograms("Multiplicity");

  TH1* mcHist = mult->GetMultiplicityVtx(etaRange)->ProjectionY("mymc", 1, 1);

  mult->DrawComparison("MinuitChi2", etaRange, kFALSE, kTRUE, mcHist, kTRUE);

  TCanvas* canvas = (TCanvas*) gROOT->FindObject("MinuitChi2_DrawComparison_1");
  ((TH1*) canvas->FindObject("proj"))->GetXaxis()->SetRangeUser(0, displayRange);
  ((TH1*) canvas->FindObject("fCurrentESD"))->GetXaxis()->SetRangeUser(0, displayRange);
  //((TH1*) canvas->FindObject("proj"))->GetXaxis()->SetTitle(GetMultTitle());
  //((TH1*) canvas->FindObject("fCurrentESD"))->GetXaxis()->->SetTitle(GetMultTitle(etaRange, kFALSE));
  canvas->SaveAs("chi2FluctuationResult.eps");
}

void DrawUnfolded(const char* fileName, const char* eps)
{
  loadlibs();
  
  TFile::Open(fileName);

  AliMultiplicityCorrection* mult = new AliMultiplicityCorrection("Multiplicity", "Multiplicity");
  mult->LoadHistograms("Multiplicity");

  TH1* mcHist = mult->GetMultiplicityVtx(etaRange)->ProjectionY("mymc", 1, mult->GetMultiplicityVtx(etaRange)->GetNbinsX());
  TH1* result = mult->GetMultiplicityESDCorrected(etaRange);

  DrawResultRatio(mcHist, result, eps);
}

void minimizationInfluenceAlpha()
{
  loadlibs();

  TFile::Open(measuredFile);
  AliMultiplicityCorrection* mult2 = new AliMultiplicityCorrection("Multiplicity2", "Multiplicity2");
  mult2->LoadHistograms("Multiplicity");

  TH1* mcHist = mult2->GetMultiplicityVtx(etaRange)->ProjectionY("mymc", 1, mult2->GetMultiplicityVtx(etaRange)->GetNbinsX());
  mcHist->Scale(1.0 / mcHist->Integral());
  mcHist->SetStats(kFALSE);
  mcHist->SetTitle(";True multiplicity n in |#eta| < 1.0;P(n)");

  TCanvas* canvas = new TCanvas("minimizationInfluenceAlpha", "minimizationInfluenceAlpha", 1000, 350);
  canvas->Divide(3, 1, 0.005);

  TFile::Open("chi2compare-influencealpha/EvaluateChi2Method.root");
  
  TH1* hist1 = (TH1*) gFile->Get("MinuitChi2_00_2_3.162278");
  TH1* hist2 = (TH1*) gFile->Get("MinuitChi2_07_2_10000.000000");
  TH1* hist3 = (TH1*) gFile->Get("MinuitChi2_13_2_10000000.000000");

  /*mcHist->Rebin(2);  mcHist->Scale(0.5);
  hist1->Rebin(2);   hist1->Scale(0.5);
  hist2->Rebin(2);   hist2->Scale(0.5);
  hist3->Rebin(2);   hist3->Scale(0.5);*/

  mcHist->GetXaxis()->SetRangeUser(0, displayRange);
  mcHist->SetLabelSize(0.06, "xy");
  mcHist->SetTitleSize(0.06, "xy");
  mcHist->GetYaxis()->SetTitleOffset(1.5);
  
  canvas->cd(1);
  
  gPad->SetLogy();
  gPad->SetRightMargin(0.03);
  gPad->SetLeftMargin(0.19);
  gPad->SetTopMargin(0.05);
  gPad->SetBottomMargin(0.13);
  gPad->SetGridx();
  gPad->SetGridy();
  mcHist->Draw();
  hist1->SetMarkerStyle(5);
  hist1->SetMarkerColor(2);
  hist1->Draw("SAME PE");

  canvas->cd(2);
  gPad->SetRightMargin(0.03);
  gPad->SetLeftMargin(0.19);
  gPad->SetTopMargin(0.05);
  gPad->SetBottomMargin(0.13);
  gPad->SetLogy();
  gPad->SetGridx();
  gPad->SetGridy();
  mcHist->Draw();
  hist2->SetMarkerStyle(5);
  hist2->SetMarkerColor(2);
  hist2->Draw("SAME PE");

  canvas->cd(3);
  gPad->SetRightMargin(0.03);
  gPad->SetLeftMargin(0.19);
  gPad->SetTopMargin(0.05);
  gPad->SetBottomMargin(0.13);
  gPad->SetLogy();
  gPad->SetGridx();
  gPad->SetGridy();
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

void StartingConditions()
{
  // data generated by runMultiplicitySelector.C StartingConditions

  const char* name = "StartingConditions";

  TFile* file = TFile::Open(Form("%s.root", name));

  TCanvas* canvas = new TCanvas(name, name, 1200, 600);
  canvas->Divide(2, 1);

  TH1* mc = (TH1*) file->Get("mc");
  mc->Sumw2();
  mc->Scale(1.0 / mc->Integral(2, displayRange));

  //Int_t marker[] = {24, 25, 26, 27, 28, 2, 3, 4, 5};

  TLegend* legend = new TLegend(0.6, 0.7, 0.99, 0.99);
  legend->SetFillColor(0);
  legend->SetTextSize(0.04);
  
  const char* names[] = { "True", "Measured 1", "Measured 2", "Measured 3", "NBD", "Flat" };

  for (Int_t i=0; i<6; ++i)
  {
    Int_t id = i;
    if (id > 2)
      id += 2;

    TH1* chi2Result = (TH1*) file->Get(Form("chi2Result_%d", id));
    TH1* bayesResult = (TH1*) file->Get(Form("bayesResult_%d", id));
    
    chi2Result->Scale(1.0 / chi2Result->Integral(2, displayRange));
    bayesResult->Scale(1.0 / bayesResult->Integral(2, displayRange));

    chi2Result->Divide(mc, chi2Result, 1, 1, "");
    bayesResult->Divide(mc, bayesResult, 1, 1, "");

    chi2Result->SetTitle(Form("a) #chi^{2}-minimization;%s;MC / unfolded", GetMultLabel()));
    chi2Result->GetXaxis()->SetRangeUser(0, displayRange);
    chi2Result->GetYaxis()->SetRangeUser(0.7, 1.3);
    chi2Result->GetYaxis()->SetTitleOffset(1.7);
    //chi2Result->SetMarkerStyle(marker[i]);
    chi2Result->SetLineColor(i+1);
    chi2Result->SetMarkerColor(i+1);
    chi2Result->SetStats(kFALSE);

    bayesResult->SetTitle(Form("b) Bayesian unfolding;%s;MC / unfolded", GetMultLabel()));
    bayesResult->GetXaxis()->SetRangeUser(0, displayRange);
    bayesResult->GetYaxis()->SetRangeUser(0.7, 1.3);
    bayesResult->GetYaxis()->SetTitleOffset(1.7);
    bayesResult->SetStats(kFALSE);
    //bayesResult->SetLineColor(2);
    bayesResult->SetLineColor(i+1);

    canvas->cd(1);
    gPad->SetRightMargin(0.05);
    gPad->SetLeftMargin(0.12);
    gPad->SetGridx();
    gPad->SetGridy();
    chi2Result->DrawCopy((i == 0) ? "HIST" : "HIST SAME");

    canvas->cd(2);
    gPad->SetRightMargin(0.05);
    gPad->SetLeftMargin(0.12);
    gPad->SetGridx();
    gPad->SetGridy();
    bayesResult->DrawCopy((i == 0) ? "HIST" : "HIST SAME");

    //TLine* line = new TLine(0, 1, 150, 1);
    //line->Draw();

    legend->AddEntry(chi2Result, names[i]);
  }

  canvas->cd(1);
  legend->Draw();
  canvas->cd(2);
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

void Draw2Unfolded(const char* file1, const char* file2, const char* output)
{
  loadlibs();
  
  TFile::Open(file1);
  AliMultiplicityCorrection* mult = new AliMultiplicityCorrection("Multiplicity", "Multiplicity");
  mult->LoadHistograms("Multiplicity");

  // result with systematic effect
  TFile::Open(file2);
  AliMultiplicityCorrection* mult2 = new AliMultiplicityCorrection("Multiplicity2", "Multiplicity2");
  mult2->LoadHistograms("Multiplicity");

  TH1* mcHist = mult2->GetMultiplicityVtx(etaRange)->ProjectionY("mymc", 1, 1);
  TH1* result1 = (TH1*) mult2->GetMultiplicityESDCorrected(etaRange)->Clone("result1"); // from file2 (with syst)
  TH1* result2 = (TH1*) mult->GetMultiplicityESDCorrected(etaRange)->Clone("result2");  // from file1 (without syst)

  DrawResultRatio(mcHist, result1, "tmp1.eps");
  DrawResultRatio(mcHist, result2, "tmp2.eps");
  Draw2ResultRatio(mcHist, result1, result2, output);
}

void PythiaPhojet()
{
  loadlibs();
  
  displayRange = 55;
  Draw2Unfolded("self.root", "pythia.root", "test.eps");
  
  canvas = (TCanvas*) gROOT->GetListOfCanvases()->Last();
  pad1 = (TCanvas*)canvas->GetListOfPrimitives()->First();
  pad2 = (TCanvas*)canvas->GetListOfPrimitives()->Last();
  legend = (TLegend*)pad1->GetListOfPrimitives()->Last();
  
  ((TH1*)pad2->GetListOfPrimitives()->At(1))->GetYaxis()->SetTitle("Ratio (Pythia / Phojet)");
  ((TLegendEntry*)legend->GetListOfPrimitives()->At(1))->SetLabel("Unfolded distribution (Pythia)");
  ((TLegendEntry*)legend->GetListOfPrimitives()->At(2))->SetLabel("Unfolded distribution (Phojet)");
  canvas->SaveAs("PythiaPhojet.eps");
}

void Misalignment()
{
  loadlibs();
  
  Draw2Unfolded("chi2_ideal.root", "chi2_misaligned.root", "test.eps");
  
  canvas = (TCanvas*) gROOT->GetListOfCanvases()->Last();
  pad1 = (TCanvas*)canvas->GetListOfPrimitives()->First();
  pad2 = (TCanvas*)canvas->GetListOfPrimitives()->Last();
  legend = (TLegend*)pad1->GetListOfPrimitives()->Last();

  ((TH1*)pad2->GetListOfPrimitives()->At(1))->GetYaxis()->SetTitle("Ratio (misaligned / realigned)");
  ((TLegendEntry*)legend->GetListOfPrimitives()->At(1))->SetLabel("Unfolded distribution (misaligned)");
  ((TLegendEntry*)legend->GetListOfPrimitives()->At(2))->SetLabel("Unfolded distribution (realigned)");
  canvas->SaveAs("SystematicMisalignment.eps");
}

void SystematicLowEfficiency()
{
  Draw2Unfolded("chi2.root", "chi2_loweff_5.root", "SystematicLowEfficiency.eps");
}

void SystematicLowEfficiency2()
{
  loadlibs();
  
  AliMultiplicityCorrection* mult = AliMultiplicityCorrection::Open("chi2.root");
  TH1* result2 = (TH1*) mult->GetMultiplicityESDCorrected(etaRange)->Clone("result2");
  TH1* mcHist = mult->GetMultiplicityVtx(etaRange)->ProjectionY("mymc", 1, 1);
  result2->Scale(1.0 / result2->Integral(2, displayRange));
  result2->Scale(mcHist->Integral(2, displayRange));
  
  mult = AliMultiplicityCorrection::Open("chi2_loweff_5.root");
  TH1* result1 = (TH1*) mult->GetMultiplicityESDCorrected(etaRange)->Clone("result1");
  
  mult = AliMultiplicityCorrection::Open("chi2_loweff_2.root");
  TH1* ratio2 = (TH1*) mult->GetMultiplicityESDCorrected(etaRange)->Clone("ratio2");
  ratio2->Scale(1.0 / ratio2->Integral(2, displayRange));
  ratio2->Scale(mcHist->Integral(2, displayRange));
  ratio2->Divide(result2);
  
  mult = AliMultiplicityCorrection::Open("chi2_loweff_1.root");
  TH1* ratio3 = (TH1*) mult->GetMultiplicityESDCorrected(etaRange)->Clone("ratio3");
  ratio3->Scale(1.0 / ratio3->Integral(2, displayRange));
  ratio3->Scale(mcHist->Integral(2, displayRange));
  ratio3->Divide(result2);
  
  Draw2ResultRatios(mcHist, result1, result2, ratio2, ratio3, "SystematicLowEfficiency2.eps");
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

void LowMomentumEffectSPD()
{
  // this function increases/reduces the correction as function of pt between 0 and 0.2 gev by +-50% to 0% and checks the effect on the overall correction factor
  // only a normal acceptance region is considered to not get a bias by the edges
  
  loadlibs();
  TFile::Open("multiplicity.root");
  
    
  AliCorrection* correction[8];
  Float_t values[3];
  
  for (Int_t loop=0; loop<3; loop++)
  {
    Float_t sumGen = 0;
    Float_t sumMeas = 0;
    
    Printf("loop %d", loop);
    for (Int_t i=0; i<4; ++i)
    {
      Printf("correction %d", i);
  
      TString name; name.Form("correction_%d", i);
      correction[i] = new AliCorrection(name, name);
      correction[i]->LoadHistograms();
      
      TH3* gene = correction[i]->GetTrackCorrection()->GetGeneratedHistogram();
      TH3* meas = correction[i]->GetTrackCorrection()->GetMeasuredHistogram();
  
      Float_t vtxRange = 5.9;
      gene->GetXaxis()->SetRangeUser(-vtxRange, vtxRange);
      meas->GetXaxis()->SetRangeUser(-vtxRange, vtxRange);
      
      Float_t etaRange = 0.99;
      gene->GetYaxis()->SetRangeUser(-etaRange, etaRange);
      meas->GetYaxis()->SetRangeUser(-etaRange, etaRange);
  
      TH1* genePt = gene->Project3D(Form("z_%d", i));
      TH1* measPt = meas->Project3D(Form("z_%d", i));
  
      if (loop > 0)
      {
        for (Int_t x=1; x<=genePt->GetNbinsX(); x++)
        {
          Float_t pt = genePt->GetXaxis()->GetBinCenter(x);
          //Printf("%f", pt);
          if (pt < 0.2)
          {
            Float_t factor = 1;
            if (loop == 1)
              factor = 1.5 - pt / 0.2 * 0.5;
            if (loop == 2)
              factor = 0.5 + pt / 0.2 * 0.5;
            //Printf("%f", factor);
            genePt->SetBinContent(x, genePt->GetBinContent(x) * factor);
            measPt->SetBinContent(x, measPt->GetBinContent(x) * factor);
          }
        }
      }
      
      //new TCanvas; genePt->DrawCopy(); measPt->DrawCopy("SAME");
  
      sumGen += genePt->Integral();
      sumMeas += measPt->Integral();  
      
      Float_t average = measPt->Integral() / genePt->Integral();
      
      Printf("The average efficiency of this correction is %f", average);
    }
    
    Float_t average = sumMeas / sumGen;
      
    Printf("The average efficiency of all corrections is %f", average);
    values[loop] = average;
  }
  
  Printf("relative is %f and %f", values[1] / values[0], values[2] / values[0]);
}
  

void EfficiencySpecies(Bool_t addDecayStopped = kFALSE)
{
  loadlibs();

  Int_t marker[] = {24, 25, 26, 27};
  Int_t color[] = {1, 2, 4, 3};

  // SPD TPC
  //const char* fileName[] = { "multiplicityMC_400k_syst.root", "multiplicityMC_TPC_4kfiles_syst.root" };
  //const char* fileName[] = { "LHC09b9_0.9TeV_0T/mb1/spd/multiplicity.root", "LHC09b8_0.9TeV_0.5T/mb1/tpc/multiplicity.root" };
  //const char* fileName[] = { "spd/multiplicity.root", "tpc/multiplicity.root" };
  const char* fileName[] = { "multiplicityparticle-efficiency.root", "multiplicity.root" };
  Float_t etaRangeArr[] = {0.49, 0.9};
  const char* titles[] = { "SPD Tracklets", "TPC Tracks" };

  TCanvas* canvas = new TCanvas("EfficiencySpecies", "EfficiencySpecies", 1000, 500);
  canvas->Divide(2, 1);

  TCanvas* canvas3 = new TCanvas("EfficiencySpecies_comb", "EfficiencySpecies_comb", 600, 600);
  gPad->SetGridx();
  gPad->SetGridy();
  gPad->SetRightMargin(0.05);
  gPad->SetTopMargin(0.05);
  
  TLegend* legends[2];
  
  for (Int_t loop=0; loop<1; ++loop)
  {
    Printf("%s", fileName[loop]);

    TCanvas* canvas2 = new TCanvas(Form("EfficiencySpecies_%d", loop), Form("EfficiencySpecies_%d", loop), 600, 600);
    gPad->SetGridx();
    gPad->SetGridy();
    gPad->SetRightMargin(0.05);
    gPad->SetTopMargin(0.05);
    
    AliCorrection* correction[8];

    canvas->cd(loop+1);

    gPad->SetGridx();
    gPad->SetGridy();
    gPad->SetRightMargin(0.05);

    TLegend* legend = new TLegend(0.6, 0.4, 0.85, 0.6);
    legend->SetFillColor(0);
    legend->SetEntrySeparation(0.2);
    legend->SetTextSize(gStyle->GetTextSize());
    
    legends[loop] = new TLegend(0.4+loop*0.3, 0.2, 0.6+loop*0.3, 0.5);
    legends[loop]->SetFillColor(0);
    legends[loop]->SetEntrySeparation(0.2);
    legends[loop]->SetTextSize(gStyle->GetTextSize());
    legends[loop]->SetHeader((loop == 0) ? "SPD" : "TPC");

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

    Float_t sumGenAbove02 = 0;
    Float_t sumMeasAbove02 = 0;
    
    for (Int_t i=0; i<3; ++i)
    {
      Printf("correction %d", i);

      TString name; name.Form("correction_%d", i);
      correction[i] = new AliCorrection(name, name);
      correction[i]->LoadHistograms();
      
      TH3* gene = correction[i]->GetTrackCorrection()->GetGeneratedHistogram();
      TH3* meas = correction[i]->GetTrackCorrection()->GetMeasuredHistogram();

      // limit vtx axis
      Float_t vtxRange = 3.9;
      gene->GetXaxis()->SetRangeUser(-vtxRange, vtxRange);
      meas->GetXaxis()->SetRangeUser(-vtxRange, vtxRange);

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
      Float_t etaBegin = -etaRangeArr[loop];
      Float_t etaEnd   = etaRangeArr[loop];
      //etaBegin = 0.01;
      //etaEnd = -0.01;
      gene->GetYaxis()->SetRangeUser(etaBegin, etaEnd);
      meas->GetYaxis()->SetRangeUser(etaBegin, etaEnd);

      TH1* genePt = gene->Project3D(Form("z_%d", i));
      TH1* measPt = meas->Project3D(Form("z_%d", i));
      
//       if (i == 2)
//       {
//         Printf("WARNING: Rebinning for protons!");
//         genePt->Rebin(2);
//         measPt->Rebin(2);
//       }

      genePt->Sumw2();
      measPt->Sumw2();
      
      for (Int_t x=0; x<=genePt->GetNbinsX()+1; x++)
      {
      	genePt->SetBinError(x, TMath::Sqrt(genePt->GetBinContent(x)));
      	measPt->SetBinError(x, TMath::Sqrt(measPt->GetBinContent(x)));
      }
      
      sumGen += genePt->Integral();
      sumMeas += measPt->Integral();

      sumGenAbove02 += genePt->Integral(genePt->GetXaxis()->FindBin(0.21), genePt->GetNbinsX());
      sumMeasAbove02 += measPt->Integral(genePt->GetXaxis()->FindBin(0.21), genePt->GetNbinsX());
      
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
      effPt->SetMarkerSize(2);

      effPt->GetXaxis()->SetRangeUser(0, 1);
      effPt->GetYaxis()->SetRangeUser(0.001, 1);

      effPt->GetXaxis()->SetTitleOffset(1.1);
      effPt->GetYaxis()->SetTitleOffset(1.2);
      
      effPt->GetXaxis()->SetTitle("p_{T} (GeV/c)");

      effPt->SetStats(kFALSE);
      effPt->SetTitle(titles[loop]);
      effPt->GetYaxis()->SetTitle("Efficiency");

      canvas->cd(loop+1);
      effPt->DrawCopy((i == 0) ? "" : "SAME");
 
      canvas2->cd();
      effPt->SetTitle("");
      effPt->DrawCopy((i == 0) ? "" : "SAME");
      
      canvas3->cd();
      effPtClone = (TH1*) effPt->Clone("effPtClone");
      effPtClone->SetMarkerStyle(marker[i]-4*loop);
      effPtClone->DrawCopy((i == 0 && loop == 0) ? "" : "SAME");

      legend->AddEntry(effPt, ((i == 0) ? "#pi^{#pm}" : ((i == 1) ? "K^{#pm}" : "p,#bar{p}")), "P");
      legends[loop]->AddEntry(effPtClone, ((i == 0) ? "#pi^{#pm}" : ((i == 1) ? "K^{#pm}" : "p,#bar{p}")), "P");
      //legend2->AddEntry(effPt, Form("%s %s", (loop == 0) ? "SPD" : "TPC", ((i == 0) ? "#pi^{#pm}" : ((i == 1) ? "K^{#pm}" : "p,#bar{p}"))), "P");

      if (addDecayStopped)
      {
        name.Form("correction_%d", i+4);
        corr = new AliCorrection(name, name);
        corr->LoadHistograms();
        
        TH3* gene = corr->GetTrackCorrection()->GetGeneratedHistogram();
        TH3* meas = corr->GetTrackCorrection()->GetMeasuredHistogram();
        
        // limit axes
        gene->GetXaxis()->SetRangeUser(-vtxRange, vtxRange);
        meas->GetXaxis()->SetRangeUser(-vtxRange, vtxRange);
        gene->GetYaxis()->SetRangeUser(-etaRangeArr[loop], etaRangeArr[loop]);
        meas->GetYaxis()->SetRangeUser(-etaRangeArr[loop], etaRangeArr[loop]);
        
        TH1* decayed = gene->Project3D(Form("z_%d", i+4));
        TH1* stopped = meas->Project3D(Form("z_%d", i+4));
        
        Printf("%d: %d decayed, %d stopped, out of %d", i, (Int_t) decayed->Integral(), (Int_t) stopped->Integral(), (Int_t) genePt->Integral());
        
        decayed->Divide(decayed, genePt, 1, 1, "B");
        stopped->Divide(stopped, genePt, 1, 1, "B");
        
        decayed->SetMarkerStyle(20);
        stopped->SetMarkerStyle(21);
        stopped->SetMarkerColor(2);
        
        new TCanvas(Form("all_%d_%d", loop, i), Form("all_%d_%d", loop, i), 600, 600);
        effPt->DrawCopy();
        decayed->DrawCopy("SAME");
        stopped->DrawCopy("SAME");
        
        decayed->Add(stopped);
        decayed->Add(effPt);
        decayed->SetMarkerStyle(22);
        decayed->SetMarkerColor(4);
        decayed->DrawCopy("SAME");
      }
      
    }

    Printf("In total %.4f of the particles are below their effective pt cut off", (Float_t) below / total);

    Printf("%f measured, %f generated, effiency: %f", sumGen, sumMeas, sumMeas / sumGen);
    Printf("Above 0.2 GeV/c: %f measured, %f generated, effiency: %f", sumGenAbove02, sumMeasAbove02, sumMeasAbove02 / sumGenAbove02);

    canvas->cd(loop+1);
    legend->Draw();
  
    canvas2->cd();
    legend->Draw();
    canvas2->SaveAs(Form("%s.eps", canvas2->GetName()));  
  }
  
  return;

  canvas->cd();
  canvas->SaveAs(Form("%s.eps", canvas->GetName()));

  canvas3->cd();
  legends[0]->Draw();
  legends[1]->Draw();
  canvas3->SaveAs(Form("%s.eps", canvas3->GetName()));
}

void DrawpTCutOff()
{
/*
aliroot -b -q $ALICE_ROOT/PWG0/multiplicity/correct.C'("multiplicityMC_pt_ref.root", "Multiplicity", "multiplicityESD.root", kTRUE)'
mv unfolded.root chi2_ptref.root
aliroot -b -q $ALICE_ROOT/PWG0/multiplicity/correct.C'("multiplicityMC_pt0.root", "Multiplicity", "multiplicityESD.root", kTRUE)'
mv unfolded.root chi2_pt0.root
aliroot -b -q $ALICE_ROOT/PWG0/multiplicity/correct.C'("multiplicityMC_pt1.root", "Multiplicity", "multiplicityESD.root", kTRUE)'
mv unfolded.root chi2_pt1.root
aliroot -b -q $ALICE_ROOT/PWG0/multiplicity/correct.C'("multiplicityMC_pt0_25.root", "Multiplicity", "multiplicityESD.root", kTRUE)'
mv unfolded.root chi2_pt0_25.root
aliroot -b -q $ALICE_ROOT/PWG0/multiplicity/correct.C'("multiplicityMC_pt1_25.root", "Multiplicity", "multiplicityESD.root", kTRUE)'
mv unfolded.root chi2_pt1_25.root
*/

  loadlibs();
  
  TH1* results[10];
  const char* files[] = { "chi2_ptref.root", "chi2_pt0.root", "chi2_pt1.root", "chi2_pt0_25.root", "chi2_pt1_25.root"};
  
  Int_t nMax = 5;
  for (Int_t i = 0; i<nMax; ++i)
  {
    AliMultiplicityCorrection* mult = AliMultiplicityCorrection::Open(files[i]);
    results[i] = (TH1*) mult->GetMultiplicityESDCorrected(etaRange)->Clone(Form("result_%d", i));
  }
  
  const char* legendStrings[] = { "Reduced 50%", "Enhanced 50%", "Reduced 25%", "Enhanced 25%" };
  DrawRatio(results[0], nMax-1, results+1, "LowMomentumSyst.eps", kFALSE, legendStrings);
}

void ParticleSpeciesComparison()
{
  loadlibs();

  TH1* results[10];
  TH1* mc = 0;
  
  // loop over cases (normal, enhanced/reduced ratios)
  Int_t nMax = 7;
  for (Int_t i = 0; i<nMax; ++i)
  {
    AliMultiplicityCorrection* mult = AliMultiplicityCorrection::Open(Form("chi2a_inel_species_%d.root", i), "Multiplicity");
    results[i] = (TH1*) mult->GetMultiplicityESDCorrected(etaRange)->Clone(Form("result_%d", i));
  }

  for (Int_t i=1; i<=results[0]->GetNbinsX(); i++)
    results[0]->SetBinError(i, 0);

  const char* legendStrings[] = { "K #times 0.7", "K #times 1.3", "p #times 0.7", "p #times 1.3", "others #times 0.7", "others #times 1.3",  };

  DrawRatio(results[0], nMax-1, results+1, "ParticleSpeciesComparison1_2.eps", kFALSE, legendStrings);
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

void CompareVertexRecoEfficiencies()
{
  loadlibs();
  
  const char* file1 = "LHC10a12_run10482X";
  const char* file2 = "LHC10a14_run10482X_Phojet";
  
  const char* suffix = "/all/spd/multiplicityMC_xsection.root";
  
  hist1 = AliMultiplicityCorrection::Open(Form("%s%s", file1, suffix))->GetEfficiency(etaRange, AliMultiplicityCorrection::kMB);
  hist2 = AliMultiplicityCorrection::Open(Form("%s%s", file2, suffix))->GetEfficiency(etaRange, AliMultiplicityCorrection::kMB);
  
  ratio = (TH1*) hist1->Clone("ratio");
  ratio->Divide(hist2);
  
  new TCanvas;
  hist1->Draw();
  hist2->SetLineColor(2);
  hist2->Draw("SAME");
  
  new TCanvas;
  ratio->Draw();
}

void TriggerVertexCorrection()
{
  //
  // plots the correction performed on the unfolded spectrum to gain the spectrum for the full inelastic sample
  //

  loadlibs();

  TFile::Open(correctionFile);
  AliMultiplicityCorrection* mult = new AliMultiplicityCorrection("Multiplicity", "Multiplicity");
  mult->LoadHistograms("Multiplicity");

  TH1* corrINEL = Invert(mult->GetEfficiency(etaRange, AliMultiplicityCorrection::kINEL));
  TH1* corrNSD = Invert(mult->GetEfficiency(etaRange, AliMultiplicityCorrection::kNSD));
  TH1* corrMB   = Invert(mult->GetEfficiency(etaRange, AliMultiplicityCorrection::kMB));
  
  TH1* triggerEff = Invert(mult->GetTriggerEfficiency(etaRange, AliMultiplicityCorrection::kNSD));

  TCanvas* canvas = new TCanvas("TriggerVertexCorrection", "TriggerVertexCorrection", 800, 500);
  gPad->SetGridx();
  gPad->SetGridy();
  gPad->SetTopMargin(0.05);
  gPad->SetRightMargin(0.05);

  corrINEL->SetStats(kFALSE);
  corrINEL->GetXaxis()->SetRangeUser(0, 12);
  corrINEL->GetYaxis()->SetRangeUser(0, 8);
  corrINEL->SetTitle(Form(";%s;Correction factor", GetMultLabel()));
  corrINEL->SetMarkerStyle(22);
  corrINEL->Draw("PE");

  corrMB->SetStats(kFALSE);
  corrMB->SetLineColor(2);
  corrMB->SetMarkerStyle(25);
  corrMB->SetMarkerColor(2);
  corrMB->Draw("SAME PE");
  
  corrNSD->SetLineColor(4);
  corrNSD->SetMarkerStyle(24);
  corrNSD->SetMarkerColor(4);
  corrNSD->Draw("SAME PE");
  
  triggerEff->SetLineColor(6);
  triggerEff->SetMarkerStyle(25);
  triggerEff->SetMarkerColor(6);
  //triggerEff->Draw("SAME PE");
  
  Printf("       MB  INEL  NSD  TRIGINEL");
  Printf("bin 0: %.2f %.2f %.2f %.2f", corrMB->GetBinContent(1), corrINEL->GetBinContent(1), corrNSD->GetBinContent(1), triggerEff->GetBinContent(1));
  Printf("bin 1: %.2f %.2f %.2f %.2f", corrMB->GetBinContent(2), corrINEL->GetBinContent(2), corrNSD->GetBinContent(2), triggerEff->GetBinContent(2));

  TLegend* legend = new TLegend(0.3, 0.6, 0.85, 0.85);
  legend->SetFillColor(0);
  legend->AddEntry(corrINEL, "Correction to inelastic sample");
  legend->AddEntry(corrNSD, "Correction to NSD sample");
  legend->AddEntry(corrMB, "Correction to triggered sample");
  legend->SetTextSize(0.04);

  legend->Draw();

  canvas->SaveAs(Form("%s.eps", canvas->GetName()));
}

void DrawStatisticalUncertainty(const char* fileName = "StatisticalUncertainty.root")
{
  TFile::Open(fileName);
  
  errorResponse = (TH1*) gFile->Get("errorResponse");
  errorMeasured = (TH1*) gFile->Get("errorMeasured");
  errorBoth = (TH1*) gFile->Get("errorBoth");
  
  TCanvas* canvas = new TCanvas("StatisticalUncertainty", "StatisticalUncertainty", 600, 400);
  canvas->SetGridx();
  canvas->SetGridy();
  canvas->SetRightMargin(0.05);
  canvas->SetTopMargin(0.05);

  errorResponse->SetLineColor(1);
  errorResponse->GetXaxis()->SetRangeUser(0, longDisplayRange);
  errorResponse->GetYaxis()->SetRangeUser(0, 1);
  errorResponse->SetStats(kFALSE);
  errorResponse->SetTitle(";true multiplicity;Uncertainty");

  errorResponse->Draw();

  errorMeasured->SetLineColor(2);
  errorMeasured->Draw("SAME");

  errorBoth->SetLineColor(4);
  errorBoth->Draw("SAME");

  errorResponse->Scale(1.0 / sqrt(2));
  errorMeasured->Scale(1.0 / sqrt(2));
  errorBoth->Scale(1.0 / sqrt(2));
  
  Printf("Average errorResponse: %f", errorResponse->Integral(2, displayRange) / (displayRange - 1));
  Printf("Average errorMeasured: %f", errorMeasured->Integral(2, displayRange) /  (displayRange - 1));
  Printf("Average errorBoth: %f", errorBoth->Integral(2, displayRange) /  (displayRange - 1));

  canvas->SaveAs(Form("%s.eps", canvas->GetName()));
}

void StatisticalUncertaintyCompare(const char* det = "SPD")
{
  TFile* file1 = TFile::Open(Form("StatisticalUncertainty%sBayesian.root", det));
  TH1* errorResponse = (TH1*) file1->Get("errorResponse");
  TH1* errorMeasured = (TH1*) file1->Get("errorMeasured");
  TH1* errorBoth = (TH1*) file1->Get("errorBoth");

  TString str;
  str.Form("StatisticalUncertaintyCompare%s", det);

  TCanvas* canvas = new TCanvas(str, str, 800, 500);
  canvas->SetGridx();
  canvas->SetGridy();
  canvas->SetRightMargin(0.05);
  canvas->SetTopMargin(0.05);
  
  errorResponse->Scale(1.0 / sqrt(2));
  errorMeasured->Scale(1.0 / sqrt(2));
  errorBoth->Scale(1.0 / sqrt(2));

  errorResponse->SetLineColor(1);
  errorResponse->GetXaxis()->SetRangeUser(0, displayRange);
  errorResponse->GetYaxis()->SetRangeUser(0, 0.18);
  errorResponse->SetStats(kFALSE);
  errorResponse->GetYaxis()->SetTitleOffset(1.2);
  errorResponse->SetTitle(Form(";%s;#sqrt{2}^{-1} #sigma(unfolded - unfolded_{0}) / unfolded_{0}", GetMultLabel()));

  errorResponse->Draw();

  errorMeasured->SetLineColor(2);
  errorMeasured->Draw("SAME");

  errorBoth->SetLineColor(4);
  errorBoth->Draw("SAME");

  TFile* file2 = TFile::Open(Form("StatisticalUncertainty%sChi2.root", det));
  TH1* errorResponse2 = (TH1*) file2->Get("errorResponse");
  TH1* errorMeasured2 = (TH1*) file2->Get("errorMeasured");
  TH1* errorBoth2 = (TH1*) file2->Get("errorBoth");

  errorResponse2->Scale(1.0 / sqrt(2));
  errorMeasured2->Scale(1.0 / sqrt(2));
  errorBoth2->Scale(1.0 / sqrt(2));
  
  errorResponse2->SetLineStyle(2);
  errorResponse2->Draw("SAME");
  
  errorMeasured2->SetLineColor(2);
  errorMeasured2->SetLineStyle(2);
  errorMeasured2->Draw("SAME");
  
  errorBoth2->SetLineColor(4);
  errorBoth2->SetLineStyle(2);
  errorBoth2->Draw("SAME");

  TLegend* legend = new TLegend(0.2, 0.5, 0.8, 0.9);
  legend->SetFillColor(0);
  legend->SetTextSize(0.04);
  legend->AddEntry(errorBoth, "Both (Bayesian unfolding)");
  legend->AddEntry(errorMeasured, "Measured (Bayesian unfolding)");
  legend->AddEntry(errorResponse, "Response matrix (Bayesian unfolding)");
  legend->AddEntry(errorBoth2, "Both (#chi^{2}-minimization)");
  legend->AddEntry(errorMeasured2, "Measured (#chi^{2}-minimization)");
  legend->AddEntry(errorResponse2, "Response matrix (#chi^{2}-minimization)");
  legend->Draw();

  canvas->SaveAs(Form("%s.eps", canvas->GetName()));
}

void EfficiencyComparison(Int_t eventType = 2, Bool_t uncertainty = kTRUE)
{
 const char* files[] = { "multiplicityMC_nd.root", "multiplicityMC_sd.root", "multiplicityMC_dd.root", "multiplicityMC_xsection.root" };

  loadlibs();

  TCanvas* canvas = new TCanvas("EfficiencyComparison", "EfficiencyComparison", 800, 600);
  canvas->SetGridx();
  canvas->SetGridy();
  canvas->SetRightMargin(0.05);
  canvas->SetTopMargin(0.05);

  AliMultiplicityCorrection* data[4];
  TH1* effArray[4];
  TH1* effErrorArray[2];

  Int_t markers[] = { 24, 25, 26, 5 };
  //Int_t markers[] = { 2, 25, 24, 5 };
  Int_t colors[] = { 1, 2, 4, 6 };
  //Int_t colors[] = { 1, 1, 1, 1 };

  //TLegend* legend = new TLegend(0.45, 0.45, 0.9, 0.7);
  TLegend* legend = new TLegend(0.3, 0.3, 0.9, 0.6);
  legend->SetTextSize(0.04);
  legend->SetFillColor(0);
 
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

    TH1* eff = 0;
    if (eventType == -1)
    {
      eff = (TH1*) data[i]->GetTriggerEfficiency(etaRange)->Clone(Form("eff_%d", i));
    }
    else
      eff = (TH1*) data[i]->GetEfficiency(etaRange, (AliMultiplicityCorrection::EventType) eventType)->Clone(Form("eff_%d", i));
    effArray[i] = eff;

    eff->GetXaxis()->SetRangeUser(0, 15);
    eff->GetYaxis()->SetRangeUser(0, 1.19);
    eff->SetStats(kFALSE);
    eff->GetXaxis()->SetTitle(GetMultLabel());
    eff->GetYaxis()->SetTitle("Efficiency");
    eff->SetTitle("");
    eff->SetLineColor(colors[i]);
    eff->SetMarkerColor(colors[i]);
    eff->SetMarkerStyle(markers[i]);

    if (i == 3)
    {
      // once for INEL, once for NSD
      for (AliMultiplicityCorrection::EventType eventType2 = AliMultiplicityCorrection::kINEL; eventType2 <= AliMultiplicityCorrection::kNSD; eventType2++)
      {
        effDiff = (TH1*) data[i]->GetEfficiency(etaRange, eventType2)->Clone(Form("effDiff_%d", i));
        
        for (Int_t bin=1; bin<=effDiff->GetNbinsX(); bin++)
          effDiff->SetBinError(bin, 0);
  
        // loop over cross section combinations
        for (Int_t j=1; j<7; ++j)
        {
          AliMultiplicityCorrection* mult = new AliMultiplicityCorrection("Multtmp", "Multtmp");
          mult->LoadHistograms(Form("Multiplicity_%d", j));
  
          TH1* eff2 = mult->GetEfficiency(etaRange, eventType2);
  
          for (Int_t bin=1; bin<=effDiff->GetNbinsX(); bin++)
          {
            // TODO we could also do asymmetric errors here
            Float_t deviation = TMath::Abs(effDiff->GetBinContent(bin) - eff2->GetBinContent(bin));
  
            effDiff->SetBinError(bin, TMath::Max(effDiff->GetBinError(bin), (Double_t) deviation));
          }
        }
  
        for (Int_t bin=1; bin<=effDiff->GetNbinsX(); bin++)
        {
          //if (eventType2 == AliMultiplicityCorrection::kINEL)
            //eff->SetBinError(bin, 0);
            //eff->SetBinError(bin, effDiff->GetBinError(bin));
          if (bin < 20 && effDiff->GetBinContent(bin) > 0)
            Printf("Bin %d: Error: %.2f", bin, 100.0 * effDiff->GetBinError(bin) / effDiff->GetBinContent(bin));
        }
        
        if (uncertainty) {
	        TH1* effError = (TH1*) effDiff->Clone(Form("effError_%s", (eventType2 == AliMultiplicityCorrection::kINEL) ? "INEL" : "NSD"));
          effErrorArray[eventType2 - AliMultiplicityCorrection::kINEL] = effError;
	        effError->Reset();
        
	        for (Int_t bin=1; bin<=effDiff->GetNbinsX(); bin++)
	          if (effDiff->GetBinContent(bin) > 0)
	            effError->SetBinContent(bin, 1.0 * effDiff->GetBinError(bin) / effDiff->GetBinContent(bin));
        
	        effError->SetLineColor(1);
	        effError->SetMarkerStyle(1);
	        
          if (eventType2 == AliMultiplicityCorrection::kNSD)
            effError->SetLineStyle(2);
          
          effError->DrawCopy("SAME HIST");
        }
      }
    }

    canvas->cd();
    if (i == 0)
    {
      eff->DrawCopy("P");
    }
    else
      eff->DrawCopy("SAME P");

    legend->AddEntry(eff, (((i == 0) ? "Non-diffractive" : ((i == 1) ? "Single-diffractive" : ((i == 2) ? "Double-diffractive" : "Pythia combined")))));
  }

  if (uncertainty)
  {
    legend->AddEntry(effErrorArray[0], "Relative syst. uncertainty: inelastic");
    legend->AddEntry(effErrorArray[1], "Relative syst. uncertainty: NSD");
  
    file = TFile::Open("uncertainty_xsection.root", "RECREATE");
    effErrorArray[0]->Write();
    effErrorArray[1]->Write();
    file->Close();
  }

  legend->Draw();

  canvas->SaveAs(Form("%s.eps", canvas->GetName()));
}

void ModelDependencyPlot()
{
  loadlibs();

  TFile::Open("multiplicityMC.root");
  AliMultiplicityCorrection* mult = new AliMultiplicityCorrection("Multiplicity", "Multiplicity");
  mult->LoadHistograms("Multiplicity");
  
  hist = mult->GetCorrelation(etaRange);
  
  for (Int_t y=0; y<=hist->GetYaxis()->GetNbins()+1; ++y)
  {
    for (Int_t z=0; z<=hist->GetZaxis()->GetNbins()+1; ++z)
    {
      hist->SetBinContent(0, y, z, 0);
      hist->SetBinContent(hist->GetXaxis()->GetNbins()+1, y, z, 0);
    }
  }

  TH2* proj = (TH2*) hist->Project3D("zy");

  TCanvas* canvas = new TCanvas("ModelDependencyPlot", "ModelDependencyPlot", 1200, 600);

  canvas->Divide(2, 1);

  canvas->cd(2);
  gPad->SetLogy();
  gPad->SetGridx();
  gPad->SetGridy();
  gPad->SetTopMargin(0.05);
  gPad->SetRightMargin(0.05);
 
  Int_t selectedMult = 30;
  Int_t yMax = 9e4;

  TH1* full = proj->ProjectionX("full");
  TH1* selected = proj->ProjectionY("selected", proj->GetXaxis()->FindBin(selectedMult), proj->GetXaxis()->FindBin(selectedMult)); 

  full->SetStats(kFALSE);
  full->GetXaxis()->SetRangeUser(0, displayRange);
  full->GetYaxis()->SetRangeUser(5, yMax);
  full->SetTitle(";Multiplicity;Entries");

  selected->SetLineColor(0);
  selected->SetMarkerColor(2);
  selected->SetMarkerStyle(5);

  full->Draw();
  selected->Draw("SAME P");

  TLegend* legend = new TLegend(0.5, 0.65, 0.85, 0.85);
  legend->SetFillColor(0);
  legend->AddEntry(full, "True");
  legend->AddEntry(selected, "Measured");
  legend->Draw();
 
  TLine* line = new TLine(selectedMult, 5, selectedMult, yMax);
  line->SetLineWidth(2);
  line->Draw();

  canvas->cd(1);
  gPad->SetLogy();
  gPad->SetGridx();
  gPad->SetGridy();
  gPad->SetTopMargin(0.05);
  gPad->SetRightMargin(0.05);

  full = proj->ProjectionY("full2");
  selected = proj->ProjectionX("selected2", proj->GetYaxis()->FindBin(selectedMult), proj->GetYaxis()->FindBin(selectedMult));

  full->SetStats(kFALSE);
  full->GetXaxis()->SetRangeUser(0, displayRange);
  full->GetYaxis()->SetRangeUser(5, yMax);
  full->SetTitle(";Multiplicity;Entries");

  full->SetLineColor(0);
  full->SetMarkerColor(2);
  full->SetMarkerStyle(5);

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

void FitPt(const char* fileName)
{
  // needs a MC file from the dNdEta analysis

  TFile::Open(fileName);

  TH1* genePt = (TH1*) gFile->Get("fdNdpT");
  
  genePt->SetTitle(";p_{T} (GeV/c);1/p_{T} dN_{ch}/dp_{T} (GeV/c)^{-2}");
  // number of events
  genePt->Scale(1.0 / 287800);
  // bin width
  genePt->Scale(1.0 / genePt->GetXaxis()->GetBinWidth(1));
  
  genePt->GetXaxis()->SetRangeUser(0, 0.4);

  TF1* func = new TF1("func", "[1]*x*exp(x*[0])");
  func->SetParameters(-1, 1);

  genePt->SetMarkerStyle(25);
  genePt->SetTitle("");
  genePt->SetStats(kFALSE);

  new TCanvas;
  genePt->DrawCopy("P");
  func->DrawCopy("SAME");
  gPad->SetLogy();

  genePt->Fit(func, "0", "", 0, 0.25);
  genePt->Fit(func, "0", "", 0, 0.25);

  TCanvas* canvas = new TCanvas("FitPt", "FitPt", 600, 600);

  gPad->SetGridx();
  gPad->SetGridy();
  gPad->SetLeftMargin(0.13);
  gPad->SetRightMargin(0.05);
  gPad->SetTopMargin(0.05);

  //genePt->GetXaxis()->SetRangeUser(0, 1);
  genePt->GetYaxis()->SetRangeUser(2, 200);
  genePt->GetYaxis()->SetTitleOffset(1.4);
  genePt->GetXaxis()->SetTitleOffset(1.1);
  genePt->DrawCopy("P");
  //func->SetRange(0, 0.3);
  func->DrawCopy("SAME");
  gPad->SetLogy();

  canvas->SaveAs(Form("%s.eps", canvas->GetName()));
  
  TH1* first = (TH1*) func->GetHistogram()->Clone("first");

  TCanvas* canvas2 = new TCanvas("FitPt2", "FitPt2", 600, 400);

  TFile* file = TFile::Open("ptspectrum_fit.root", "RECREATE");

  for (Int_t param=0; param<2; param++)
  {
    for (Int_t sign=0; sign<2; sign++)
    {
      TF1* func2 = (TF1*) func->Clone(Form("func_%d_%d", param, sign));  // new TF1(Form("func_%d_%d", param, sign), &FitPtFunc, 0, 2, 6);
      func2->SetParameters(func->GetParameters());
      //TF1* func2 = (TF1*) func->Clone(); // SetParameter after this does not work

      Float_t factor = ((sign == 0) ? 0.75 : 1.25);
      func2->SetParameter(param, func2->GetParameter(param) * factor);
      //func2->Print();

      canvas->cd();
      func2->SetLineWidth(2);
      func2->SetLineColor(2);
      func2->SetLineStyle(2);
      func2->DrawCopy("SAME");

      canvas2->cd();
      TH1* second = func2->GetHistogram();
      second->Divide(first);
      second->SetLineColor(param + 1);
      // set to 1 above 0.2 GeV
      for (Int_t bin=second->GetXaxis()->FindBin(0.20001); bin<=second->GetNbinsX(); bin++)
        second->SetBinContent(bin, 1);
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

TH1* SystematicsSummary(Bool_t tpc = 0, Bool_t nsd = kTRUE)
{
  Int_t nEffects = 7;

  TH1* effects[10];
  const char** names = 0;
  Int_t colors[] = { 1, 2, 4, 1, 2, 4 };
  Int_t styles[] = { 1, 2, 3, 1, 2, 3 };
  Int_t widths[] = { 1, 1, 1, 2, 2, 2 };
  Int_t markers[] = { 20, 21, 22, 23, 24, 25, 26 };
  
  TH1* dummy = new TH2F("dummy", Form(";%s;Uncertainty", GetMultLabel()), 202, -1.5, 200.5, 100, 0, 0.4);
  dummy->SetStats(0);

  for (Int_t i=0; i<nEffects; ++i)
    effects[i] = new TH1F("SystematicsSummary", Form(";%s;Uncertainty", GetMultLabel()), 201, -0.5, 200.5);

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
    nEffects = 5;

    //const char* namesSPD[] = { "Particle composition",  "p_{t} cut-off", "Unfolding Method (#chi^{2})", "Relative cross-sections (INEL)", "Relative cross-sections (NSD)"};
    const char* namesSPD[] = { "Unfolding Method (#chi^{2})", "Relative cross-sections (INEL)", "Relative cross-sections (NSD)", "Particle composition",  "p_{t} cut-off"};
    names = namesSPD;

    currentEffect = 0;

    // method
    TFile* file = TFile::Open("StatisticalUncertaintySPDChi2.root");
    TH1* hist = (TH1*) file->Get("errorBoth");
    hist->Scale(1.0 / sqrt(2));

    // smooth a bit, but skip 0 bin
    /*effects[currentEffect]->SetBinContent(1, hist->GetBinContent(1));
    for (Int_t i=2; i<=201; ++i)
      effects[currentEffect]->SetBinContent(i, (hist->GetBinContent(i) + hist->GetBinContent(i+1)) / 2);*/
    effects[currentEffect] = hist;

    currentEffect++;

    // relative x-section
    file = TFile::Open("uncertainty_xsection.root");
    effects[currentEffect++] = (TH1*) file->Get("effError_INEL");
    effects[currentEffect] = (TH1*) file->Get("effError_NSD");
    effects[currentEffect]->SetLineStyle(1);
    //effects[2]->SetBinContent(1, 0.20);
    //effects[2]->SetBinContent(2, 0.01);
    //effects[2]->SetBinContent(3, 0.002);
    
    currentEffect++;
    
    // particle composition
    effects[currentEffect]->SetBinContent(1, 0.16);
    for (Int_t i=2; i<=81; ++i)
    {
      effects[currentEffect]->SetBinContent(i, 0.01 + 0.05 * i / 81);
    }
    
    currentEffect++;

    // pt spectrum
    effects[currentEffect]->SetBinContent(1, 0.06);
    effects[currentEffect]->SetBinContent(2, 0.03);
    for (Int_t i=3; i<=81; ++i)
    {
      if (i <= 61)
      {
        effects[currentEffect]->SetBinContent(i, 0.01);
      }
      else if (i <= 81)
      {
        effects[currentEffect]->SetBinContent(i, 0.01 + 0.05 * (i - 61) / 20);
      }
    }
    
//     currentEffect++;
//         
//     // material budget
//     for (Int_t i=1; i<=81; ++i)
//     {
//       if (i < 5)
//         effects[currentEffect]->SetBinContent(i, 0.05 - 0.01 * i);
//       if (i > 51)
//         effects[currentEffect]->SetBinContent(i, 0.05 * (i - 50) / 30);
//     }
//         
    currentEffect++;
    
  }

  TCanvas* canvas = new TCanvas("SystematicsSummary.eps", "SystematicsSummary.eps", 800, 500);
  canvas->SetRightMargin(0.05);
  canvas->SetTopMargin(0.05);
  //canvas->SetGridx();
  canvas->SetGridy();
  TLegend* legend = new TLegend(0.2, 0.4, 0.7, 0.4 + 0.5 * nEffects / 7);
  legend->SetFillColor(0);
  legend->SetTextSize(0.04);
  dummy->Draw();
  dummy->GetXaxis()->SetRangeUser(0, displayRange);

  for (Int_t i=0; i<nEffects; ++i)
  {
    TH1* current = (TH1*) effects[i]->Clone(Form("current_%d", i));
    /*current->Reset();
    for (Int_t j=0; j<nEffects-i; ++j)
      current->Add(effects[j]);*/

    current->SetLineColor(colors[i]);
    current->SetLineStyle(styles[i]);
    current->SetLineWidth(widths[i]);
    //current->SetFillColor(colors[i]);
    current->SetMarkerColor(colors[i]);
    //current->SetMarkerStyle(markers[i]);

    current->SetStats(kFALSE);
    current->GetYaxis()->SetRangeUser(0, 0.4);
    current->DrawCopy("SAME");
    legend->AddEntry(current, names[i]);

    //TLatex* text = new TLatex(displayRange+2, current->GetBinContent(displayRange+1), names[i]);
    TLatex* text = new TLatex(displayRange+2, 0.1 - i * 0.02, names[i]);
    text->SetTextSize(0.04);
    text->SetTextColor(colors[i]);
    //text->Draw();
  }

  // add total in square
  TH1* totalINEL = (TH1*) effects[0]->Clone("totalINEL");
  totalINEL->Reset();
  TH1* totalNSD = (TH1*) totalINEL->Clone("totalNSD");

  for (Int_t i=0; i<nEffects; ++i)
  {
    //Printf("%d %f", i, effects[i]->GetBinContent(20));
    effects[i]->Multiply(effects[i]);
    
    if (i != 2)
      totalINEL->Add(effects[i]);
    if (i != 1)
      totalNSD->Add(effects[i]);
  }
  
  for (Int_t i=1; i<=totalINEL->GetNbinsX(); ++i)
  {
    totalINEL->SetBinError(i, 0);
    if (totalINEL->GetBinContent(i) > 0)
      totalINEL->SetBinContent(i, TMath::Min(sqrt(totalINEL->GetBinContent(i)), 1.0));
    totalNSD->SetBinError(i, 0);
    if (totalNSD->GetBinContent(i) > 0)
      totalNSD->SetBinContent(i, TMath::Min(sqrt(totalNSD->GetBinContent(i)), 1.0));
  }

  //Printf("%f", total->GetBinContent(20));

  totalINEL->SetMarkerStyle(5);
  totalINEL->SetMarkerColor(1);
  legend->AddEntry(totalINEL, "Total (INEL)", "P");
  
  totalNSD->SetMarkerStyle(24);
  totalNSD->SetMarkerColor(2);
  legend->AddEntry(totalNSD, "Total (NSD)", "P");
  
  Printf("total in bin 0 is INEL: %f NSD: %f", totalINEL->GetBinContent(1), totalNSD->GetBinContent(1));
  totalINEL->DrawCopy("SAME P"); //->SetBinContent(1, 0);
  totalNSD->DrawCopy("SAME P"); //->SetBinContent(1, 0);

  legend->Draw();

  canvas->SaveAs(canvas->GetName());
  
  file = TFile::Open(Form("%s_syst_error.root", (tpc) ? "tpc" : "spd"), "RECREATE");
  totalINEL->Write("inel_1");
  totalNSD->Write("nsd_1");
  file->Close();

  return (nsd) ? totalNSD : totalINEL;
}

void finalPlot(Bool_t tpc = 0, Bool_t small = kFALSE)
{
  loadlibs();

  if (tpc)
    SetTPC();

  //TH1* errorNSD = SystematicsSummary(tpc, 1);

  TCanvas* canvas = new TCanvas("finalPlot.eps", "finalPlot.eps", 800, 500); 
  canvas->SetRightMargin(0.05);
  canvas->SetTopMargin(0.05);
  canvas->SetGridx();
  canvas->SetGridy();
  
  legend = new TLegend(0.5, 0.6, 0.9, 0.8);
  legend->SetFillColor(0);
  legend->SetTextSize(0.04);
  
  for (AliMultiplicityCorrection::EventType eventType = AliMultiplicityCorrection::kINEL; eventType <= AliMultiplicityCorrection::kNSD; eventType++)
  {
    AliMultiplicityCorrection* mult = AliMultiplicityCorrection::Open((eventType == AliMultiplicityCorrection::kINEL) ? "chi2_inel.root" : "chi2_nsd.root");
    TH1* mcHist = mult->GetMultiplicityMC(etaRange, eventType)->ProjectionY("mymc");
    TH1* result = mult->GetMultiplicityESDCorrected(etaRange);
  
    DrawResultRatio(mcHist, result, Form("finalPlotCheck_%d.eps", (Int_t) eventType));

    // normalize result
    //result->Scale(1.0 / result->Integral(2, displayRange));
  
    result->GetXaxis()->SetRangeUser(0, displayRange);
    //result->SetBinContent(1, 0); result->SetBinError(1, 0);
    result->SetTitle(Form(";True multiplicity in |#eta| < %.1f;Entries", (etaRange+1) * 0.5));
    result->SetMarkerStyle(0);
    result->SetLineColor(1);
    result->SetStats(kFALSE);
  
    // systematic error
    TH1* error = SystematicsSummary(tpc, (eventType == AliMultiplicityCorrection::kNSD));
    
    TH1* systError = (TH1*) result->Clone("systError");
    for (Int_t i=1; i<=systError->GetNbinsX(); ++i)
      systError->SetBinError(i, systError->GetBinContent(i) * error->GetBinContent(i));
  
    // change error drawing style
    systError->SetFillColor(15);
    
    if (eventType == AliMultiplicityCorrection::kNSD)
    {
      result->SetLineColor(2);
      result->SetMarkerColor(2);
      result->SetMarkerStyle(5);
    }
    
    canvas->cd();
    systError->DrawCopy(Form("E2 ][ %s", (eventType == AliMultiplicityCorrection::kINEL) ? "" : "SAME"));
    result->DrawCopy("SAME E ][");
    canvas->SetLogy();
    
    legend->AddEntry(result, (eventType == AliMultiplicityCorrection::kINEL) ? "Inelastic cross-section" : "NSD cross-section", (eventType == AliMultiplicityCorrection::kINEL) ? "L" : "P");
  }
  
  legend->Draw();
  /*
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
  */

/*  TText* text = new TText(10, 1e-4, "Systematic errors summed quadratically");
  text->SetTextSize(0.04);
  text->DrawText(10, 5e-5, "0.6 #cdot 10^{6} minimum bias events");
  text->DrawText(10, 3e-5, "TPC tracks in |#eta| < 0.9");
  text->DrawText(10, 1e-5, "corrected to ineleastic events in |#eta| < 0.9");
  text->Draw();*/


  canvas->SaveAs(canvas->GetName());
}

void finalPlot2(Bool_t tpc = 0)
{
  loadlibs();

  if (tpc)
    SetTPC();

  //TH1* errorNSD = SystematicsSummary(tpc, 1);

  TCanvas* canvas = new TCanvas("finalPlot2.eps", "finalPlot2.eps", 800, 600); 
  canvas->SetRightMargin(0.05);
  canvas->SetTopMargin(0.05);
  canvas->SetGridx();
  canvas->SetGridy();
  canvas->SetLogy();
  
  legend = new TLegend(0.5, 0.7, 0.9, 0.9);
  legend->SetFillColor(0);
  legend->SetTextSize(0.04);
  
  Int_t displayRanges[] = { 30, 45, 65 };
  
  TH2* dummy = new TH2F("dummy", ";True multiplicity N_{ch};P(N_{ch})", 100, -0.5, displayRanges[2]+10, 1000, 5e-6, 5);
  dummy->SetStats(0);
  dummy->Draw();
  
  TList list;
  
  // get syst error
  file = TFile::Open(Form("%s_syst_error.root", (tpc) ? "tpc" : "spd"));
  TH1* totalINEL = (TH1*) file->Get("inel_1");
  TH1* totalNSD = (TH1*) file->Get("nsd_1");
  
  Int_t count = 0;
  for (AliMultiplicityCorrection::EventType eventType = AliMultiplicityCorrection::kINEL; eventType <= AliMultiplicityCorrection::kNSD; eventType++)
  {
    AliMultiplicityCorrection* mult = AliMultiplicityCorrection::Open((eventType == AliMultiplicityCorrection::kINEL) ? "chi2_inel.root" : "chi2_nsd.root");
    //TH1* mcHist = mult->GetMultiplicityMC(etaRange, eventType)->ProjectionY("mymc");
    for (Int_t etaR = 0; etaR <= 2; etaR++)
    {
      TH1* result = mult->GetMultiplicityESDCorrected(etaR);
    
      //DrawResultRatio(mcHist, result, Form("finalPlotCheck_%d.eps", (Int_t) eventType));
  
      // normalize result
      result->Scale(TMath::Power(5, etaR) / result->Integral(1, displayRanges[etaR]));
    
      result->GetXaxis()->SetRangeUser(0, displayRanges[etaR]);
      //result->SetBinContent(1, 0); result->SetBinError(1, 0);
      //result->SetTitle(Form(""));
      result->SetMarkerStyle(0);
      result->SetLineColor(1);
      result->SetLineWidth(2);
      //result->SetMarkerStyle(4);
      //result->SetStats(kFALSE);
    
      // systematic error
      //TH1* error = SystematicsSummary(tpc, (eventType == AliMultiplicityCorrection::kNSD));
      TH1* error = (eventType == AliMultiplicityCorrection::kNSD) ? totalNSD : totalINEL;
      
      TH1* systError = (TH1*) result->Clone("systError");
      // small hack until we have syst. errors for all eta ranges
      Float_t factor = 80.0 / displayRanges[etaR];
      for (Int_t i=1; i<=systError->GetNbinsX(); ++i)
      {
        systError->SetBinError(i, systError->GetBinContent(i) * error->GetBinContent(factor * i));
      }
    
      // change error drawing style
      systError->SetFillColor(15);
      
      if (eventType == AliMultiplicityCorrection::kNSD)
      {
        result->SetLineColor(2);
        result->SetMarkerColor(2);
        result->SetMarkerStyle(5);
      }
      
      canvas->cd();
      systError->GetXaxis()->SetRangeUser(0, displayRanges[etaR]);
      systError->DrawCopy("E2 ][ SAME");
      list.Add(result);
      
      if (eventType == AliMultiplicityCorrection::kINEL)
      {
        TLatex* Tl = new TLatex;
        Tl->SetTextSize(0.04);
        //Tl->SetBit(TLatex::kTextNDC);
        Tl->SetTextAlign(32);

        Float_t etaRangeArr[] = { 0.5, 1.0, 1.4 };
        TString tmpStr;
        tmpStr.Form("|#eta| < %.1f (x %d)", etaRangeArr[etaR], (Int_t) TMath::Power(5, etaR));
        if (etaR == 0)
          Tl->DrawLatex(15, result->GetBinContent(20), tmpStr);
        if (etaR == 1)
        {
          Tl->SetTextAlign(12);
          Tl->DrawLatex(40, result->GetBinContent(40), tmpStr);
        }
        if (etaR == 2)
        {
          Tl->SetTextAlign(12);
          Tl->DrawLatex(60, result->GetBinContent(50), tmpStr);
        }
      }
      
      if (etaR == 0)
        legend->AddEntry(result, (eventType == AliMultiplicityCorrection::kINEL) ? "Inelastic events" : "NSD events", (eventType == AliMultiplicityCorrection::kINEL) ? "L" : "P");
      
      count++;
    }
  }
  
  for (Int_t i=0; i<list.GetEntries(); i++)
    ((TH1*) list.At(i))->DrawCopy("SAME E ][");
  
  legend->Draw();

  canvas->SaveAs(canvas->GetName());
}

TMatrixD* NonInvertable()
{
  const Int_t kSize = 5;

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

  matrix.Print();

  //TMatrixD inverted(matrix);
  //inverted.Invert();
  
  //inverted.Print();
  
  return new TMatrixD(matrix);
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

  TH1F* inputDist = new TH1F("inputDist", ";t;Entries", kSize, -0.5, (Float_t) kSize - 0.5);
  TVectorD inputDistVector(kSize);
  TH1F* unfolded = inputDist->Clone("unfolded");
  TH1F* measuredIdealDist = inputDist->Clone("measuredIdealDist");
  measuredIdealDist->SetTitle(";m;Entries");
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

  TCanvas* canvas = new TCanvas("BlobelUnfoldingExample", "BlobelUnfoldingExample", 1200, 600);
  canvas->SetTopMargin(0.05);
  canvas->Divide(2, 1);

  canvas->cd(1);
  canvas->cd(1)->SetLeftMargin(0.15);
  canvas->cd(1)->SetRightMargin(0.05);
  canvas->cd(1)->SetTopMargin(0.05);
  gPad->SetGridx();
  gPad->SetGridy();
  measuredDist->GetYaxis()->SetRangeUser(-600, 2799);
  measuredDist->GetYaxis()->SetTitleOffset(1.9);
  measuredDist->SetStats(0);
  measuredDist->DrawCopy();
  gaus->Draw("SAME");

  canvas->cd(2);
  canvas->cd(2)->SetLeftMargin(0.15);
  canvas->cd(2)->SetRightMargin(0.05);
  canvas->cd(2)->SetTopMargin(0.05);
  gPad->SetGridx();
  gPad->SetGridy();
  unfolded->GetYaxis()->SetRangeUser(-600, 2799);
  unfolded->GetYaxis()->SetTitleOffset(1.9);
  unfolded->SetStats(0);
  unfolded->DrawCopy();
  gaus->Draw("SAME");

  canvas->SaveAs("BlobelUnfoldingExample.eps");
  
  return;
  
  // now unfold this with Bayesian
  loadlibs();
  
  // fill a multiplicity object
  mult = new AliMultiplicityCorrection("mult", "mult");
  for (Int_t x=0; x<kSize; x++)
  {
    mult->GetMultiplicityVtx(0)->SetBinContent(1, x+1, inputDistVector(x));
    mult->GetMultiplicityESD(0)->SetBinContent(1, x+1, measuredDistVector(x)*10000);
    for (Int_t y=0; y<kSize; y++)
      mult->GetCorrelation(0)->SetBinContent(1, x+1, y+1, matrix(x, y));
  }
  
  //mult->DrawHistograms();
  
  mult->SetRegularizationParameters(AliMultiplicityCorrection::kPol0, 0);
  //mult->SetCreateBigBin(kFALSE);
  mult->ApplyMinuitFit(0, kFALSE, AliMultiplicityCorrection::kTrVtx, kFALSE); //hist2->ProjectionY("mymchist"));
  
  //mult->ApplyBayesianMethod(0, kFALSE, AliMultiplicityCorrection::kTrVtx, 0, -1, 0, kFALSE);
  
  mult->DrawComparison("BlobelExample", 0, kFALSE, kTRUE, mult->GetMultiplicityVtx(0)->ProjectionY("mcmchist", 1, mult->GetMultiplicityVtx(0)->GetNbinsX()));
}

void TestWithDiagonalMatrix()
{
  const Int_t kSize = 20;

  TMatrixD matrix(kSize, kSize);
  for (Int_t x=0; x<kSize; x++)
    matrix(x, x) = 1;
 
  if (1)
  { 
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
  }

  //matrix.Print();

  TH1F* inputDist = new TH1F("inputDist", ";t;Entries", kSize, -0.5, (Float_t) kSize - 0.5);
  TVectorD inputDistVector(kSize);
  
  TH1F* measuredIdealDist = inputDist->Clone("measuredIdealDist");
  measuredIdealDist->SetTitle(";m;Entries");
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

  // randomize
  //measuredDist->FillRandom(measuredIdealDist, 10000);
  measuredDist = measuredIdealDist;
  measuredDist->Sumw2();
  
  for (Int_t x=0; x<kSize; x++)
    Printf("bin %d %.2f +- %.2f", x+1, measuredDist->GetBinContent(x+1), measuredDist->GetBinError(x+1));

  // now unfold this 
  loadlibs();
  
  // fill a multiplicity object
  mult = new AliMultiplicityCorrection("mult", "mult");
  for (Int_t x=0; x<kSize; x++)
  {
    mult->GetMultiplicityVtx(0)->SetBinContent(1, x+1, inputDistVector(x));
    mult->GetMultiplicityESD(0)->SetBinContent(1, x+1, measuredDist->GetBinContent(x+1));
    for (Int_t y=0; y<kSize; y++)
      mult->GetCorrelation(0)->SetBinContent(1, x+1, y+1, matrix(x, y));
  }
  
  //mult->DrawHistograms();
  
  AliUnfolding::SetNbins(20, 20);
  AliUnfolding::SetChi2Regularization(AliUnfolding::kNone, 0);
  AliUnfolding::SetChi2Regularization(AliUnfolding::kPol1, 1e-14);
  //mult->SetCreateBigBin(kFALSE);
  mult->ApplyMinuitFit(0, kFALSE, AliMultiplicityCorrection::kTrVtx, kFALSE); //hist2->ProjectionY("mymchist"));
  
  //mult->ApplyBayesianMethod(0, kFALSE, AliMultiplicityCorrection::kTrVtx, 0, -1, 0, kFALSE);
  
  mult->DrawComparison("TestWithDiagonalMatrix", 0, kFALSE, kTRUE, mult->GetMultiplicityVtx(0)->ProjectionY("mcmchist", 1, mult->GetMultiplicityVtx(0)->GetNbinsX()));
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
  fCurrentESD->SetMarkerStyle(0);
  fCurrentESD->SetLineColor(1);
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

void DifferentSamples()
{
  loadlibs();

  Int_t n = 2;
  const char* filesChi2[] = { "chi2_100k_1.root", "chi2_100k_2.root" };
  const char* filesBayesian[] = { "bayesian_100k_1.root", "bayesian_100k_2.root" };

  TCanvas* canvas = new TCanvas("DifferentSamples", "DifferentSamples", 1200, 600);
  canvas->Divide(2, 1);
  
  legend = new TLegend(0.15, 0.7, 0.65, 0.9);
  legend->SetFillColor(0);
  legend->SetTextSize(0.04);

  for (Int_t i=0; i<n; i++)
  {
    AliMultiplicityCorrection* chi2 = new AliMultiplicityCorrection("Multiplicity", "Multiplicity");
    TFile::Open(filesChi2[i]);
    chi2->LoadHistograms("Multiplicity");

    AliMultiplicityCorrection* bayesian = new AliMultiplicityCorrection("Multiplicity", "Multiplicity");
    TFile::Open(filesBayesian[i]);
    bayesian->LoadHistograms("Multiplicity");
    
    chi2Hist = chi2->GetMultiplicityESDCorrected(etaRange);
    bayesianHist = bayesian->GetMultiplicityESDCorrected(etaRange);
    
    mc = chi2->GetMultiplicityVtx(etaRange)->ProjectionY("mymc", 1, 1);
    
    // normalize and divide
    chi2Hist->Scale(1.0 / chi2Hist->Integral(2, displayRange+1) * mc->Integral(2, displayRange));
    bayesianHist->Scale(1.0 / bayesianHist->Integral(2, displayRange+1) * mc->Integral(2, displayRange));
    
    chi2Hist->Divide(mc, chi2Hist);
    bayesianHist->Divide(mc, bayesianHist);
    
    canvas->cd(i+1);
    gPad->SetTopMargin(0.05);
    gPad->SetRightMargin(0.05);
    //gPad->SetLeftMargin(0.12);
    gPad->SetGridx();
    gPad->SetGridy(); 
    
    chi2Hist->GetXaxis()->SetRangeUser(0, displayRange);
    chi2Hist->GetYaxis()->SetTitleOffset(1.3);
    chi2Hist->SetStats(0);
    chi2Hist->SetTitle(Form(";%s;MC / unfolded", GetMultLabel()));
    chi2Hist->GetYaxis()->SetRangeUser(0.2, 1.8);
    chi2Hist->Draw("HIST");
    
    for (Int_t x=1; x<=bayesianHist->GetNbinsX(); x++)
      bayesianHist->SetBinError(x, 1e-6);
    
    bayesianHist->SetLineColor(2);
    bayesianHist->SetMarkerColor(2);
    bayesianHist->SetMarkerStyle(5);
    bayesianHist->Draw("HIST E SAME");
    
    if (i == 0)
    {
      legend->AddEntry(chi2Hist, "#chi^{2}-minimization", "L");
      legend->AddEntry(bayesianHist, "Bayesian unfolding", "LP");
    }
    legend->Draw();
  }
  
  canvas->SaveAs("DifferentSamples.eps");
}

void PileUp()
{
  loadlibs();
  
  AliMultiplicityCorrection* mult = AliMultiplicityCorrection::Open("multiplicityMC.root");
  hist2d = mult->GetMultiplicityMC(etaRange, AliMultiplicityCorrection::kINEL);
  mult1 = hist2d->ProjectionY("mult1", 1, hist2d->GetNbinsX());
  
  conv = (TH1*) mult1->Clone("conv");
  conv->Reset();
  
  mult1->Scale(1.0 / mult1->Integral());
  
  for (Int_t i=1; i<=mult1->GetNbinsX(); i++)
    for (Int_t j=1; j<=mult1->GetNbinsX(); j++)
      conv->Fill(mult1->GetBinCenter(i)+mult1->GetBinCenter(j), mult1->GetBinContent(i) * mult1->GetBinContent(j));
  
  conv->Scale(1.0 / conv->Integral());
  
  c = new TCanvas("c", "c", 800, 500);
  gPad->SetLogy();
  gPad->SetTopMargin(0.05);
  gPad->SetRightMargin(0.05);
  mult1->SetTitle(Form(";%s;Probability", GetMultLabel()));
  mult1->SetStats(0);
  gPad->SetGridx();
  gPad->SetGridy();
  mult1->Draw();
  mult1->GetYaxis()->SetRangeUser(1e-7, 2 * mult1->GetMaximum());
  mult1->GetXaxis()->SetRangeUser(0, displayRange);
  mult1->GetXaxis()->SetTitleOffset(1.15);
  conv->SetLineColor(2);
  conv->SetMarkerColor(2);
  conv->SetMarkerStyle(5);
  conv->DrawCopy("SAME P");
  
  conv->Scale(0.00058);
  conv->DrawCopy("SAME P");
  
  legend = new TLegend(0.73, 0.73, 0.93, 0.93);
  legend->SetFillColor(0);
  legend->SetTextSize(0.04);
  legend->AddEntry(mult1, "1 collision");
  legend->AddEntry(conv, "2 collisions", "P");
  legend->Draw();
  
  c->SaveAs("pileup.eps");

  new TCanvas;
  conv->Divide(mult1);
  conv->Draw();
}

void TestErrorDetermination(Int_t nRandomizations)
{
  TF1* func = new TF1("nbd", "exp(log([0]) + TMath::LnGamma([2]+x) - TMath::LnGamma([2]) - TMath::LnGamma(x+1) + log([1] / ([1]+[2])) * x + log(1.0 + [1]/[2]) * -[2])", 0, 100);
  func->SetParNames("scaling", "averagen", "k");
  func->SetParameters(1, 15, 2);
  
  TF1* func2 = new TF1("nbd2", "exp(log([0]) + TMath::LnGamma([2]+x) - TMath::LnGamma([2]) - TMath::LnGamma(x+1) + log([1] / ([1]+[2])) * x + log(1.0 + [1]/[2]) * -[2])", 0, 100);
  func2->SetParNames("scaling", "averagen", "k");
  func2->SetParLimits(0, 0.5, 2);
  func2->SetParLimits(1, 1, 50);
  func2->SetParLimits(2, 1, 10);
  func2->SetParameters(1, 15, 2);
  //func2->FixParameter(0, 1);
  
  //new TCanvas; func->Draw("L");
  
  hist1 = new TH1F("hist1", "", 100, 0.5, 100.5);
  hist2 = new TH1F("hist2", "", 100, 0.5, 100.5);
  hist1->Sumw2();
  
  TH1* params[3];
  params[0] = new TH1F("param_0", Form("param_%d", 0), 100, 0.95, 1.05);
  params[1] = new TH1F("param_1", Form("param_%d", 1), 100, 14, 16);
  params[2] = new TH1F("param_2", Form("param_%d", 2), 100, 1.8, 2.2);
  
  const Int_t nTries = 1000;
  for (Int_t i=0; i<nTries; i++)
  {
    hist1->Reset();
    
    if (nRandomizations == 1)
    {
      hist1->FillRandom("nbd", 10000);
    }
    else if (nRandomizations == 2)
    {
      hist2->Reset();
      hist2->FillRandom("nbd", 10000);
      hist1->FillRandom(hist2, 10000);
    }
    else if (nRandomizations == 3)
    {
      hist2->Reset();
      hist1->FillRandom("nbd", 10000);
      hist2->FillRandom(hist1, 10000);
      hist1->Reset();
      hist1->FillRandom(hist2, 10000);
    }
    else
      return;
  
    //new TCanvas; hist1->Draw();
  
    hist1->Scale(1.0 / hist1->Integral());
    hist1->Fit(func2, "NQ");
    hist1->Fit(func2, "NQ");
    for (Int_t j=0; j<3; j++)
      params[j]->Fill(func2->GetParameter(j));
  }
  
  for (Int_t j=0; j<3; j++)
  {
    new TCanvas; params[j]->Draw();
    params[j]->Fit("gaus");
    Printf("sigma of param %d if %f", j, ((TF1*) params[j]->FindObject("gaus"))->GetParameter(2));
  }
}

void DrawRawDistributions(const char* fileName = "multiplicityESD.root")
{
  loadlibs();

  AliMultiplicityCorrection* mult = AliMultiplicityCorrection::Open(fileName);
  
  c = new TCanvas("c", "c", 600, 600);
  
  dummy = new TH2F("dummy", ";measured multiplicity", 100, -0.5, 149.5, 100, 0.5, 4e4);
  dummy->SetStats(0);
  dummy->Draw();
  gPad->SetGridx();
  gPad->SetGridy();
  gPad->SetLogy();
  
  Int_t colors[] = { 1, 2, 4 };
  
  for (Int_t i=2; i>=0; i--)
  {
    hist = mult->GetMultiplicityESD(i)->ProjectionY();
    
    hist->SetLineColor(colors[i]);
    hist->DrawCopy("SAME");
  }
  
  
}

void FindUnfoldedLimit()
{
  loadlibs();
  
  AliMultiplicityCorrection* mult = AliMultiplicityCorrection::Open("multiplicityMC.root");
  AliMultiplicityCorrection* esd = AliMultiplicityCorrection::Open("multiplicityESD.root");
  
  TH1* hist = mult->GetCorrelation(etaRange);
  for (Int_t y=0; y<=hist->GetYaxis()->GetNbins()+1; ++y)
  {
    for (Int_t z=0; z<=hist->GetZaxis()->GetNbins()+1; ++z)
    {
      hist->SetBinContent(0, y, z, 0);
      hist->SetBinContent(hist->GetXaxis()->GetNbins()+1, y, z, 0);
    }
  }
  TH2* corr = (TH2*) ((TH3*) mult->GetCorrelation(etaRange))->Project3D("zy");

  TH1* esd_proj = esd->GetMultiplicityESD(etaRange)->ProjectionY("esd_proj", 1, esd->GetMultiplicityESD(etaRange)->GetNbinsX());
  //esd_proj = corr->ProjectionY("esd_proj", 1, corr->GetNbinsX());
  
  new TCanvas; corr->Draw("COLZ");
  new TCanvas; esd_proj->DrawCopy();
  
  TH1* percentage = (TH1*) (esd_proj->Clone("percentage"));
  percentage->SetTitle("percentage");
  percentage->Reset();
  
  for (Int_t i=1; i<=esd_proj->GetNbinsX(); i++)
    if (esd_proj->GetBinContent(i) > 0)
      esd_proj->SetBinContent(i, 1);
  
  Int_t limit = -1;
  for (Int_t i=1; i<=percentage->GetNbinsX(); i++)
  {
    TH1* binResponse = corr->ProjectionY("proj", i, i);
    if (binResponse->Integral() <= 0)
      continue;
    binResponse->Scale(1.0 / binResponse->Integral());
    binResponse->Multiply(esd_proj);
    //new TCanvas; binResponse->Draw();
    Float_t value = binResponse->Integral();
    percentage->SetBinContent(i, value);
    if (limit == -1 && value < 0.9)
      limit = percentage->GetXaxis()->GetBinCenter(i-1);
    //return;
  }
  
  Printf("Limit is %d", limit);
  
  new TCanvas; percentage->Draw();
  
  new TCanvas;
  mc = esd->GetMultiplicityVtx(etaRange)->ProjectionY("mc", 1, esd->GetMultiplicityVtx(etaRange)->GetNbinsX());
  mc->SetLineColor(2);
  mc->Draw("");
}
   
void FindUnfoldedLimitAll()
{
  loadlibs();
  for (etaRange = 0; etaRange <= 2; etaRange++)
    FindUnfoldedLimit();
}
   
void CompareUnfoldedWithMC()
{
  loadlibs();
  
  AliMultiplicityCorrection* mult = AliMultiplicityCorrection::Open("unfolded.root");

  //mult->GetMultiplicityESDCorrected(etaRange)->Scale(1.0 / mult->GetMultiplicityESDCorrected(etaRange)->Integral());
  mult->GetMultiplicityESDCorrected(etaRange)->SetTitle(";multiplicity;events");
  mult->GetMultiplicityESDCorrected(etaRange)->SetStats(0);
  mult->GetMultiplicityESDCorrected(etaRange)->Draw();
  //mcHist->Scale(1.0 / mcHist->Integral()); 
  mcHist = mult->GetMultiplicityVtx(etaRange)->ProjectionY("mymc", 1, 1);
  mcHist->SetLineColor(2);
  mcHist->Draw("SAME");
  gPad->SetLogy();
}

void PaperNumbers()
{
  const char* label[] = { "SD", "DD", "ND" };
  const char* files[] = { "multiplicity_syst_sd.root", "multiplicity_syst_dd.root", "multiplicity_syst_nd.root" };
  
  loadlibs();
  
  Printf("vertex reco");
  
  Float_t oneParticle[3];
  for (Int_t i=0; i<3; i++)
  {
    AliMultiplicityCorrection* mult = AliMultiplicityCorrection::Open(files[i]);
    
    eff = mult->GetEfficiency(2, AliMultiplicityCorrection::kMB);
    //eff = mult->GetTriggerEfficiency(2);
    
    oneParticle[i] = 100.0 * eff->GetBinContent(2);
    
    Printf("%s: %.2f", label[i], oneParticle[i]);
  }
  
  Float_t ua5_SD = 0.153;
  Float_t ua5_DD = 0.080;
  Float_t ua5_ND = 0.767;
  
  Float_t vtxINELUA5 = ua5_SD * oneParticle[0] + ua5_DD * oneParticle[1] + ua5_ND * oneParticle[2];
  Float_t vtxNSDUA5  = (ua5_DD * oneParticle[1] + ua5_ND * oneParticle[2]) / (ua5_DD + ua5_ND);
  
  Printf("INEL (UA5)  = %.1f", vtxINELUA5);
  Printf("NSD (UA5)  = %.1f", vtxNSDUA5);
  
  Printf("total for >= 2 charged tracks");
  
  AliMultiplicityCorrection* mult = AliMultiplicityCorrection::Open("multiplicity.root");
  
  vtx = mult->GetMultiplicityVtx(2)->ProjectionY("vtx", 1, mult->GetMultiplicityVtx(2)->GetNbinsX(), "e");
  all = mult->GetMultiplicityINEL(2)->ProjectionY("all", 1, mult->GetMultiplicityINEL(2)->GetNbinsX(), "e"); 
  
  Int_t begin = vtx->GetXaxis()->FindBin(2);
  Int_t end = vtx->GetNbinsX();
  
  Float_t above2 = vtx->Integral(begin, end) / all->Integral(begin, end);
  
  Printf("%.1f", 100.0 * above2);
}

void DrawRawEta()
{
  TFile::Open(measuredFile);
  
  Int_t colors[] = {1, 2, 4, 6, 7, 8, 9, 10};
  
  for (Int_t i=0; i<3; i++)
  {
    hist = dynamic_cast<TH1*> (gFile->Get(Form("fEta_%d", i)));
    
    hist->GetYaxis()->SetRangeUser(0, 8);
    hist->SetLineColor(colors[i]);
    hist->SetStats(kFALSE);
    hist->Draw((i == 0) ? "HIST" : "HIST SAME");
  }   
  
  gPad->SetGridx();
  gPad->SetGridy();
}

void CrossSectionUncertainties(Int_t xsectionID, Int_t eventType)
{
  const char* files[] = { "multiplicitySD.root", "multiplicityDD.root", "multiplicityND.root" };

  loadlibs();

  AliMultiplicityCorrection* data[3];

  Float_t ref_SD = -1;
  Float_t ref_DD = -1;
  Float_t ref_ND = -1;

  Float_t error_SD = -1;
  Float_t error_DD = -1;
  Float_t error_ND = -1;

  gROOT->ProcessLine(gSystem->ExpandPathName(".L $ALICE_ROOT/PWG0/dNdEta/drawSystematics.C"));
  GetRelativeFractions(xsectionID, ref_SD, ref_DD, ref_ND, error_SD, error_DD, error_ND);
  
  TH1* unc[9*3];
  TH1* unc2[9*3];
  const char* legendStrings[8];
  
  for (Int_t i=0; i<9; i++)
  {
    Float_t factorSD = 0;
    Float_t factorDD = 0;
    
    if (i > 0 && i < 3)
      factorSD = (i % 2 == 0) ? 1 : -1;
    else if (i >= 3 && i < 5)
      factorDD = (i % 2 == 0) ? 1 : -1;
    else if (i >= 5 && i < 9)
    {
      factorSD = ((i % 2 == 0) ? 1.0 : -1.0) / TMath::Sqrt(2);
      if (i == 5 || i == 6)
        factorDD = factorSD;
      else
        factorDD = -factorSD;
    }
    
    Float_t scalesSD = ref_SD + factorSD * error_SD;
    Float_t scalesDD = ref_DD + factorDD * error_DD;
    Float_t scalesND = 1.0 - scalesDD[i] - scalesSD[i];
    
    str = new TString;
    str->Form("Case %d: SD: %.2f DD: %.2f ND: %.2f", i, scalesSD, scalesDD, scalesND);
    Printf("%s", str->Data());
    if (i > 0)
      legendStrings[i-1] = str->Data();
  
    for (Int_t j=0; j<3; ++j)
    {
      TString name;
      name.Form("Multiplicity_%d", j);
  
      TFile::Open(files[j]);
      data[j] = new AliMultiplicityCorrection(name, name);
      data[j]->LoadHistograms("Multiplicity");
    }
    
    // calculating relative
    Float_t sd = data[0]->GetMultiplicityINEL(0)->Integral(0, data[0]->GetMultiplicityINEL(0)->GetNbinsX()+1);
    Float_t dd = data[1]->GetMultiplicityINEL(0)->Integral(0, data[1]->GetMultiplicityINEL(0)->GetNbinsX()+1);
    Float_t nd = data[2]->GetMultiplicityINEL(0)->Integral(0, data[2]->GetMultiplicityINEL(0)->GetNbinsX()+1);
    Float_t total = nd + dd + sd;
    
    nd /= total;
    sd /= total;
    dd /= total;
    
    Printf("Ratios in the correction map are: ND=%f, DD=%f, SD=%f", nd, dd, sd);
    
    Float_t ratio[3];
    ratio[0] = scalesSD / sd;
    ratio[1] = scalesDD / dd;
    ratio[2] = scalesND / nd;
    
    Printf("SD=%.2f, DD=%.2f, ND=%.2f", ratio[0], ratio[1], ratio[2]);
    
    TList list;
    for (Int_t k=0; k<3; ++k)
    {
      // modify x-section
      for (Int_t j=0; j<AliMultiplicityCorrection::kMCHists; j++)
      {
        data[k]->GetMultiplicityVtx(j)->Scale(ratio[k]);
        data[k]->GetMultiplicityMB(j)->Scale(ratio[k]);
        data[k]->GetMultiplicityINEL(j)->Scale(ratio[k]);
        data[k]->GetMultiplicityNSD(j)->Scale(ratio[k]);
      }
  
      if (k > 0)
        list.Add(data[k]);
    }
  
    printf("Case %d, %s: Entries in response matrix: SD: %.2f DD: %.2f ND: %.2f", xsectionID, data[0]->GetName(), data[0]->GetCorrelation(0)->Integral(), data[1]->GetCorrelation(0)->Integral(), data[2]->GetCorrelation(0)->Integral());
  
    data[0]->Merge(&list);
    data[0]->FixTriggerEfficiencies(20);
  
    Printf(" Total: %.2f", data[0]->GetCorrelation(0)->Integral());
  
    for (Int_t etaR = 0; etaR < 3; etaR++)
    {
      unc[i+(etaR*9)] = (TH1*) data[0]->GetTriggerEfficiency(etaR, eventType)->Clone(Form("eff_%d_%d", etaR, i));
      unc2[i+(etaR*9)] = (TH1*) data[0]->GetEfficiency(etaR, eventType)->Clone(Form("eff_%d_%d", etaR, i));
    }
  }
      
  file = TFile::Open(Form("systunc_fractions_%s.root", (eventType == 2) ? "inel" : "nsd"), "RECREATE");
  
  for (Int_t etaR = 0; etaR < 3; etaR++)
  {
    largestError1 = (TH1*) unc[0]->Clone(Form("fractions_trig_%d", etaR));
    largestError1->Reset();
    largestError2 = (TH1*) unc[0]->Clone(Form("fractions_trigvtx_%d", etaR));
    largestError2->Reset();
  
    etaRange = etaR; // correct labels in DrawRatio
    DrawRatio(unc[etaR*9],  8, unc+1+etaR*9,  Form("xsections_trig%d.png", etaR),    kFALSE, legendStrings, kFALSE, largestError1);
    DrawRatio(unc2[etaR*9], 8, unc2+1+etaR*9, Form("xsections_trigvtx%d.png", etaR), kFALSE, legendStrings, kFALSE, largestError2);
    
    largestError2->SetBinContent(1, largestError1->GetBinContent(1));
    
    largestError2->GetYaxis()->UnZoom();
    largestError2->Write();
    
    new TCanvas;
    largestError2->DrawCopy();
  }
  file->Close();
}

void PythiaPhojetUncertainties()
{
  PythiaPhojetUncertainties(900, 2);
  PythiaPhojetUncertainties(900, 3);
  PythiaPhojetUncertainties(2360, 2);
  PythiaPhojetUncertainties(2360, 3);
}

void PythiaPhojetUncertainties(Int_t energy, Int_t eventType)
{
  loadlibs();
  
  AliMultiplicityCorrection* mult[2];

  const char* trigger = "all";
  if (eventType == 3)
    trigger = "v0and";
    
  if (energy == 2360)
    trigger = "spdgfobits";

  if (energy == 7000)
    trigger = "all-onepart";
  
  if (energy == 900)
    TFile::Open(Form("LHC10a12_run10482X/%s/spd/multiplicityMC_xsection.root", trigger));
  else if (energy == 2360)
    TFile::Open(Form("LHC10a10_run105054_7/%s/spd/multiplicityMC_xsection.root", trigger));
  else
    TFile::Open(Form("LHC10b2_run114931/%s/spd/multiplicity.root", trigger));
  mult[0] = new AliMultiplicityCorrection("Multiplicity", "Multiplicity");
  mult[0]->LoadHistograms("Multiplicity");
  
  if (energy == 900)
    TFile::Open(Form("LHC10a14_run10482X_Phojet/%s/spd/multiplicityMC_xsection.root", trigger));
  else if (energy == 2360)
    TFile::Open(Form("LHC10a15_run10505X_Phojet/%s/spd/multiplicityMC_xsection.root", trigger));
  else
    TFile::Open(Form("LHC10b1_run114931/%s/spd/multiplicity.root", trigger));
  mult[1] = new AliMultiplicityCorrection("Multiplicity", "Multiplicity");
  mult[1]->LoadHistograms("Multiplicity");
  
  TH1* unc[6];
  TH1* unc2[6];
  
  for (Int_t i=0; i<2; i++)
  {
    mult[i]->FixTriggerEfficiencies(20);
    for (Int_t etaR = 0; etaR < 3; etaR++)
    {
      unc[i+(etaR*2)] = (TH1*) mult[i]->GetTriggerEfficiency(etaR, eventType)->Clone(Form("eff_%d_%d", etaR, i));
      unc2[i+(etaR*2)] = (TH1*) mult[i]->GetEfficiency(etaR, eventType)->Clone(Form("eff_%d_%d", etaR, i));
    }
  }
  
  file = TFile::Open(Form("systunc_pythiaphojet_%d_%s.root", energy, (eventType == 2) ? "inel" : "nsd"), "RECREATE");
  
  for (Int_t etaR = 0; etaR < 3; etaR++)
  {
    largestError1Low = (TH1*) unc[0]->Clone(Form("pythiaphojet_trig_%d_low", etaR));
    largestError1Low->Reset();
    largestError1High = (TH1*) largestError1Low->Clone(Form("pythiaphojet_trig_%d_high", etaR));
    largestError2Low  = (TH1*) largestError1Low->Clone(Form("pythiaphojet_trigvtx_%d_low", etaR));
    largestError2High = (TH1*) largestError1Low->Clone(Form("pythiaphojet_trigvtx_%d_high", etaR));
  
    DrawRatio(unc[etaR*2], 1, unc+1+etaR*2, Form("pythiaphojet_trig%d.png", etaR), kFALSE, 0, kFALSE, largestError1Low, largestError1High);
    DrawRatio(unc2[etaR*2], 1, unc2+1+etaR*2, Form("pythiaphojet_trigvtx%d.png", etaR), kFALSE, 0, kFALSE, largestError2Low, largestError2High);
    
    largestError2Low->SetBinContent(1, largestError1Low->GetBinContent(1));
    largestError2High->SetBinContent(1, largestError1High->GetBinContent(1));
    
    largestError2Low->GetYaxis()->UnZoom();
    largestError2High->GetYaxis()->UnZoom();
    largestError2Low->Write();
    largestError2High->Write();
    
    new TCanvas;
    largestError2Low->DrawCopy();
    
    new TCanvas;
    largestError2High->DrawCopy();
  }
  file->Close();
}

void CombinatoricsAsFunctionOfN(Int_t etaR)
{
  loadlibs();
  
  mult = AliMultiplicityCorrection::Open("multiplicityESD.root");
  
  //Float_t integratedCombinatorics[] = { 0.00062, 0.00074, 0.001 }; // 900 GeV
  Float_t integratedCombinatorics[] = { 7.5e-4, 1.1e-3, 1.2e-3 };   // 2.36 TeV
    
  multHist = mult->GetMultiplicityESD(etaR)->ProjectionY();
  multHist->Scale(1.0 / multHist->Integral());
  
  Float_t integral = 0;
  Float_t tracks = 0;
  for (Int_t i=1; i<= multHist->GetNbinsX(); i++)
  {
    Float_t x = multHist->GetXaxis()->GetBinCenter(i);
    integral += x*x * multHist->GetBinContent(i+1);
    tracks += x * multHist->GetBinContent(i+1);
    Printf("%f %f %f", x, multHist->GetBinContent(i+1), integral);
  }
    
  Printf("Integral is %f; %f", integral, integral / tracks);
  Printf("alpha is %e", integratedCombinatorics[etaR] / (integral / tracks));
}

void TryCombinatoricsCorrection(Float_t etaR, Float_t alpha)
{
  loadlibs();
  
  mult = AliMultiplicityCorrection::Open("multiplicity.root");
  
  multHist = mult->GetMultiplicityESD(etaR)->ProjectionY();
  
  target = (TH1*) multHist->Clone("target");
  target->Reset();
  
  Int_t totalTracks1 = 0;
  Int_t totalTracks2 = 0;
  
  for (Int_t i=1; i<=multHist->GetNbinsX(); i++)
  {
    for (Int_t j=0; j<multHist->GetBinContent(i); j++)
    {
      Int_t origTracks = (Int_t) multHist->GetXaxis()->GetBinCenter(i);
      tracks = origTracks;
      //tracks -= gRandom->Poisson(1.3e-4 * origTracks * origTracks);
      if (gRandom->Uniform() < alpha * origTracks * origTracks)
        tracks--;
      target->Fill(tracks);
      
      totalTracks1 += origTracks;
      totalTracks2 += tracks;
    }
  }
  
  new TCanvas; multHist->Draw(); target->Draw("SAME"); target->SetLineColor(2);

  Printf("%d %d %f %f", totalTracks1, totalTracks2, 1. - (Float_t) totalTracks2 / totalTracks1, 1. - target->GetMean() / multHist->GetMean());
}

void ApplyCombinatoricsCorrection()
{
  loadlibs();
  
  mult = AliMultiplicityCorrection::Open("multiplicity.root");
  
  fractionMC = new TF1("fractionMC", "[0]+[1]*x", 0, 200);
  fractionMC->SetParameters(0.0004832, 5.82e-5); // mariella's presentation 16.04.10
  
  fractionData = new TF1("fractionData", "[0]+[1]*x", 0, 200);
  fractionData->SetParameters(0.00052, 0.0001241); // mariella's presentation 16.04.10
  
  Int_t etaR = 1;
  esd = mult->GetMultiplicityESD(etaR);
  
  for (Int_t x=1; x<=esd->GetNbinsX(); x++)
  {
    for (Int_t y=2; y<=esd->GetNbinsY(); y++)
    {
      printf("Before: %f  %f ", esd->GetBinContent(x, y-1), esd->GetBinContent(x, y));
    
      Float_t n = esd->GetYaxis()->GetBinCenter(y);
      Int_t addComb = (fractionData->Eval(n) - fractionMC->Eval(n)) * n * esd->GetBinContent(x, y);
      
      // shift addComb down one bin
      esd->SetBinContent(x, y-1, esd->GetBinContent(x, y-1) + addComb);
      esd->SetBinContent(x, y,   esd->GetBinContent(x, y) - addComb);
      
      // recalc errors
      esd->SetBinError(x, y-1, TMath::Sqrt(esd->GetBinContent(x, y-1)));
      esd->SetBinError(x, y, TMath::Sqrt(esd->GetBinContent(x, y)));
      
      Printf("comb: %d; after: %f +- %f  %f +- %f", addComb, esd->GetBinContent(x, y-1), esd->GetBinError(x, y-1), esd->GetBinContent(x, y), esd->GetBinError(x, y));
    }
  }
  
  TFile::Open("out.root", "RECREATE");
  mult->SaveHistograms();
}
