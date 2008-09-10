/* $Id$ */

#if !defined(__CINT__) || defined(__MAKECINT__)

#include "AliPWG0Helper.h"
#include "dNdEtaAnalysis.h"
#include "AlidNdEtaCorrection.h"

#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TLine.h>
#include <TSystem.h>

#endif

Int_t gMax = 5;

extern TSystem* gSystem;

void loadlibs()
{
  gSystem->Load("libANALYSIS");
  gSystem->Load("libPWG0base");
}

void SetRanges(TAxis* axis)
{
  if (strcmp(axis->GetTitle(), "#eta") == 0)
    axis->SetRangeUser(-1.7999, 1.7999);
  if (strcmp(axis->GetTitle(), "p_{T} [GeV/c]") == 0)
    axis->SetRangeUser(0, 4.9999);
  if (strcmp(axis->GetTitle(), "vtx z [cm]") == 0)
    axis->SetRangeUser(-15, 14.9999);
  if (strcmp(axis->GetTitle(), "Ntracks") == 0)
    axis->SetRangeUser(0, 99.9999);
}

void SetRanges(TH1* hist)
{
  SetRanges(hist->GetXaxis());
  SetRanges(hist->GetYaxis());
  SetRanges(hist->GetZaxis());
}


void Prepare3DPlot(TH3* hist)
{
  hist->GetXaxis()->SetTitleOffset(1.5);
  hist->GetYaxis()->SetTitleOffset(1.5);
  hist->GetZaxis()->SetTitleOffset(1.5);

  hist->SetStats(kFALSE);
}

void Prepare2DPlot(TH2* hist)
{
  hist->SetStats(kFALSE);
  hist->GetYaxis()->SetTitleOffset(1.4);

  hist->SetMinimum(0);
  hist->SetMaximum(gMax);

  SetRanges(hist);
}

void Prepare1DPlot(TH1* hist)
{
  hist->SetLineWidth(2);
  hist->SetStats(kFALSE);

  hist->GetXaxis()->SetLabelOffset(0.02);
  hist->GetXaxis()->SetTitleOffset(1.3);
  hist->GetYaxis()->SetTitleOffset(1.3);

  SetRanges(hist);
}

void InitPad()
{
  gPad->Range(0, 0, 1, 1);
  gPad->SetLeftMargin(0.15);
  //gPad->SetRightMargin(0.05);
  //gPad->SetTopMargin(0.13);
  gPad->SetBottomMargin(0.12);

  gPad->SetGridx();
  gPad->SetGridy();
}

void InitPadCOLZ()
{
  gPad->Range(0, 0, 1, 1);
  gPad->SetRightMargin(0.15);
  gPad->SetLeftMargin(0.12);
  gPad->SetTopMargin(0.05);

  gPad->SetGridx();
  gPad->SetGridy();
}

// --- end of helpers --- begin functions ---

void DrawOverview(const char* fileName = "correction_map.root", const char* dirName = "dndeta_correction")
{
  loadlibs();
  if (!TFile::Open(fileName))
    return;

  AlidNdEtaCorrection* dNdEtaCorrection = new AlidNdEtaCorrection(dirName, dirName);
  if (!dNdEtaCorrection->LoadHistograms())
    return;

  dNdEtaCorrection->Finish();

  dNdEtaCorrection->DrawOverview();
}

void PrintInfo(const char* fileName = "correction_map.root", const char* dirName = "dndeta_correction")
{
  loadlibs();
  if (!TFile::Open(fileName))
    return;

  AlidNdEtaCorrection* dNdEtaCorrection = new AlidNdEtaCorrection(dirName, dirName);
  if (!dNdEtaCorrection->LoadHistograms())
    return;

  dNdEtaCorrection->Finish();

  for (Int_t i=AlidNdEtaCorrection::kTrack2Particle; i<=AlidNdEtaCorrection::kND; i++)
  {
    Printf("Correction %d", i);
    dNdEtaCorrection->GetCorrection(i)->PrintInfo(0.3);
  }
}

void PrintAllInfos()
{
  PrintInfo();

  Printf("RAW ESD");
  TFile::Open("analysis_esd_raw.root");
  dNdEtaAnalysis* fdNdEtaAnalysis = new dNdEtaAnalysis("dndeta", "dndeta");
  fdNdEtaAnalysis->LoadHistograms("fdNdEtaAnalysisESD");
  fdNdEtaAnalysis->GetData()->PrintInfo(0.3);

  const Int_t num = 3;
  const char* results[] = { "dndeta", "dndetaTr", "dndetaTrVtx" };

  TFile::Open("analysis_esd.root");
  for (Int_t i=0; i<num; i++)
  {
    Printf("CORRECTED %s", results[i]);
    dNdEtaAnalysis* fdNdEtaAnalysis = new dNdEtaAnalysis("dndeta", "dndeta");
    fdNdEtaAnalysis->LoadHistograms(results[i]);
    fdNdEtaAnalysis->GetData()->PrintInfo(0.3);
  }
}  

void ComparedNdEta(const char* ESDfolder = "dndeta", const char* MCfolder = "dndeta", const char* esdFile = "analysis_esd.root", const char* mcFile = "analysis_mc.root")
{
  gSystem->Load("libPWG0base");

  TFile::Open(esdFile);
  dNdEtaAnalysis* fdNdEtaAnalysisESD = new dNdEtaAnalysis(ESDfolder, ESDfolder);
  fdNdEtaAnalysisESD->LoadHistograms();

  TFile::Open(mcFile);
  dNdEtaAnalysis* fdNdEtaAnalysisMC = new dNdEtaAnalysis(MCfolder, MCfolder);
  fdNdEtaAnalysisMC->LoadHistograms();
  //fdNdEtaAnalysisMC->Finish(0, 0.3, AlidNdEtaCorrection::kNone);

  for (Int_t i=0; i<dNdEtaAnalysis::kVertexBinning; ++i)
    fdNdEtaAnalysisESD->GetdNdEtaPtCutOffCorrectedHistogram(i)->Divide(fdNdEtaAnalysisMC->GetdNdEtaPtCutOffCorrectedHistogram(i));

  fdNdEtaAnalysisESD->DrawHistograms();
}

void CompareVertexDist(Int_t plot = 0, const char* correctionMapFile = "correction_map.root", const char* correctionMapFolder = "dndeta_correction")
{
  gSystem->Load("libPWG0base");

  const char* ESDfolder = 0;

  if (plot == 0) // all
    ESDfolder = "dndeta";
  else if (plot == 1) // mb
    ESDfolder = "dndeta_mb";
  else if (plot == 2) // mb vtx
    ESDfolder = "dndeta_mbvtx";

  TFile::Open("analysis_esd.root");
  dNdEtaAnalysis* fdNdEtaAnalysisESD = new dNdEtaAnalysis(ESDfolder, ESDfolder);
  fdNdEtaAnalysisESD->LoadHistograms();

  AlidNdEtaCorrection* dNdEtaCorrection = new AlidNdEtaCorrection(correctionMapFolder, correctionMapFolder);
  dNdEtaCorrection->LoadHistograms(correctionMapFile, correctionMapFolder);

  TH2F* hist = 0;

  if (plot == 0) // all
    hist = dNdEtaCorrection->GetTriggerBiasCorrection()->GetGeneratedHistogram();
  else if (plot == 1) // mb
    hist = dNdEtaCorrection->GetTriggerBiasCorrection()->GetMeasuredHistogram();
  else if (plot == 2) // mb vtx
    hist = dNdEtaCorrection->GetVertexRecoCorrection()->GetMeasuredHistogram();

  TH1* proj = hist->ProjectionX();

  TH1* vertex = fdNdEtaAnalysisESD->GetVtxHistogram();
  for (Int_t i=1; i<=vertex->GetNbinsX(); ++i)
  {
    Float_t value = proj->GetBinContent(proj->FindBin(vertex->GetBinCenter(i)));
    if (value != 0)
    {
      printf("vtx = %f, esd = %f, corr = %f, ratio = %f\n", vertex->GetBinCenter(i), vertex->GetBinContent(i), value, vertex->GetBinContent(i) / value);
      vertex->SetBinContent(i, vertex->GetBinContent(i) / value);
    }
  }

  new TCanvas;
  vertex->DrawCopy();
}

void CompareTrack2ParticleWithAnalysisData(const char* correctionMapFile = "correction_map.root", const char* correctionMapFolder = "dndeta_correction")
{
  gSystem->Load("libPWG0base");

  TFile::Open("analysis_esd.root");
  dNdEtaAnalysis* fdNdEtaAnalysisESD = new dNdEtaAnalysis("dndeta_mbvtx", "dndeta_mbvtx");
  fdNdEtaAnalysisESD->LoadHistograms();

  AlidNdEtaCorrection* dNdEtaCorrection = new AlidNdEtaCorrection(correctionMapFolder, correctionMapFolder);
  dNdEtaCorrection->LoadHistograms(correctionMapFile, correctionMapFolder);

  //TH1* histESD = fdNdEtaAnalysisESD->GetUncorrectedHistogram();
  //TH1* histCorr = dNdEtaCorrection->GetTrack2ParticleCorrection()->GetMeasuredHistogram();

  TH1* histESD = fdNdEtaAnalysisESD->GetHistogram();
  TH1* histCorr = dNdEtaCorrection->GetTrack2ParticleCorrection()->GetGeneratedHistogram();

  TH1F* diff = new TH1F("diff", "diff", 100, 0.98, 1.02);

  new TCanvas;
  histESD->Draw();

  new TCanvas;
  histCorr->Draw();

  for (Int_t x=1; x<=histESD->GetNbinsX(); ++x)
    for (Int_t y=1; y<=histESD->GetNbinsY(); ++y)
      for (Int_t z=1; z<=histESD->GetNbinsZ(); ++z)
      {
        Float_t value1 = histESD->GetBinContent(x, y, z);
        Float_t value2 = histCorr->GetBinContent(histCorr->FindBin(histESD->GetXaxis()->GetBinCenter(x), histESD->GetYaxis()->GetBinCenter(y), histESD->GetZaxis()->GetBinCenter(z)));

        if (value2 > 0 && value1 > 0)
        {
          printf("%f %f %f\n", value1, value2, value1 / value2);
          diff->Fill(value1 / value2);
        }
      }

  new TCanvas;
  diff->Draw();
}

Double_t PrintIntegratedDeviation(TH1* histMC, TH1* histESD, const char* desc = "")
{
  Double_t avgMC = 0;
  Double_t avgESD = 0;
  for (Int_t bin = histMC->FindBin(-0.7999); bin <= histMC->FindBin(0.7999); bin++)
  {
    avgMC += histMC->GetBinContent(bin);
    avgESD += histESD->GetBinContent(bin);
  }
  Int_t nBins = histMC->FindBin(0.7999) - histMC->FindBin(-0.7999) + 1;

  avgMC /= nBins;
  avgESD /= nBins;

  // deviation when integrate in |eta| < 0.8 between mc and esd
  Double_t diffFullRange = (avgMC - avgESD) / avgMC;

  Printf("%s: Integrated deviation in |eta| < 0.8 is %.2f %%", desc, diffFullRange * 100);

  return diffFullRange;
}

void dNdEtaNoResolution()
{
  loadlibs();

  TFile::Open("correction_map.root");

  const char* correctionMapFolder = "dndeta_correction";
  AlidNdEtaCorrection* dNdEtaCorrection = new AlidNdEtaCorrection(correctionMapFolder, correctionMapFolder);
  dNdEtaCorrection->LoadHistograms();
  dNdEtaCorrection->GetTrack2ParticleCorrection()->PrintInfo(0.3);

  TFile::Open("analysis_mc.root");
  fdNdEtaAnalysis = new dNdEtaAnalysis("dndetaTracks", "dndetaTracks");
  fdNdEtaAnalysis->LoadHistograms();
  fdNdEtaAnalysis->Finish(dNdEtaCorrection, 0.3, AlidNdEtaCorrection::kTrack2Particle, "ESD (no resolution effect) -> MB with vertex");
  fdNdEtaAnalysis->GetdNdEtaPtCutOffCorrectedHistogram(0)->SetMarkerStyle(21);

  TFile::Open("analysis_mc.root");
  dNdEtaAnalysis* fdNdEtaAnalysisMC = new dNdEtaAnalysis("dndetaTrVtx", "dndetaTrVtx");
  fdNdEtaAnalysisMC->LoadHistograms();
  fdNdEtaAnalysisMC->Finish(0, 0, AlidNdEtaCorrection::kNone, "MC: MB with vertex");

  DrawdNdEtaRatio(fdNdEtaAnalysis->GetdNdEtaPtCutOffCorrectedHistogram(0), fdNdEtaAnalysisMC->GetdNdEtaPtCutOffCorrectedHistogram(0), "MB with vertex (no resolution effect)", 3);
}

TH1* GetMCHist(const char* folder, Float_t ptCut, const char* tag)
{
  dNdEtaAnalysis* fdNdEtaAnalysis = new dNdEtaAnalysis(folder, folder);
  fdNdEtaAnalysis->LoadHistograms();
  fdNdEtaAnalysis->Finish(0, ptCut, AlidNdEtaCorrection::kNone, tag);
  return fdNdEtaAnalysis->GetdNdEtaHistogram(0);
}

void dNdEta(Bool_t onlyESD = kFALSE, Bool_t save = kTRUE)
{
  TFile* file = TFile::Open("analysis_esd.root");
  TH1* histESD = (TH1*) file->Get("dndeta/dNdEta_corrected");
  TH1* histESDnsd = (TH1*) file->Get("dndetaNSD/dNdEta_corrected");
  TH1* histESDNoPt = (TH1*) file->Get("dndeta/dNdEta");
  TH1* histESDMB = (TH1*) file->Get("dndetaTr/dNdEta_corrected");
  TH1* histESDMBNoPt = (TH1*) file->Get("dndetaTr/dNdEta");
  TH1* histESDMBVtx = (TH1*) file->Get("dndetaTrVtx/dNdEta_corrected");
  TH1* histESDMBVtxNoPt = (TH1*) file->Get("dndetaTrVtx/dNdEta");
  TH1* histESDMBTracksNoPt = (TH1*) file->Get("dndetaTracks/dNdEta");

  Prepare1DPlot(histESD);
  Prepare1DPlot(histESDnsd);
  Prepare1DPlot(histESDMB);
  Prepare1DPlot(histESDMBVtx);

  Prepare1DPlot(histESDNoPt);
  Prepare1DPlot(histESDMBNoPt);
  Prepare1DPlot(histESDMBVtxNoPt);
  Prepare1DPlot(histESDMBTracksNoPt);

  histESD->SetLineWidth(0);
  histESDnsd->SetLineWidth(0);
  histESDMB->SetLineWidth(0);
  histESDMBVtx->SetLineWidth(0);

  histESDNoPt->SetLineWidth(0);
  histESDMBNoPt->SetLineWidth(0);
  histESDMBVtxNoPt->SetLineWidth(0);

  histESD->SetMarkerColor(1);
  histESDnsd->SetMarkerColor(6);
  histESDMB->SetMarkerColor(2);
  histESDMBVtx->SetMarkerColor(3);

  histESD->SetLineColor(1);
  histESDnsd->SetLineColor(6);
  histESDMB->SetLineColor(2);
  histESDMBVtx->SetLineColor(3);

  histESDNoPt->SetMarkerColor(1);
  histESDMBNoPt->SetMarkerColor(2);
  histESDMBVtxNoPt->SetMarkerColor(3);
  histESDMBTracksNoPt->SetMarkerColor(4);

  histESD->SetMarkerStyle(20);
  histESDnsd->SetMarkerStyle(29);
  histESDMB->SetMarkerStyle(21);
  histESDMBVtx->SetMarkerStyle(22);

  histESDNoPt->SetMarkerStyle(20);
  histESDMBNoPt->SetMarkerStyle(21);
  histESDMBVtxNoPt->SetMarkerStyle(22);
  histESDMBTracksNoPt->SetMarkerStyle(23);
  
  //Float_t etaLimit = 1.2999;
  Float_t etaLimit = 2.41;
  Float_t etaPlotLimit = 2.6;

  histESDMBVtx->GetXaxis()->SetRangeUser(-etaLimit, etaLimit);
  histESDMB->GetXaxis()->SetRangeUser(-etaLimit, etaLimit);
  histESD->GetXaxis()->SetRangeUser(-etaLimit, etaLimit);
  histESDnsd->GetXaxis()->SetRangeUser(-etaLimit, etaLimit);

  histESDNoPt->GetXaxis()->SetRangeUser(-etaLimit, etaLimit);
  histESDMBNoPt->GetXaxis()->SetRangeUser(-etaLimit, etaLimit);
  histESDMBVtxNoPt->GetXaxis()->SetRangeUser(-etaLimit, etaLimit);
  histESDMBTracksNoPt->GetXaxis()->SetRangeUser(-etaLimit, etaLimit);

  Float_t max = TMath::Max(histESDMBVtx->GetMaximum(), histESDMB->GetMaximum());
  max = TMath::Max(max, histESD->GetMaximum());

  TH2F* dummy = new TH2F("dummy", "", 100, -etaPlotLimit, etaPlotLimit, 1000, 0, max * 1.1);
  Prepare1DPlot(dummy);
  dummy->SetStats(kFALSE);
  dummy->SetXTitle("#eta");
  dummy->SetYTitle("dN_{ch}/d#eta");
  dummy->GetYaxis()->SetTitleOffset(1);

  TCanvas* canvas = new TCanvas("dNdEta1", "dNdEta1", 500, 500);

  dummy->DrawCopy();
  histESDMBVtx->Draw("SAME");
  histESDMB->Draw("SAME");
  histESD->Draw("SAME");

  if (save)
  {
    canvas->SaveAs("dNdEta1.gif");
    canvas->SaveAs("dNdEta1.eps");
  }

  if (onlyESD)
    return;

  loadlibs();

  TFile* file2 = TFile::Open("analysis_mc.root");

  TH1* histMC =            (TH1*) GetMCHist("dndeta", -1, "MC: full inelastic")->Clone("histMC");
  TH1* histMCTr =          (TH1*) GetMCHist("dndetaTr", -1, "MC: minimum bias")->Clone("histMCTr");
  TH1* histMCTrVtx =       (TH1*) GetMCHist("dndetaTrVtx", -1, "MC: MB with trigger")->Clone("histMCTrVtx");
  TH1* histMCnsd =         (TH1*) GetMCHist("dndetaNSD", -1, "MC: NSD")->Clone("histMCnsd");

  TH1* histMCPtCut =       (TH1*) GetMCHist("dndeta", 0.3, "MC: full inelastic, pt cut")->Clone("histMCPtCut");
  TH1* histMCTrPtCut =     (TH1*) GetMCHist("dndetaTr", 0.3, "MC: minimum bias, pt cut")->Clone("histMCTrPtCut");
  TH1* histMCTrVtxPtCut =  (TH1*) GetMCHist("dndetaTrVtx", 0.3, "MC: MB with trigger, pt cut")->Clone("histMCTrVtxPtCut");
  TH1* histMCTracksPtCut = (TH1*) GetMCHist("dndetaTracks", 0.3, "MC: Tracks w/o resolution effect, pt cut")->Clone("histMCTracksPtCut");

  Prepare1DPlot(histMC);
  Prepare1DPlot(histMCnsd);
  Prepare1DPlot(histMCTr);
  Prepare1DPlot(histMCTrVtx);

  Prepare1DPlot(histMCPtCut);
  Prepare1DPlot(histMCTrPtCut);
  Prepare1DPlot(histMCTrVtxPtCut);
  Prepare1DPlot(histMCTracksPtCut);

  histMC->GetXaxis()->SetRangeUser(-etaLimit, etaLimit);
  histMCnsd->GetXaxis()->SetRangeUser(-etaLimit, etaLimit);
  histMCTr->GetXaxis()->SetRangeUser(-etaLimit, etaLimit);
  histMCTrVtx->GetXaxis()->SetRangeUser(-etaLimit, etaLimit);

  histMCPtCut->GetXaxis()->SetRangeUser(-etaLimit, etaLimit);
  histMCTrPtCut->GetXaxis()->SetRangeUser(-etaLimit, etaLimit);
  histMCTrVtxPtCut->GetXaxis()->SetRangeUser(-etaLimit, etaLimit);
  histMCTracksPtCut->GetXaxis()->SetRangeUser(-etaLimit, etaLimit);

  histMC->SetLineColor(1);
  histMCnsd->SetLineColor(6);
  histMCTr->SetLineColor(2);
  histMCTrVtx->SetLineColor(3);

  histMCPtCut->SetLineColor(1);
  histMCTrPtCut->SetLineColor(2);
  histMCTrVtxPtCut->SetLineColor(3);
  if (histMCTracksPtCut)
    histMCTracksPtCut->SetLineColor(4);

  TCanvas* canvas2 = new TCanvas("dNdEta2", "dNdEta2", 500, 500);

  TH2* dummy2 = (TH2F*) dummy->Clone("dummy2");
  dummy2->GetYaxis()->SetRangeUser(0, max * 1.1);

  dummy2->DrawCopy();
  histMC->Draw("SAME");
  histMCnsd->Draw("SAME");
  histMCTr->Draw("SAME");
  histMCTrVtx->Draw("SAME");
  histESD->Draw("SAME");
  histESDnsd->Draw("SAME");
  histESDMB->Draw("SAME");
  histESDMBVtx->Draw("SAME");
  histESDNoPt->Draw("SAME");
  histESDMBNoPt->Draw("SAME");
  histESDMBVtxNoPt->Draw("SAME");
  histESDMBTracksNoPt->Draw("SAME");
  histMCPtCut->Draw("SAME");
  histMCTrPtCut->Draw("SAME");
  histMCTrVtxPtCut->Draw("SAME");
  if (histMCTracksPtCut)
    histMCTracksPtCut->Draw("SAME");

  if (save)
  {
    canvas2->SaveAs("dNdEta2.gif");
    canvas2->SaveAs("dNdEta2.eps");
  }

  DrawdNdEtaRatio(histESD, histMC, "full_inelastic", etaPlotLimit);
  DrawdNdEtaRatio(histESDMB, histMCTr, "triggered", etaPlotLimit);
  DrawdNdEtaRatio(histESDMBVtx, histMCTrVtx, "triggered_vertex", etaPlotLimit);
  DrawdNdEtaRatio(histESDnsd, histMCnsd, "NSD", etaPlotLimit);

  new TCanvas;
  dummy2->DrawCopy();
  histMCnsd->Draw("SAME");
  histESDnsd->Draw("SAME");

  TH1* ratio = (TH1*) histMC->Clone("ratio");
  TH1* ratioNoPt = (TH1*) histMCPtCut->Clone("ratioNoPt");

  ratio->Divide(histESD);
  ratioNoPt->Divide(histESDNoPt);

  ratio->GetXaxis()->SetRangeUser(-etaLimit, etaLimit);

  ratio->SetLineColor(1);
  ratioNoPt->SetLineColor(2);

  Double_t average = 0;       // average deviation from 1 in ratio (depends on the number of bins if statistical)
  for (Int_t bin = ratio->FindBin(-0.7999); bin <= ratio->FindBin(0.7999); bin++)
    average += TMath::Abs(ratio->GetBinContent(bin) - 1);
  Int_t nBins = ratio->FindBin(0.7999) - ratio->FindBin(-0.7999) + 1;
  average /= nBins;
  Printf("Average deviation in |eta| < 0.8 is %.2f %%", average * 100);

  PrintIntegratedDeviation(histMC, histESD, "all events");
  PrintIntegratedDeviation(histMCnsd, histESDnsd, "all events (NSD)");
  PrintIntegratedDeviation(histMCTr, histESDMB, "triggered");
  PrintIntegratedDeviation(histMCTrVtx, histESDMBVtx, "trigger, vertex");
  PrintIntegratedDeviation(histMCPtCut, histESDNoPt, "all events (no pt corr)");
  PrintIntegratedDeviation(histMCTrPtCut, histESDMBNoPt, "triggered (no pt corr)");
  PrintIntegratedDeviation(histMCTrVtxPtCut, histESDMBVtxNoPt, "trigger, vertex (no pt corr)");

  TCanvas* canvas3 = new TCanvas("dNdEta", "dNdEta", 700, 600);
  canvas3->Range(0, 0, 1, 1);
  //canvas3->Divide(1, 2, 0, 0);

  //canvas3->cd(1);
  TPad* pad1 = new TPad("dNdEta_1", "", 0, 0.5, 0.98, 0.98);
  pad1->Draw();

  TPad* pad2 = new TPad("dNdEta_2", "", 0, 0.02, 0.98, 0.5);
  pad2->Draw();

  pad1->SetRightMargin(0.05);
  pad2->SetRightMargin(0.05);

  // no border between them
  pad1->SetBottomMargin(0);
  pad2->SetTopMargin(0);

  pad1->cd();

  TLegend* legend = new TLegend(0.4, 0.05, 0.65, 0.3);
  legend->SetFillColor(0);
  legend->AddEntry(histESDMBVtx, "triggered, vertex");
  legend->AddEntry(histESDMB, "triggered");
  legend->AddEntry(histESD, "all events");
  legend->AddEntry(histMC, "MC prediction");

  dummy->GetXaxis()->SetLabelSize(0.06);
  dummy->GetYaxis()->SetLabelSize(0.06);
  dummy->GetXaxis()->SetTitleSize(0.06);
  dummy->GetYaxis()->SetTitleSize(0.06);
  dummy->GetYaxis()->SetTitleOffset(0.7);
  dummy->DrawCopy();
  histESDMBVtx->Draw("SAME");
  histESDMB->Draw("SAME");
  histESD->Draw("SAME");
  histMC->Draw("SAME");

  legend->Draw();

  pad2->cd();
  pad2->SetBottomMargin(0.15);

  Float_t minR = 0.9; //TMath::Min(0.961, ratio->GetMinimum() * 0.95);
  Float_t maxR = 1.1; //TMath::Max(1.049, ratio->GetMaximum() * 1.05);

  TH1F dummy3("dummy3", ";#eta;Ratio: MC / ESD", 100, -etaPlotLimit, etaPlotLimit);
  dummy3.SetStats(kFALSE);
  for (Int_t i=1; i<=100; ++i)
    dummy3.SetBinContent(i, 1);
  dummy3.GetYaxis()->SetRangeUser(minR, maxR);
  dummy3.SetLineWidth(2);
  dummy3.GetXaxis()->SetLabelSize(0.06);
  dummy3.GetYaxis()->SetLabelSize(0.06);
  dummy3.GetXaxis()->SetTitleSize(0.06);
  dummy3.GetYaxis()->SetTitleSize(0.06);
  dummy3.GetYaxis()->SetTitleOffset(0.7);
  dummy3.DrawCopy();

  ratio->Draw("SAME");

  //pad2->Draw();

  canvas3->Modified();

  if (save)
  {
    canvas3->SaveAs("dNdEta.gif");
    canvas3->SaveAs("dNdEta.eps");
  }

  TCanvas* canvas4 = new TCanvas("ratio", "ratio", 700, 500);

  ratio->Draw();
  ratioNoPt->Draw("SAME");

  TLegend* legend = new TLegend(0.6, 0.7, 0.95, 0.9);
  legend->SetFillColor(0);
  legend->AddEntry(ratio, "mc/esd");
  legend->AddEntry(ratioNoPt, "mc/esd, not pt cut off corrected");
  legend->Draw();
}

void DrawdNdEtaRatio(TH1* corr, TH1* mc, const char* name, Float_t etaPlotLimit)
{
  TCanvas* canvas3 = new TCanvas(name, name, 700, 600);
  canvas3->Range(0, 0, 1, 1);

  TPad* pad1 = new TPad(Form("%s_1", name), "", 0, 0.5, 0.98, 0.98);
  pad1->Draw();

  TPad* pad2 = new TPad(Form("%s_2", name), "", 0, 0.02, 0.98, 0.5);
  pad2->Draw();

  pad1->SetRightMargin(0.05);
  pad2->SetRightMargin(0.05);

  // no border between them
  pad1->SetBottomMargin(0);
  pad2->SetTopMargin(0);

  pad1->cd();

  TLegend* legend = new TLegend(0.4, 0.05, 0.65, 0.3);
  legend->SetFillColor(0);
  legend->AddEntry(corr, "corrected");
  legend->AddEntry(mc, "MC prediction");

  TH2F* dummy = new TH2F("dummy", "", 100, -etaPlotLimit, etaPlotLimit, 1000, 0, corr->GetMaximum() * 1.1);
  Prepare1DPlot(dummy);
  dummy->SetStats(kFALSE);
  dummy->SetXTitle("#eta");
  dummy->SetYTitle("dN_{ch}/d#eta");
  dummy->GetYaxis()->SetTitleOffset(1);

  dummy->GetXaxis()->SetLabelSize(0.06);
  dummy->GetYaxis()->SetLabelSize(0.06);
  dummy->GetXaxis()->SetTitleSize(0.06);
  dummy->GetYaxis()->SetTitleSize(0.06);
  dummy->GetYaxis()->SetTitleOffset(0.7);
  dummy->DrawCopy();

  corr->Draw("SAME");
  mc->Draw("SAME");

  legend->Draw();

  pad2->cd();
  pad2->SetBottomMargin(0.15);

  TH1* ratio = (TH1*) mc->Clone("ratio");
  ratio->Divide(corr);

  Float_t minR = TMath::Min(0.961, ratio->GetMinimum() * 0.95);
  Float_t maxR = TMath::Max(1.049, ratio->GetMaximum() * 1.05);

  TH1F dummy3("dummy3", ";#eta;Ratio: MC / corr", 100, -etaPlotLimit, etaPlotLimit);
  dummy3.SetStats(kFALSE);
  for (Int_t i=1; i<=100; ++i)
  	dummy3.SetBinContent(i, 1);
  dummy3.GetYaxis()->SetRangeUser(minR, maxR);
  dummy3.SetLineWidth(2);
  dummy3.GetXaxis()->SetLabelSize(0.06);
  dummy3.GetYaxis()->SetLabelSize(0.06);
  dummy3.GetXaxis()->SetTitleSize(0.06);
  dummy3.GetYaxis()->SetTitleSize(0.06);
  dummy3.GetYaxis()->SetTitleOffset(0.7);
  dummy3.DrawCopy();

  ratio->Draw("SAME");

  canvas3->Modified();
}

void ptSpectrum()
{
  TFile* file = TFile::Open("analysis_esd.root");
  TH1* histESD = (TH1*) file->Get("dndeta/dndeta_pt");

  TFile* file2 = TFile::Open("analysis_mc.root");
  TH1* histMC = (TH1*) file2->Get("dndeta/dndeta_pt");

  TCanvas* canvas = new TCanvas("ptSpectrum", "ptSpectrum", 500, 500);
  InitPad();
  gPad->SetLogy();

  Prepare1DPlot(histMC);
  Prepare1DPlot(histESD);

  histESD->SetTitle("");
  histESD->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  histESD->GetYaxis()->SetTitle("#frac{dN}{d#eta dp_{T}} [c/GeV]");

  histMC->SetLineColor(kBlue);
  histESD->SetLineColor(kRed);

  histESD->GetYaxis()->SetTitleOffset(1.5);
  histESD->GetXaxis()->SetRangeUser(0, 4.9999);

  histESD->SetMaximum(TMath::Max(histESD->GetMaximum(), histMC->GetMaximum()) * 2);

  histESD->Draw();
  histMC->Draw("SAME");

  canvas->SaveAs("ptSpectrum.gif");
  canvas->SaveAs("ptSpectrum.eps");
}

void TriggerBiasVtxRecon(const char* fileName = "correction_map.root", const char* folder = "dndeta_correction")
{
  gSystem->Load("libPWG0base");

  TFile::Open(fileName);
  AlidNdEtaCorrection* dNdEtaCorrection = new AlidNdEtaCorrection("dndeta_correction", "dndeta_correction");
  dNdEtaCorrection->LoadHistograms();

  TH2* corrTrigger = dNdEtaCorrection->GetTriggerBiasCorrectionINEL()->GetEventCorrection()->GetCorrectionHistogram();
  TH2* corrVtx = dNdEtaCorrection->GetVertexRecoCorrection()->GetEventCorrection()->GetCorrectionHistogram();

  Prepare2DPlot(corrTrigger);
  corrTrigger->SetTitle("b) Trigger bias correction");

  Prepare2DPlot(corrVtx);
  corrVtx->SetTitle("a) Vertex reconstruction correction");

  corrTrigger->GetYaxis()->SetTitle("Multiplicity");
  corrVtx->GetYaxis()->SetTitle("Multiplicity");

  TCanvas* canvas = new TCanvas("TriggerBiasVtxRecon", "TriggerBiasVtxRecon", 1000, 500);
  canvas->Divide(2, 1);

  canvas->cd(1);
  InitPadCOLZ();
  corrVtx->DrawCopy("COLZ");

  canvas->cd(2);
  InitPadCOLZ();
  corrTrigger->DrawCopy("COLZ");

  canvas->SaveAs(Form("TriggerBiasVtxRecon_%d.gif", gMax));
  canvas->SaveAs(Form("TriggerBiasVtxRecon_%d.eps", gMax));

  canvas = new TCanvas("TriggerBiasVtxReconZoom", "TriggerBiasVtxReconZoom", 1000, 500);
  canvas->Divide(2, 1);

  corrTrigger->GetYaxis()->SetRangeUser(0, 5);
  corrVtx->GetYaxis()->SetRangeUser(0, 5);

  canvas->cd(1);
  InitPadCOLZ();
  corrVtx->DrawCopy("COLZ");

  canvas->cd(2);
  InitPadCOLZ();
  corrTrigger->DrawCopy("COLZ");

  canvas->SaveAs(Form("TriggerBiasVtxReconZoom_%d.gif", gMax));
  canvas->SaveAs(Form("TriggerBiasVtxReconZoom_%d.eps", gMax));
}

void TriggerBias(const char* fileName = "correction_map.root")
{
  TFile* file = TFile::Open(fileName);

  TH2* corr = dynamic_cast<TH2*> (file->Get("dndeta_correction/corr_dndeta_correction_trigger"));

  Prepare2DPlot(corr);
  corr->SetTitle("Trigger bias correction");

  TCanvas* canvas = new TCanvas("TriggerBias", "TriggerBias", 500, 500);
  InitPadCOLZ();
  corr->DrawCopy("COLZ");

  canvas->SaveAs(Form("TriggerBias_%d.gif", gMax));
  canvas->SaveAs(Form("TriggerBias_%d.eps", gMax));

  corr->GetYaxis()->SetRangeUser(0, 5);

  canvas = new TCanvas("TriggerBiasZoom", "TriggerBiasZoom", 500, 500);
  InitPadCOLZ();
  corr->DrawCopy("COLZ");

  canvas->SaveAs(Form("TriggerBiasZoom_%d.gif", gMax));
  canvas->SaveAs(Form("TriggerBiasZoom_%d.eps", gMax));
}

void TriggerBias1D(const char* fileName = "correction_map.root", const char* folderName = "dndeta_correction")
{
  gSystem->Load("libPWG0base");

  TFile* file = TFile::Open(fileName);
  AlidNdEtaCorrection* dNdEtaCorrection = new AlidNdEtaCorrection(folderName, folderName);
  dNdEtaCorrection->LoadHistograms();

  TH1* hist = dNdEtaCorrection->GetTriggerBiasCorrectionINEL()->GetEventCorrection()->Get1DCorrection("x");
  TH1* hist2 = dNdEtaCorrection->GetTriggerBiasCorrectionINEL()->GetEventCorrection()->Get1DCorrection("y", -10, 10);

  TCanvas* canvas = new TCanvas("TriggerBias1D", "TriggerBias1D", 1000, 500);
  canvas->Divide(2, 1);

  canvas->cd(1);
  InitPad();

  Prepare1DPlot(hist);
  hist->SetTitle("");
  hist->GetYaxis()->SetTitle("correction factor");
  hist->GetYaxis()->SetRangeUser(1, 1.5);
  hist->GetYaxis()->SetTitleOffset(1.6);
  hist->Draw();

  canvas->cd(2);
  InitPad();

  Prepare1DPlot(hist2);
  hist2->SetTitle("");
  hist2->GetYaxis()->SetTitle("correction factor");
  hist2->GetXaxis()->SetRangeUser(0, 5);
  hist2->GetYaxis()->SetTitleOffset(1.6);
  hist2->GetXaxis()->SetTitle("multiplicity");
  hist2->Draw();

  TPaveText* pave = new TPaveText(0.6, 0.8, 0.8, 0.85, "NDC");
  pave->SetFillColor(0);
  pave->AddText("|z| < 10 cm");
  pave->Draw();

  canvas->SaveAs("TriggerBias1D.eps");
}

void VtxRecon()
{
  TFile* file = TFile::Open("correction_map.root");

  TH2* corr = dynamic_cast<TH2*> (file->Get("dndeta_correction/corr_dndeta_correction_vtxReco"));

  Prepare2DPlot(corr);
  corr->SetTitle("Vertex reconstruction correction");

  TCanvas* canvas = new TCanvas("VtxRecon", "VtxRecon", 500, 500);
  InitPadCOLZ();
  corr->DrawCopy("COLZ");

  canvas->SaveAs(Form("VtxRecon_%d.eps", gMax));
  canvas->SaveAs(Form("VtxRecon_%d.eps", gMax));

  corr->GetYaxis()->SetRangeUser(0, 5);

  canvas = new TCanvas("VtxReconZoom", "VtxReconZoom", 500, 500);
  InitPadCOLZ();
  corr->DrawCopy("COLZ");

  canvas->SaveAs(Form("VtxReconZoom_%d.gif", gMax));
  canvas->SaveAs(Form("VtxReconZoom_%d.eps", gMax));
}

void VtxRecon1D(const char* fileName = "correction_map.root", const char* folderName = "dndeta_correction")
{
  gSystem->Load("libPWG0base");

  TFile* file = TFile::Open(fileName);
  AlidNdEtaCorrection* dNdEtaCorrection = new AlidNdEtaCorrection(folderName, folderName);
  dNdEtaCorrection->LoadHistograms();

  TH1* hist = dNdEtaCorrection->GetVertexRecoCorrection()->GetEventCorrection()->Get1DCorrection("x");
  TH1* hist2 = dNdEtaCorrection->GetVertexRecoCorrection()->GetEventCorrection()->Get1DCorrection("y", -10, 10);

  TCanvas* canvas = new TCanvas("VtxRecon1D", "VtxRecon1D", 1000, 500);
  canvas->Divide(2, 1);

  canvas->cd(1);
  InitPad();

  Prepare1DPlot(hist);
  hist->SetTitle("");
  hist->GetYaxis()->SetTitle("correction factor");
  hist->GetYaxis()->SetRangeUser(1, 1.8);
  hist->GetYaxis()->SetTitleOffset(1.6);
  hist->DrawCopy();

  canvas->cd(2);
  InitPad();

  Prepare1DPlot(hist2);
  hist2->SetTitle("");
  hist2->GetYaxis()->SetTitle("correction factor");
  hist2->GetXaxis()->SetRangeUser(0, 20);
  hist2->GetYaxis()->SetTitleOffset(1.6);
  hist2->GetXaxis()->SetTitle("multiplicity");
  hist2->Draw();

  TPaveText* pave = new TPaveText(0.6, 0.8, 0.8, 0.85, "NDC");
  pave->SetFillColor(0);
  pave->AddText("|z| < 10 cm");
  pave->Draw();

  canvas->SaveAs("VtxRecon1D.eps");

  Correction1DCreatePlots(fileName, folderName, 9.9, 2);

  TH1* corrX = dynamic_cast<TH1*> (gROOT->FindObject("generated_x_div_measured_x"));
  TH1* corrZ = dynamic_cast<TH1*> (gROOT->FindObject("generated_z_div_measured_z"));

  Prepare1DPlot(corrX);
  Prepare1DPlot(corrZ);

  corrX->GetYaxis()->SetTitleOffset(1.5);
  corrZ->GetYaxis()->SetTitleOffset(1.5);

  corrX->SetTitle("a) z projection");
  corrZ->SetTitle("b) p_{T} projection");

  corrX->GetYaxis()->SetTitle("correction factor");
  corrZ->GetYaxis()->SetTitle("correction factor");

  corrZ->GetXaxis()->SetRangeUser(0.11, 9.9);

  TString canvasName;
  canvasName.Form("VtxRecon1D_Track");
  TCanvas* canvas = new TCanvas(canvasName, canvasName, 800, 400);
  canvas->Divide(2, 1);

  canvas->cd(1);
  InitPad();
  corrX->DrawCopy();

  canvas->cd(2);
  InitPad();
  gPad->SetLogx();
  corrZ->Draw();

  canvas->SaveAs("VtxRecon1D_Track.eps");
  canvas->SaveAs("VtxRecon1D_Track.gif");
}

void Track2ParticleAsNumber(const char* fileName = "correction_map.root")
{
  gSystem->Load("libPWG0base");

  TFile::Open(fileName);
  AlidNdEtaCorrection* dNdEtaCorrection = new AlidNdEtaCorrection("dndeta_correction", "dndeta_correction");
  dNdEtaCorrection->LoadHistograms(fileName, "dndeta_correction");

  TH3F* gene = dNdEtaCorrection->GetTrack2ParticleCorrection()->GetGeneratedHistogram();
  TH3F* meas = dNdEtaCorrection->GetTrack2ParticleCorrection()->GetMeasuredHistogram();

  TH3F* gene = dNdEtaCorrection->GetTrack2ParticleCorrection()->GetGeneratedHistogram();
  TH3F* meas = dNdEtaCorrection->GetTrack2ParticleCorrection()->GetMeasuredHistogram();

  gene->GetYaxis()->SetRangeUser(-0.8, 0.8);
  meas->GetYaxis()->SetRangeUser(-0.8, 0.8);
  gene->GetXaxis()->SetRangeUser(-10, 10);
  meas->GetXaxis()->SetRangeUser(-10, 10);

  Float_t eff1 = gene->Integral() / meas->Integral();
  Float_t error1 = TMath::Sqrt(gene->Integral()) / meas->Integral();

  printf("Correction without pT cut: %f +- %f\n", eff1, error1);

  gene->GetZaxis()->SetRangeUser(0.3, 10);
  meas->GetZaxis()->SetRangeUser(0.3, 10);

  Float_t eff2 = gene->Integral() / meas->Integral();
  Float_t error2 = TMath::Sqrt(gene->Integral()) / meas->Integral();

  printf("Correction with pT cut: %f +- %f\n", eff2, error2);

  gene->GetZaxis()->SetRangeUser(0.3, 1);
  meas->GetZaxis()->SetRangeUser(0.3, 1);

  Float_t eff3 = gene->Integral() / meas->Integral();
  Float_t error3 = TMath::Sqrt(gene->Integral()) / meas->Integral();

  printf("Correction with 0.3 < pT < 0.5: %f +- %f\n", eff3, error3);
}

void Correction1DCreatePlots(const char* fileName = "correction_map.root", const char* folderName = "dndeta_correction", Float_t upperPtLimit = 9.9, Int_t correctionType = 0)
{
  TFile::Open(fileName);
  AlidNdEtaCorrection* dNdEtaCorrection = new AlidNdEtaCorrection(folderName, folderName);
  dNdEtaCorrection->LoadHistograms();

  TH3F* gene = dNdEtaCorrection->GetCorrection(correctionType)->GetTrackCorrection()->GetGeneratedHistogram();
  TH3F* meas = dNdEtaCorrection->GetCorrection(correctionType)->GetTrackCorrection()->GetMeasuredHistogram();

  gene->GetZaxis()->SetRangeUser(0.3, upperPtLimit);
  meas->GetZaxis()->SetRangeUser(0.3, upperPtLimit);
  gene->GetYaxis()->SetRangeUser(-0.8, 0.8);
  meas->GetYaxis()->SetRangeUser(-0.8, 0.8);
  AliPWG0Helper::CreateDividedProjections(gene, meas, "x", kFALSE);
  gene->GetYaxis()->SetRange(0, 0);
  meas->GetYaxis()->SetRange(0, 0);

  gene->GetXaxis()->SetRangeUser(-9.9, 9.9);
  meas->GetXaxis()->SetRangeUser(-9.9, 9.9);
  AliPWG0Helper::CreateDividedProjections(gene, meas, "y", kFALSE);
  gene->GetZaxis()->SetRange(0, 0);
  meas->GetZaxis()->SetRange(0, 0);

  gene->GetYaxis()->SetRangeUser(-0.8, 0.8);
  meas->GetYaxis()->SetRangeUser(-0.8, 0.8);
  AliPWG0Helper::CreateDividedProjections(gene, meas, "z", kFALSE);
}

TCanvas* Correction1D(Int_t correctionType = 0, const char* fileName = "correction_map.root", const char* folder = "dndeta_correction", Float_t upperPtLimit = 9.9)
{
  gSystem->Load("libPWG0base");

  Correction1DCreatePlots(fileName, folder, upperPtLimit, correctionType);

  TH1* corrX = dynamic_cast<TH1*> (gROOT->FindObject(Form("generated_x_div_measured_x", folder, folder)));
  TH1* corrY = dynamic_cast<TH1*> (gROOT->FindObject(Form("generated_y_div_measured_y", folder, folder)));
  TH1* corrZ = dynamic_cast<TH1*> (gROOT->FindObject(Form("generated_z_div_measured_z", folder, folder)));

  Prepare1DPlot(corrX);
  Prepare1DPlot(corrY);
  Prepare1DPlot(corrZ);

  corrX->SetTitle("a) z projection");
  corrY->SetTitle("b) #eta projection");
  corrZ->SetTitle("c) p_{T} projection");

  corrX->GetYaxis()->SetTitle("correction factor");
  corrY->GetYaxis()->SetTitle("correction factor");
  corrZ->GetYaxis()->SetTitle("correction factor");
  corrX->GetYaxis()->SetTitleOffset(1.7);
  corrY->GetYaxis()->SetTitleOffset(1.7);
  corrZ->GetYaxis()->SetTitleOffset(1.7);
  corrX->GetYaxis()->SetRangeUser(0.8, 1.5);
  corrY->GetYaxis()->SetRangeUser(0.8, 1.5);
  corrZ->GetYaxis()->SetRangeUser(0.8, 1.5);

  corrZ->GetXaxis()->SetRangeUser(0.11, upperPtLimit);

  TString canvasName;
  canvasName.Form(Form("Correction1D_%d_%s_%f", correctionType, fileName, upperPtLimit));
  TCanvas* canvas = new TCanvas(canvasName, canvasName, 1200, 400);
  canvas->Divide(3, 1);

  TLatex* Tl = new TLatex;
  Tl->SetTextSize(0.04);
  Tl->SetBit(TLatex::kTextNDC);

  canvas->cd(1);
  InitPad();
  corrX->DrawCopy();
  Tl->DrawLatex(0.6, 0.8, "0.3 < p_{T} < 10");
  Tl->DrawLatex(0.6, 0.75, "|#eta| < 0.8");

  canvas->cd(2);
  InitPad();
  corrY->Draw();
  Tl->DrawLatex(0.6, 0.8, "0.3 < p_{T} < 10");
  Tl->DrawLatex(0.6, 0.75, "|vtx-z| < 10");

  canvas->cd(3);
  InitPad();
  gPad->SetLogx();
  corrZ->Draw();
  Tl->DrawLatex(0.6, 0.8, "|vtx-z| < 10");
  Tl->DrawLatex(0.6, 0.75, "|#eta| < 0.8");

  return canvas;
}

void Track2Particle1D(const char* fileName = "correction_map.root", const char* folder = "dndeta_correction", Float_t upperPtLimit = 9.9)
{
  gSystem->Load("libPWG0base");

  Correction1DCreatePlots(fileName, folder, upperPtLimit, AlidNdEtaCorrection::kTrack2Particle);

  TH1* corrX = dynamic_cast<TH1*> (gROOT->FindObject(Form("generated_x_div_measured_x", folder, folder)));
  TH1* corrY = dynamic_cast<TH1*> (gROOT->FindObject(Form("generated_y_div_measured_y", folder, folder)));
  TH1* corrZ = dynamic_cast<TH1*> (gROOT->FindObject(Form("generated_z_div_measured_z", folder, folder)));

  Prepare1DPlot(corrX);
  Prepare1DPlot(corrY);
  Prepare1DPlot(corrZ);

  corrX->SetTitle("a) z projection");
  corrY->SetTitle("a) #eta projection");
  corrZ->SetTitle("b) p_{T} projection");

  corrY->GetYaxis()->SetTitle("correction factor");
  corrZ->GetYaxis()->SetTitle("correction factor");

  corrZ->GetXaxis()->SetRangeUser(0, upperPtLimit);

  TString canvasName;
  canvasName.Form("Track2Particle1D_%s", folder);
  TCanvas* canvas = new TCanvas(canvasName, canvasName, 1200, 400);
  canvas->Divide(3, 1);

  canvas->cd(1);
  InitPad();
  corrX->DrawCopy();

  canvas->cd(2);
  InitPad();
  corrY->Draw();

  canvas->cd(3);
  InitPad();
  corrZ->Draw();

  canvas->SaveAs(Form("Track2Particle1D_%s_%f.gif", fileName, upperPtLimit));
  canvas->SaveAs(Form("Track2Particle1D_%s_%f.eps", fileName, upperPtLimit));

  //TPaveText* pave = new TPaveText(-0.4, 1.35, 0.4, 1.45);

  canvasName.Form("Track2Particle1D_%s_etapt", folder);
  TCanvas* canvas = new TCanvas(canvasName, canvasName, 1000, 500);
  canvas->Divide(2, 1);

  canvas->cd(1);
  InitPad();
  corrY->GetXaxis()->SetRangeUser(-0.99, 0.99);
  corrY->GetYaxis()->SetRangeUser(1, 1.5);
  corrY->GetYaxis()->SetTitleOffset(1.5);
  corrY->DrawCopy();
  TPaveText* pave = new TPaveText(0.3, 0.7, 0.7, 0.8, "NDC");
  pave->AddText("|z| < 10 cm");
  pave->AddText("0.3 GeV/c < p_{T} < 10 GeV/c");
  pave->Draw();

  canvas->cd(2);
  InitPad();
  gPad->SetLogx();
  corrZ->GetYaxis()->SetRangeUser(1, 2.5);
  corrZ->GetXaxis()->SetRangeUser(0.101, upperPtLimit);
  corrZ->GetYaxis()->SetTitleOffset(1.5);
  corrZ->DrawCopy();
  pave = new TPaveText(0.5, 0.7, 0.8, 0.8, "NDC");
  pave->AddText("|z| < 10 cm");
  pave->AddText("|#eta| < 0.8");
  pave->Draw();

  canvas->SaveAs(Form("Track2Particle1D_etapt_%s_%f.eps", fileName, upperPtLimit));
  canvas->SaveAs(Form("Track2Particle1D_etapt_%s_%f.gif", fileName, upperPtLimit));
}

void CompareTrack2Particle1D(Float_t upperPtLimit = 9.9)
{
  gSystem->Load("libPWG0base");

  // particle type
  for (Int_t particle=0; particle<4; ++particle)
  {
    TString dirName;
    dirName.Form("correction_%d", particle);
    Track2Particle1DCreatePlots("systematics-detail-only-positive.root", dirName, upperPtLimit);

    TString tmpx, tmpy, tmpz;
    tmpx.Form("gene_%s_nTrackToNPart_x_div_meas_%s_nTrackToNPart_x", dirName.Data(), dirName.Data());
    tmpy.Form("gene_%s_nTrackToNPart_y_div_meas_%s_nTrackToNPart_y", dirName.Data(), dirName.Data());
    tmpz.Form("gene_%s_nTrackToNPart_z_div_meas_%s_nTrackToNPart_z", dirName.Data(), dirName.Data());

    TH1* posX = dynamic_cast<TH1*> (gROOT->FindObject(tmpx)->Clone("pos_x"));
    TH1* posY = dynamic_cast<TH1*> (gROOT->FindObject(tmpy)->Clone("pos_y"));
    TH1* posZ = dynamic_cast<TH1*> (gROOT->FindObject(tmpz)->Clone("pos_z"));

    Track2Particle1DCreatePlots("systematics-detail-only-negative.root", dirName, upperPtLimit);

    TH1* negX = dynamic_cast<TH1*> (gROOT->FindObject(tmpx)->Clone("neg_x"));
    TH1* negY = dynamic_cast<TH1*> (gROOT->FindObject(tmpy)->Clone("neg_y"));
    TH1* negZ = dynamic_cast<TH1*> (gROOT->FindObject(tmpz)->Clone("neg_z"));

    posX->Divide(negX);
    posY->Divide(negY);
    posZ->Divide(negZ);

    Prepare1DPlot(posX);
    Prepare1DPlot(posY);
    Prepare1DPlot(posZ);

    Float_t min = 0.8;
    Float_t max = 1.2;

    posX->SetMinimum(min);
    posX->SetMaximum(max);
    posY->SetMinimum(min);
    posY->SetMaximum(max);
    posZ->SetMinimum(min);
    posZ->SetMaximum(max);

    posZ->GetXaxis()->SetRangeUser(0, upperPtLimit);

    posX->GetYaxis()->SetTitleOffset(1.7);
    posX->GetYaxis()->SetTitle("C_{+} / C_{-}");
    posY->GetYaxis()->SetTitleOffset(1.7);
    posY->GetYaxis()->SetTitle("C_{+} / C_{-}");
    posZ->GetYaxis()->SetTitleOffset(1.7);
    posZ->GetYaxis()->SetTitle("C_{+} / C_{-}");

    posZ->GetXaxis()->SetRangeUser(0, 1);

    TString canvasName;
    canvasName.Form("PosNegRatios_%s_%f", ((particle == 0) ? "Pi" : ((particle == 1) ? "K" : ((particle == 2) ? "p" : "other"))), upperPtLimit);

    TCanvas* canvas = new TCanvas(canvasName, canvasName, 1200, 400);
    canvas->Divide(3, 1);

    canvas->cd(1);
    InitPad();
    posX->DrawCopy();

    canvas->cd(2);
    InitPad();
    posY->DrawCopy();

    canvas->cd(3);
    InitPad();
    posZ->DrawCopy();

    canvas->SaveAs(Form("%s.gif", canvas->GetName()));
    canvas->SaveAs(Form("%s.eps", canvas->GetName()));
  }
}

void CompareTrack2Particle1D(const char* file1, const char* file2, Float_t upperPtLimit = 9.9)
{
  loadlibs();

  const char* folderName = "dndeta_correction";

  c = new TCanvas("CompareTrack2Particle1D", "CompareTrack2Particle1D", 1200, 400);
  c->Divide(3, 1);

  for (Int_t fileId = 0; fileId < 2; fileId++)
  {
    const char* file = ((fileId == 0) ? file1 : file2);
    Correction1DCreatePlots(file, folderName, upperPtLimit, 1);

    TH1* corr[3];
    corr[0] = dynamic_cast<TH1*> (gROOT->FindObject("generated_x_div_measured_x"));
    corr[1] = dynamic_cast<TH1*> (gROOT->FindObject("generated_y_div_measured_y"));
    corr[2] = dynamic_cast<TH1*> (gROOT->FindObject("generated_z_div_measured_z"));
    /*corr[0] = dynamic_cast<TH1*> (gROOT->FindObject("generated_x"))->Clone(Form("hist_x_%d", fileId));
    corr[1] = dynamic_cast<TH1*> (gROOT->FindObject("generated_y"))->Clone(Form("hist_y_%d", fileId));
    corr[2] = dynamic_cast<TH1*> (gROOT->FindObject("generated_z"))->Clone(Form("hist_z_%d", fileId));*/

    for (Int_t i=0; i<3; i++)
    {
      c->cd(i+1);
      InitPad();
      corr[i]->GetYaxis()->SetRangeUser(0.8, 2);
      corr[i]->SetLineColor(fileId+1);
      corr[i]->DrawCopy((fileId == 0) ? "" : "SAME");
    }
  }

  return;

  c->SaveAs(Form("%s.gif", canvas->GetName()));
  c->SaveAs(Form("%s.eps", canvas->GetName()));
}

void Track2Particle2DCreatePlots(const char* fileName = "correction_map.root")
{
  TFile::Open(fileName);
  AlidNdEtaCorrection* dNdEtaCorrection = new AlidNdEtaCorrection("dndeta_correction", "dndeta_correction");
  dNdEtaCorrection->LoadHistograms();

  TH3F* gene = dNdEtaCorrection->GetTrack2ParticleCorrection()->GetTrackCorrection()->GetGeneratedHistogram();
  TH3F* meas = dNdEtaCorrection->GetTrack2ParticleCorrection()->GetTrackCorrection()->GetMeasuredHistogram();

  gene->GetZaxis()->SetRangeUser(0.3, 10);
  meas->GetZaxis()->SetRangeUser(0.3, 10);
  AliPWG0Helper::CreateDividedProjections(gene, meas, "yx");
  gene->GetZaxis()->SetRange(0, 0);
  meas->GetZaxis()->SetRange(0, 0);

  gene->GetYaxis()->SetRangeUser(-0.8, 0.8);
  meas->GetYaxis()->SetRangeUser(-0.8, 0.8);
  AliPWG0Helper::CreateDividedProjections(gene, meas, "zx");
  gene->GetYaxis()->SetRange(0, 0);
  meas->GetYaxis()->SetRange(0, 0);

  gene->GetXaxis()->SetRangeUser(-10, 10);
  meas->GetXaxis()->SetRangeUser(-10, 10);
  AliPWG0Helper::CreateDividedProjections(gene, meas, "zy");
  gene->GetXaxis()->SetRange(0, 0);
  meas->GetXaxis()->SetRange(0, 0);
}

TCanvas* Track2Particle2D(const char* fileName = "correction_map.root", const char* folder = "dndeta_correction")
{
  gSystem->Load("libPWG0base");

  Track2Particle2DCreatePlots(fileName);

  TH2* corrYX = dynamic_cast<TH2*> (gROOT->FindObject("generated_yx_div_measured_yx"));
  TH2* corrZX = dynamic_cast<TH2*> (gROOT->FindObject("generated_zx_div_measured_zx"));
  TH2* corrZY = dynamic_cast<TH2*> (gROOT->FindObject("generated_zy_div_measured_zy"));

  Prepare2DPlot(corrYX);
  Prepare2DPlot(corrZX);
  Prepare2DPlot(corrZY);

  const char* title = "";
  corrYX->SetTitle(title);
  corrZX->SetTitle(title);
  corrZY->SetTitle(title);

  TCanvas* canvas = new TCanvas("Track2Particle2D", "Track2Particle2D", 1200, 400);
  canvas->Divide(3, 1);

  canvas->cd(1);
  InitPadCOLZ();
  corrYX->Draw("COLZ");

  canvas->cd(2);
  InitPadCOLZ();
  corrZX->Draw("COLZ");

  canvas->cd(3);
  InitPadCOLZ();
  corrZY->Draw("COLZ");

  canvas->SaveAs(Form("corr_track2particle_%d.gif", gMax));
  canvas->SaveAs(Form("corr_track2particle_%d.eps", gMax));

  return canvas;
}

void CompareTrack2Particle2D()
{
  gSystem->Load("libPWG0base");

  Track2Particle2DCreatePlots("correction_maponly-positive.root");

  TH2* posYX = dynamic_cast<TH2*> (gROOT->FindObject("generated_yx_div_measured_yx")->Clone("pos_yx"));
  TH2* posZX = dynamic_cast<TH2*> (gROOT->FindObject("generated_zx_div_measured_zx")->Clone("pos_zx"));
  TH2* posZY = dynamic_cast<TH2*> (gROOT->FindObject("generated_zy_div_measured_zy")->Clone("pos_zy"));

  Track2Particle2DCreatePlots("correction_maponly-negative.root");

  TH2* negYX = dynamic_cast<TH2*> (gROOT->FindObject("generated_yx_div_measured_yx")->Clone("neg_yx"));
  TH2* negZX = dynamic_cast<TH2*> (gROOT->FindObject("generated_zx_div_measured_zx")->Clone("neg_zx"));
  TH2* negZY = dynamic_cast<TH2*> (gROOT->FindObject("generated_zy_div_measured_zy")->Clone("neg_zy"));

  posYX->Divide(negYX);
  posZX->Divide(negZX);
  posZY->Divide(negZY);

  Prepare2DPlot(posYX);
  Prepare2DPlot(posZX);
  Prepare2DPlot(posZY);

  Float_t min = 0.8;
  Float_t max = 1.2;

  posYX->SetMinimum(min);
  posYX->SetMaximum(max);
  posZX->SetMinimum(min);
  posZX->SetMaximum(max);
  posZY->SetMinimum(min);
  posZY->SetMaximum(max);

  TCanvas* canvas = new TCanvas("CompareTrack2Particle2D", "CompareTrack2Particle2D", 1200, 400);
  canvas->Divide(3, 1);

  canvas->cd(1);
  InitPadCOLZ();
  posYX->Draw("COLZ");

  canvas->cd(2);
  InitPadCOLZ();
  posZX->Draw("COLZ");

  canvas->cd(3);
  InitPadCOLZ();
  posZY->Draw("COLZ");

  canvas->SaveAs("CompareTrack2Particle2D.gif");
  canvas->SaveAs("CompareTrack2Particle2D.eps");
}

void Track2Particle3D()
{
  // get left margin proper

  TFile* file = TFile::Open("correction_map.root");

  TH3* corr = dynamic_cast<TH3*> (file->Get("dndeta_correction/corr_nTrackToNPart"));

  corr->SetTitle("Correction Factor");
  SetRanges(corr->GetZaxis());

  Prepare3DPlot(corr);

  TCanvas* canvas = new TCanvas("Track2Particle3D", "Track2Particle3D", 500, 500);
  canvas->SetTheta(29.428);
  canvas->SetPhi(16.5726);

  corr->Draw();

  canvas->SaveAs("Track2Particle3D.gif");
  canvas->SaveAs("Track2Particle3D.eps");
}

void Track2Particle3DAll()
{
  TFile* file = TFile::Open("correction_map.root");

  TH3* gene = dynamic_cast<TH3*> (file->Get("dndeta_correction/gene_nTrackToNPart"));
  TH3* meas = dynamic_cast<TH3*> (file->Get("dndeta_correction/meas_nTrackToNPart"));
  TH3* corr = dynamic_cast<TH3*> (file->Get("dndeta_correction/corr_nTrackToNPart"));

  gene->SetTitle("Generated Particles");
  meas->SetTitle("Measured Tracks");
  corr->SetTitle("Correction Factor");

  Prepare3DPlot(gene);
  Prepare3DPlot(meas);
  Prepare3DPlot(corr);

  TCanvas* canvas = new TCanvas("Track2Particle3DAll", "Track2Particle3DAll", 1200, 400);
  canvas->Divide(3, 1);

  canvas->cd(1);
  InitPad();
  gene->Draw();

  canvas->cd(2);
  meas->Draw();

  canvas->cd(3);
  corr->Draw();

  canvas->SaveAs("Track2Particle3DAll.gif");
  canvas->SaveAs("Track2Particle3DAll.eps");
}

void MultiplicityMC(Int_t xRangeMax = 50)
{
  TFile* file = TFile::Open("multiplicityMC.root");

  if (!file)
  {
    printf("multiplicityMC.root could not be opened.\n");
    return;
  }

  TH1F* fMultiplicityESD = dynamic_cast<TH1F*> (file->Get("fMultiplicityESD"));
  TH1F* fMultiplicityMC = dynamic_cast<TH1F*> (file->Get("fMultiplicityMC"));
  TH2F* fCorrelation = dynamic_cast<TH2F*> (file->Get("fCorrelation"));

  TH1F* correction = new TH1F("MultiplicityMC_correction", "MultiplicityMC_correction;Ntracks;Npart", 76, -0.5, 75.5);
  TH1F* correctionWidth = new TH1F("MultiplicityMC_correctionwidth", "MultiplicityMC_correctionwidth;Ntracks;Npart", 76, -0.5, 75.5);
  //fMultiplicityMC->GetNbinsX(), fMultiplicityMC->GetXaxis()->GetXmin(), fMultiplicityMC->GetXaxis()->GetXmax());
  for (Int_t i=1; i<=correction->GetNbinsX(); ++i)
  {
    TH1D* proj = fCorrelation->ProjectionX("_px", i, i+1);
    proj->Fit("gaus", "0");
    correction->SetBinContent(i, proj->GetFunction("gaus")->GetParameter(1));
    correctionWidth->SetBinContent(i, proj->GetFunction("gaus")->GetParameter(2));

    continue;

    // draw for debugging
    new TCanvas;
    proj->DrawCopy();
    proj->GetFunction("gaus")->DrawCopy("SAME");
  }

  TH1F* fMultiplicityESDCorrected = new TH1F("fMultiplicityESDCorrected", "fMultiplicityESDCorrected", 2010, -0.5, 200.5);

  for (Int_t i=1; i<=correction->GetNbinsX(); ++i)
  {
    Float_t mean = correction->GetBinContent(i);
    Float_t width = correctionWidth->GetBinContent(i);

    Int_t fillBegin = fMultiplicityESDCorrected->FindBin(mean - width * 3);
    Int_t fillEnd   = fMultiplicityESDCorrected->FindBin(mean + width * 3);
    printf("bin %d mean %f width %f, filling from %d to %d\n", i, mean, width, fillBegin, fillEnd);

    for (Int_t j=fillBegin; j <= fillEnd; ++j)
    {
      fMultiplicityESDCorrected->AddBinContent(j, TMath::Gaus(fMultiplicityESDCorrected->GetXaxis()->GetBinCenter(j), mean, width, kTRUE) * fMultiplicityESD->GetBinContent(i));
    }
  }

  TH1F* fMultiplicityESDCorrectedRebinned = dynamic_cast<TH1F*> (fMultiplicityESDCorrected->Clone("fMultiplicityESDCorrectedRebinned"));
  fMultiplicityESDCorrectedRebinned->Rebin(10);
  fMultiplicityESDCorrectedRebinned->Scale(0.1);

  TH1F* ratio = dynamic_cast<TH1F*> (fMultiplicityESD->Clone("multiplicity_ratio"));
  ratio->SetTitle("ratio;Ntracks;Nreco/Ngene");
  ratio->Divide(fMultiplicityMC);

  TH1F* ratio2 = dynamic_cast<TH1F*> (fMultiplicityESDCorrectedRebinned->Clone("multiplicity_ratio_corrected"));
  ratio2->Divide(fMultiplicityMC);

  TCanvas* canvas = new TCanvas("MultiplicityMC", "MultiplicityMC", 1500, 1000);
  canvas->Divide(3, 2);

  fMultiplicityESD->GetXaxis()->SetRangeUser(0, xRangeMax);
  ratio->GetXaxis()->SetRangeUser(0, xRangeMax);
  fCorrelation->GetXaxis()->SetRangeUser(0, xRangeMax);
  fCorrelation->GetYaxis()->SetRangeUser(0, xRangeMax);
  correction->GetXaxis()->SetRangeUser(0, xRangeMax);
  fMultiplicityESDCorrected->GetXaxis()->SetRangeUser(0, xRangeMax);
  fMultiplicityESDCorrectedRebinned->GetXaxis()->SetRangeUser(0, xRangeMax);

  canvas->cd(1); //InitPad();
  fMultiplicityESD->Draw();
  fMultiplicityMC->SetLineColor(2);
  fMultiplicityMC->Draw("SAME");

  TLegend* legend = new TLegend(0.6, 0.7, 0.85, 0.85);
  legend->AddEntry(fMultiplicityESD, "ESD");
  legend->AddEntry(fMultiplicityMC, "MC");
  legend->Draw();

  canvas->cd(2);
  fCorrelation->Draw("COLZ");

  canvas->cd(3);
  correction->Draw();
  //correction->Fit("pol1");
  correctionWidth->SetLineColor(2);
  correctionWidth->Draw("SAME");

  legend = new TLegend(0.2, 0.7, 0.45, 0.85);
  legend->AddEntry(correction, "#bar{x}");
  legend->AddEntry(correctionWidth, "#sigma");
  legend->Draw();

  canvas->cd(4);
  ratio->Draw();

  ratio2->SetLineColor(2);
  ratio2->Draw("SAME");

  legend = new TLegend(0.6, 0.7, 0.85, 0.85);
  legend->AddEntry(ratio, "uncorrected");
  legend->AddEntry(ratio2, "corrected");
  legend->Draw();

  canvas->cd(5);
  fMultiplicityESDCorrected->SetLineColor(kBlue);
  fMultiplicityESDCorrected->Draw();
  fMultiplicityMC->Draw("SAME");
  fMultiplicityESD->Draw("SAME");

  legend = new TLegend(0.6, 0.7, 0.85, 0.85);
  legend->AddEntry(fMultiplicityESDCorrected, "ESD corrected");
  legend->AddEntry(fMultiplicityMC, "MC");
  legend->AddEntry(fMultiplicityESD, "ESD");
  legend->Draw();

  canvas->cd(6);
  fMultiplicityESDCorrectedRebinned->SetLineColor(kBlue);
  fMultiplicityESDCorrectedRebinned->Draw();
  fMultiplicityMC->Draw("SAME");

  legend = new TLegend(0.6, 0.7, 0.85, 0.85);
  legend->AddEntry(fMultiplicityESDCorrectedRebinned, "ESD corrected");
  legend->AddEntry(fMultiplicityMC, "MC");
  legend->Draw();

  canvas->SaveAs("MultiplicityMC.gif");
}

void MultiplicityESD()
{
  TFile* file = TFile::Open("multiplicityESD.root");

  if (!file)
  {
    printf("multiplicityESD.root could not be opened.\n");
    return;
  }

  TH1F* fMultiplicityESD = dynamic_cast<TH1F*> (file->Get("fMultiplicity"));

  TCanvas* canvas = new TCanvas("MultiplicityESD", "MultiplicityESD", 500, 500);

  fMultiplicityESD->Draw();
}

void drawPlots(Int_t max)
{
  gMax = max;

  ptCutoff();
  TriggerBias();
  VtxRecon();
  Track2Particle2D();
  Track2Particle3D();
  ptSpectrum();
  dNdEta();
}

void drawPlots()
{
  drawPlots(5);
  drawPlots(2);
}

void CompareCorrection2Measured(const char* dataInput = "analysis_esd_raw.root", const char* correctionMapFile = "correction_map.root", const char* correctionMapFolder = "dndeta_correction")
{
  loadlibs();

  AlidNdEtaCorrection* dNdEtaCorrection = new AlidNdEtaCorrection(correctionMapFolder, correctionMapFolder);
  TFile::Open(correctionMapFile);
  dNdEtaCorrection->LoadHistograms();

  TFile* file = TFile::Open(dataInput);

  if (!file)
  {
    cout << "Error. File not found" << endl;
    return;
  }

  dNdEtaAnalysis* fdNdEtaAnalysis = new dNdEtaAnalysis("dndeta", "dndeta");
  fdNdEtaAnalysis->LoadHistograms("fdNdEtaAnalysisESD");

  gROOT->cd();
  
  TH3* hist1 = (TH3*) dNdEtaCorrection->GetTrack2ParticleCorrection()->GetTrackCorrection()->GetMeasuredHistogram()->Clone("mc");
  hist1->SetTitle("mc");
  Printf("mc contains %f entries", hist1->Integral());
  Printf("mc contains %f entries in |vtx-z| < 10, pt > 0.3", hist1->Integral(hist1->GetXaxis()->FindBin(-9.9), hist1->GetXaxis()->FindBin(9.9), 1, hist1->GetNbinsY(), hist1->GetZaxis()->FindBin(0.301), hist1->GetNbinsZ()));

  TH3* hist2 = (TH3*) fdNdEtaAnalysis->GetData()->GetTrackCorrection()->GetMeasuredHistogram()->Clone("esd");
  hist2->SetTitle("esd");
  Printf("esd contains %f entries", hist2->Integral());
  Printf("esd contains %f entries in |vtx-z| < 10, pt > 0.3", hist2->Integral(hist2->GetXaxis()->FindBin(-9.9), hist2->GetXaxis()->FindBin(9.9), 1, hist2->GetNbinsY(), hist2->GetZaxis()->FindBin(0.301), hist2->GetNbinsZ()));

  AliPWG0Helper::CreateDividedProjections(hist1, hist2);
  AliPWG0Helper::CreateDividedProjections(hist1, hist2, "x");

  hist1->GetXaxis()->SetRange(hist1->GetXaxis()->FindBin(-10), hist2->GetXaxis()->FindBin(10));
  hist2->GetXaxis()->SetRange(hist1->GetXaxis()->FindBin(-10), hist2->GetXaxis()->FindBin(10));
  AliPWG0Helper::CreateDividedProjections(hist1, hist2, "y");

  new TCanvas; gROOT->FindObject("mc_yx_div_esd_yx")->Draw("COLZ");
  new TCanvas; gROOT->FindObject("mc_zx_div_esd_zx")->Draw("COLZ");
  new TCanvas; gROOT->FindObject("mc_zy_div_esd_zy")->Draw("COLZ");
  new TCanvas; gROOT->FindObject("mc_x_div_esd_x")->Draw("COLZ");
  new TCanvas; gROOT->FindObject("mc_y_div_esd_y")->Draw("COLZ");
}

void CompareMeasured2Measured(const char* dataInput = "analysis_esd_raw.root", const char* dataInput2 = "analysis_esd_raw.root")
{
  loadlibs();

  TFile* file = TFile::Open(dataInput);

  if (!file)
  {
    cout << "Error. File not found" << endl;
    return;
  }

  dNdEtaAnalysis* fdNdEtaAnalysis = new dNdEtaAnalysis("dndeta", "dndeta");
  fdNdEtaAnalysis->LoadHistograms("fdNdEtaAnalysisESD");

  TFile* file = TFile::Open(dataInput2);

  if (!file)
  {
    cout << "Error. File not found" << endl;
    return;
  }

  dNdEtaAnalysis* fdNdEtaAnalysis2 = new dNdEtaAnalysis("dndeta2", "dndeta2");
  fdNdEtaAnalysis2->LoadHistograms("fdNdEtaAnalysisESD");

  gROOT->cd();

  TH3* hist1 = (TH3*) fdNdEtaAnalysis->GetData()->GetTrackCorrection()->GetMeasuredHistogram()->Clone("esd1");
  hist1->SetTitle("esd1");
  Printf("esd1 contains %f entries", hist1->GetEntries());
  Printf("esd1 contains %f entries in |vtx-z| < 10, pt > 0.3", hist1->Integral(hist1->GetXaxis()->FindBin(-9.9), hist1->GetXaxis()->FindBin(9.9), 1, hist1->GetNbinsY(), hist1->GetZaxis()->FindBin(0.301), hist1->GetNbinsZ()));

  TH3* hist2 = (TH3*) fdNdEtaAnalysis2->GetData()->GetTrackCorrection()->GetMeasuredHistogram()->Clone("esd2");
  hist2->SetTitle("esd2");
  Printf("esd2 contains %f entries", hist2->GetEntries());
  Printf("esd2 contains %f entries in |vtx-z| < 10, pt > 0.3", hist2->Integral(hist2->GetXaxis()->FindBin(-9.9), hist2->GetXaxis()->FindBin(9.9), 1, hist2->GetNbinsY(), hist2->GetZaxis()->FindBin(0.301), hist2->GetNbinsZ()));

  AliPWG0Helper::CreateDividedProjections(hist1, hist2);

  new TCanvas; gROOT->FindObject("esd1_yx_div_esd2_yx")->Draw("COLZ");
  new TCanvas; gROOT->FindObject("esd1_zx_div_esd2_zx")->Draw("COLZ");
  new TCanvas; gROOT->FindObject("esd1_zy_div_esd2_zy")->Draw("COLZ");

  TH2* event1 = (TH2*) fdNdEtaAnalysis->GetData()->GetEventCorrection()->GetMeasuredHistogram()->Clone("event1");
  TH2* event2 = (TH2*) fdNdEtaAnalysis2->GetData()->GetEventCorrection()->GetMeasuredHistogram()->Clone("event2");

  Printf("event1 contains %f entries", event1->GetEntries());
  Printf("event2 contains %f entries", event2->GetEntries());
  Printf("event1 integral is %f", event1->Integral());
  Printf("event2 integral is %f", event2->Integral());
  Printf("event1 contains %f entries in |vtx-z| < 10", event1->Integral(event1->GetXaxis()->FindBin(-9.9), event1->GetXaxis()->FindBin(9.9), 1, event1->GetNbinsY()));
  Printf("event2 contains %f entries in |vtx-z| < 10", event2->Integral(event2->GetXaxis()->FindBin(-9.9), event2->GetXaxis()->FindBin(9.9), 1, event2->GetNbinsY()));

  projx1 = event1->ProjectionX();
  projx2 = event2->ProjectionX();

  new TCanvas; projx1->DrawCopy(); projx2->SetLineColor(2); projx2->DrawCopy("SAME");

  projx1->Divide(projx2);
  new TCanvas; projx1->Draw();

  event1->Divide(event2);
  new TCanvas; event1->Draw("COLZ");

}

void DrawTrackletOrigin()
{
  TFile::Open("correction_map.root");

  Int_t colors[]  = {1,2,3,4,6,7,8,102};

  Int_t maxHists = 8;
  TH1* hist[8];

  const char* titles[] = { "PP", "SS", "PP'", "PS", "PS*", "SP", "SS'", "" };

  TLegend* legend = new TLegend(0.75, 0.6, 0.95, 0.95);

  Int_t total = 0;
  for (Int_t i=0; i<maxHists; i++)
  {
    hist[i] = (TH1*) gFile->Get(Form("fDeltaPhi_%d", i));
    //hist[i]->Rebin(20);
    hist[i]->SetStats(kFALSE);
    hist[i]->SetLineColor(colors[i]);
    hist[i]->GetXaxis()->SetRangeUser(-0.2, 0.2);
    hist[i]->Draw(((i == 0) ? "" : "SAME"));

    total += hist[i]->GetEntries();

    if (i != 7)
      legend->AddEntry(hist[i], titles[i]);
  }

  legend->Draw();
  gPad->SetLogy();

  Printf("Total: %d", total);
  for (Int_t i=0; i<maxHists; i++)
    Printf("Histogram %d (%s) containts %.2f %% of the entries", i, titles[i], 100.0 * hist[i]->GetEntries() / total);

  printf("|  Delta phi  |  Acc. %%  |  ");
  for (Int_t i=0; i<maxHists; i++)
    printf("%3s %%   |  ", titles[i]);
  Printf("");

  for (Float_t f = 0.01; f < 0.09; f += 0.01)
  {
    Int_t integralBegin = hist[0]->GetXaxis()->FindBin(-f);
    Int_t integralEnd = hist[0]->GetXaxis()->FindBin(f);

    Int_t total2 = 0;
    for (Int_t i=0; i<maxHists; i++)
      total2 += (Int_t) hist[i]->Integral(integralBegin, integralEnd);

    printf("|    %.2f     |  %6.2f  |  ", f, 100.0 * total2 / total);

    for (Int_t i=0; i<maxHists; i++)
      printf("%6.2f  |  ", (hist[i]->GetEntries() > 0) ? (100.0 * hist[i]->Integral(integralBegin, integralEnd) / hist[i]->GetEntries()) : -1.0);
    Printf("");
  }
}

TH2* GetCorrection(const char* fileName = "correction_map.root", const char* dirName = "dndeta_correction", Double_t ptmin=0.2)
{
  // returns the correction factor with pt integrated out

  loadlibs();

  TFile::Open(fileName);

  AlidNdEtaCorrection* dNdEtaCorrection = new AlidNdEtaCorrection(dirName, dirName);
  if (!dNdEtaCorrection->LoadHistograms())
    return;

  //  hist = dNdEtaCorrection->GetTrack2ParticleCorrection()->GetTrackCorrection()->GetCorrectionHistogram();

  gener = dNdEtaCorrection->GetTrack2ParticleCorrection()->GetTrackCorrection()->GetGeneratedHistogram();
  measu = dNdEtaCorrection->GetTrack2ParticleCorrection()->GetTrackCorrection()->GetMeasuredHistogram();

  gener->GetZaxis()->SetRange(gener->GetZaxis()->FindBin(ptmin), gener->GetNbinsZ()+1);
  TH2D *gener_xy = gener->Project3D("yx");

  measu->GetZaxis()->SetRange(measu->GetZaxis()->FindBin(ptmin), measu->GetNbinsZ()+1);
  TH2D *measu_xy = measu->Project3D("yx");

  cout << measu->GetZaxis()->FindBin(ptmin) << " " << measu->GetNbinsZ()+1 << endl;

  TCanvas *canp = new TCanvas("canp","canp",600,1000);
  canp->Divide(1,2,0.0001,0.0001);
  canp->cd(1);
  gener_xy->Draw("COLZ");
  canp->cd(2);
  measu_xy->Draw("COLZ");


  TCanvas *canpr = new TCanvas("canpr","canpr",700,500);
  canpr->cd();
  TH2D *proj = new TH2D(*gener_xy);
  proj->Divide(measu_xy);

//   proj = hist->Project3D("yx");
  proj->Draw("COLZ");

  return proj;
}

void DetermineAcceptance(const char* fileName = "correction_map.root", const char* dirName = "dndeta_correction", Double_t ptmin=0.2)
{
  TH2* proj = GetCorrection(fileName, dirName, ptmin);

  const Float_t limit = 5;

  TString array = "{";
  TString arrayEnd = "}";

  for (Int_t y=1; y<=proj->GetNbinsY(); ++y)
  {
    Int_t begin = -1;
    Int_t end = -1;
    for (Int_t x=1; x<=proj->GetNbinsX(); ++x)
    {
      if (begin == -1 && proj->GetBinContent(x, y) > 0 && proj->GetBinContent(x, y) < limit)
        begin = x;
      if (begin != -1 && proj->GetBinContent(x, y) > 0 && proj->GetBinContent(x, y) < limit)
        end = x;
    }
    Printf("Limits for y = %d are %d to %d", y, begin, end);

    if (y > 1)
      array += ", ";
    array += Form("%d", begin);

    if (y > 1)
      arrayEnd.Prepend(", ");
    arrayEnd.Prepend(Form("%d", (end == -1) ? -1 : proj->GetNbinsX() + 1 - end));
  }
  array += "}";
  arrayEnd.Prepend("{");

  Printf("Begin array:");
  Printf("%s", array.Data());

  Printf("End array (mirrored) (should be the same):");
  Printf("%s", arrayEnd.Data());
}

void AverageMultiplicity(const char* fileName = "correction_map.root", const char* correctionMapFolder = "dndeta_correction")
{
  loadlibs();

  TFile::Open(fileName);

  AlidNdEtaCorrection* dNdEtaCorrection = new AlidNdEtaCorrection(correctionMapFolder, correctionMapFolder);
  dNdEtaCorrection->LoadHistograms();
  TH2* events = dNdEtaCorrection->GetTriggerBiasCorrectionINEL()->GetEventCorrection()->GetGeneratedHistogram();
  TH3* tracks = dNdEtaCorrection->GetTriggerBiasCorrectionINEL()->GetTrackCorrection()->GetGeneratedHistogram();

  Float_t nEvents = events->Integral(events->GetXaxis()->FindBin(-1), events->GetXaxis()->FindBin(1), 0, events->GetNbinsY()+1);
  Float_t nTracks = tracks->Integral(tracks->GetXaxis()->FindBin(-1), tracks->GetXaxis()->FindBin(1), tracks->GetYaxis()->FindBin(-0.39), tracks->GetYaxis()->FindBin(0.59), 0, tracks->GetNbinsZ()+1);

  Printf("%f %f --> %f", nEvents, nTracks, nTracks / nEvents);
}

void GetAverageCorrectionFactor(Float_t etaRange = 1.5, Float_t vertexRange = 9.9, const char* rawFile = "analysis_esd_raw.root", const char* mcFile = "analysis_mc.root")
{
  loadlibs();

  TFile::Open(rawFile);
  dNdEtaAnalysis* raw = new dNdEtaAnalysis("dndeta", "dndeta");
  raw->LoadHistograms("fdNdEtaAnalysisESD");
  raw->GetData()->GetTrackCorrection()->GetMeasuredHistogram()->GetXaxis()->SetRangeUser(-vertexRange, vertexRange);
  tracks = raw->GetData()->GetTrackCorrection()->GetMeasuredHistogram()->Project3D("y");
  events = raw->GetData()->GetEventCorrection()->GetMeasuredHistogram()->ProjectionX("events", 0, raw->GetData()->GetEventCorrection()->GetMeasuredHistogram()->GetNbinsY() + 1);
  Float_t nEvents = events->Integral(events->FindBin(-vertexRange), events->FindBin(vertexRange));
  tracks->Scale(1.0 / nEvents / tracks->GetBinWidth(1));

  TFile::Open(mcFile);
  dNdEtaAnalysis* mc = new dNdEtaAnalysis("dndeta", "dndeta");
  mc->LoadHistograms("dndetaTrVtx");
  mcH = mc->GetdNdEtaPtCutOffCorrectedHistogram(0);

  new TCanvas;
  mcH->SetLineColor(2);
  mcH->DrawCopy();
  tracks->DrawCopy("SAME");

  new TCanvas;
  mcH->GetYaxis()->SetRangeUser(0, 5);
  mcH->Divide(tracks);
  mcH->DrawCopy();
  mcH->Fit("pol0", "", "", -etaRange, etaRange);
}

void TrackCuts_Comparison(char* histName, Int_t plotWhich = 0, const char* fileName = "correction_map.root")
{
  // for the nsigmaplot it is needed to run with all cuts except the nsigmatovertex
  //    --> manually disable it in the run.C
  //
  // plotWhich: 0 = only before
  //            1 = both
  //            2 = only after

  file = TFile::Open(fileName);

  Int_t count = 0;
  Int_t colors[] = { 1, 2, 3, 4, 5, 6 };

  TLegend* legend = new TLegend(0.5, 0.7, 1, 1);
  TLegend* legend2 = new TLegend(0.4, 0.6, 1, 1);
  TLegend* legend3 = new TLegend(0.6, 0.5, 1, 0.7);

  TCanvas* c1 = new TCanvas("c1", "c1", 800, 1200);
  c1->Divide(1, 2);
  //TCanvas* c2 = new TCanvas("c2", "c2", 800, 600);
  //TCanvas* c3 = new TCanvas("c3", "c3", 800, 600);

  const char* folders2[] = { "before_cuts", "after_cuts" };
  Bool_t first = kTRUE;
  for (Int_t j = ((plotWhich < 2) ? 0 : 1); j < ((plotWhich > 0) ? 2 : 1); j++)
  {
    const char* folders1[] = { "esd_track_cuts", "esd_track_cuts_primaries", "esd_track_cuts_secondaries" };
    const char* names[] =    { "all", "primaries", "secondaries" };
    TH1* base = 0;
    TH1* prim = 0;
    TH1* sec = 0;
    for (Int_t i = 0; i < 3; i++)
    {
      TString folder;
      folder.Form("%s/%s/%s", folders1[i], folders2[j], histName);
      TH1* hist = (TH1*) file->Get(folder);
      legend->AddEntry(hist, Form("%s %s", names[i], folders2[j]));

      c1->cd(1);
      hist->SetLineColor(colors[count]);
      hist->DrawCopy((count == 0) ? "" : "SAME");

      switch (i)
      {
        case 0: base = hist; break;
        case 1: prim = hist; break;
        case 2: sec = hist; break;
      }

      count++;
    }

    TH1* eff    = (TH1*) prim->Clone("eff"); eff->Reset();
    TH1* purity = (TH1*) prim->Clone("purity"); purity->Reset();

    for (Int_t bin = 1; bin <= prim->GetNbinsX(); bin++)
    {
      eff->SetBinContent(bin, prim->Integral(1, bin) / prim->Integral(1, prim->GetNbinsX() + 1));
      if (prim->Integral(1, bin) + sec->Integral(1, bin) > 0)
        purity->SetBinContent(bin, sec->Integral(1, bin) / (prim->Integral(1, bin) + sec->Integral(1, bin)));
    }

    eff->GetYaxis()->SetRangeUser(0, 1);
    eff->SetLineColor(colors[0+j*2]);
    eff->SetStats(kFALSE);
    purity->SetLineColor(colors[1+j*2]);

    legend3->AddEntry(eff, Form("%s: efficiency", folders2[j]));
    legend3->AddEntry(purity, Form("%s: contamination", folders2[j]));

    c1->cd(2);
    eff->DrawCopy((first) ? "" : "SAME");
    first = kFALSE;
    purity->DrawCopy("SAME");
  }

  c1->cd(1)->SetLogy();
  c1->cd(1)->SetGridx();
  c1->cd(1)->SetGridy();
  legend->Draw();

  //c2->cd();
 // c2->SetGridx();
 // c2->SetGridy();
  //legend2->Draw();

  c1->cd(2)->SetGridx();
  c1->cd(2)->SetGridy();
  legend3->Draw();

  c1->SaveAs(Form("%s.png", histName));
}

void TrackCuts_DCA()
{
  file = TFile::Open("correction_map.root");
  hist = (TH2*) file->Get("esd_track_cuts/before_cuts/dXYvsDZ");

  TCanvas* c1 = new TCanvas("c1", "c1", 600, 600);
  c1->SetLogz();
  c1->SetRightMargin(0.12);
  c1->SetBottomMargin(0.12);

  hist->SetStats(kFALSE);
  hist->Draw("COLZ");

  ellipse = new TEllipse(0, 0, 4);
  ellipse->SetLineWidth(2);
  ellipse->SetLineStyle(2);
  ellipse->SetFillStyle(0);
  ellipse->Draw();

  c1->SaveAs("trackcuts_dca_2d.eps");
}

void FindNSigma(TH2* hist, Int_t nSigma = 1)
{
  TH1* proj = hist->ProjectionY();
  proj->Reset();

  for (Int_t bin=1; bin<=proj->GetNbinsX(); bin++)
  {
    if (hist->Integral(1, hist->GetNbinsX(), bin, bin) == 0)
      continue;

    Int_t limit = -1;
    for (limit = 1; limit<=hist->GetNbinsX(); limit++)
  }
}

void ShowOnlyAccepted(TH2* input, const char* fileName = "correction_map.root", const char* dirName = "dndeta_correction", Double_t ptmin=0.2)
{
  TH2* proj = GetCorrection(fileName, dirName, ptmin);

  for (Int_t y=1; y<=proj->GetNbinsY(); ++y)
    for (Int_t x=1; x<=proj->GetNbinsX(); ++x)
      if (proj->GetBinContent(x, y) > 5 || proj->GetBinContent(x, y) == 0)
      {
        proj->SetBinContent(x, y, 0);
      }
      else
        proj->SetBinContent(x, y, 1);


  input->Multiply(proj);
}

void MakeGaussianProfile(const char* histName = "fVertexCorrelation", Bool_t subtractMean = kFALSE)
{
    TFile::Open("correction_map.root");

    TH2* hist2d = (TH2*) gFile->Get(histName);
    hist2d->Sumw2();

    TH1* result = hist2d->ProjectionX("result");
    result->GetYaxis()->SetTitle(hist2d->GetYaxis()->GetTitle());
    result->Reset();

    for (Int_t x=1; x<hist2d->GetNbinsX(); ++x)
    {
        hist = hist2d->ProjectionY(Form("temp_%d", x), x, x);
        if (hist->GetEntries() == 0)
            continue;
        if (hist->Fit("gaus") == 0)
        {
            func = hist->GetFunction("gaus");
            mean = func->GetParameter(1);
            error = func->GetParError(1);

            if (subtractMean)
                mean = hist2d->GetXaxis()->GetBinCenter(x) - mean;

            result->SetBinContent(x, mean);
            result->SetBinError(x, error);

            if (x % 10 == 0)
            {
                new TCanvas;
                ((TH1*) hist->Clone())->DrawCopy();
            }
        }
        //break;
    }

    new TCanvas;
    result->GetYaxis()->SetRangeUser(-0.2, 0.2);
    result->Draw();
}
