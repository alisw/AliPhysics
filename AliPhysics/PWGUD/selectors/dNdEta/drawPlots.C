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
  gSystem->Load("libVMC");
  gSystem->Load("libTree");
  gSystem->Load("libSTEERBase");
  gSystem->Load("libESD");
  gSystem->Load("libAOD");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libPWG0base");
}

void SetRanges(TAxis* axis)
{
  if (strcmp(axis->GetTitle(), "#eta") == 0)
    axis->SetRangeUser(-1.4999, 1.4999);
  if (strcmp(axis->GetTitle(), "p_{T} [GeV/c]") == 0 || strcmp(axis->GetTitle(), "p_{T} (GeV/c)") == 0)
  {
    axis->SetRangeUser(0, 4.9999);
    axis->SetTitle("p_{T} (GeV/c)");
  }
  if (strcmp(axis->GetTitle(), "vtx z [cm]") == 0 || strcmp(axis->GetTitle(), "vtx z (cm)") == 0 || strcmp(axis->GetTitle(), "vtx-z [cm]") == 0 || strcmp(axis->GetTitle(), "vtx-z (cm)") == 0)
  {
    axis->SetRangeUser(-15, 14.9999);
    axis->SetTitle("vtx-z (cm)");
  }
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
    dNdEtaCorrection->GetCorrection(i)->PrintInfo(0.2);
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
  loadlibs();

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
  loadlibs();
  
  dNdEtaAnalysis* fdNdEtaAnalysis = new dNdEtaAnalysis(folder, folder);
  fdNdEtaAnalysis->LoadHistograms();
  fdNdEtaAnalysis->Finish(0, ptCut, AlidNdEtaCorrection::kNone, tag);
  return fdNdEtaAnalysis->GetdNdEtaHistogram(0);
}

void dNdEtaFinal(Bool_t spd = kTRUE)
{
  TFile* file = TFile::Open("analysis_esd.root");
  TH1* histESD = (TH1*) file->Get("dndeta/dNdEta_corrected");
  TH1* histESDnsd = (TH1*) file->Get("dndetaNSD/dNdEta_corrected");
  Prepare1DPlot(histESD);
  Prepare1DPlot(histESDnsd);
  
  TCanvas* canvas = new TCanvas("dNdEtaFinal", "dNdEtaFinal", 600, 600);
  gPad->SetTopMargin(0.05);
  gPad->SetRightMargin(0.05);
  gPad->SetLeftMargin(0.12);
  gPad->SetBottomMargin(0.12);
  gPad->SetGridx();
  gPad->SetGridy();
  
  Float_t etaMax = 1.9;
  Float_t histMax = 1.39;
  Float_t systErrorValue = 0.023;
  Float_t systErrorNSDValue = 0.081;
  if (!spd)
  {
    //etaMax = 1.5;
    histMax = 0.99;
    systErrorValue = 0.043;
    systErrorNSDValue = 0.088;
  }
  
  dummy = new TH2F("dummy", ";#eta;dN_{ch}/d#eta", 100, -etaMax, etaMax, 100, 3, 8);
  dummy->SetStats(0);
  dummy->GetYaxis()->SetTitleOffset(1.3);
  
  histESD->SetMarkerStyle(20);
  histESDnsd->SetMarkerStyle(21);
  histESDnsd->SetMarkerColor(4);
  histESDnsd->SetLineColor(4);
  histESD->SetMarkerSize(1.5);
  histESDnsd->SetMarkerSize(1.5);
  
  histESD->GetXaxis()->SetRangeUser(-histMax, histMax);
  histESDnsd->GetXaxis()->SetRangeUser(-histMax, histMax);
  
  legend = new TLegend(0.3, 0.2, 0.78, 0.4);
  legend->SetFillColor(0);
  legend->SetTextSize(0.04);
  legend->AddEntry(histESD, "Inelastic events", "P");
  legend->AddEntry(histESDnsd, "NSD events", "P");
  
  dummy->Draw();
  
  // syst errors.
  TH1* systError = (TH1*) histESD->Clone("systError");
  for (Int_t i=1; i<=systError->GetNbinsX(); ++i)
    systError->SetBinError(i, systError->GetBinContent(i) * systErrorValue);
  // change error drawing style
  systError->SetFillColor(15);    
  systError->DrawCopy("SAME E2 ][");
  
  // syst error NSD
  for (Int_t i=1; i<=systError->GetNbinsX(); ++i)
  {
    systError->SetBinContent(i, histESDnsd->GetBinContent(i));
    systError->SetBinError(i, systError->GetBinContent(i) * systErrorNSDValue);
  }
  systError->DrawCopy("SAME E2 ][");
  
  histESD->Draw("SAME");
  histESDnsd->Draw("SAME");
  legend->Draw();  
  
  canvas->SaveAs(Form("%s_dndeta_final.eps", (spd) ? "spd" : "tpc"));
}

void dNdEtaPythiaPhojet()
{
  // evtl. deactivate acceptance maps in dNdEtaAnalysis.cxx

  loadlibs();

  TH1* hist[4];
  
  TFile::Open("LHC08c11_10TeV_0.5T/mb1/spd/analysis_mc.root");
  hist[0] =         (TH1*) GetMCHist("dndeta", -1, "MC: full inelastic")->Clone("histMC");
  hist[1] =         (TH1*) GetMCHist("dndetaNSD", -1, "MC: NSD")->Clone("histMCnsd");

  TFile::Open("LHC08c15_10TeV_0.5T_Phojet/mb1/spd/analysis_mc.root");
  hist[2] =         (TH1*) GetMCHist("dndeta", -1, "MC: full inelastic")->Clone("histMCPhojet");
  hist[3] =         (TH1*) GetMCHist("dndetaNSD", -1, "MC: NSD")->Clone("histMCnsdPhojet");
  
  file = TFile::Open("pythia_phojet_dndeta.root", "RECREATE");
  for (Int_t i=0; i<4; i++)
    hist[i]->Write();
  file->Close();
}
 
void dNdEta(Bool_t onlyESD = kFALSE, Bool_t save = kTRUE)
{
  loadlibs();

  TFile* file = TFile::Open("analysis_esd.root");
  
  dNdEtaAnalysis* fdNdEtaAnalysis = new dNdEtaAnalysis("dndeta", "dndeta");
  fdNdEtaAnalysis->LoadHistograms("dndeta");
  
  TH1* histESD = (TH1*) file->Get("dndeta/dNdEta_corrected");
  TH1* histESD1 = (TH1*) file->Get("dndeta/dNdEta_corrected_1");
  TH1* histESD2 = (TH1*) file->Get("dndeta/dNdEta_corrected_2");
  TH1* histESDnsd = (TH1*) file->Get("dndetaNSD/dNdEta_corrected");
  TH1* histESDnsdNoPt = (TH1*) file->Get("dndetaNSD/dNdEta");
  TH1* histESDonePart = (TH1*) file->Get("dndetaOnePart/dNdEta_corrected");
  if (!histESDonePart)
    histESDonePart = new TH1F;
  TH1* histESDNoPt = (TH1*) file->Get("dndeta/dNdEta");
  TH1* histESDMB = (TH1*) file->Get("dndetaTr/dNdEta_corrected");
  TH1* histESDMBNoPt = (TH1*) file->Get("dndetaTr/dNdEta");
  TH1* histESDMBVtx = (TH1*) file->Get("dndetaTrVtx/dNdEta_corrected");
  TH1* histESDMBVtxNoPt = (TH1*) file->Get("dndetaTrVtx/dNdEta");
  TH1* histESDMBTracksNoPt = (TH1*) file->Get("dndetaTracks/dNdEta");

  Prepare1DPlot(histESD);
  Prepare1DPlot(histESD1);
  Prepare1DPlot(histESD2);
  Prepare1DPlot(histESDnsd);
  Prepare1DPlot(histESDonePart);
  Prepare1DPlot(histESDMB);
  Prepare1DPlot(histESDMBVtx);

  Prepare1DPlot(histESDNoPt);
  Prepare1DPlot(histESDMBNoPt);
  Prepare1DPlot(histESDMBVtxNoPt);
  Prepare1DPlot(histESDMBTracksNoPt);

  histESD->SetLineWidth(0);
  histESDnsd->SetLineWidth(0);
  histESDonePart->SetLineWidth(0);
  histESDMB->SetLineWidth(0);
  histESDMBVtx->SetLineWidth(0);

  histESDNoPt->SetLineWidth(0);
  histESDMBNoPt->SetLineWidth(0);
  histESDMBVtxNoPt->SetLineWidth(0);

  histESD->SetMarkerColor(1);
  histESDnsd->SetMarkerColor(6);
  histESDonePart->SetMarkerColor(3);
  histESDMB->SetMarkerColor(2);
  histESDMBVtx->SetMarkerColor(4);

  histESD->SetLineColor(1);
  histESDnsd->SetLineColor(6);
  histESDonePart->SetLineColor(3);
  histESDMB->SetLineColor(2);
  histESDMBVtx->SetLineColor(4);

  histESDNoPt->SetMarkerColor(1);
  histESDMBNoPt->SetMarkerColor(2);
  histESDMBVtxNoPt->SetMarkerColor(4);
  histESDMBTracksNoPt->SetMarkerColor(3);

  histESD->SetMarkerStyle(20);
  histESDnsd->SetMarkerStyle(29);
  histESDonePart->SetMarkerStyle(24);
  histESDMB->SetMarkerStyle(21);
  histESDMBVtx->SetMarkerStyle(22);

  histESDNoPt->SetMarkerStyle(20);
  histESDMBNoPt->SetMarkerStyle(21);
  histESDMBVtxNoPt->SetMarkerStyle(22);
  histESDMBTracksNoPt->SetMarkerStyle(23);
  
  Float_t etaLimit = (fdNdEtaAnalysis->GetAnalysisMode() & AliPWG0Helper::kTPC) ? 0.89 : 1.99;
  Float_t etaPlotLimit = (fdNdEtaAnalysis->GetAnalysisMode() & AliPWG0Helper::kTPC) ? 1.2 : 2.3;
  //Float_t etaLimit = (fdNdEtaAnalysis->GetAnalysisMode() == AliPWG0Helper::kTPC) ? 0.89 : 1.39;
  //Float_t etaPlotLimit = (fdNdEtaAnalysis->GetAnalysisMode() == AliPWG0Helper::kTPC) ? 1.2 : 1.9;

  histESDMBVtx->GetXaxis()->SetRangeUser(-etaLimit, etaLimit);
  histESDMB->GetXaxis()->SetRangeUser(-etaLimit, etaLimit);
  histESD->GetXaxis()->SetRangeUser(-etaLimit, etaLimit);
  histESDnsd->GetXaxis()->SetRangeUser(-etaLimit, etaLimit);
  histESDonePart->GetXaxis()->SetRangeUser(-etaLimit, etaLimit);

  histESDNoPt->GetXaxis()->SetRangeUser(-etaLimit, etaLimit);
  histESDMBNoPt->GetXaxis()->SetRangeUser(-etaLimit, etaLimit);
  histESDMBVtxNoPt->GetXaxis()->SetRangeUser(-etaLimit, etaLimit);
  histESDMBTracksNoPt->GetXaxis()->SetRangeUser(-etaLimit, etaLimit);

  Float_t max = TMath::Max(histESDMBVtx->GetMaximum(), histESDMB->GetMaximum());
  max = TMath::Max(max, histESD->GetMaximum());
  max = TMath::Max(max, histESDMBTracksNoPt->GetMaximum());

  TLegend* legend = new TLegend(0.35, 0.05, 0.75, 0.4);
  legend->SetFillColor(0);
  legend->AddEntry(histESDMBVtx, "Triggered, vertex");
  legend->AddEntry(histESDMB, "Triggered");
  legend->AddEntry(histESD, "All INEL events");
  legend->AddEntry(histESDnsd, "All NSD events");
  legend->AddEntry(histESDonePart, "One Particle");

  TH2F* dummy = new TH2F("dummy", "", 100, -etaPlotLimit, etaPlotLimit, 1000, 0, max * 1.1);
  dummy->GetYaxis()->SetRangeUser(2.1, max * 1.1);
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
  histESDnsd->Draw("SAME");
  histESDonePart->Draw("SAME");
  legend->Draw();
  
  if (save)
  {
    canvas->SaveAs("dNdEta1.png");
    //canvas->SaveAs("dNdEta1.eps");
  }
  
  histESD->Fit("pol0", "0", "", -0.49, 0.49);
  histESDnsd->Fit("pol0", "0", "", -0.49, 0.49);
  histESDonePart->Fit("pol0", "0", "", -0.49, 0.49);
  histESDonePart->Fit("pol0", "0", "", -0.99, 0.99);

  canvas = new TCanvas("dNdEta1_mirrored", "dNdEta1_mirrored", 500, 500);
  canvas->SetGridx();
  canvas->SetGridy();

  dummy->DrawCopy()->GetXaxis()->SetRangeUser(0, 100);
  histESD->DrawCopy("SAME")->SetMarkerStyle(24);
  histESDnsd->DrawCopy("SAME")->SetMarkerStyle(24);
  
  graph = new TGraphErrors(histESD);
  for (Int_t i=0; i<graph->GetN(); i++)
    graph->GetX()[i] *= -1;
  graph->SetMarkerStyle(5);
  graph->Draw("P SAME");

  graph = new TGraphErrors(histESDnsd);
  for (Int_t i=0; i<graph->GetN(); i++)
    graph->GetX()[i] *= -1;
  graph->SetMarkerStyle(5);
  graph->SetMarkerColor(histESDnsd->GetMarkerColor());
  graph->Draw("P SAME");
  
  canvas = new TCanvas("dNdEta1_ratio", "dNdEta1_ratio", 500, 500);
  canvas->SetGridx();
  canvas->SetGridy();
  
  dummy_clone = dummy->DrawCopy();
  dummy_clone->GetXaxis()->SetRangeUser(0, 100);
  dummy_clone->GetYaxis()->SetRangeUser(0.5, 1.5);
  
  graph = new TGraphErrors(histESD);
  for (Int_t i=0; i<graph->GetN(); i++)
  {
    Int_t bin = histESD->GetXaxis()->FindBin(-graph->GetX()[i]);
    if (histESD->GetBinContent(bin) > 0 && graph->GetY()[i] > 0)
    {
      graph->GetEY()[i] = sqrt(graph->GetEY()[i] * graph->GetEY()[i] / graph->GetY()[i] / graph->GetY()[i] 
        + histESD->GetBinError(bin) * histESD->GetBinError(bin) / histESD->GetBinContent(bin) / histESD->GetBinContent(bin));
      graph->GetY()[i] /= histESD->GetBinContent(bin);
      graph->GetEY()[i] *= graph->GetY()[i];
    }
    else
      graph->GetY()[i] = 0;
  }
  graph->SetMarkerStyle(5);
  graph->Draw("P SAME");
  
  graph = new TGraphErrors(histESDnsd);
  for (Int_t i=0; i<graph->GetN(); i++)
  {
    Int_t bin = histESDnsd->GetXaxis()->FindBin(-graph->GetX()[i]);
    if (histESDnsd->GetBinContent(bin) > 0 && graph->GetY()[i] > 0)
    {
      graph->GetEY()[i] = sqrt(graph->GetEY()[i] * graph->GetEY()[i] / graph->GetY()[i] / graph->GetY()[i] 
        + histESDnsd->GetBinError(bin) * histESDnsd->GetBinError(bin) / histESDnsd->GetBinContent(bin) / histESDnsd->GetBinContent(bin));
      graph->GetY()[i] /= histESDnsd->GetBinContent(bin);
      graph->GetEY()[i] *= graph->GetY()[i];
      graph->GetY()[i] += 0.2;
    }
  }
  graph->SetMarkerStyle(5);
  graph->SetMarkerColor(histESDnsd->GetMarkerColor());
  graph->Draw("P SAME");
  
  canvas = new TCanvas("dNdEta1_vertex", "dNdEta1_vertex", 500, 500);
  dummy->DrawCopy();
  histESD->DrawCopy("SAME");
  histESD1->SetLineColor(2);
  histESD1->DrawCopy("SAME");
  histESD2->SetLineColor(4);
  histESD2->DrawCopy("SAME");
  
  PrintIntegratedDeviation(histESDnsd, histESDMB, "factor MB / NSD");
  
  if (onlyESD)
    return;

  loadlibs();

  TFile* file2 = TFile::Open("analysis_mc.root");

  TH1* histMCTrVtx =       (TH1*) GetMCHist("dndetaTrVtx", -1, "MC: MB with vertex")->Clone("histMCTrVtx");
  TH1* ratioTrVtx = (TH1*) DrawdNdEtaRatio(histESDMBVtx, histMCTrVtx, "triggered_vertex", etaPlotLimit)->Clone();
  
  TH1* histMC =            (TH1*) GetMCHist("dndeta", -1, "MC: full inelastic")->Clone("histMC");
  TH1* histMCTr =          (TH1*) GetMCHist("dndetaTr", -1, "MC: minimum bias")->Clone("histMCTr");
  TH1* histMCnsd =         (TH1*) GetMCHist("dndetaNSD", -1, "MC: NSD")->Clone("histMCnsd");
  TH1* histMConePart =     (TH1*) GetMCHist("dndetaOnePart", -1, "MC: OnePart")->Clone("histMConePart");

  TH1* histMCPtCut =       (TH1*) GetMCHist("dndeta", 0.151, "MC: full inelastic, pt cut")->Clone("histMCPtCut");
  TH1* histMCTrPtCut =     (TH1*) GetMCHist("dndetaTr", 0.151, "MC: minimum bias, pt cut")->Clone("histMCTrPtCut");
  TH1* histMCTrVtxPtCut =  (TH1*) GetMCHist("dndetaTrVtx", 0.151, "MC: MB with vertex, pt cut")->Clone("histMCTrVtxPtCut");
  TH1* histMCnsdNoPt =     (TH1*) GetMCHist("dndetaNSD", 0.151, "MC: NSD, put cut")->Clone("histMCnsdNoPt");
  TH1* histMCTracksPtCut = (TH1*) GetMCHist("dndetaTracks", 0.151, "MC: Tracks w/o resolution effect, pt cut")->Clone("histMCTracksPtCut");

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
  histMCTrVtx->SetLineColor(4);

  histMCPtCut->SetLineColor(1);
  histMCTrPtCut->SetLineColor(2);
  histMCTrVtxPtCut->SetLineColor(4);
  if (histMCTracksPtCut)
    histMCTracksPtCut->SetLineColor(3);

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

  TH1* ratio = (TH1*) DrawdNdEtaRatio(histESD, histMC, "full_inelastic", etaPlotLimit)->Clone();
  TH1* ratioTr = (TH1*) DrawdNdEtaRatio(histESDMB, histMCTr, "triggered", etaPlotLimit)->Clone();
  TH1* ratioTrVtx = (TH1*) DrawdNdEtaRatio(histESDMBVtx, histMCTrVtx, "triggered_vertex", etaPlotLimit)->Clone();
  TH1* ratioTrVtxNoPt = (TH1*) DrawdNdEtaRatio(histESDMBVtxNoPt, histMCTrVtxPtCut, "triggered_vertex_nopt", etaPlotLimit)->Clone();
  TH1* ratioNSD = (TH1*) DrawdNdEtaRatio(histESDnsd, histMCnsd, "NSD", etaPlotLimit)->Clone();
  TH1* ratioOnePart = (TH1*) DrawdNdEtaRatio(histESDonePart, histMConePart, "OnePart", etaPlotLimit)->Clone();

  // draw ratios of single steps
  c7 = new TCanvas("all_ratios", "all_ratios", 600, 600);
  c7->SetRightMargin(0.05);
  c7->SetTopMargin(0.05);
  c7->SetGridx();
  c7->SetGridy();
  
  ratioTrVtxNoPt->SetMarkerStyle(20);
  ratioTrVtx->SetMarkerStyle(21);
  ratioTr->SetMarkerStyle(23);
  ratio->SetMarkerStyle(22);
  ratioNSD->SetMarkerStyle(26);
  
  ratioTrVtxNoPt->SetMarkerSize(2);
  ratioTrVtx->SetMarkerSize(2);
  ratioTr->SetMarkerSize(2);
  ratio->SetMarkerSize(2);
  ratioNSD->SetMarkerSize(2);
  
  ratioTrVtxNoPt->SetMarkerColor(1);
  ratioTrVtx->SetMarkerColor(2);
  ratioTr->SetMarkerColor(4);
  ratio->SetMarkerColor(2);
  ratioNSD->SetMarkerColor(1);
  
  ratioTrVtxNoPt->SetLineColor(1);
  ratioTrVtx->SetLineColor(2);
  ratioTr->SetLineColor(4);
  ratio->SetLineColor(2);
  ratioNSD->SetLineColor(1);
  
  legend7 = new TLegend(0.13, 0.7, 0.94, 0.9);
  legend7->SetFillColor(0);
  legend7->SetTextSize(0.035);
  legend7->SetNColumns(2);
  
  flat = new TF1("flat", "-1", -5, 5);
  ratioTrVtxNoPt->Add(flat);
  ratioTrVtx->Add(flat);
  ratioTr->Add(flat);
  ratio->Add(flat);
  ratioNSD->Add(flat);
  
  ratioTrVtxNoPt->Scale(100);
  ratioTrVtx->Scale(100);
  ratioTr->Scale(100);
  ratio->Scale(100);
  ratioNSD->Scale(100);
  
  ratio->Add(ratioTr, -1);
  ratioNSD->Add(ratioTr, -1);
  ratioTr->Add(ratioTrVtx, -1);
  ratioTrVtx->Add(ratioTrVtxNoPt, -1);
  
  legend7->AddEntry(ratioTrVtxNoPt, "Track-to-particle", "P");
  legend7->AddEntry(ratio, "Trigger-bias INEL", "P");
  legend7->AddEntry(ratioTr, "Vertex-reconstruction", "P");
  legend7->AddEntry(ratioNSD, "Trigger-bias NSD", "P");
  if ((fdNdEtaAnalysis->GetAnalysisMode() & AliPWG0Helper::kFieldOn) && (fdNdEtaAnalysis->GetAnalysisMode() & AliPWG0Helper::kTPC || fdNdEtaAnalysis->GetAnalysisMode() & AliPWG0Helper::kTPCITS))
    legend7->AddEntry(ratioTrVtx, "p_{T} cut-off", "P");
  
  TH1* dummy7 = new TH2F("dummy7", ";#eta;Deviation in %", 100, -etaPlotLimit, etaPlotLimit, 100, -5, 7);
  dummy7->SetStats(0);
  dummy7->Draw();
  
  ratio->GetXaxis()->SetRangeUser(-etaLimit, etaLimit);
  ratioTr->GetXaxis()->SetRangeUser(-etaLimit, etaLimit);
  ratioTrVtx->GetXaxis()->SetRangeUser(-etaLimit, etaLimit);
  ratioTrVtxNoPt->GetXaxis()->SetRangeUser(-etaLimit, etaLimit);
  ratioNSD->GetXaxis()->SetRangeUser(-etaLimit, etaLimit);
  
  ratio->Draw("HIST EP SAME");
  ratioTr->Draw("HIST EP SAME");
  if ((fdNdEtaAnalysis->GetAnalysisMode() & AliPWG0Helper::kFieldOn) && (fdNdEtaAnalysis->GetAnalysisMode() & AliPWG0Helper::kTPC || fdNdEtaAnalysis->GetAnalysisMode() & AliPWG0Helper::kTPCITS))
    ratioTrVtx->Draw("HIST EP SAME");
  ratioTrVtxNoPt->Draw("HIST EP SAME");
  ratioNSD->Draw("HIST EP SAME");
  legend7->Draw();
  
  //c7->SaveAs("ratios.eps");

  new TCanvas;
  dummy2->DrawCopy();
  histMCnsd->Draw("SAME");
  histESDnsd->Draw("SAME");

  ratio = (TH1*) histMC->Clone("ratio");
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
  PrintIntegratedDeviation(histMConePart, histESDonePart, "all events (INEL>0)");
  PrintIntegratedDeviation(histMCTr, histESDMB, "triggered");
  PrintIntegratedDeviation(histMCTrVtx, histESDMBVtx, "trigger, vertex");
  PrintIntegratedDeviation(histMCPtCut, histESDNoPt, "all events (no pt corr)");
  PrintIntegratedDeviation(histMCnsdNoPt, histESDnsdNoPt, "all events (NSD) (no pt corr)");
  PrintIntegratedDeviation(histMCTrPtCut, histESDMBNoPt, "triggered (no pt corr)");
  PrintIntegratedDeviation(histMCTrVtxPtCut, histESDMBVtxNoPt, "trigger, vertex (no pt corr)");
  PrintIntegratedDeviation(histESD, histESDNoPt, "pt cut off correction");

  TCanvas* canvas3 = new TCanvas("dNdEta", "dNdEta", 600, 600);
  canvas3->Range(0, 0, 1, 1);
  //canvas3->Divide(1, 2, 0, 0);

  //canvas3->cd(1);
  TPad* pad1 = new TPad("dNdEta_1", "", 0, 0.5, 0.98, 0.98);
  pad1->SetTopMargin(0.05);
  pad1->SetLeftMargin(0.13);
  pad1->Draw();

  TPad* pad2 = new TPad("dNdEta_2", "", 0, 0.02, 0.98, 0.5);
  pad2->SetLeftMargin(0.13);
  pad2->Draw();

  pad1->SetRightMargin(0.01);
  pad2->SetRightMargin(0.01);

  // no border between them
  pad1->SetBottomMargin(0);
  pad2->SetTopMargin(0);

  pad1->cd();
  pad1->SetGridx();
  pad1->SetGridy();

  legend->AddEntry(histMC, "MC prediction");

  dummy->GetXaxis()->SetLabelSize(0.08);
  dummy->GetYaxis()->SetLabelSize(0.08);
  dummy->GetXaxis()->SetTitleSize(0.08);
  dummy->GetYaxis()->SetTitleSize(0.08);
  dummy->GetYaxis()->SetTitleOffset(0.8);
  dummy->DrawCopy();
  histESDMBVtx->Draw("SAME");
  histESDMB->Draw("SAME");
  histESD->Draw("SAME");
  histMC->Draw("SAME");

  legend->SetTextSize(0.08);
  legend->Draw();

  pad2->cd();
  pad2->SetBottomMargin(0.15);
  //pad2->SetGridx();
  //pad2->SetGridy();

  Float_t minR = 0.91; //TMath::Min(0.961, ratio->GetMinimum() * 0.95);
  Float_t maxR = 1.09; //TMath::Max(1.049, ratio->GetMaximum() * 1.05);

  TH1F dummy3("dummy3", ";#eta;Ratio: MC / corr", 100, -etaPlotLimit, etaPlotLimit);
  dummy3.SetStats(kFALSE);
  for (Int_t i=1; i<=100; ++i)
    dummy3.SetBinContent(i, 1);
  dummy3.GetYaxis()->SetRangeUser(minR, maxR);
  dummy3.SetLineWidth(2);
  dummy3.GetXaxis()->SetLabelSize(0.08);
  dummy3.GetYaxis()->SetLabelSize(0.08);
  dummy3.GetXaxis()->SetTitleSize(0.08);
  dummy3.GetYaxis()->SetTitleSize(0.08);
  dummy3.GetYaxis()->SetTitleOffset(0.8);
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

void CompareTwodNdEta(const char* fileName1, const char* fileName2, Bool_t errorsCorrelated = kFALSE)
{
  c = new TCanvas;
  
  c->SetGridx();
  c->SetGridy();

  hist = new TH2F("dummy", ";#eta;dN_{ch}/d#eta", 100, -2.5, 2.5, 100, 0, 8);
  hist->SetStats(0);
  hist->DrawCopy();//->GetYaxis()->SetRangeUser(2, 4.5);
  
  l = new TLegend(0.2, 0.13, 0.8, 0.35);
  l->SetNColumns(2);
  l->SetFillColor(0);
  
  TH1* histESD[2];
  TH1* histESDnsd[2];
  
  for (Int_t i=0; i<2; i++)
  {
    if (i == 0)
      file = TFile::Open(fileName1);
    if (i == 1)
    {
      if (fileName2 == 0)
        break;
      file = TFile::Open(fileName2);
    }
  
    histESD[i] = (TH1*) file->Get("dndeta/dNdEta_corrected");
    histESDnsd[i] = (TH1*) file->Get("dndetaNSD/dNdEta_corrected");
    
    histESD[i]->SetMarkerStyle(20 + i*4);
    histESDnsd[i]->SetMarkerStyle(21 + i*4);
    
    histESD[i]->SetMarkerColor(i+1);
    histESD[i]->SetLineColor(i+1);
    histESDnsd[i]->SetMarkerColor(i+1);
    histESDnsd[i]->SetLineColor(i+1);
    
    histESD[i]->DrawCopy("SAME");
    histESDnsd[i]->DrawCopy("SAME");
    
    l->AddEntry(histESD[i], Form("Data %d INEL", i), "P");
    l->AddEntry(histESDnsd[i], Form("Data %d NSD", i), "P");
  }

  if (0)
  {
    TGraphErrors *gre = new TGraphErrors(16);
    gre->SetFillColor(4);
    gre->SetMarkerColor(4);
    gre->SetMarkerStyle(26);
    gre->SetPoint(0,0.125,3.14);
    gre->SetPointError(0,0,0.07);
    gre->SetPoint(1,0.375,3.04);
    gre->SetPointError(1,0,0.07);
    gre->SetPoint(2,0.625,3.17);
    gre->SetPointError(2,0,0.07);
    gre->SetPoint(3,0.875,3.33);
    gre->SetPointError(3,0,0.07);
    gre->SetPoint(4,1.125,3.33);
    gre->SetPointError(4,0,0.07);
    gre->SetPoint(5,1.375,3.53);
    gre->SetPointError(5,0,0.07);
    gre->SetPoint(6,1.625,3.46);
    gre->SetPointError(6,0,0.07);
    gre->SetPoint(7,1.875,3.41);
    gre->SetPointError(7,0,0.07);
    gre->SetPoint(8,-0.125,3.14);
    gre->SetPointError(8,0,0.07);
    gre->SetPoint(9,-0.375,3.04);
    gre->SetPointError(9,0,0.07);
    gre->SetPoint(10,-0.625,3.17);
    gre->SetPointError(10,0,0.07);
    gre->SetPoint(11,-0.875,3.33);
    gre->SetPointError(11,0,0.07);
    gre->SetPoint(12,-1.125,3.33);
    gre->SetPointError(12,0,0.07);
    gre->SetPoint(13,-1.375,3.53);
    gre->SetPointError(13,0,0.07);
    gre->SetPoint(14,-1.625,3.46);
    gre->SetPointError(14,0,0.07);
    gre->SetPoint(15,-1.875,3.41);
    gre->SetPointError(15,0,0.07);
    gre->Draw("p");
    
    l->AddEntry(gre, "UA5 INEL", "P");
    
    gre = new TGraphErrors(16);
    gre->SetMarkerColor(4);
    gre->SetFillColor(4);
    gre->SetMarkerStyle(22);
    gre->SetPoint(0,0.125,3.48);
    gre->SetPointError(0,0,0.07);
    gre->SetPoint(1,0.375,3.38);
    gre->SetPointError(1,0,0.07);
    gre->SetPoint(2,0.625,3.52);
    gre->SetPointError(2,0,0.07);
    gre->SetPoint(3,0.875,3.68);
    gre->SetPointError(3,0,0.07);
    gre->SetPoint(4,1.125,3.71);
    gre->SetPointError(4,0,0.07);
    gre->SetPoint(5,1.375,3.86);
    gre->SetPointError(5,0,0.07);
    gre->SetPoint(6,1.625,3.76);
    gre->SetPointError(6,0,0.07);
    gre->SetPoint(7,1.875,3.66);
    gre->SetPointError(7,0,0.07);
    gre->SetPoint(8,-0.125,3.48);
    gre->SetPointError(8,0,0.07);
    gre->SetPoint(9,-0.375,3.38);
    gre->SetPointError(9,0,0.07);
    gre->SetPoint(10,-0.625,3.52);
    gre->SetPointError(10,0,0.07);
    gre->SetPoint(11,-0.875,3.68);
    gre->SetPointError(11,0,0.07);
    gre->SetPoint(12,-1.125,3.71);
    gre->SetPointError(12,0,0.07);
    gre->SetPoint(13,-1.375,3.86);
    gre->SetPointError(13,0,0.07);
    gre->SetPoint(14,-1.625,3.76);
    gre->SetPointError(14,0,0.07);
    gre->SetPoint(15,-1.875,3.66);
    gre->SetPointError(15,0,0.07);
    gre->Draw("p");
    
    l->AddEntry(gre, "UA5 NSD", "P");
  }

  l->Draw();
  
  if (fileName2 == 0)
    return;
  
  new TCanvas;
  gPad->SetGridx();
  gPad->SetGridy();
  
  if (errorsCorrelated)
  {
    for (Int_t i=1; i<=histESD[1]->GetNbinsX(); i++)
    {
      histESD[1]->SetBinError(i, 0);
      histESDnsd[1]->SetBinError(i, 0);
    }
  }
  
  histESD[0]->Divide(histESD[0], histESD[1]);
  histESDnsd[0]->Divide(histESDnsd[0], histESDnsd[1]);
  
  for (Int_t i=1; i<=histESD[1]->GetNbinsX(); i++)
    histESDnsd[0]->SetBinContent(i, histESDnsd[0]->GetBinContent(i) + 0.2);
  
  hist->DrawCopy()->GetYaxis()->SetRangeUser(0.8, 1.4);
  histESD[0]->Draw("SAME");
  histESDnsd[0]->Draw("SAME");
}

TH1* DrawdNdEtaRatio(TH1* corr, TH1* mc, const char* name, Float_t etaPlotLimit)
{
  TCanvas* canvas3 = new TCanvas(name, name, 600, 600);
  canvas3->Range(0, 0, 1, 1);

  TPad* pad1 = new TPad(Form("%s_1", name), "", 0, 0.5, 0.98, 0.98);
  pad1->Draw();

  TPad* pad2 = new TPad(Form("%s_2", name), "", 0, 0.02, 0.98, 0.5);
  pad2->Draw();

  pad1->SetRightMargin(0.01);
  pad2->SetRightMargin(0.01);
  pad1->SetTopMargin(0.05);
  pad1->SetLeftMargin(0.13);
  pad2->SetLeftMargin(0.13);
  pad2->SetBottomMargin(0.15);
  
  // no border between them
  pad1->SetBottomMargin(0);
  pad2->SetTopMargin(0);

  pad1->cd();
  pad1->SetGridx();
  pad1->SetGridy();

  TLegend* legend = new TLegend(0.35, 0.05, 0.75, 0.3);
  legend->SetFillColor(0);
  legend->AddEntry(corr, "Corrected");
  legend->AddEntry(mc, "MC prediction");
  legend->SetTextSize(0.08);

  TH2F* dummy = new TH2F("dummy", "", 100, -etaPlotLimit, etaPlotLimit, 1000, 2.7, corr->GetMaximum() * 1.1);
  Prepare1DPlot(dummy);
  dummy->SetStats(kFALSE);
  dummy->SetXTitle("#eta");
  dummy->SetYTitle("dN_{ch}/d#eta");
  dummy->GetYaxis()->SetTitleOffset(1);

  dummy->GetXaxis()->SetLabelSize(0.08);
  dummy->GetYaxis()->SetLabelSize(0.08);
  dummy->GetXaxis()->SetTitleSize(0.08);
  dummy->GetYaxis()->SetTitleSize(0.08);
  dummy->GetYaxis()->SetTitleOffset(0.8);
  dummy->DrawCopy();

  corr->Draw("SAME");
  mc->Draw("SAME");

  legend->Draw();

  pad2->cd();
  pad2->SetBottomMargin(0.15);
  //pad2->SetGridx();
  //pad2->SetGridy();

  TH1* ratio = (TH1*) mc->Clone("ratio");
  ratio->Divide(corr);

  Float_t minR = TMath::Min(0.91, ratio->GetMinimum() * 0.95);
  Float_t maxR = TMath::Max(1.09, ratio->GetMaximum() * 1.05);

  TH1F dummy3("dummy3", ";#eta;Ratio: MC / corr", 100, -etaPlotLimit, etaPlotLimit);
  dummy3.SetStats(kFALSE);
  for (Int_t i=1; i<=100; ++i)
  	dummy3.SetBinContent(i, 1);
  dummy3.GetYaxis()->SetRangeUser(minR, maxR);
  dummy3.SetLineWidth(2);
  dummy3.GetXaxis()->SetLabelSize(0.08);
  dummy3.GetYaxis()->SetLabelSize(0.08);
  dummy3.GetXaxis()->SetTitleSize(0.08);
  dummy3.GetYaxis()->SetTitleSize(0.08);
  dummy3.GetYaxis()->SetTitleOffset(0.8);
  dummy3.DrawCopy();

  ratio->Draw("SAME");

  canvas3->Modified();

  return ratio;
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

  TH1* hist = dNdEtaCorrection->GetTriggerBiasCorrectionOnePart()->GetEventCorrection()->Get1DCorrection("x");
  TH1* hist2 = dNdEtaCorrection->GetTriggerBiasCorrectionOnePart()->GetEventCorrection()->Get1DCorrection("y", -5, 5);

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
  pave->AddText("|z| < 5 cm");
  pave->Draw();
  
  Float_t triggerEff = 100.0 / hist2->GetBinContent(1);
  Printf("trigger eff in 0 bin is: %.2f %%", triggerEff);
  
  return;

  TH1* hist2 = dNdEtaCorrection->GetVertexRecoCorrection()->GetEventCorrection()->Get1DCorrection("y", -10, 10);
  //new TCanvas;
  //hist2->Draw();

  Printf("vertex reco eff in 0 bin is: %.2f %%", 100.0 / hist2->GetBinContent(1));
  
  Printf("combined efficiency is %.2f %%", triggerEff / hist2->GetBinContent(1));
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

  corrX->GetYaxis()->SetTitle("Correction factor");
  corrZ->GetYaxis()->SetTitle("Correction factor");

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

void Correction1DCreatePlots(const char* fileName = "correction_map.root", const char* folderName = "dndeta_correction", Float_t upperPtLimit = 9.9, Int_t correctionType = 0, Int_t correctionType2 = -1)
{
  if (correctionType2 == -1)
    correctionType2 = correctionType;

  TFile::Open(fileName);
  AlidNdEtaCorrection* dNdEtaCorrection = new AlidNdEtaCorrection(folderName, folderName);
  dNdEtaCorrection->LoadHistograms();

  TH3F* gene = dNdEtaCorrection->GetCorrection(correctionType)->GetTrackCorrection()->GetGeneratedHistogram();
  TH3F* meas = dNdEtaCorrection->GetCorrection(correctionType2)->GetTrackCorrection()->GetMeasuredHistogram();

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

TCanvas* Correction1D(Int_t correctionType = 0, const char* fileName = "correction_map.root", const char* folder = "dndeta_correction", Float_t upperPtLimit = 9.9, Int_t correctionType2 = -1)
{
  gSystem->Load("libPWG0base");

  Correction1DCreatePlots(fileName, folder, upperPtLimit, correctionType, correctionType2);

  TH1* corrX = dynamic_cast<TH1*> (gROOT->FindObject(Form("generated_x_div_measured_x", folder, folder)));
  TH1* corrY = dynamic_cast<TH1*> (gROOT->FindObject(Form("generated_y_div_measured_y", folder, folder)));
  TH1* corrZ = dynamic_cast<TH1*> (gROOT->FindObject(Form("generated_z_div_measured_z", folder, folder)));

  Prepare1DPlot(corrX);
  Prepare1DPlot(corrY);
  Prepare1DPlot(corrZ);

  /*
  corrX->SetTitle("a) z projection");
  corrY->SetTitle("b) #eta projection");
  corrZ->SetTitle("c) p_{T} projection");
  */
  
  corrX->SetTitle("");
  corrY->SetTitle("");
  corrZ->SetTitle("");

  corrX->SetTitleSize(0.06, "xyz");
  corrX->SetLabelSize(0.06, "xyz");
  corrY->SetTitleSize(0.06, "xyz");
  corrY->SetLabelSize(0.06, "xyz");
  corrZ->SetTitleSize(0.06, "xyz");
  corrZ->SetLabelSize(0.06, "xyz");

  corrX->GetYaxis()->SetTitle("Correction factor");
  corrY->GetYaxis()->SetTitle("Correction factor");
  corrZ->GetYaxis()->SetTitle("Correction factor");
  //corrX->GetYaxis()->SetTitleOffset(1.7);
  //corrY->GetYaxis()->SetTitleOffset(1.7);
  //corrZ->GetYaxis()->SetTitleOffset(1.7);
  corrX->GetYaxis()->SetRangeUser(0.8, 1.5);
  corrY->GetYaxis()->SetRangeUser(0.8, 1.5);
  corrZ->GetYaxis()->SetRangeUser(0.8, 1.5);

  corrZ->GetXaxis()->SetRangeUser(0.11, upperPtLimit);

  TString canvasName;
  canvasName.Form(Form("Correction1D_%d_%s_%f", correctionType, fileName, upperPtLimit));
  TCanvas* canvas = new TCanvas(canvasName, canvasName, 1200, 400);
  canvas->Divide(3, 1);

  TLatex* Tl = new TLatex;
  Tl->SetTextSize(0.06);
  Tl->SetBit(TLatex::kTextNDC);

  canvas->cd(1);
  InitPad();
  gPad->SetTopMargin(0.05);
  gPad->SetBottomMargin(0.15);
  corrX->DrawCopy();
  Tl->DrawLatex(0.5, 0.88, "0.3 < p_{T} < 10");
  Tl->DrawLatex(0.5, 0.8, "|#eta| < 0.8");

  canvas->cd(2);
  InitPad();
  gPad->SetTopMargin(0.05);
  gPad->SetBottomMargin(0.15);
  corrY->Draw();
  Tl->DrawLatex(0.5, 0.88, "0.3 < p_{T} < 10");
  Tl->DrawLatex(0.5, 0.8, "|vtx-z| < 10 cm");

  canvas->cd(3);
  InitPad();
  gPad->SetTopMargin(0.05);
  gPad->SetBottomMargin(0.15);
  gPad->SetLogx();
  corrZ->Draw();
  corrZ->GetXaxis()->SetLabelOffset(0.005);
  corrZ->GetXaxis()->SetTitleOffset(1.2);
  Tl->DrawLatex(0.5, 0.88, "|vtx-z| < 10 cm");
  Tl->DrawLatex(0.5, 0.8, "|#eta| < 0.8");

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

/*
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
*/

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

void Track2Particle2DCreatePlots(const char* fileName = "correction_map.root", const char* folder = "dndeta_correction")
{
  TFile::Open(fileName);
  AlidNdEtaCorrection* dNdEtaCorrection = new AlidNdEtaCorrection(folder, folder);
  dNdEtaCorrection->LoadHistograms();

  TH3F* gene = dNdEtaCorrection->GetTrack2ParticleCorrection()->GetTrackCorrection()->GetGeneratedHistogram();
  TH3F* meas = dNdEtaCorrection->GetTrack2ParticleCorrection()->GetTrackCorrection()->GetMeasuredHistogram();

  gene->GetZaxis()->SetRangeUser(0.2, 10);
  meas->GetZaxis()->SetRangeUser(0.2, 10);
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

  Track2Particle2DCreatePlots(fileName, folder);

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

void CompareCorrection2Measured(Float_t ptMin = 0.301, const char* dataInput = "analysis_esd_raw.root", const char* correctionMapFile = "correction_map.root", const char* correctionMapFolder = "dndeta_correction")
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
  Printf("mc contains %f entries in |vtx-z| < 10, |eta| < 1, pt > 0.3", hist1->Integral(hist1->GetXaxis()->FindBin(-9.9), hist1->GetXaxis()->FindBin(9.9), hist1->GetYaxis()->FindBin(-0.99), hist1->GetYaxis()->FindBin(0.99), hist1->GetZaxis()->FindBin(ptMin), hist1->GetNbinsZ()));

  TH3* hist2 = (TH3*) fdNdEtaAnalysis->GetData()->GetTrackCorrection()->GetMeasuredHistogram()->Clone("esd");
  hist2->SetTitle("esd");
  Printf("esd contains %f entries", hist2->Integral());
  Printf("esd contains %f entries in |vtx-z| < 10, |eta| < 1, pt > 0.3", hist2->Integral(hist2->GetXaxis()->FindBin(-9.9), hist2->GetXaxis()->FindBin(9.9), hist2->GetYaxis()->FindBin(-0.99), hist2->GetYaxis()->FindBin(0.99), hist2->GetZaxis()->FindBin(ptMin), hist2->GetNbinsZ()));

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

  TH2* hist3 = (TH2*) dNdEtaCorrection->GetVertexRecoCorrection()->GetEventCorrection()->GetMeasuredHistogram()->Clone("mc2");
  hist3->SetTitle("mc2");
  Printf("mc event contains %f entries", hist3->Integral());
  Printf("mc event contains %f entries in |vtx-z| < 10", hist3->Integral(hist3->GetXaxis()->FindBin(-9.9), hist3->GetXaxis()->FindBin(9.9), 1, hist3->GetNbinsY()));

  TH2* hist4 = (TH2*) fdNdEtaAnalysis->GetData()->GetEventCorrection()->GetMeasuredHistogram()->Clone("esd2");
  hist4->SetTitle("esd2");
  Printf("esd event contains %f entries", hist4->Integral());
  Printf("esd event contains %f entries in |vtx-z| < 10", hist4->Integral(hist4->GetXaxis()->FindBin(-9.9), hist4->GetXaxis()->FindBin(9.9), 1, hist4->GetNbinsY()));
  
  ratio = (TH2*) hist3->Clone("ratio");
  ratio->Divide(hist4);
  
  new TCanvas; ratio->Draw("COLZ");
}

void CompareCorrection2Generated(Float_t ptMin = 0.301, const char* dataInput = "analysis_mc.root", const char* correctionMapFile = "correction_map.root", const char* correctionMapFolder = "dndeta_correction")
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

  dNdEtaAnalysis* fdNdEtaAnalysis = new dNdEtaAnalysis("dndetaTrVtx", "dndetaTrVtx");
  fdNdEtaAnalysis->LoadHistograms("dndetaTrVtx");

  gROOT->cd();
  
  TH3* hist1 = (TH3*) dNdEtaCorrection->GetTrack2ParticleCorrection()->GetTrackCorrection()->GetGeneratedHistogram()->Clone("mc");
  hist1->SetTitle("mc");
  Printf("mc contains %f entries", hist1->Integral());
  Printf("mc contains %f entries in |vtx-z| < 10, pt > 0.3", hist1->Integral(hist1->GetXaxis()->FindBin(-9.9), hist1->GetXaxis()->FindBin(9.9), hist1->GetYaxis()->FindBin(-0.99), hist1->GetYaxis()->FindBin(0.99), hist1->GetZaxis()->FindBin(ptMin), hist1->GetNbinsZ()));

  TH3* hist2 = (TH3*) fdNdEtaAnalysis->GetData()->GetTrackCorrection()->GetGeneratedHistogram()->Clone("esd");
  hist2->SetTitle("esd");
  Printf("esd contains %f entries", hist2->Integral());
  Printf("esd contains %f entries in |vtx-z| < 10, pt > 0.3", hist2->Integral(hist2->GetXaxis()->FindBin(-9.9), hist2->GetXaxis()->FindBin(9.9), hist2->GetYaxis()->FindBin(-0.99), hist2->GetYaxis()->FindBin(0.99), hist2->GetZaxis()->FindBin(ptMin), hist2->GetNbinsZ()));

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

void DrawTrackletOrigin(const char* fileName = "correction_map.root", Bool_t myFile = kTRUE)
{
  TFile::Open(fileName);

  Int_t maxHists = 8;
  TH1* hist[8];
  
  const Int_t kRebin = 8;

  const char* titles[] = { "PP", "SS", "PP'", "PS'", "PS", "SP'", "SS'", "" };

  if (myFile)
  {
    for (Int_t i=0; i<maxHists; i++)
    {
      hist[i] = (TH1*) gFile->Get(Form("fDeltaPhi_%d", i));
      if (hist[i]->GetDimension() == 2)
        hist[i] = ((TH2*) hist[i])->ProjectionX(Form("fDeltaPhi_clone_%d", i));
    }
  }
  else
  {
    maxHists = 6;
    const char* names[] = { "DePhiPPTracklets", "DePhiSecTracklets", "DePhiPpTracklets", "DePhiPSTracklets", "DePhiPSdaugTracklets", "DePhiSPTracklets" }; 
    for (Int_t i=0; i<maxHists; i++)
      hist[i] = (TH1*) gFile->Get(names[i]);
  }
  
  // clone before rebinning
  good = (TH1*) hist[0]->Clone("good");
  good->Add(hist[4]);
  
  bad = (TH1*) hist[1]->Clone("bad");
  bad->Add(hist[2]);
  bad->Add(hist[3]);
  bad->Add(hist[5]);
  if (myFile)
    bad->Add(hist[6]);
  
  c = (TCanvas*) gROOT->GetListOfCanvases()->FindObject("c");
  TH1* ref = 0;
  Bool_t nw = kFALSE;
  if (!c)
  {
    c = new TCanvas("c", "c", 600, 600);
    nw = kTRUE;
    ref = (TH1*) c->GetListOfPrimitives()->At(1);
  }  
  c->cd();
  c->SetRightMargin(0.05);
  c->SetTopMargin(0.05);
  c->SetLogy();
  c->SetGridx();
  c->SetGridy();
  
  Int_t order[] = { 0, 4, 1, 2, 3, 5, 6, 7 };
  //Int_t colors[]  = {1,2,4,1,2,4,1,2,4};
  Int_t colors[]  = {1,2,3,4,6,7,8,102};
  Int_t markers[]  = {20, 21, 22, 23, 24, 25, 26, 27, 28};
  
  TLegend* legend = new TLegend(0.75, 0.6, 0.93, 0.93);
  legend->SetFillColor(0);
  legend->SetTextSize(0.04);

  Int_t total = 0;
  for (Int_t ii=0; ii<maxHists; ii++)
  {
    i = order[ii];
    
    hist[i]->Rebin(kRebin);
    hist[i]->SetStats(kFALSE);
    hist[i]->SetLineColor(colors[i]);
    hist[i]->SetLineWidth(2);
    //hist[i]->SetMarkerStyle(markers[i]);
    //hist[i]->SetMarkerColor(colors[i]);
    //hist[i]->SetLineStyle(ii+1);
    hist[i]->GetXaxis()->SetRangeUser(-0.09, 0.09);
    hist[i]->GetYaxis()->SetRangeUser(5, hist[i]->GetMaximum() * 2);
    hist[i]->GetYaxis()->SetTitleOffset(1.3);
    hist[i]->GetXaxis()->SetTitle("#Delta#varphi (rad.)");
    
    if (i == 0 && ref)
      hist[i]->Scale(1.0 / hist[i]->GetMaximum() * ref->GetMaximum());
    
    hist[i]->DrawCopy(((i == 0 && nw) ? "" : "SAME"));

    total += hist[i]->GetEntries();

    if (i != 7)
      legend->AddEntry(hist[i], titles[i], "L");
  }

  legend->Draw();
  c->SaveAs("spd_tracklets_deltaphi_detailed.eps");

  Printf("Total: %d", total);
  for (Int_t i=0; i<maxHists; i++)
    Printf("Histogram %d (%s) contains %.2f %% of the entries", i, titles[i], 100.0 * hist[i]->GetEntries() / total);

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
  
  eff = new TH1F("eff", ";#Delta#varphi cut (rad.)", 101,-0.0005, 0.1005);
  cont = new TH1F("cont", "cont", 101,-0.0005, 0.1005);
  signalOverBg = new TH1F("signalOverBg", "signalOverBg", 101,-0.0005, 0.1005);
  for (Float_t cut=0.000; cut<0.10; cut += 0.001)
  {
    Float_t accGood = good->Integral(good->GetXaxis()->FindBin(-cut), good->GetXaxis()->FindBin(cut));
    Float_t accBad = bad->Integral(bad->GetXaxis()->FindBin(-cut), bad->GetXaxis()->FindBin(cut));
    Float_t sB = accGood / accBad;
    eff->Fill(cut, 100.0 * accGood / good->Integral());
    cont->Fill(cut, 100.0 * accBad / (accGood + accBad));
    signalOverBg->Fill(cut, sB);
  }
  
  //new TCanvas; signalOverBg->Draw();
  
  c = (TCanvas*) gROOT->GetListOfCanvases()->FindObject("c2");
  Bool_t nw = kFALSE;
  if (!c)
  {
    c = new TCanvas("c2", "c2", 600, 600);
    nw = kTRUE;
  }
  c->cd();
  c->SetRightMargin(0.05);
  c->SetTopMargin(0.05);
  c->SetGridx();
  c->SetGridy();
  gPad->SetLogy();
  good->Rebin(kRebin);
  bad->Rebin(kRebin);
  good->GetXaxis()->SetRangeUser(-0.09, 0.09);
  good->GetYaxis()->SetTitleOffset(1.3);
  good->SetStats(0);
  good->GetXaxis()->SetTitle("#Delta#varphi (rad.)");  
  good->DrawCopy((nw) ? "" : "SAME");
  
  bad->SetLineColor(2);
  bad->SetLineStyle(2);
  bad->SetLineWidth(2);
  //bad->SetMarkerColor(2);
  //bad->SetMarkerStyle(7);
  bad->DrawCopy("SAME");
  
  TLegend* legend = new TLegend(0.2, 0.13, 0.85, 0.25);
  legend->SetFillColor(0);
  legend->SetTextSize(0.04);
  legend->AddEntry(good, "Primaries", "L");
  legend->AddEntry(bad, "Secondaries + Background", "L");
  legend->Draw();
  
  c->SaveAs("spd_tracklets_deltaphi.eps");
  
  c = (TCanvas*) gROOT->GetListOfCanvases()->FindObject("c3");
  Bool_t nw = kFALSE;
  if (!c)
  {
    c = new TCanvas("c3", "c3", 600, 600);
    nw = kTRUE;
  }
  c->cd();
  c->SetRightMargin(0.05);
  c->SetTopMargin(0.05);
  c->SetGridx();
  c->SetGridy();
  
  TLegend* legend = new TLegend(0.5, 0.6, 0.93, 0.75);
  legend->SetFillColor(0);
  legend->SetTextSize(0.04);
  legend->AddEntry(eff, "Efficiency (%)", "L");
  legend->AddEntry(cont, "Contamination (%)", "L");
  
  eff->SetStats(0);
  eff->GetXaxis()->SetRangeUser(0, 0.08);
  eff->GetYaxis()->SetRangeUser(1e-3, 105);
  eff->SetLineWidth(2);
  eff->DrawCopy((nw) ? "" : "SAME");
  cont->SetLineStyle(2);
  cont->SetLineWidth(2);
  cont->SetLineColor(2);
  cont->DrawCopy("SAME");
  legend->Draw();
  
  c->SaveAs("spd_tracklets_efficiency.eps");
}

void DrawTrackletOrigin_Compare(const char* file1, const char* file2)
{
  DrawTrackletOrigin(file1);
  good1 = (TH1*) gROOT->FindObject("good")->Clone("good1");
  bad1 = (TH1*) gROOT->FindObject("bad")->Clone("bad1");

  DrawTrackletOrigin(file2);
  good2 = (TH1*) gROOT->FindObject("good")->Clone("good2");
  bad2 = (TH1*) gROOT->FindObject("bad")->Clone("bad2");
     
  c = new TCanvas("c4", "c4", 600, 600);
  c->SetRightMargin(0.05);
  c->SetTopMargin(0.05);
  c->SetGridx();
  c->SetGridy();
  gPad->SetLogy();
  
  good1->Draw();
  bad1->SetLineColor(1);
  bad1->SetMarkerColor(1);
  bad1->Draw("SAME");
  
  Float_t factor = (good1->Integral() + bad1->Integral()) / (good2->Integral() + bad2->Integral());
  
  good2->Scale(factor);
  bad2->Scale(factor);
  
  good2->SetLineColor(2);
  bad2->SetMarkerColor(2);
  
  good2->Draw("SAME");
  bad2->Draw("SAME");
  
  good1->GetYaxis()->SetRangeUser(1, TMath::Max(good1->GetMaximum(), good2->GetMaximum()) * 1.1);
}
  
void Tracklets_Asymmetry()
{
  TFile::Open("correction_map.root");

  Int_t maxHists = 7;
  TH1* hist[8];

  Int_t colors[]  = {1,2,3,4,6,7,8,102};
  const char* titles[] = { "PP", "SS", "PP'", "PS'", "PS", "SP'", "SS'", "" };

  TLegend* legend = new TLegend(0.75, 0.6, 0.93, 0.93);
  
  for (Int_t i=0; i<maxHists; i++)
  {
    hist[i] = (TH1*) gFile->Get(Form("fDeltaPhi_%d", i));
    hist[i]->Rebin(10);
    
    for (Int_t j=hist[i]->GetNbinsX()/2; j<=hist[i]->GetNbinsX(); j++)
      if (hist[i]->GetBinContent(j) > 0)
        hist[i]->SetBinContent(j, (hist[i]->GetBinContent(j) -  hist[i]->GetBinContent(hist[i]->GetXaxis()->FindBin(-hist[i]->GetXaxis()->GetBinCenter(j)))) / hist[i]->GetBinContent(j));
      
    hist[i]->SetStats(kFALSE);
    hist[i]->SetLineColor(colors[i]);
    hist[i]->GetXaxis()->SetRangeUser(0.001, 0.09);
    //hist[i]->GetYaxis()->SetRangeUser(5, hist[i]->GetMaximum() * 2);
    hist[i]->GetYaxis()->SetTitleOffset(1.3);
    hist[i]->GetXaxis()->SetTitle("#Delta#varphi (rad.)");
    hist[i]->Draw(((i == 0) ? "" : "SAME"));
    
    legend->AddEntry(hist[i], titles[i], "L");
  }
  
  legend->Draw();
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

void TrackCuts_Comparison_MC(char* histName, Int_t plotWhich = 0, const char* fileName = "correction_map.root", Bool_t mirror = kFALSE)
{
  // for the nsigmaplot it is needed to run with all cuts except the nsigmatovertex
  //    --> manually disable it in the run.C
  //
  // plotWhich: 0 = only before
  //            1 = both
  //            2 = only after
  //
  // mirror: kTRUE --> project negative values on the positive side
  

  file = TFile::Open(fileName);

  Int_t count = 0;
  Int_t colors[] = { 1, 2, 3, 4, 5, 6 };

  TLegend* legend = new TLegend(0.5, 0.7, 1, 1);
  TLegend* legend2 = new TLegend(0.4, 0.6, 1, 1);
  TLegend* legend3 = new TLegend(0.6, 0.5, 1, 0.7);

  TCanvas* c1 = new TCanvas(histName, histName, 800, 1200);
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
      
      if (mirror)
      {
        for (Int_t bin=1; bin<=hist->GetXaxis()->FindBin(-0.0001); bin++)
        {
          Int_t newBin = hist->GetXaxis()->FindBin(-hist->GetXaxis()->GetBinCenter(bin));
          if (bin != newBin)
          {
            hist->Fill(-hist->GetXaxis()->GetBinCenter(bin), hist->GetBinContent(bin));
            hist->SetBinContent(bin, 0);
          }
        }
      }
      
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

  //c1->SaveAs(Form("%s.png", histName));
}

void TrackCuts_Comparison_Data(char* histName, Int_t plotWhich, const char* fileName1, const char* fileName2, Bool_t mirror = kFALSE, const char* label1 = "file1", const char* label2 = "file2")
{
  // for the nsigmaplot it is needed to run with all cuts except the nsigmatovertex
  //    --> manually disable it in the run.C
  //
  // plotWhich: 0 = only before
  //            1 = both
  //            2 = only after
  //
  // mirror: kTRUE --> project negative values on the positive side
  

  Int_t count = 0;
  Int_t colors[] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };

  TLegend* legend = new TLegend(0.5, 0.7, 1, 1);
  legend->SetFillColor(0);
  TLegend* legend2 = new TLegend(0.4, 0.6, 1, 1);
  TLegend* legend3 = new TLegend(0.6, 0.5, 1, 0.7);

  TCanvas* c1 = new TCanvas(histName, histName, 600, 600);
  //TCanvas* c2 = new TCanvas("c2", "c2", 800, 600);
  //TCanvas* c3 = new TCanvas("c3", "c3", 800, 600);

  const char* folders2[] = { "before_cuts", "after_cuts" };
  Bool_t first = kTRUE;
  for (Int_t j = ((plotWhich < 2) ? 0 : 1); j < ((plotWhich > 0) ? 2 : 1); j++)
  {
    const char* folders1[] = { "esd_track_cuts", "esd_track_cuts_primaries", "esd_track_cuts_secondaries" };
    const char* names[] =    { "all", "primaries", "secondaries" };
    
    Float_t normalize[3];
    
    for (Int_t i = 0; i < 2; i++)
    {
      file = TFile::Open((i == 0) ? fileName1 : fileName2);
      
      for (Int_t k = 1; k < 3; k++)
      {
        TString folder;
        folder.Form("%s/%s/%s", folders1[k], folders2[j], histName);
        Printf("%s", folder.Data());
        TH1* hist = (TH1*) file->Get(folder);
        
        if (mirror)
        {
          for (Int_t bin=1; bin<=hist->GetXaxis()->FindBin(-0.0001); bin++)
          {
            Int_t newBin = hist->GetXaxis()->FindBin(-hist->GetXaxis()->GetBinCenter(bin));
            if (bin != newBin)
            {
              hist->Fill(-hist->GetXaxis()->GetBinCenter(bin), hist->GetBinContent(bin));
              hist->SetBinContent(bin, 0);
            }
          }
        }
      
        if (i == 0)
        {
          normalize[k] = hist->Integral();
        }
        else
          hist->Scale(normalize[k] / hist->Integral());
        
        legend->AddEntry(hist, Form("%s %s %s", (i == 0) ? label1 : label2, (k == 1) ? "primaries" : "secondaries", folders2[j]));
  
        c1->cd();
        hist->SetStats(0);
        hist->SetLineColor(colors[count]);
        hist->DrawCopy((count == 0) ? "" : "SAME");
  
        count++;
      }
    }

  }

  //c1->SetLogy();
  c1->SetGridx();
  c1->SetGridy();
  legend->Draw();
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

TH2* GetAcceptance(void* corr2d_void)
{
        corr2d = (AliCorrectionMatrix2D*) corr2d_void;
        corr_xy = (TH2*) corr2d->GetCorrectionHistogram()->Clone("acceptance");

        // fold in acceptance
        for (Int_t x=1; x<=corr_xy->GetNbinsX(); ++x)
                for (Int_t y=1; y<=corr_xy->GetNbinsY(); ++y)
                {
                        if (corr_xy->GetBinContent(x, y) > 1.5)
                                corr_xy->SetBinContent(x, y, 0);

                        if (corr_xy->GetBinContent(x, y) > 0)
                                corr_xy->SetBinContent(x, y, 1);

                        corr_xy->SetBinError(x, y, 0);
                }

        return corr_xy;
}

void ZeroOutsideAcceptance(TH2* acc, TH3* hist)
{
  for (Int_t x=0; x<=acc->GetNbinsX()+1; ++x)
    for (Int_t y=0; y<=acc->GetNbinsY()+1; ++y)
    {
      if (acc->GetBinContent(x, y) > 2 || acc->GetBinContent(x, y) == 0)
      {
        for (Int_t z=0; z<=hist->GetNbinsZ()+1; ++z)
        {
          hist->SetBinContent(x, y, z, 0);
          hist->SetBinError(x, y, z, 0);
        }
      }
    }
}

void DrawPhi()
{
  loadlibs();

  TFile::Open("correction_map.root");
  AlidNdEtaCorrection* dNdEtaCorrection = new AlidNdEtaCorrection("dndeta_correction", "dndeta_correction");
  if (!dNdEtaCorrection->LoadHistograms())
    return 0;

  TFile::Open("analysis_esd.root");
  dNdEtaAnalysis* fdNdEtaAnalysis = new dNdEtaAnalysis("dndetaTrVtx", "dndetaTrVtx");
  fdNdEtaAnalysis->LoadHistograms();

  // acc. map!
  //acc = GetAcceptance(dNdEtaCorrection->GetCorrection(1)->GetTrackCorrection()->Get2DCorrection("yx", 0, 1000));
  acc = dNdEtaCorrection->GetCorrection(1)->GetTrackCorrection()->Get2DCorrection("yx", 0, 1000)->GetCorrectionHistogram();
  //new TCanvas; acc->Draw("COLZ");

  histG = fdNdEtaAnalysis->GetData()->GetTrackCorrection()->GetGeneratedHistogram();
  ZeroOutsideAcceptance(acc, histG);
  //new TCanvas; histG->Project3D("yx")->Draw("COLZ");
  //histG->GetYaxis()->SetRangeUser(-0.9, 0.9);

  histM = fdNdEtaAnalysis->GetData()->GetTrackCorrection()->GetMeasuredHistogram();
  ZeroOutsideAcceptance(acc, histM);
  //histM->GetYaxis()->SetRangeUser(-0.9, 0.9);

  TFile::Open("analysis_mc.root");
  dNdEtaAnalysis* fdNdEtaAnalysis2 = new dNdEtaAnalysis("dndetaTrVtxMC", "dndetaTrVtxMC");
  fdNdEtaAnalysis2->LoadHistograms("dndetaTrVtx");

  histMC = fdNdEtaAnalysis2->GetData()->GetTrackCorrection()->GetMeasuredHistogram();
  ZeroOutsideAcceptance(acc, histMC);
  //new TCanvas; histMC->Project3D("yx2")->Draw("COLZ");

  //histG->GetZaxis()->SetRangeUser(1,2); histMC->GetZaxis()->SetRangeUser(1,2);
  new TCanvas; a = histG->Project3D("yx3"); a->Add(histMC->Project3D("yx4"), -1); a->Draw("COLZ");

  //histMC->GetYaxis()->SetRangeUser(-0.9, 0.9);

  c = new TCanvas;

  histG->GetXaxis()->SetRangeUser(-9.9, 9.9);
  histG->Project3D("z")->DrawCopy();

  histM->GetXaxis()->SetRangeUser(-9.9, 9.9);
  proj = histM->Project3D("z2");
  proj->SetLineColor(2);
  proj->DrawCopy("SAME");

  histMC->GetXaxis()->SetRangeUser(-9.9, 9.9);
  projMC = histMC->Project3D("z3");
  projMC->SetLineColor(4);
  projMC->DrawCopy("SAME");
}

void PrintEventStats(Int_t corrID = 3, const char* fileName = "correction_map.root", const char* dir = "dndeta_correction")
{
  loadlibs();

  /*
  TFile::Open("analysis_mc.root");
  fdNdEtaAnalysis = new dNdEtaAnalysis("dndetaNSD", "dndetaNSD");
  fdNdEtaAnalysis->LoadHistograms();
  trackHist = fdNdEtaAnalysis->GetData()->GetTrackCorrection()->GetGeneratedHistogram();
  eventHist = fdNdEtaAnalysis->GetData()->GetEventCorrection()->GetGeneratedHistogram();
  */

  TFile::Open(fileName);
  AlidNdEtaCorrection* dNdEtaCorrection = new AlidNdEtaCorrection(dir, dir);
  if (!dNdEtaCorrection->LoadHistograms())
    return;
  trackHist = dNdEtaCorrection->GetCorrection(corrID)->GetTrackCorrection()->GetGeneratedHistogram();
  eventHist = dNdEtaCorrection->GetCorrection(corrID)->GetEventCorrection()->GetGeneratedHistogram();

  trackHist->GetXaxis()->SetRange(0, trackHist->GetNbinsX()+1);
  trackHist->GetZaxis()->SetRange(0, trackHist->GetNbinsZ()+1);
  eta = trackHist->Project3D("y");

  events = eventHist->Integral(0, eventHist->GetNbinsX()+1, 0, eventHist->GetNbinsY()+1);

  eta->Scale(1.0 / events);

  Float_t avgN = eta->Integral(0, eta->GetNbinsX()+1);
  Printf("<N> = %f", avgN);

  eta->Scale(1.0 / eta->GetXaxis()->GetBinWidth(1));

  Printf("dndeta | eta = 0 is %f", (eta->GetBinContent(eta->FindBin(0.01)) + eta->GetBinContent(eta->FindBin(-0.01))) / 2);
  Printf("dndeta in |eta| < 0.5 is %f", eta->Integral(eta->FindBin(-0.49), eta->FindBin(0.49)) / (eta->FindBin(0.49) - eta->FindBin(-0.49) + 1));
  Printf("dndeta in |eta| < 1 is %f", eta->Integral(eta->FindBin(-0.99), eta->FindBin(0.99)) / (eta->FindBin(0.99) - eta->FindBin(-0.99) + 1));
  Printf("dndeta in |eta| < 1.5 is %f", eta->Integral(eta->FindBin(-1.49), eta->FindBin(1.49)) / (eta->FindBin(1.49) - eta->FindBin(-1.49) + 1));

  stats = (TH2*) gFile->Get("fEventStats");
  proj = stats->ProjectionX();
  gROOT->ProcessLine(".L PrintHist.C");
  PrintHist2D(stats);
  PrintHist(proj);
  
  Float_t ua5_SD = 0.153;
  Float_t ua5_DD = 0.080;
  Float_t ua5_ND = 0.767;
  
  Printf("+++ FRACTIONS +++");
  
  Printf("ND: %f", proj->GetBinContent(3) / proj->GetBinContent(1));
  Printf("SD: %f", proj->GetBinContent(4) / proj->GetBinContent(1));
  Printf("DD: %f", proj->GetBinContent(5) / proj->GetBinContent(1));
  
  Printf("+++ TRIGGER EFFICIENCIES +++");
  
  Printf("INEL = %.1f", 100. * (proj->GetBinContent(1) - stats->GetBinContent(1, 1) - stats->GetBinContent(1, 3)) / proj->GetBinContent(1));
  Printf("NSD  = %.1f", 100. * (proj->GetBinContent(2) - stats->GetBinContent(2, 1) - stats->GetBinContent(2, 3)) / proj->GetBinContent(2));
  
  
  Float_t trigND = 100. * (proj->GetBinContent(3) - stats->GetBinContent(3, 1) - stats->GetBinContent(3, 3)) / proj->GetBinContent(3);
  Float_t trigSD = 100. * (proj->GetBinContent(4) - stats->GetBinContent(4, 1) - stats->GetBinContent(4, 3)) / proj->GetBinContent(4);
  Float_t trigDD = 100. * (proj->GetBinContent(5) - stats->GetBinContent(5, 1) - stats->GetBinContent(5, 3)) / proj->GetBinContent(5);
  
  Printf("ND  = %.1f", trigND);
  Printf("SD  = %.1f", trigSD);
  Printf("DD  = %.1f", trigDD);
  
  Float_t trigINELUA5 = ua5_SD * trigSD + ua5_DD * trigDD + ua5_ND * trigND;
  Float_t trigNSDUA5  = (ua5_DD * trigDD + ua5_ND * trigND) / (ua5_DD + ua5_ND);
  Printf("INEL (UA5)  = %.1f", trigINELUA5);
  Printf("NSD (UA5)  = %.1f", trigNSDUA5);
  
  Printf("+++ VERTEX EFFICIENCIES +++");
  
  Printf("INEL = %.1f", 100. * (stats->GetBinContent(1, 3) + stats->GetBinContent(1, 4)) / proj->GetBinContent(1));
  Printf("NSD  = %.1f", 100. * (stats->GetBinContent(2, 3) + stats->GetBinContent(2, 4)) / proj->GetBinContent(2));
  
  Float_t vtxND = 100. * (stats->GetBinContent(3, 3) + stats->GetBinContent(3, 4)) / proj->GetBinContent(3);
  Float_t vtxSD = 100. * (stats->GetBinContent(4, 3) + stats->GetBinContent(4, 4)) / proj->GetBinContent(4);
  Float_t vtxDD = 100. * (stats->GetBinContent(5, 3) + stats->GetBinContent(5, 4)) / proj->GetBinContent(5);
  Printf("ND  = %.1f", vtxND);
  Printf("SD  = %.1f", vtxSD);
  Printf("DD  = %.1f", vtxDD);
  
  Float_t vtxINELUA5 = ua5_SD * vtxSD + ua5_DD * vtxDD + ua5_ND * vtxND;
  Float_t vtxNSDUA5  = (ua5_DD * vtxDD + ua5_ND * vtxND) / (ua5_DD + ua5_ND);
  Printf("INEL (UA5)  = %.1f", vtxINELUA5);
  Printf("NSD (UA5)  = %.1f", vtxNSDUA5);
  
  Printf("+++ TRIGGER + VERTEX EFFICIENCIES +++");
  
  Printf("INEL = %.1f", 100. * stats->GetBinContent(1, 4) / proj->GetBinContent(1));
  Printf("NSD  = %.1f", 100. * stats->GetBinContent(2, 4) / proj->GetBinContent(2));
  Printf("ND  = %.1f",  100. * stats->GetBinContent(3, 4) / proj->GetBinContent(3));
  Printf("SD  = %.1f",  100. * stats->GetBinContent(4, 4) / proj->GetBinContent(4));
  Printf("DD  = %.1f",  100. * stats->GetBinContent(5, 4) / proj->GetBinContent(5));
  
  
  
  for (Int_t i=7; i<=proj->GetNbinsX(); i++)
    if (proj->GetBinContent(i) > 0)
      Printf("bin %d (process type %d) = %.2f", i, (Int_t) proj->GetXaxis()->GetBinCenter(i), 100.0 * (proj->GetBinContent(i) - stats->GetBinContent(i, 1)) / proj->GetBinContent(i));
  
  //eta->Draw();
}

void TestAsymmetry()
{
  loadlibs();

  TFile* file2 = TFile::Open("analysis_mc.root");
  
  dNdEtaAnalysis* fdNdEtaAnalysis = new dNdEtaAnalysis("dndeta", "dndeta");
  fdNdEtaAnalysis->LoadHistograms();
  fdNdEtaAnalysis->Finish(0, 0, AlidNdEtaCorrection::kNone, "...");
  
  fdNdEtaAnalysis->GetdNdEtaHistogram(0)->DrawCopy();
  
  hist = (TH1*) fdNdEtaAnalysis->GetData()->GetTrackCorrection()->GetMeasuredHistogram();
  hist2 = (TH1*) hist->Clone("hist2");
  for (Int_t x=1; x<=hist->GetNbinsX(); x++)
    for (Int_t y=1; y<=hist->GetNbinsY(); y++)
      for (Int_t z=1; z<=hist->GetNbinsZ(); z++)
      {
        Printf("%d %d %d %d", x, y, z, hist->GetNbinsY() + 1 - y);
        hist->SetBinContent(x, y, z, hist2->GetBinContent(hist->GetNbinsX() / 2, TMath::Min(y, hist->GetNbinsY() + 1 - y), z));
      }
  
  hist = fdNdEtaAnalysis->GetData()->GetEventCorrection()->GetMeasuredHistogram();
  for (Int_t x=1; x<=hist->GetNbinsX(); x++)
    for (Int_t y=1; y<=hist->GetNbinsY(); y++)
      {
        //Printf("%d %d %d %d", x, y, z, hist->GetNbinsY() + 1 - y);
        hist->SetBinContent(x, y, hist->GetBinContent(hist->GetNbinsX() / 2, y));
      }
  
  fdNdEtaAnalysis->Finish(0, 0, AlidNdEtaCorrection::kNone, "...");
  fdNdEtaAnalysis->GetdNdEtaHistogram(0)->SetMarkerColor(2);
  fdNdEtaAnalysis->GetdNdEtaHistogram(0)->SetLineColor(2);
  fdNdEtaAnalysis->GetdNdEtaHistogram(0)->SetMarkerStyle(5);
  fdNdEtaAnalysis->GetdNdEtaHistogram(0)->DrawCopy("SAMEP");
}

void DeltaPhiFromPt(Float_t smearing = 0.005)
{
  loadlibs();

  TFile::Open("analysis_mc.root");
  hist = (TH1*) gFile->Get("dndeta_check_pt");
  
  dPhiHist = new TH1F("dPhiHist", ";#Delta phi", 400, -0.1, 0.1);
  dPhiHist2 = new TH1F("dPhiHist2", ";#Delta phi", 400, -0.1, 0.1);
  
  for (Int_t i=1; i<=hist->GetNbinsX(); i++)
  {
    Float_t pt = hist->GetBinCenter(i);
    Float_t deltaPhi = (0.076 - 0.039) / 2 / (pt / 0.15);
    
    if (smearing > 0)
    {
      gaus = new TF1("mygaus", "gaus(0)", -0.1, 0.1);
      gaus->SetParameters(1, -deltaPhi, smearing);
    
      dPhiHist->FillRandom("mygaus", hist->GetBinContent(i) / 2 * 1000);
    
      dPhiHist2->FillRandom("mygaus", hist->GetBinContent(i) / 2 * 1000);
      gaus->SetParameters(1, deltaPhi, smearing);
      dPhiHist2->FillRandom("mygaus", hist->GetBinContent(i) / 2 * 1000);
    }
    else
{
dPhiHist->Fill(deltaPhi, hist->GetBinContent(i) / 2);
dPhiHist2->Fill(deltaPhi, hist->GetBinContent(i) / 2);
dPhiHist2->Fill(-deltaPhi, hist->GetBinContent(i) / 2);
}
  }
  
  new TCanvas;
  dPhiHist->Draw();
  dPhiHist2->SetLineColor(2);
  dPhiHist2->Draw("SAME");
  gPad->SetLogy();
  
  TFile::Open("trackletsDePhi.root");
  //TFile::Open("tmp/correction_maponly-positive.root");
  //TFile::Open("tmp/correction_map.root");
  //tracklets = (TH1*) gFile->Get(Form("fDeltaPhi_%d", 0));
  tracklets = (TH1*) gFile->Get("DePhiPPTracklets");
  tracklets->Scale(1.0 / tracklets->GetMaximum() * dPhiHist->GetMaximum());
  tracklets->SetLineColor(4);
  tracklets->Draw("SAME");
}

void VertexDistributions()
{
  loadlibs();
  
  TFile::Open("correction_map.root");
  AlidNdEtaCorrection* dNdEtaCorrection = new AlidNdEtaCorrection("dndeta_correction", "dndeta_correction");
  if (!dNdEtaCorrection->LoadHistograms())
    return;
  
  all = dNdEtaCorrection->GetCorrection(AlidNdEtaCorrection::kINEL)->GetEventCorrection()->GetGeneratedHistogram()->ProjectionX("all");
  trigger = dNdEtaCorrection->GetCorrection(AlidNdEtaCorrection::kINEL)->GetEventCorrection()->GetMeasuredHistogram()->ProjectionX("trigger");
  vtx = dNdEtaCorrection->GetCorrection(AlidNdEtaCorrection::kVertexReco)->GetEventCorrection()->GetMeasuredHistogram()->ProjectionX("vtx");
 
  nottriggered = (TH1*) all->Clone("nottriggered");
  nottriggered->Add(trigger, -1);

  novertex = (TH1*) trigger->Clone("novertex");
  novertex->Add(vtx, -1);
  
  temphist = dNdEtaCorrection->GetCorrection(AlidNdEtaCorrection::kVertexReco)->GetEventCorrection()->GetMeasuredHistogram();
  highmult = temphist->ProjectionX("highmult", temphist->GetYaxis()->FindBin(10), temphist->GetNbinsY());
  //all = dNdEtaCorrection->GetCorrection(AlidNdEtaCorrection::kINEL)->GetEventCorrection()->GetGeneratedHistogram()->ProjectionX("all", temphist->GetYaxis()->FindBin(10), temphist->GetNbinsY());
 
  for (Int_t i=1; i<=trigger->GetNbinsX(); i++)
  {
    all->SetBinContent(i, all->GetBinContent(i) / all->GetBinWidth(i));
    trigger->SetBinContent(i, trigger->GetBinContent(i) / trigger->GetBinWidth(i));
    vtx->SetBinContent(i, vtx->GetBinContent(i) / vtx->GetBinWidth(i));
    nottriggered->SetBinContent(i, nottriggered->GetBinContent(i) / nottriggered->GetBinWidth(i));
    novertex->SetBinContent(i, novertex->GetBinContent(i) / novertex->GetBinWidth(i));
    highmult->SetBinContent(i, highmult->GetBinContent(i) / highmult->GetBinWidth(i));
  }

  new TCanvas;
  vtx->SetTitle("");
  vtx->SetStats(0);
  vtx->DrawCopy("HIST");

  all->Scale(1.0 / all->Integral());
  nottriggered->Scale(1.0 / nottriggered->Integral());
  novertex->Scale(1.0 / novertex->Integral());
  highmult->Scale(1.0 / highmult->Integral());

  new TCanvas;
  all->Draw("HIST");
  novertex->SetLineColor(2);
  novertex->Draw("HISTSAME");
  highmult->SetLineColor(3);
  highmult->Draw("HISTSAME");

  legend = new TLegend(0.5, 0.5, 0.8, 0.8);
  legend->SetFillColor(0);
  legend->AddEntry(all, "all");
  legend->AddEntry(novertex, "no vertex");
  legend->AddEntry(highmult, "mult > 10");
  legend->Draw();
  
  new TCanvas;
  trigger->Scale(1.0 / trigger->Integral());
  vtx->Scale(1.0 / vtx->Integral());
  
  trigger->Divide(vtx);
  
  trigger->Draw();
  //vtx->SetLineColor(2);
  //vtx->Draw("SAME");

  //temphist = dNdEtaCorrection->GetCorrection(AlidNdEtaCorrection::kVertexReco)->GetEventCorrection()->GetMeasuredHistogram();
  temphist = dNdEtaCorrection->GetCorrection(AlidNdEtaCorrection::kINEL)->GetEventCorrection()->GetGeneratedHistogram();
  //temphist = dNdEtaCorrection->GetCorrection(AlidNdEtaCorrection::kINEL)->GetEventCorrection()->GetMeasuredHistogram();
  
  temphist = (TH2*) gFile->Get("fTemp1");
  
  new TCanvas;
  legend = new TLegend(0.7, 0.7, 0.99, 0.99);
  legend->SetFillColor(0);
 
  Bool_t first = kTRUE; 
  for (Int_t i=0; i<20; i+=5)
  {
    highmult = temphist->ProjectionX("highmult", i+1, i+1+4);
    highmult->Rebin(10);
    Printf("%f", highmult->Integral());
    if (highmult->Integral() <= 0)
      continue;
  
    for (Int_t j=1; j<=trigger->GetNbinsX(); j++)
      highmult->SetBinContent(j, highmult->GetBinContent(j) / highmult->GetBinWidth(j));

    highmult->Scale(1.0 / highmult->Integral());
    highmult->SetLineColor((i/5)+1);
    highmult->GetYaxis()->SetRangeUser(0, 0.15);
    if (first)
    {
      highmult->DrawCopy();
      first = kFALSE;
    }
    else
      highmult->DrawCopy("SAME");
    legend->AddEntry(highmult->Clone(), Form("%d <= N <= %d", i, i+4));
  }
  legend->Draw();
 
}

void PlotPt1DCorrection()
{
  const char* files[] = { "field.root", "field_onlyprim.root", "nofield.root", "nofield_onlyprim.root" };
  const char* names[] = { "B: all", "B: primaries", "No B: all", "No B: primaries" };
  Int_t colors[] = { 1, 2, 3, 4 };
  
  loadlibs();
  
  dummy = new TH2F("dummy", ";p_{T};correction", 100, 0, 1.4, 100, 0.5, 3);
  dummy->SetStats(0);
  //dummy->GetYaxis()->SetTitleOffset(1.3);
  dummy->Draw();
  
  legend = new TLegend(0.48, 0.57, 0.88, 0.88);
  legend->SetFillColor(0);
  
  for (Int_t i=0; i<4; i++)
  {
    TFile::Open(files[i]);
    AlidNdEtaCorrection* dNdEtaCorrection = new AlidNdEtaCorrection("dndeta_correction", "dndeta_correction");
    if (!dNdEtaCorrection->LoadHistograms())
      return;
      
    hist = dNdEtaCorrection->GetTrack2ParticleCorrection()->GetTrackCorrection()->Get1DCorrectionHistogram("z", -9.9, 9.9, -0.79, 0.79);
    hist->SetLineColor(colors[i]);
    hist->SetLineWidth(2);
    hist->SetMarkerColor(colors[i]);
    hist->Draw("SAME");
    
    legend->AddEntry(hist, names[i], "L");
  }
  
  legend->Draw();
}

void FitDiamond()
{
  TFile::Open("analysis_esd_raw.root");
  
  hist = (TH3*) gFile->Get("vertex_check");
  
  gStyle->SetOptFit(1);
  
  TH1* proj[3];
  proj[0] = hist->ProjectionX();
  proj[1] = hist->ProjectionY();
  proj[2] = hist->ProjectionZ();
  
  for (Int_t i=0; i<3; i++)
  {
    c = new TCanvas;
    proj[i]->Draw();
    proj[i]->Fit("gaus");
    
    c->SaveAs(Form("FitDiamond_%d.png", i));
  }
}

void CompareDiamond(const char* mc, const char* data)
{
  TFile::Open(mc);
  
  hist = (TH3*) gFile->Get("vertex_check");
  
  gStyle->SetOptFit(1);
  
  TH1* proj[3];
  proj[0] = hist->ProjectionX("vertex_check_px");
  proj[1] = hist->ProjectionY("vertex_check_py");
  proj[2] = hist->ProjectionZ("vertex_check_pz");
  
  TFile::Open(data);
  
  hist = (TH3*) gFile->Get("vertex_check");
  
  TH1* proj2[3];
  proj2[0] = hist->ProjectionX("vertex_check_px2");
  proj2[1] = hist->ProjectionY("vertex_check_py2");
  proj2[2] = hist->ProjectionZ("vertex_check_pz2");

  for (Int_t i=0; i<3; i++)
  {
    CompareQualityHists(proj[i], proj2[i], 1, 1);
  }
}

void FitDiamondVsMult()
{
  TFile::Open("analysis_esd_raw.root");
  
  fVertexVsMult = (TH3*) gFile->Get("fVertexVsMult");
  fVertexVsMult->GetZaxis()->SetTitle("multiplicity");
  
  TH2* proj[2];
  proj[0] = (TH2*) fVertexVsMult->Project3D("xz");
  proj[1] = (TH2*) fVertexVsMult->Project3D("yz");
  
  gStyle->SetPadGridX(kTRUE);
  gStyle->SetPadGridY(kTRUE);
  
  Int_t max = 40;
  
  for (Int_t i=0; i<2; i++)
  {
    proj[i]->Rebin2D(4, 1);
    proj[i]->FitSlicesY();
    
    c = new TCanvas(Form("c_%d", i), Form("c_%d", i), 800, 400);
    c->Divide(2, 1);
    
    c->cd(1);
    hist = (TH1*) gROOT->FindObject(Form("fVertexVsMult_%sz_1", (i == 0) ? "x" : "y"));
    hist->GetXaxis()->SetRangeUser(0, max);
    hist->GetYaxis()->SetRangeUser(-0.4, 0.4);
    hist->Draw();
    
    c->cd(2);
    hist = (TH1*) gROOT->FindObject(Form("fVertexVsMult_%sz_2", (i == 0) ? "x" : "y"));
    hist->GetXaxis()->SetRangeUser(0, max);
    hist->GetYaxis()->SetRangeUser(0, 0.2);
    hist->Draw();
    
    c->SaveAs(Form("FitDiamondVsMult_%d.png", i));
  }
}

void CompareQualityHists(const char* fileName1, const char* fileName2, const char* plotName, Int_t rebin1 = 1, Int_t rebin2 = 1, const char* exec = 0)
{
  file1 = TFile::Open(fileName1);
  hist1 = (TH1*) file1->Get(plotName);
  
  file2 = TFile::Open(fileName2);
  hist2 = (TH1*) file2->Get(plotName);
  
  hist1->SetStats(0);
  
  Printf("Entries in histograms: %f %f", hist1->Integral(), hist2->Integral());
  
  if (exec)
  {
    hist1 = (TH1*) gROOT->ProcessLine(Form(exec, hist1, "hist1a"));
    hist2 = (TH1*) gROOT->ProcessLine(Form(exec, hist2, "hist2a"));
    hist1->Sumw2();
    hist2->Sumw2();
    Printf("Entries in histograms: %f %f", hist1->Integral(), hist2->Integral());
  }
  
  CompareQualityHists(hist1, hist2, rebin1, rebin2);
}

void CompareQualityHists(TH1* hist1, TH1* hist2, Int_t rebin1 = 1, Int_t rebin2 = 1)
{
  hist1->SetLineColor(1);
  hist2->SetLineColor(2);
 
  if (rebin1 != 0 && rebin2 != 0)
  { 
    hist1->Rebin(TMath::Abs(rebin1));
    hist2->Rebin(TMath::Abs(rebin2));
  }
  
  //hist2 = hist2->Rebin(hist1->GetNbinsX(), Form("%s_rebinned", hist2->GetName()), hist1->GetXaxis()->GetXbins()->GetArray());
      
      //hist1->Scale(1.0 / 0.83);

//hist1->GetXaxis()->SetRangeUser(0, 50);
/*  hist1->GetYaxis()->SetRangeUser(0.9, 1.2);
  hist1->Scale(1.0 / 0.808751);*/
  
  //hist1->Scale(1.0 / 1.24632);
  //hist1->Scale(1.0 / 1.23821);
  //hist1->Scale(1.0 / 1.26213);
  
  if (rebin1 > 0 && rebin2 > 0)
  {
    hist1->Scale(hist2->Integral() / hist1->Integral() * hist2->GetXaxis()->GetBinWidth(1) / hist1->GetXaxis()->GetBinWidth(1) / rebin1 * rebin2);
    
    //hist1->Scale(0.5);
    //hist2->Scale(0.5);
  }

  c = new TCanvas;
  if (strcmp(hist1->GetName(), "fDeltaTheta") == 0 || strcmp(hist1->GetName(), "fDeltaPhi") == 0 || strcmp(hist1->GetName(), "fMultVtx") == 0 || TString(hist1->GetName()).BeginsWith("vertex_check"))
    c->SetLogy();
  
  if (TString(hist1->GetName()).BeginsWith("fMultiplicityESD"))
  {
    c->SetLogy();
    loadlibs();
    AliPWG0Helper::NormalizeToBinWidth(hist1);
    AliPWG0Helper::NormalizeToBinWidth(hist2);
  }
  
  Printf("Means: %f %f %e", hist1->GetMean(), hist2->GetMean(), 1.0 - hist2->GetMean() / hist1->GetMean());
  
  //hist1->GetYaxis()->SetRangeUser(0.01, hist1->GetMaximum() * 1.3);
  hist1->DrawCopy("HISTE");
  hist2->DrawCopy("HISTE SAME");
  gPad->SetGridx();
  gPad->SetGridy();
  //gPad->SetLogy();
  c->SaveAs(Form("%s_1.png", hist1->GetName()));
  
  if (rebin1 == rebin2)
  {
    for (Int_t i=1; i<=hist1->GetNbinsX(); i++)
      if (hist1->GetBinContent(i) == 0 && hist2->GetBinContent(i) > 0 || hist1->GetBinContent(i) > 0 && hist2->GetBinContent(i) == 0)
        Printf("Inconsistent bin %d: %f %f", i, hist1->GetBinContent(i), hist2->GetBinContent(i));
  
    c2 = new TCanvas;
    hist1->GetYaxis()->SetRangeUser(0.5, 1.5);
    hist1->Divide(hist2);
    hist1->DrawCopy("HIST");
    gPad->SetGridx();
    gPad->SetGridy();
    c2->SaveAs(Form("%s_2.png", hist1->GetName()));
    
    /*
    for (Int_t i=1; i<=hist1->GetNbinsX(); i++)
      if (hist1->GetBinContent(i) > 0.9 && hist1->GetBinContent(i) < 1.1)
        hist1->SetBinContent(i, 0);
        
    new TCanvas;
    hist1->SetMarkerStyle(20);
    hist1->DrawCopy("P");
    */
  }
}

void DrawClustersVsTracklets()
{
  TFile::Open("analysis_esd_raw.root");
  
  hist = (TH2*) gFile->Get("fTrackletsVsClusters");
  
  c = new TCanvas("c", "c", 600, 600);
  c->SetRightMargin(0.05);
  c->SetTopMargin(0.05);
  
  hist->SetStats(0);
  hist->GetYaxis()->SetRangeUser(0, 400);
  hist->GetYaxis()->SetTitleOffset(1.3);
  hist->GetXaxis()->SetRangeUser(0, 30);
  hist->Draw("BOX");
  
  func = new TF1("func", "80 + x * 11", 0, 30);
  func->Draw("SAME");
  
  c->SaveAs("clusters_vs_tracklets.eps");
}

void VertexPlotBackgroundNote()
{
  TFile::Open("all.root");
  
  hist = (TH3*) gFile->Get("vertex_check");
  proj = (TH1*) hist->ProjectionZ()->Clone("all");
  proj->Rebin(2);
  
  proj->Draw();
  
  TFile::Open("analysis_esd_raw.root");
  hist = (TH3*) gFile->Get("vertex_check");
  proj = (TH1*) hist->ProjectionZ()->Clone("afterbg");
  proj->Rebin(2);
  
  proj->SetLineColor(2);
  proj->Draw("SAME");
}

void BackgroundAnalysis(const char* signal, const char* background)
{
  TFile::Open(signal);
  signalHist = (TH2*) gFile->Get("fTrackletsVsClusters");
  
  TFile::Open(background);
  backgroundHist = (TH2*) gFile->Get("fTrackletsVsClusters");
  
  Printf("For events with >= 1 tracklet:");
  
  func = new TF1("func", "[0] + x * 11", 0, 30);
  for (Int_t a = 50; a <= 100; a += 10)
  {
    func->SetParameter(0, a);
    
    Float_t signalCount = 0;
    Float_t backgroundCount = 0;
    for (Int_t x = 2; x <= signalHist->GetNbinsX(); x++)
    {
      signalCount += signalHist->Integral(x, x, signalHist->GetYaxis()->FindBin(func->Eval(signalHist->GetXaxis()->GetBinCenter(x))), signalHist->GetNbinsY());
      backgroundCount += backgroundHist->Integral(x, x, signalHist->GetYaxis()->FindBin(func->Eval(signalHist->GetXaxis()->GetBinCenter(x))), signalHist->GetNbinsY());
    }
    
    Float_t signalFraction = 100.0 * signalCount / signalHist->Integral(2, signalHist->GetNbinsX(), 1, signalHist->GetNbinsY());
    Float_t backgroundFraction = 100.0 * backgroundCount / backgroundHist->Integral(2, signalHist->GetNbinsX(), 1, signalHist->GetNbinsY());
    
    Printf("Cut at a = %d; Removed %.2f %% of the background (%.0f events); Removed %.2f %% of the signal", a, backgroundFraction, backgroundCount, signalFraction);
  }
}

void ZPhiPlots()
{
  TFile::Open("analysis_esd_raw.root");
  
  for (Int_t i=0; i<2; i++)
  {  
    hist = (TH2*) gFile->Get(Form("fZPhi_%d", i));
    
    c = new TCanvas;
    hist->SetStats(0);
    hist->Draw("COLZ");
    c->SaveAs(Form("ZPhi_%d.png", i));
  }
}

void DrawStats(Bool_t all = kFALSE)
{
  if (all)
  {
    Int_t count = 4;
    const char* list[] = { "CINT1B-ABCE-NOPF-ALL/spd", "CINT1A-ABCE-NOPF-ALL/spd", "CINT1C-ABCE-NOPF-ALL/spd", "CINT1-E-NOPF-ALL/spd" };
  }
  else
  {
    Int_t count = 1;
    const char* list[] = { "." };
  }
  
  for (Int_t i=0; i<count; i++)
  {
    TFile::Open(Form("%s/analysis_esd_raw.root", list[i]));
  
    hist = (TH2*) gFile->Get("fStats2");
    
    c = new TCanvas(list[i], list[i], 800, 600);
    gPad->SetBottomMargin(0.2);
    gPad->SetLeftMargin(0.2);
    gPad->SetRightMargin(0.2);
    hist->Draw("TEXT");
    hist->SetMarkerSize(2);
    //hist->GetYaxis()->SetRangeUser(0, 0);
    
    gROOT->Macro("increaseFonts.C");
  
    c->SaveAs(Form("%s/stats.png", list[i]));
  }
}

void CompareMCDataTrigger(const char* mcFile, const char* dataFile)
{
  TH1* stat[2];

  TFile::Open(mcFile);
  mc = (TH1*) gFile->Get("trigger_histograms_/fHistFiredBitsSPD");
  stat[0] = (TH1*) gFile->Get("fHistStatistics");
  
  TFile::Open(dataFile);
  data = (TH1*) gFile->Get("trigger_histograms_+CINT1B-ABCE-NOPF-ALL/fHistFiredBitsSPD");
  if (!data)
    data = (TH1*) gFile->Get("trigger_histograms_+CSMBB-ABCE-NOPF-ALL/fHistFiredBitsSPD");

  stat[1] = (TH1*) gFile->Get("fHistStatistics");
  
  CompareQualityHists(mc, data);
  
  for (Int_t i=0; i<2; i++)
  {
    Float_t total = stat[i]->GetBinContent(stat[i]->GetXaxis()->FindBin("Trigger class"), 1);
    Float_t spd = stat[i]->GetBinContent(stat[i]->GetXaxis()->FindBin("FO >= 2"), 1);
    Float_t v0A = stat[i]->GetBinContent(stat[i]->GetXaxis()->FindBin("V0A"), 1);
    Float_t v0C = stat[i]->GetBinContent(stat[i]->GetXaxis()->FindBin("V0C"), 1);
    
    Printf("%s:\nSPD / V0A: %.3f\nSPD / V0C: %.3f\nV0A / V0C: %.3f", (i == 0) ? "MC  " : "Data", spd / v0A, spd / v0C, v0A / v0C);
    Printf("SPD / Total: %.3f\nV0A / Total: %.3f\nV0C / Total: %.3f\n", spd / total, v0A / total, v0C / total);
  }
}

void CompareMCDatadNdEta(const char* mcFile, const char* dataFile)
{
  //CompareQualityHists(mcFile, dataFile, "fEtaPhi", 4, 4, "((TH2*)%p)->ProjectionY(\"%s\", 1, 40)");
  //CompareQualityHists(mcFile, dataFile, "fEtaPhi", 4, 4, "((TH2*)%p)->ProjectionY(\"%s\", 41, 80)");

  CompareQualityHists(mcFile, dataFile, "fEtaPhi", 1, 1, "((TH2*)%p)->ProjectionX(\"%s\", 271, 360)");
}

void TrigVsTrigVtx(const char* fileName = "correction_map.root", const char* dirName = "dndeta_correction")
{
  loadlibs();
  if (!TFile::Open(fileName))
    return;

  AlidNdEtaCorrection* dNdEtaCorrection = new AlidNdEtaCorrection(dirName, dirName);
  if (!dNdEtaCorrection->LoadHistograms())
    return;
  
  TH2* eTrig =    dNdEtaCorrection->GetVertexRecoCorrection()->GetEventCorrection()->GetGeneratedHistogram();
  TH2* eTrigVtx = dNdEtaCorrection->GetVertexRecoCorrection()->GetEventCorrection()->GetMeasuredHistogram();
  
  eTrig_proj = eTrig->ProjectionY("y1", eTrig->GetYaxis()->FindBin(-9.9), eTrig->GetYaxis()->FindBin(9.9));
  eTrigVtx_proj = eTrigVtx->ProjectionY("y2", eTrig->GetYaxis()->FindBin(-9.9), eTrig->GetYaxis()->FindBin(9.9));
  
  new TCanvas;
  eTrig_proj->Draw();
  eTrig_proj->GetXaxis()->SetRangeUser(0, 20);
  eTrigVtx_proj->SetLineColor(2);
  eTrigVtx_proj->Draw("SAME");
  
  gPad->SetLogy();
}

void PrintAverageNSDCorrectionFactors()
{
  // factors estimated from MC, can be slighly different with data b/c correction is applies as function of measured multiplicity

  loadlibs();

  if (!TFile::Open("correction_mapprocess-types.root"))
    return;

  AlidNdEtaCorrection* dNdEtaCorrection = new AlidNdEtaCorrection("dndeta_correction_ND", "dndeta_correction_ND");
  if (!dNdEtaCorrection->LoadHistograms())
    return;

  AlidNdEtaCorrection* dNdEtaCorrection2 = new AlidNdEtaCorrection("dndeta_correction_DD", "dndeta_correction_DD");
  if (!dNdEtaCorrection2->LoadHistograms())
    return;
    
  // for scaling factors see drawSystematics.C; GetRelativeFractions()
  // 900 GeV
  //dNdEtaCorrection->Scale(1.06);
  //dNdEtaCorrection->Add(dNdEtaCorrection2, 9.5 / 12.3);
  // 2.36 TeV
  dNdEtaCorrection->Scale(1.036);
  dNdEtaCorrection->Add(dNdEtaCorrection2, 0.075 * 1.43 / 0.127);
  
  Printf("event adding: %f", dNdEtaCorrection->GetTriggerBiasCorrectionINEL()->GetEventCorrection()->GetGeneratedHistogram()->Integral() / dNdEtaCorrection->GetTriggerBiasCorrectionINEL()->GetEventCorrection()->GetMeasuredHistogram()->Integral());
  
  Printf("track adding: %f", dNdEtaCorrection->GetTriggerBiasCorrectionINEL()->GetTrackCorrection()->GetGeneratedHistogram()->Integral() / dNdEtaCorrection->GetTriggerBiasCorrectionINEL()->GetTrackCorrection()->GetMeasuredHistogram()->Integral());
  
  AlidNdEtaCorrection* dNdEtaCorrection3 = new AlidNdEtaCorrection("dndeta_correction_SD", "dndeta_correction_SD");
  if (!dNdEtaCorrection3->LoadHistograms())
    return;

  // 900 GeV
  //dNdEtaCorrection3->Scale(0.153 / 0.189);
  // 2.36 TeV
  dNdEtaCorrection3->Scale(0.159 / 0.166);
  dNdEtaCorrection->Add(dNdEtaCorrection3);

  Printf("event subtraction: %f", dNdEtaCorrection3->GetTriggerBiasCorrectionINEL()->GetEventCorrection()->GetMeasuredHistogram()->Integral() / dNdEtaCorrection->GetTriggerBiasCorrectionINEL()->GetEventCorrection()->GetMeasuredHistogram()->Integral());
  
  Printf("track subtraction: %f", dNdEtaCorrection3->GetTriggerBiasCorrectionINEL()->GetTrackCorrection()->GetMeasuredHistogram()->Integral() / dNdEtaCorrection->GetTriggerBiasCorrectionINEL()->GetTrackCorrection()->GetMeasuredHistogram()->Integral());
  
  dNdEtaCorrection->GetTriggerBiasCorrectionNSD()->PrintInfo(0.0);
}
    
