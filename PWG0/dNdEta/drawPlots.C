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

void SetRanges(TAxis* axis)
{
  if (strcmp(axis->GetTitle(), "#eta") == 0)
    axis->SetRangeUser(-1.7999, 1.7999);
  if (strcmp(axis->GetTitle(), "p_{T} [GeV/c]") == 0)
    axis->SetRangeUser(0, 9.9999);
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

  gPad->SetGridx();
  gPad->SetGridy();
}

void dNdEta(Bool_t onlyESD = kFALSE)
{
  TFile* file = TFile::Open("analysis_esd.root");
  TH1* histESD = (TH1*) file->Get("dndeta/dndeta_dNdEta_corrected_2");
  TH1* histESDNoPt = (TH1*) file->Get("dndeta/dndeta_dNdEta_2");
  TH1* histESDMB = (TH1*) file->Get("dndeta_mb/dndeta_mb_dNdEta_corrected_2");
  TH1* histESDMBVtx = (TH1*) file->Get("dndeta_mbvtx/dndeta_mbvtx_dNdEta_corrected_2");

  TCanvas* canvas = new TCanvas("dNdEta1", "dNdEta1", 500, 500);

  Prepare1DPlot(histESD);
  Prepare1DPlot(histESDNoPt);
  Prepare1DPlot(histESDMB);
  Prepare1DPlot(histESDMBVtx);

  histESD->SetLineColor(0);
  histESDMB->SetLineColor(0);
  histESDMBVtx->SetLineColor(0);

  histESD->SetMarkerColor(kRed);
  histESDMB->SetMarkerColor(kBlue);
  histESDMBVtx->SetMarkerColor(103);

  histESD->SetMarkerStyle(20);
  histESDMB->SetMarkerStyle(21);
  histESDMBVtx->SetMarkerStyle(22);

  TH2F* dummy = new TH2F("dummy", "", 100, -1.5, 1.5, 100, 0, histESDMBVtx->GetMaximum() * 1.1);
  Prepare1DPlot(dummy);
  dummy->SetStats(kFALSE);
  dummy->SetXTitle("#eta");
  dummy->SetYTitle("dN_{ch}/d#eta");
  dummy->GetYaxis()->SetTitleOffset(1);

  Float_t etaLimit = 1.1999;

  histESDMBVtx->GetXaxis()->SetRangeUser(-etaLimit, etaLimit);
  histESDMB->GetXaxis()->SetRangeUser(-etaLimit, etaLimit);
  histESD->GetXaxis()->SetRangeUser(-etaLimit, etaLimit);
  histESDNoPt->GetXaxis()->SetRangeUser(-etaLimit, etaLimit);

  dummy->DrawCopy();
  histESDMBVtx->Draw("SAME");
  histESDMB->Draw("SAME");
  histESD->Draw("SAME");

  canvas->SaveAs("dNdEta1.gif");
  canvas->SaveAs("dNdEta1.eps");

  if (onlyESD)
    return;

  TFile* file2 = TFile::Open("analysis_mc.root");
  TH1* histMC = (TH1*) file2->Get("dndeta/dndeta_dNdEta_corrected_2")->Clone("cloned");

  gSystem->Load("libPWG0base");
  dNdEtaAnalysis* fdNdEtaAnalysis = new dNdEtaAnalysis("dndeta", "dndeta");
  fdNdEtaAnalysis->LoadHistograms();
  fdNdEtaAnalysis->Finish(0, 0.3);
  TH1* histMCPtCut = fdNdEtaAnalysis->GetdNdEtaHistogram(2);

  TCanvas* canvas2 = new TCanvas("dNdEta2", "dNdEta2", 500, 500);

  Prepare1DPlot(histMC);
  Prepare1DPlot(histMCPtCut);

  histMC->SetLineColor(kBlue);
  histMCPtCut->SetLineColor(104);
  histESDNoPt->SetLineColor(102);

  TH2* dummy2 = (TH2F*) dummy->Clone("dummy2");
  Prepare1DPlot(dummy2);
  dummy2->GetYaxis()->SetRangeUser(0, histESD->GetMaximum() * 1.1);

  dummy2->DrawCopy();
  histMC->Draw("SAME");
//  histMC->Draw();
  histESD->Draw("SAME");
  histESDNoPt->Draw("SAME");
  histMCPtCut->Draw("SAME");

  canvas2->SaveAs("dNdEta2.gif");
  canvas2->SaveAs("dNdEta2.eps");

  TCanvas* canvas3 = new TCanvas("dNdEta", "dNdEta", 700, 500);
  //InitPad();
  gPad->SetBottomMargin(0.12);

  TLegend* legend = new TLegend(0.35, 0.2, 0.6, 0.4);
  legend->SetFillColor(0);
  legend->AddEntry(histESDMBVtx, "triggered, vertex");
  legend->AddEntry(histESDMB, "triggered");
  legend->AddEntry(histESD, "all events");
  legend->AddEntry(histMC, "MC prediction");

  dummy->DrawCopy();
  histESDMBVtx->Draw("SAME");
  histESDMB->Draw("SAME");
  histESD->Draw("SAME");
  histMC->Draw("SAME");

  legend->Draw();

  canvas3->SaveAs("dNdEta.gif");
  canvas3->SaveAs("dNdEta.eps");

  TH1* ratio = (TH1*) histMC->Clone("ratio");
  TH1* ratioNoPt = (TH1*) histMCPtCut->Clone("ratioNoPt");

  ratio->Divide(histESD);
  ratioNoPt->Divide(histESDNoPt);

  TCanvas* canvas4 = new TCanvas("ratio", "ratio", 700, 500);

  ratio->GetXaxis()->SetRangeUser(-0.7999, 0.7999);

  ratio->SetLineColor(1);
  ratioNoPt->SetLineColor(2);

  ratio->Draw();
  ratioNoPt->Draw("SAME");

  TLegend* legend = new TLegend(0.6, 0.7, 0.95, 0.9);
  legend->SetFillColor(0);
  legend->AddEntry(ratio, "mc/esd");
  legend->AddEntry(ratioNoPt, "mc/esd, not pt cut off corrected");
  legend->Draw();
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

void ptCutoff()
{
  gSystem->Load("libPWG0base");

  AlidNdEtaCorrection* dNdEtaCorrection = new AlidNdEtaCorrection("dndeta_correction", "dndeta_correction");
  dNdEtaCorrection->LoadHistograms("correction_map.root","dndeta_correction");

  dNdEtaCorrection->GetMeasuredFraction(0.3, -100, kTRUE);

  TH1* hist = dynamic_cast<TH1*> (gROOT->FindObject("gene_dndeta_correction_nTrackToNPart_pt")->Clone("ptcutoff"));

  hist->GetXaxis()->SetRangeUser(0, 0.9999);
  hist->SetMinimum(0);

  hist->SetTitle("Generated Particles");
  Prepare1DPlot(hist);

  TCanvas* canvas = new TCanvas("ptCutoff", "ptCutoff", 700, 500);
  hist->DrawCopy();

  TLine* line = new TLine(0.3, 0 - hist->GetMaximum() * 0, 0.3, hist->GetMaximum() * 1.1);
  line->SetLineWidth(3);
  line->SetLineColor(kRed);
  line->Draw();

  canvas->SaveAs("ptCutoff.gif");
  canvas->SaveAs("ptCutoff.eps");

  TH1F* factor = new TH1F("factor", ";#eta;correction factor", 10, -1, 1.000001);
  for (Float_t eta = -0.9; eta<1; eta += 0.2)
    factor->Fill(eta, 1.0 / dNdEtaCorrection->GetMeasuredFraction(0.3, eta, kFALSE));

  TCanvas* canvas = new TCanvas("ptCutoff_factor", "ptCutoff_factor", 700, 500);
  InitPad();

  Prepare1DPlot(factor);
  factor->GetYaxis()->SetRangeUser(1, 2);
  factor->GetYaxis()->SetTitleOffset(1);
  factor->Draw();

  canvas->SaveAs("ptCutoff_factor.eps");
}

void TriggerBiasVtxRecon(const char* fileName = "correction_map.root", const char* folder = "dndeta_correction")
{
  TFile* file = TFile::Open(fileName);

  TH2* corrTrigger = dynamic_cast<TH2*> (file->Get(Form("%s/corr_%s_trigger", folder, folder)));
  TH2* corrVtx = dynamic_cast<TH2*> (file->Get(Form("%s/corr_%s_vtxReco", folder, folder)));

  Prepare2DPlot(corrTrigger);
  corrTrigger->SetTitle("a) Trigger bias correction");

  Prepare2DPlot(corrVtx);
  corrVtx->SetTitle("b) Vertex reconstruction correction");

  corrTrigger->GetYaxis()->SetTitle("Multiplicity");
  corrVtx->GetYaxis()->SetTitle("Multiplicity");

  TCanvas* canvas = new TCanvas("TriggerBiasVtxRecon", "TriggerBiasVtxRecon", 1000, 500);
  canvas->Divide(2, 1);

  canvas->cd(1);
  InitPadCOLZ();
  corrTrigger->DrawCopy("COLZ");

  canvas->cd(2);
  InitPadCOLZ();
  corrVtx->DrawCopy("COLZ");

  canvas->SaveAs(Form("TriggerBiasVtxRecon_%d.gif", gMax));
  canvas->SaveAs(Form("TriggerBiasVtxRecon_%d.eps", gMax));

  canvas = new TCanvas("TriggerBiasVtxReconZoom", "TriggerBiasVtxReconZoom", 1000, 500);
  canvas->Divide(2, 1);

  corrTrigger->GetYaxis()->SetRangeUser(0, 5);
  corrVtx->GetYaxis()->SetRangeUser(0, 5);

  canvas->cd(1);
  InitPadCOLZ();
  corrTrigger->DrawCopy("COLZ");

  canvas->cd(2);
  InitPadCOLZ();
  corrVtx->DrawCopy("COLZ");

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
  dNdEtaCorrection->LoadHistograms(fileName, folderName);

  TH1* hist = dNdEtaCorrection->GetTriggerCorrection()->Get1DCorrection("x");

  TCanvas* canvas = new TCanvas("TriggerBias1D", "TriggerBias1D", 500, 500);
  InitPad();

  Prepare1DPlot(hist);
  hist->SetTitle("");
  hist->GetYaxis()->SetTitle("correction factor");
  hist->GetYaxis()->SetRangeUser(1, 1.5);
  hist->GetYaxis()->SetTitleOffset(1.6);
  hist->Draw();

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
  dNdEtaCorrection->LoadHistograms(fileName, folderName);

  TH1* hist = dNdEtaCorrection->GetVertexRecoCorrection()->Get1DCorrection("x");

  TCanvas* canvas = new TCanvas("VtxRecon1D", "VtxRecon1D", 500, 500);
  InitPad();

  Prepare1DPlot(hist);
  hist->SetTitle("");
  hist->GetYaxis()->SetTitle("correction factor");
  hist->GetYaxis()->SetRangeUser(1, 1.8);
  hist->GetYaxis()->SetTitleOffset(1.6);
  hist->Draw();

  canvas->SaveAs("VtxRecon1D.eps");
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

void Track2Particle1DCreatePlots(const char* fileName = "correction_map.root", const char* folderName = "dndeta_correction", Float_t upperPtLimit = 10)
{
  TFile::Open(fileName);
  AlidNdEtaCorrection* dNdEtaCorrection = new AlidNdEtaCorrection(folderName, folderName);
  dNdEtaCorrection->LoadHistograms(fileName, folderName);

  TH3F* gene = dNdEtaCorrection->GetTrack2ParticleCorrection()->GetGeneratedHistogram();
  TH3F* meas = dNdEtaCorrection->GetTrack2ParticleCorrection()->GetMeasuredHistogram();

  gene->GetZaxis()->SetRangeUser(0.3, upperPtLimit);
  meas->GetZaxis()->SetRangeUser(0.3, upperPtLimit);
  gene->GetYaxis()->SetRangeUser(-0.8, 0.8);
  meas->GetYaxis()->SetRangeUser(-0.8, 0.8);
  AliPWG0Helper::CreateDividedProjections(gene, meas, "x", kTRUE);
  gene->GetYaxis()->SetRange(0, 0);
  meas->GetYaxis()->SetRange(0, 0);

  gene->GetXaxis()->SetRangeUser(-10, 10);
  meas->GetXaxis()->SetRangeUser(-10, 10);
  AliPWG0Helper::CreateDividedProjections(gene, meas, "y", kTRUE);
  gene->GetZaxis()->SetRange(0, 0);
  meas->GetZaxis()->SetRange(0, 0);

  gene->GetYaxis()->SetRangeUser(-0.8, 0.8);
  meas->GetYaxis()->SetRangeUser(-0.8, 0.8);
  AliPWG0Helper::CreateDividedProjections(gene, meas, "z", kTRUE);
}

void Track2Particle1D(const char* fileName = "correction_map.root", const char* folder = "dndeta_correction", Float_t upperPtLimit = 9.9)
{
  gSystem->Load("libPWG0base");

  Track2Particle1DCreatePlots(fileName, folder, upperPtLimit);

  TH1* corrX = dynamic_cast<TH1*> (gROOT->FindObject(Form("gene_%s_nTrackToNPart_x_div_meas_%s_nTrackToNPart_x", folder, folder)));
  TH1* corrY = dynamic_cast<TH1*> (gROOT->FindObject(Form("gene_%s_nTrackToNPart_y_div_meas_%s_nTrackToNPart_y", folder, folder)));
  TH1* corrZ = dynamic_cast<TH1*> (gROOT->FindObject(Form("gene_%s_nTrackToNPart_z_div_meas_%s_nTrackToNPart_z", folder, folder)));

  Prepare1DPlot(corrX);
  Prepare1DPlot(corrY);
  Prepare1DPlot(corrZ);

  //corrX->SetTitle("a) z projection");
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
  corrY->DrawCopy();

  canvas->cd(3);
  InitPad();
  corrZ->DrawCopy();

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

  Track2Particle1DCreatePlots("correction_maponly-positive.root", "dndeta_correction", upperPtLimit);

  TH1* posX = dynamic_cast<TH1*> (gROOT->FindObject("gene_nTrackToNPart_x_div_meas_nTrackToNPart_x")->Clone("pos_x"));
  TH1* posY = dynamic_cast<TH1*> (gROOT->FindObject("gene_nTrackToNPart_y_div_meas_nTrackToNPart_y")->Clone("pos_y"));
  TH1* posZ = dynamic_cast<TH1*> (gROOT->FindObject("gene_nTrackToNPart_z_div_meas_nTrackToNPart_z")->Clone("pos_z"));

  Track2Particle1DCreatePlots("correction_maponly-negative.root", "dndeta_correction", upperPtLimit);

  TH1* negX = dynamic_cast<TH1*> (gROOT->FindObject("gene_nTrackToNPart_x_div_meas_nTrackToNPart_x")->Clone("neg_x"));
  TH1* negY = dynamic_cast<TH1*> (gROOT->FindObject("gene_nTrackToNPart_y_div_meas_nTrackToNPart_y")->Clone("neg_y"));
  TH1* negZ = dynamic_cast<TH1*> (gROOT->FindObject("gene_nTrackToNPart_z_div_meas_nTrackToNPart_z")->Clone("neg_z"));

  //printf("%f %f %f %f\n", posX->GetBinContent(20), posX->GetBinError(20), negX->GetBinContent(20), negX->GetBinError(20));

  posX->Divide(negX);
  posY->Divide(negY);
  posZ->Divide(negZ);

  //printf("%f %f\n", posX->GetBinContent(20), posX->GetBinError(20));

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

  TCanvas* canvas = new TCanvas("CompareTrack2Particle1D", "CompareTrack2Particle1D", 1200, 400);
  canvas->Divide(3, 1);

  canvas->cd(1);
  InitPad();
  posX->Draw();

  canvas->cd(2);
  InitPad();
  posY->Draw();

  canvas->cd(3);
  InitPad();
  posZ->Draw();

  canvas->SaveAs(Form("CompareTrack2Particle1D_%f.gif", upperPtLimit));
  canvas->SaveAs(Form("CompareTrack2Particle1D_%f.eps", upperPtLimit));
}

void Track2Particle2DCreatePlots(const char* fileName = "correction_map.root")
{
  TFile::Open(fileName);
  AlidNdEtaCorrection* dNdEtaCorrection = new AlidNdEtaCorrection("dndeta_correction", "dndeta_correction");
  dNdEtaCorrection->LoadHistograms(fileName, "dndeta_correction");

  TH3F* gene = dNdEtaCorrection->GetTrack2ParticleCorrection()->GetGeneratedHistogram();
  TH3F* meas = dNdEtaCorrection->GetTrack2ParticleCorrection()->GetMeasuredHistogram();

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

void Track2Particle2D(const char* fileName = "correction_map.root", const char* folder = "dndeta_correction")
{
  gSystem->Load("libPWG0base");

  Track2Particle2DCreatePlots(fileName);

  TH2* corrYX = dynamic_cast<TH2*> (gROOT->FindObject(Form("gene_%s_nTrackToNPart_yx_div_meas_%s_nTrackToNPart_yx", folder, folder)));
  TH2* corrZX = dynamic_cast<TH2*> (gROOT->FindObject(Form("gene_%s_nTrackToNPart_zx_div_meas_%s_nTrackToNPart_zx", folder, folder)));
  TH2* corrZY = dynamic_cast<TH2*> (gROOT->FindObject(Form("gene_%s_nTrackToNPart_zy_div_meas_%s_nTrackToNPart_zy", folder, folder)));

  /* this reads them from the file
  TH2* corrYX = dynamic_cast<TH2*> (file->Get("dndeta_correction/gene_nTrackToNPart_yx_div_meas_nTrackToNPart_yx"));
  TH2* corrZX = dynamic_cast<TH2*> (file->Get("dndeta_correction/gene_nTrackToNPart_zx_div_meas_nTrackToNPart_zx"));
  TH2* corrZY = dynamic_cast<TH2*> (file->Get("dndeta_correction/gene_nTrackToNPart_zy_div_meas_nTrackToNPart_zy"));*/

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

  canvas->SaveAs(Form("Track2Particle2D_%d.gif", gMax));
  canvas->SaveAs(Form("Track2Particle2D_%d.eps", gMax));
}

void CompareTrack2Particle2D()
{
  gSystem->Load("libPWG0base");

  Track2Particle2DCreatePlots("correction_maponly-positive.root");

  TH2* posYX = dynamic_cast<TH2*> (gROOT->FindObject("gene_nTrackToNPart_yx_div_meas_nTrackToNPart_yx")->Clone("pos_yx"));
  TH2* posZX = dynamic_cast<TH2*> (gROOT->FindObject("gene_nTrackToNPart_zx_div_meas_nTrackToNPart_zx")->Clone("pos_zx"));
  TH2* posZY = dynamic_cast<TH2*> (gROOT->FindObject("gene_nTrackToNPart_zy_div_meas_nTrackToNPart_zy")->Clone("pos_zy"));

  Track2Particle2DCreatePlots("correction_maponly-negative.root");

  TH2* negYX = dynamic_cast<TH2*> (gROOT->FindObject("gene_nTrackToNPart_yx_div_meas_nTrackToNPart_yx")->Clone("neg_yx"));
  TH2* negZX = dynamic_cast<TH2*> (gROOT->FindObject("gene_nTrackToNPart_zx_div_meas_nTrackToNPart_zx")->Clone("neg_zx"));
  TH2* negZY = dynamic_cast<TH2*> (gROOT->FindObject("gene_nTrackToNPart_zy_div_meas_nTrackToNPart_zy")->Clone("neg_zy"));

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
