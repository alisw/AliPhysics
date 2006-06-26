/* $Id$ */

Int_t gMax = 5;

void drawPlots()
{
  drawPlots(5);
  drawPlots(2);
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

void dNdEta()
{
  TFile* file = TFile::Open("analysis_esd.root");
  TH1* histESD = file->Get("dndeta/dndeta_dNdEta_corrected_2");
  TH1* histESDNoPt = file->Get("dndeta/dndeta_dNdEta_2");
  TH1* histESDMB = file->Get("dndeta_mb/dndeta_mb_dNdEta_corrected_2");
  TH1* histESDMBVtx = file->Get("dndeta_mbvtx/dndeta_mbvtx_dNdEta_corrected_2");

  TCanvas* canvas = new TCanvas("dNdEta1", "dNdEta1", 500, 500);

  Prepare1DPlot(histESD);
  Prepare1DPlot(histESDNoPt);
  Prepare1DPlot(histESDMB);
  Prepare1DPlot(histESDMBVtx);

  histESD->SetLineColor(kRed);
  histESDMB->SetLineColor(kBlue);
  histESDMBVtx->SetLineColor(103);

  TH2F* dummy = new TH2F("dummy", "", 100, -1.5, 1.5, 100, 0, histESDMBVtx->GetMaximum() * 1.1);
  dummy->SetStats(kFALSE);
  dummy->SetXTitle("#eta");
  dummy->SetYTitle("dN/d#eta");

  histESDMBVtx->GetXaxis()->SetRangeUser(-0.7999, 0.7999);
  histESDMB->GetXaxis()->SetRangeUser(-0.7999, 0.7999);
  histESD->GetXaxis()->SetRangeUser(-0.7999, 0.7999);
  histESDNoPt->GetXaxis()->SetRangeUser(-0.7999, 0.7999);

  dummy->Draw();
  histESDMBVtx->Draw("SAME");
  histESDMB->Draw("SAME");
  histESD->Draw("SAME");

  canvas->SaveAs("dNdEta1.gif");

  TFile* file2 = TFile::Open("analysis_mc.root");
  TH1* histMC = file2->Get("dndeta/dndeta_dNdEta_corrected_2")->Clone("cloned");

  gSystem->Load("libPWG0base");
  fdNdEtaAnalysis = new dNdEtaAnalysis("dndeta", "dndeta");
  fdNdEtaAnalysis->LoadHistograms();
  fdNdEtaAnalysis->Finish(0, 0.3);
  TH1* histMCPtCut = fdNdEtaAnalysis->GetdNdEtaHistogram(2);

  TCanvas* canvas2 = new TCanvas("dNdEta2", "dNdEta2", 500, 500);

  Prepare1DPlot(histMC);
  Prepare1DPlot(histMCPtCut);

  histMC->SetLineColor(kBlue);
  histMCPtCut->SetLineColor(104);
  histESDNoPt->SetLineColor(102);

  TH2* dummy2 = dummy->Clone("dummy2");
  dummy2->GetYaxis()->SetRangeUser(0, histESD->GetMaximum() * 1.1);

  dummy2->Draw();
  histMC->Draw("SAME");
//  histMC->Draw();
  histESD->Draw("SAME");
  histESDNoPt->Draw("SAME");
  histMCPtCut->Draw("SAME");

  canvas2->SaveAs("dNdEta2.gif");
}

void ptSpectrum()
{
  TFile* file = TFile::Open("analysis_esd.root");
  TH1* histESD = file->Get("dndeta/dndeta_pt");

  TFile* file2 = TFile::Open("analysis_mc.root");
  TH1* histMC = file2->Get("dndeta/dndeta_pt");

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
}

void ptCutoff()
{
  gSystem->Load("libPWG0base");

  AlidNdEtaCorrection* dNdEtaCorrection = new AlidNdEtaCorrection();
  dNdEtaCorrection->LoadHistograms("correction_map.root","dndeta_correction");

  dNdEtaCorrection->GetMeasuredFraction(0.3, -1, kTRUE);

  TH1* hist = dynamic_cast<TH1*> (gROOT->FindObject("gene_nTrackToNPart_pt"));

  hist->GetXaxis()->SetRangeUser(0, 0.9999);
  hist->SetMinimum(0);

  hist->SetTitle("Generated Particles");
  Prepare1DPlot(hist);

  TCanvas* canvas = new TCanvas("ptCutoff", "ptCutoff", 500, 500);
  hist->Draw();

  TLine* line = new TLine(0.3, 0 - hist->GetMaximum() * 0, 0.3, hist->GetMaximum() * 1.1);
  line->SetLineWidth(3);
  line->SetLineColor(kRed);
  line->Draw();

  canvas->SaveAs("ptCutoff.gif");
}

void TriggerBias()
{
  TFile* file = TFile::Open("correction_map.root");

  TH2* corr = dynamic_cast<TH2*> (file->Get("dndeta_correction/corr_trigger"));

  Prepare2DPlot(corr);
  corr->SetTitle("Trigger bias correction");

  TCanvas* canvas = new TCanvas("TriggerBias", "TriggerBias", 500, 500);
  InitPadCOLZ();
  corr->DrawCopy("COLZ");

  canvas->SaveAs(Form("TriggerBias_%d.gif", gMax));

  corr->GetYaxis()->SetRangeUser(0, 5);

  canvas = new TCanvas("TriggerBiasZoom", "TriggerBiasZoom", 500, 500);
  InitPadCOLZ();
  corr->DrawCopy("COLZ");

  canvas->SaveAs(Form("TriggerBiasZoom_%d.gif", gMax));
}

void VtxRecon()
{
  TFile* file = TFile::Open("correction_map.root");

  TH2* corr = dynamic_cast<TH2*> (file->Get("dndeta_correction/corr_vtxReco"));

  Prepare2DPlot(corr);
  corr->SetTitle("Vertex reconstruction correction");

  TCanvas* canvas = new TCanvas("VtxRecon", "VtxRecon", 500, 500);
  InitPadCOLZ();
  corr->DrawCopy("COLZ");

  canvas->SaveAs(Form("VtxRecon_%d.gif", gMax));


  corr->GetYaxis()->SetRangeUser(0, 5);

  TCanvas* canvas = new TCanvas("VtxReconZoom", "VtxReconZoom", 500, 500);
  InitPadCOLZ();
  corr->DrawCopy("COLZ");

  canvas->SaveAs(Form("VtxReconZoom_%d.gif", gMax));
}

void Track2Particle2D()
{
  gSystem->Load("libPWG0base");

  TFile* file = TFile::Open("correction_map.root");
  AlidNdEtaCorrection* dNdEtaCorrection = new AlidNdEtaCorrection();
  dNdEtaCorrection->LoadHistograms("correction_map.root","dndeta_correction");

  TH3* gene = dNdEtaCorrection->GetTrack2ParticleCorrection()->GetGeneratedHistogram();
  TH3* meas = dNdEtaCorrection->GetTrack2ParticleCorrection()->GetMeasuredHistogram();

  gene->GetZaxis()->SetRangeUser(0.3, 10);
  meas->GetZaxis()->SetRangeUser(0.3, 10);
  AliPWG0Helper::CreateDividedProjections(gene, meas, "yx");
  gene->GetZaxis()->UnZoom();
  meas->GetZaxis()->UnZoom();

  gene->GetYaxis()->SetRangeUser(-0.8, 0.8);
  meas->GetYaxis()->SetRangeUser(-0.8, 0.8);
  AliPWG0Helper::CreateDividedProjections(gene, meas, "zx");
  gene->GetYaxis()->UnZoom();
  meas->GetYaxis()->UnZoom();

  gene->GetXaxis()->SetRangeUser(-10, 10);
  meas->GetXaxis()->SetRangeUser(-10, 10);
  AliPWG0Helper::CreateDividedProjections(gene, meas, "zy");
  gene->GetXaxis()->UnZoom();
  meas->GetXaxis()->UnZoom();

  TH2* corrYX = dynamic_cast<TH2*> (gROOT->FindObject("gene_nTrackToNPart_yx_div_meas_nTrackToNPart_yx"));
  TH2* corrZX = dynamic_cast<TH2*> (gROOT->FindObject("gene_nTrackToNPart_zx_div_meas_nTrackToNPart_zx"));
  TH2* corrZY = dynamic_cast<TH2*> (gROOT->FindObject("gene_nTrackToNPart_zy_div_meas_nTrackToNPart_zy"));

  /* this reads them from the file
  TH2* corrYX = dynamic_cast<TH2*> (file->Get("dndeta_correction/gene_nTrackToNPart_yx_div_meas_nTrackToNPart_yx"));
  TH2* corrZX = dynamic_cast<TH2*> (file->Get("dndeta_correction/gene_nTrackToNPart_zx_div_meas_nTrackToNPart_zx"));
  TH2* corrZY = dynamic_cast<TH2*> (file->Get("dndeta_correction/gene_nTrackToNPart_zy_div_meas_nTrackToNPart_zy"));*/

  Prepare2DPlot(corrYX);
  Prepare2DPlot(corrZX);
  Prepare2DPlot(corrZY);

  const char* title = "Track2Particle Correction";
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
}

void SetRanges(TH1* hist)
{
  SetRanges(hist->GetXaxis());
  SetRanges(hist->GetYaxis());
  SetRanges(hist->GetZaxis());
}

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
}

void InitPad()
{
  gPad->Range(0, 0, 1, 1);
  gPad->SetLeftMargin(0.15);
  //gPad->SetRightMargin(0.05);
  //gPad->SetTopMargin(0.13);
  //gPad->SetBottomMargin(0.1);

  //gPad->SetGridx();
  //gPad->SetGridy();
}

void InitPadCOLZ()
{
  gPad->Range(0, 0, 1, 1);
  gPad->SetRightMargin(0.15);
  gPad->SetLeftMargin(0.12);

  gPad->SetGridx();
  gPad->SetGridy();
}
