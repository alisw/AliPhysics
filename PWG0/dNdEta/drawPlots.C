/* $Id$ */

void drawPlots()
//void Track2Particle2D()
{
  TFile* file = TFile::Open("correction_map.root");

  TH2* corrYX = dynamic_cast<TH2*> (file->Get("dndeta_correction/corr_nTrackToNPart_yx"));
  TH2* corrZX = dynamic_cast<TH2*> (file->Get("dndeta_correction/corr_nTrackToNPart_zx"));
  TH2* corrZY = dynamic_cast<TH2*> (file->Get("dndeta_correction/corr_nTrackToNPart_zy"));

  Prepare2DPlot(corrYX);
  Prepare2DPlot(corrZX);
  Prepare2DPlot(corrZY);

  SetRanges(corrYX);
  SetRanges(corrZX);
  SetRanges(corrZY);

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
}

void Track2Particle3D()
{
  // get left margin proper

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

  TCanvas* canvas = new TCanvas("Track2Particle3D", "Track2Particle3D", 1200, 400);
  canvas->Divide(3, 1);

  canvas->cd(1);
  InitPad();
  gene->Draw();

  canvas->cd(2);
  meas->Draw();

  canvas->cd(3);
  corr->Draw();
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
    axis->SetRangeUser(-0.8, 0.79999);
  if (strcmp(axis->GetTitle(), "p_{T}") == 0)
    axis->SetRangeUser(0, 9.9999);
  if (strcmp(axis->GetTitle(), "vtx z [cm]") == 0)
    axis->SetRangeUser(-10, 9.9999);
}

void Prepare3DPlot(TH3* hist)
{
  hist->GetXaxis()->SetTitleOffset(1.5);
  hist->GetYaxis()->SetTitleOffset(1.5);
  hist->GetZaxis()->SetTitleOffset(1.5);
}

void Prepare2DPlot(TH2* hist)
{
  hist->SetStats(kFALSE);
}

void InitPad()
{
  gPad->Range(0, 0, 1, 1);
  gPad->SetLeftMargin(0);
  gPad->SetRightMargin(0.05);
  gPad->SetTopMargin(0.13);
  gPad->SetBottomMargin(0.1);

  gPad->SetGridx();
  gPad->SetGridy();
}

void InitPadCOLZ()
{  
  gPad->Range(0, 0, 1, 1);
  gPad->SetRightMargin(0.15);
  gPad->SetLeftMargin(0.12);
}
