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

  SetRanges(hist);
}

void Prepare1DPlot(TH1* hist)
{
  hist->SetLineWidth(2);
  hist->SetStats(kFALSE);

  SetRanges(hist);
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

void Secondaries()
{
  TFile* file = TFile::Open("systematics.root");

  TH3F* secondaries = dynamic_cast<TH3F*> (file->Get("fSecondaries"));
  if (!secondaries)
  {
    printf("Could not read histogram\n");
    return;
  }

  for (Int_t ptBin=1; ptBin<=secondaries->GetNbinsZ(); ptBin++)
  //for (Int_t ptBin = 1; ptBin<=2; ptBin++)
  {
    TGraph* graph = new TGraph;
    graph->Clear();
    graph->SetTitle(Form("%f < p_{T} < %f", secondaries->GetZaxis()->GetBinLowEdge(ptBin), secondaries->GetZaxis()->GetBinUpEdge(ptBin)));

    for (Int_t cBin=1; cBin<=secondaries->GetNbinsY(); ++cBin)
    {
      if (secondaries->GetBinContent(0, cBin, ptBin) > 0)
        printf("WARNING: Underflow bin not empty!");
      if (secondaries->GetBinContent(secondaries->GetNbinsX()+1, cBin, ptBin) > 0)
        printf("WARNING: Overflow bin not empty!");

      Double_t sum = 0;
      Double_t count = 0;
      for (Int_t nBin=1; nBin<=secondaries->GetNbinsX(); ++nBin)
      {
        //printf("%f %f\n", secondaries->GetXaxis()->GetBinCenter(nBin), secondaries->GetBinContent(nBin, cBin, ptBin));
        sum += secondaries->GetXaxis()->GetBinCenter(nBin) * secondaries->GetBinContent(nBin, cBin, ptBin);
        count += secondaries->GetBinContent(nBin, cBin, ptBin);
      }

      printf("%f %f\n", sum, count);

      if (count > 0)
        graph->SetPoint(graph->GetN(), secondaries->GetYaxis()->GetBinCenter(cBin), sum / count);
    }

    new TCanvas;
    graph->SetMarkerStyle(21);
    graph->Draw("AP");
    graph->Print();
  }
}

void Composition()
{
  gSystem->Load("libPWG0base");

  AlidNdEtaCorrection* fdNdEtaCorrection[4];

  TFile::Open("systematics.root");

  for (Int_t i=0; i<4; ++i)
  {
    TString name;
    name.Form("correction_%d", i);
    fdNdEtaCorrection[i] = new AlidNdEtaCorrection(name, name);
    fdNdEtaCorrection[i]->LoadHistograms("systematics.root", name);
  }

  //fdNdEtaCorrection[0]->GetTrack2ParticleCorrection()->GetMeasuredHistogram()->Scale(2);
  //fdNdEtaCorrection[0]->GetTrack2ParticleCorrection()->GetGeneratedHistogram()->Scale(2);

  AlidNdEtaCorrection* finalCorrection = new AlidNdEtaCorrection("dndeta_correction", "dndeta_correction");

  TList* collection = new TList;

  for (Int_t i=0; i<4; ++i)
    collection->Add(fdNdEtaCorrection[i]);

  finalCorrection->Merge(collection);

  delete collection;

  finalCorrection->Finish();

  TFile* file = TFile::Open("temp.root", "RECREATE");
  finalCorrection->SaveHistograms();
  file->Write();
  file->Close();

  gROOT->ProcessLine(".L drawPlots.C");
  Track2Particle1D("temp.root");
}

void drawSystematics()
{
  Composition();
}
