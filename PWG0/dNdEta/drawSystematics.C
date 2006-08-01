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

extern TPad* gPad;

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
  if (!gPad)
    return;

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
    TH1F* hist = new TH1F(Form("secondaries_%d", ptBin), Form("secondaries_%d", ptBin), secondaries->GetNbinsY(), secondaries->GetYaxis()->GetXmin(), secondaries->GetYaxis()->GetXmax());

    hist->SetTitle(Form("%f < p_{T} < %f", secondaries->GetZaxis()->GetBinLowEdge(ptBin), secondaries->GetZaxis()->GetBinUpEdge(ptBin)));

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
      {
        hist->SetBinContent(cBin, sum / count);
        hist->SetBinError(cBin, TMath::Sqrt(sum) / count);
      }
    }

    hist->Scale(1.0 / hist->GetBinContent(hist->GetXaxis()->FindBin(1)));
    hist->Add(new TF1("flat", "-1", 0, 2));

    new TCanvas;
    hist->SetMarkerStyle(21);
    hist->Draw("");
  }
}

void Track2Particle1DComposition(const char* fileName = "correction_map.root", Int_t folderCount, const char** folderNames, Float_t upperPtLimit = 9.9)
{
  gSystem->Load("libPWG0base");

  TString canvasName;
  canvasName.Form("Track2Particle1DComposition");
  TCanvas* canvas = new TCanvas(canvasName, canvasName, 1200, 400);
  canvas->Divide(3, 1);

  TLegend* legend = new TLegend(0.7, 0.7, 0.95, 0.95);

  for (Int_t i=0; i<folderCount; ++i)
  {
    Track2Particle1DCreatePlots(fileName, folderNames[i], upperPtLimit);

    TH1* corrX = dynamic_cast<TH1*> (gROOT->FindObject("gene_nTrackToNPart_x_div_meas_nTrackToNPart_x"));
    TH1* corrY = dynamic_cast<TH1*> (gROOT->FindObject("gene_nTrackToNPart_y_div_meas_nTrackToNPart_y"));
    TH1* corrZ = dynamic_cast<TH1*> (gROOT->FindObject("gene_nTrackToNPart_z_div_meas_nTrackToNPart_z"));

    Prepare1DPlot(corrX);
    Prepare1DPlot(corrY);
    Prepare1DPlot(corrZ);

    const char* title = "Track2Particle Correction";
    corrX->SetTitle(title);
    corrY->SetTitle(title);
    corrZ->SetTitle(title);

    corrZ->GetXaxis()->SetRangeUser(0, upperPtLimit);

    corrX->SetLineColor(i+1);
    corrY->SetLineColor(i+1);
    corrZ->SetLineColor(i+1);

    canvas->cd(1);
    InitPad();
    corrX->DrawCopy(((i>0) ? "SAME" : ""));

    canvas->cd(2);
    InitPad();
    corrY->DrawCopy(((i>0) ? "SAME" : ""));

    canvas->cd(3);
    InitPad();
    corrZ->DrawCopy(((i>0) ? "SAME" : ""));

    legend->AddEntry(corrZ, folderNames[i]);
  }

  legend->Draw();

  //canvas->SaveAs(Form("Track2Particle1D_%s_%d_%f.gif", fileName, gMax, upperPtLimit));
  //canvas->SaveAs(Form("Track2Particle1D_%s_%d_%f.eps", fileName, gMax, upperPtLimit));
}

void DrawDifferentSpecies()
{
  gROOT->ProcessLine(".L drawPlots.C");

  const char* folderNames[] = { "correction_0", "correction_1", "correction_2", "correction_3" };

  Track2Particle1DComposition("systematics.root", 4, folderNames);
}

void ScalePtDependent(TH3F* hist, TF1* function)
{
  // assumes that pt is the third dimension of hist
  // scales with function(pt)

  for (Int_t z=1; z<=hist->GetNbinsZ(); ++z)
  {
    Double_t factor = function->Eval(hist->GetZaxis()->GetBinCenter(z));
    printf("z = %d, pt = %f, scaling with %f\n", z, hist->GetZaxis()->GetBinCenter(z), factor);

    for (Int_t x=1; x<=hist->GetNbinsX(); ++x)
      for (Int_t y=1; y<=hist->GetNbinsY(); ++y)
        hist->SetBinContent(x, y, z, hist->GetBinContent(x, y, z) * factor);
  }
}

void ChangeComposition(void** correctionPointer, Int_t index)
{
  AlidNdEtaCorrection** fdNdEtaCorrection = (AlidNdEtaCorrection**) correctionPointer;

  switch (index)
  {
    case 0:
      break;

    case 1:
      fdNdEtaCorrection[0]->GetTrack2ParticleCorrection()->GetMeasuredHistogram()->Scale(10);
      fdNdEtaCorrection[0]->GetTrack2ParticleCorrection()->GetGeneratedHistogram()->Scale(10);
      break;

    case 2:
      TF1* ptDependence = new TF1("simple", "x", 0, 100);
      ScalePtDependent(fdNdEtaCorrection[0]->GetTrack2ParticleCorrection()->GetGeneratedHistogram(), ptDependence);
      ScalePtDependent(fdNdEtaCorrection[0]->GetTrack2ParticleCorrection()->GetMeasuredHistogram(), ptDependence);
      break;

  }
}

void Composition()
{
  gSystem->Load("libPWG0base");

  gSystem->Unlink("new_compositions.root");

  for (Int_t comp = 0; comp < 3; ++comp)
  {
    AlidNdEtaCorrection* fdNdEtaCorrection[4];

    TFile::Open("systematics.root");

    for (Int_t i=0; i<4; ++i)
    {
      TString name;
      name.Form("correction_%d", i);
      fdNdEtaCorrection[i] = new AlidNdEtaCorrection(name, name);
      fdNdEtaCorrection[i]->LoadHistograms("systematics.root", name);
    }

    ChangeComposition(fdNdEtaCorrection, comp);

    Double_t geneCount[5];
    Double_t measCount[5];
    geneCount[4] = 0;
    measCount[4] = 0;

    for (Int_t i=0; i<4; ++i)
    {
      geneCount[i] = fdNdEtaCorrection[i]->GetTrack2ParticleCorrection()->GetGeneratedHistogram()->Integral();
      measCount[i] = fdNdEtaCorrection[i]->GetTrack2ParticleCorrection()->GetMeasuredHistogram()->Integral();

      geneCount[4] += geneCount[i];
      measCount[4] += measCount[i];

      printf("Particle %s: %d gene, %d meas\n", ((i == 0) ? "pi" : (i == 1) ? "K" : (i == 2) ? "p" : "others"), (Int_t) geneCount[i], (Int_t) measCount[i]);
    }

    printf("Generated ratios are:     %f pi, %f K, %f p, %f others\n", geneCount[0] / geneCount[4], geneCount[1] / geneCount[4], geneCount[2] / geneCount[4], geneCount[3] / geneCount[4]);

    printf("Reconstructed ratios are: %f pi, %f K, %f p, %f others\n", measCount[0] / measCount[4], measCount[1] / measCount[4], measCount[2] / measCount[4], measCount[3] / measCount[4]);

    TList* collection = new TList;

    for (Int_t i=0; i<4; ++i)
      collection->Add(fdNdEtaCorrection[i]);

    TString correctionName;
    correctionName.Form("new_composition_%d", comp);

    AlidNdEtaCorrection* newComposition = new AlidNdEtaCorrection(correctionName, correctionName);
    newComposition->Merge(collection);
    newComposition->Finish();

    delete collection;

    TFile* file = TFile::Open("new_compositions.root", "UPDATE");
    newComposition->SaveHistograms();
    //file->Write();
    file->Close();
  }

  gROOT->ProcessLine(".L drawPlots.C");

  const char* folderNames[] = { "new_composition_0", "new_composition_1", "new_composition_2" };

  Track2Particle1DComposition("new_compositions.root", 3, folderNames);
}

Double_t ConvSigma1To2D(Double_t sigma)
{
  return TMath::Sqrt( - TMath::Log( 1 - TMath::Erf(sigma / TMath::Sqrt(2)) ) * 2);
}

Double_t ConvDistance1DTo2D(Double_t distance)
{
  return TMath::ErfInverse(1 - TMath::Exp(-distance * distance / 2)) * TMath::Sqrt(2);
}

Double_t Sigma2VertexCount(TH2F* tracks, Double_t nSigma)
{
  Double_t count = 0;

  //nSigma = ConvSigma1To2D(nSigma);

  for (Int_t x=1; x<=tracks->GetNbinsX(); ++x)
    for (Int_t y=1; y<=tracks->GetNbinsY(); ++y)
    {
      Double_t impactX = tracks->GetXaxis()->GetBinCenter(x);
      Double_t impactY = tracks->GetYaxis()->GetBinCenter(y);

      Float_t d = TMath::Sqrt(impactX*impactX + impactY*impactY);

      d = ConvDistance1DTo2D(d);

      if (d < nSigma)
        count += tracks->GetBinContent(x, y);
    }

  return count;
}

TH2F* Sigma2VertexGaussianTracksHist()
{
  TH2F* tracks = new TH2F("Sigma2Vertex_tracks", "Sigma2Vertex_tracks", 200, -5, 5, 200, -5, 5);

  TF2* gaussian2D = new TF2("gaussian2D", "xgausn(0) * ygausn(3)", -5, 5, -5, 5);
  gaussian2D->SetParameters(1, 0, 1, 1, 0, 1);

  for (Int_t x=1; x<=tracks->GetNbinsX(); ++x)
    for (Int_t y=1; y<=tracks->GetNbinsY(); ++y)
      tracks->SetBinContent(x, y, gaussian2D->Eval(tracks->GetXaxis()->GetBinCenter(x), tracks->GetYaxis()->GetBinCenter(y)));

  //normalize
  tracks->Scale(1.0 / tracks->Integral());

  return tracks;
}

TH1F* Sigma2VertexGaussian()
{
  TH2F* tracks = Sigma2VertexGaussianTracksHist();

  TCanvas* canvas = new TCanvas("Sigma2VertexGaussian", "Sigma2VertexGaussian", 1000, 1000);
  canvas->Divide(2, 2);

  canvas->cd(1);
  tracks->Draw("COLZ");

  TH1F* ratio = new TH1F("Sigma2Vertex_ratio", "Sigma2Vertex_ratio;n sigma;included", 10, 0.25, 5.25);
  for (Double_t nSigma = 0.5; nSigma < 5.1; nSigma += 0.5)
    ratio->Fill(nSigma, Sigma2VertexCount(tracks, nSigma));
  ratio->SetMarkerStyle(21);

  canvas->cd(2);
  ratio->DrawCopy("P");

  TH1F* ratio2 = new TH1F("Sigma2Vertex_ratio2", "Sigma2Vertex_ratio2;nSigma;% included 3 sigma / % included n sigma", 10, 0.25, 5.25);
  Double_t sigma3 = Sigma2VertexCount(tracks, 3);
  for (Double_t nSigma = 0.5; nSigma < 5.1; nSigma += 0.5)
    ratio2->Fill(nSigma, sigma3 / Sigma2VertexCount(tracks, nSigma));
  ratio2->SetMarkerStyle(21);

  canvas->cd(3);
  ratio2->DrawCopy("P");

  canvas->SaveAs("Sigma2Vertex.eps");

  return ratio2;
}

TH1F* Sigma2VertexSimulation()
{
  TFile* file = TFile::Open("systematics.root");

  TH1F* sigmavertex = dynamic_cast<TH1F*> (file->Get("fSigmaVertex"));
  if (!sigmavertex)
  {
    printf("Could not read histogram\n");
    return;
  }

  TH1F* ratio = new TH1F("sigmavertexsimulation_ratio", "sigmavertexsimulation_ratio;Nsigma;% included 3 sigma / % included n sigma", sigmavertex->GetNbinsX(), sigmavertex->GetXaxis()->GetXmin(), sigmavertex->GetXaxis()->GetXmax());

  for (Int_t i=1; i<=sigmavertex->GetNbinsX(); ++i)
    ratio->SetBinContent(i, sigmavertex->GetBinContent(sigmavertex->GetXaxis()->FindBin(3)) / sigmavertex->GetBinContent(i));

  TCanvas* canvas = new TCanvas("Sigma2VertexSimulation", "Sigma2VertexSimulation", 1000, 500);
  canvas->Divide(2, 1);

  canvas->cd(1);
  sigmavertex->SetMarkerStyle(21);
  sigmavertex->Draw("P");

  canvas->cd(2);
  ratio->SetMarkerStyle(21);
  ratio->DrawCopy("P");

  return ratio;
}

void Sigma2VertexCompare()
{
  TH1F* ratio1 = Sigma2VertexGaussian();
  TH1F* ratio2 = Sigma2VertexSimulation();

  ratio1->SetStats(kFALSE);
  ratio2->SetStats(kFALSE);

  ratio1->SetMarkerStyle(0);
  ratio2->SetMarkerStyle(0);

  TLegend* legend = new TLegend(0.647177,0.775424,0.961694,0.966102);
  legend->AddEntry(ratio1, "Gaussian");
  legend->AddEntry(ratio2, "Simulation");

  ratio1->GetXaxis()->SetTitleOffset(1.5);

  TCanvas* canvas = new TCanvas("Sigma2VertexCompare", "Sigma2VertexCompare", 500, 500);
  InitPad();

  ratio1->Draw();
  ratio2->SetLineColor(kRed);
  ratio2->Draw("SAME");

  legend->Draw();
}

void drawSystematics()
{
  //Secondaries();
  //DrawDifferentSpecies();
  //Composition();

  Sigma2VertexSimulation();
}
