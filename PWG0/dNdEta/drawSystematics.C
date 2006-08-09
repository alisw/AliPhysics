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
#include <TLegend.h>
#include <TPad.h>
#include <TF1.h>

extern TPad* gPad;

void Track2Particle1DCreatePlots(const char* fileName = "correction_map.root", const char* folderName = "dndeta_correction", Float_t upperPtLimit = 10);

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

void Track2Particle1DComposition(const char** fileNames, Int_t folderCount, const char** folderNames, Float_t upperPtLimit = 9.9)
{
  gSystem->Load("libPWG0base");

  TString canvasName;
  canvasName.Form("Track2Particle1DComposition");
  TCanvas* canvas = new TCanvas(canvasName, canvasName, 1200, 400);
  canvas->Divide(3, 1);

  TLegend* legend = new TLegend(0.7, 0.7, 0.95, 0.95);

  for (Int_t i=0; i<folderCount; ++i)
  {
    Track2Particle1DCreatePlots(fileNames[i], folderNames[i], upperPtLimit);

    TH1* corrX = dynamic_cast<TH1*> (gROOT->FindObject(Form("gene_%s_nTrackToNPart_x_div_meas_%s_nTrackToNPart_x", folderNames[i], folderNames[i])));
    TH1* corrY = dynamic_cast<TH1*> (gROOT->FindObject(Form("gene_%s_nTrackToNPart_y_div_meas_%s_nTrackToNPart_y", folderNames[i], folderNames[i])));
    TH1* corrZ = dynamic_cast<TH1*> (gROOT->FindObject(Form("gene_%s_nTrackToNPart_z_div_meas_%s_nTrackToNPart_z", folderNames[i], folderNames[i])));

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

TH1** DrawRatios(const char* fileName = "systematics.root")
{
  gSystem->Load("libPWG0base");

  TFile* file = TFile::Open(fileName);

  TString canvasName;
  canvasName.Form("DrawRatios");
  TCanvas* canvas = new TCanvas(canvasName, canvasName, 800, 400);
  canvas->Divide(2, 1);
  canvas->cd(1);

  TH1** ptDists = new TH1*[4];

  TLegend* legend = new TLegend(0.7, 0.7, 0.95, 0.95);

  const char* folderNames[] = { "correction_0", "correction_1", "correction_2", "correction_3" };
  const char* particleNames[] = { "#pi", "K", "p", "other" };
  for (Int_t i=0; i<4; ++i)
  {
    AlidNdEtaCorrection* dNdEtaCorrection = new AlidNdEtaCorrection(folderNames[i], folderNames[i]);
    dNdEtaCorrection->LoadHistograms(fileName, folderNames[i]);

    TH3F* gene = dNdEtaCorrection->GetTrack2ParticleCorrection()->GetGeneratedHistogram();

    gene->GetYaxis()->SetRangeUser(-0.8, 0.8);
    gene->GetXaxis()->SetRangeUser(-10, 10);

    ptDists[i] = dynamic_cast<TH1*> (gene->Project3D("z")->Clone(Form("%s_z", folderNames[i])));
    ptDists[i]->SetTitle(";p_{T};Count");
    if (!ptDists[i])
    {
      printf("Problem getting distribution %d.\n", i);
      return 0;
    }

    ptDists[i]->GetYaxis()->SetRangeUser(1, ptDists[i]->GetMaximum()*1.1);
    ptDists[i]->GetXaxis()->SetRangeUser(0, 9.9);
    ptDists[i]->SetLineColor(i+1);
    ptDists[i]->DrawCopy((i == 0) ? "" : "SAME");
    ptDists[i]->GetYaxis()->SetRange(0, 0);

    legend->AddEntry(ptDists[i], particleNames[i]);
  }
  gPad->SetLogy();

  TH1* total = dynamic_cast<TH1*> (ptDists[0]->Clone("total"));

  for (Int_t i=1; i<4; ++i)
    total->Add(ptDists[i]);

  canvas->cd(2);
  for (Int_t i=0; i<4; ++i)
  {
    ptDists[i]->Divide(total);
    ptDists[i]->SetTitle(";p_{T};Ratio");
    ptDists[i]->GetYaxis()->SetRangeUser(0, 1);
    ptDists[i]->Draw((i == 0) ? "" : "SAME");
  }
  legend->Draw();

  canvas->SaveAs("DrawRatios.gif");

  file->Close();

  return ptDists;
}

void DrawDifferentSpecies()
{
  gROOT->ProcessLine(".L drawPlots.C");

  const char* fileNames[] = { "systematics.root", "systematics.root", "systematics.root", "systematics.root" };
  const char* folderNames[] = { "correction_0", "correction_1", "correction_2", "correction_3" };

  Track2Particle1DComposition(fileNames, 4, folderNames);
}

void DrawSpeciesAndCombination()
{
  gROOT->ProcessLine(".L drawPlots.C");

  const char* fileNames[] = { "systematics.root", "systematics.root", "systematics.root", "new_compositions.root" };
  const char* folderNames[] = { "correction_0", "correction_1", "correction_2", "PythiaRatios" };

  Track2Particle1DComposition(fileNames, 4, folderNames);
}

void DrawSimulatedVsCombined()
{
  gROOT->ProcessLine(".L drawPlots.C");

  const char* fileNames[] = { "new_compositions.root", "new_compositions.root" };
  const char* folderNames[] = { "Pythia", "PythiaRatios" };

  Track2Particle1DComposition(fileNames, 2, folderNames);
}

void DrawBoosts()
{
  gROOT->ProcessLine(".L drawPlots.C");

  const char* fileNames[] = { "new_compositions.root", "new_compositions.root", "new_compositions.root", "new_compositions.root" };
  const char* folderNames[] = { "PythiaRatios", "PiBoosted", "KBoosted", "pBoosted" };

  Track2Particle1DComposition(fileNames, 4, folderNames);
}

TH1F* HistogramDifferences(const char* filename, const char* folder1, const char* folder2)
{
  gSystem->Load("libPWG0base");

  AlidNdEtaCorrection* fdNdEtaCorrection[2];

  TFile::Open(filename);

  fdNdEtaCorrection[0] = new AlidNdEtaCorrection(folder1, folder1);
  fdNdEtaCorrection[0]->LoadHistograms(filename, folder1);

  fdNdEtaCorrection[1] = new AlidNdEtaCorrection(folder2, folder2);
  fdNdEtaCorrection[1]->LoadHistograms(filename, folder2);

  TH1F* difference = new TH1F("difference", Form(";#DeltaC_{pT, z, #eta} %s - %s;Count", folder2, folder1), 1000, -0.5, 0.5);

  TH3F* hist1 = fdNdEtaCorrection[0]->GetTrack2ParticleCorrection()->GetCorrectionHistogram();
  TH3F* hist2 = fdNdEtaCorrection[1]->GetTrack2ParticleCorrection()->GetCorrectionHistogram();

  for (Int_t x=hist1->GetXaxis()->FindBin(-10); x<=hist1->GetXaxis()->FindBin(10); ++x)
    for (Int_t y=hist1->GetYaxis()->FindBin(-0.8); y<=hist1->GetYaxis()->FindBin(0.8); ++y)
      for (Int_t z=hist1->GetZaxis()->FindBin(0.3); z<=hist1->GetZaxis()->FindBin(9.9); ++z)
        difference->Fill(hist2->GetBinContent(x, y, z) - hist1->GetBinContent(x, y, z));

  printf("Over-/Underflow bins: %d %d\n", difference->GetBinContent(0), difference->GetBinContent(difference->GetNbinsX()+1));

  return difference;
}

void HistogramDifferences()
{
  TH1F* PiBoosted = HistogramDifferences("new_compositions.root", "PythiaRatios", "PiBoosted");
  TH1F* KBoosted = HistogramDifferences("new_compositions.root", "PythiaRatios", "KBoosted");
  TH1F* pBoosted = HistogramDifferences("new_compositions.root", "PythiaRatios", "pBoosted");

  TCanvas* canvas = new TCanvas("HistogramDifferences", "HistogramDifferences", 1200, 400);
  canvas->Divide(3, 1);

  canvas->cd(1);
  PiBoosted->GetXaxis()->SetRangeUser(-0.01, 0.01);
  PiBoosted->Draw();

  canvas->cd(2);
  KBoosted->GetXaxis()->SetRangeUser(-0.01, 0.01);
  KBoosted->Draw();

  canvas->cd(3);
  pBoosted->GetXaxis()->SetRangeUser(-0.01, 0.01);
  pBoosted->Draw();

  canvas->SaveAs("HistogramDifferences.gif");
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

void ScalePtDependent(TH3F* hist, TH1* function)
{
  // assumes that pt is the third dimension of hist
  // scales with histogram(pt)

  for (Int_t z=1; z<=hist->GetNbinsZ(); ++z)
  {
    Double_t factor = function->GetBinContent(function->GetXaxis()->FindBin(hist->GetZaxis()->GetBinCenter(z)));
    printf("z = %d, pt = %f, scaling with %f\n", z, hist->GetZaxis()->GetBinCenter(z), factor);

    for (Int_t x=1; x<=hist->GetNbinsX(); ++x)
      for (Int_t y=1; y<=hist->GetNbinsY(); ++y)
        hist->SetBinContent(x, y, z, hist->GetBinContent(x, y, z) * factor);
  }
}

const char* ChangeComposition(void** correctionPointer, Int_t index)
{
  AlidNdEtaCorrection** fdNdEtaCorrection = (AlidNdEtaCorrection**) correctionPointer;

  switch (index)
  {
    case 0: // result from pp events
      {
        TFile::Open("pythiaratios.root");

        for (Int_t i=0; i<4; ++i)
        {
          TString name;
          name.Form("correction_%d", i);
          fdNdEtaCorrection[i] = new AlidNdEtaCorrection(name, name);
          fdNdEtaCorrection[i]->LoadHistograms("pythiaratios.root", name);
        }
      }
      return "Pythia";
      break;

    case 1: // each species rated with pythia ratios
      TH1** ptDists = DrawRatios("pythiaratios.root");

      for (Int_t i=0; i<3; ++i)
      {
        ScalePtDependent(fdNdEtaCorrection[i]->GetTrack2ParticleCorrection()->GetMeasuredHistogram(), ptDists[i]);
        ScalePtDependent(fdNdEtaCorrection[i]->GetTrack2ParticleCorrection()->GetGeneratedHistogram(), ptDists[i]);
      }
      return "PythiaRatios";
      break;

      // each species rated with pythia ratios
    case 2: // + 10% pions
    case 3: // - 10% pions
    case 4: // + 10% kaons
    case 5: // - 10% kaons
    case 6: // + 10% protons
    case 7: // - 10% protons
      TH1** ptDists = DrawRatios("pythiaratios.root");
      Int_t functionIndex = (index - 2) / 2;
      Double_t scaleFactor = (index % 2 == 0) ? 1.1 : 0.9;
      ptDists[functionIndex]->Scale(scaleFactor);

      for (Int_t i=0; i<3; ++i)
      {
        ScalePtDependent(fdNdEtaCorrection[i]->GetTrack2ParticleCorrection()->GetMeasuredHistogram(), ptDists[i]);
        ScalePtDependent(fdNdEtaCorrection[i]->GetTrack2ParticleCorrection()->GetGeneratedHistogram(), ptDists[i]);
      }
      TString* str = new TString;
      str->Form("%s%s", (functionIndex == 0) ? "Pi" : ((functionIndex == 1) ? "K" : "p"), (index % 2 == 0) ? "Boosted" : "Reduced");
      return str->Data();
      break;

    case 999:
      TF1* ptDependence = new TF1("simple", "x", 0, 100);
      ScalePtDependent(fdNdEtaCorrection[0]->GetTrack2ParticleCorrection()->GetGeneratedHistogram(), ptDependence);
      ScalePtDependent(fdNdEtaCorrection[0]->GetTrack2ParticleCorrection()->GetMeasuredHistogram(), ptDependence);
      break;

  }

  return "noname";
}

void Composition()
{
  gSystem->Load("libPWG0base");

  gSystem->Unlink("new_compositions.root");

  Int_t nCompositions = 8;
  for (Int_t comp = 0; comp < nCompositions; ++comp)
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

    const char* newName = ChangeComposition(fdNdEtaCorrection, comp);

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

    AlidNdEtaCorrection* newComposition = new AlidNdEtaCorrection(newName, newName);
    newComposition->Merge(collection);
    newComposition->Finish();

    delete collection;

    TFile* file = TFile::Open("new_compositions.root", "UPDATE");
    newComposition->SaveHistograms();
    //file->Write();
    file->Close();
  }

  gROOT->ProcessLine(".L drawPlots.C");

  const char* fileNames[] = {"new_compositions.root", "new_compositions.root", "new_compositions.root", "new_compositions.root", "new_compositions.root", "new_compositions.root", "new_compositions.root", "new_compositions.root" };
  const char* folderNames[] = { "Pythia", "PythiaRatios", "PiBoosted", "PiReduced", "KBoosted", "KReduced", "pBoosted", "pReduced" };

  Track2Particle1DComposition(fileNames, nCompositions, folderNames);
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
