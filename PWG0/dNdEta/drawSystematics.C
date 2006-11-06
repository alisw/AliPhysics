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

Int_t markers[] = {20,20,21,22,23,28,29};
Int_t colors[]  = {1,2,3,4,6,8,102};

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

void Prepare1DPlot(TH1* hist, Bool_t setRanges = kTRUE)
{
  hist->SetLineWidth(2);
  hist->SetStats(kFALSE);

  hist->GetXaxis()->SetTitleOffset(1.2);
  hist->GetYaxis()->SetTitleOffset(1.2);

  if (setRanges)
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

void Secondaries()
{
  TFile* file = TFile::Open("systematics.root");

  TH2F* secondaries = dynamic_cast<TH2F*> (file->Get("fSecondaries"));
  if (!secondaries)
  {
    printf("Could not read histogram\n");
    return;
  }

  TCanvas* canvas = new TCanvas("Secondaries", "Secondaries", 1000, 1000);
  canvas->Divide(3, 3);
  for (Int_t i=1; i<=8; i++)
  {
    TH1D* hist = secondaries->ProjectionY(Form("proj_%d", i), i, i);
    hist->SetTitle(secondaries->GetXaxis()->GetBinLabel(i));

    canvas->cd(i);
    hist->Draw();
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

  TLegend* legend = new TLegend(0.73, 0.73, 0.98, 0.98);

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
    ptDists[i]->SetStats(kFALSE);
    ptDists[i]->SetTitle(";p_{T};Fraction of total");
    ptDists[i]->GetYaxis()->SetRangeUser(0, 1);
    ptDists[i]->Draw((i == 0) ? "" : "SAME");
  }
  legend->SetFillColor(0);
  legend->Draw();

  canvas->SaveAs("DrawRatios.gif");


  canvas = new TCanvas("PythiaRatios", "PythiaRatios", 500, 500);
  for (Int_t i=0; i<4; ++i)
  {
    TH1* hist = ptDists[i]->Clone();
    hist->GetXaxis()->SetRangeUser(0, 1.9);
    hist->Draw((i == 0) ? "" : "SAME");
  }
  legend->Draw();

  canvas->SaveAs("pythiaratios.eps");

  file->Close();

  return ptDists;
}

void DrawpiKpAndCombinedZOnly(Float_t upperPtLimit=0.99)
{
  gROOT->ProcessLine(".L drawPlots.C");
  gSystem->Load("libPWG0base");

  const char* fileNames[] = { "systematics.root", "systematics.root", "systematics.root", "correction_map.root" };
  const char* folderNames[] = { "correction_0", "correction_1", "correction_2", "dndeta_correction" };
  const char* legendNames[] = { "#pi", "K", "p", "standard" };
  Int_t folderCount = 3;

  /*const char* fileNames[] = { "systematics.root", "systematics.root", "systematics.root", "systematics.root", "correction_map.root" };
  const char* folderNames[] = { "correction_0", "correction_1", "correction_2", "correction_3", "dndeta_correction" };
  const char* legendNames[] = { "#pi", "K", "p", "others", "standard" };
  Int_t folderCount = 5;*/

  TString canvasName;
  canvasName.Form("Track2Particle1DComposition");
  TCanvas* canvas = new TCanvas(canvasName, canvasName, 700, 500);
  canvas->SetGridx();
  canvas->SetGridy();
  canvas->SetBottomMargin(0.12);
  //InitPad();

  TLegend* legend = new TLegend(0.8, 0.7, 0.95, 0.95);
  legend->SetFillColor(0);

  Int_t mycolors[] = {1, 2, 4};

  for (Int_t i=0; i<folderCount; ++i)
  {
    Track2Particle1DCreatePlots(fileNames[i], folderNames[i], upperPtLimit);

    TH1* corrZ = dynamic_cast<TH1*> (gROOT->FindObject(Form("gene_%s_nTrackToNPart_z_div_meas_%s_nTrackToNPart_z", folderNames[i], folderNames[i])));

    Prepare1DPlot(corrZ);

    corrZ->SetTitle("");
    corrZ->GetXaxis()->SetRangeUser(0, upperPtLimit);
    corrZ->GetYaxis()->SetRangeUser(0.51, 6);
    corrZ->SetMarkerColor(mycolors[i]);
    corrZ->SetLineColor(mycolors[i]);
    corrZ->SetMarkerStyle(markers[i+1]);
    corrZ->GetYaxis()->SetTitle("correction factor");

    corrZ->DrawCopy(((i>0) ? "SAMEP" : "P"));

    legend->AddEntry(corrZ, legendNames[i]);
  }

  legend->Draw();

  canvas->SaveAs("ptcutoff_species.eps");
}

void DrawCompareToReal()
{
  gROOT->ProcessLine(".L drawPlots.C");

  const char* fileNames[] = { "correction_map.root", "new_compositions.root" };
  const char* folderNames[] = { "dndeta_correction", "PythiaRatios" };

  Track2Particle1DComposition(fileNames, 2, folderNames);
}

void DrawDifferentSpecies()
{
  gROOT->ProcessLine(".L drawPlots.C");

  const char* fileNames[] = { "systematics.root", "systematics.root", "systematics.root", "systematics.root" };
  const char* folderNames[] = { "correction_0", "correction_1", "correction_2", "correction_3" };

  Track2Particle1DComposition(fileNames, 4, folderNames);
}

void DrawpiKpAndCombined()
{
  gROOT->ProcessLine(".L drawPlots.C");

  const char* fileNames[] = { "systematics.root", "systematics.root", "systematics.root", "correction_map.root" };
  const char* folderNames[] = { "correction_0", "correction_1", "correction_2", "dndeta_correction" };

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

TH2F* HistogramDifferences(const char* filename, const char* folder1, const char* folder2)
{
  gSystem->Load("libPWG0base");

  AlidNdEtaCorrection* fdNdEtaCorrection[2];

  TFile::Open(filename);

  fdNdEtaCorrection[0] = new AlidNdEtaCorrection(folder1, folder1);
  fdNdEtaCorrection[0]->LoadHistograms(filename, folder1);

  fdNdEtaCorrection[1] = new AlidNdEtaCorrection(folder2, folder2);
  fdNdEtaCorrection[1]->LoadHistograms(filename, folder2);

  TH3F* hist1 = fdNdEtaCorrection[0]->GetTrack2ParticleCorrection()->GetCorrectionHistogram();
  TH3F* hist2 = fdNdEtaCorrection[1]->GetTrack2ParticleCorrection()->GetCorrectionHistogram();

  //TH1F* difference = new TH1F("difference", Form(";#DeltaC_{pT, z, #eta} %s / %s;Count", folder2, folder1), 1000, 0.9, 1.1);
  TH2F* difference = new TH2F(Form("difference_%s_%s", folder1, folder2), Form(";#Sigma (C_{pT, z} %s / C_{pT, z} %s);#eta;Count", folder2, folder1), 100, 0.9, 1.1, hist1->GetYaxis()->GetNbins(), hist1->GetYaxis()->GetXmin(), hist1->GetYaxis()->GetXmax());

  for (Int_t x=hist1->GetXaxis()->FindBin(-10); x<=hist1->GetXaxis()->FindBin(10); ++x)
    for (Int_t y=hist1->GetYaxis()->FindBin(-0.8); y<=hist1->GetYaxis()->FindBin(0.8); ++y)
      for (Int_t z=hist1->GetZaxis()->FindBin(0.3); z<=hist1->GetZaxis()->FindBin(9.9); ++z)
        if (hist1->GetBinContent(x, y, z) != 0)
          difference->Fill(hist2->GetBinContent(x, y, z) / hist1->GetBinContent(x, y, z), hist1->GetYaxis()->GetBinCenter(y));

  difference->GetYaxis()->SetRangeUser(-0.8, 0.8);

  printf("Over-/Underflow bins: %d %d\n", difference->GetBinContent(0), difference->GetBinContent(difference->GetNbinsX()+1));

  return difference;
}

void HistogramDifferences()
{
  TH2F* KBoosted = HistogramDifferences("new_compositions.root", "PythiaRatios", "KBoosted");
  TH2F* pBoosted = HistogramDifferences("new_compositions.root", "PythiaRatios", "pBoosted");
  TH2F* KReduced = HistogramDifferences("new_compositions.root", "PythiaRatios", "KReduced");
  TH2F* pReduced = HistogramDifferences("new_compositions.root", "PythiaRatios", "pReduced");

  TCanvas* canvas = new TCanvas("HistogramDifferences", "HistogramDifferences", 1000, 1000);
  canvas->Divide(2, 2);

  canvas->cd(1);
  KBoosted->GetXaxis()->SetRangeUser(-0.05, 0.05);
  KBoosted->Draw("COLZ");

  canvas->cd(2);
  KReduced->GetXaxis()->SetRangeUser(-0.05, 0.05);
  KReduced->Draw("COLZ");

  canvas->cd(3);
  pBoosted->GetXaxis()->SetRangeUser(-0.02, 0.02);
  pBoosted->Draw("COLZ");

  canvas->cd(4);
  pReduced->GetXaxis()->SetRangeUser(-0.02, 0.02);
  pReduced->Draw("COLZ");

  canvas->SaveAs("HistogramDifferences.gif");

  canvas = new TCanvas("HistogramDifferencesProfile", "HistogramDifferencesProfile", 1000, 500);
  canvas->Divide(2, 1);

  canvas->cd(1);
  gPad->SetBottomMargin(0.13);
  KBoosted->SetTitle("a) Ratio of correction maps");
  KBoosted->SetStats(kFALSE);
  KBoosted->GetXaxis()->SetTitleOffset(1.4);
  KBoosted->GetXaxis()->SetLabelOffset(0.02);
  KBoosted->Draw("COLZ");

  canvas->cd(2);
  gPad->SetGridx();
  gPad->SetGridy();

  TLegend* legend = new TLegend(0.73, 0.73, 0.98, 0.98);

  for (Int_t i=0; i<4; ++i)
  {
    TH2F* hist = 0;
    TString name;
    switch (i)
    {
      case 0: hist = KBoosted; name = "K enhanced"; break;
      case 1: hist = KReduced; name = "K reduced"; break;
      case 2: hist = pBoosted; name = "p enhanced"; break;
      case 3: hist = pReduced; name = "p reduced"; break;
    }

    TProfile* profile = hist->ProfileY();
    profile->SetTitle("b) Mean and RMS");
    profile->GetXaxis()->SetRange(hist->GetYaxis()->GetFirst(), hist->GetYaxis()->GetLast());
    profile->GetXaxis()->SetTitleOffset(1.2);
    profile->GetXaxis()->SetLabelOffset(0.02);
    profile->GetYaxis()->SetRangeUser(0.98, 1.02);
    profile->SetStats(kFALSE);
    profile->SetLineColor(i+1);
    profile->SetMarkerColor(i+1);
    profile->DrawCopy(((i > 0) ? "SAME" : ""));


    legend->AddEntry(profile, name);
  }

  legend->Draw();
  canvas->SaveAs("particlecomposition_result.eps");
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
      /*TH1** ptDists = DrawRatios("pythiaratios.root");

      for (Int_t i=0; i<3; ++i)
      {
        ScalePtDependent(fdNdEtaCorrection[i]->GetTrack2ParticleCorrection()->GetMeasuredHistogram(), ptDists[i]);
        ScalePtDependent(fdNdEtaCorrection[i]->GetTrack2ParticleCorrection()->GetGeneratedHistogram(), ptDists[i]);
      }*/
      return "PythiaRatios";
      break;

      // one species enhanced / reduced
    case 2: // + 50% pions
    case 3: // - 50% pions
    case 4: // + 50% kaons
    case 5: // - 50% kaons
    case 6: // + 50% protons
    case 7: // - 50% protons
      Int_t correctionIndex = (index - 2) / 2;
      Double_t scaleFactor = (index % 2 == 0) ? 1.5 : 0.5;

      fdNdEtaCorrection[correctionIndex]->GetTrack2ParticleCorrection()->GetMeasuredHistogram()->Scale(scaleFactor);
      fdNdEtaCorrection[correctionIndex]->GetTrack2ParticleCorrection()->GetGeneratedHistogram()->Scale(scaleFactor);

      TString* str = new TString;
      str->Form("%s%s", (correctionIndex == 0) ? "Pi" : ((correctionIndex == 1) ? "K" : "p"), (index % 2 == 0) ? "Boosted" : "Reduced");
      return str->Data();
      break;

      // each species rated with pythia ratios
    case 12: // + 50% pions
    case 13: // - 50% pions
    case 14: // + 50% kaons
    case 15: // - 50% kaons
    case 16: // + 50% protons
    case 17: // - 50% protons
      TH1** ptDists = DrawRatios("pythiaratios.root");
      Int_t functionIndex = (index - 2) / 2;
      Double_t scaleFactor = (index % 2 == 0) ? 1.5 : 0.5;
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

  TH1F* ratio = new TH1F("Sigma2Vertex_ratio", "Sigma2Vertex_ratio;n sigma;included", 50, 0.05, 5.05);
  for (Double_t nSigma = 0.1; nSigma < 5.05; nSigma += 0.1)
    ratio->Fill(nSigma, Sigma2VertexCount(tracks, nSigma));
  ratio->SetMarkerStyle(21);

  canvas->cd(2);
  ratio->DrawCopy("P");

  TH1F* ratio2 = new TH1F("Sigma2Vertex_ratio2", "Sigma2Vertex_ratio2;nSigma;% included 3 sigma / % included n sigma", 50, 0.05, 5.05);
  Double_t sigma3 = Sigma2VertexCount(tracks, 3);
  for (Double_t nSigma = 0.1; nSigma < 5.05; nSigma += 0.1)
    ratio2->Fill(nSigma, sigma3 / ratio->GetBinContent(ratio->FindBin(nSigma)));
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

  TH1F* ratio = new TH1F("sigmavertexsimulation_ratio", "sigmavertexsimulation_ratio;N#sigma;% included in 3 #sigma / % included in N#sigma", sigmavertex->GetNbinsX(), sigmavertex->GetXaxis()->GetXmin(), sigmavertex->GetXaxis()->GetXmax());

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

  ratio1->SetLineWidth(2);
  ratio2->SetLineWidth(2);

  TLegend* legend = new TLegend(0.7, 0.8, 0.95, 0.95);
  legend->SetFillColor(0);
  legend->AddEntry(ratio1, "Gaussian");
  legend->AddEntry(ratio2, "Simulation");

  ratio2->SetTitle("");
  ratio2->GetYaxis()->SetTitleOffset(1.5);
  ratio2->GetXaxis()->SetRangeUser(2, 4);

  TCanvas* canvas = new TCanvas("Sigma2VertexCompare", "Sigma2VertexCompare", 500, 500);
  InitPad();

  ratio2->SetMarkerStyle(21);
  ratio1->SetMarkerStyle(22);

  ratio2->GetYaxis()->SetRangeUser(0.8, 1.2);
  ratio2->SetLineColor(kRed);
  ratio2->SetMarkerColor(kRed);
  ratio2->Draw("PL");
  ratio1->Draw("SAMEPL");

  legend->Draw();

  canvas->SaveAs("Sigma2VertexCompare.eps");
}

void drawSystematics()
{
  //Secondaries();
  //DrawDifferentSpecies();
  //Composition();

  Sigma2VertexSimulation();

}

void DrawdNdEtaDifferences()
{
  TH1* hists[5];

  TLegend* legend = new TLegend(0.3, 0.73, 0.70, 0.98);
  legend->SetFillColor(0);

  TCanvas* canvas = new TCanvas("DrawdNdEtaDifferences", "DrawdNdEtaDifferences", 1000, 500);
  canvas->Divide(2, 1);

  canvas->cd(1);

  for (Int_t i=0; i<5; ++i)
  {
    hists[i] = 0;
    TFile* file = 0;
    TString title;

    switch(i)
    {
      case 0 : file = TFile::Open("systematics_dndeta_reference.root"); title = "standard composition"; break;
      case 1 : file = TFile::Open("systematics_dndeta_KBoosted.root"); title = "+ 50% kaons"; break;
      case 2 : file = TFile::Open("systematics_dndeta_KReduced.root"); title = "- 50% kaons"; break;
      case 3 : file = TFile::Open("systematics_dndeta_pBoosted.root"); title = "+ 50% protons"; break;
      case 4 : file = TFile::Open("systematics_dndeta_pReduced.root"); title = "- 50% protons"; break;
      default: return;
    }

    if (file)
    {
      hists[i] = (TH1*) file->Get("dndeta/dndeta_dNdEta_corrected_2");
      hists[i]->SetTitle("a)");

      Prepare1DPlot(hists[i], kFALSE);
      hists[i]->GetXaxis()->SetRangeUser(-0.7999, 0.7999);
      hists[i]->GetYaxis()->SetRangeUser(6, 7);
      hists[i]->SetLineColor(colors[i]);
      hists[i]->SetMarkerColor(colors[i]);
      hists[i]->SetMarkerStyle(markers[i]);
      hists[i]->GetXaxis()->SetLabelOffset(0.015);
      hists[i]->GetYaxis()->SetTitleOffset(1.5);
      gPad->SetLeftMargin(0.12);
      hists[i]->DrawCopy(((i > 0) ? "SAME" : ""));

      legend->AddEntry(hists[i], title);
      hists[i]->SetTitle(title);
    }
  }
  legend->Draw();

  canvas->cd(2);
  gPad->SetLeftMargin(0.14);

  TLegend* legend2 = new TLegend(0.73, 0.73, 0.98, 0.98);
  legend2->SetFillColor(0);

  for (Int_t i=1; i<5; ++i)
  {
    if (hists[i])
    {
      legend2->AddEntry(hists[i]);

      hists[i]->Divide(hists[0]);
      hists[i]->SetTitle("b)");
      hists[i]->SetLineColor(colors[i-1]);
      hists[i]->SetMarkerColor(colors[i-1]);
      hists[i]->GetYaxis()->SetRangeUser(0.95, 1.05);
      hists[i]->GetYaxis()->SetTitle("Ratio to standard composition");
      hists[i]->GetYaxis()->SetTitleOffset(1.8);
      hists[i]->DrawCopy(((i > 1) ? "SAME" : ""));
    }
  }

  legend2->Draw();

  canvas->SaveAs("particlecomposition_result_detail.gif");

  TCanvas* canvas2 = new TCanvas("DrawdNdEtaDifferences2", "DrawdNdEtaDifferences2", 700, 500);

  for (Int_t i=1; i<5; ++i)
  {
    if (hists[i])
    {
      hists[i]->SetTitle("");
      hists[i]->GetYaxis()->SetTitleOffset(1.1);
      hists[i]->DrawCopy(((i > 1) ? "SAME" : ""));
    }
  }

  legend2->Draw();

  canvas2->SaveAs("particlecomposition_result.gif");
  canvas2->SaveAs("particlecomposition_result.eps");
}

mergeCorrectionsWithDifferentCrosssections(Char_t* standardCorrectionFileName="correction_map.root",
					     Char_t* systematicCorrectionFileName="systematics.root",
					     Char_t* outputFileName="systematics_vtxtrigger_compositions.root") {
  //
  // Function used to merge standard corrections with vertex
  // reconstruction corrections obtained by a certain mix of ND, DD
  // and SD events.
  // 

  gSystem->Load("libPWG0base");

  const Char_t* typeName[] = { "vertexreco", "trigger", "vtxtrigger" };
  const Char_t* changes[]  = {"pythia","ddmore","ddless","sdmore","sdless", "dmore", "dless"};
  Float_t scalesDD[] = {1.0, 1.5, 0.5, 1.0, 1.0, 1.5, 0.5};
  Float_t scalesSD[] = {1.0, 1.0, 1.0, 1.5, 0.5, 1.5, 0.5};

  // cross section from Pythia
  Float_t sigmaND = 55.2;
  Float_t sigmaDD = 9.78;
  Float_t sigmaSD = 14.30;

  AlidNdEtaCorrection* corrections[21];
  for (Int_t j=0; j<3; j++) { // j = 0 (change vtx), j = 1 (change trg), j = 2 (change both)
    AlidNdEtaCorrection* correctionStandard = new AlidNdEtaCorrection("dndeta_correction","dndeta_correction");
    correctionStandard->LoadHistograms(standardCorrectionFileName);

    // dont take vertexreco from this one
    correctionStandard->GetVertexRecoCorrection()->Reset();
    // dont take triggerbias from this one
    correctionStandard->GetTriggerBiasCorrectionINEL()->Reset();

    for (Int_t i=0; i<7; i++) {
      TString name;
      name.Form("dndeta_correction_syst_%s_%s", typeName[j], changes[i]);
      AlidNdEtaCorrection* current = new AlidNdEtaCorrection(name, name);

      name.Form("vertexRecoND");
      AlidNdEtaCorrection* dNdEtaCorrectionVtxND = new AlidNdEtaCorrection(name,name);
      dNdEtaCorrectionVtxND->LoadHistograms(systematicCorrectionFileName, name);
      name.Form("vertexRecoDD");
      AlidNdEtaCorrection* dNdEtaCorrectionVtxDD = new AlidNdEtaCorrection(name,name);
      dNdEtaCorrectionVtxDD->LoadHistograms(systematicCorrectionFileName, name);
      name.Form("vertexRecoSD");
      AlidNdEtaCorrection* dNdEtaCorrectionVtxSD = new AlidNdEtaCorrection(name,name);
      dNdEtaCorrectionVtxSD->LoadHistograms(systematicCorrectionFileName, name);

      name.Form("triggerBiasND");
      AlidNdEtaCorrection* dNdEtaCorrectionTriggerND = new AlidNdEtaCorrection(name,name);
      dNdEtaCorrectionTriggerND->LoadHistograms(systematicCorrectionFileName, name);
      name.Form("triggerBiasDD");
      AlidNdEtaCorrection* dNdEtaCorrectionTriggerDD = new AlidNdEtaCorrection(name,name);
      dNdEtaCorrectionTriggerDD->LoadHistograms(systematicCorrectionFileName, name);
      name.Form("triggerBiasSD");
      AlidNdEtaCorrection* dNdEtaCorrectionTriggerSD = new AlidNdEtaCorrection(name,name);
      dNdEtaCorrectionTriggerSD->LoadHistograms(systematicCorrectionFileName, name);

      // calculating relative
      Float_t nd = 100 * sigmaND/(sigmaND + (scalesDD[i]*sigmaDD) + (scalesDD[i]*sigmaSD));
      Float_t dd = 100 * (scalesDD[i]*sigmaDD)/(sigmaND + (scalesDD[i]*sigmaDD) + (scalesDD[i]*sigmaSD));
      Float_t sd = 100 * (scalesSD[i]*sigmaSD)/(sigmaND + (scalesDD[i]*sigmaDD) + (scalesDD[i]*sigmaSD));

      printf(Form("%s : ND=%.1f\%, DD=%.1f\%, SD=%.1f\% \n",changes[i],nd,dd,sd));
      current->SetTitle(Form("ND=%.2f\%,DD=%.2f\%,SD=%.2f\%",nd,dd,sd));
      current->SetTitle(name);

      // scale
      if (j == 0 || j == 2)
      {
        dNdEtaCorrectionVtxDD->GetVertexRecoCorrection()->GetMeasuredHistogram()->Scale(scalesDD[i]);
        dNdEtaCorrectionVtxSD->GetVertexRecoCorrection()->GetMeasuredHistogram()->Scale(scalesSD[i]);
        dNdEtaCorrectionVtxDD->GetVertexRecoCorrection()->GetGeneratedHistogram()->Scale(scalesDD[i]);
        dNdEtaCorrectionVtxSD->GetVertexRecoCorrection()->GetGeneratedHistogram()->Scale(scalesSD[i]);
      }
      if (j == 1 || j == 2)
      {
        dNdEtaCorrectionTriggerDD->GetTriggerBiasCorrectionINEL()->GetMeasuredHistogram()->Scale(scalesDD[i]);
        dNdEtaCorrectionTriggerSD->GetTriggerBiasCorrectionINEL()->GetMeasuredHistogram()->Scale(scalesSD[i]);
        dNdEtaCorrectionTriggerDD->GetTriggerBiasCorrectionINEL()->GetGeneratedHistogram()->Scale(scalesDD[i]);
        dNdEtaCorrectionTriggerSD->GetTriggerBiasCorrectionINEL()->GetGeneratedHistogram()->Scale(scalesSD[i]);
      }

      //clear track, trigger in Vtx correction
      dNdEtaCorrectionVtxND->GetTrack2ParticleCorrection()->Reset();
      dNdEtaCorrectionVtxND->GetTriggerBiasCorrectionNSD()->Reset();
      dNdEtaCorrectionVtxND->GetTriggerBiasCorrectionND()->Reset();
      dNdEtaCorrectionVtxND->GetTriggerBiasCorrectionINEL()->Reset();
      dNdEtaCorrectionVtxSD->GetTrack2ParticleCorrection()->Reset();
      dNdEtaCorrectionVtxSD->GetTriggerBiasCorrectionNSD()->Reset();
      dNdEtaCorrectionVtxSD->GetTriggerBiasCorrectionND()->Reset();
      dNdEtaCorrectionVtxSD->GetTriggerBiasCorrectionINEL()->Reset();
      dNdEtaCorrectionVtxDD->GetTrack2ParticleCorrection()->Reset();
      dNdEtaCorrectionVtxDD->GetTriggerBiasCorrectionNSD()->Reset();
      dNdEtaCorrectionVtxDD->GetTriggerBiasCorrectionND()->Reset();
      dNdEtaCorrectionVtxDD->GetTriggerBiasCorrectionINEL()->Reset();

      //clear track, vertexreco in trigger correction
      dNdEtaCorrectionTriggerND->GetTrack2ParticleCorrection()->Reset();
      dNdEtaCorrectionTriggerND->GetTriggerBiasCorrectionNSD()->Reset();
      dNdEtaCorrectionTriggerND->GetTriggerBiasCorrectionND()->Reset();
      dNdEtaCorrectionTriggerND->GetVertexRecoCorrection()->Reset();
      dNdEtaCorrectionTriggerSD->GetTrack2ParticleCorrection()->Reset();
      dNdEtaCorrectionTriggerSD->GetTriggerBiasCorrectionNSD()->Reset();
      dNdEtaCorrectionTriggerSD->GetTriggerBiasCorrectionND()->Reset();
      dNdEtaCorrectionTriggerSD->GetVertexRecoCorrection()->Reset();
      dNdEtaCorrectionTriggerDD->GetTrack2ParticleCorrection()->Reset();
      dNdEtaCorrectionTriggerDD->GetTriggerBiasCorrectionNSD()->Reset();
      dNdEtaCorrectionTriggerDD->GetTriggerBiasCorrectionND()->Reset();
      dNdEtaCorrectionTriggerDD->GetVertexRecoCorrection()->Reset();

      TList collection;
      collection.Add(correctionStandard);
      collection.Add(dNdEtaCorrectionVtxND);
      collection.Add(dNdEtaCorrectionVtxDD);
      collection.Add(dNdEtaCorrectionVtxSD);
      collection.Add(dNdEtaCorrectionTriggerND);
      collection.Add(dNdEtaCorrectionTriggerDD);
      collection.Add(dNdEtaCorrectionTriggerSD);

      current->Merge(&collection);
      current->Finish();

      corrections[i+j*7] = current;
    }
  }

  TFile* fout = new TFile(outputFileName,"RECREATE");

  for (Int_t i=0; i<21; i++)
    corrections[i]->SaveHistograms();

  fout->Write();
  fout->Close();
}


DrawVertexRecoSyst() {
  // Draws the ratio of the dN/dEta obtained with changed SD and DD
  // cross-sections vertex reco corrections to the dN/dEta obtained
  // from the standard pythia cross-sections 
  //
  // The files with the vertex reco corrections for different
  // processes (and with the other standard corrections) are expected
  // to have the names "analysis_esd_X.root", where the Xs are defined
  // in the array changes below.

  Char_t* changes[]  = {"pythia","ddmore","ddless","sdmore","sdless", "dmore", "dless"};
  Char_t* descr[]  =   {"",
			"#sigma_{DD}' = 1.5#times#sigma_{DD}",
			"#sigma_{DD}' = 0.5#times#sigma_{DD}",
			"#sigma_{SD}' = 1.5#times#sigma_{SD}",
			"#sigma_{SD}' = 0.5#times#sigma_{SD}",
			"#sigma_{D}' = 1.5#times#sigma_{D}",
			"#sigma_{D}' = 0.5#times#sigma_{D}"};

  Float_t scalesDD[] = {1.0, 1.5, 0.5, 1.5, 0.5, 1.5, 0.5};
  Float_t scalesSD[] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.5, 0.5};

  Int_t markers[] = {20,20,21,22,23,28,29};
  Int_t colors[]  = {1,2,3,4,6,8,102};

  // cross section from Pythia
  Float_t sigmaND = 55.2;
  Float_t sigmaDD = 9.78;
  Float_t sigmaSD = 14.30;

  TH1F* dNdEta[7];
  TH1F* ratios[7];

  TFile* fin;

  for (Int_t i=0; i<7; i++) {
    // calculating relative
    fin = TFile::Open(Form("analysis_esd_%s.root",changes[i]));

    dNdEta[i] = (TH1F*)(fin->Get("dndeta/dndeta_dNdEta_corrected_2"))->Clone();

    for (Int_t b=0; b<dNdEta[i]->GetNbinsX(); b++) {
      if (TMath::Abs(dNdEta[i]->GetBinCenter(b))>0.9) {
	dNdEta[i]->SetBinContent(b,0);
	dNdEta[i]->SetBinError(b,0);
      }
    }

    dNdEta[i]->Rebin(4);

    dNdEta[i]->SetLineWidth(2);
    dNdEta[i]->SetLineColor(colors[i]);
    dNdEta[i]->SetMarkerStyle(markers[i]);
    dNdEta[i]->SetMarkerSize(0.9);
    dNdEta[i]->SetMarkerColor(colors[i]);

    ratios[i] = (TH1F*)dNdEta[i]->Clone("ratio");
    ratios[i]->Divide(ratios[i],dNdEta[0],1,1,"B");
    
    ratios[i]->SetName(changes[i]);
    ratios[i]->SetTitle(changes[i]);
  }
  
  //##########################################################
  
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetOptFit(0);

  gStyle->SetTextSize(0.2);
  gStyle->SetTitleSize(0.05,"xyz");
  gStyle->SetLabelOffset(0.01, "xyz");


  gStyle->SetTitleOffset(1.2, "y");
  gStyle->SetTitleOffset(1.2, "x");
  gStyle->SetEndErrorSize(0.0);

  //##############################################

  //making canvas and pads
  TCanvas *c = new TCanvas(Form("vertex_reco_syst_%s", plot), Form("vertex_reco_syst_%s", plot),600,500);

  TPad* p1 = new TPad("pad1","", 0, 0.0, 1.0, 1.0, 0, 0, 0);

  p1->SetBottomMargin(0.15);
  p1->SetTopMargin(0.03);
  p1->SetLeftMargin(0.15);
  p1->SetRightMargin(0.03);

  p1->SetGridx();
  p1->SetGridy();

  p1->Draw();
  p1->cd();
  
  
  TH2F* null = new TH2F("","",100,-1.5,1.5,100,0.9601,1.099);
  null->SetXTitle("#eta");
  null->GetXaxis()->CenterTitle(kTRUE);
  null->SetYTitle("dN/d#eta / dN/d#eta_{pythia}");
  null->GetYaxis()->CenterTitle(kTRUE);
  null->Draw();

  for (Int_t i=1; i<7; i++) 
    ratios[i]->Draw("same");

  TLegend* legend = new TLegend(0.6, 0.6, 0.95, 0.95);
  legend->SetFillColor(0);

  TLatex* text[7];
  for (Int_t i=1; i<7; i++) {
    legend->AddEntry(dNdEta[i], descr[i]);
  }

  legend->Draw();
  //text(0.2,0.88,"Effect of changing",0.045,1,kTRUE);
  //text(0.2,0.83,"relative cross-sections",0.045,1,kTRUE);
  //text(0.2,0.78,"(vertex reconstruction corr.)",0.043,13,kTRUE);

  c->SaveAs(Form("%s.gif", c->GetName()));
  c->SaveAs(Form("%s.eps", c->GetName()));
}



DrawTriggerEfficiency(Char_t* fileName) {

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetOptFit(0);

  gStyle->SetTextSize(0.04);
  gStyle->SetTitleSize(0.05,"xyz");
  //gStyle->SetTitleFont(133, "xyz");
  //gStyle->SetLabelFont(133, "xyz");
  //gStyle->SetLabelSize(17, "xyz");
  gStyle->SetLabelOffset(0.01, "xyz");

  gStyle->SetTitleOffset(1.1, "y");
  gStyle->SetTitleOffset(1.1, "x");
  gStyle->SetEndErrorSize(0.0);

  //##############################################

  //making canvas and pads
  TCanvas *c = new TCanvas("trigger_eff", "",600,500);

  TPad* p1 = new TPad("pad1","", 0, 0.0, 1.0, 1.0, 0, 0, 0);

  p1->SetBottomMargin(0.15);
  p1->SetTopMargin(0.03);
  p1->SetLeftMargin(0.15);
  p1->SetRightMargin(0.03);
  
  p1->SetGridx();
  p1->SetGridy();

  p1->Draw();
  p1->cd();

  gSystem->Load("libPWG0base");

  AlidNdEtaCorrection* corrections[4];
  AliCorrectionMatrix2D* triggerBiasCorrection[4];

  TH1F* hTriggerEffInv[4];
  TH1F* hTriggerEff[4];

  Char_t* names[] = {"triggerBiasND", "triggerBiasDD", "triggerBiasSD", "triggerBiasINEL"};
  
  for (Int_t i=0; i<4; i++) {
    corrections[i] = new AlidNdEtaCorrection(names[i], names[i]);
    corrections[i]->LoadHistograms(fileName, names[i]);    

    triggerBiasCorrection[i] = corrections[i]->GetTriggerBiasCorrectionINEL();

    
    hTriggerEffInv[i] = triggerBiasCorrection[i]->Get1DCorrection();
    hTriggerEff[i]    = (TH1F*)hTriggerEffInv[i]->Clone();
    
    for (Int_t b=0; b<=hTriggerEff[i]->GetNbinsX(); b++) {
      hTriggerEff[i]->SetBinContent(b,1);
      hTriggerEff[i]->SetBinError(b,0);
    }
    hTriggerEff[i]->Divide(hTriggerEff[i],hTriggerEffInv[i]);
    hTriggerEff[i]->Scale(100);
  }

  Int_t colors[] = {2,3,4,1};
  Float_t effs[4];
  for (Int_t i=0; i<4; i++) {
    hTriggerEff[i]->Fit("pol0","","",-20,20);
    effs[i] = ((TF1*)hTriggerEff[i]->GetListOfFunctions()->At(0))->GetParameter(0);
    ((TF1*)hTriggerEff[i]->GetListOfFunctions()->At(0))->SetLineWidth(1);
    ((TF1*)hTriggerEff[i]->GetListOfFunctions()->At(0))->SetLineStyle(2);
    ((TF1*)hTriggerEff[i]->GetListOfFunctions()->At(0))->SetLineColor(colors[i]);
    cout << effs[i] << endl;
  }


  Char_t* text[] = {"ND", "DD", "SD", "INEL"};
  TLatex* latex[4];

  TH2F* null = new TH2F("","",100,-25,35,100,0,110);
  null->SetXTitle("Vertex z [cm]");
  null->GetXaxis()->CenterTitle(kTRUE);
  null->SetYTitle("Trigger efficiency [%]");
  null->GetYaxis()->CenterTitle(kTRUE);
  null->Draw();


  for (Int_t i=0; i<4; i++) {
    hTriggerEff[i]->SetLineWidth(2);
    hTriggerEff[i]->SetLineColor(colors[i]);

    hTriggerEff[i]->Draw("same");

    latex[i] = new TLatex(22,effs[i]-1.5, Form("%s (%0.1f)",text[i],effs[i]));
    latex[i]->SetTextColor(colors[i]);
    latex[i]->Draw();
  }
  
}


DrawSpectraPID(Char_t* fileName) {

  gSystem->Load("libPWG0base");

  Char_t* names[] = {"correction_0", "correction_1", "correction_2", "correction_3"};
  AlidNdEtaCorrection* corrections[4];
  AliCorrectionMatrix3D* trackToPartCorrection[4];

  TH1F* measuredPt[4];
  TH1F* generatedPt[4];
  TH1F* ratioPt[4];

  for (Int_t i=0; i<4; i++) {
    corrections[i] = new AlidNdEtaCorrection(names[i], names[i]);
    corrections[i]->LoadHistograms(fileName, names[i]);    

    trackToPartCorrection[i] = corrections[i]->GetTrack2ParticleCorrection();

    Int_t binX1 = (TH1F*)trackToPartCorrection[i]->GetMeasuredHistogram()->GetXaxis()->FindBin(-10);
    Int_t binX2 = (TH1F*)trackToPartCorrection[i]->GetMeasuredHistogram()->GetXaxis()->FindBin(10);
    Int_t binY1 = (TH1F*)trackToPartCorrection[i]->GetMeasuredHistogram()->GetYaxis()->FindBin(-1);
    Int_t binY2 = (TH1F*)trackToPartCorrection[i]->GetMeasuredHistogram()->GetYaxis()->FindBin(1);

    measuredPt[i]  = (TH1F*)trackToPartCorrection[i]->GetMeasuredHistogram()->ProjectionZ(Form("m_%d",i),binX1,binX2,binY1,binY2);
    generatedPt[i] = (TH1F*)trackToPartCorrection[i]->GetGeneratedHistogram()->ProjectionZ(Form("g_%d",i),binX1,binX2,binY1,binY2);
    ratioPt[i] = (TH1F*)generatedPt[i]->Clone(Form("r_%d",i));
    ratioPt[i]->Divide(measuredPt[i], generatedPt[i], 1,1,"B");
  }
  
  ratioPt[0]->Draw();
  
  for (Int_t i=0; i<3; i++) {
    ratioPt[i]->SetLineColor(i+1);
    ratioPt[i]->SetLineWidth(2);
    
    ratioPt[i]->Draw("same");
    
  }

  return;
  measuredPt[0]->SetLineColor(2);
  measuredPt[0]->SetLineWidth(5);

  measuredPt[0]->Draw();
  generatedPt[0]->Draw("same");
}

void changePtSpectrum()
{
  TFile* file = TFile::Open("analysis_mc.root");
  TH1F* hist = dynamic_cast<TH1F*> (file->Get("dndeta_check_pt"));

  //hist->Rebin(3);
  //hist->Scale(1.0/3);

  TH1F* clone1 = dynamic_cast<TH1F*> (hist->Clone("clone1"));
  TH1F* clone2 = dynamic_cast<TH1F*> (hist->Clone("clone2"));

  TH1F* scale1 =  dynamic_cast<TH1F*> (hist->Clone("scale1"));
  TH1F* scale2 =  dynamic_cast<TH1F*> (hist->Clone("scale2"));

  Float_t ptCutOff = 0.3;

  for (Int_t i=1; i <= hist->GetNbinsX(); ++i)
  {
    if (hist->GetBinCenter(i) > ptCutOff)
    {
      scale1->SetBinContent(i, 1);
      scale2->SetBinContent(i, 1);
    }
    else
    {
      // 90 % at pt = 0, 0% at pt = ptcutoff
      scale1->SetBinContent(i, 1 - (ptCutOff - hist->GetBinCenter(i)) / ptCutOff * 0.3);

      // 110% at pt = 0, ...
      scale2->SetBinContent(i, 1 + (ptCutOff - hist->GetBinCenter(i)) / ptCutOff * 0.3);
    }
    scale1->SetBinError(i, 0);
    scale2->SetBinError(i, 0);
  }

  new TCanvas;

  scale1->Draw();
  scale2->SetLineColor(kRed);
  scale2->Draw("SAME");

  clone1->Multiply(scale1);
  clone2->Multiply(scale2);

  Prepare1DPlot(hist);
  Prepare1DPlot(clone1);
  Prepare1DPlot(clone2);

  /*hist->SetMarkerStyle(markers[0]);
  clone1->SetMarkerStyle(markers[0]);
  clone2->SetMarkerStyle(markers[0]);*/

  hist->SetTitle(";p_{T} in GeV/c;dN/dp_{T} in c/GeV");
  hist->GetXaxis()->SetRangeUser(0, 0.7);
  hist->GetYaxis()->SetRangeUser(0.01, clone2->GetMaximum() * 1.1);
  hist->GetYaxis()->SetTitleOffset(1);

  TCanvas* canvas = new TCanvas;
  gPad->SetBottomMargin(0.12);
  hist->Draw("H");
  clone1->SetLineColor(kRed);
  clone1->Draw("HSAME");
  clone2->SetLineColor(kBlue);
  clone2->Draw("HSAME");
  hist->Draw("HSAME");

  Float_t fraction =  hist->Integral(hist->GetXaxis()->FindBin(ptCutOff), hist->GetNbinsX()) / hist->Integral(1, hist->GetNbinsX());
  Float_t fraction1 = clone1->Integral(clone1->GetXaxis()->FindBin(ptCutOff), clone1->GetNbinsX()) / clone1->Integral(1, clone1->GetNbinsX());
  Float_t fraction2 = clone2->Integral(clone2->GetXaxis()->FindBin(ptCutOff), clone2->GetNbinsX()) / clone2->Integral(1, clone2->GetNbinsX());

  printf("%f %f %f\n", fraction, fraction1, fraction2);
  printf("Rel. %f %f\n", fraction1 / fraction, fraction2 / fraction);

  canvas->SaveAs("changePtSpectrum.gif");
  canvas->SaveAs("changePtSpectrum.eps");
}

void FractionBelowPt()
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

  Double_t geneCount[5];
  Double_t measCount[5];
  geneCount[4] = 0;
  measCount[4] = 0;

  for (Int_t i=0; i<4; ++i)
  {
    TH3F* hist = (TH3F*) fdNdEtaCorrection[i]->GetTrack2ParticleCorrection()->GetGeneratedHistogram();
    geneCount[i] = hist->Integral(hist->GetXaxis()->FindBin(-10), hist->GetXaxis()->FindBin(10),
                                  hist->GetYaxis()->FindBin(-0.8), hist->GetYaxis()->FindBin(0.8),
                                  1, hist->GetZaxis()->FindBin(0.3));

    hist = (TH3F*) fdNdEtaCorrection[i]->GetTrack2ParticleCorrection()->GetMeasuredHistogram();
    measCount[i] = hist->Integral(hist->GetXaxis()->FindBin(-10), hist->GetXaxis()->FindBin(10), hist->GetYaxis()->FindBin(-0.8), hist->GetYaxis()->FindBin(0.8), 1, hist->GetZaxis()->FindBin(0.3));

    geneCount[4] += geneCount[i];
    measCount[4] += measCount[i];

    printf("Particle %s: %d gene, %d meas\n", ((i == 0) ? "pi" : (i == 1) ? "K" : (i == 2) ? "p" : "others"), (Int_t) geneCount[i], (Int_t) measCount[i]);
  }

  printf("Generated ratios are:     %f pi, %f K, %f p, %f others\n", geneCount[0] / geneCount[4], geneCount[1] / geneCount[4], geneCount[2] / geneCount[4], geneCount[3] / geneCount[4]);

  printf("Reconstructed ratios are: %f pi, %f K, %f p, %f others\n", measCount[0] / measCount[4], measCount[1] / measCount[4], measCount[2] / measCount[4], measCount[3] / measCount[4]);
}


mergeCorrectionsMisalignment(Char_t* alignedFile = "correction_map_aligned.root",
					     Char_t* misalignedFile = "correction_map_misaligned.root",
					     Char_t* outputFileName="correction_map_misaligned_single.root")
{
  //
  // from the aligned and misaligned corrections, 3 new corrections are created
  // in these new corrections only one of the corrections (track2particle, vertex, trigger)
  // is taken from the misaligned input to allow study of the effect on the different
  // corrections

  gSystem->Load("libPWG0base");

  const Char_t* typeName[] = { "track2particle", "vertex", "trigger" };

  AlidNdEtaCorrection* corrections[3];
  for (Int_t j=0; j<3; j++) { // j = 0 (track2particle), j = 1 (vertex), j = 2 (trigger)
    AlidNdEtaCorrection* alignedCorrection = new AlidNdEtaCorrection("dndeta_correction","dndeta_correction");
    alignedCorrection->LoadHistograms(alignedFile);

    AlidNdEtaCorrection* misalignedCorrection = new AlidNdEtaCorrection("dndeta_correction","dndeta_correction");
    misalignedCorrection->LoadHistograms(misalignedFile);

    TString name;
    name.Form("dndeta_correction_alignment_%s", typeName[j]);
    AlidNdEtaCorrection* current = new AlidNdEtaCorrection(name, name);

    switch (j)
    {
      case 0:
        alignedCorrection->GetTrack2ParticleCorrection()->Reset();
        misalignedCorrection->GetTriggerBiasCorrectionNSD()->Reset();
        misalignedCorrection->GetTriggerBiasCorrectionND()->Reset();
        misalignedCorrection->GetTriggerBiasCorrectionINEL()->Reset();
        misalignedCorrection->GetVertexRecoCorrection()->Reset();
        break;

      case 1:
        misalignedCorrection->GetTrack2ParticleCorrection()->Reset();
        misalignedCorrection->GetTriggerBiasCorrectionNSD()->Reset();
        misalignedCorrection->GetTriggerBiasCorrectionND()->Reset();
        misalignedCorrection->GetTriggerBiasCorrectionINEL()->Reset();
        alignedCorrection->GetVertexRecoCorrection()->Reset();
        break;

      case 2:
        misalignedCorrection->GetTrack2ParticleCorrection()->Reset();
        alignedCorrection->GetTriggerBiasCorrectionNSD()->Reset();
        alignedCorrection->GetTriggerBiasCorrectionND()->Reset();
        alignedCorrection->GetTriggerBiasCorrectionINEL()->Reset();
        misalignedCorrection->GetVertexRecoCorrection()->Reset();
        break;

      default:
        return;
    }

    TList collection;
    collection.Add(misalignedCorrection);
    collection.Add(alignedCorrection);

    current->Merge(&collection);
    current->Finish();

    corrections[j] = current;
  }

  TFile* fout = new TFile(outputFileName, "RECREATE");

  for (Int_t i=0; i<3; i++)
    corrections[i]->SaveHistograms();

  fout->Write();
  fout->Close();
}


void DrawdNdEtaDifferencesAlignment()
{
  TH1* hists[5];

  TLegend* legend = new TLegend(0.3, 0.73, 0.70, 0.98);
  legend->SetFillColor(0);

  TCanvas* canvas = new TCanvas("DrawdNdEtaDifferences", "DrawdNdEtaDifferences", 1000, 500);
  canvas->Divide(2, 1);

  canvas->cd(1);

  for (Int_t i=0; i<5; ++i)
  {
    hists[i] = 0;
    TFile* file = 0;
    TString title;

    switch(i)
    {
      case 0 : file = TFile::Open("systematics_misalignment_aligned.root"); title = "aligned"; break;
      case 1 : file = TFile::Open("systematics_misalignment_misaligned.root"); title = "fully misaligned"; break;
      case 2 : file = TFile::Open("systematics_misalignment_track2particle.root"); title = "only track2particle"; break;
      case 3 : file = TFile::Open("systematics_misalignment_vertex.root"); title = "only vertex rec."; break;
      case 4 : file = TFile::Open("systematics_misalignment_trigger.root"); title = "only trigger bias"; break;
      default: return;
    }

    if (file)
    {
      hists[i] = (TH1*) file->Get("dndeta/dndeta_dNdEta_corrected_2");
      hists[i]->SetTitle("");

      Prepare1DPlot(hists[i], kFALSE);
      hists[i]->GetXaxis()->SetRangeUser(-0.7999, 0.7999);
      hists[i]->GetYaxis()->SetRangeUser(6, 7);
      hists[i]->SetLineWidth(1);
      hists[i]->SetLineColor(colors[i]);
      hists[i]->SetMarkerColor(colors[i]);
      hists[i]->SetMarkerStyle(markers[i]);
      hists[i]->GetXaxis()->SetLabelOffset(0.015);
      hists[i]->GetYaxis()->SetTitleOffset(1.5);
      gPad->SetLeftMargin(0.12);
      hists[i]->DrawCopy(((i > 0) ? "SAME" : ""));

      legend->AddEntry(hists[i], title);
      hists[i]->SetTitle(title);
    }
  }
  legend->Draw();

  canvas->cd(2);
  gPad->SetLeftMargin(0.14);

  TLegend* legend2 = new TLegend(0.63, 0.73, 0.98, 0.98);
  legend2->SetFillColor(0);

  for (Int_t i=1; i<5; ++i)
  {
    if (hists[i])
    {
      legend2->AddEntry(hists[i]);

      hists[i]->Divide(hists[0]);
      hists[i]->SetTitle("b)");
      hists[i]->GetYaxis()->SetRangeUser(0.9, 1.1);
      hists[i]->GetYaxis()->SetTitle("Ratio to standard composition");
      hists[i]->GetYaxis()->SetTitleOffset(1.8);
      hists[i]->DrawCopy(((i > 1) ? "SAME" : ""));
    }
  }

  legend2->Draw();

  canvas->SaveAs("misalignment_result.eps");
  canvas->SaveAs("misalignment_result.gif");
}

void EvaluateMultiplicityEffect()
{
  gSystem->Load("libPWG0base");

  AlidNdEtaCorrection* fdNdEtaCorrectionLow[4];
  AlidNdEtaCorrection* fdNdEtaCorrectionHigh[4];

  TFile::Open("systematics-low-multiplicity.root");

  for (Int_t i=0; i<4; ++i)
  {
    TString name;
    name.Form("correction_%d", i);
    fdNdEtaCorrectionLow[i] = new AlidNdEtaCorrection(name, name);
    fdNdEtaCorrectionLow[i]->LoadHistograms("systematics-low-multiplicity.root", name);
  }

  TList list;
  for (Int_t i=1; i<4; ++i)
    list.Add(fdNdEtaCorrectionLow[i]);

  fdNdEtaCorrectionLow[0]->Merge(&list);
  fdNdEtaCorrectionLow[0]->Finish();

  TFile::Open("systematics-high-multiplicity.root");

  for (Int_t i=0; i<4; ++i)
  {
    TString name;
    name.Form("correction_%d", i);
    fdNdEtaCorrectionHigh[i] = new AlidNdEtaCorrection(name, name);
    fdNdEtaCorrectionHigh[i]->LoadHistograms("systematics-high-multiplicity.root", name);
  }

  TList list2;
  for (Int_t i=1; i<4; ++i)
    list2.Add(fdNdEtaCorrectionHigh[i]);

  fdNdEtaCorrectionHigh[0]->Merge(&list2);
  fdNdEtaCorrectionHigh[0]->Finish();

  TH1F* outputLow = new TH1F("Track2ParticleLow", "Track2Particle at low multiplicity", 200, 0, 2);
  TH1F* outputHigh = new TH1F("Track2ParticleHigh", "Track2Particle at high multiplicity", 200, 0, 2);

  TH3F* hist = fdNdEtaCorrectionLow[0]->GetTrack2ParticleCorrection()->GetCorrectionHistogram();
  TH3F* hist2 = fdNdEtaCorrectionHigh[0]->GetTrack2ParticleCorrection()->GetCorrectionHistogram();
  for (Int_t x=hist->GetXaxis()->FindBin(-10); x<=hist->GetXaxis()->FindBin(10); ++x)
    for (Int_t y=hist->GetYaxis()->FindBin(-0.8); y<=hist->GetYaxis()->FindBin(0.8); ++y)
      for (Int_t z=hist->GetZaxis()->FindBin(0.3); z<=hist->GetZaxis()->FindBin(9.9); ++z)
      //for (Int_t z=1; z<=hist->GetNbinsZ(); ++z)
      {
        if (hist->GetBinContent(x, y, z) > 0)
          outputLow->Fill(hist->GetBinContent(x, y, z));
        //if (hist->GetBinContent(x, y, z) == 1)
        //  printf("z = %f, eta = %f, pt = %f: %f %f %f\n", hist->GetXaxis()->GetBinCenter(x), hist->GetYaxis()->GetBinCenter(y), hist->GetZaxis()->GetBinCenter(z), hist->GetBinContent(x, y, z), fdNdEtaCorrectionLow[0]->GetTrack2ParticleCorrection()->GetGeneratedHistogram()->GetBinContent(x, y, z), fdNdEtaCorrectionLow[0]->GetTrack2ParticleCorrection()->GetMeasuredHistogram()->GetBinContent(x, y, z));

        if (hist2->GetBinContent(x, y, z) > 0)
          outputHigh->Fill(hist2->GetBinContent(x, y, z));
      }

  TCanvas* canvas = new TCanvas("EvaluateMultiplicityEffect", "EvaluateMultiplicityEffect", 1000, 500);
  canvas->Divide(2, 1);

  canvas->cd(1);
  outputLow->Draw();
  outputLow->Fit("gaus", "0");
  outputLow->GetFunction("gaus")->SetLineColor(2);
  outputLow->GetFunction("gaus")->DrawCopy("SAME");

  canvas->cd(2);
  outputHigh->Draw();
  outputHigh->Fit("gaus", "0");
  outputHigh->GetFunction("gaus")->DrawCopy("SAME");

  canvas->SaveAs(Form("%s.gif", canvas->GetName()));
}

void PlotErrors(Float_t xmin, Float_t xmax, Float_t ymin, Float_t ymax, Float_t zmin, Float_t zmax, const char* correctionMapFile = "correction_map.root", const char* correctionMapFolder = "dndeta_correction")
{
    gSystem->Load("libPWG0base");
	
	AlidNdEtaCorrection* dNdEtaCorrection = new AlidNdEtaCorrection(correctionMapFolder, correctionMapFolder);
    dNdEtaCorrection->LoadHistograms(correctionMapFile, correctionMapFolder);
    
    dNdEtaCorrection->GetTrack2ParticleCorrection()->PlotBinErrors(xmin, xmax, ymin, ymax, zmin, zmax)->Draw();
}
