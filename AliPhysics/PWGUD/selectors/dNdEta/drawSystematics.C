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

void loadlibs()
{
  gSystem->Load("libTree");
  gSystem->Load("libVMC");

  gSystem->Load("libSTEERBase");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libPWG0base");
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
    Correction1DCreatePlots(fileNames[i], folderNames[i], upperPtLimit);

    TH1* corrX = dynamic_cast<TH1*> (gROOT->FindObject("generated_x_div_measured_x"));
    TH1* corrY = dynamic_cast<TH1*> (gROOT->FindObject("generated_y_div_measured_y"));
    TH1* corrZ = dynamic_cast<TH1*> (gROOT->FindObject("generated_z_div_measured_z"));

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
    case 2: // + 30% kaons
    case 3: // - 30% kaons
    case 4: // + 30% protons
    case 5: // - 30% protons
    case 6: // + 30% kaons + 30% protons
    case 7: // - 30% kaons - 30% protons
    case 8: // + 30% kaons - 30% protons
    case 9: // - 30% kaons + 30% protons
    case 10: // + 30% others
    case 11: // - 30% others
      TString* str = new TString;
      if (index < 6)
      {
        Int_t correctionIndex = index / 2;
        Double_t scaleFactor = (index % 2 == 0) ? 1.3 : 0.7;
  
        fdNdEtaCorrection[correctionIndex]->GetTrack2ParticleCorrection()->GetTrackCorrection()->Scale(scaleFactor);
        str->Form("%s%s", (correctionIndex == 0) ? "Pi" : ((correctionIndex == 1) ? "K" : (correctionIndex == 2) ? "p" : "others"), (index % 2 == 0) ? "Boosted" : "Reduced");
      }
      else if (index < 10)
      {
        Double_t scaleFactor = (index % 2 == 0) ? 1.3 : 0.7;
        fdNdEtaCorrection[1]->GetTrack2ParticleCorrection()->GetTrackCorrection()->Scale(scaleFactor);
        str->Form("%s%s", "K", (scaleFactor > 1) ? "Boosted" : "Reduced");
        
        if (index >= 8)
          scaleFactor = (index % 2 == 0) ? 0.3 : 1.7;
        fdNdEtaCorrection[2]->GetTrack2ParticleCorrection()->GetTrackCorrection()->Scale(scaleFactor);
        *str += Form("%s%s", "p", (scaleFactor > 1) ? "Boosted" : "Reduced");
      }
      else
      {
        Double_t scaleFactor = (index % 2 == 0) ? 1.3 : 0.7;
        fdNdEtaCorrection[3]->GetTrack2ParticleCorrection()->GetTrackCorrection()->Scale(scaleFactor);
        str->Form("%s%s", "others", (scaleFactor > 1) ? "Boosted" : "Reduced");
      }

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
        ScalePtDependent(fdNdEtaCorrection[i]->GetTrack2ParticleCorrection()->GetTrackCorrection()->GetMeasuredHistogram(), ptDists[i]);
        ScalePtDependent(fdNdEtaCorrection[i]->GetTrack2ParticleCorrection()->GetTrackCorrection()->GetGeneratedHistogram(), ptDists[i]);
      }
      TString* str = new TString;
      str->Form("%s%s", (functionIndex == 0) ? "Pi" : ((functionIndex == 1) ? "K" : "p"), (index % 2 == 0) ? "Boosted" : "Reduced");
      return str->Data();
      break;

    case 999:
      TF1* ptDependence = new TF1("simple", "x", 0, 100);
      ScalePtDependent(fdNdEtaCorrection[0]->GetTrack2ParticleCorrection()->GetTrackCorrection()->GetGeneratedHistogram(), ptDependence);
      ScalePtDependent(fdNdEtaCorrection[0]->GetTrack2ParticleCorrection()->GetTrackCorrection()->GetMeasuredHistogram(), ptDependence);
      break;

  }

  return "noname";
}

void Composition()
{
  loadlibs();

  gSystem->Unlink("new_compositions.root");
  gSystem->Unlink("new_compositions_analysis.root");
  
  const char* names[] = { "pi", "K", "p", "other" };
  TH1* hRatios[20];

  //backgroundEvents = 1162+434; // Michele for MB1, run 104892, 15.02.10
  backgroundEvents = -1;    // use 0 bin from MC! for 2.36 TeV
  
  Printf("Subtracting %d background events!!!", backgroundEvents);
  gSystem->Sleep(1000);
  
  Int_t nCompositions = 12;
  Int_t counter = 0;
  for (Int_t comp = 1; comp < nCompositions; ++comp)
  {
    AlidNdEtaCorrection* fdNdEtaCorrection[4];

    TFile::Open("correction_mapparticle-species.root");

    for (Int_t i=0; i<4; ++i)
    {
      TString name;
      name.Form("dndeta_correction_%s", names[i]);
      fdNdEtaCorrection[i] = new AlidNdEtaCorrection(name, name);
      fdNdEtaCorrection[i]->LoadHistograms();
    }

    const char* newName = ChangeComposition(fdNdEtaCorrection, comp);

    Double_t geneCount[5];
    Double_t measCount[5];
    geneCount[4] = 0;
    measCount[4] = 0;

    for (Int_t i=0; i<4; ++i)
    {
      geneCount[i] = fdNdEtaCorrection[i]->GetTrack2ParticleCorrection()->GetTrackCorrection()->GetGeneratedHistogram()->Integral();
      measCount[i] = fdNdEtaCorrection[i]->GetTrack2ParticleCorrection()->GetTrackCorrection()->GetMeasuredHistogram()->Integral();

      geneCount[4] += geneCount[i];
      measCount[4] += measCount[i];

      printf("Particle %s: %d gene, %d meas\n", ((i == 0) ? "pi" : (i == 1) ? "K" : (i == 2) ? "p" : "others"), (Int_t) geneCount[i], (Int_t) measCount[i]);
    }

    printf("Generated ratios are:     %f pi, %f K, %f p, %f others\n", geneCount[0] / geneCount[4], geneCount[1] / geneCount[4], geneCount[2] / geneCount[4], geneCount[3] / geneCount[4]);

    printf("Reconstructed ratios are: %f pi, %f K, %f p, %f others\n", measCount[0] / measCount[4], measCount[1] / measCount[4], measCount[2] / measCount[4], measCount[3] / measCount[4]);

    TList* collection = new TList;

    // skip "other" particle correction here
    // with them has to be dealt differently, maybe just increasing the neutral particles...
    for (Int_t i=1; i<4; ++i)
      collection->Add(fdNdEtaCorrection[i]);

    fdNdEtaCorrection[0]->Merge(collection);
    fdNdEtaCorrection[0]->Finish();

    delete collection;

    // save everything
    TFile* file = TFile::Open("new_compositions.root", "UPDATE");
    fdNdEtaCorrection[0]->SetName(newName);
    fdNdEtaCorrection[0]->SaveHistograms();
    //file->Write();
    file->Close();
    
    // correct dNdeta distribution with modified correction map
    TFile::Open("analysis_esd_raw.root");

    dNdEtaAnalysis* fdNdEtaAnalysis = new dNdEtaAnalysis("fdNdEtaAnalysisESD", "fdNdEtaAnalysisESD");
    fdNdEtaAnalysis->LoadHistograms();

    fdNdEtaAnalysis->Finish(fdNdEtaCorrection[0], 0.2, 3, newName);
    
    hRatios[counter] = (TH1F*) fdNdEtaAnalysis->GetdNdEtaHistogram()->Clone(newName);
    hRatios[counter]->SetTitle(newName);
    hRatios[counter]->SetYTitle("dN_{ch}/d#eta ratio #frac{default composition}{modified composition}");

    if (counter > 0)
      hRatios[counter]->Divide(hRatios[0],hRatios[counter],1,1);

    file = TFile::Open("new_compositions_analysis.root", "UPDATE");
    hRatios[counter]->Write();
    file->Close();
    
    delete fdNdEtaAnalysis;

    counter++;
  }

  /*
  gROOT->ProcessLine(".L drawPlots.C");

  const char* fileNames[] = {"new_compositions.root", "new_compositions.root", "new_compositions.root", "new_compositions.root", "new_compositions.root", "new_compositions.root", "new_compositions.root", "new_compositions.root" };
  const char* folderNames[] = { "PythiaRatios", "PiBoosted", "PiReduced", "KBoosted", "KReduced", "pBoosted", "pReduced" };

  Track2Particle1DComposition(fileNames, nCompositions, folderNames);
  */
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

void mergeCorrectionsWithDifferentCrosssections(Int_t origin, Int_t correctionTarget = 3, Char_t* correctionFileName="correction_mapprocess-types.root", const char* analysisFileName = "analysis_esd_raw.root", const Char_t* outputFileName=0) {
  //
  // Function used to merge standard corrections with vertex
  // reconstruction corrections obtained by a certain mix of ND, DD
  // and SD events.
  //
  // the dn/deta spectrum is corrected and the ratios
  // (standard to changed x-section) of the different dN/deta
  // distributions are saved to a file.
  //
  // correctionTarget is of type AlidNdEtaCorrection::CorrectionType
  //    kINEL = 3
  //    kNSD = 4
  //    kOnePart = 6

  if (outputFileName == 0)
  {
    if (correctionTarget == 3)
      outputFileName = "systematics_vtxtrigger_compositions_inel.root";
    if (correctionTarget == 4)
      outputFileName = "systematics_vtxtrigger_compositions_nsd.root";
    if (correctionTarget == 6)
      outputFileName = "systematics_vtxtrigger_compositions_onepart.root";
  }

  loadlibs();

  const Char_t* typeName[] = { "vertexreco", "trigger", "vtxtrigger" };

  //Karel:
//     fsd = 0.153 +- 0.031 (0.050 to take into account SD definition) --> change
//     fdd = 0.080 +- 0.050 --> change 
//     fnd = 0.767 +- 0.059 --> keep (error small)

//  const Char_t* changes[]  = { "pythia","ddmore","ddless","sdmore","sdless", "dmore", "dless", "sdlessddmore", "sdmoreddless", "ddmore25","ddless25","sdmore25","sdless25", "dmore25", "dless25", "sdlessddmore25", "sdmoreddless25"};
  //Float_t scalesND[] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0 };
  //Float_t scalesDD[] = {1.0, 1.5, 0.5, 1.0, 1.0, 1.5, 0.5, 1.5, 0.5, 1.25, 0.75, 1.0,  1.0,  1.25, 0.75, 1.25, 0.75};
  //Float_t scalesSD[] = {1.0, 1.0, 1.0, 1.5, 0.5, 1.5, 0.5, 0.5, 1.5, 1.0,  1.0,  1.25, 0.75, 1.25, 0.75, 0.75, 1.25};
  //Float_t scalesDD[] = {1.0, 1.4, 0.6, 1.0, 1.0, 1.4, 0.6, 1.4, 0.6, 1.25, 0.75, 1.0,  1.0,  1.25, 0.75, 1.25, 0.75};
  //Float_t scalesSD[] = {1.0, 1.0, 1.0, 1.4, 0.6, 1.4, 0.6, 0.4, 1.6, 1.0,  1.0,  1.25, 0.75, 1.25, 0.75, 0.75, 1.25};
/*  Int_t nChanges = 9;

  const Char_t* changes[]  = { "ua5","ddmore","ddless","sdmore","sdless", "dmore", "dless", "sdlessddmore", "sdmoreddless" };
  Float_t scalesND[] = {0.767, 0.767, 0.767, 0.767, 0.767, 0.767, 0.767, 0.767, 0.767};
  Float_t scalesDD[] = {0.080, 0.130, 0.030, 0.080, 0.080, 0.130, 0.030, 0.130, 0.030};
  Float_t scalesSD[] = {0.153, 0.153, 0.153, 0.203, 0.103, 0.203, 0.103, 0.103, 0.203};*/
  
  Float_t ref_SD = -1;
  Float_t ref_DD = -1;
  Float_t ref_ND = -1;
  
  Float_t error_SD = -1;
  Float_t error_DD = -1;
  Float_t error_ND = -1;
  
  GetRelativeFractions(origin, ref_SD, ref_DD, ref_ND, error_SD, error_DD, error_ND);
  
  Printf("Using x-sections:\n SD: %f +- %f\n DD: %f +- %f\n ND: %f +- %f", ref_SD, error_SD, ref_DD, error_DD, ref_ND, error_ND);
  
  const Char_t* changes[]  = { "default","sdless","sdmore","ddless","ddmore", "dless", "dmore", "sdlessddmore", "sdmoreddless" };
  Int_t nChanges = 9;
  Float_t scalesSD[9];
  Float_t scalesDD[9];
  Float_t scalesND[9];
  
  if (1)
  {
    // sample 8 points on the error ellipse
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
      
      scalesSD[i] = ref_SD + factorSD * error_SD;
      scalesDD[i] = ref_DD + factorDD * error_DD;
      scalesND[i] = 1.0 - scalesDD[i] - scalesSD[i];
      
      Printf("Case %d: SD: %f DD: %f ND: %f", i, scalesSD[i], scalesDD[i], scalesND[i]);
    }
  }
  else
  {
    Printf("WARNING: Special treatment for ratios active");
    gSystem->Sleep(1000);
    
    // constrained values by allowed changing of cross-sections
    Float_t pythiaScaling = 0.224 / 0.189;

    if (origin == 10)
    {
      // 900 GeV
      for (Int_t i=0; i<9; i++)
      {
        scalesSD[i] = 15.3;
        scalesDD[i] = 9.5;
      }

      scalesSD[1] = 15.7;
      scalesSD[2] = 17.6;
      scalesSD[3] = 13.5;
      scalesSD[4] = 17.6;

      scalesDD[5] = 15.5;
      scalesDD[6] = 8.8;
      scalesDD[7] = 13.8;
      scalesDD[8] = 7.6;
    }
    else if (origin == 20)
    {
      // 2.36 TeV
      pythiaScaling = 0.217 / 0.167;
      
      for (Int_t i=0; i<9; i++)
      {
        scalesSD[i] = 15.9;
        scalesDD[i] = 10.7;
      }

      scalesSD[1] = 13.5;
      scalesSD[2] = 15.2;
      scalesSD[3] = 13.5;
      scalesSD[4] = 17.6;

      scalesDD[5] = 13.8;
      scalesDD[6] = 7.6;
      scalesDD[7] = 13.8;
      scalesDD[8] = 7.6;
    }
    else
      AliFatal("Not supported");

    for (Int_t i=0; i<9; i++)
    {
      scalesSD[i] /= 100;
      scalesSD[i] *= pythiaScaling;
      scalesDD[i] /= 100;
      scalesND[i] = 1.0 - scalesDD[i] - scalesSD[i];
      Printf("Case %d: SD: %f DD: %f ND: %f", i, scalesSD[i], scalesDD[i], scalesND[i]);
    }
  }
  
  Int_t backgroundEvents = 0;
  
  //backgroundEvents = 1162+434; // Michele for MB1, run 104892, 15.02.10
  //backgroundEvents = 6;          // Michele for V0AND, run 104892, 15.02.10
  
  //backgroundEvents = 4398+961;   // Michele for MB1, run 104824-52, 16.02.10
  //backgroundEvents = 19;         // Michele for V0AND, run 104824-52, 16.02.10
  
  backgroundEvents = -1;    // use 0 bin from MC! for 2.36 TeV
  
  Printf("Subtracting %d background events!!!", backgroundEvents);
  gSystem->Sleep(1000);
  
  /*
  const Char_t* changes[]  = { "pythia", "qgsm", "phojet"};
  Float_t scalesND[] = {1.0, 1.10, 1.11};
  Float_t scalesSD[] = {1.0, 0.69, 0.86};
  Float_t scalesDD[] = {1.0, 0.98, 0.61};
  Int_t nChanges = 3;
  */
  
  // cross section from Pythia
  // 14 TeV!
//   Float_t sigmaND = 55.2;
//   Float_t sigmaDD = 9.78;
//   Float_t sigmaSD = 14.30;

  // standard correction
  TFile::Open(correctionFileName);
  AlidNdEtaCorrection* correctionStandard = new AlidNdEtaCorrection("dndeta_correction","dndeta_correction");
  correctionStandard->LoadHistograms();

  // dont take vertexreco from this one
  correctionStandard->GetVertexRecoCorrection()->Reset();
  // dont take triggerbias from this one
  correctionStandard->GetTriggerBiasCorrectionINEL()->Reset();
  correctionStandard->GetTriggerBiasCorrectionNSD()->Reset();
  correctionStandard->GetTriggerBiasCorrectionND()->Reset();
  correctionStandard->GetTriggerBiasCorrectionOnePart()->Reset();

  AlidNdEtaCorrection* corrections[100];
  TH1F* hRatios[100];

  Int_t counter = 0;
  for (Int_t j=2; j<3; j++) { // j = 0 (change vtx), j = 1 (change trg), j = 2 (change both)

    for (Int_t i=0; i<nChanges; i++) {
      TFile::Open(correctionFileName);

      TString name;
      name.Form("dndeta_correction_syst_%s_%s", typeName[j], changes[i]);
      AlidNdEtaCorrection* current = new AlidNdEtaCorrection(name, name);
      current->LoadHistograms("dndeta_correction");
      current->Reset();

      name.Form("dndeta_correction_ND");
      AlidNdEtaCorrection* dNdEtaCorrectionND = new AlidNdEtaCorrection(name,name);
      dNdEtaCorrectionND->LoadHistograms();
      name.Form("dndeta_correction_DD");
      AlidNdEtaCorrection* dNdEtaCorrectionDD = new AlidNdEtaCorrection(name,name);
      dNdEtaCorrectionDD->LoadHistograms();
      name.Form("dndeta_correction_SD");
      AlidNdEtaCorrection* dNdEtaCorrectionSD = new AlidNdEtaCorrection(name,name);
      dNdEtaCorrectionSD->LoadHistograms();

      // calculating relative
      Float_t nd = dNdEtaCorrectionND->GetTriggerBiasCorrectionINEL()->GetEventCorrection()->GetGeneratedHistogram()->Integral();
      Float_t dd = dNdEtaCorrectionDD->GetTriggerBiasCorrectionINEL()->GetEventCorrection()->GetGeneratedHistogram()->Integral();
      Float_t sd = dNdEtaCorrectionSD->GetTriggerBiasCorrectionINEL()->GetEventCorrection()->GetGeneratedHistogram()->Integral();
      Float_t total = nd + dd + sd;
      
      nd /= total;
      sd /= total;
      dd /= total;
      
      Printf("Ratios in the correction map are: ND=%f, DD=%f, SD=%f", nd, dd, sd);
      
      Float_t scaleND = scalesND[i] / nd;
      Float_t scaleDD = scalesDD[i] / dd;
      Float_t scaleSD = scalesSD[i] / sd;
      
      Printf("ND=%.2f, DD=%.2f, SD=%.2f",scaleND, scaleDD, scaleSD);      
      
/*      Float_t nd = 100 * sigmaND/(sigmaND + (scalesDD[i]*sigmaDD) + (scalesDD[i]*sigmaSD));
      Float_t dd = 100 * (scalesDD[i]*sigmaDD)/(sigmaND + (scalesDD[i]*sigmaDD) + (scalesDD[i]*sigmaSD));
      Float_t sd = 100 * (scalesSD[i]*sigmaSD)/(sigmaND + (scalesDD[i]*sigmaDD) + (scalesDD[i]*sigmaSD));

      printf(Form("%s : ND=%.2f\%, DD=%.2f\%, SD=%.2f\% \n",changes[i],nd,dd,sd));*/
      current->SetTitle(Form("ND=%.2f\%,DD=%.2f\%,SD=%.2f\%",scaleND,scaleDD,scaleSD));
      current->SetTitle(name);

      // scale
      if (j == 0 || j == 2)
      {
        dNdEtaCorrectionND->GetVertexRecoCorrection()->Scale(scaleND);
        dNdEtaCorrectionDD->GetVertexRecoCorrection()->Scale(scaleDD);
        dNdEtaCorrectionSD->GetVertexRecoCorrection()->Scale(scaleSD);
      }
      if (j == 1 || j == 2)
      {
        dNdEtaCorrectionND->GetTriggerBiasCorrectionINEL()->Scale(scaleND);
        dNdEtaCorrectionDD->GetTriggerBiasCorrectionINEL()->Scale(scaleDD);
        dNdEtaCorrectionSD->GetTriggerBiasCorrectionINEL()->Scale(scaleSD);

        dNdEtaCorrectionND->GetTriggerBiasCorrectionNSD()->Scale(scaleND);
        dNdEtaCorrectionDD->GetTriggerBiasCorrectionNSD()->Scale(scaleDD);
        dNdEtaCorrectionSD->GetTriggerBiasCorrectionNSD()->Scale(scaleSD);

        dNdEtaCorrectionND->GetTriggerBiasCorrectionND()->Scale(scaleND);
        dNdEtaCorrectionDD->GetTriggerBiasCorrectionND()->Scale(scaleDD);
        dNdEtaCorrectionSD->GetTriggerBiasCorrectionND()->Scale(scaleSD);

        dNdEtaCorrectionND->GetTriggerBiasCorrectionOnePart()->Scale(scaleND);
        dNdEtaCorrectionDD->GetTriggerBiasCorrectionOnePart()->Scale(scaleDD);
        dNdEtaCorrectionSD->GetTriggerBiasCorrectionOnePart()->Scale(scaleSD);
      }

      //clear track in correction
      dNdEtaCorrectionND->GetTrack2ParticleCorrection()->Reset();
      dNdEtaCorrectionDD->GetTrack2ParticleCorrection()->Reset();
      dNdEtaCorrectionSD->GetTrack2ParticleCorrection()->Reset();

      TList collection;
      collection.Add(correctionStandard);
      collection.Add(dNdEtaCorrectionND);
      collection.Add(dNdEtaCorrectionDD);
      collection.Add(dNdEtaCorrectionSD);

      current->Merge(&collection);
      current->Finish();
      
      // print 0 bin efficiency
      TH1* hist2 = current->GetTriggerBiasCorrectionINEL()->GetEventCorrection()->Get1DCorrection("y", -10, 10);
      if (hist2->GetBinContent(1))
      {
        Float_t triggerEff = 100.0 / hist2->GetBinContent(1);
        Printf("trigger eff in 0 bin is: %.2f %%", triggerEff);
      }

      corrections[counter] = current;

      // now correct dNdeta distribution with modified correction map
      TFile::Open(analysisFileName);

      dNdEtaAnalysis* fdNdEtaAnalysis = new dNdEtaAnalysis("fdNdEtaAnalysisESD", "fdNdEtaAnalysisESD");
      fdNdEtaAnalysis->LoadHistograms();

      fdNdEtaAnalysis->Finish(current, 0.151, correctionTarget, Form("%d %d", j, i), backgroundEvents);

      name = "ratio";
      if (j==0) name.Append("_vetexReco_");
      if (j==1) name.Append("_triggerBias_");
      if (j==2) name.Append("_vertexReco_triggerBias_");
      name.Append(changes[i]);

      hRatios[counter] = (TH1F*) fdNdEtaAnalysis->GetdNdEtaHistogram()->Clone(name);

      name.Form("ND #times %0.2f DD #times %0.2f, SD #times %0.2f",scaleND,scaleDD,scaleSD);
      hRatios[counter]->SetTitle(name.Data());
      hRatios[counter]->SetYTitle("dN_{ch}/d#eta ratio #frac{default cross-section}{modified cross-sections}");
      
      TF1* pol0 = new TF1("pol0", "[0]", -0.49, 0.49);
      pol0->SetParameter(0, 0);
      hRatios[counter]->Fit(pol0, "RN");
      Printf("Case %d: %f", i, pol0->GetParameter(0));
      
      if (counter > 0)
        hRatios[counter]->Divide(hRatios[0],hRatios[counter],1,1);

      delete fdNdEtaAnalysis;

      counter++;
    }
  }

  TFile* fout = new TFile(outputFileName,"RECREATE");

  // to make everything consistent
  hRatios[0]->Divide(hRatios[0],hRatios[0],1,1);

  for (Int_t i=0; i<counter; i++)
  {
    corrections[i]->SaveHistograms();
    hRatios[i]->Write();
  }

  fout->Write();
  fout->Close();
}

void GetRelativeFractions(Int_t origin, Float_t& ref_SD, Float_t& ref_DD, Float_t& ref_ND, Float_t& error_SD, Float_t& error_DD, Float_t& error_ND)
{
  // origin: 
  //   -1 = Pythia (test)
  //   0 = UA5
  //   1 = Data 1.8 TeV
  //   2 = Tel-Aviv
  //   3 = Durham
  //

  switch (origin)
  {
    case -10: // Pythia default at 7 GeV, 50% error
      Printf("PYTHIA x-sections");
      ref_SD = 0.192637; error_SD = ref_SD * 0.5;
      ref_DD = 0.129877; error_DD = ref_DD * 0.5;
      ref_ND = 0.677486; error_ND = 0;
      break;

    case -1: // Pythia default at 900 GeV, as test
      Printf("PYTHIA x-sections");
      ref_SD = 0.223788;
      ref_DD = 0.123315;
      ref_ND = 0.652897;
      break;
      
    case 0: // UA5
      Printf("UA5 x-sections a la first paper");
      ref_SD = 0.153; error_SD = 0.05;
      ref_DD = 0.080; error_DD = 0.05;
      ref_ND = 0.767; error_ND = 0;
      break;
      
    case 10: // UA5
      Printf("UA5 x-sections hadron level definition for Pythia"); 
      // Fractions in Pythia with UA5 cuts selection for SD
      // ND: 0.688662
      // SD: 0.188588 --> this should be 15.3
      // DD: 0.122750
      ref_SD = 0.224 * 0.153 / 0.189; error_SD = 0.023 * 0.224 / 0.189;
      ref_DD = 0.095;                 error_DD = 0.06; 
      ref_ND = 1.0 - ref_SD - ref_DD; error_ND = 0;
      break;
    
    case 11: // UA5
      Printf("UA5 x-sections hadron level definition for Phojet"); 
      // Fractions in Phojet with UA5 cuts selection for SD
      // ND: 0.783573
      // SD: 0.151601 --> this should be 15.3
      // DD: 0.064827
      ref_SD = 0.191 * 0.153 / 0.152; error_SD = 0.023 * 0.191 / 0.152;
      ref_DD = 0.095;                 error_DD = 0.06; 
      ref_ND = 1.0 - ref_SD - ref_DD; error_ND = 0;
      break;
      
    case 20: // E710, 1.8 TeV
      Printf("E710 x-sections hadron level definition for Pythia");
      // ND: 0.705709
      // SD: 0.166590 --> this should be 15.9
      // DD: 0.127701
      ref_SD = 0.217 * 0.159 / 0.167; error_SD = 0.024 * 0.217 / 0.167;
      ref_DD = 0.075 * 1.43;          error_DD = 0.02 * 1.43; 
      ref_ND = 1.0 - ref_SD - ref_DD; error_ND = 0;
      break;
    
    case 21: // E710, 1.8 TeV
      Printf("E710 x-sections hadron level definition for Phojet"); 
      // ND: 0.817462
      // SD: 0.125506 --> this should be 15.9
      // DD: 0.057032
      ref_SD = 0.161 * 0.159 / 0.126; error_SD = 0.024 * 0.161 / 0.126;
      ref_DD = 0.075 * 1.43;         error_DD = 0.02 * 1.43;
      ref_ND = 1.0 - ref_SD - ref_DD; error_ND = 0;
      break;
    
    case 1: // data 1.8 TeV
      Printf("??? x-sections");
      ref_SD = 0.152;
      ref_DD = 0.092;
      ref_ND = 1 - ref_SD - ref_DD;
      break;
      
    case 2: // tel-aviv model
      Printf("Tel-aviv model x-sections");
      ref_SD = 0.171;
      ref_DD = 0.094;
      ref_ND = 1 - ref_SD - ref_DD;
      break;
    
    case 3: // durham model
      Printf("Durham model x-sections");
      ref_SD = 0.190;
      ref_DD = 0.125;
      ref_ND = 1 - ref_SD - ref_DD;
      break;
    
    default:
      AliFatal(Form("Unknown origin %d", origin));
  }
}

void CreateCorrectionsWithUA5CrossSections(Int_t origin, const Char_t* correctionFileName="correction_mapprocess-types.root", const Char_t* outputFileName="correction_map2.root") {
  //
  // Function used to merge standard corrections with vertex
  // reconstruction corrections obtained by a certain mix of ND, DD
  // and SD events.
  //
  loadlibs();

  const Char_t* typeName[] = { "vertexreco", "trigger", "vtxtrigger" };
  
  Float_t ref_SD = -1;
  Float_t ref_DD = -1;
  Float_t ref_ND = -1;
  
  Float_t error_SD = -1;
  Float_t error_DD = -1;
  Float_t error_ND = -1;
  
  GetRelativeFractions(origin, ref_SD, ref_DD, ref_ND, error_SD, error_DD, error_ND);
  
  Printf("Using x-sections:\n SD: %f +- %f\n DD: %f +- %f\n ND: %f +- %f", ref_SD, error_SD, ref_DD, error_DD, ref_ND, error_ND);
  
//Karel (UA5):
//     fsd = 0.153 +- 0.031
//     fdd = 0.080 +- 0.050
//     fnd = 0.767 +- 0.059

//       Karel (1.8 TeV):
//       
//       Tel-Aviv model Sd/Inel = 0.171           Dd/Inel = 0.094
//       Durham model   Sd/Inel = 0.190           Dd/Inel = 0.125
//       Data           Sd/Inel = 0.152 +- 0.030  Dd/Inel = 0.092 +- 0.45

  // standard correction
  TFile::Open(correctionFileName);
  AlidNdEtaCorrection* correctionStandard = new AlidNdEtaCorrection("dndeta_correction","dndeta_correction");
  correctionStandard->LoadHistograms();

  // dont take vertexreco from this one
  correctionStandard->GetVertexRecoCorrection()->Reset();
  // dont take triggerbias from this one
  correctionStandard->GetTriggerBiasCorrectionINEL()->Reset();
  correctionStandard->GetTriggerBiasCorrectionNSD()->Reset();
  correctionStandard->GetTriggerBiasCorrectionND()->Reset();
  correctionStandard->GetTriggerBiasCorrectionOnePart()->Reset();

  AlidNdEtaCorrection* corrections[100];
  TH1F* hRatios[100];

  Int_t counter = 0;
      
  TFile::Open(correctionFileName);

  AlidNdEtaCorrection* current = new AlidNdEtaCorrection("dndeta_correction_ua5", "dndeta_correction_ua5");
  current->LoadHistograms("dndeta_correction");
  current->Reset();

  TString name;
  name.Form("dndeta_correction_ND");
  AlidNdEtaCorrection* dNdEtaCorrectionND = new AlidNdEtaCorrection(name,name);
  dNdEtaCorrectionND->LoadHistograms();
  name.Form("dndeta_correction_DD");
  AlidNdEtaCorrection* dNdEtaCorrectionDD = new AlidNdEtaCorrection(name,name);
  dNdEtaCorrectionDD->LoadHistograms();
  name.Form("dndeta_correction_SD");
  AlidNdEtaCorrection* dNdEtaCorrectionSD = new AlidNdEtaCorrection(name,name);
  dNdEtaCorrectionSD->LoadHistograms();

  // calculating relative
  Float_t nd = dNdEtaCorrectionND->GetTriggerBiasCorrectionINEL()->GetEventCorrection()->GetGeneratedHistogram()->Integral();
  Float_t dd = dNdEtaCorrectionDD->GetTriggerBiasCorrectionINEL()->GetEventCorrection()->GetGeneratedHistogram()->Integral();
  Float_t sd = dNdEtaCorrectionSD->GetTriggerBiasCorrectionINEL()->GetEventCorrection()->GetGeneratedHistogram()->Integral();
  Float_t total = nd + dd + sd;
  
  nd /= total;
  sd /= total;
  dd /= total;
  
  Printf("Ratios in the correction map are: ND=%f, DD=%f, SD=%f", nd, dd, sd);
  
  Float_t scaleND = ref_ND / nd;
  Float_t scaleDD = ref_DD / dd;
  Float_t scaleSD = ref_SD / sd;
  
  Printf("ND=%.2f, DD=%.2f, SD=%.2f",scaleND, scaleDD, scaleSD);

  // scale
  dNdEtaCorrectionND->GetVertexRecoCorrection()->Scale(scaleND);
  dNdEtaCorrectionDD->GetVertexRecoCorrection()->Scale(scaleDD);
  dNdEtaCorrectionSD->GetVertexRecoCorrection()->Scale(scaleSD);
    
  dNdEtaCorrectionND->GetTriggerBiasCorrectionINEL()->Scale(scaleND);
  dNdEtaCorrectionDD->GetTriggerBiasCorrectionINEL()->Scale(scaleDD);
  dNdEtaCorrectionSD->GetTriggerBiasCorrectionINEL()->Scale(scaleSD);

  dNdEtaCorrectionND->GetTriggerBiasCorrectionNSD()->Scale(scaleND);
  dNdEtaCorrectionDD->GetTriggerBiasCorrectionNSD()->Scale(scaleDD);
  dNdEtaCorrectionSD->GetTriggerBiasCorrectionNSD()->Scale(scaleSD);

  dNdEtaCorrectionND->GetTriggerBiasCorrectionND()->Scale(scaleND);
  dNdEtaCorrectionDD->GetTriggerBiasCorrectionND()->Scale(scaleDD);
  dNdEtaCorrectionSD->GetTriggerBiasCorrectionND()->Scale(scaleSD);

  dNdEtaCorrectionND->GetTriggerBiasCorrectionOnePart()->Scale(scaleND);
  dNdEtaCorrectionDD->GetTriggerBiasCorrectionOnePart()->Scale(scaleDD);
  dNdEtaCorrectionSD->GetTriggerBiasCorrectionOnePart()->Scale(scaleSD);

  //clear track in correction
  dNdEtaCorrectionND->GetTrack2ParticleCorrection()->Reset();
  dNdEtaCorrectionDD->GetTrack2ParticleCorrection()->Reset();
  dNdEtaCorrectionSD->GetTrack2ParticleCorrection()->Reset();

  TList collection;
  collection.Add(correctionStandard);
  collection.Add(dNdEtaCorrectionND);
  collection.Add(dNdEtaCorrectionDD);
  collection.Add(dNdEtaCorrectionSD);

  current->Merge(&collection);
  current->Finish();

  // print 0 bin efficiency
  TH1* hist2 = current->GetTriggerBiasCorrectionINEL()->GetEventCorrection()->Get1DCorrection("y", -10, 10);
  if (hist2->GetBinContent(1) > 0)
  {
    Float_t triggerEff = 100.0 / hist2->GetBinContent(1);
    Printf("trigger eff in 0 bin is: %.2f %%", triggerEff);
  }
  
  TFile* fout = new TFile(outputFileName,"RECREATE");
  current->SaveHistograms();

  fout->Write();
  fout->Close();

  Printf("Trigger efficiencies:");
  Printf("ND: %.2f %%", 100.0 * dNdEtaCorrectionND->GetTriggerBiasCorrectionINEL()->GetEventCorrection()->GetMeasuredHistogram()->Integral() / dNdEtaCorrectionND->GetTriggerBiasCorrectionINEL()->GetEventCorrection()->GetGeneratedHistogram()->Integral());
  Printf("SD: %.2f %%", 100.0 * dNdEtaCorrectionSD->GetTriggerBiasCorrectionINEL()->GetEventCorrection()->GetMeasuredHistogram()->Integral() / dNdEtaCorrectionSD->GetTriggerBiasCorrectionINEL()->GetEventCorrection()->GetGeneratedHistogram()->Integral());
  Printf("DD: %.2f %%", 100.0 * dNdEtaCorrectionDD->GetTriggerBiasCorrectionINEL()->GetEventCorrection()->GetMeasuredHistogram()->Integral() / dNdEtaCorrectionDD->GetTriggerBiasCorrectionINEL()->GetEventCorrection()->GetGeneratedHistogram()->Integral());
  Printf("INEL: %.2f %%", 100.0 * current->GetTriggerBiasCorrectionINEL()->GetEventCorrection()->GetMeasuredHistogram()->Integral() / current->GetTriggerBiasCorrectionINEL()->GetEventCorrection()->GetGeneratedHistogram()->Integral());
  Printf("NSD: %.2f %%", 100.0 * (dNdEtaCorrectionND->GetTriggerBiasCorrectionINEL()->GetEventCorrection()->GetMeasuredHistogram()->Integral() + dNdEtaCorrectionDD->GetTriggerBiasCorrectionINEL()->GetEventCorrection()->GetMeasuredHistogram()->Integral()) / (dNdEtaCorrectionND->GetTriggerBiasCorrectionINEL()->GetEventCorrection()->GetGeneratedHistogram()->Integral() + dNdEtaCorrectionDD->GetTriggerBiasCorrectionINEL()->GetEventCorrection()->GetGeneratedHistogram()->Integral()));
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

void changePtSpectrum(const char* fileName = "analysis_mc.root", Float_t ptCutOff = 0.2, const char* fileName2 = 0)
{
  Float_t factor = 0.5;

  TFile* file = TFile::Open(fileName);
  TH1F* hist = dynamic_cast<TH1F*> (file->Get("dndeta_check_pt")->Clone());
  
  TH1* hist2 = 0;
  if (fileName2)
  {
    file2 = TFile::Open(fileName2);
    hist2 = dynamic_cast<TH1*> (file2->Get("dndeta_check_pt")->Clone());
    hist2->Scale(hist->GetBinContent(hist->FindBin(ptCutOff)) / hist2->GetBinContent(hist2->FindBin(ptCutOff)));
  }
  
  //hist->Scale(1.0 / hist->Integral());

  //hist->Rebin(3);
  //hist->Scale(1.0/3);

  TH1F* clone1 = dynamic_cast<TH1F*> (hist->Clone("clone1"));
  TH1F* clone2 = dynamic_cast<TH1F*> (hist->Clone("clone2"));

  TH1F* scale1 =  dynamic_cast<TH1F*> (hist->Clone("scale1"));
  TH1F* scale2 =  dynamic_cast<TH1F*> (hist->Clone("scale2"));

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
      scale1->SetBinContent(i, 1 - (ptCutOff - hist->GetBinCenter(i)) / ptCutOff * factor);

      // 110% at pt = 0, ...
      scale2->SetBinContent(i, 1 + (ptCutOff - hist->GetBinCenter(i)) / ptCutOff * factor);
    }
    scale1->SetBinError(i, 0);
    scale2->SetBinError(i, 0);
  }

  /*
  new TCanvas;
  scale1->Draw();
  scale2->SetLineColor(kRed);
  scale2->Draw("SAME");
  */

  clone1->Multiply(scale1);
  clone2->Multiply(scale2);

  Prepare1DPlot(hist);
  Prepare1DPlot(clone1);
  Prepare1DPlot(clone2);

  /*hist->SetMarkerStyle(markers[0]);
  clone1->SetMarkerStyle(markers[0]);
  clone2->SetMarkerStyle(markers[0]);*/

  hist->SetTitle(";p_{T} in GeV/c;dN_{ch}/dp_{T} in c/GeV");
  hist->GetXaxis()->SetRangeUser(0, 0.5);
  hist->GetYaxis()->SetRangeUser(0.01, clone2->GetMaximum() * 1.1);
  hist->GetYaxis()->SetTitleOffset(1);

  TCanvas* canvas = new TCanvas("c", "c", 600, 600);
  gPad->SetGridx();
  gPad->SetGridy();
  gPad->SetTopMargin(0.05);
  gPad->SetRightMargin(0.05);
  gPad->SetBottomMargin(0.12);
  hist->Draw("H");
  clone1->SetLineColor(kRed);
  clone1->Draw("HSAME");
  clone2->SetLineColor(kBlue);
  clone2->Draw("HSAME");
  hist->Draw("HSAME");
  
  if (hist2)
  {
    Prepare1DPlot(hist2);
    hist2->SetLineStyle(2);
    hist2->Draw("HSAME");
  }

  Float_t fraction =  hist->Integral(hist->GetXaxis()->FindBin(ptCutOff), hist->GetNbinsX()) / hist->Integral(1, hist->GetNbinsX());
  Float_t fraction1 = clone1->Integral(clone1->GetXaxis()->FindBin(ptCutOff), clone1->GetNbinsX()) / clone1->Integral(1, clone1->GetNbinsX());
  Float_t fraction2 = clone2->Integral(clone2->GetXaxis()->FindBin(ptCutOff), clone2->GetNbinsX()) / clone2->Integral(1, clone2->GetNbinsX());

  printf("%f %f %f\n", fraction, fraction1, fraction2);
  printf("Rel. %f %f\n", fraction1 / fraction, fraction2 / fraction);

  //canvas->SaveAs("changePtSpectrum.gif");
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


void runStudy(const char* baseCorrectionMapFile = "correction_map.root", const char* baseCorrectionMapFolder = "dndeta_correction", const char* changedCorrectionMapFile = "correction_map.root", const char* changedCorrectionMapFolder = "dndeta_correction", const char* dataFile = "analysis_esd_raw.root", const char* output = "analysis_esd_syst.root")
{
  gSystem->Load("libPWG0base");

  TFile* file = TFile::Open(output, "RECREATE");

  const Int_t max = 5;
  dNdEtaAnalysis* fdNdEtaAnalysis[5];

  new TCanvas;
  TLegend* legend = new TLegend(0.63, 0.73, 0.98, 0.98);
  legend->SetFillColor(0);

  for (Int_t i = 0; i < max; ++i)
  {
    TFile::Open(baseCorrectionMapFile);
    AlidNdEtaCorrection* baseCorrection = new AlidNdEtaCorrection(baseCorrectionMapFolder, baseCorrectionMapFolder);
    baseCorrection->LoadHistograms();

    AlidNdEtaCorrection::CorrectionType correctionType = AlidNdEtaCorrection::kNone;
    const char* name = 0;

    TFile::Open(changedCorrectionMapFile);
    switch (i)
    {
      case 0 :
        name = "default";
        break;

      case 1 :
        baseCorrection->GetTrack2ParticleCorrection()->LoadHistograms(Form("%s/Track2Particle", changedCorrectionMapFolder));
        name = "Track2Particle";
        break;

      case 2 :
        baseCorrection->GetVertexRecoCorrection()->LoadHistograms(Form("%s/VertexReconstruction", changedCorrectionMapFolder));
        name = "VertexReco";
        break;

      case 3 :
        baseCorrection->GetTriggerBiasCorrectionINEL()->LoadHistograms(Form("%s/TriggerBias_MBToINEL", changedCorrectionMapFolder));
        name = "TriggerBias_MBToINEL";
        break;

      case 4 :
        baseCorrection->LoadHistograms(changedCorrectionMapFolder);
        name = "all";
        break;

      default: return;
    }

    TFile::Open(dataFile);
    fdNdEtaAnalysis[i] = new dNdEtaAnalysis(name, name);
    fdNdEtaAnalysis[i]->LoadHistograms("dndeta");

    fdNdEtaAnalysis[i]->Finish(baseCorrection, 0.3, AlidNdEtaCorrection::kINEL);
    file->cd();
    fdNdEtaAnalysis[i]->SaveHistograms();

    TH1* hist = fdNdEtaAnalysis[i]->GetdNdEtaHistogram(0);
    hist->SetLineColor(colors[i]);
    hist->SetMarkerColor(colors[i]);
    hist->SetMarkerStyle(markers[i]+1);
    hist->DrawCopy((i == 0) ? "" : "SAME");
    legend->AddEntry(hist, name);
  }

  legend->Draw();
}

void ChangePtInCorrection(const char* fileName = "correction_map.root", const char* dirName = "dndeta_correction")
{
  loadlibs();
  if (!TFile::Open(fileName))
    return;

  AlidNdEtaCorrection* dNdEtaCorrection = new AlidNdEtaCorrection(dirName, dirName);
  if (!dNdEtaCorrection->LoadHistograms())
    return;

  dNdEtaCorrection->Finish();

  AliCorrection* corr = dNdEtaCorrection->GetCorrection(AlidNdEtaCorrection::kTrack2Particle);
 
  Printf(">>>> Before");
  corr->PrintInfo(0);

  Float_t factor = 0.5;
  Float_t ptCutOff = 0.2;
  
  TH3* gene = corr->GetTrackCorrection()->GetGeneratedHistogram();
  TH3* meas = corr->GetTrackCorrection()->GetMeasuredHistogram();
  
  for (Int_t z = 1; z <= gene->GetZaxis()->FindBin(ptCutOff - 0.01); z++)
  {
    Float_t localFactor = 1 - (ptCutOff - gene->GetZaxis()->GetBinCenter(z)) / ptCutOff * factor;
    Printf("%f %f", gene->GetZaxis()->GetBinCenter(z), localFactor);
    for (Int_t x = 1; x <= gene->GetNbinsX(); x++)
      for (Int_t y = 1; y <= gene->GetNbinsY(); y++)
      {
        gene->SetBinContent(x, y, z, gene->GetBinContent(x, y, z) * localFactor);
        meas->SetBinContent(x, y, z, meas->GetBinContent(x, y, z) * localFactor);
      }
  }
  
  dNdEtaCorrection->Finish();
  
  Printf(">>>> After");
  corr->PrintInfo(0);
}

Float_t FitAverage(TH1* hist, Int_t n, Float_t* begin, Float_t *end, Int_t color, Int_t& totalBins)
{
  Float_t average = 0;
  totalBins = 0;
  
  for (Int_t i=0; i<n; i++)
  {
    func = new TF1("func", "[0]", hist->GetXaxis()->GetBinLowEdge(hist->GetXaxis()->FindBin(begin[i])), hist->GetXaxis()->GetBinUpEdge(hist->GetXaxis()->FindBin(end[i])));
    Int_t bins = hist->GetXaxis()->FindBin(end[i]) - hist->GetXaxis()->FindBin(begin[i]) + 1;
    func->SetParameter(0, 1);
    func->SetLineColor(color);

    hist->Fit(func, "RNQ");
    func->Draw("SAME");
    
    average += func->GetParameter(0) * bins;
    totalBins += bins;
  }
  
  return average / totalBins;
}

void SPDIntegrateGaps(Bool_t all, const char* mcFile = "../../../LHC10b2/v0or/spd/analysis_esd_raw.root")
{
  Float_t eta = 1.29;
  Int_t binBegin = ((TH2*) gFile->Get("fEtaPhi"))->GetXaxis()->FindBin(-eta);
  Int_t binEnd   = ((TH2*) gFile->Get("fEtaPhi"))->GetXaxis()->FindBin(eta);
  
  Printf("eta range: %f bins: %d %d", eta, binBegin, binEnd);
  
  if (!all)
    Printf("Eta smaller than 0 side");
  
  c = new TCanvas;
  TFile::Open("analysis_esd_raw.root");
  hist = ((TH2*) gFile->Get("fEtaPhi"))->ProjectionY("hist", binBegin, (all) ? binEnd : 40);
  hist->Rebin(2);
  hist->SetStats(0);
  hist->Sumw2();
  hist->Draw("HIST");
  gPad->SetGridx();
  gPad->SetGridy();
  
  TFile::Open(mcFile);  
  mcHist = ((TH2*) gFile->Get("fEtaPhi"))->ProjectionY("mcHist", binBegin, (all) ? binEnd : 40);
  mcHist->Rebin(2);
  mcHist->SetLineColor(2);
  mcHist->Scale(hist->Integral() / mcHist->Integral());
  mcHist->Draw("SAME");
  
  Float_t add = 0;
  Int_t bins;
  
  Float_t okRangeBegin[] = { 0.04, 0.67, 1.34 };
  Float_t okRangeEnd[] =   { 0.55, 1.24, 1.63 };
  Float_t gapRangeBegin[] = { 0.6, 1.27  };
  Float_t gapRangeEnd[] =   { 0.65, 1.32 };
  Float_t averageOK  = FitAverage(hist, 3, okRangeBegin, okRangeEnd, 1, bins);
  Float_t averageGap = FitAverage(hist, 2, gapRangeBegin, gapRangeEnd, 2, bins);
  add += bins * (averageOK - averageGap);
  Printf("Average OK: %f %f %d: %f", averageOK, averageGap, bins, add);
  
  Float_t okRangeBegin2[] = { 2.4, 2.65 };
  Float_t okRangeEnd2[] =   { 2.55, 3.2 };
  Float_t gapRangeBegin2[] = { 2.59, 3.3 };
  Float_t gapRangeEnd2[] =   { 2.61, 3.3 };
  averageOK  = FitAverage(hist, 2, okRangeBegin2, okRangeEnd2, 1, bins);
  averageGap = FitAverage(hist, 2, gapRangeBegin2, gapRangeEnd2, 2, bins);
  add += bins * (averageOK - averageGap);
  Printf("Average OK: %f %f %d: %f", averageOK, averageGap, bins, add);
  
  Float_t okRangeBegin3[] = { 3.55, 3.9 };
  Float_t okRangeEnd3[] =   { 3.8, 4.15 };
  Float_t gapRangeBegin3[] = { 3.83  };
  Float_t gapRangeEnd3[] =   { 3.86 };
  averageOK  = FitAverage(hist, 2, okRangeBegin3, okRangeEnd3, 1, bins);
  averageGap = FitAverage(hist, 1, gapRangeBegin3, gapRangeEnd3, 2, bins);
  add += bins * (averageOK - averageGap);
  Printf("Average OK: %f %f %d: %f", averageOK, averageGap, bins, add);
  
  Float_t okRangeBegin4[] = { 4.2, 4.5 };
  Float_t okRangeEnd4[] =   { 4.4, 4.7 };
  Float_t gapRangeBegin4[] = { 4.45  };
  Float_t gapRangeEnd4[] =   { 4.45 };
  averageOK  = FitAverage(hist, 2, okRangeBegin4, okRangeEnd4, 1, bins);
  averageGap = FitAverage(hist, 1, gapRangeBegin4, gapRangeEnd4, 2, bins);
  add += bins * (averageOK - averageGap);
  Printf("Average OK: %f %f %d: %f", averageOK, averageGap, bins, add);
  
  Float_t okRangeBegin5[] = { 5.4, 5.7 };
  Float_t okRangeEnd5[] =   { 5.6, 5.8 };
  Float_t gapRangeBegin5[] = { 5.63  };
  Float_t gapRangeEnd5[] =   { 5.67 };
  averageOK  = FitAverage(hist, 2, okRangeBegin5, okRangeEnd5, 1, bins);
  averageGap = FitAverage(hist, 1, gapRangeBegin5, gapRangeEnd5, 2, bins);
  add += bins * (averageOK - averageGap);
  Printf("Average OK: %f %f %d: %f", averageOK, averageGap, bins, add);
  
  Printf("This adds %.2f %% to the total number of tracklets (%f)", 100.0 * add / hist->Integral(), hist->Integral());
  c->SaveAs("gap1.png");
  
  add1 = add;
  total1 = hist->Integral();

  if (all)
    return;

  Printf("\nEta larger than 0 side");
  
  c = new TCanvas;
  TFile::Open("analysis_esd_raw.root");
  hist = ((TH2*) gFile->Get("fEtaPhi"))->ProjectionY("hist2", 41, binEnd);
  hist->Rebin(2);
  hist->SetStats(0);
  hist->Sumw2();
  hist->Draw("HIST");
  gPad->SetGridx();
  gPad->SetGridy();
  
  TFile::Open(mcFile);  
  mcHist = ((TH2*) gFile->Get("fEtaPhi"))->ProjectionY("mcHist", 41, binEnd);
  mcHist->Rebin(2);
  mcHist->SetLineColor(2);
  mcHist->Scale(hist->Integral() / mcHist->Integral());
  mcHist->Draw("SAME");
  
  add = 0;
  
  Float_t okRangeBegin[] = { 0.04, 0.67, 1.34 };
  Float_t okRangeEnd[] =   { 0.55, 1.24, 1.63 };
  Float_t gapRangeBegin[] = { 0.6, 1.27  };
  Float_t gapRangeEnd[] =   { 0.65, 1.32 };
  Float_t averageOK  = FitAverage(hist, 3, okRangeBegin, okRangeEnd, 1, bins);
  Float_t averageGap = FitAverage(hist, 2, gapRangeBegin, gapRangeEnd, 2, bins);
  add += bins * (averageOK - averageGap);
  Printf("Average OK: %f %f %d: %f", averageOK, averageGap, bins, add);
  
  Float_t okRangeBegin2[] = { 2.32, 2.65 };
  Float_t okRangeEnd2[] =   { 2.55, 3.2 };
  Float_t gapRangeBegin2[] = { 2.59, 3.3 };
  Float_t gapRangeEnd2[] =   { 2.61, 3.3 };
  averageOK  = FitAverage(hist, 2, okRangeBegin2, okRangeEnd2, 1, bins);
  averageGap = FitAverage(hist, 2, gapRangeBegin2, gapRangeEnd2, 2, bins);
  add += bins * (averageOK - averageGap);
  Printf("Average OK: %f %f %d: %f", averageOK, averageGap, bins, add);
  
  Float_t okRangeBegin3[] = { 3.55, 3.9 };
  Float_t okRangeEnd3[] =   { 3.8, 4.15 };
  Float_t gapRangeBegin3[] = { 3.83  };
  Float_t gapRangeEnd3[] =   { 3.86 };
  averageOK  = FitAverage(hist, 2, okRangeBegin3, okRangeEnd3, 1, bins);
  averageGap = FitAverage(hist, 1, gapRangeBegin3, gapRangeEnd3, 2, bins);
  add += bins * (averageOK - averageGap);
  Printf("Average OK: %f %f %d: %f", averageOK, averageGap, bins, add);
  
  Float_t okRangeBegin4[] = { 4.2, 4.5 };
  Float_t okRangeEnd4[] =   { 4.4, 4.7 };
  Float_t gapRangeBegin4[] = { 4.45  };
  Float_t gapRangeEnd4[] =   { 4.45 };
  averageOK  = FitAverage(hist, 2, okRangeBegin4, okRangeEnd4, 1, bins);
  averageGap = FitAverage(hist, 1, gapRangeBegin4, gapRangeEnd4, 2, bins);
  add += bins * (averageOK - averageGap);
  Printf("Average OK: %f %f %d: %f", averageOK, averageGap, bins, add);

  Printf("This adds %.2f %% to the total number of tracklets (%f)", 100.0 * add / hist->Integral(), hist->Integral());
  c->SaveAs("gap2.png");
  
  Printf("In total we add %.2f %%.", 100.0 * (add1 + add) / (total1 + hist->Integral()));
}
