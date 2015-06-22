#include "TH2D.h"
#include "THnSparse.h"
#include "TCanvas.h"
#include "TList.h"
#include "TString.h"
#include "TFile.h"
#include "TGraph.h"
#include "TMath.h"
#include "TObjArray.h"
#include "THnSparse.h"
#include "TSpline.h"

#include <iostream>
#include <iomanip>

//__________________________________________________________________________________________
TSpline3* loadSplines(TFile* file, TObjArray* arr, TString name, Bool_t useOADB)
{
  if (useOADB) {
    if (!arr)
      return 0x0;
    
    return (TSpline3*)arr->FindObject(name.Data());
  }
  
  if (!file)
    return 0x0;
  
  return (TSpline3*)file->Get(name.Data());
}



//__________________________________________________________________________________________
Int_t compareSplines(TString pathNameSplines1, TString pathNameSplines2,
                     TString period = "LHC10D_PASS2", TString dataType1 = "DATA", TString beamType1 = "PP",
                     TString period2 = "", TString dataType2 = "", TString beamType2 = "",
                     Bool_t useOADBforFirstSplines = kFALSE, Bool_t useOADBforSecondSplines = kFALSE,
                     TString displayNameSplines1 = "", TString displayNameSplines2 = "",  Bool_t onlyPlotFirstSplines = kFALSE)
{ 
  if (period2.IsNull())
    period2 = period;
  
  if (dataType2.IsNull())
    dataType2 = dataType1;
  
  if (beamType2.IsNull())
    beamType2 = beamType1;
  
  if (displayNameSplines1.IsNull())
    displayNameSplines1 = "Spline 1";
  
  if (displayNameSplines2.IsNull())
    displayNameSplines2 = "Spline 2";
  
  TFile* f1 = 0x0;
  f1 = TFile::Open(pathNameSplines1.Data());
  if (!f1)  {
    std::cout << "Failed to open file \"" << pathNameSplines1.Data() << "\"!" << std::endl;
    return -1;
  }
  
  TObjArray* arr = 0x0;
  
  if (useOADBforFirstSplines) {
    arr = (TObjArray*)f1->Get("TPCPIDResponse");
    
    if (!arr)  {
      std::cout << "Failed to load first array \"TPCPIDResponse\"!" << std::endl;
      return -1;
    }
  }

  TSpline3* splPion = loadSplines(f1, arr, Form("TSPLINE3_%s_PION_%s_%s_MEAN", dataType1.Data(), period.Data(), beamType1.Data()), 
                                  useOADBforFirstSplines);
  
  TSpline3* splProton = loadSplines(f1, arr, Form("TSPLINE3_%s_PROTON_%s_%s_MEAN", dataType1.Data(), period.Data(), beamType1.Data()), 
                                    useOADBforFirstSplines);
  
  TSpline3* splKaon = loadSplines(f1, arr, Form("TSPLINE3_%s_KAON_%s_%s_MEAN", dataType1.Data(), period.Data(), beamType1.Data()), 
                                  useOADBforFirstSplines);
  
  TSpline3* splElectron = loadSplines(f1, arr, Form("TSPLINE3_%s_ELECTRON_%s_%s_MEAN", dataType1.Data(), period.Data(), beamType1.Data()), 
                                      useOADBforFirstSplines);
  
  
  if (onlyPlotFirstSplines) {
    TCanvas* c = new TCanvas("c", "Splines 1",  100,10,1380,800);
    c->SetTopMargin(0.02);
    c->SetBottomMargin(0.14);
    c->SetLeftMargin(0.1);
    c->SetRightMargin(0.025);
    c->SetLogx(kTRUE);
    
    const Int_t nBins = 1000;
    const Double_t fromLow = 0.15;
    const Double_t toHigh = 100.;
    const Double_t factor = TMath::Power(toHigh/fromLow, 1./nBins);

    // Log binning 
    Double_t pBins[nBins + 1];
    pBins[0] = fromLow;
    for (Int_t i = 0 + 1; i <= nBins; i++)
      pBins[i] = factor * pBins[i - 1];
    
    TH2D* hDummy = new TH2D("hDummy", ";#it{p} (GeV/#it{c});<d#it{E}/d#it{x}> (arb. unit)", nBins, pBins, 1520, 42, 200);
    hDummy->GetYaxis()->SetLabelSize(0.05);
    hDummy->GetYaxis()->SetTitleSize(0.07);
    hDummy->GetYaxis()->SetTitleOffset(0.7);
    hDummy->GetXaxis()->SetNoExponent(kTRUE);
    hDummy->GetXaxis()->SetMoreLogLabels(kTRUE);
    hDummy->GetXaxis()->SetTitleSize(0.07);
    hDummy->GetXaxis()->SetLabelSize(0.05);
    hDummy->GetXaxis()->SetTitleOffset(0.85);
    hDummy->SetStats(kFALSE);
    hDummy->Draw("colz");
    
    TGraph* gPion = new TGraph(nBins + 1);
    gPion->SetTitle("#pi");
    gPion->SetFillStyle(0);
    gPion->SetFillColor(kWhite);
    gPion->SetLineWidth(2);
    for (Int_t i = 0; i < nBins + 1; i++)
      gPion->SetPoint(i, pBins[i], 50.*(splPion->Eval(pBins[i] / AliPID::ParticleMass(AliPID::kPion))));
    gPion->SetLineColor(kRed);
    gPion->SetMarkerColor(kRed);
    gPion->Draw("same");
    
    
    TGraph* gProton = new TGraph(nBins + 1);
    gProton->SetTitle("p");
    gProton->SetFillStyle(0);
    gProton->SetFillColor(kWhite);
    gProton->SetLineWidth(2);
    for (Int_t i = 0; i < nBins + 1; i++) {
      gProton->SetPoint(i, pBins[i], 50.*(splProton->Eval(pBins[i] / AliPID::ParticleMass(AliPID::kProton))));
    }
    gProton->SetLineColor(kBlue);
    gProton->SetMarkerColor(kBlue);
    gProton->Draw("same");
    
    
    TGraph* gKaon = new TGraph(nBins + 1);
    gKaon->SetTitle("K");
    gKaon->SetFillStyle(0);
    gKaon->SetFillColor(kWhite);
    gKaon->SetLineWidth(2);
    for (Int_t i = 0; i < nBins + 1; i++) {
      gKaon->SetPoint(i, pBins[i], 50.*(splKaon->Eval(pBins[i] / AliPID::ParticleMass(AliPID::kKaon))));
    }
    gKaon->SetLineColor(kGreen);
    gKaon->SetMarkerColor(kGreen);
    gKaon->Draw("same");
    
    
    TGraph* gElectron = new TGraph(nBins + 1);
    gElectron->SetTitle("e");
    gElectron->SetFillStyle(0);
    gElectron->SetFillColor(kWhite);
    gElectron->SetLineWidth(2);
    for (Int_t i = 0; i < nBins + 1; i++) {
      gElectron->SetPoint(i, pBins[i], 50.*(splElectron->Eval(pBins[i] / AliPID::ParticleMass(AliPID::kElectron))));
    }
    gElectron->SetLineColor(kMagenta);
    gElectron->SetMarkerColor(kMagenta);
    gElectron->Draw("same");
    
    TLegend* leg = new TLegend(0.65, 0.65, 0.85, 0.88);
    leg->SetNColumns(2);
    leg->SetTextSize(0.06);
    leg->SetHeader(displayNameSplines1.Data());
    leg->AddEntry(gPion, Form(" %s", gPion->GetTitle()), "l");
    leg->AddEntry(gKaon, Form(" %s", gKaon->GetTitle()), "l");
    leg->AddEntry(gProton, Form(" %s", gProton->GetTitle()), "l");
    leg->AddEntry(gElectron, Form(" %s", gElectron->GetTitle()), "l");
    leg->SetFillColor(kWhite);
    leg->SetBorderSize(0.);
    c->SetGridx(kFALSE);
    c->SetGridy(kFALSE);
    leg->Draw();
    
    c->RedrawAxis(""); // To get the ticks visible again
    //c->RedrawAxis("G"); // To get the grid visible again
    
    c->Update();
  
    return 0;
  }
  
  
  TFile* f2 = 0x0;
  f2 = TFile::Open(pathNameSplines2.Data());
  if (!f2)  {
    std::cout << "Failed to open file \"" << pathNameSplines2.Data() << "\"!" << std::endl;
    return -1;
  }

  TCanvas* c = new TCanvas("c", "",  100,10,1380,800);
  c->SetLogx(kTRUE);LHC10E_PASS2
  
  const Int_t nPoints = 20000;
  const Float_t stepSize = 0.02;
  const Float_t stepSizeEl = stepSize * 100.;
  const Int_t nPointsEl = nPoints * 10.;
  
  TH2D* hDummy = new TH2D("hDummy", Form("; #it{p} (GeV/#it{c}); %s / %s", displayNameSplines1.Data(), displayNameSplines2.Data()),
                          1000, 0.15, 60, 1000, 0.9, 1.1);
  hDummy->GetYaxis()->SetLabelSize(0.03);
  hDummy->GetYaxis()->SetTitleSize(0.05);
  hDummy->GetYaxis()->SetTitleOffset(1.);
  hDummy->GetXaxis()->SetNoExponent(kTRUE);
  hDummy->GetXaxis()->SetMoreLogLabels(kTRUE);
  hDummy->GetXaxis()->SetTitleSize(0.05);
  hDummy->GetXaxis()->SetLabelSize(0.04);
  hDummy->GetXaxis()->SetTitleOffset(1.1);
  hDummy->SetStats(kFALSE);
  hDummy->Draw("colz");
  
  
  
  
  
  TSpline3* splPion2 = 0x0;
  TSpline3* splKaon2 = 0x0;
  TSpline3* splElectron2 = 0x0;
  TSpline3* splProton2 = 0x0;
  
  TObjArray* arr2 = 0x0;
  
  if (useOADBforSecondSplines) {
    arr2 = (TObjArray*)f2->Get("TPCPIDResponse");
    
    if (!arr2)  {
      std::cout << "Failed to second array \"TPCPIDResponse\"!" << std::endl;
      return -1;
    }
  }
  
  splPion2 = loadSplines(f2, arr2, Form("TSPLINE3_%s_PION_%s_%s_MEAN", dataType2.Data(), period2.Data(), beamType2.Data()), useOADBforSecondSplines);
  TGraph* gPion = new TGraph(nPoints);
  gPion->SetTitle("#pi");
  gPion->SetFillStyle(0);
  gPion->SetFillColor(kWhite);
  for (Int_t i = 0; i < nPoints; i++) {
    gPion->SetPoint(i, (0. + i * stepSize) * 0.1396, (splPion->Eval((0. + i * stepSize)) / splPion2->Eval((0. + i * stepSize))));
  }
  gPion->SetLineColor(kRed);
  gPion->SetMarkerColor(kRed);
  gPion->Draw("same");


  splProton2 = loadSplines(f2, arr2, Form("TSPLINE3_%s_PROTON_%s_%s_MEAN", dataType2.Data(), period2.Data(), beamType2.Data()), useOADBforSecondSplines);
  TGraph* gProton = new TGraph(nPoints);
  gProton->SetTitle("p");
  gProton->SetFillStyle(0);
  gProton->SetFillColor(kWhite);
  for (Int_t i = 0; i < nPoints; i++) {
    gProton->SetPoint(i, (0. + i * stepSize) * 0.938, (splProton->Eval((0. + i * stepSize)) / splProton2->Eval((0. + i * stepSize))));
  }
  gProton->SetLineColor(kBlue);
  gProton->SetMarkerColor(kBlue);
  gProton->Draw("same");


  splKaon2 = loadSplines(f2, arr2, Form("TSPLINE3_%s_KAON_%s_%s_MEAN", dataType2.Data(), period2.Data(), beamType2.Data()), useOADBforSecondSplines);
  TGraph* gKaon = new TGraph(nPoints);
  gKaon->SetTitle("K");
  gKaon->SetFillStyle(0);
  gKaon->SetFillColor(kWhite);
  for (Int_t i = 0; i < nPoints; i++) {
    gKaon->SetPoint(i, (0. + i * stepSize) * 0.49368, (splKaon->Eval((0. + i * stepSize)) / splKaon2->Eval((0. + i * stepSize))));
  }
  gKaon->SetLineColor(kGreen);
  gKaon->SetMarkerColor(kGreen);
  gKaon->Draw("same");

  
  splElectron2 = loadSplines(f2, arr2, Form("TSPLINE3_%s_ELECTRON_%s_%s_MEAN", dataType2.Data(), period2.Data(), beamType2.Data()),
                             useOADBforSecondSplines);
  TGraph* gElectron = new TGraph(nPoints);
  gElectron->SetTitle("e");
  gElectron->SetFillStyle(0);
  gElectron->SetFillColor(kWhite);
  for (Int_t i = 0; i < nPointsEl; i++) {
    gElectron->SetPoint(i, (0. + i * stepSizeEl) * 0.000511, 
                        (splElectron->Eval((0. + i * stepSizeEl)) / splElectron2->Eval((0. + i * stepSizeEl))));
  }
  gElectron->SetLineColor(kMagenta);
  gElectron->SetMarkerColor(kMagenta);
  gElectron->Draw("same");
  
  TLegend* leg = c->BuildLegend();
  leg->GetListOfPrimitives()->RemoveAt(0);
  leg->SetFillColor(kWhite);
  c->SetGridx(kTRUE);
  c->SetGridy(kTRUE);
  
  c->Update();
  
  return 0;
}