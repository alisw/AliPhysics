#include <TFile.h>
#include <TH1.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TList.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TMath.h>
#include <TSystem.h>
#include <TLatex.h>
#include <TASImage.h>
#include <TLine.h>
#include <TBox.h>

#include "my_tools.C"


/*
  Info:
  Multiplies the fraction with the charged spectrum. 
  Applies a rapidity correction in case of K and p.

  To run:

  gSystem->AddIncludePath("-I../macros")
  gSystem->AddIncludePath("-I../lib")
  gROOT->SetMacroPath(".:../macros")
  .L my_tools.C+
  .L makeSpectra.C+

  MakeSpectra("../corrected_fraction/fraction_7tev_b_all.root", "NOPID_pp_7TeV.root", "final_spectra_7tev_b.root");


 */

Double_t rap_correction(Double_t* x, Double_t* par);

void MakeSpectra(const Char_t* fileName,
		 const Char_t* fileNameCharged,
		 const Char_t* outFileName,
		 const Char_t* endName="PP")
{
  TFile* fileData = FindFileFresh(fileName);

  if(!fileData)
    return;

  TFile* fileDataCharged = FindFileFresh(fileNameCharged);
  
  if(!fileDataCharged)
    return;
  
  const Int_t nPid = 3;
  const Char_t* pidNames[nPid] = {"Pion", "Kaon", "Proton"}; 
  const Char_t* titleNames[nPid] = {"#pi^{+} + #pi^{-}", "K^{+} + K^{-}", "p + #bar{p}"}; 
  const Double_t mass[nPid] = {0.140, 0.494, 0.938};

  TH1D* hPtSpectrum[nPid] = {0, 0, 0};
  TH1D* hPtSpectrumSyst[nPid] = {0, 0, 0};

  TH1D* hPtCharged = (TH1D*)fileDataCharged->Get(Form("hDnDptCharged_%s", endName));
  TH1D* hPtChargedSyst = (TH1D*)fileDataCharged->Get(Form("hDnDptChargedSyst_%s", endName));

  TF1* fRap = new TF1("fRap", rap_correction, 0.0, 50.0, 2);

  for(Int_t pid = 0; pid < nPid; pid++) {
    
    // Pions
    TH1D* hFraction     = (TH1D*)fileData->Get(Form("h%sFraction_%s", pidNames[pid], endName));
    TH1D* hFractionSyst = (TH1D*)fileData->Get(Form("h%sFractionSyst_%s", pidNames[pid], endName));
    
    hPtSpectrum[pid] = (TH1D*)hFraction->Clone(Form("h%sSpectrum_%s", pidNames[pid], endName));
    hPtSpectrum[pid]->GetYaxis()->SetTitle(Form("1/N_{evt} 1/(2#pi p_{T}) (d^{2}N_{%s})/(dy dp_{T}) (GeV/c)^{-2}", titleNames[pid]));
    hPtSpectrum[pid]->Multiply(hPtCharged);

    hPtSpectrumSyst[pid] = (TH1D*)hFractionSyst->Clone(Form("h%sSpectrumSyst_%s", pidNames[pid], endName));
    hPtSpectrumSyst[pid]->SetMarkerStyle(1);
    hPtSpectrum[pid]->GetYaxis()->SetTitle(Form("1/N_{evt} 1/(2#pi p_{T}) (d^{2}N_{%s})/(dy dp_{T}) (GeV/c)^{-2}", titleNames[pid]));
    hPtSpectrumSyst[pid]->Multiply(hPtChargedSyst);
    hPtSpectrumSyst[pid]->SetFillStyle(1001);
    hPtSpectrumSyst[pid]->SetFillColor(kGray);

    //
    // Rapidity correction
    //
    fRap->SetParameters(0.8, mass[pid]);
    hPtSpectrum[pid]->Divide(fRap);
    hPtSpectrumSyst[pid]->Divide(fRap);
  }

  TFile* fileOut = new TFile(outFileName, "RECREATE");
  for(Int_t pid = 0; pid < nPid; pid++) {
    hPtSpectrumSyst[pid]->Write();
    hPtSpectrum[pid]->Write();
  }
  fileOut->Close();
}

//______________________________________________________________________________
Double_t rap_correction(Double_t* x, Double_t* par)
{
  Double_t pt = x[0];  

  Double_t eta  = par[0];
  Double_t mass = par[1];

  const Double_t mt = TMath::Sqrt(pt*pt + mass*mass);
  
  const Double_t rap = TMath::ASinH(pt/mt*TMath::SinH(eta));
  
  return rap/eta;
}
