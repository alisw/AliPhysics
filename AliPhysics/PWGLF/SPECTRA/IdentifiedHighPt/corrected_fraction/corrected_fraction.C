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

#include "my_tools.C"


/*
  Info:
  * Efficiencies are currently set by hand!!!!!
  * Systematic errors needs an update also.
  * Need to accomodate k+p at soem point!

  To run:
  
  gSystem->AddIncludePath("-I../macros")
  gSystem->AddIncludePath("-I../lib")
  gROOT->SetMacroPath(".:../macros")
  .L my_tools.C++
  .L corrected_fraction.C+
  
  
  CorrectedFraction("../ratios_7tevb/fit_yields_results/7tev_b.root", "fraction_7tev_b_pi.root", 1, "PP", 0, 10.0);
  CorrectedFraction("../ratios_7tevb/fit_yields_results/7tev_b.root", "fraction_7tev_b_k.root", 2, "PP", 0, 10.0);
  CorrectedFraction("../ratios_7tevb/fit_yields_results/7tev_b.root", "fraction_7tev_b_p.root", 3, "PP", 0, 10.0);


  And in the shell one can add files like:
  
  hadd fraction_7tev_b_all.root fraction_7tev_b_pi.root fraction_7tev_b_k.root fraction_7tev_b_p.root

 */

void CorrectedFraction(const Char_t* inFileName, const Char_t* outFileName, Int_t pid, 
		       const Char_t* endName="PP", Int_t centBin=0, 
		       Double_t ptMax=20.0, const Char_t* outDir="plots");
TH1D* GetSystErrorHist(TH1D* hRatio, Int_t centBin, TF1* electronFraction);


//________________________________________________________________________________________
void CorrectedFraction(const Char_t* inFileName, const Char_t* outFileName, Int_t pid, 
		   const Char_t* endName, Int_t centBin, 
		   Double_t ptMax, const Char_t* outDir)
{
  gStyle->SetOptStat(0);
  
  
  TFile* fileData = FindFileFresh(inFileName);

  if(!fileData)
    return;
  
  TF1*  fElectronFraction   = (TF1*)fileData->Get("fElectronFraction");
  cout << "Electron fraction found!" << endl;
  TH1D* hPionRatio   = 0;
  if(pid==1) {
    hPionRatio = (TH1D*)fileData->Get("hPionRatio");
    CutHistogram(hPionRatio, 2.0, ptMax);
    hPionRatio->SetMarkerStyle(24);
  } else if(pid==2) {
    hPionRatio = (TH1D*)fileData->Get("hKaonRatio");
    CutHistogram(hPionRatio, 3.0, ptMax);
    hPionRatio->SetMarkerStyle(25);
  } else if(pid==3) {
    hPionRatio = (TH1D*)fileData->Get("hProtonRatio");
    CutHistogram(hPionRatio, 3.0, ptMax);
    hPionRatio->SetMarkerStyle(25);
  }

  // Global variable
  TH1D* hSystFraction = 0;
  if(pid==1)
    hSystFraction =  GetSystErrorHist(hPionRatio, centBin, fElectronFraction);
  else
    hSystFraction =  GetSystErrorHist(hPionRatio, centBin, 0);
  hSystFraction->SetLineColor(1);
  hSystFraction->SetMarkerStyle(1);
  hSystFraction->SetFillStyle(1001);
  hSystFraction->SetFillColor(kGray);

  TH1D* hPionFraction = (TH1D*)hPionRatio->Clone("hPionFraction");
  if(pid==1 && !fElectronFraction) {
    TF1 f1("f1", "1.0", 0.0, 50.0);
    cout << "NO ELECTRON FRACTION!!!" << endl;
    hPionFraction->Add(&f1, -0.01); // correct for muons and electrons
    CutHistogram(hPionFraction, 3.0, ptMax);
    hSystFraction->Add(&f1, -0.01); // correct for muons and electrons
    CutHistogram(hSystFraction, 3.0, ptMax);
  }

  if(pid==1 || pid ==3) {

    hPionFraction->Scale(0.94); // correct for efficiency ratio
    hSystFraction->Scale(0.94); // correct for efficiency ratio
  } else {

    TF1* effRatioK = new TF1("effRatioK", "exp(-[1]*x)+[0]", 0.0, 50.0);
    effRatioK->SetParameters(9.82065e-01, 1.28157e+00);
    hPionFraction->Multiply(effRatioK); // correct for efficiency ratio
    hSystFraction->Multiply(effRatioK); // correct for efficiency ratio
  }

  if(pid==1)
    hPionFraction->SetMarkerStyle(20);
  else
    hPionFraction->SetMarkerStyle(20);
  TLatex latex;
  latex.SetNDC();
  latex.SetTextSize(0.05);

  TCanvas* cRatios = new TCanvas("cRatios", "Particle fractions");
  cRatios->Clear();
  hSystFraction->GetXaxis()->SetRangeUser(0, ptMax);
  hSystFraction->GetYaxis()->SetRangeUser(0, 1);
  hSystFraction->Draw("E5");
  hPionRatio->Draw("SAME");
  hPionFraction->Draw("SAME");
  latex.DrawLatex(0.6, 0.95, Form("%s", endName));
  
  CreateDir(outDir);
  if(pid==1) {
    hPionFraction->SetName(Form("hPionFraction_%s", endName));
    hSystFraction->SetName(Form("hPionFractionSyst_%s", endName));
    cRatios->SaveAs(Form("%s/fractions_%s_pions.gif", outDir, endName));
  } else if (pid==2) {
    hPionFraction->SetName(Form("hKaonFraction_%s", endName));
    hSystFraction->SetName(Form("hKaonFractionSyst_%s", endName));
    cRatios->SaveAs(Form("%s/fractions_%s_kaons.gif", outDir, endName));
  } else if (pid==3) {
    hPionFraction->SetName(Form("hProtonFraction_%s", endName));
    hSystFraction->SetName(Form("hProtonFractionSyst_%s", endName));
    cRatios->SaveAs(Form("%s/fractions_%s_protons.gif", outDir, endName));
  }

  TFile* outFile = new TFile(outFileName, "RECREATE");
  hPionFraction->Write();
  hSystFraction->Write();
  outFile->Close();
}

//__________________________________________________________________________________
TH1D* GetSystErrorHist(TH1D* hRatio, Int_t centBin, TF1* electronFraction)
{
  TFile* fileSyst = FindFile("syst_vs_pt.root");
  
  TF1* fSyst = (TF1*)fileSyst->Get(Form("systHigh%d", centBin));
  fSyst->SetRange(2.0, 20.0);
  fSyst->Print();
  
  Double_t syst_error_mc = 0.03; // pp
  if (centBin>0)
    syst_error_mc = 0.05; // Pb+Pb
  
  Double_t syst_error_correction = 0.01; // pp
  if (centBin>0)
    syst_error_correction = 0.03; // Pb+Pb

  TH1D* hSystError = (TH1D*)hRatio->Clone("hPionFractionSyst");
  hSystError->Multiply(fSyst);

  const Int_t nBins = hSystError->GetNbinsX();
  for(Int_t bin = 1; bin <= nBins; bin++) {
  
    Double_t value      = hRatio->GetBinContent(bin);
    Double_t stat_error = hSystError->GetBinContent(bin);

    if(value==0)
      continue;

    Double_t syst_error = stat_error*stat_error;
    syst_error += value*value*syst_error_mc*syst_error_mc;
    syst_error += value*value*syst_error_correction*syst_error_correction;
    
    
    if (electronFraction) {
      
      Double_t systEandMu = electronFraction->Eval(hRatio->GetXaxis()->GetBinCenter(bin));
      syst_error += systEandMu*systEandMu;
    }

    syst_error = TMath::Sqrt(syst_error);
    hSystError->SetBinContent(bin, value);
    hSystError->SetBinError(bin, syst_error);
  }

  return hSystError;
}

