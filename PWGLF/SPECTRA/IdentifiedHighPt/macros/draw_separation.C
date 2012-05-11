#include <TFile.h>
#include <TList.h>
#include <TTree.h>
#include <TChain.h>
#include <TClonesArray.h>
#include <TObjString.h>
#include <TString.h>
#include <TROOT.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TCanvas.h>

#include "my_tools.C"
#include "my_functions.C"

#include <iostream>
#include <fstream>
#include <string>

using namespace std;

/*
  To run calibrations:
  ====================

  Use root:

  .L libMyDeDxAnalysis.so 
  .L my_functions.C+
  .L my_tools.C+
  .L draw_separation.C+

  DrawSeparation("fitparameters/lhc10h_aod_all.root", 0, 50)
 */


TF1* piFunc = 0;
TF1* kFunc  = 0;
TF1* pFunc = 0;
TF1* sigmaFunc = 0;

Double_t Sep(Double_t* xx, Double_t* par);

//____________________________________________________________________________
void DrawSeparation(const Char_t* fitFileName, Double_t pLow, Double_t pHigh)
{
  gStyle->SetOptStat(0);
  
  TFile* fitFile = FindFileFresh(fitFileName);
  if(!fitFile)
    return;
  DeDxFitInfo* fitPar = (DeDxFitInfo*)fitFile->Get("fitInfo");
  fitPar->Print();
  
  fixMIP      = fitPar->MIP;
  fixPlateau  = fitPar->plateau;
  
  Double_t dedxPar[6]  = {0, 0, 0, 0, 0, 0};
  Double_t sigmaPar[6] = {0, 0, 0, 0, 0, 0};
  
  dedxPar[0] = fitPar->optionDeDx;
  for(Int_t i = 0; i < fitPar->nDeDxPar; i++) {
    dedxPar[i+1] = fitPar->parDeDx[i];
  }
  
  sigmaPar[0] = fitPar->optionSigma;
  for(Int_t i = 0; i < fitPar->nSigmaPar; i++) {
    sigmaPar[i+1] = fitPar->parSigma[i];
  }
  
  piFunc = new TF1("piFunc", FitFunc, 0, 100, fitPar->nDeDxPar+1);
  piFunc->SetParameters(dedxPar);
  
  kFunc = new TF1("kFunc", FitFunc, 0, 100, fitPar->nDeDxPar+1);
  kFunc->SetParameters(dedxPar);
  kFunc->SetParameter(0, kFunc->GetParameter(0)+10);
  
  pFunc = new TF1("pFunc", FitFunc, 0, 100, fitPar->nDeDxPar+1);
  pFunc->SetParameters(dedxPar);
  pFunc->SetParameter(0, pFunc->GetParameter(0)+20);
  
  sigmaFunc = new TF1("sigmaFunc", SigmaFunc, 0, 100, fitPar->nSigmaPar+1); 
  sigmaFunc->SetParameters(sigmaPar);

  TCanvas* c1 = new TCanvas("c1", "c1"); 
  
  TH1F* hist = new TH1F("hist", "Separation in pp vs p; p [GeV/c]; Separation",
			100, 0, pHigh);
  hist->SetMinimum(0.0);
  hist->SetMaximum(6.0);
  hist->Draw();
  
  TLegend* legend = new TLegend(0.74, 0.64, 0.89, 0.89);    
  legend->SetBorderSize(0);
  legend->SetFillColor(0);

  TF1* pipFunc = new TF1("pipFunc", Sep,
			 pLow, pHigh, 1);
  pipFunc->SetParameter(0, 0);
  pipFunc->SetLineColor(2);
  pipFunc->Draw("SAME");

  TF1* pikFunc = new TF1("pikFunc", Sep,
			 pLow, pHigh, 1);
  pikFunc->SetParameter(0, 1);
  pikFunc->SetLineColor(3);
  pikFunc->Draw("SAME");

  TF1* kpFunc = new TF1("kpFunc", Sep,
			 pLow, pHigh, 1);
  kpFunc->SetParameter(0, 2);
  kpFunc->SetLineColor(4);
  kpFunc->Draw("SAME");

  legend->AddEntry(pipFunc, "#pi-p", "L");
  legend->AddEntry(pikFunc, "#pi-K", "L");
  legend->AddEntry(kpFunc, "K-p", "L");
  legend->Draw();
  gROOT->ProcessLine(".x drawText.C");
  c1->SaveAs("separation.gif");
  c1->SaveAs("separation.pdf");
}
  
//______________________________________________________________________________
Double_t Sep(Double_t* xx, Double_t* par)
{
  //
  // Could speed up fit by forcing it to use <p>. In that way the parameters
  // could be amde statis cand only changed when going to a new p bin
  //
  Double_t p = xx[0];

  Int_t option = Int_t(par[0]);

  TF1* f1 = 0;
  TF1* f2 = 0;
  switch (option) {
    
  case 0: // pi - p
    f1 = piFunc;
    f2 = pFunc;
    break;
  case 1: // pi - k
    f1 = piFunc;
    f2 = kFunc;
    break;
  case 2: // k - p
    f1 = kFunc;
    f2 = pFunc;
    break;
  default:
    cout << "Error in Sep: option " << option << " not supported!!!!!" << endl;
    return 0;
    break;
  }

  Double_t dedx1  = f1->Eval(p);
  Double_t dedx2  = f2->Eval(p);
  Double_t sigma1 = sigmaFunc->Eval(dedx1);
  Double_t sigma2 = sigmaFunc->Eval(dedx2);
  
  return (dedx1-dedx2)/TMath::Sqrt(sigma1*sigma2);
}
