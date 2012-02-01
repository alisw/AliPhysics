/**
 * Test script to fit the energy loss spectra 
 * 
 * @deprecated
 * This is a simple test script 
 *
 * @ingroup pwg2_forward_scripts_tests
 */
#ifndef __CINT__
# include "AliForwardUtil.h"
# include <TFile.h>
# include <TList.h>
# include <TH1.h>
# include <TF1.h>
# include <TFitResult.h>
# include <TMath.h>
# include <TStyle.h>
# include <TArrow.h>
# include <TCanvas.h>
#else
class TCanvas;
class TFile;
class TH1;
class TList;
class TF1;
#endif

//____________________________________________________________________
/** 
 * 
 * 
 * @param ef 
 * @param d 
 * @param r 
 * @param etabin 
 * 
 * @return 
 *
 * @deprecated
 * This is a simple test script 
 *
 * @ingroup pwg2_forward_scripts_tests
 */
TH1* GetEDist(TList* ef, UShort_t d, Char_t r, UShort_t etabin)
{
  TList* dL = static_cast<TList*>(ef->FindObject(Form("FMD%d%c",d,r)));
  if (!dL) {
    Error("GetEDist", "Couldn't find list FMD%d%c",d,r);
    ef->ls();
    return 0;
  }
  if (etabin > 999) {
    TH1* hist = static_cast<TH1*>(dL->FindObject(Form("FMD%d%c_edist",d,r)));
    if (hist) {
      Error("GetEDist", "Couldn't find EDists histogram for FMD%d%c",d,r);
      return 0;
    }
  }
      
  TList* edL = static_cast<TList*>(dL->FindObject("EDists"));
  if (!edL) {
    Error("GetEDist", "Couldn't find list EDists for FMD%d%c",d,r);
    return 0;
  }
  
  TH1* hist = static_cast<TH1*>(edL->FindObject(Form("FMD%d%c_etabin%03d", 
						     d, r, etabin)));
  if (!hist) {
    Error("GetEDist", "Couldn't find histogra FMD%d%c_etabin%03d",
	  d,r, etabin);
    return 0;
  }
  
  return hist;
}

//____________________________________________________________________
/** 
 * 
 * 
 * @param ef 
 * @param d 
 * @param r 
 * @param eta 
 * 
 * @return  
 *
 * @deprecated
 * This is a simple test script 
 *
 * @ingroup pwg2_forward_scripts_tests
*/
TH1* GetEDist(TList* ef, UShort_t d, Char_t r, Float_t eta)
{
  if (!ef) { 
    Error("GetEDist", "EF not set");
    return 0;
  }
  TAxis* etaAxis = static_cast<TAxis*>(ef->FindObject("etaAxis"));
  if (!etaAxis) { 
    Error("GetEDist", "Couldn't find eta axis in list");
    return 0;
  }

  UShort_t bin = etaAxis->FindBin(eta);
  if (bin <= 0 || bin > etaAxis->GetNbins()) { 
    Error("GetEDist", "eta=%f out of range [%f,%f] - getting ring histo", 
	  eta, etaAxis->GetXmin(), etaAxis->GetXmax());
    return GetEDist(ef, d, r, UShort_t(1000));
  }

  return GetEDist(ef, d, r, bin);
}

//____________________________________________________________________
/** 
 * 
 * 
 * @param file 
 * 
 * @return 
 *
 * @deprecated
 * This is a simple test script 
 *
 * @ingroup pwg2_forward_scripts_tests
 */
TList* GetEF(TFile* file) 
{
  TList* forward = static_cast<TList*>(file->Get("PWG2forwardDnDeta/Forward"));
  if (!forward) {
    Error("GetEF", "Failed to get list PWG2forwardDnDeta/Forward from file");
    return 0;
  }
  TList* ef = static_cast<TList*>(forward->FindObject("fmdEnergyFitter"));
  if (!ef) {
    Error("GetEF", "Failed to get energy fitter list");
    return 0;
  }
  
  return ef;
}

//____________________________________________________________________
TList* ef = 0;

//____________________________________________________________________
/** 
 * 
 * 
 * 
 * @return 
 *
 * @deprecated
 * This is a simple test script 
 *
 * @ingroup pwg2_forward_scripts_tests
 */
TList*  CheckEF()
{
  if (ef) return ef;
  
  TFile* file = TFile::Open("AnalysisResults.root", "READ");
  if (!file) { 
    Error("Fit1", "Failed to open file");
    return 0;
  }
  return GetEF(file);
}

//____________________________________________________________________
TCanvas* c = 0;

//____________________________________________________________________
/** 
 * 
 * 
 * 
 * @return 
 *
 * @deprecated
 * This is a simple test script 
 *
 * @ingroup pwg2_forward_scripts_tests
 */
TCanvas* CheckC()
{
  if (c) return c;

  gStyle->SetOptFit(1111);
  gStyle->SetCanvasColor(0);
  gStyle->SetCanvasBorderSize(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasDefW(800);
  gStyle->SetCanvasDefH(800);
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetTitleStyle(0);
  gStyle->SetStatBorderSize(1);
  gStyle->SetStatColor(0);
  gStyle->SetStatStyle(0);
  gStyle->SetStatX(0.95);
  gStyle->SetStatY(0.95);
  gStyle->SetStatW(0.15);
  gStyle->SetStatH(0.15);

  c = new TCanvas("c", "c");
  c->SetLogy();

  return c;
}

//____________________________________________________________________
/** 
 * 
 * 
 * @param f 
 *
 * @deprecated
 * This is a simple test script 
 *
 * @ingroup pwg2_forward_scripts_tests
 */
void PrintFit(TF1* f)
{
  Int_t    nu     = f->GetNDF();
  Double_t chi2   = f->GetChisquare();
  Double_t chi2nu = (nu > 0 ? chi2/nu : 0);
  Printf("%s: chi^2/nu=%f/%d=%f [%f,%f]", 
	 f->GetName(), chi2, nu, chi2nu,
	 f->GetXmin(), f->GetXmax());
  for (Int_t i = 0; i < f->GetNpar(); i++) { 
    Double_t v = f->GetParameter(i);
    Double_t e = f->GetParError(i);
    Double_t r = 100*(v == 0 ? 1 : e / v);
    Printf("%32s = %14.7f +/- %14.7f (%5.1f%%)",f->GetParName(i),v,e,r);
  }
}
  
//____________________________________________________________________
/** 
 * 
 * 
 * @param n 
 * @param d 
 * @param r 
 * @param eta 
 *
 * @deprecated
 * This is a simple test script 
 *
 * @ingroup pwg2_forward_scripts_tests
 */
void TestFitELoss(Int_t n, UShort_t d, Char_t r, Float_t eta)
{
  TList* ef1 = CheckEF();
  TCanvas* c1 = CheckC();
  if (!ef1) return;
  if (!c1)  return;

  TH1* dist = GetEDist(ef1, d, r, eta);
  if (!dist) return;

  AliForwardUtil::ELossFitter f(0.4, 10, 4);
  TF1* landau1 = new TF1(*f.Fit1Particle(dist, 0));
  landau1->SetName("Landau1");
  PrintFit(landau1);
  TF1* landaun = new TF1(*f.FitNParticle(dist, n, 0)); 
  landau1->SetName(Form("Landau%d", n));
  PrintFit(landaun);
  landau1->SetRange(0,10);
  landaun->SetRange(0,10);
  landau1->SetLineWidth(4);
  landaun->SetLineWidth(4);
  landau1->SetLineColor(kBlack);
  landaun->SetLineColor(kBlack);

  dist->GetListOfFunctions()->Clear();
  dist->GetListOfFunctions()->Add(landau1);
  dist->GetListOfFunctions()->Add(landaun);
  dist->GetListOfFunctions()->ls();
  dist->Draw();
  landau1->Draw("same");
  landaun->Draw("same");

  Double_t mp = landaun->GetParameter(1);
  Double_t xi = landaun->GetParameter(2);
  for (Int_t i  = 1; i <= n; i++) { 
    Double_t x  = i * (mp + xi * TMath::Log(i));
    Double_t y  = landaun->Eval(x);
    Double_t y1 = y < 0.05 ? 1 : 0.01;
    TArrow* a = new TArrow(x,y1,x,y,0.03,"|>");
    Info("FitSteer", "Delta_{p,%d}=%f", i, x);
    a->SetLineWidth(2);
    a->SetAngle(30);
    a->Draw();
  }

  c1->cd();
}
//
// EOF
//
