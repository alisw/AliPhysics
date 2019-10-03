
#include <TFile.h>
#include <TH1.h>
#include <TProfile.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TGraphAsymmErrors.h>
#include <TList.h>
#include <TMath.h>
#include <TString.h>
#include <TSystem.h>
#include <TLatex.h>
#include <TF1.h>
#include <TF2.h>
#include <TFitResultPtr.h>
#include <TFitResult.h>

#include "AliHighPtDeDxCalib.h"

#include "my_tools.C"
#include "my_functions.C"

#include <iostream>
#include <fstream>
#include <string>

using namespace std;

/*
  Ideas to improve:
  =================

  Use the real mean p -> Might give some effect
  Effect: Push the <dE/dx> down (let'a see if it is still needed in the fits)

  Use the real <nCl> 
  Effect: Increase sigma for p<15. For p>15 seems to saturate.
  
  To use:

  =======
  root is enough

  gSystem->AddIncludePath("-I$ALICE_ROOT/TPC")
  gSystem->AddIncludePath("-I../lib")
  gSystem->AddIncludePath("-I../grid")
  gSystem->AddIncludePath("-I../macros")
  gROOT->SetMacroPath(".:../macros:../grid:../lib/")
  .L libMyDeDxAnalysis.so 
  .L my_functions.C+
  .L my_tools.C+
  .L DebugClasses.C+
  .L calibrate_de_dx.C+
  
  // Step 2: test calibrations done in step 1
  DrawStep2("calib_eta/7tev_b.root", kFALSE);
  DrawStep2("calib_eta/7tev_b.root", kTRUE);
  TestResolutionVsDeDx("calib_eta/7tev_b.root");
  TestResolutionVsEta("calib_eta/7tev_b.root", kTRUE);

  // Here we want to see that the data is symmetric in eta and roughly symmetric for high vs low eta
  CompareYields("calib_eta/7tev_b.root", 0, 3.0, 30.0, "eta-80", "eta08")
  CompareYields("calib_eta/7tev_b.root", 0, 3.0, 30.0, "etaabs04", "etaabs48") 

  // Step 3: extract
  FitDeDxVsP("calib_eta/7tev_b.root", 3.0, 30.0, 4, 2, kTRUE)


  // Step 4: inspect the results in results/calibdedx and results/calibplots
  // using e.g. gwenview. Event the fit did not converge it is still ok if it
  // describes the data!



  //
  // Test
  //
  // Step 2: test calibrations done in step 1
  DrawStep2("calib_eta/7tev_b_test.root", kFALSE);
  DrawStep2("calib_eta/7tev_b_test.root", kTRUE);
  TestResolutionVsDeDx("calib_eta/7tev_b_test.root");
  TestResolutionVsEta("calib_eta/7tev_b_test.root", kTRUE);
  // Step 3: extract
  FitDeDxVsP("calib_eta/7tev_b_test.root", 2.0, 10.0, 4, 2, kTRUE)


 */


//____________________________________________________________________________
void FitDeDxVsP(const Char_t* calibFileName,
		Double_t pStart, Double_t pStop,
		Int_t optionDeDx, Int_t optionSigma,
		Bool_t performFit = kFALSE,
		Int_t run    = 0,
		Int_t filter = 1,
		Bool_t usePhiCut = kTRUE,
		const Char_t* v0FileName=0,
		Bool_t fixAllPar=kFALSE)
{
  gStyle->SetOptStat(0);

  TString baseName(gSystem->BaseName(calibFileName));
  baseName.ReplaceAll(".root", "");

  TFile* calibFile = FindFileFresh(calibFileName);
  if(!calibFile)
    return;
  AliHighPtDeDxCalib* calib = (AliHighPtDeDxCalib*)GetObject(calibFile, filter, usePhiCut, run);
  calib->Print();

  fixMIP      = calib->GetHistDeDx(kTRUE, 0)->GetMean();
  fixPlateau  = calib->GetHistDeDx(kFALSE, 0)->GetMean();
  hDeDxVsP = calib->GetHistDeDxVsP(0);
  hMeanP = calib->GetHistMeanP();

  AliHighPtDeDxCalib* lambdaCalib = 0;
  AliHighPtDeDxCalib* kaonCalib   = 0;
  TH2D* hDeDxVsPV0Lambda = 0;
  TH2D* hDeDxVsPV0Kaon   = 0;
  if(v0FileName) {
    TFile* v0File = FindFileFresh(v0FileName);
    if(!v0File)
      return;
    lambdaCalib = (AliHighPtDeDxCalib*)GetObject(v0File, 0, usePhiCut, run, "lambda");
    lambdaCalib->Print();
    hDeDxVsPV0Lambda = lambdaCalib->GetHistDeDxVsP(0);
    kaonCalib   = (AliHighPtDeDxCalib*)GetObject(v0File, 0, usePhiCut, run, "kaon");
    kaonCalib->Print();
    hDeDxVsPV0Kaon = kaonCalib->GetHistDeDxVsP(0);
  }

  CreateDir("old");
  gSystem->Exec("mv results/calibdedx/* old/");
  if(calib->IsMc())
    gSystem->Exec("mv debugmc/* old/");

  if(hDeDxVsPV0Lambda) {
    gSystem->Exec("mv debuglambda/* old/");
    gSystem->Exec("mv debugkaon/* old/");
  }

  TH2D* hDeDxVsPPi = 0;
  TH2D* hDeDxVsPK  = 0;
  TH2D* hDeDxVsPP  = 0;
  TH2D* hDeDxVsPMu = 0;

  if(calib->IsMc()) {

    hDeDxVsPPi = calib->GetHistDeDxVsP(1);
    hDeDxVsPK  = calib->GetHistDeDxVsP(2);
    hDeDxVsPP  = calib->GetHistDeDxVsP(3);
    hDeDxVsPMu = calib->GetHistDeDxVsP(5);
  }

  TCanvas* cDeDxVsP = new TCanvas("cDeDxVsP", "dE/dx vs p", 400, 300);
  cDeDxVsP->Clear();
  cDeDxVsP->cd();
  cDeDxVsP->SetLogz();
  hDeDxVsP->SetTitle(0);
  hDeDxVsP->Draw("COLZ");

  TCanvas* cDeDxVsPLogX = new TCanvas("cDeDxVsPLogX", "dE/dx vs p", 400, 300);
  cDeDxVsPLogX->Clear();
  cDeDxVsPLogX->cd();
  cDeDxVsPLogX->SetLogz();
  cDeDxVsPLogX->SetLogx();
  hDeDxVsP->Draw("COLZ");

  // Root is a bit stupid with finidng bins so we have to add and subtract a
  // little to be sure we get the right bin as we typically put edges as
  // limits
  const Int_t binStart = hDeDxVsP->GetXaxis()->FindBin(pStart+0.01);
  pStart = hDeDxVsP->GetXaxis()->GetBinLowEdge(binStart);
  const Int_t binStop  = hDeDxVsP->GetXaxis()->FindBin(pStop-0.01);
  pStop = hDeDxVsP->GetXaxis()->GetBinUpEdge(binStop);
  const Int_t nBins    = binStop - binStart + 1;

  cout << "Doing 2d fit from pTlow = " << pStart << " (bin: " << binStart
       << ") to pThigh = " << pStop << " (bin: " << binStop << ")" << endl;
  
  // Double_t binSize = (histdEdxvsP->GetXaxis()->GetXmax() - histdEdxvsP->GetXaxis()->GetXmin())/ histdEdxvsP->GetXaxis()->GetNbins();
  
  Double_t parDeDx[3]  = {0, 0, 0};
  Double_t parSigma[3] = {0, 0, 0};
  
  const Int_t nLocalParams  = 3; // pi, K, p yields
  Int_t nDeDxPar      = 0;
  Int_t nSigmaPar     = 0;

  switch(optionDeDx) {
    
  case 1:
    nDeDxPar = 2;
    parDeDx[0] = 39.7;
    parDeDx[1] =  6.3;
    break;
  case 2:
    nDeDxPar = 1;
    parDeDx[0] =  7.3;
    break;
  case 3:
    nDeDxPar = 2;
    parDeDx[0] =  6.85097;
    parDeDx[1] =  29.5933;
    break;
  // case 4:
  //   nDeDxPar = 3;
  //   parDeDx[0] = 35.5471;
  //   parDeDx[1] =  6.85097;
  //   parDeDx[2] =  29.5933;
  //   break;
  case 4:
    nDeDxPar = 3;
    parDeDx[0] =  35.6694;
    parDeDx[1] =  6.80835;
    parDeDx[2] =  7.6542;
    break;
  default:

    cout << "optionDeDx does not support option: " << optionDeDx << endl;
    return;
    break;
  }

  switch(optionSigma) {
    
  case 1:
    nSigmaPar = 1;
    parSigma[0] = 3.0;
    break;
  case 2:
    nSigmaPar = 1;
    parSigma[0] = 0.0655688;
    break;
  case 3:
    nSigmaPar = 2;
    parSigma[0] = 0.06;
    parSigma[1] = -0.001;
  case 4:
    nSigmaPar = 2;
    parSigma[0] = 0.06;
    parSigma[1] = 1.0;
    break;
  case 5:
    nSigmaPar = 2;
    parSigma[0] = 0.06;
    parSigma[1] = 1.0;
    break;
  default:

    cout << "optionSigma does not support option: " << optionSigma << endl;
    return;
    break;
  }

  const Int_t nGlobalParams = 2  // binStart, binStop, 
    + 2 + nDeDxPar               // optionDeDx, nDeDxPar, dedxpar0 ....
    + 2 + nSigmaPar;             // nSigmaPar, sigmapar0 .....
  
  const Int_t nParams = nBins*nLocalParams + nGlobalParams;
  
  
  TF2* fitAll = new TF2("fitAll", fitf3G, pStart, pStop, 30, 90, nParams);
  Double_t parametersIn[nParams]; 
  
  parametersIn[0] = binStart;
  fitAll->SetParName(0,"binStart");
  fitAll->FixParameter(0, parametersIn[0]);

  parametersIn[1] = binStop;
  fitAll->SetParName(1,"binStop");
  fitAll->FixParameter(1, parametersIn[1]);

  // dE/dx parameters
  parametersIn[2] = nDeDxPar;
  fitAll->SetParName(2,"nDeDxPar");
  fitAll->FixParameter(2, parametersIn[2]);

  parametersIn[3] = optionDeDx;
  fitAll->SetParName(3,"optionDeDx");
  fitAll->FixParameter(3, parametersIn[3]);

  for(Int_t i = 0; i < nDeDxPar; i++) {

    parametersIn[4+i] = parDeDx[i];
    fitAll->SetParName(4+i,Form("dEdXpar%d", i));
    // if(optionDeDx == 4 && i <2) {
      
    //   fitAll->FixParameter(4+i, parametersIn[4+i]);
    // }
    if(fixAllPar) {

      fitAll->FixParameter(4+i, parametersIn[4+i]);
    }
  }

  // sigma parameters
  parametersIn[4+nDeDxPar] = nSigmaPar;
  fitAll->SetParName(4+nDeDxPar,"nSigmaPar");
  fitAll->FixParameter(4+nDeDxPar, parametersIn[4+nDeDxPar]);

  parametersIn[5+nDeDxPar] = optionSigma;
  fitAll->SetParName(5+nDeDxPar,"optionSigma");
  fitAll->FixParameter(5+nDeDxPar, parametersIn[5+nDeDxPar]);

  for(Int_t i = 0; i < nSigmaPar; i++) {

    parametersIn[6+nDeDxPar+i] = parSigma[i];
    fitAll->SetParName(6+nDeDxPar+i,Form("sigmapar%d", i));
    if(fixAllPar) {

      fitAll->FixParameter(6+nDeDxPar+i, parametersIn[6+nDeDxPar+i]);
    }
  }
  
  // Initial yields

  for(Int_t bin = binStart; bin <= binStop; bin++) {
    
    TH1D* hDeDxVsPProj =(TH1D*)hDeDxVsP->ProjectionY("hDeDxVsPProj", bin, bin);
    
    const Int_t offset = nGlobalParams + (bin - binStart)*nLocalParams;
    const Double_t all = hDeDxVsPProj->Integral();
    parametersIn[offset + 0] = (all)*0.6;
    parametersIn[offset + 1] = (all)*0.3;
    parametersIn[offset + 2] = (all)*0.1;
    
    fitAll->SetParName(offset + 0, Form("piYield%d", bin));
    fitAll->SetParLimits(offset + 0, 0, 10*all);
    fitAll->SetParName(offset + 1, Form("kYield%d", bin));
    fitAll->SetParLimits(offset + 1, 0, 10*all);
    fitAll->SetParName(offset + 2, Form("pYield%d", bin));
    fitAll->SetParLimits(offset + 2, 0, 10*all);
    // fitAll->SetParLimits(offset + 0, 0., histdEdxvsPpy->GetEntries());
    // fitAll->SetParLimits(offset + 1, 0., histdEdxvsPpy->GetEntries());
    // fitAll->SetParLimits(offset + 2, 0., histdEdxvsPpy->GetEntries());    
  }
  
  fitAll->SetParameters(parametersIn);
  fitAll->Print();

  Bool_t converged = kFALSE;

  if(performFit) {
    for(Int_t i = 0; i < 4; i++) {
      TFitResultPtr fitResultPtr =  hDeDxVsP->Fit(fitAll, "0S", "", pStart+0.01, pStop-0.01);
      if(!fitResultPtr->Status()) {
	//      if(!fitResultPtr->Status() && !fitResultPtr->CovMatrixStatus()) {
	
	converged = kTRUE;
	break;
      }
    }
  }
  // else we just draw how the results looks with the input parameters

  Double_t parametersOut[nParams];
  fitAll->GetParameters(parametersOut);
  const Double_t* parameterErrorsOut = fitAll->GetParErrors();

  // overlay the fitfunction
  

  TF1* fDeDxPi = new TF1("fDeDxPi", FitFunc, 0, 50, nDeDxPar+1); // +1 because of option! 
  fDeDxPi->SetParameters(&parametersOut[3]);
  fDeDxPi->SetLineWidth(2);
  cDeDxVsP->cd();
  fDeDxPi->Draw("SAME");
  cDeDxVsPLogX->cd();
  fDeDxPi->Draw("SAME");

  TF1* fDeDxK = new TF1("fDeDxK", FitFunc, 0, 50, nDeDxPar+1); 
  fDeDxK->SetParameters(&parametersOut[3]);
  fDeDxK->SetParameter(0, fDeDxK->GetParameter(0)+10);
  fDeDxK->SetLineWidth(2);
  cDeDxVsP->cd();
  fDeDxK->Draw("SAME");
  cDeDxVsPLogX->cd();
  fDeDxK->Draw("SAME");

  TF1* fDeDxP = new TF1("fDeDxP", FitFunc, 0, 50, nDeDxPar+1); 
  fDeDxP->SetParameters(&parametersOut[3]);
  fDeDxP->SetParameter(0, fDeDxP->GetParameter(0)+20);
  fDeDxP->SetLineWidth(2);
  cDeDxVsP->cd();
  fDeDxP->Draw("SAME");
  cDeDxVsPLogX->cd();
  fDeDxP->Draw("SAME");

  TF1* fDeDxE = new TF1("fDeDxE", FitFunc, 0, 50, nDeDxPar+1); 
  fDeDxE->SetParameters(&parametersOut[3]);
  fDeDxE->SetParameter(0, fDeDxE->GetParameter(0)+30);
  fDeDxE->SetLineWidth(2);
  cDeDxVsP->cd();
  fDeDxE->Draw("SAME");
  cDeDxVsPLogX->cd();
  fDeDxE->Draw("SAME");

  TF1* fSigma = new TF1("fSigma", SigmaFunc, 0, 50, nSigmaPar+1); 
  fSigma->SetParameters(&parametersOut[5 + nDeDxPar]);

  //  fitAll->Draw("same cont3"); 

  CreateDir("results/calibdedx");

  cDeDxVsP->cd();
  cDeDxVsP->Modified();
  cDeDxVsP->Update();
  gROOT->ProcessLine(".x drawText.C");
  cDeDxVsP->SaveAs("results/calibdedx/dedx_vs_p.gif");
  cDeDxVsP->SaveAs("results/calibdedx/dedx_vs_p.pdf");

  cDeDxVsPLogX->cd();
  cDeDxVsPLogX->Modified();
  cDeDxVsPLogX->Update();
  gROOT->ProcessLine(".x drawText.C");
  cDeDxVsPLogX->SaveAs("results/calibdedx/dedx_vs_p_logx.gif");
  cDeDxVsPLogX->SaveAs("results/calibdedx/dedx_vs_p_logx.pdf");

  //cross check
  TCanvas* cFits = new TCanvas("cFits", "Fit comparison to data", 1200, 800);
  cFits->Clear();
  cFits->Divide(7, 5);

  TF1* pion = new TF1("pion", "gausn", 30, 90);
  pion->SetLineWidth(2);
  pion->SetLineColor(kRed);
  TF1* kaon = new TF1("kaon", "gausn", 30, 90);
  kaon->SetLineWidth(2);
  kaon->SetLineColor(kGreen);
  TF1* proton = new TF1("proton", "gausn", 30, 90);
  proton->SetLineWidth(2);
  proton->SetLineColor(kBlue);
  TF1* total = new TF1("total", "gausn(0)+gausn(3)+gausn(6)", 30, 90);
  total->SetLineColor(kBlack);
  total->SetLineWidth(2);
  total->SetLineStyle(2);

  TLegend* legend = new TLegend(0.11, 0.6, 0.35, 0.85);    
  legend->SetBorderSize(0);
  legend->SetFillColor(0);
  legend->AddEntry(total, "3-Gauss fit", "L");
  legend->AddEntry(pion, "#pi", "L");
  legend->AddEntry(kaon, "K", "L");
  legend->AddEntry(proton, "p", "L");


  TCanvas* cSingleFit = new TCanvas("cSingleFit", "single fit");

  TH1D* hPionRatio =(TH1D*)hDeDxVsP->ProjectionX("hPionRatio", 1, 1);
  hPionRatio->Reset();
  hPionRatio->GetXaxis()->SetRangeUser(pStart+0.001, pStop-0.001);
  hPionRatio->GetYaxis()->SetRangeUser(0.0, 1.0);
  hPionRatio->SetTitle("particle fractions; p [GeV/c]; particle fraction");
  TH1D* hKaonRatio   = (TH1D*)hPionRatio->Clone("hKaonRatio");
  TH1D* hProtonRatio = (TH1D*)hPionRatio->Clone("hProtonRatio");

  TH1D* hPionRatioMc = (TH1D*)hPionRatio->Clone("hPionRatioMc");
  TH1D* hKaonRatioMc = (TH1D*)hPionRatio->Clone("hKaonRatioMc");
  TH1D* hProtonRatioMc = (TH1D*)hPionRatio->Clone("hProtonRatioMc");
  TH1D* hMuonRatioMc = (TH1D*)hPionRatio->Clone("hMuonRatioMc");
  
  for(Int_t bin = binStart; bin <= binStop; bin++){
    
    cout << "Making projection for bin: " << bin << endl;
    
    const Int_t j = bin-binStart;
    cFits->cd();
    cFits->cd(j + 1);
    
    TH1D* hDeDxVsPProj =(TH1D*)hDeDxVsP->ProjectionY(Form("hDeDxVsPProj%d", bin), bin, bin);
    hDeDxVsPProj->GetXaxis()->SetRangeUser(40, 90);
    hDeDxVsPProj->SetTitle(Form("%.1f<p<%.1f GeV/c", 
				hDeDxVsP->GetXaxis()->GetBinLowEdge(bin),
				hDeDxVsP->GetXaxis()->GetBinUpEdge(bin)));
    hDeDxVsPProj->Draw();
    
    const Int_t offset = nGlobalParams + j*nLocalParams; 
    const Double_t p = hDeDxVsP->GetXaxis()->GetBinCenter(bin);
    const Double_t pKeff = p*piMass/kMass; // corresponding p of a pion with same dE/dx
    const Double_t pPeff = p*piMass/pMass; // corresponding p of a pion with same dE/dx
    const Double_t meanDeDxPi = fDeDxPi->Eval(p);
    const Double_t meanDeDxK  = fDeDxPi->Eval(pKeff);
    const Double_t meanDeDxP  = fDeDxPi->Eval(pPeff);
    Double_t gausParams[9] = { 
      parametersOut[offset + 0], 
      meanDeDxPi, 
      fSigma->Eval(meanDeDxPi) ,
      parametersOut[offset + 1], 
      meanDeDxK, 
      fSigma->Eval(meanDeDxK) ,
      parametersOut[offset + 2], 
      meanDeDxP, 
      fSigma->Eval(meanDeDxP) ,
    };

    for(Int_t i = 0; i < 9; i++) 
      cout << gausParams[i] << ", ";

    cout << endl;
    
    pion->SetParameters(&gausParams[0]);
    pion->DrawCopy("same");
    Double_t all =  hDeDxVsPProj->Integral();
    hPionRatio->SetBinContent(bin, parametersOut[offset + 0]/all);
    hPionRatio->SetBinError(bin, parameterErrorsOut[offset + 0]/all);

    kaon->SetParameters(&gausParams[3]);
    kaon->DrawCopy("same");
    hKaonRatio->SetBinContent(bin, parametersOut[offset + 1]/all);
    hKaonRatio->SetBinError(bin, parameterErrorsOut[offset + 1]/all);
    
    proton->SetParameters(&gausParams[6]);
    proton->DrawCopy("same");
    hProtonRatio->SetBinContent(bin, parametersOut[offset + 2]/all);
    hProtonRatio->SetBinError(bin, parameterErrorsOut[offset + 2]/all);
    
    total->SetParameters(gausParams);
    total->DrawCopy("same");

    cSingleFit->cd();
    cSingleFit->Clear();
    //    cSingleFit->SetLogy();
    hDeDxVsPProj->Draw();
    pion->DrawCopy("same");
    kaon->DrawCopy("same");
    proton->DrawCopy("same");
    total->DrawCopy("same");
    
    gROOT->ProcessLine(".x drawText.C(2)");
    cSingleFit->SaveAs(Form("results/calibdedx/ptspectrum_bin%d.gif", bin));
    cSingleFit->SaveAs(Form("results/calibdedx/ptspectrum_bin%d.pdf", bin));
    //    legend->Draw();

    if(calib->IsMc()) {

      cSingleFit->cd();
      cSingleFit->Clear();
      TH1D* hDeDxVsPPiProj =(TH1D*)hDeDxVsPPi->ProjectionY(Form("hDeDxVsPPiProj%d", bin), bin, bin);
      hDeDxVsPPiProj->SetMarkerStyle(20);
      hDeDxVsPPiProj->SetMarkerColor(2);
      hDeDxVsPPiProj->GetXaxis()->SetRangeUser(40, 90);
      hDeDxVsPPiProj->SetTitle(Form("%.1f<p<%.1f GeV/c", 
				    hDeDxVsP->GetXaxis()->GetBinLowEdge(bin),
				    hDeDxVsP->GetXaxis()->GetBinUpEdge(bin)));
      hDeDxVsPPiProj->Draw("P");
      hPionRatioMc->SetBinContent(bin, hDeDxVsPPiProj->Integral()/all);
      TH1D* hDeDxVsPKProj =(TH1D*)hDeDxVsPK->ProjectionY(Form("hDeDxVsPKProj%d", bin), bin, bin);
      hDeDxVsPKProj->SetMarkerStyle(21);
      hDeDxVsPKProj->SetMarkerColor(3);
      hDeDxVsPKProj->Draw("SAME P");
      hKaonRatioMc->SetBinContent(bin, hDeDxVsPKProj->Integral()/all);
      TH1D* hDeDxVsPPProj =(TH1D*)hDeDxVsPP->ProjectionY(Form("hDeDxVsPPProj%d", bin), bin, bin);
      hDeDxVsPPProj->SetMarkerStyle(22);
      hDeDxVsPPProj->SetMarkerColor(4);
      hDeDxVsPPProj->Draw("SAME P");
      hProtonRatioMc->SetBinContent(bin, hDeDxVsPPProj->Integral()/all);
      TH1D* hDeDxVsPMuProj =(TH1D*)hDeDxVsPMu->ProjectionY(Form("hDeDxVsPMuProj%d", bin), bin, bin);
      hDeDxVsPMuProj->SetMarkerStyle(22);
      hDeDxVsPMuProj->SetMarkerColor(6);
      //      hDeDxVsPMuProj->Draw("SAME P");
      hMuonRatioMc->SetBinContent(bin, hDeDxVsPMuProj->Integral()/all);
      //    cSingleFit->SetLogy();
      pion->DrawCopy("same");
      kaon->DrawCopy("same");
      proton->DrawCopy("same");
      CreateDir("results/calibdedx/debugmc");
      cSingleFit->SaveAs(Form("results/calibdedx/debugmc/ptspectrum_bin%d.gif", bin));
    }

    if(hDeDxVsPV0Lambda) {
      cSingleFit->cd();
      cSingleFit->Clear();
      TH1D* hDeDxVsPV0LambdaProj =(TH1D*)hDeDxVsPV0Lambda->ProjectionY(Form("hDeDxVsPV0LambdaProj%d", bin), bin, bin);
      hDeDxVsPV0LambdaProj->SetMarkerStyle(20);
      hDeDxVsPV0LambdaProj->SetMarkerColor(4);
      hDeDxVsPV0LambdaProj->GetXaxis()->SetRangeUser(40, 90);
      hDeDxVsPV0LambdaProj->SetTitle(Form("%.1f<p<%.1f GeV/c", 
					  hDeDxVsP->GetXaxis()->GetBinLowEdge(bin),
					  hDeDxVsP->GetXaxis()->GetBinUpEdge(bin)));
      hDeDxVsPV0LambdaProj->Draw("P");
      proton->SetParameter(0, hDeDxVsPV0LambdaProj->Integral());
      proton->DrawCopy("same");

      CreateDir("results/calibdedx/debuglambda");
      cSingleFit->SaveAs(Form("results/calibdedx/debuglambda/ptspectrum_bin%d.gif", bin));
    }

    if(hDeDxVsPV0Kaon) {
      cSingleFit->cd();
      cSingleFit->Clear();
      TH1D* hDeDxVsPV0KaonProj =(TH1D*)hDeDxVsPV0Kaon->ProjectionY(Form("hDeDxVsPV0KaonProj%d", bin), bin, bin);
      hDeDxVsPV0KaonProj->SetMarkerStyle(20);
      hDeDxVsPV0KaonProj->SetMarkerColor(2);
      hDeDxVsPV0KaonProj->GetXaxis()->SetRangeUser(40, 90);
      hDeDxVsPV0KaonProj->SetTitle(Form("%.1f<p<%.1f GeV/c", 
					  hDeDxVsP->GetXaxis()->GetBinLowEdge(bin),
					  hDeDxVsP->GetXaxis()->GetBinUpEdge(bin)));
      hDeDxVsPV0KaonProj->Draw("P");
      pion->SetParameter(0, hDeDxVsPV0KaonProj->Integral());
      pion->DrawCopy("same");
      CreateDir("results/calibdedx/debugkaon");
      cSingleFit->SaveAs(Form("results/calibdedx/debugkaon/ptspectrum_bin%d.gif", bin));
    }
  }

  TCanvas* cRatio = new TCanvas("cRatio", "ratios/all vs p", 600, 400);
  cRatio->Clear();
  hPionRatio->SetMaximum(0.8);
  hPionRatio->SetMarkerStyle(20);
  hPionRatio->SetMarkerColor(2);
  hPionRatio->Draw("P E");

  hKaonRatio->SetMarkerStyle(20);
  hKaonRatio->SetMarkerColor(3);
  hKaonRatio->Draw("SAME P E");

  hProtonRatio->SetMarkerStyle(20);
  hProtonRatio->SetMarkerColor(4);
  hProtonRatio->Draw("SAME P E");
  gROOT->ProcessLine(".x drawText.C(2)");
  cRatio->SaveAs("results/calibdedx/particle_ratios.gif");
  cRatio->SaveAs("results/calibdedx/particle_ratios.pdf");

  if(calib->IsMc()) {
    
    hPionRatioMc->SetMarkerStyle(24);
    hPionRatioMc->SetMarkerColor(2);
    hPionRatioMc->Draw("SAME P");
    
    hKaonRatioMc->SetMarkerStyle(24);
    hKaonRatioMc->SetMarkerColor(3);
    hKaonRatioMc->Draw("SAME P");
    
    hProtonRatioMc->SetMarkerStyle(24);
    hProtonRatioMc->SetMarkerColor(4);
    hProtonRatioMc->Draw("SAME P");

    hMuonRatioMc->SetMarkerStyle(24);
    hMuonRatioMc->SetMarkerColor(6);
    //    hMuonRatioMc->Draw("SAME P");
    cRatio->SaveAs("results/calibdedx/debugmc/particle_ratios_mc.gif");
  }


  cout << "MIP was fixed to: " << fixMIP << endl
       << "Plateau was fixed to: " << fixPlateau << endl;

  //
  // Store the <dE/dx> parameters in a file that we can get them back to use for the Delta-pi!
  //
  DeDxFitInfo* fitInfo = new DeDxFitInfo();
  fitInfo->MIP     = fixMIP;
  fitInfo->plateau = fixPlateau; 
  fitInfo->optionDeDx = optionDeDx; 
  fitInfo->nDeDxPar = nDeDxPar; 
  for(Int_t i = 0; i < nDeDxPar; i++) {
    fitInfo->parDeDx[i] = fDeDxPi->GetParameter(i+1); // 1st parameter is option
  }
  fitInfo->optionSigma = optionSigma; 
  fitInfo->nSigmaPar = nSigmaPar; 
  for(Int_t i = 0; i < nSigmaPar; i++) {
    fitInfo->parSigma[i] = fSigma->GetParameter(i+1); // 1st parameter is option 
  }
  if(!strstr(calibFileName, "no_calib_eta"))
    fitInfo->calibFileName = calibFileName;

  fitInfo->Print();

  CreateDir("fitparameters");
  TFile* outFile = new TFile(Form("fitparameters/%s", gSystem->BaseName(calibFileName)), "RECREATE");
  fitInfo->Write("fitInfo");
  outFile->Close();

  if(converged) {

    cout << "Fit converged and error matrix was ok" << endl;
  } else {

    cout << "WARNING!!!!!!!!!!!!!!! Fit did not converge, or error matrix was  not ok!" << endl;
  }
}



//____________________________________________________________________________
void DrawStep2(const Char_t* calibFileName,
	       Bool_t forMIP = kTRUE,
	       Int_t run    = 0,
	       Int_t filter = 1,
	       Bool_t usePhiCut = kTRUE)
{
  // option:
  // 0 - compare to ALICE INEL
  // 1 - compare to MC input spectrum
  // 2 - compare to MC correct spectrum
  gStyle->SetOptStat(0);

  CreateDir("results/calibplots");
  TString baseName(gSystem->BaseName(calibFileName));
  baseName.ReplaceAll(".root", "");

  TFile* calibFile = FindFileFresh(calibFileName);
  if(!calibFile)
    return;
  AliHighPtDeDxCalib* calib = (AliHighPtDeDxCalib*)GetObject(calibFile, filter, usePhiCut, run);
  calib->Print();
 
  TCanvas* cPhiCut = calib->DrawPhiCutHistograms();
  gROOT->ProcessLine(".x drawText.C(2)");
  cPhiCut->SaveAs(Form("results/calibplots/phicut_%s.gif", baseName.Data())); 
  cPhiCut->SaveAs(Form("results/calibplots/phicut_%s.pdf", baseName.Data())); 
  TCanvas* cEta = calib->DrawEta(forMIP);
  gROOT->ProcessLine(".x drawText.C(2)");
  TCanvas* cEtaCal = calib->DrawEtaCalibrated(forMIP);
  gROOT->ProcessLine(".x drawText.C(2)");
  if(forMIP) {
    cEta->cd();
    cEta->SaveAs(Form("results/calibplots/MIP_dedx_vs_eta_nocal_%s.gif", baseName.Data())); 
    cEta->SaveAs(Form("results/calibplots/MIP_dedx_vs_eta_nocal_%s.pdf", baseName.Data())); 
    cEtaCal->cd();
    cEtaCal->SaveAs(Form("results/calibplots/MIP_dedx_vs_eta_final_%s.gif", baseName.Data())); 
    cEtaCal->SaveAs(Form("results/calibplots/MIP_dedx_vs_eta_final_%s.pdf", baseName.Data())); 
  } else {
    cEta->cd();
    cEta->SaveAs(Form("results/calibplots/Electrons_dedx_vs_eta_nocal_%s.gif", baseName.Data())); 
    cEta->SaveAs(Form("results/calibplots/Electrons_dedx_vs_eta_nocal_%s.pdf", baseName.Data())); 
    cEtaCal->cd();
    cEtaCal->SaveAs(Form("results/calibplots/Electrons_dedx_vs_eta_final_%s.gif", baseName.Data())); 
    cEtaCal->SaveAs(Form("results/calibplots/Electrons_dedx_vs_eta_final_%s.pdf", baseName.Data())); 
  }
  TCanvas* cEtaSel = calib->DrawSelectionHistograms();
  gROOT->ProcessLine(".x drawText.C(2)");
  cEtaSel->SaveAs(Form("results/calibplots/selection_%s.gif", baseName.Data())); 
  cEtaSel->SaveAs(Form("results/calibplots/selection_%s.pdf", baseName.Data())); 

  // TCanvas* cNcl = calib->DrawNclCal();
  // gROOT->ProcessLine(".x drawText.C");
  // cNcl->SaveAs(Form("results/calibplots/dedx_vs_ncl_calib__%s.gif", baseName.Data())); 
  // cNcl->SaveAs(Form("results/calibplots/dedx_vs_ncl_calib__%s.pdf", baseName.Data())); 
}

//____________________________________________________________________________
void TestResolutionVsDeDx(const Char_t* calibFileName,
			  Int_t run    = 0,
			  Int_t filter = 1,
			  Bool_t usePhiCut = kTRUE)
{
  gStyle->SetOptStat(0);
  
  TFile* calibFile = FindFileFresh(calibFileName);
  if(!calibFile)
    return;
  AliHighPtDeDxCalib* calib = (AliHighPtDeDxCalib*)GetObject(calibFile, filter, usePhiCut, run);
  calib->Print();
 
  TH2D* hDeDxVsNcl  = calib->GetHistDeDxVsNcl(kTRUE, 0);
  TH2D* hDeDxVsNclE = calib->GetHistDeDxVsNcl(kFALSE, 0);
  TH1D* hDeDx       = calib->GetHistDeDx(kTRUE, 0);
  TH1D* hDeDxE      = calib->GetHistDeDx(kFALSE, 0);

  Double_t mdedx  = hDeDx->GetMean();
  Double_t mdedxE = hDeDxE->GetMean();

  TObjArray* helpArray = new TObjArray();

  hDeDxVsNcl->FitSlicesY(0, 0, -1, 0, "QNR", helpArray);
  TH1D* hMeanVsNcl =  (TH1D*)helpArray->At(1);
  hMeanVsNcl->SetName("hMeanVsNcl");
  hMeanVsNcl->SetMarkerStyle(29);
  TH1D* hSigmaVsNcl =  (TH1D*)helpArray->At(2);
  hSigmaVsNcl->SetName("hSigmaVsNcl");
  hSigmaVsNcl->SetTitle("#sigma vs Ncl; Ncl; #sigma");
  hSigmaVsNcl->SetMarkerStyle(29);
  hSigmaVsNcl->SetMarkerColor(2);
  hSigmaVsNcl->Scale(1.0/mdedx);

  hDeDxVsNclE->FitSlicesY(0, 0, -1, 0, "QNR", helpArray);
  TH1D* hMeanVsNclE =  (TH1D*)helpArray->At(1);
  hMeanVsNclE->SetName("hSigmaVsNclE");
  hMeanVsNclE->SetMarkerStyle(29);
  TH1D* hSigmaVsNclE =  (TH1D*)helpArray->At(2);
  hSigmaVsNclE->SetName("hSigmaVsNclE");
  hSigmaVsNclE->SetTitle("#sigma vs Ncl; Ncl; #sigma");
  hSigmaVsNclE->SetMarkerStyle(29);
  hSigmaVsNclE->SetMarkerColor(kMagenta);
  hSigmaVsNclE->Scale(1.0/mdedxE);

  TLegend* legend = new TLegend(0.54, 0.67, 0.84, 0.87);
  legend->SetFillColor(0);
  legend->SetBorderSize(0);
  legend->SetTextSize(0.05);
  legend->AddEntry(hSigmaVsNcl, "MIP pions", "P");
  legend->AddEntry(hSigmaVsNclE, "electrons", "P");

  TCanvas* cTestDeDx = new TCanvas("cTestDeDx", "Comparing resolution for MIPs and electrons", 
				   1200, 400);
  cTestDeDx->Clear();
  cTestDeDx->Divide(3, 1);
  
  cTestDeDx->cd(1);
  hDeDxVsNcl->Draw("COL");
  hMeanVsNcl->Draw("SAME");

  cTestDeDx->cd(2);
  hDeDxVsNclE->Draw("COL");
  hMeanVsNclE->Draw("SAME");

  cTestDeDx->cd(3);
  TF1* fSDedxVsNclSqrt = new TF1("fSDedxVsNclSqrt", "[0]*sqrt(159.0/x)",  0, 160);
  fSDedxVsNclSqrt->SetParameter(0, 0.05);
  hSigmaVsNcl->SetMinimum(0.0);
  hSigmaVsNcl->SetMaximum(0.15);
  hSigmaVsNcl->SetTitle("#sigma/<dE/dx> vs Ncl; Ncl; rel. #sigma");
  //  hSigmaVsNcl->Fit(fSDedxVsNclSqrt);
  hSigmaVsNcl->Draw();
  hSigmaVsNclE->Draw("SAME");
  legend->Draw();
  gROOT->ProcessLine(".x drawText.C(2)");
  TString baseName(gSystem->BaseName(calibFileName));
  baseName.ReplaceAll(".root", "");
  cTestDeDx->SaveAs(Form("results/calibplots/test_resolution_vs_dedx_%s.gif", baseName.Data())); 
  cTestDeDx->SaveAs(Form("results/calibplots/test_resolution_vs_dedx_%s.pdf", baseName.Data()));   
}  

//___________________________________________________________________________________________
void TestResolutionVsEta(const Char_t* calibFileName,
			 Bool_t forMIP,
			 Int_t run    = 0,
			 Int_t filter = 1,
			 Bool_t usePhiCut = kTRUE)
{
  gStyle->SetOptStat(0);
  
  TFile* calibFile = FindFileFresh(calibFileName);
  if(!calibFile)
    return;
  AliHighPtDeDxCalib* calib = (AliHighPtDeDxCalib*)GetObject(calibFile, filter, usePhiCut, run);
  calib->Print();
 
  const Int_t nEta = 4;
  Int_t color[nEta] = {kBlack, kBlue, kGreen+2, kRed};
  TH1D* hDeDx[nEta] = {0, 0, 0, 0};
  TH1D* hMeanVsNcl[nEta] = {0, 0, 0, 0};
  TH1D* hSigmaVsNcl[nEta] = {0, 0, 0, 0};

  const Double_t mdedx  = (calib->GetHistDeDx(forMIP, 0))->GetMean();
  
  TObjArray* helpArray = new TObjArray();
  
  for(Int_t bin = 1; bin <= nEta; bin++) {
    
    const Int_t i = bin - 1;
    hDeDx[i]  = calib->GetHistDeDx(forMIP, bin);
    hDeDx[i]->Scale(1.0/hDeDx[i]->GetEntries());
    hDeDx[i]->SetMarkerStyle(24+i);
    hDeDx[i]->SetMarkerColor(color[i]);

    TH2D* hDeDxVsNcl  = calib->GetHistDeDxVsNcl(forMIP, bin);
    hDeDxVsNcl->FitSlicesY(0, 0, -1, 0, "QNR", helpArray);

    hMeanVsNcl[i] =  (TH1D*)helpArray->At(1);
    hMeanVsNcl[i]->SetTitle("<dE/dx> vs Ncl; Ncl; #sigma");
    hMeanVsNcl[i]->SetName(Form("hMeanVsNcl%d", bin));
    hMeanVsNcl[i]->SetMarkerStyle(24+i);
    hMeanVsNcl[i]->SetMarkerColor(color[i]);
    hMeanVsNcl[i]->SetMinimum(mdedx*0.9);
    hMeanVsNcl[i]->SetMaximum(mdedx*1.1);

    hSigmaVsNcl[i] =  (TH1D*)helpArray->At(2);
    hSigmaVsNcl[i]->SetName("hSigmaVsNcl");
    hSigmaVsNcl[i]->SetTitle("#sigma/<dE/dx> vs Ncl; Ncl; rel. #sigma");
    hSigmaVsNcl[i]->SetMarkerStyle(24+i);
    hSigmaVsNcl[i]->SetMarkerColor(color[i]);
    hSigmaVsNcl[i]->Scale(1.0/mdedx);
    hSigmaVsNcl[i]->SetMinimum(0.0);
    hSigmaVsNcl[i]->SetMaximum(0.15);
  }

  TLegend* legend = new TLegend(0.54, 0.67, 0.84, 0.87);
  legend->SetFillColor(0);
  legend->SetBorderSize(0);
  legend->SetTextSize(0.05);
  
  legend->AddEntry(hDeDx[0], "|#eta|<0.2", "P");
  legend->AddEntry(hDeDx[1], "0.2<|#eta|<0.4", "P");
  legend->AddEntry(hDeDx[2], "0.4<|#eta|<0.6", "P");
  legend->AddEntry(hDeDx[3], "0.6<|#eta|<0.8", "P");

  TCanvas* cTestEta = new TCanvas("cTestEta", "Comparing resolution for MIPs and electrons", 
				   1200, 400);
  cTestEta->Clear();
  cTestEta->Divide(3, 1);
  
  for(Int_t i = 0; i < nEta; i++) {

    cTestEta->cd(1);
    if(i==0)
      hDeDx[i]->Draw("P");
    else
      hDeDx[i]->Draw("SAME P");
      	
    cTestEta->cd(2);
    if(i==0)
      hMeanVsNcl[i]->Draw("P");
    else
      hMeanVsNcl[i]->Draw("SAME P");

    cTestEta->cd(3);
    if(i==0) {
      hSigmaVsNcl[i]->Draw("P");
      legend->Draw();
    } else
      hSigmaVsNcl[i]->Draw("SAME P");
  }

  gROOT->ProcessLine(".x drawText.C(2)");
  TString baseName(gSystem->BaseName(calibFileName));
  baseName.ReplaceAll(".root", "");
  cTestEta->SaveAs(Form("results/calibplots/test_resolution_vs_eta_%s.gif", baseName.Data())); 
  cTestEta->SaveAs(Form("results/calibplots/test_resolution_vs_eta_%s.pdf", baseName.Data()));   
}  
  
  // hMeanVsNcl =  (TH1D*)helpArray->At(1);
  // hMeanVsNcl->SetName("hMeanVsNcl");
  // hMeanVsNcl->SetTitle("Mean vs Ncl; Ncl; Mean");
  // hMeanVsNcl->SetMarkerStyle(29);
  // hMeanVsNcl->SetMinimum(45);
  // hMeanVsNcl->SetMaximum(55);
  
  // if(isMC) {
  // 	hDeDxVsNclElectrons->FitSlicesY(0, 0, -1, 0, "QNR", helpArray);
  // 	hSigmaVsNclElectrons =  (TH1D*)helpArray->At(2);
  // 	hSigmaVsNclElectrons->SetName("hSigmaVsNcl");
  // 	hSigmaVsNclElectrons->SetTitle("#sigma vs Ncl; Ncl; #sigma");
  //     }

  //     TCanvas* cNcl2 = new TCanvas("cNcl2", "sigma dE/dx vs Ncl", 600, 400);
  //     cNcl2->Clear();
  //     cNcl2->cd();
      
  //     hSigmaVsNcl->DrawCopy();
  //     hSigmaVsNcl->Fit(fSDeDxVsNcl, "0");
  //     hSigmaVsNcl->Fit(fSDeDxVsNclSqrt, "0");
      
  //     fSDeDxVsNcl->DrawCopy("SAME");
  //     fSDeDxVsNclSqrt->SetLineStyle(2);
  //     fSDeDxVsNclSqrt->DrawCopy("SAME");

  //     TLegend* legFit = new TLegend(0.64, 0.67, 0.84, 0.87);
  //     legFit->SetFillColor(0);
  //     legFit->SetBorderSize(0);
  //     legFit->SetTextSize(0.05);
  //     legFit->AddEntry(fSDeDxVsNcl, "Linear fit", "L");
  //     legFit->AddEntry(fSDeDxVsNclSqrt, "1/Sqrt fit", "L");
  //     legFit->Draw();

  //     cNcl2->SaveAs(Form("%s/sigma_vs_ncl.gif", baseName.Data()));

  //     if(isMC) {
  // 	TCanvas* cNcl2Electrons = new TCanvas("cNcl2Electrons", "sigma dE/dx vs Ncl (comparison)", 600, 400);
  // 	cNcl2Electrons->Clear();
  // 	cNcl2Electrons->cd();
	
  // 	hSigmaVsNcl->SetMarkerStyle(29);
  // 	hSigmaVsNcl->DrawCopy();
  // 	hSigmaVsNclElectrons->Scale(50.0 / hDeDxVsEtaCalibratedElectrons->GetMean(2));
  // 	hSigmaVsNclElectrons->SetMarkerStyle(29);
  // 	hSigmaVsNclElectrons->SetMarkerColor(6);
  // 	hSigmaVsNclElectrons->DrawCopy("SAME");
  // 	fSDeDxVsNcl->DrawCopy("SAME");
  // 	fSDeDxVsNclSqrt->DrawCopy("SAME");
	
  // 	TLegend* legSigma = new TLegend(0.64, 0.67, 0.84, 0.87);
  // 	legSigma->SetFillColor(0);
  // 	legSigma->SetBorderSize(0);
  // 	legSigma->SetTextSize(0.05);
  // 	legSigma->AddEntry(hSigmaVsNcl, "MIP pions", "P");
  // 	legSigma->AddEntry(fSDeDxVsNcl, "Linear fit", "L");
  // 	legSigma->AddEntry(fSDeDxVsNclSqrt, "1/Sqrt fit", "L");
  // 	legSigma->AddEntry(hSigmaVsNclElectrons, "Scaled electrons", "P");
  // 	legSigma->Draw();
	
  // 	cNcl2Electrons->SaveAs(Form("%s/electrons_sigma_vs_ncl_comparison.gif", baseName.Data()));
  //     }

  //     TCanvas* cDeDx = new TCanvas("cDeDx", "dE/dx", 600, 400);
  //     cDeDx->Clear();
  //     hDeDx->Scale(1.0/hDeDx->GetEntries());
  //     hDeDx->SetMarkerStyle(29);
  //     hDeDx->DrawCopy();
  //     TCanvas* cNcl3 = new TCanvas("cNcl3", "#sigma vs Ncl", 600, 400);
  //     cNcl3->Clear();
  //     hSigmaVsNcl->DrawCopy();
  //     TCanvas* cNcl3Mean = new TCanvas("cNcl3Mean", "Mean vs Ncl", 600, 400);
  //     cNcl3Mean->Clear();
  //     hMeanVsNcl->DrawCopy();

  //     TLegend* legDeDx = new TLegend(0.7, 0.67, 0.9, 0.87);
  //     legDeDx->SetFillColor(0);
  //     legDeDx->SetBorderSize(0);
  //     legDeDx->SetTextSize(0.05);
  //     legDeDx->AddEntry(hDeDx, "No #eta cut", "P");

  //     TLegend* legScan = new TLegend(0.44, 0.67, 0.64, 0.87);
  //     legScan->SetFillColor(0);
  //     legScan->SetBorderSize(0);
  //     legScan->SetTextSize(0.05);
  //     legScan->AddEntry(hSigmaVsNcl, "No #eta cut", "P");

  //     TH2D* hArray[4] = {hDeDxVsNcl0, hDeDxVsNcl1, hDeDxVsNcl2, hDeDxVsNcl3};
  //     TH1D* hArrayDeDx[4] = {hDeDx0, hDeDx1, hDeDx3, hDeDx4};

  //     Int_t color[4] = {2, kGreen+1, 4, 6};
  //     Double_t x[4];
  //     Double_t x_err[4];
  //     Double_t slope[4];
  //     Double_t slope_err[4];
  //     Double_t offset[4];
  //     Double_t offset_err[4];
  //     for(Int_t i = 0; i < 4; i++) {

  // 	hArrayDeDx[i]->SetMarkerStyle(20);
  // 	hArrayDeDx[i]->SetMarkerColor(color[i]);
  // 	hArrayDeDx[i]->Scale(1.0/hArrayDeDx[i]->GetEntries());
  // 	cDeDx->cd();
  // 	hArrayDeDx[i]->DrawCopy("SAME");

  // 	x[i] = hMeanEta->GetBinContent(i+1);
  // 	x_err[i] = hMeanEta->GetBinError(i+1);

  // 	hArray[i]->FitSlicesY(0, 0, -1, 0, "QNR", helpArray);
  // 	hSigmaVsNclHelp =  (TH1D*)helpArray->At(2);
  // 	hSigmaVsNclHelp->SetName(Form("hSigmaVsNcl%d", i));
  // 	hSigmaVsNclHelp->SetTitle("#sigma vs Ncl; Ncl; #sigma");
  // 	hSigmaVsNclHelp->SetMarkerStyle(20);
  // 	hSigmaVsNclHelp->SetMarkerColor(color[i]);
  // 	hSigmaVsNclHelp->Fit(fSDeDxVsNclHelp, "0");
  // 	slope[i] = fSDeDxVsNclHelp->GetParameter(0);
  // 	slope_err[i] = fSDeDxVsNclHelp->GetParError(0);
  // 	offset[i] = fSDeDxVsNclHelp->GetParameter(1);
  // 	offset_err[i] = fSDeDxVsNclHelp->GetParError(1);
  // 	cNcl3->cd();
  // 	hSigmaVsNclHelp->DrawCopy("SAME");
  // 	//	fSDeDxVsNclHelp->DrawCopy("SAME");

  // 	if(i==0) {
  // 	  legDeDx->AddEntry(hArrayDeDx[i], "|#eta|<0.2", "P");
  // 	  legScan->AddEntry(hSigmaVsNclHelp, "|#eta|<0.2", "P");
  // 	} else if(i==1) {
  // 	  legDeDx->AddEntry(hArrayDeDx[i], "0.2<|#eta|<0.4", "P");
  // 	  legScan->AddEntry(hSigmaVsNclHelp, "0.2<|#eta|<0.4", "P");
  // 	} else if(i==2) {
  // 	  legDeDx->AddEntry(hArrayDeDx[i], "0.4<|#eta|<0.6", "P");
  // 	  legScan->AddEntry(hSigmaVsNclHelp, "0.4<|#eta|<0.6", "P");
  // 	} else if(i==3) {
  // 	  legDeDx->AddEntry(hArrayDeDx[i], "0.6<|#eta|<0.8", "P");
  // 	  legScan->AddEntry(hSigmaVsNclHelp, "0.6<|#eta|<0.8", "P");
  // 	}

  // 	cNcl3Mean->cd();
  // 	hMeanVsNclHelp =  (TH1D*)helpArray->At(1);
  // 	hMeanVsNclHelp->SetName(Form("hMeanVsNcl%d", i));
  // 	hMeanVsNclHelp->SetTitle("Mean vs Ncl; Ncl; Mean");
  // 	hMeanVsNclHelp->SetMarkerStyle(20);
  // 	hMeanVsNclHelp->SetMarkerColor(color[i]);
  // 	hMeanVsNclHelp->DrawCopy("SAME");
  //     }

  //     cDeDx->cd();
  //     hDeDx->DrawCopy("SAME");
  //     legDeDx->Draw();
  //     cDeDx->SaveAs(Form("%s/dedx_resolution_vs_eta.gif", baseName.Data()));

  //     cNcl3->cd();
  //     hSigmaVsNcl->DrawCopy("SAME");
  //     legScan->Draw();
  //     cNcl3->SaveAs(Form("%s/sigma_vs_ncl_vs_eta.gif", baseName.Data()));

  //     cNcl3Mean->cd();
  //     hMeanVsNcl->DrawCopy("SAME");
  //     legScan->Draw();
  //     cNcl3Mean->SaveAs(Form("%s/mean_vs_ncl_vs_eta.gif", baseName.Data()));

  //     TGraphErrors* graphSlope = new TGraphErrors(4, x, slope, x_err, slope_err);
  //     graphSlope->SetTitle("slope vs <#eta>; <#eta>; slope");
  //     graphSlope->SetMarkerStyle(29);
  //     TCanvas* cNcl4 = new TCanvas("cNcl4", "Slope vs <eta>", 600, 400);
  //     cNcl4->Clear();
  //     graphSlope->Draw("AP");
  //     cNcl4->SaveAs(Form("%s/fit_slope_vs_eta.gif", baseName.Data()));

  //     TGraphErrors* graphOffset = new TGraphErrors(4, x, offset, x_err, offset_err);
  //     graphOffset->SetMarkerStyle(29);
  //     graphOffset->SetTitle("offset vs <#eta>; <#eta>; offset");
  //     TCanvas* cNcl5 = new TCanvas("cNcl5", "Offset vs <eta>", 600, 400);
  //     cNcl5->Clear();
  //     graphOffset->Draw("AP");
  //     cNcl5->SaveAs(Form("%s/fit_offset_vs_eta.gif", baseName.Data()));
      
  //   }

//____________________________________________________________________________
void CompareYields(const Char_t* dataFileName1,
		   const Char_t* dataFileName2,
		   Double_t ptStart, Double_t ptStop,
		   const Char_t* endName1=0,
		   const Char_t* endName2=0,
		   Int_t run    = 0,
		   Int_t filter = 1,
		   Bool_t usePhiCut = kTRUE)
{
  gStyle->SetOptStat(0);

  
  TFile* dataFile1 = FindFileFresh(dataFileName1);
  if(!dataFile1)
    return;
  AliHighPtDeDxCalib* data1 = (AliHighPtDeDxCalib*)GetObject(dataFile1, filter, usePhiCut, run, "filter", endName1);
  data1->Print();

  //  gSystem->Exec("mv debug/* olddebug/");


  TH2D* hDeltaPiVsPt1 = data1->GetHistDeDxVsP(0);

  AliHighPtDeDxCalib* data2 = data1;
  if(dataFileName2) {
    TFile* dataFile2 = FindFileFresh(dataFileName2);
    if(!dataFile2)
      return;
    data2 = (AliHighPtDeDxCalib*)GetObject(dataFile2, filter, usePhiCut, run, "filter", endName2);
    data2->Print();
  } else if (endName2) {

    data2 = (AliHighPtDeDxCalib*)GetObject(dataFile1, filter, usePhiCut, run, "filter", endName2);
  }

  TH2D* hDeltaPiVsPt2 = data2->GetHistDeDxVsP(0);

  // Root is a bit stupid with finidng bins so we have to add and subtract a
  // little to be sure we get the right bin as we typically put edges as
  // limits
  const Int_t binStart = hDeltaPiVsPt1->GetXaxis()->FindBin(ptStart+0.01);
  ptStart = hDeltaPiVsPt1->GetXaxis()->GetBinLowEdge(binStart);
  const Int_t binStop  = hDeltaPiVsPt1->GetXaxis()->FindBin(ptStop-0.01);
  ptStop = hDeltaPiVsPt1->GetXaxis()->GetBinUpEdge(binStop);
  //  const Int_t nBins    = binStop - binStart + 1;

  cout << "Doing fits from pTlow = " << ptStart << " (bin: " << binStart
       << ") to pThigh = " << ptStop << " (bin: " << binStop << ")" << endl;
  

  //cross check
  TCanvas* cFits = new TCanvas("cFits", "Fit comparison to data", 1200, 800);
  cFits->Clear();
  cFits->Divide(7, 4);

  TLegend* legend = new TLegend(0.11, 0.6, 0.35, 0.85);    
  legend->SetBorderSize(0);
  legend->SetFillColor(0);

  TCanvas* cSingleFit = new TCanvas("cSingleFit", "single fit");

  for(Int_t bin = binStart; bin <= binStop; bin++){
    
    cout << "Making projection for bin: " << bin << endl;
    
    const Int_t j = bin-binStart;
    
    TH1D* hDeltaPiVsPtProj1 =(TH1D*)hDeltaPiVsPt1->ProjectionY(Form("hDeltaPiVsPtProj1%d", bin), bin, bin);
    //    hDeltaPiVsPtProj1->GetXaxis()->SetRangeUser(-20, 20);
    hDeltaPiVsPtProj1->SetTitle(Form("%.1f<p<%.1f GeV/c", 
  				hDeltaPiVsPt1->GetXaxis()->GetBinLowEdge(bin),
  				hDeltaPiVsPt1->GetXaxis()->GetBinUpEdge(bin)));

    TH1D* hDeltaPiVsPtProj2 =(TH1D*)hDeltaPiVsPt2->ProjectionY(Form("hDeltaPiVsPtProj2%d", bin), bin, bin);

    hDeltaPiVsPtProj1->SetLineColor(4);
    hDeltaPiVsPtProj2->SetLineColor(2);
    hDeltaPiVsPtProj1->SetNormFactor();
    hDeltaPiVsPtProj2->SetNormFactor();

    cFits->cd();
    cFits->cd(j + 1);
    hDeltaPiVsPtProj1->DrawCopy();
    hDeltaPiVsPtProj2->DrawCopy("SAME");

    cSingleFit->cd();
    cSingleFit->Clear();
    hDeltaPiVsPtProj1->DrawCopy();
    hDeltaPiVsPtProj2->DrawCopy("SAME");

    CreateDir("results/debug");
    cSingleFit->SaveAs(Form("results/debug/ptspectrum_bin%d.gif", bin));
  }
}
