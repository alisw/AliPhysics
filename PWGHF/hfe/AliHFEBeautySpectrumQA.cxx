
/**************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
*                                                                        *
* Author: The ALICE Off-line Project.                                    *
* Contributors are mentioned in the code where appropriate.              *
*                                                                        *
* Permission to use, copy, modify and distribute this software and its   *
* documentation strictly for non-commercial purposes is hereby granted   *
* without fee, provided that the above copyright notice appears in all   *
* copies and that both the copyright notice and this permission notice   *
* appear in the supporting documentation. The authors make no claims     *
* about the suitability of this software for any purpose. It is          *
* provided "as is" without express or implied warranty.                  *
**************************************************************************/
//
// Class for spectrum correction
// Subtraction of hadronic background, Unfolding of the data and
// Renormalization done here
// The following containers have to be set:
//  - Correction framework container for real data
//  - Correction framework container for MC (Efficiency Map)
//  - Correction framework container for background coming from data
//  - Correction framework container for background coming from MC
//
//  Author: 
//            Raphaelle Bailhache <R.Bailhache@gsi.de>
//

#include <TArrayD.h>
#include <TH1.h>
#include <TList.h>
#include <TObjArray.h>
#include <TObject.h>
#include <TROOT.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TMath.h>
#include <TAxis.h>
#include <TGraphErrors.h>
#include <TFile.h>
#include <TPad.h>
#include <TH2D.h>
#include <TF1.h>

#include "AliPID.h"
#include "AliCFContainer.h"
#include "AliCFDataGrid.h"
#include "AliCFEffGrid.h"
#include "AliCFGridSparse.h"
#include "AliCFUnfolding.h"
#include "AliLog.h"

#include "AliHFEBeautySpectrumQA.h"
#include "AliHFECorrectSpectrumBase.h"
#include "AliHFEcuts.h"
#include "AliHFEcontainer.h"
#include "AliHFEtools.h"

ClassImp(AliHFEBeautySpectrumQA)

const Char_t *AliHFEBeautySpectrumQA::fgkNameCanvas[AliHFEBeautySpectrumQA::kNTypeEfficiency] = {
  "MCEfficiency",
  "ParametrizedEfficiency"
};

//____________________________________________________________
AliHFEBeautySpectrumQA::AliHFEBeautySpectrumQA():
  TNamed(),
  fPtMax(8.0),
  fListOfResult(),
  fWriteToFile(kTRUE)
  {
  //
  // Default constructor
  //

  fListOfResult = new TObjArray(kNResults);
  fListOfResult->SetName("ListOfResults");
 

}
//____________________________________________________________
AliHFEBeautySpectrumQA::AliHFEBeautySpectrumQA(const char *name):
  TNamed(name, ""),
  fPtMax(8.0),
  fListOfResult(),
  fWriteToFile(kTRUE)
  {
  //
  // Default constructor
  //

  fListOfResult = new TObjArray(kNResults);
  fListOfResult->SetName("ListOfResults");
 

}

//____________________________________________________________
AliHFEBeautySpectrumQA::~AliHFEBeautySpectrumQA(){
  //
  // Destructor
  //
  if(fListOfResult) delete fListOfResult;
 
}
//____________________________________________________________
void AliHFEBeautySpectrumQA::AddResultAt(TObject *obj,Int_t index)
{
  //
  // Init what we need for the correction:
  //

  if(fListOfResult) fListOfResult->AddAt(obj,index);

}
//____________________________________________________________
TObject *AliHFEBeautySpectrumQA::GetResult(Int_t index)
{
  //
  // Get result
  //

  if(fListOfResult) return fListOfResult->UncheckedAt(index);
  else return 0x0;

}
//____________________________________________________________
void AliHFEBeautySpectrumQA::DrawProjections() const
{
  //
  // get spectrum for beauty 2nd method
  //
  //
  AliCFContainer *data = (AliCFContainer *) fListOfResult->UncheckedAt(kDataProjection);
  THnSparseF *correlation = (THnSparseF *) fListOfResult->UncheckedAt(kCMProjection);
  if(!data || !correlation) return;

  Int_t ndimcont = data->GetNVar();
  Int_t ndimcor = correlation->GetNdimensions();
  Int_t charge = 3;
  Int_t centrality = 5;
  Int_t eta = 1;

  TCanvas * canvas = new TCanvas("Projections","Projections",1000,700);
  Int_t n = 0;
  if(charge < ndimcont) n++;
  if(centrality < ndimcont) n++;
  if(eta < ndimcont) n++;
  canvas->Divide(2,n);
  Int_t counter = 1;

  if(charge < ndimcont) {
   
    canvas->cd(counter);
    TH1 *checkcharge = (TH1 *) data->Project(data->GetNStep()-1,charge);
    checkcharge->Draw();
    counter++;
    canvas->cd(counter);
    TH2F* projectioncharge = (TH2F *) correlation->Projection(charge,charge+((Int_t)(ndimcor/2.)));
    projectioncharge->Draw("colz");
    counter++;  

  }

  if(centrality < ndimcont) {
    canvas->cd(counter);
    TH1 *checkcentrality = (TH1 *) data->Project(data->GetNStep()-1,centrality);
    checkcentrality->Draw();
    counter++;
    canvas->cd(counter);
    TH2F *projectioncentrality = (TH2F *) correlation->Projection(centrality,centrality+((Int_t)(ndimcor/2.)));
    projectioncentrality->Draw("colz");
    counter++;  
  }

  if(eta < ndimcont) {
    canvas->cd(counter);
    TH1 *checketa = (TH1 *) data->Project(data->GetNStep()-1,eta);
    checketa->Draw();
    counter++;
    canvas->cd(counter);
    TH2D* projectioneta = (TH2D *) correlation->Projection(eta,eta+((Int_t)(ndimcor/2.)));
    projectioneta->Draw("colz");
  }


}
//____________________________________________________________
void AliHFEBeautySpectrumQA::DrawSubtractContamination() const
{
  //
  // get spectrum for beauty 2nd method
  //
  //
  TH1D *measuredTH1Daftersubstraction = (TH1D *) fListOfResult->UncheckedAt(kAfterSC);
  TH1D *measuredTH1Dbeforesubstraction = (TH1D *) fListOfResult->UncheckedAt(kBeforeSC);
  TH1D *measuredTH1background = (TH1D *) fListOfResult->UncheckedAt(kMeasBG);
  if(!measuredTH1Daftersubstraction || !measuredTH1Dbeforesubstraction) return;

  SetStyle();

  TCanvas * cbackgroundsubtraction = new TCanvas("backgroundsubtraction","backgroundsubtraction",1000,700);
  cbackgroundsubtraction->Divide(3,1);
  cbackgroundsubtraction->cd(1);
  gPad->SetLogy();
  gPad->SetTicks();
  measuredTH1Daftersubstraction->SetStats(0);
  measuredTH1Daftersubstraction->SetTitle("");
  measuredTH1Daftersubstraction->GetYaxis()->SetTitleOffset(1.5);
  measuredTH1Daftersubstraction->GetYaxis()->SetTitle("dN/dp_{T} [(GeV/c)^{-1}]");
  measuredTH1Daftersubstraction->GetXaxis()->SetTitle("p^{rec}_{T} [GeV/c]");
  //measuredTH1Daftersubstraction->GetXaxis()->SetRangeUser(0.0,fPtMax);
  //measuredTH1Daftersubstraction->SetMarkerStyle(25);
  //measuredTH1Daftersubstraction->SetMarkerColor(kBlack);
  //measuredTH1Daftersubstraction->SetLineColor(kBlack);
  measuredTH1Dbeforesubstraction->SetStats(0);
  measuredTH1Dbeforesubstraction->SetTitle("");
  measuredTH1Dbeforesubstraction->GetYaxis()->SetTitle("dN/dp_{T} [(GeV/c)^{-1}]");
  measuredTH1Dbeforesubstraction->GetXaxis()->SetTitle("p^{rec}_{T} [GeV/c]");
  //measuredTH1Dbeforesubstraction->GetXaxis()->SetRangeUser(0.0,fPtMax);
  //measuredTH1Dbeforesubstraction->SetMarkerStyle(24);
  //measuredTH1Dbeforesubstraction->SetMarkerColor(kBlue);
  //measuredTH1Dbeforesubstraction->SetLineColor(kBlue);
  measuredTH1Daftersubstraction->Draw();
  measuredTH1Dbeforesubstraction->Draw("same");
  TLegend *legsubstraction = new TLegend(0.4,0.6,0.89,0.89);
  legsubstraction->AddEntry(measuredTH1Dbeforesubstraction,"With hadron contamination","p");
  legsubstraction->AddEntry(measuredTH1Daftersubstraction,"Without hadron contamination ","p");
  legsubstraction->SetFillStyle(0);
  legsubstraction->SetLineStyle(0);
  legsubstraction->SetLineColor(0);
  legsubstraction->Draw("same");
  cbackgroundsubtraction->cd(2);
  gPad->SetLogy();
  gPad->SetTicks();
  TH1D* ratiomeasuredcontamination = (TH1D*)measuredTH1Dbeforesubstraction->Clone();
  ratiomeasuredcontamination->SetName("ratiomeasuredcontamination");
  ratiomeasuredcontamination->SetTitle("");
  ratiomeasuredcontamination->GetYaxis()->SetTitleOffset(1.5);
  ratiomeasuredcontamination->GetYaxis()->SetTitle("(with contamination - without contamination) / with contamination");
  ratiomeasuredcontamination->GetXaxis()->SetTitle("p^{rec}_{T} [GeV/c]");
  ratiomeasuredcontamination->GetYaxis()->SetRangeUser(0.8,1.2);
  ratiomeasuredcontamination->GetXaxis()->SetRangeUser(0.0,fPtMax);
  ratiomeasuredcontamination->Sumw2();
  ratiomeasuredcontamination->Add(measuredTH1Daftersubstraction,-1.0);
  ratiomeasuredcontamination->Divide(measuredTH1Dbeforesubstraction);
  ratiomeasuredcontamination->SetStats(0);
  ratiomeasuredcontamination->SetMarkerStyle(26);
  ratiomeasuredcontamination->SetMarkerColor(kBlack);
  ratiomeasuredcontamination->SetLineColor(kBlack);
  for(Int_t k=0; k < ratiomeasuredcontamination->GetNbinsX(); k++){
    ratiomeasuredcontamination->SetBinError(k+1,0.0);
  }
  ratiomeasuredcontamination->Draw("P");
  cbackgroundsubtraction->cd(3);
  measuredTH1background->SetStats(0);
  measuredTH1background->SetTitle("");
  measuredTH1background->GetYaxis()->SetTitle("dN/dp_{T} [(GeV/c)^{-1}]");
  measuredTH1background->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  measuredTH1background->SetMarkerStyle(26);
  measuredTH1background->SetMarkerColor(kBlack);
  measuredTH1background->SetLineColor(kBlack);
  measuredTH1background->Draw();
  if(fWriteToFile) cbackgroundsubtraction->SaveAs("BackgroundSubtracted.png");

}


//____________________________________________________________
void AliHFEBeautySpectrumQA::DrawCorrectWithEfficiency(Int_t typeeff) const
{
  //
  // Correct the spectrum for efficiency and unfolding
  // with both method and compare
  //
  
  TH1D *afterE = 0x0;
  TH1D *beforeE = 0x0;
  TH1D *efficiencyDproj = 0x0;
  TF1 *efficiencyparametrized = 0x0;
  
  if(typeeff== kMC) {
    afterE = (TH1D *) fListOfResult->UncheckedAt(kAfterMCE);
    beforeE = (TH1D *) fListOfResult->UncheckedAt(kBeforeMCE);
    efficiencyDproj = (TH1D *) fListOfResult->UncheckedAt(kMCEfficiency);
  }
  if(typeeff== kParametrized) {
    afterE = (TH1D *) fListOfResult->UncheckedAt(kAfterPE);
    beforeE = (TH1D *) fListOfResult->UncheckedAt(kBeforePE);
    efficiencyparametrized = (TF1 *) fListOfResult->UncheckedAt(kPEfficiency);
  }
  
  if(typeeff==kMC && (!afterE || !beforeE || !efficiencyDproj)) return;
  if(typeeff==kParametrized && (!afterE || !beforeE || !efficiencyparametrized)) return;
  
  SetStyle();
  
  TCanvas * cEfficiency = new TCanvas(AliHFEBeautySpectrumQA::fgkNameCanvas[typeeff],AliHFEBeautySpectrumQA::fgkNameCanvas[typeeff],1000,700);
  cEfficiency->Divide(2,1);
  cEfficiency->cd(1);
  gPad->SetLogy();
  gPad->SetTicks();
  afterE->SetStats(0);
  afterE->SetTitle("");
  afterE->GetYaxis()->SetTitleOffset(1.5);
  afterE->GetYaxis()->SetTitle("dN/dp_{T} [(GeV/c)^{-1}]");
  afterE->GetXaxis()->SetTitle("p^{rec}_{T} [GeV/c]");
  afterE->GetXaxis()->SetRangeUser(0.0,fPtMax);
  afterE->SetMarkerStyle(25);
  afterE->SetMarkerColor(kBlack);
  afterE->SetLineColor(kBlack);
  beforeE->SetStats(0);
  beforeE->SetTitle("");
  beforeE->GetYaxis()->SetTitle("dN/dp_{T} [(GeV/c)^{-1}]");
  beforeE->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  beforeE->GetXaxis()->SetRangeUser(0.0,fPtMax);
  beforeE->SetMarkerStyle(24);
  beforeE->SetMarkerColor(kBlue);
  beforeE->SetLineColor(kBlue);
  afterE->Draw();
  beforeE->Draw("same");
  TLegend *legefficiency = new TLegend(0.4,0.6,0.89,0.89);
  legefficiency->AddEntry(beforeE,"Before Efficiency correction","p");
  legefficiency->AddEntry(afterE,"After Efficiency correction","p");
  legefficiency->SetFillStyle(0);
  legefficiency->SetLineStyle(0);
  legefficiency->SetLineColor(0);
  legefficiency->Draw("same");
  cEfficiency->cd(2);
  gPad->SetTicks();
  if(typeeff==kMC) {
    if(efficiencyDproj) {
      efficiencyDproj->SetTitle("");
      efficiencyDproj->SetStats(0);
      efficiencyDproj->GetYaxis()->SetTitleOffset(1.5);
      efficiencyDproj->GetYaxis()->SetRangeUser(0.0,1.0);
      efficiencyDproj->GetYaxis()->SetTitle("Efficiency");
      efficiencyDproj->GetXaxis()->SetTitle("p^{rec}_{T} [GeV/c]");
      efficiencyDproj->GetXaxis()->SetRangeUser(0.0,fPtMax);
      efficiencyDproj->SetMarkerStyle(25);
      efficiencyDproj->Draw();
    }
  }
  if(typeeff==kParametrized) {
    if(efficiencyparametrized) efficiencyparametrized->Draw();
  }
  
  if(fWriteToFile) {
    if(typeeff==kMC) cEfficiency->SaveAs("EfficiencyMC.png");
    if(typeeff==kParametrized) cEfficiency->SaveAs("EfficiencyParametrized.png");
  }

}

//____________________________________________________________
void AliHFEBeautySpectrumQA::DrawUnfolding() const
{
  //
  // Draw unfolding
  //
  TH1D *measuredspectrumD = (TH1D *) fListOfResult->UncheckedAt(kBeforeU);
  TH1D *residualspectrumD = (TH1D *) fListOfResult->UncheckedAt(kResidualU);
  TH1D *efficiencyDproj = (TH1D *) fListOfResult->UncheckedAt(kUEfficiency);
  THnSparseF *correlation = (THnSparseF *) fListOfResult->UncheckedAt(kCMProjection);
  
  if(!measuredspectrumD || !residualspectrumD || !efficiencyDproj || !correlation) return;
  
  Int_t ndimcor = (Int_t) correlation->GetNdimensions()/2.;

  SetStyle();

  TCanvas * cunfolding = new TCanvas("unfolding","unfolding",1000,700);
  cunfolding->Divide(2,2);
  cunfolding->cd(1);
  gPad->SetLogy();
  gPad->SetTicks();
  residualspectrumD->GetYaxis()->SetTitle("dN/dp_{T} [(GeV/c)^{-1}]");
  residualspectrumD->GetXaxis()->SetTitle("p^{rec}_{T} [GeV/c]");
  residualspectrumD->GetXaxis()->SetRangeUser(0.0,fPtMax);
  residualspectrumD->SetStats(0);
  residualspectrumD->SetTitle("");
  residualspectrumD->GetYaxis()->SetTitleOffset(1.5);
  residualspectrumD->SetMarkerStyle(26);
  residualspectrumD->SetMarkerColor(kBlue);
  residualspectrumD->SetLineColor(kBlue);
  residualspectrumD->Sumw2();
  residualspectrumD->Draw("P");
  measuredspectrumD->SetStats(0);
  measuredspectrumD->SetTitle("");  
  measuredspectrumD->GetYaxis()->SetTitleOffset(1.5);
  measuredspectrumD->SetMarkerStyle(25);
  measuredspectrumD->SetMarkerColor(kBlack);
  measuredspectrumD->SetLineColor(kBlack);
  measuredspectrumD->Draw("same");
  TLegend *legres = new TLegend(0.4,0.6,0.89,0.89);
  legres->AddEntry(residualspectrumD,"Residual","p");
  legres->AddEntry(measuredspectrumD,"Measured","p");
  legres->SetFillStyle(0);
  legres->SetLineStyle(0);
  legres->SetLineColor(0);
  legres->Draw("same");
  cunfolding->cd(2);
  gPad->SetTicks();
  TH1D* ratioresidual = (TH1D*)residualspectrumD->Clone();
  ratioresidual->SetName("ratioresidual");
  ratioresidual->SetTitle("");
  ratioresidual->GetYaxis()->SetRangeUser(0.6,1.4);
  ratioresidual->GetYaxis()->SetTitle("Folded/Measured");
  ratioresidual->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  ratioresidual->Divide(measuredspectrumD);
  ratioresidual->SetStats(0);
  ratioresidual->Draw();
  cunfolding->cd(3);
  gPad->SetTicks();
  efficiencyDproj->GetYaxis()->SetTitle("Efficiency");
  efficiencyDproj->GetXaxis()->SetTitle("p^{MC}_{T} [GeV/c]");
  efficiencyDproj->GetXaxis()->SetRangeUser(0.0,fPtMax);
  efficiencyDproj->GetYaxis()->SetRangeUser(0.0,1.0);
  efficiencyDproj->SetStats(0);
  efficiencyDproj->SetTitle("");
  efficiencyDproj->GetYaxis()->SetTitleOffset(1.5);
  efficiencyDproj->SetMarkerStyle(26);
  efficiencyDproj->SetMarkerColor(kBlue);
  efficiencyDproj->SetLineColor(kBlue);
  efficiencyDproj->Sumw2();
  efficiencyDproj->Draw("P");
  cunfolding->cd(4);
  TH2F *projectioncorr = (TH2F *) correlation->Projection(0,ndimcor);
  projectioncorr->GetYaxis()->SetTitle("p^{ESD}_{T} [GeV/c]");
  projectioncorr->GetXaxis()->SetTitle("p^{MC}_{T} [GeV/c]");
  projectioncorr->GetXaxis()->SetRangeUser(0.0,fPtMax);
  projectioncorr->GetYaxis()->SetRangeUser(0.0,fPtMax);
  projectioncorr->SetStats(0);
  projectioncorr->SetTitle("");
  projectioncorr->Draw("colz");

  if(fWriteToFile){
    cunfolding->SaveAs("Unfolding.png");
  }
  
}
//____________________________________________________________
void AliHFEBeautySpectrumQA::DrawResult()
{
  //
  // Draw Results
  //
  TGraphErrors* correctedspectrumD = (TGraphErrors *) fListOfResult->UncheckedAt(kFinalResultUnfolded);
  TGraphErrors* alltogetherspectrumD = (TGraphErrors *) fListOfResult->UncheckedAt(kFinalResultDirectEfficiency);
  THnSparse* correctedspectrum = (THnSparse *) fListOfResult->UncheckedAt(kFinalResultUnfSparse);
  AliCFDataGrid* alltogetherCorrection = (AliCFDataGrid *) fListOfResult->UncheckedAt( kFinalResultDirectEffSparse);
  if(!correctedspectrumD || !alltogetherspectrumD) return;
  
  SetStyle();

  TCanvas * ccorrected = new TCanvas("corrected","corrected",1000,700);
  ccorrected->Divide(2,1);
  ccorrected->cd(1);
  gPad->SetLogy();
  correctedspectrumD->SetTitle("");
  correctedspectrumD->GetYaxis()->SetTitleOffset(1.5);
  correctedspectrumD->GetYaxis()->SetRangeUser(0.000000001,1.0);
  correctedspectrumD->SetMarkerStyle(26);
  correctedspectrumD->SetMarkerColor(kBlue);
  correctedspectrumD->SetLineColor(kBlue);
  correctedspectrumD->Draw("AP");
  alltogetherspectrumD->SetTitle("");
  alltogetherspectrumD->GetYaxis()->SetTitleOffset(1.5);
  alltogetherspectrumD->GetYaxis()->SetRangeUser(0.000000001,1.0);
  alltogetherspectrumD->SetMarkerStyle(25);
  alltogetherspectrumD->SetMarkerColor(kBlack);
  alltogetherspectrumD->SetLineColor(kBlack);
  alltogetherspectrumD->Draw("P");
  TLegend *legcorrected = new TLegend(0.4,0.6,0.89,0.89);
  legcorrected->AddEntry(correctedspectrumD,"Corrected","p");
  legcorrected->AddEntry(alltogetherspectrumD,"Alltogether","p");
  legcorrected->Draw("same");
  ccorrected->cd(2);
  TH1D *correctedTH1D = correctedspectrum->Projection(0);
  TH1D *alltogetherTH1D = (TH1D *) alltogetherCorrection->Project(0);
  TH1D* ratiocorrected = (TH1D*)correctedTH1D->Clone();
  ratiocorrected->SetName("ratiocorrected");
  ratiocorrected->SetTitle("");
  ratiocorrected->GetYaxis()->SetTitle("Unfolded/DirectCorrected");
  ratiocorrected->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  ratiocorrected->Divide(correctedTH1D,alltogetherTH1D,1,1);
  ratiocorrected->SetStats(0);
  ratiocorrected->Draw();
  if(fWriteToFile){ 
    ccorrected->SaveAs("CorrectedBeauty.eps");    
    TFile *out;
    out = new TFile("finalSpectrum.root","recreate");
    out->cd();
    //
    correctedspectrumD->SetName("UnfoldingCorrectedSpectrum");
    correctedspectrumD->Write();
    alltogetherspectrumD->SetName("AlltogetherSpectrum");
    alltogetherspectrumD->Write();
    ratiocorrected->SetName("RatioUnfoldingAlltogetherSpectrum");
    ratiocorrected->Write();
    //
 
    out->Close(); 
    delete out;
  }
}
//____________________________________________________________
TH1D *AliHFEBeautySpectrumQA::DivideSpectra(TGraphErrors *ga, TGraphErrors *gb) 
{
  //
  // Divide Spectra
  //

  TH1D *afterE = (TH1D *) fListOfResult->UncheckedAt(kAfterMCE);
  if(!afterE) return 0x0;

  TH1D *histoB = (TH1D*) afterE->Clone();
  histoB->Sumw2();
  histoB->SetName("ratio");
  TH1D *histoa = (TH1D*) afterE->Clone();
  histoa->Sumw2();
  histoa->SetName("a");
  TH1D *histob = (TH1D*) afterE->Clone();
  histob->Sumw2();
  histob->SetName("b");
  
  double xa,ya,xb,yb,eya,eyb;
  Int_t npointsa = ga->GetN();
  Int_t npointsb = gb->GetN();
  if(npointsa != npointsb) {
    printf("Problem the two spectra have not the same number of points\n");
    return 0x0;
  }
  for(Int_t k = 0; k < npointsa; k++){
    ga->GetPoint(k,xa,ya);
    gb->GetPoint(k,xb,yb);
    //
    Double_t centerhisto = histoa->GetBinCenter(k+1);
    //
    if((TMath::Abs(xa-xb) > 0.0001) || (TMath::Abs(xa-centerhisto) > 0.0001)) {
      printf("Problem not the same x axis\n");
      return 0x0;
    }
    histoa->SetBinContent(k+1,ya);
    histob->SetBinContent(k+1,yb);
    //
    eya = ga->GetErrorY(k);
    eyb = gb->GetErrorY(k);
    //
    histoa->SetBinError(k+1,eya);
    histob->SetBinError(k+1,eyb);
   
  }
  
  histoB->Sumw2();
  histoB->Divide(histoa,histob,1.0,1.0,"B");

  return histoB;  

}
//__________________________________________
void AliHFEBeautySpectrumQA::SetStyle() const
{
  //
  // Set style
  //

  gStyle->SetPalette(1);
  gStyle->SetOptStat(1111);
  gStyle->SetPadBorderMode(0);
  gStyle->SetCanvasColor(10);
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPadRightMargin(0.13);

}
