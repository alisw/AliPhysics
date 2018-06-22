//
// Macro to extract 1D or 2D D-meson efficiency maps from D2H CF output file
// It is possible to set the variables to be used on the axes and a cut-off value for the y-axis
// 
// Usage:
// gSystem->SetIncludePath("-I. -I$ALICE_ROOT/include -I$ALICE_PHYSICS/include -I$ROOTSYS/include");
// .L DrawEfficiency_HFCJ.C++
// and then call the needed methods
// Call steps 9,0 for RecoPID/LimAcc; otherwise 9,2 for RecoPID/GenAcc + Apply_GenAccLimAcc_Factor method to go to LimAcc
//
// Author: F.Colamaria (fabiocolamaria@cern.ch)
//

#include <iostream>
#include "TSystem.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TDirectoryFile.h"
#include "TAxis.h"
#include "TH1D.h"
#include "TH2D.h"

#include "AliCFContainer.h"

//#include "ComputeAcceptance.C"

using namespace std;

/* 
LIST OF VARIABLES FROM CF (in Cheetah mode) - Always cross-check with master!
 ipTFast = 0; 
 iyFast = 1;
 icTFast = 2;
 iphiFast = 3;
 izvtxFast = 4;
 icentFast = 5;
 ifakeFast = 6;
 imultFast = 7;

LIST OF STEPS FROM CF - Always cross-check with master!
 kStepGeneratedLimAcc = 0;
 kStepGenerated = 1;
 kStepAcceptance = 2;
 kStepVertex = 3;
 kStepRefit = 4;
 kStepReconstructed = 5;
 kStepRecoAcceptance = 6;
 kStepRecoITSClusters = 7;
 kStepRecoPPR = 8;
 kStepRecoPID = 9;
*/

/*
INFOS: 
numStep/denStep = CF steps for which the numerator/denominator shall be extracted;
xVar/yVar = Cf variables to plot on the x and y axis (x axis only for 1D maps);
lastVal_y = cutoff value for the y-axis (to remove low-stat fluctuating bins). After that value the maps are flattened
*/

void DrawEfficiency_2D(TString fileName="AnalysisResults.root", TString dirFileName="PWG3_D2H_CFtaskD0toKpi_c", TString contName="CFHFccontainer0_c", TString outFileName="EfficiencyMap_2D_Dzero_c.root", Int_t numStep=9, Int_t denStep=0, Int_t xVar=0, Int_t yVar=7, Double_t lastVal_y=-1);
void DrawEfficiency_1D(TString fileName="AnalysisResults.root", TString dirFileName="PWG3_D2H_CFtaskD0toKpi_c", TString contName="CFHFccontainer0_c", TString outFileName="EfficiencyMap_1D_Dzero_c.root", Int_t numStep=9, Int_t denStep=0, Int_t xVar=0);
void PlotEfficiency_2D(TString fileName="EfficiencyMap_2D_Dzero_c.root", TString fileNameOut="EfficiencyMap_2D_Dzero_c_Plot", Double_t xMin=0, Double_t xMax=24, Double_t yMin=0, Double_t yMax=90);
void PlotEfficiency_1D(TString fileName="EfficiencyMap_1D_Dzero_c.root", TString fileNameOut="EfficiencyMap_1D_Dzero_c_Plot", Double_t xMin=0, Double_t xMax=24);
void Apply_GenAccLimAcc_Factor_2D(TString inputEffFile="EfficiencyMap_2D_Dzero_c.root", TString MCtoyfile="Acceptance_Toy_D0Kpi_yfidPtDep_etaDau08_ptDau300_FONLL7ptshape.root", TString fileNameOut="EfficiencyMap_2D_Dzero_c_wLimAcc.root");
void Apply_GenAccLimAcc_Factor_1D(TString inputEffFile="EfficiencyMap_1D_Dzero_c.root", TString MCtoyfile="Acceptance_Toy_D0Kpi_yfidPtDep_etaDau08_ptDau300_FONLL7ptshape.root", TString fileNameOut="EfficiencyMap_1D_Dzero_c_wLimAcc.root");

void DrawEfficiency_2D(TString fileName, TString dirFileName, TString contName, TString outFileName, Int_t numStep, Int_t denStep, Int_t xVar, Int_t yVar, Double_t lastVal_y) {

  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);
  gStyle->SetCanvasColor(0);
  gStyle->SetFrameFillColor(0);
  gStyle->SetOptTitle(0);

  // Get the container
  TFile* f = new TFile(fileName.Data());
  TDirectoryFile* d = (TDirectoryFile*)f->Get(dirFileName.Data());
  AliCFContainer *data = (AliCFContainer*) (d->Get(contName.Data()));

  // Extract numerator and denominator
  TH2D *hNum = data->ShowProjection(xVar, yVar, numStep);
  hNum->Sumw2();
  TH2D *hDen = data->ShowProjection(xVar, yVar, denStep);
  hDen->Sumw2();

  // Rebin of multiplicity axis
  if(xVar==7) {
    hNum->RebinX(4);
    hDen->RebinX(4);
  } else if (yVar==7) {
    hNum->RebinY(4);
    hDen->RebinY(4);
  }

  if(lastVal_y>0) {
    Int_t binLim = hNum->GetYaxis()->FindBin(lastVal_y);
    printf("Flattening the efficiency on y axis from %f (bin %d) onwards\n",lastVal_y,binLim);
    for(Int_t i = 1; i <= hNum->GetNbinsX(); i++) {
      Double_t numSum = 0;
      Double_t denSum = 0;
      for(Int_t j = binLim+1; j <= hNum->GetNbinsY(); j++) { //from the bin following the one containing 'lastVal_y', onwards: sum entries in num and den
        numSum += hNum->GetBinContent(i,j,hNum->GetBinContent(i,j));
        denSum += hDen->GetBinContent(i,j,hDen->GetBinContent(i,j));
      }
      for(Int_t j = binLim+1; j <= hNum->GetNbinsY(); j++) { //from the bin following the one containing 'lastVal_y', onwards: set to each bin num_sum and den_sum (so in the division in each bin you get the average eff of those bins)
        hNum->SetBinContent(i,j,numSum);
        hNum->SetBinError(i,j,TMath::Sqrt(numSum));
        hDen->SetBinContent(i,j,denSum);
        hDen->SetBinError(i,j,TMath::Sqrt(denSum));
      }
    }
  } //end if
  
  //Evaluate the efficiency
  TH2D* hEff = (TH2D*)hNum->Clone("h_Eff");  
  hEff->Divide(hNum,hDen,1,1,"b"); //"b2 = apply binomial error

  for(Int_t i = 1; i <= hEff->GetNbinsX(); i++) {
    for(Int_t j = 1; j <= hEff->GetNbinsY(); j++) {
      if(hDen->GetBinContent(i,j) == 0) {
        if (hNum->GetBinContent(i,j) != 0) printf("Warning! Den. set at 1 with non-0 num, in (%d,%d) bin! Values: num = %1.3f, den = %1.3f\n",i,j,hNum->GetBinContent(i,j),hDen->GetBinContent(i,j));
        hDen->SetBinContent(i,j,1);  //if you have some Inf/NaN propagated...
      }
      printf("Num, den bin %d,%d (vals x: %.1f-%.1f, y: %.1f-%.1f) = %1.3f, %1.3f\n",i,j,hNum->GetXaxis()->GetBinLowEdge(i),hNum->GetXaxis()->GetBinUpEdge(i),hNum->GetYaxis()->GetBinLowEdge(j),hNum->GetYaxis()->GetBinUpEdge(j),hNum->GetBinContent(i,j),hDen->GetBinContent(i,j));
      printf("hEff bin %d,%d  (vals x: %.1f-%.1f, y: %.1f-%.1f) = %1.3f +- %1.3f\n",i,j,hNum->GetXaxis()->GetBinLowEdge(i),hNum->GetXaxis()->GetBinUpEdge(i),hNum->GetYaxis()->GetBinLowEdge(j),hNum->GetYaxis()->GetBinUpEdge(j),hEff->GetBinContent(i,j),hEff->GetBinError(i,j));
    }
  }

  // Save the map
  TFile* effMapFile = new TFile(Form("%s",outFileName.Data()),"RECREATE");
  hEff->SetDrawOption("lego2");
  hEff->Write("h_Eff");
  effMapFile->Close();

  return;
}


void DrawEfficiency_1D(TString fileName, TString dirFileName, TString contName, TString outFileName, Int_t numStep, Int_t denStep, Int_t xVar) {

  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);
  gStyle->SetCanvasColor(0);
  gStyle->SetFrameFillColor(0);
  gStyle->SetOptTitle(0);

  // Get the container
  TFile* f = new TFile(fileName.Data());
  TDirectoryFile* d = (TDirectoryFile*)f->Get(dirFileName.Data());
  AliCFContainer *data = (AliCFContainer*) (d->Get(contName.Data()));

  // Extract numerator and denominator
  TH1D *hNum = data->ShowProjection(xVar, numStep);
  hNum->Sumw2();
  TH1D *hDen = data->ShowProjection(xVar, denStep);
  hDen->Sumw2();

  // Rebin of multiplicity axis
  if(xVar==7) {
    hNum->RebinX(4);
    hDen->RebinX(4);
  } 
  
  //Evaluate the efficiency
  TH1D* hEff = (TH1D*)hNum->Clone("h_Eff");  
  hEff->Divide(hNum,hDen,1,1,"b"); //"b2 = apply binomial error

  for(Int_t i = 1; i <= hEff->GetNbinsX(); i++) {
    if(hDen->GetBinContent(i) == 0) {
      if (hDen->GetBinContent(i) != 0) printf("Warning! Den. set at 1 with non-0 num, in %d bin! Values: num = %1.3f, den = %1.3f\n",i,hNum->GetBinContent(i),hDen->GetBinContent(i));
      hDen->SetBinContent(i,1);  //if you have some Inf/NaN propagated...
    }
    printf("Num, den %d = %1.2f, %1.2f\n",i,hNum->GetBinContent(i),hDen->GetBinContent(i));
    printf("hEff %d = %1.2f +- %1.2f\n",i,hEff->GetBinContent(i),hEff->GetBinError(i));
  }

  // Save the map
  TFile* effMapFile = new TFile(Form("%s",outFileName.Data()),"RECREATE");
  hEff->SetDrawOption("lego2");
  hEff->Write("h_Eff");
  effMapFile->Close();

  return;
}


void PlotEfficiency_2D(TString fileName, TString fileNameOut, Double_t xMin, Double_t xMax, Double_t yMin, Double_t yMax) {

  TFile* effMapFile = new TFile(fileName.Data(),"READ");  
  TH2D *effMap = (TH2D*)effMapFile->Get("h_Eff");

  TCanvas *c = new TCanvas("c","2D Efficiency map",150,150,900,600);
  c->cd();
  effMap->GetXaxis()->SetRangeUser(xMin,xMax);
  effMap->GetYaxis()->SetRangeUser(yMin,yMax);
  effMap->GetXaxis()->SetTitleOffset(1.5);
  effMap->GetYaxis()->SetTitleOffset(1.5);
  effMap->DrawCopy("lego2");

  c->SaveAs(Form("%s.png",fileNameOut.Data()));
  c->SaveAs(Form("%s.root",fileNameOut.Data()));
  effMapFile->Close();

  return;
}


void PlotEfficiency_1D(TString fileName, TString fileNameOut, Double_t xMin, Double_t xMax) {

  TFile *effMapFile = new TFile(fileName.Data(),"READ");  
  TH1D *effMap = (TH1D*)effMapFile->Get("h_Eff");

  TCanvas *c = new TCanvas("c","1D Efficiency map",150,150,900,600);
  c->cd();
  effMap->GetXaxis()->SetRangeUser(xMin,xMax);
  effMap->GetXaxis()->SetTitleOffset(1.5);
  effMap->SetLineColor(kBlack);
  effMap->SetLineWidth(2);
  effMap->SetMarkerColor(kBlue);
  effMap->SetMarkerStyle(21);
  effMap->SetMarkerSize(0.9);
  effMap->DrawCopy();

  c->SaveAs(Form("%s.png",fileNameOut.Data()));
  c->SaveAs(Form("%s.root",fileNameOut.Data()));
  effMapFile->Close();

  return;
}

void Apply_GenAccLimAcc_Factor_2D(TString inputEffFile, TString MCtoyfile, TString fileNameOut) {

/* //NOT WORKING, COMMENTED OUT	
  // Run MC toy to generate GenAcc/LimAcc factor
  if(runMCtoy) {
  	printf("Running MC toy to generate GenAcc/LimAcc factor...\n");
  	ComputeAcceptance();
  } */
  
  // Load GenAcc and LimAcc plots  
  printf("Apply GenAcc/LimAcc factor to the efficiency...\n");
  TFile *fMap = new TFile(inputEffFile.Data(),"read");
  TH2D *hMap2D = (TH2D*)fMap->Get("h_Eff");
  
  TFile *fFact = new TFile(MCtoyfile.Data(),"read");
  TH1D *hGenAcc = (TH1D*)fFact->Get("hPtGenAcc");
  TH1D *hLimAcc = (TH1D*)fFact->Get("hPtGenLimAcc");  
  
  // Rebin to map x-axis binning and produce ratio GenAcc/LimAcc
  Int_t xBins = hMap2D->GetNbinsX();
  Double_t xVals[xBins+1];
  for(Int_t i=0; i<xBins; i++) xVals[i] = hMap2D->GetXaxis()->GetBinLowEdge(i+1);
  xVals[xBins] = hMap2D->GetXaxis()->GetBinUpEdge(xBins);

  TH1D *hGenAcc_Reb = (TH1D*)hGenAcc->Rebin(xBins,"hGenAcc_Reb",xVals);
  TH1D *hLimAcc_Reb = (TH1D*)hLimAcc->Rebin(xBins,"hLimAcc_Reb",xVals);  
  
  TH1D* hRatio=(TH1D*)hGenAcc_Reb->Clone("hRatio");  
  hRatio->Divide(hGenAcc_Reb,hLimAcc_Reb,1,1,"B");
  
  // Go 2-dimensional
  Int_t yBins = hMap2D->GetNbinsY();
  Double_t yVals[yBins+1];
  for(Int_t j=0; j<yBins; j++) yVals[j] = hMap2D->GetYaxis()->GetBinLowEdge(j+1);
  yVals[yBins] = hMap2D->GetYaxis()->GetBinUpEdge(yBins);  

  TH2D *hRatio2D = new TH2D("hRatio2D","hRatio2D",xBins,xVals,yBins,yVals);
  for(Int_t i=1; i<=xBins; i++) {
  	for(Int_t j=1; j<=yBins; j++) {
  	  hRatio2D->SetBinContent(i,j,hRatio->GetBinContent(i));	
  	  hRatio2D->SetBinError(i,j,0); //treat this correction as a pure number, with no uncertainty
	}
  }
  
//for(Int_t i=0; i<=xBins; i++) printf("xaxis %i = %f\n",i,xVals[i]); //DEBUG
//for(Int_t i=0; i<=yBins; i++) printf("yaxis %i = %f\n",i,yVals[i]);
    
  //Correct the efficiency maps and save it on new file
  hMap2D->Multiply(hRatio2D);
  
  TFile *fOut = new TFile(fileNameOut.Data(),"recreate");
  hMap2D->Write();
  // hRatio->Write(); //DEBUG
  // hRatio2D->Write(); //DEBUG 
  
  fMap->Close();
  fFact->Close();
  fOut->Close();

  return;
}

void Apply_GenAccLimAcc_Factor_1D(TString inputEffFile, TString MCtoyfile, TString fileNameOut) {

/* //NOT WORKING, COMMENTED OUT	
  // Run MC toy to generate GenAcc/LimAcc factor
  if(runMCtoy) {
  	printf("Running MC toy to generate GenAcc/LimAcc factor...\n");
  	ComputeAcceptance();
  } */
  
  // Load GenAcc and LimAcc plots  
  printf("Apply GenAcc/LimAcc factor to the efficiency...\n");
  TFile *fMap = new TFile(inputEffFile.Data(),"read");
  TH1D *hMap1D = (TH1D*)fMap->Get("h_Eff");
  
  TFile *fFact = new TFile(MCtoyfile.Data(),"read");
  TH1D *hGenAcc = (TH1D*)fFact->Get("hPtGenAcc");
  TH1D *hLimAcc = (TH1D*)fFact->Get("hPtGenLimAcc");  
  
  // Rebin to map x-axis binning and produce ratio GenAcc/LimAcc
  Int_t xBins = hMap1D->GetNbinsX();
  Double_t xVals[xBins+1];
  for(Int_t i=0; i<xBins; i++) xVals[i] = hMap1D->GetXaxis()->GetBinLowEdge(i+1);
  xVals[xBins] = hMap1D->GetXaxis()->GetBinUpEdge(xBins);

  TH1D *hGenAcc_Reb = (TH1D*)hGenAcc->Rebin(xBins,"hGenAcc_Reb",xVals);
  TH1D *hLimAcc_Reb = (TH1D*)hLimAcc->Rebin(xBins,"hLimAcc_Reb",xVals);  
  
  TH1D* hRatio=(TH1D*)hGenAcc_Reb->Clone("hRatio");  
  hRatio->Divide(hGenAcc_Reb,hLimAcc_Reb,1,1,"B");
  
//for(Int_t i=0; i<=xBins; i++) printf("xaxis %i = %f\n",i,xVals[i]); //DEBUG
    
  //Correct the efficiency maps and save it on new file
  for(Int_t i=1; i<=xBins; i++) {
  	hRatio->SetBinError(i,0); //treat this correction as a pure number, with no uncertainty
  }  
  hMap1D->Multiply(hRatio);
  
  TFile *fOut = new TFile(fileNameOut.Data(),"recreate");
  hMap1D->Write();
  // hRatio->Write(); //DEBUG
  
  fMap->Close();
  fFact->Close();
  fOut->Close();

  return;
}
