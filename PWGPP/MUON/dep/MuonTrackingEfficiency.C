/**
 *
 * \file MuonTrackingEfficiency.C
 * \author Philippe Pillot, Antoine Lardeux, Lizardo Valencia Palomo, Javier Martin Blanco
 * \brief Compute trk efficiency at DE, chamber, station and spectro levels vs various variables from the output of the efficiency task
 *
 * 1) efficiency estimator and error calculation at chamber and DE level:
 *  2 options:
 *    - using bayesian method with uniform prior
 *    - using Clopper-Pearson or other frequentist method
 *  if n ≠ 0: use above methods
 *  if n = 0: eff = -1 ± 0
 *
 * 2) efficiency and error propagation at station and spectrometer level:
 *  if eff = -1 for one or several ch/DE:
 *    - assume eff_ch = 1 ± 0 to compute eff_up and err_up with std error propagation at nth order
 *    - assume eff_ch = 0 ± 0 to compute eff_low and err_low with std error propagation at nth order
 *    - eff_spectro = eff_up + err_up - (eff_up-eff_low + err_low)
 *  otherwise: std efficiency and error propagation at nth order
 *
 * 3) efficiency and error integration over runs
 *  - performed from efficiency plots versus run
 *  - compute weighted average and apply standard error propagation to err_up and err_low
 *  - do it for both extreme cases above (with eff_up ± err–up and eff_low ± err_low)
 *  - int_eff = int_eff_up + int_err_up - (int_eff_up-int_eff_low + int_err_low)
 *
 */

#include <TROOT.h>
#include <TH1.h>
#include <TH2.h>
#include <THnSparse.h>
#include <TAxis.h>
#include <TString.h>
#include <TObjString.h>
#include <Riostream.h>
#include <TFile.h>
#include <TList.h>
#include <TCanvas.h>
#include <TGraphAsymmErrors.h>
#include <TMath.h>
#include <TArrayD.h>
#include <TStyle.h>
#include <THashList.h>
#include <TParameter.h>

//const Char_t *effErrMode = "cp"; // Clopper-Pearson
const Char_t *effErrMode = "b(1,1)mode"; // Bayesian with uniform prior

Double_t centMin = -999.;
Double_t centMax = 999.;
Double_t ptMin = 0.;
Double_t ptMax = -1.;
Double_t yMin = -99.;
Double_t yMax = 99.;
Double_t phiMin = -99.;
Double_t phiMax = 99.;
Int_t charge = 0; // 0 selects + and -, -1 and +1 selects - or + muons respectively

Bool_t moreTrackCandidates = kFALSE;

THashList *runWeights = 0x0;
TString gObjNameExtension = "";


void PlotMuonEfficiencyVsX(TString var, TString fileNameData, TString fileNameSave, Bool_t saveEdges, Bool_t print, Bool_t draw);
void PlotIntegratedMuonEfficiencyVsX(TString var, TString runList, TString fileNameWeights,
                                     TString fileNameData, TString fileNameSave, Bool_t print, Bool_t draw);

void PlotMuonEfficiencyVsXY(TString xVar, TString yVar, TString fileNameData, TString fileNameSave, Bool_t draw, Bool_t rap = kFALSE);

void PlotMuonEfficiency(TString fileNameData, TString fileNameSave, Bool_t saveEdges, Bool_t print, Bool_t draw);
void PlotMuonEfficiencyVsRun(TString runList, TString fileNameData, TString fileNameSave, Bool_t print, Bool_t draw);
void PlotIntegratedMuonEfficiency(TString fileNameWeights, TString fileNameSave, Bool_t print, Bool_t draw);

void PlotMuonEfficiencyPerDE(TString fileNameData, TString fileNameSave, Bool_t saveEdges);
void PlotMuonEfficiencyPerDEVsRun(TString runList, TString fileNameData, TString fileNameSave);
void PlotIntegratedMuonEfficiencyPerDE(TString fileNameWeights, TString fileNameSave);

Bool_t GetChamberEfficiency(THnSparse &TT, THnSparse &TD, TArrayD &chEff, TArrayD chEffErr[2], Bool_t printError = kFALSE);
void GetDEEfficiency(THnSparse &TT, THnSparse &TD, TGraphAsymmErrors &effVsDE);
void ComputeStationEfficiency(TArrayD &chEff, TArrayD chEffErr[2], Int_t iSt, Double_t &stEff, Double_t stEffErr[2]);
void GetStationEfficiency(TArrayD &chEff, TArrayD chEffErr[2], Int_t iSt, TGraphAsymmErrors *effVsX[3], Int_t ip, Double_t xp);
void ComputeStation45Efficiency(TArrayD &chEff, TArrayD chEffErr[2], Double_t &st45Eff, Double_t st45EffErr[2]);
void GetStation45Efficiency(TArrayD &chEff, TArrayD chEffErr[2], TGraphAsymmErrors *effVsX[3], Int_t ip, Double_t xp);
void ComputeTrackingEfficiency(Double_t stEff[6], Double_t stEffErr[6][2], Double_t &spectroEff, Double_t spectroEffErr[2]);
void GetTrackingEfficiency(TArrayD &chEff, TArrayD chEffErr[2], TGraphAsymmErrors *effVsSt[3],
                           TGraphAsymmErrors *effVsX[3], Int_t ip, Double_t xp, Bool_t print = kFALSE);
void IntegrateMuonEfficiency(TGraphAsymmErrors &effVsRunLow, TGraphAsymmErrors &effVsRunUp,
                             TGraphAsymmErrors &effVsX, Int_t ip, Double_t xp);

void LoadRunWeights(TString fileName);
void SetCentPtCh(THnSparse& SparseData);
TGraphAsymmErrors* CreateGraph(const char* name, const char* title, int value=-1);
void BeautifyGraph(TGraphAsymmErrors &g, const char* xAxisName, const char* yAxisName);
void BeautifyGraphs(TObjArray& array, const char* xAxisName, const char* yAxisName);
void SetRunLabel(TGraphAsymmErrors &g, Int_t irun, const TList& runs);
void SetRunLabel(TObjArray& array, Int_t irun, const TList& runs);


//---------------------------------------------------------------------------
void MuonTrackingEfficiency(TString runList = "runList.txt",
                            TString fileNameWeights = "",
                            TString objNameExtension = "",
                            TString fileNameData ="AnalysisResults.root",
                            TString fileNameSave = "efficiency_new.root")
{
  /// main function to compute, print and plot efficiencies
  
  gObjNameExtension = objNameExtension;
  
  PlotMuonEfficiencyVsX("centrality", fileNameData, fileNameSave, kFALSE, kFALSE, kTRUE);
  PlotMuonEfficiencyVsX("pt", fileNameData, fileNameSave, kFALSE, kFALSE, kTRUE);
  PlotMuonEfficiencyVsX("y", fileNameData, fileNameSave, kFALSE, kFALSE, kTRUE);
  PlotMuonEfficiencyVsX("phi", fileNameData, fileNameSave, kFALSE, kFALSE, kTRUE);
  PlotMuonEfficiencyVsX("charge", fileNameData, fileNameSave, kFALSE, kFALSE, kTRUE);
  
  PlotIntegratedMuonEfficiencyVsX("centrality", runList, fileNameWeights, fileNameData, fileNameSave, kFALSE, kTRUE);
  PlotIntegratedMuonEfficiencyVsX("pt", runList, fileNameWeights, fileNameData, fileNameSave, kFALSE, kTRUE);
  PlotIntegratedMuonEfficiencyVsX("y", runList, fileNameWeights, fileNameData, fileNameSave, kFALSE, kTRUE);
  PlotIntegratedMuonEfficiencyVsX("phi", runList, fileNameWeights, fileNameData, fileNameSave, kFALSE, kTRUE);
  PlotIntegratedMuonEfficiencyVsX("charge", runList, fileNameWeights, fileNameData, fileNameSave, kFALSE, kTRUE);
  
  PlotMuonEfficiencyVsXY("pt", "centrality", fileNameData, fileNameSave, kTRUE);
  PlotMuonEfficiencyVsXY("y", "centrality", fileNameData, fileNameSave, kTRUE);
  PlotMuonEfficiencyVsXY("pt", "y", fileNameData, fileNameSave, kTRUE);
  PlotMuonEfficiencyVsXY("phi", "y", fileNameData, fileNameSave, kTRUE, kTRUE);
  
  PlotMuonEfficiency(fileNameData, fileNameSave, kFALSE, kTRUE, kTRUE);
  PlotMuonEfficiencyVsRun(runList, fileNameData, fileNameSave, kFALSE, kTRUE);
  PlotIntegratedMuonEfficiency(fileNameWeights, fileNameSave, kTRUE, kTRUE);
  
  PlotMuonEfficiencyPerDE(fileNameData, fileNameSave, kFALSE);
  PlotMuonEfficiencyPerDEVsRun(runList, fileNameData, fileNameSave);
  PlotIntegratedMuonEfficiencyPerDE(fileNameWeights, fileNameSave);
  
}


//---------------------------------------------------------------------------
void PlotMuonEfficiencyVsX(TString var, TString fileNameData, TString fileNameSave, Bool_t saveEdges, Bool_t print, Bool_t draw)
{
  /// plot the tracking efficiency versus X
  
  printf("plotting efficiency versus %s...\n", var.Data());
  
  Int_t xDim = -1;
  if (var == "centrality") xDim = 1;
  else if (var == "pt") xDim = 2;
  else if (var == "y") xDim = 3;
  else if (var == "phi") xDim = 4;
  else if (var == "charge") xDim = 5;
  else {
    printf("incorrect variable. Choices are centrality, pt, y, phi and charge.\n");
    return;
  }
  
  // get input hists
  TFile *file = new TFile(fileNameData.Data(), "read");
  if (!file || !file->IsOpen()) {
    printf("cannot open file %s \n",fileNameData.Data());
    return;
  }
  TList *listTT = static_cast<TList*>(file->FindObjectAny(Form("TotalTracksPerChamber%s", gObjNameExtension.Data())));
  TList *listTD = static_cast<TList*>(file->FindObjectAny(Form("TracksDetectedPerChamber%s", gObjNameExtension.Data())));
  THnSparse *TT = static_cast<THnSparse*>(listTT->At(10));
  THnSparse *TD = static_cast<THnSparse*>(listTD->At(10));
  
  // output graph
  TGraphAsymmErrors *effVsX[3] = {0x0, 0x0, 0x0};
  TString nameAdd[3] = {"", "Low", "Up"};
  TString titleAdd[3] = {"", " - lower limit", " - upper limit"};
  for (Int_t i = 0; i < 1 || (saveEdges && i < 3); ++i) {
    effVsX[i] = CreateGraph(Form("trackingEffVs%s%s",var.Data(),nameAdd[i].Data()),
                            Form("Measured tracking efficiency versus %s%s",var.Data(),titleAdd[i].Data()));
  }
  
  // set the centrality and pT range for integration
  SetCentPtCh(*TT);
  SetCentPtCh(*TD);
  
  // loop over X bins
  TArrayD chEff(11);
  TArrayD chEffErr[2];
  chEffErr[0].Set(11);
  chEffErr[1].Set(11);
  for (Int_t ix = 1; ix <= TT->GetAxis(xDim)->GetNbins(); ix++) {
    
    if (print) cout << var.Data() << " " << TT->GetAxis(xDim)->GetBinLowEdge(ix) << "-" << TT->GetAxis(xDim)->GetBinUpEdge(ix) << ":" << endl;
    
    // set the var range to the current bin
    TT->GetAxis(xDim)->SetRange(ix, ix);
    TD->GetAxis(xDim)->SetRange(ix, ix);
    
    // compute chamber efficency and errors
    if (GetChamberEfficiency(*TT, *TD, chEff, chEffErr, print)) {
      
      // compute the overall tracking efficiency
      TGraphAsymmErrors *dummy[3] = {0x0, 0x0, 0x0};
      GetTrackingEfficiency(chEff, chEffErr, dummy, effVsX, ix-1, TT->GetAxis(xDim)->GetBinCenter(ix), print);
      
    } else {
      
      // fill graph with 1 +0 -1
      effVsX[0]->SetPoint(ix-1,TT->GetAxis(xDim)->GetBinCenter(ix),1.);
      effVsX[0]->SetPointError(ix-1,0.,0.,1.,0.);
      
      if (saveEdges) {
        
        // lower = 0 ± 0
        effVsX[1]->SetPoint(ix-1,TT->GetAxis(xDim)->GetBinCenter(ix),0.);
        effVsX[1]->SetPointError(ix-1,0.,0.,0.,0.);
        
        // upper = 1 ± 0
        effVsX[2]->SetPoint(ix-1,TT->GetAxis(xDim)->GetBinCenter(ix),1.);
        effVsX[2]->SetPointError(ix-1,0.,0.,0.,0.);
        
      }
      
    }
    
  }
  
  // close input file
  file->Close();
  
  // display
  for (Int_t i = 0; i < 1 || (saveEdges && i < 3); ++i) BeautifyGraph(*effVsX[i], var.Data(), "efficiency");
  
  if (draw) {
    new TCanvas(Form("cTrackingEffVs%s",var.Data()), Form("Measured tracking efficiency versus %s",var.Data()),1000,400);
    effVsX[0]->DrawClone("ap");
  }
  
  // save output
  file = new TFile(fileNameSave.Data(),"update");
  for (Int_t i = 0; i < 1 || (saveEdges && i < 3); ++i) effVsX[i]->Write(0x0, TObject::kOverwrite);
  file->Close();
  
  // clean memory
  for (Int_t i = 0; i < 1 || (saveEdges && i < 3); ++i) delete effVsX[i];
  
}


//---------------------------------------------------------------------------
void PlotIntegratedMuonEfficiencyVsX(TString var, TString runList, TString fileNameWeights,
                                     TString fileNameData, TString fileNameSave, Bool_t print, Bool_t draw)
{
  /// plot the tracking efficiency versus X for each run and integrated
  
  printf("plotting integrated efficiency versus %s...\n", var.Data());
  
  // load run weights
  if (!fileNameWeights.IsNull()) LoadRunWeights(fileNameWeights);
  if (!runWeights) {
    printf("Cannot compute integrated efficiency without run-by-run weights\n");
    return;
  }
  
  // open run list
  ifstream inFile(runList.Data());
  if (!inFile.is_open()) {
    printf("cannot open file %s\n",runList.Data());
    return;
  }
  
  // output graph
  TGraphAsymmErrors *intEffVsX = CreateGraph(Form("integratedTrackingEffVs%s",var.Data()),
                                             Form("Integrated tracking efficiency versus %s",var.Data()));
  
  Int_t n = -1;
  TArrayD x;
  TArrayD rec[2];
  Double_t gen = 0.;
  TArrayD effErr[2];
  
  // loop over runs
  while (!inFile.eof()) {
    
    // get current run number
    TString currRun;
    currRun.ReadLine(inFile,kTRUE);
    if(currRun.IsNull() || !currRun.IsDec()) continue;
    Int_t run = currRun.Atoi();
    
    printf("run %d: ", run);
    
    // compute efficiency vs var
    TString dataFile = Form("runs/%d/%s", run, fileNameData.Data());
    TString outFile = Form("runs/%d/%s", run, fileNameSave.Data());
    PlotMuonEfficiencyVsX(var, dataFile, outFile, kTRUE, print, kFALSE);
    
    // get run weight
    TParameter<Double_t> *weight = static_cast<TParameter<Double_t>*>(runWeights->FindObject(currRun.Data()));
    if (!weight) {
      printf("weight not found for run %s\n", currRun.Data());
      continue;
    }
    Double_t w = weight->GetVal();
    Double_t w2 = w*w;
    
    // get result
    TFile *file = new TFile(outFile.Data(), "read");
    if (!file || !file->IsOpen()) {
      printf("cannot open file\n");
      continue;
    }
    TGraphAsymmErrors *effVsX[2];
    effVsX[0] = static_cast<TGraphAsymmErrors*>(file->FindObjectAny(Form("trackingEffVs%sLow",var.Data())));
    effVsX[1] = static_cast<TGraphAsymmErrors*>(file->FindObjectAny(Form("trackingEffVs%sUp",var.Data())));
    if (!effVsX[0] || !effVsX[1]) {
      printf("trackingEffVs%sLow(Up) object not found\n", var.Data());
      continue;
    }
    
    // prepare the arrays if not already done
    if (n < 0) {
      n = effVsX[0]->GetN();
      x.Set(n, effVsX[0]->GetX());
      for (Int_t i = 0; i < 2; ++i) {
        rec[i].Set(n);
        effErr[i].Set(n);
      }
    } else if (n != effVsX[0]->GetN()) {
      printf("number of points in graph trackingEffVs%sLow(Up) for run %d is different than from previous runs\n", var.Data(), run);
      continue;
    }
    
    // integrate for all bins
    gen += w;
    for (Int_t ix = 0; ix < n; ++ix) {
      Double_t ieffErr[2] = {effVsX[0]->GetErrorYlow(ix), effVsX[1]->GetErrorYhigh(ix)};
      for (Int_t i = 0; i < 2; ++i) {
        rec[i][ix] += w*effVsX[i]->GetY()[ix];
        effErr[i][ix] += w2*ieffErr[i]*ieffErr[i];
      }
    }
    
    // close input file
    file->Close();
    
  }
  
  inFile.close();
  
  // fill output graph
  if (gen > 0.) {
    
    for (Int_t ix = 0; ix < n; ++ix) {
      
      intEffVsX->SetPoint(ix, x[ix], rec[1][ix]/gen);
      intEffVsX->SetPointError(ix, 0., 0., (rec[1][ix]-rec[0][ix]+TMath::Sqrt(effErr[0][ix]))/gen, TMath::Sqrt(effErr[1][ix])/gen);
      
    }
    
  } else {
    
    for (Int_t ix = 0; ix < n; ++ix) {
      
      printf("impossible to integrate, all weights = 0 or unknown ?!?\n");
      
      intEffVsX->SetPoint(ix, x[ix], -1.);
      intEffVsX->SetPointError(ix, 0., 0., 0., 0.);
      
    }
    
  }
  
  // display
  BeautifyGraph(*intEffVsX, var.Data(), "efficiency");
  
  if (draw) {
    new TCanvas(Form("cIntegratedTrackingEffVs%s",var.Data()), Form("Integrated tracking efficiency versus %s",var.Data()),1000,400);
    intEffVsX->DrawClone("ap");
  }
  
  // save output
  TFile *file = new TFile(fileNameSave.Data(),"update");
  intEffVsX->Write(0x0, TObject::kOverwrite);
  file->Close();
  
  // clean memory
  delete intEffVsX;
  
}


//---------------------------------------------------------------------------
void PlotMuonEfficiencyVsXY(TString xVar, TString yVar, TString fileNameData, TString fileNameSave, Bool_t draw, Bool_t rap)
{
  /// plot the tracking efficiency versus X,Y
  
  printf("plotting efficiency versus %s/%s...\n", xVar.Data(), yVar.Data());
  
  Int_t xDim = -1;
  if (xVar == "centrality") xDim = 1;
  else if (xVar == "pt") xDim = 2;
  else if (xVar == "y") xDim = 3;
  else if (xVar == "phi") xDim = 4;
  else if (xVar == "charge") xDim = 5;
  else {
    printf("incorrect variable. Choices are centrality, pt, y, phi and charge.\n");
    return;
  }
  Int_t yDim = -1;
  if (yVar == "centrality") yDim = 1;
  else if (yVar == "pt") yDim = 2;
  else if (yVar == "y") yDim = 3;
  else if (yVar == "phi") yDim = 4;
  else if (yVar == "charge") yDim = 5;
  else {
    printf("incorrect variable. Choices are centrality, pt, y, phi and charge.\n");
    return;
  }
  
  // get input hists
  TFile *file = new TFile(fileNameData.Data(), "read");
  if (!file || !file->IsOpen()) {
    printf("cannot open file %s\n",fileNameData.Data());
    return;
  }
  TList *listTT = static_cast<TList*>(file->FindObjectAny(Form("TotalTracksPerChamber%s", gObjNameExtension.Data())));
  TList *listTD = static_cast<TList*>(file->FindObjectAny(Form("TracksDetectedPerChamber%s", gObjNameExtension.Data())));
  THnSparse *TT = static_cast<THnSparse*>(listTT->At(10));
  THnSparse *TD = static_cast<THnSparse*>(listTD->At(10));
  
  // output map
  Int_t nxBins = TT->GetAxis(xDim)->GetNbins();
  Int_t nyBins = TT->GetAxis(yDim)->GetNbins();
  TH2F *effVsXY = new TH2F(Form("trackingEffVs%s-%s",xVar.Data(),yVar.Data()),
			   Form("Measured tracking efficiency versus %s and %s",xVar.Data(),yVar.Data()),
			   nxBins, TT->GetAxis(xDim)->GetBinLowEdge(1), TT->GetAxis(xDim)->GetBinUpEdge(nxBins),
			   nyBins, TT->GetAxis(yDim)->GetBinLowEdge(1), TT->GetAxis(yDim)->GetBinUpEdge(nyBins));
  effVsXY->SetDirectory(0);
  
  // set the centrality and pT range for integration
  SetCentPtCh(*TT);
  SetCentPtCh(*TD);
  
  // loop over X/Y bins
  TArrayD chEff(11);
  TArrayD chEffErr[2];
  chEffErr[0].Set(11);
  chEffErr[1].Set(11);
  for (Int_t ix = 1; ix <= nxBins; ++ix) {
    
    // set x range
    TT->GetAxis(xDim)->SetRange(ix, ix);
    TD->GetAxis(xDim)->SetRange(ix, ix);
    
    for (Int_t iy = 1; iy <= nyBins; ++iy) {
      
      // set y range
      TT->GetAxis(yDim)->SetRange(iy, iy);
      TD->GetAxis(yDim)->SetRange(iy, iy);
      
      // compute chamber efficency and errors
      if (GetChamberEfficiency(*TT, *TD, chEff, chEffErr, kFALSE)) {
        
        // compute the overall tracking efficiency
        TGraphAsymmErrors *dummy[3] = {0x0, 0x0, 0x0};
        GetTrackingEfficiency(chEff, chEffErr, dummy, dummy, 0, 0.);
        
        // fill histo
        effVsXY->Fill(TT->GetAxis(xDim)->GetBinCenter(ix),TT->GetAxis(yDim)->GetBinCenter(iy),chEff[0]);
        effVsXY->SetBinError(ix,iy,TMath::Max(chEffErr[0][0], chEffErr[1][0]));
        
      } else {
        
        // fill histo with 0 ± 1
	effVsXY->Fill(TT->GetAxis(xDim)->GetBinCenter(ix),TT->GetAxis(yDim)->GetBinCenter(iy),0.);
	effVsXY->SetBinError(ix,iy,1.);
        
      }
      
    }
    
  }
  
  // close input file
  file->Close();
  
  // display
  effVsXY->GetXaxis()->SetTitle(xVar.Data());
  effVsXY->GetXaxis()->CenterTitle(kTRUE);
  effVsXY->GetXaxis()->SetLabelFont(22);
  effVsXY->GetXaxis()->SetTitleFont(22);
  effVsXY->GetYaxis()->SetTitle(yVar.Data());
  effVsXY->GetYaxis()->CenterTitle(kTRUE);
  effVsXY->GetYaxis()->SetLabelFont(22);
  effVsXY->GetYaxis()->SetTitleFont(22);
  effVsXY->GetZaxis()->SetTitle("efficiency");
  effVsXY->GetZaxis()->SetLabelFont(22);
  effVsXY->GetZaxis()->SetTitleFont(22);
  
  if (draw) {
    new TCanvas(Form("cTrackingEffVs%s-%s",xVar.Data(),yVar.Data()), Form("Measured tracking efficiency versus %s and %s",xVar.Data(),yVar.Data()),700,600);
    effVsXY->DrawClone("surf1");
  }
  
  // save output
  file = new TFile(fileNameSave.Data(),"update");
  effVsXY->Write(0x0, TObject::kOverwrite);
  
  // add an histo with variable size rapidity bins
  if (yDim == 3 && rap) {
    
    TH2F* effVsXYrap = new TH2F();
    TString rapName = Form("trackingEffVs%s-%sRapBins", xVar.Data(), yVar.Data());
    TString rapTitle = Form("Measured tracking efficiency versus %s and %s", xVar.Data(), yVar.Data());
    effVsXYrap->SetTitle(rapTitle.Data());
    effVsXYrap->SetName(rapName.Data());
    
    Double_t xBinEdge[nxBins+1];
    for (Int_t xbin = 0; xbin <= nxBins; ++xbin)
      xBinEdge[xbin] = effVsXY->GetXaxis()->GetBinLowEdge(xbin+1);
    
    Double_t yBinEdge[nyBins+1];
    for (Int_t ybin = 0; ybin <= nyBins; ++ybin)
      yBinEdge[ybin] = 2*TMath::ATan(TMath::Exp((effVsXY->GetYaxis()->GetBinLowEdge(ybin+1))));
    
    effVsXYrap->SetBins(nxBins, xBinEdge, nyBins, yBinEdge);
    
    for (Int_t xbin = 1; xbin <= nxBins; ++xbin)
      for (Int_t ybin = 1; ybin <= nyBins; ++ybin)
        effVsXYrap->SetBinContent(xbin, ybin, effVsXY->GetBinContent(xbin,ybin));
    
    effVsXYrap->Write(0x0, TObject::kOverwrite);
    
    delete effVsXYrap;

  }
  
  file->Close();
  
  // clean memory
  delete effVsXY;
  
}


//---------------------------------------------------------------------------
void PlotMuonEfficiency(TString fileNameData, TString fileNameSave, Bool_t saveEdges, Bool_t print, Bool_t draw)
{
  /// plot chamber, station and overall tracking efficiency
  
  printf("plotting efficiency...\n");
  
  // get input hists
  TFile *file = new TFile(fileNameData.Data(), "read");
  if (!file || !file->IsOpen()) {
    printf("cannot open file %s \n",fileNameData.Data());
    return;
  }
  TList *listTT = static_cast<TList*>(file->FindObjectAny(Form("TotalTracksPerChamber%s", gObjNameExtension.Data())));
  TList *listTD = static_cast<TList*>(file->FindObjectAny(Form("TracksDetectedPerChamber%s", gObjNameExtension.Data())));
  THnSparse *TT = static_cast<THnSparse*>(listTT->At(10));
  THnSparse *TD = static_cast<THnSparse*>(listTD->At(10));
  
  // output graphs
  TGraphAsymmErrors *effVsCh[3] = {0x0, 0x0, 0x0};
  TGraphAsymmErrors *effVsSt[3] = {0x0, 0x0, 0x0};
  TString nameAdd[3] = {"", "Low", "Up"};
  TString titleAdd[3] = {"", " - lower limit", " - upper limit"};
  for (Int_t i = 0; i < 1 || (saveEdges && i < 3); ++i) {
    effVsCh[i] = CreateGraph(Form("chamberEff%s",nameAdd[i].Data()),
                             Form("Measured efficiency per chamber (0 = spectro)%s",titleAdd[i].Data()));
    effVsSt[i] = CreateGraph(Form("stationEff%s",nameAdd[i].Data()),
                             Form("Measured efficiency per station (6 = st4&5)%s",titleAdd[i].Data()));
  }
  
  // set the centrality and pT ranges for integration
  SetCentPtCh(*TT);
  SetCentPtCh(*TD);
  
  TArrayD chEff(11);
  TArrayD chEffErr[2];
  chEffErr[0].Set(11);
  chEffErr[1].Set(11);
  
  // compute chamber efficency and errors
  if (GetChamberEfficiency(*TT, *TD, chEff, chEffErr, print)) {
    
    // compute the overall tracking efficiency
    GetTrackingEfficiency(chEff, chEffErr, effVsSt, effVsCh, 0, 0., print);
    
  } else {
    
    // set tracking efficiency to 1 +0 -1
    effVsCh[0]->SetPoint(0,0.,1.);
    effVsCh[0]->SetPointError(0,0.,0.,1.,0.);
    
    for (Int_t iSt = 0; iSt < 6; ++iSt) {
      effVsSt[0]->SetPoint(iSt,iSt+1,1.);
      effVsSt[0]->SetPointError(iSt,0.,0.,1.,0.);
    }
    
    if (saveEdges) {
      
      // lower = 0 ± 0
      effVsCh[1]->SetPoint(0,0.,0.);
      effVsCh[1]->SetPointError(0,0.,0.,0.,0.);
      
      // upper = 1 ± 0
      effVsCh[2]->SetPoint(0,0.,1.);
      effVsCh[2]->SetPointError(0,0.,0.,0.,0.);
      
      for (Int_t iSt = 0; iSt < 6; ++iSt) {
        
        effVsSt[1]->SetPoint(iSt,iSt+1,0.);
        effVsSt[1]->SetPointError(iSt,0.,0.,0.,0.);
        
        effVsSt[2]->SetPoint(iSt,iSt+1,1.);
        effVsSt[2]->SetPointError(iSt,0.,0.,0.,0.);
        
      }
      
    }
    
  }
  
  // fill graph vs chamber
  for (Int_t iCh = 1; iCh < 11; iCh++) {
    for (Int_t i = 0; i < 1 || (saveEdges && i < 3); ++i) {
      effVsCh[i]->SetPoint(iCh,iCh,chEff[iCh]);
      effVsCh[i]->SetPointError(iCh,0.,0.,chEffErr[0][iCh],chEffErr[1][iCh]);
    }
  }
  
  // close input file
  file->Close();
  
  // display
  for (Int_t i = 0; i < 1 || (saveEdges && i < 3); ++i) {
    effVsCh[i]->GetXaxis()->Set(12, -0.5, 10.5);
    effVsCh[i]->GetXaxis()->SetNdivisions(11);
    BeautifyGraph(*effVsCh[i], "chamber", "efficiency");
    effVsSt[i]->GetXaxis()->Set(7, 0.5, 6.5);
    effVsSt[i]->GetXaxis()->SetNdivisions(6);
    BeautifyGraph(*effVsSt[i], "station", "efficiency");
  }
  
  if (draw) {
    TCanvas *c = new TCanvas("cEfficiency", "Measured tracking efficiency" , 1000, 400);
    c->Divide(2,1);
    gROOT->SetSelectedPad(c->cd(1));
    effVsCh[0]->DrawClone("ap");
    gROOT->SetSelectedPad(c->cd(2));
    effVsSt[0]->DrawClone("ap");
  }
  
  // save output
  file = new TFile(fileNameSave.Data(),"update");
  for (Int_t i = 0; i < 1 || (saveEdges && i < 3); ++i) effVsCh[i]->Write(0x0, TObject::kOverwrite);
  for (Int_t i = 0; i < 1 || (saveEdges && i < 3); ++i) effVsSt[i]->Write(0x0, TObject::kOverwrite);
  file->Close();
  
  // clean memory
  for (Int_t i = 0; i < 1 || (saveEdges && i < 3); ++i) {
    delete effVsCh[i];
    delete effVsSt[i];
  }
  
}


//---------------------------------------------------------------------------
void PlotMuonEfficiencyVsRun(TString runList, TString fileNameData, TString fileNameSave, Bool_t print, Bool_t draw)
{
  /// plot chamber, station and overall tracking efficiency versus run
  
  printf("plotting efficiency versus run...\n");
  
  // open run list
  ifstream inFile(runList.Data());
  if (!inFile.is_open()) {
    printf("cannot open file %s\n",runList.Data());
    return;
  }
  
  // output graphs
  TObjArray chamberVsRunGraphs; // 10 graphs: 1 for each chamber
  chamberVsRunGraphs.SetOwner(kTRUE);
  for ( Int_t iCh = 1; iCh < 11; ++iCh)
    chamberVsRunGraphs.Add(CreateGraph("effCh%dVsRun", "Measured efficiency for chamber %d versus run",iCh));
  
  TObjArray stationVsRunGraphs[3]; // 6 graphs: 1 for each station, and 1 for station 4&5
  TGraphAsymmErrors *trkVsRun[3];
  TString nameAdd[3] = {"", "Low", "Up"};
  TString titleAdd[3] = {"", " - lower limit", " - upper limit"};
  for (Int_t i = 0; i < 3; ++i) {
    
    stationVsRunGraphs[i].SetOwner(kTRUE);
    for ( Int_t iSt = 1; iSt < 6; ++iSt)
      stationVsRunGraphs[i].Add(CreateGraph(Form("effSt%%dVsRun%s",nameAdd[i].Data()),
                                            Form("Measured efficiency for station %%d versus run%s",titleAdd[i].Data()),iSt));
    stationVsRunGraphs[i].Add(CreateGraph(Form("effSt4&5VsRun%s",nameAdd[i].Data()),
                                          Form("Measured efficiency for station 4&5 versus run%s",titleAdd[i].Data())));
    
    trkVsRun[i] = CreateGraph(Form("trackingEffVsRun%s",nameAdd[i].Data()),
                              Form("Measured tracking efficiency versus run%s", titleAdd[i].Data()));
    
  }
  
  Int_t irun = -1;
  TList runs;
  runs.SetOwner();
  
  // loop over runs
  while (!inFile.eof()) {
    
    // get current run number
    TString currRun;
    currRun.ReadLine(inFile,kTRUE);
    if(currRun.IsNull()) continue;
    runs.AddLast(new TObjString(currRun));
    irun++;
    
    Int_t run = currRun.Atoi();
    
    printf("run %d: ", run);
    
    // compute efficiencies for this run
    TString dataFile = Form("runs/%d/%s", run, fileNameData.Data());
    TString outFile = Form("runs/%d/%s", run, fileNameSave.Data());
    PlotMuonEfficiency(dataFile, outFile, kTRUE, print, kFALSE);
    
    TFile *file = new TFile(outFile.Data(), "read");
    if (!file || !file->IsOpen()) {
      printf("cannot open file\n");
      continue;
    }
    
    // loop over central value and edges
    for (Int_t i = 0; i < 3; ++i) {
      
      // get results
      TGraphAsymmErrors *effVsCh;
      TGraphAsymmErrors *effVsSt;
      effVsCh = static_cast<TGraphAsymmErrors*>(file->FindObjectAny(Form("chamberEff%s",nameAdd[i].Data())));
      effVsSt = static_cast<TGraphAsymmErrors*>(file->FindObjectAny(Form("stationEff%s",nameAdd[i].Data())));
      if (!effVsCh || !effVsSt) {
        printf("efficiency graph not found\n");
        continue;
      }
      
      // fill chamber efficiency
      if (i == 0) for ( Int_t iCh = 0; iCh < 10; ++iCh) {
        TGraphAsymmErrors *g = static_cast<TGraphAsymmErrors*>(chamberVsRunGraphs.UncheckedAt(iCh));
        g->SetPoint(irun,irun,effVsCh->GetY()[iCh+1]);
        g->SetPointError(irun,0.,0.,effVsCh->GetErrorYlow(iCh+1),effVsCh->GetErrorYhigh(iCh+1));
      }
      
      // fill station efficiency
      for ( Int_t iSt = 0; iSt < 6; ++iSt) {
        TGraphAsymmErrors *g = static_cast<TGraphAsymmErrors*>(stationVsRunGraphs[i].UncheckedAt(iSt));
        g->SetPoint(irun,irun,effVsSt->GetY()[iSt]);
        g->SetPointError(irun,0.,0.,effVsSt->GetErrorYlow(iSt),effVsSt->GetErrorYhigh(iSt));
      }
      
      // fill spectrometer efficiency
      trkVsRun[i]->SetPoint(irun,irun,effVsCh->GetY()[0]);
      trkVsRun[i]->SetPointError(irun,0.,0.,effVsCh->GetErrorYlow(0),effVsCh->GetErrorYhigh(0));
      
    }
    
    file->Close();
    
  }
  
  inFile.close();
  
  // display
  BeautifyGraphs(chamberVsRunGraphs,"run number","efficiency");
  SetRunLabel(chamberVsRunGraphs,irun,runs);
  for (Int_t i = 0; i < 3; ++i) {
    BeautifyGraph(*trkVsRun[i],"run number","efficiency");
    SetRunLabel(*trkVsRun[i],irun,runs);
    BeautifyGraphs(stationVsRunGraphs[i],"run number","efficiency");
    SetRunLabel(stationVsRunGraphs[i],irun,runs);
  }
  
  if (draw) {
    new TCanvas("cTrackingEffVsRun", "Tracking efficiency versus run",1000,400);
    trkVsRun[0]->DrawClone("ap");
  }
  
  // save output
  TFile* file = new TFile(fileNameSave.Data(),"update");
  chamberVsRunGraphs.Write("ChamberEffVsRun", TObject::kOverwrite | TObject::kSingleKey);
  for (Int_t i = 0; i < 3; ++i)
    stationVsRunGraphs[i].Write(Form("StationEffVsRun%s",nameAdd[i].Data()), TObject::kOverwrite | TObject::kSingleKey);
  for (Int_t i = 0; i < 3; ++i) trkVsRun[i]->Write(0x0, TObject::kOverwrite);
  file->Close();
  
  // clean memory
  for (Int_t i = 0; i < 3; ++i) delete trkVsRun[i];
  
}


//---------------------------------------------------------------------------
void PlotIntegratedMuonEfficiency(TString fileNameWeights, TString fileNameSave, Bool_t print, Bool_t draw)
{
  /// plot chamber, station and overall tracking efficiency integrated over runs
  
  printf("plotting integrated efficiency...\n");
  
  // load run weights
  if (!fileNameWeights.IsNull()) LoadRunWeights(fileNameWeights);
  if (!runWeights) {
    printf("Cannot compute integrated efficiency without run-by-run weights\n");
    return;
  }
  
  // get input hists
  TFile *file = new TFile(fileNameSave.Data(), "update");
  if (!file || !file->IsOpen()) {
    printf("cannot open file\n");
    return;
  }
  TObjArray *chamberVsRunGraphs = static_cast<TObjArray*>(file->FindObjectAny("ChamberEffVsRun"));
  TObjArray *stationVsRunGraphs[2];
  stationVsRunGraphs[0] = static_cast<TObjArray*>(file->FindObjectAny("StationEffVsRunLow"));
  stationVsRunGraphs[1] = static_cast<TObjArray*>(file->FindObjectAny("StationEffVsRunUp"));
  TGraphAsymmErrors *trkVsRun[2];
  trkVsRun[0] = static_cast<TGraphAsymmErrors*>(file->FindObjectAny("trackingEffVsRunLow"));
  trkVsRun[1] = static_cast<TGraphAsymmErrors*>(file->FindObjectAny("trackingEffVsRunUp"));
  if (!chamberVsRunGraphs || !stationVsRunGraphs[0] || !stationVsRunGraphs[1] || !trkVsRun[0] || !trkVsRun[1]) {
    printf("object not found --> you must first plot the efficiency versus run\n");
    return;
  }
  
  // output graphs
  TGraphAsymmErrors *effVsCh = CreateGraph("integratedChamberEff", "Integrated efficiency per chamber (0 = spectro)");
  TGraphAsymmErrors *effVsSt = CreateGraph("integratedStationEff", "Integrated efficiency per station (6 = st4&5)");
  
  // integrate spectrometer efficiency
  IntegrateMuonEfficiency(*trkVsRun[0], *trkVsRun[1], *effVsCh, 0, 0.);
  
  // integrate chamber efficiency
  for ( Int_t iCh = 0; iCh < 10; ++iCh) {
    TGraphAsymmErrors *g = static_cast<TGraphAsymmErrors*>(chamberVsRunGraphs->UncheckedAt(iCh));
    IntegrateMuonEfficiency(*g, *g, *effVsCh, iCh+1, iCh+1);
  }
  
  // integrate station efficiency
  for ( Int_t iSt = 0; iSt < 6; ++iSt) {
    TGraphAsymmErrors *gLow = static_cast<TGraphAsymmErrors*>(stationVsRunGraphs[0]->UncheckedAt(iSt));
    TGraphAsymmErrors *gUp = static_cast<TGraphAsymmErrors*>(stationVsRunGraphs[1]->UncheckedAt(iSt));
    IntegrateMuonEfficiency(*gLow, *gUp, *effVsSt, iSt, iSt+1);
  }
  
  // print results
  if (print) {
    for (Int_t iCh = 1; iCh < 11; ++iCh) {
      cout << "Efficiency chamber " << iCh << " : ";
      cout << effVsCh->GetY()[iCh] << " + " << effVsCh->GetErrorYhigh(iCh) << " - " << effVsCh->GetErrorYlow(iCh) << endl;
    }
    for (Int_t iSt = 0; iSt < 6; ++iSt) {
      if (iSt < 5) cout << "Station " << iSt+1 << " = ";
      else cout << "Station 45 = ";
      cout << effVsSt->GetY()[iSt] << " + " << effVsSt->GetErrorYhigh(iSt) << " - " << effVsSt->GetErrorYlow(iSt) << endl;
    }
    cout << "Total tracking efficiency : ";
    cout << effVsCh->GetY()[0] << " + " << effVsCh->GetErrorYhigh(0) << " - " << effVsCh->GetErrorYlow(0) << endl << endl;
  }
  
  // display
  effVsCh->GetXaxis()->Set(12, -0.5, 10.5);
  effVsCh->GetXaxis()->SetNdivisions(11);
  BeautifyGraph(*effVsCh, "chamber", "efficiency");
  effVsSt->GetXaxis()->Set(7, 0.5, 6.5);
  effVsSt->GetXaxis()->SetNdivisions(6);
  BeautifyGraph(*effVsSt, "station", "efficiency");
  
  if (draw) {
    TCanvas *c = new TCanvas("cIntegratedEfficiency", "Integrated tracking efficiency" , 1000, 400);
    c->Divide(2,1);
    gROOT->SetSelectedPad(c->cd(1));
    effVsCh->DrawClone("ap");
    gROOT->SetSelectedPad(c->cd(2));
    effVsSt->DrawClone("ap");
  }
  
  // save output
  effVsCh->Write(0x0, TObject::kOverwrite);
  effVsSt->Write(0x0, TObject::kOverwrite);
  file->Close();
  
  // clean memory
  delete effVsCh;
  delete effVsSt;
  
}


//---------------------------------------------------------------------------
void PlotMuonEfficiencyPerDE(TString fileNameData, TString fileNameSave, Bool_t saveEdges)
{
  /// plot chamber and station efficiency per DE
  
  printf("plotting efficiency per DE...\n");
  
  // get input hists
  TFile *file = new TFile(fileNameData.Data(), "read");
  if (!file || !file->IsOpen()) {
    printf("cannot open file %s \n",fileNameData.Data());
    return;
  }
  TList *listTT = static_cast<TList*>(file->FindObjectAny(Form("TotalTracksPerChamber%s", gObjNameExtension.Data())));
  TList *listTD = static_cast<TList*>(file->FindObjectAny(Form("TracksDetectedPerChamber%s", gObjNameExtension.Data())));
  
  // output graph arrays
  TObjArray chamberVsDEGraphs; // 10 graphs: 1 for each chamber
  chamberVsDEGraphs.SetOwner(kTRUE);
  TGraphAsymmErrors *gCh[10]; // shortcut
  
  TObjArray stationVsDEGraphs[3]; // 6 graphs: 1 for each station, and 1 for station 4&5
  TString nameAdd[3] = {"", "Low", "Up"};
  TString titleAdd[3] = {"", " - lower limit", " - upper limit"};
  for (Int_t i = 0; i < 3; ++i) stationVsDEGraphs[i].SetOwner(kTRUE);
  TGraphAsymmErrors *gSt[6][3]; // shortcut
  for (Int_t iSt = 0; iSt < 6; ++iSt) for (Int_t i = 0; i < 3; ++i) gSt[iSt][i] = 0x0;
  
  // loop over chambers
  for (Int_t iCh = 0; iCh < 10; ++iCh) {
    
    // output graphs
    chamberVsDEGraphs.Add(CreateGraph("effCh%dVsDE","Measured efficiency for chamber %d per DE",iCh+1));
    gCh[iCh] = static_cast<TGraphAsymmErrors*>(chamberVsDEGraphs.UncheckedAt(iCh));
    
    // get input hists
    THnSparse *TT = static_cast<THnSparse*>(listTT->At(iCh));
    THnSparse *TD = static_cast<THnSparse*>(listTD->At(iCh));
    
    // set the centrality and pT range for integration
    SetCentPtCh(*TT);
    SetCentPtCh(*TD);
    
    // compute DE efficency and errors
    GetDEEfficiency(*TT, *TD, *gCh[iCh]);
    
  }
  
  // close input file
  file->Close();
  
  TArrayD chEff(11);
  TArrayD chEffErr[2];
  chEffErr[0].Set(11);
  chEffErr[1].Set(11);
  
  // loop over the first 3 stations
  for (Int_t iSt = 0; iSt < 3; ++iSt) {
    
    // output graphs
    for (Int_t i = 0; i < 1 || (saveEdges && i < 3); ++i) {
      stationVsDEGraphs[i].Add(CreateGraph(Form("effSt%%dVsDE%s",nameAdd[i].Data()),
                                           Form("Measured efficiency for station %%d per DE%s",titleAdd[i].Data()),iSt+1));
      gSt[iSt][i] = static_cast<TGraphAsymmErrors*>(stationVsDEGraphs[i].UncheckedAt(iSt));
    }
    
    // loop over DE
    Int_t nDE = gCh[2*iSt]->GetN();
    for (Int_t iDE = 0; iDE < nDE; ++iDE) {
      
      // copy efficiency
      for (Int_t iCh = 0; iCh < 2; ++iCh) {
        chEff[2*iSt+iCh+1] = gCh[2*iSt+iCh]->GetY()[iDE];
        chEffErr[0][2*iSt+iCh+1] = gCh[2*iSt+iCh]->GetErrorYlow(iDE);
        chEffErr[1][2*iSt+iCh+1] = gCh[2*iSt+iCh]->GetErrorYhigh(iDE);
      }
      
      // compute station efficiency
      GetStationEfficiency(chEff, chEffErr, iSt, gSt[iSt], iDE, iDE);
      
    }
    
  }
  
  // output graphs for last 2 stations
  for (Int_t i = 0; i < 1 || (saveEdges && i < 3); ++i) {
    for (Int_t iSt = 3; iSt < 5; ++iSt) {
      stationVsDEGraphs[i].Add(CreateGraph(Form("effSt%%dVsDE%s",nameAdd[i].Data()),
                                           Form("Measured efficiency for station %%d per DE%s",titleAdd[i].Data()),iSt+1));
      gSt[iSt][i] = static_cast<TGraphAsymmErrors*>(stationVsDEGraphs[i].UncheckedAt(iSt));
    }
    stationVsDEGraphs[i].Add(CreateGraph(Form("effSt4&5VsDE%s",nameAdd[i].Data()),
                                         Form("Measured efficiency for station 4&5 per DE%s",titleAdd[i].Data())));
    gSt[5][i] = static_cast<TGraphAsymmErrors*>(stationVsDEGraphs[i].UncheckedAt(5));
  }
  
  // loop over DE
  Int_t nDE = gCh[6]->GetN();
  for (Int_t iDE = 0; iDE < nDE; ++iDE) {
    
    // copy efficiency
    for (Int_t iCh = 6; iCh < 10; ++iCh) {
      chEff[iCh+1] = gCh[iCh]->GetY()[iDE];
      chEffErr[0][iCh+1] = gCh[iCh]->GetErrorYlow(iDE);
      chEffErr[1][iCh+1] = gCh[iCh]->GetErrorYhigh(iDE);
    }
    
    // compute station 4&5 efficiency individually
    for (Int_t iSt = 3; iSt < 5; ++iSt) GetStationEfficiency(chEff, chEffErr, iSt, gSt[iSt], iDE, iDE);
    
    // compute station 4&5 efficiency together
    GetStation45Efficiency(chEff, chEffErr, gSt[5], iDE, iDE);
    
  }
  
  // display
  for (Int_t iCh = 0; iCh < 10; ++iCh) {
    Int_t nDE = gCh[iCh]->GetN();
    gCh[iCh]->GetXaxis()->Set(nDE+1, -0.5, nDE-0.5);
    gCh[iCh]->GetXaxis()->SetNdivisions(nDE);
  }
  BeautifyGraphs(chamberVsDEGraphs,"Detection Element","efficiency");
  for (Int_t i = 0; i < 1 || (saveEdges && i < 3); ++i) {
    for (Int_t iSt = 0; iSt < 6; ++iSt) {
      Int_t nDE = gSt[iSt][i]->GetN();
      gSt[iSt][i]->GetXaxis()->Set(nDE+1, -0.5, nDE-0.5);
      gSt[iSt][i]->GetXaxis()->SetNdivisions(nDE);
    }
    BeautifyGraphs(stationVsDEGraphs[i],"Detection Element","efficiency");
    
  }
  
  // save Output
  file = new TFile(fileNameSave.Data(),"update");
  chamberVsDEGraphs.Write("ChamberEffVsDE", TObject::kOverwrite | TObject::kSingleKey);
  for (Int_t i = 0; i < 1 || (saveEdges && i < 3); ++i)
    stationVsDEGraphs[i].Write(Form("StationEffVsDE%s",nameAdd[i].Data()), TObject::kOverwrite | TObject::kSingleKey);
  file->Close();
  
}


//---------------------------------------------------------------------------
void PlotMuonEfficiencyPerDEVsRun(TString runList, TString fileNameData, TString fileNameSave)
{
  /// plot chamber and station efficiency per DE versus run
  
  printf("plotting efficiency per DE versus run...\n");
  
  // open run list
  ifstream inFile(runList.Data());
  if (!inFile.is_open()) {
    printf("cannot open file %s\n",runList.Data());
    return;
  }
  
  // output graphs
  TObjArray deVsRunGraphs; // 1 graph per DE
  deVsRunGraphs.SetOwner(kTRUE);
  TObjArray stDEVsRunGraphs[3]; // 1 graph per pair (quartet) of DE in individual stations (stations 4&5)
  TString nameAdd[3] = {"", "Low", "Up"};
  TString titleAdd[3] = {"", " - lower limit", " - upper limit"};
  for (Int_t i = 0; i < 3; ++i) stDEVsRunGraphs[i].SetOwner(kTRUE);
  Bool_t createGraph = kTRUE;
  
  Int_t irun = -1;
  TList runs;
  runs.SetOwner();
  
  // loop over runs
  while (!inFile.eof()) {
    
    // get current run number
    TString currRun;
    currRun.ReadLine(inFile,kTRUE);
    if(currRun.IsNull()) continue;
    runs.AddLast(new TObjString(currRun));
    irun++;
    
    Int_t run = currRun.Atoi();
    
    printf("run %d: ", run);
    
    // compute efficiencies for this run
    TString dataFile = Form("runs/%d/%s", run, fileNameData.Data());
    TString outFile = Form("runs/%d/%s", run, fileNameSave.Data());
    PlotMuonEfficiencyPerDE(dataFile, outFile, kTRUE);
    
    TFile *file = new TFile(outFile.Data(), "read");
    if (!file || !file->IsOpen()) {
      printf("cannot open file\n");
      continue;
    }
    
    // get results
    TObjArray *chamberVsDEGraphs = static_cast<TObjArray*>(file->FindObjectAny("ChamberEffVsDE"));
    if (!chamberVsDEGraphs) {
      printf("efficiency graph not found\n");
      continue;
    }
    
    Int_t currentDE = 0;
    
    // loop over chambers
    for ( Int_t iCh = 0; iCh < 10; ++iCh) {
      
      // input graph
      TGraphAsymmErrors *g = static_cast<TGraphAsymmErrors*>(chamberVsDEGraphs->UncheckedAt(iCh));
      Int_t nDE = g->GetN();
      
      // loop over DE
      for (Int_t iDE = 0; iDE < nDE; ++iDE) {
        
        // output graph
        if (createGraph) deVsRunGraphs.Add(CreateGraph("effDE%dVsRun","Measured efficiency for DE %d versus run",100*(iCh+1)+iDE));
        TGraphAsymmErrors *gDE= static_cast<TGraphAsymmErrors*>(deVsRunGraphs.UncheckedAt(currentDE++));
        
        // fill DE efficiency
        gDE->SetPoint(irun,irun,g->GetY()[iDE]);
        gDE->SetPointError(irun,0.,0.,g->GetErrorYlow(iDE),g->GetErrorYhigh(iDE));
        
      }
      
    }
    
    // loop over central value and edges
    for (Int_t i = 0; i < 3; ++i) {
      
      // get results
      TObjArray *stationVsDEGraphs = static_cast<TObjArray*>(file->FindObjectAny(Form("StationEffVsDE%s",nameAdd[i].Data())));
      if (!stationVsDEGraphs) {
        printf("efficiency graph not found\n");
        continue;
      }
      
      Int_t currentStDE = 0;
      
      // loop over stations
      for ( Int_t iSt = 0; iSt < 6; ++iSt) {
        
        // input graph
        TGraphAsymmErrors *g = static_cast<TGraphAsymmErrors*>(stationVsDEGraphs->UncheckedAt(iSt));
        Int_t nDE = g->GetN();
        
        // loop over DE
        for (Int_t iDE = 0; iDE < nDE; ++iDE) {
          
          // output graph
          if (createGraph) {
            TString sSt = (iSt<5) ? Form("%d",iSt+1) : "4&5";
            stDEVsRunGraphs[i].Add(CreateGraph(Form("effSt%sDE%%dVsRun%s",sSt.Data(),nameAdd[i].Data()),
                                               Form("Measured efficiency for DE %%d in station %s versus run%s",sSt.Data(),titleAdd[i].Data()),iDE));
          }
          TGraphAsymmErrors *gDE= static_cast<TGraphAsymmErrors*>(stDEVsRunGraphs[i].UncheckedAt(currentStDE++));
          
          // fill DE efficiency
          gDE->SetPoint(irun,irun,g->GetY()[iDE]);
          gDE->SetPointError(irun,0.,0.,g->GetErrorYlow(iDE),g->GetErrorYhigh(iDE));
          
        }
        
      }
      
    }
      
    file->Close();
    
    createGraph = kFALSE;
    
  }
  
  inFile.close();
  
  // display
  BeautifyGraphs(deVsRunGraphs,"run number","efficiency");
  SetRunLabel(deVsRunGraphs,irun,runs);
  for (Int_t i = 0; i < 3; ++i) {
    BeautifyGraphs(stDEVsRunGraphs[i],"run number","efficiency");
    SetRunLabel(stDEVsRunGraphs[i],irun,runs);
  }
  
  // save output
  TFile* file = new TFile(fileNameSave.Data(),"update");
  deVsRunGraphs.Write("DEEffVsRun", TObject::kOverwrite | TObject::kSingleKey);
  for (Int_t i = 0; i < 3; ++i)
    stDEVsRunGraphs[i].Write(Form("DEEffPerStationVsRun%s",nameAdd[i].Data()), TObject::kOverwrite | TObject::kSingleKey);
  file->Close();
  
}


//---------------------------------------------------------------------------
void PlotIntegratedMuonEfficiencyPerDE(TString fileNameWeights, TString fileNameSave)
{
  /// plot chamber and station efficiency per DE integrated over runs
  
  printf("plotting integrated efficiency per DE...\n");
  
  // load run weights
  if (!fileNameWeights.IsNull()) LoadRunWeights(fileNameWeights);
  if (!runWeights) {
    printf("Cannot compute integrated efficiency without run-by-run weights\n");
    return;
  }
  
  // get input hists
  TFile *file = new TFile(fileNameSave.Data(), "update");
  if (!file || !file->IsOpen()) {
    printf("cannot open file\n");
    return;
  }
  TObjArray *deVsRunGraphs = static_cast<TObjArray*>(file->FindObjectAny("DEEffVsRun"));
  TObjArray *stDEVsRunGraphs[2];
  stDEVsRunGraphs[0] = static_cast<TObjArray*>(file->FindObjectAny("DEEffPerStationVsRunLow"));
  stDEVsRunGraphs[1] = static_cast<TObjArray*>(file->FindObjectAny("DEEffPerStationVsRunUp"));
  if (!deVsRunGraphs || !stDEVsRunGraphs[0] || !stDEVsRunGraphs[1]) {
    printf("object not found --> you must first plot the efficiency versus run\n");
    return;
  }
  
  // output graph arrays
  TObjArray chamberVsDEGraphs; // 10 graphs: 1 for each chamber
  chamberVsDEGraphs.SetOwner(kTRUE);
  for ( Int_t iCh = 0; iCh < 10; ++iCh)
    chamberVsDEGraphs.Add(CreateGraph("integratedEffCh%dVsDE","Integrated efficiency for chamber %d per DE",iCh+1));
  
  TObjArray stationVsDEGraphs; // 6 graphs: 1 for each station, and 1 for station 4&5
  stationVsDEGraphs.SetOwner(kTRUE);
  for ( Int_t iSt = 0; iSt < 5; ++iSt)
    stationVsDEGraphs.Add(CreateGraph("integratedEffSt%dVsDE","Integrated efficiency for station %d per DE",iSt+1));
  stationVsDEGraphs.Add(CreateGraph("integratedEffSt4&5VsDE","Integrated efficiency for station 4&5 per DE"));
  
  // Loop over DE
  TIter nextDE(deVsRunGraphs);
  TGraphAsymmErrors *gDE = 0x0;
  while ((gDE = static_cast<TGraphAsymmErrors*>(nextDE()))) {
    
    // get chamber and DE indices
    Int_t deId;
    sscanf(gDE->GetName(), "effDE%dVsRun", &deId);
    Int_t iCh = deId/100-1;
    Int_t iDE = deId%100;
    
    // integrate DE efficiency
    TGraphAsymmErrors *g = static_cast<TGraphAsymmErrors*>(chamberVsDEGraphs.UncheckedAt(iCh));
    IntegrateMuonEfficiency(*gDE, *gDE, *g, iDE, iDE);
    
  }
  
  // Loop over DE per station
  Int_t ng = stDEVsRunGraphs[0]->GetEntries();
  for (Int_t ig = 0; ig < ng; ++ig) {
    
    TGraphAsymmErrors *gDELow = static_cast<TGraphAsymmErrors*>(stDEVsRunGraphs[0]->UncheckedAt(ig));
    TGraphAsymmErrors *gDEUp = static_cast<TGraphAsymmErrors*>(stDEVsRunGraphs[1]->UncheckedAt(ig));
    
    // get station and DE indices
    Int_t iSt, iDE;
    if (strstr(gDELow->GetName(), "4&5")) {
      iSt = 5;
      sscanf(gDELow->GetName(), "effSt4&5DE%dVsRunLow", &iDE);
    } else {
      sscanf(gDELow->GetName(), "effSt%dDE%dVsRunLow", &iSt, &iDE);
      iSt--;
    }
    
    // Integrate DE efficiency per station
    TGraphAsymmErrors *g = static_cast<TGraphAsymmErrors*>(stationVsDEGraphs.UncheckedAt(iSt));
    IntegrateMuonEfficiency(*gDELow, *gDEUp, *g, iDE, iDE);
    
  }
  
  // display
  for ( Int_t iCh = 0; iCh < 10; ++iCh) {
    TGraphAsymmErrors *g = static_cast<TGraphAsymmErrors*>(chamberVsDEGraphs.UncheckedAt(iCh));
    Int_t nDE = g->GetN();
    g->GetXaxis()->Set(nDE+1, -0.5, nDE-0.5);
    g->GetXaxis()->SetNdivisions(nDE);
  }
  BeautifyGraphs(chamberVsDEGraphs,"Detection Element","efficiency");
  for ( Int_t iSt = 0; iSt < 6; ++iSt) {
    TGraphAsymmErrors *g = static_cast<TGraphAsymmErrors*>(stationVsDEGraphs.UncheckedAt(iSt));
    Int_t nDE = g->GetN();
    g->GetXaxis()->Set(nDE+1, -0.5, nDE-0.5);
    g->GetXaxis()->SetNdivisions(nDE);
  }
  BeautifyGraphs(stationVsDEGraphs,"Detection Element","efficiency");

  // save Output
  chamberVsDEGraphs.Write("IntegratedChamberEffVsDE", TObject::kOverwrite | TObject::kSingleKey);
  stationVsDEGraphs.Write("IntegratedStationEffVsDE", TObject::kOverwrite | TObject::kSingleKey);
  file->Close();
  
}


//---------------------------------------------------------------------------
Bool_t GetChamberEfficiency(THnSparse &TT, THnSparse &TD, TArrayD &chEff, TArrayD chEffErr[2], Bool_t printError)
{
  /// compute chamber efficiency and errors
  /// return kFALSE if efficiency unknown for all chambers
  
  // project track hists over the chamber axis
  TH1D *TTdraw = TT.Projection(0,"e");
  TH1D *TDdraw = TD.Projection(0,"e");
  
  // compute chamber efficiency and errors
  TGraphAsymmErrors *efficiency = (TTdraw->GetEntries() > 0) ? new TGraphAsymmErrors(TDdraw, TTdraw,Form("%se0",effErrMode)) : 0x0;
  Bool_t ok = (efficiency);
  
  // fill arrays
  if (ok) {
    
    Bool_t missingEff = kFALSE;
    
    for (Int_t i = 0; i < 10; i++) {
      
      if (TTdraw->GetBinContent(i+1) > 0) {
        
        chEff[i+1] = efficiency->GetY()[i];
        chEffErr[0][i+1] = efficiency->GetErrorYlow(i);
        chEffErr[1][i+1] = efficiency->GetErrorYhigh(i);
        
      } else {
        
        chEff[i+1] = -1.;
        chEffErr[0][i+1] = 0.;
        chEffErr[1][i+1] = 0.;
        
        missingEff = kTRUE;
        
      }
      
    }
    
    if (missingEff && printError) cout << "efficiency partially unknown" << endl;
    
  } else {
    
    for (Int_t i = 0; i < 10; i++) {
      
      chEff[i+1] = -1.;
      chEffErr[0][i+1] = 0.;
      chEffErr[1][i+1] = 0.;
      
    }
    
    if (printError) cout << "efficiency unknown" << endl << endl;
    
  }
  
  // clean memory
  delete TTdraw;
  delete TDdraw;
  delete efficiency;
  
  return ok;
  
}


//---------------------------------------------------------------------------
void GetDEEfficiency(THnSparse &TT, THnSparse &TD, TGraphAsymmErrors &effVsDE)
{
  /// compute DE efficiency and errors
  
  // project track hists over the chamber axis
  TH1D *TTdraw = TT.Projection(0,"e");
  TH1D *TDdraw = TD.Projection(0,"e");
  
  // compute DE efficiency and errors
  TGraphAsymmErrors *efficiency = (TTdraw->GetEntries() > 0) ? new TGraphAsymmErrors(TDdraw, TTdraw,Form("%se0",effErrMode)) : 0x0;
  Int_t nDE = TTdraw->GetNbinsX();
  
  if (efficiency) {
    
    for (Int_t iDE = 0; iDE < nDE; ++iDE) {
      
      if (TTdraw->GetBinContent(iDE+1) > 0) {
        
        effVsDE.SetPoint(iDE,iDE,efficiency->GetY()[iDE]);
        effVsDE.SetPointError(iDE,0,0,efficiency->GetErrorYlow(iDE),efficiency->GetErrorYhigh(iDE));
        
      } else {
        
        effVsDE.SetPoint(iDE,iDE,-1);
        effVsDE.SetPointError(iDE,0,0,0,0);
        
      }
      
    }
    
  } else {
    
    for (Int_t iDE = 0; iDE < nDE; ++iDE) {
      
      effVsDE.SetPoint(iDE,iDE,-1);
      effVsDE.SetPointError(iDE,0,0,0,0);
      
    }
    
  }
  
  // clean memory
  delete TTdraw;
  delete TDdraw;
  delete efficiency;
  
}


//---------------------------------------------------------------------------
void ComputeStationEfficiency(TArrayD &chEff, TArrayD chEffErr[2], Int_t iSt, Double_t &stEff, Double_t stEffErr[2])
{
  /// compute the station iSt (0...4) efficiency and errors from the individual chamber efficiencies and errors
  
  stEff = 1.-(1.-chEff[2*iSt+1])*(1.-chEff[2*iSt+2]);
  
  Double_t d1 = (1. - chEff[2*iSt+2]); d1 *= d1;
  Double_t d2 = (1. - chEff[2*iSt+1]); d2 *= d2;
  
  for (Int_t i = 0; i < 2; ++i) {
    Double_t s1 = chEffErr[i][2*iSt+1] * chEffErr[i][2*iSt+1];
    Double_t s2 = chEffErr[i][2*iSt+2] * chEffErr[i][2*iSt+2];
    stEffErr[i] = TMath::Sqrt(d1*s1 + d2*s2 + s1*s2);
  }
  
  stEffErr[0] = TMath::Min(stEff, stEffErr[0]);
  stEffErr[1] = TMath::Min(1.-stEff, stEffErr[1]);
  
}


//---------------------------------------------------------------------------
void GetStationEfficiency(TArrayD &chEff, TArrayD chEffErr[2], Int_t iSt, TGraphAsymmErrors *effVsX[3], Int_t ip, Double_t xp)
{
  /// compute the station iSt (0...4) efficiency and errors from the individual chamber efficiencies and errors
  
  if (chEff[2*iSt+1] >= 0 && chEff[2*iSt+2] >= 0) {
    
    // compute station efficiency from known chamber efficiency
    Double_t stEff, stEffErr[2];
    ComputeStationEfficiency(chEff, chEffErr, iSt, stEff, stEffErr);
    
    // fill graphs if required
    for (Int_t i = 0; i < 3; ++i) {
      if (effVsX[i]) {
        effVsX[i]->SetPoint(ip,xp,stEff);
        effVsX[i]->SetPointError(ip,0.,0.,stEffErr[0],stEffErr[1]);
      }
    }
    
  } else {
    
    Double_t edge[2] = {0., 1.};
    TArrayD chEffEdge[2];
    Double_t stEffEdge[2];
    Double_t stEffEdgeErr[2][2];
    
    for (Int_t i = 0; i < 2; ++i) {
      
      // set lower(upper) limit of chamber efficiency
      chEffEdge[i].Set(11);
      for (Int_t iCh = 1; iCh < 3; ++iCh) chEffEdge[i][2*iSt+iCh] = (chEff[2*iSt+iCh] < 0) ? edge[i] : chEff[2*iSt+iCh];
      
      // compute station efficiency
      ComputeStationEfficiency(chEffEdge[i], chEffErr, iSt, stEffEdge[i], stEffEdgeErr[i]);
      
      // fill graph if required
      if (effVsX[i+1]) {
        effVsX[i+1]->SetPoint(ip,xp,stEffEdge[i]);
        effVsX[i+1]->SetPointError(ip,0.,0.,stEffEdgeErr[i][0],stEffEdgeErr[i][1]);
      }
      
    }
    
    // combine extreme cases to get station efficiency and fill graph if required
    if (effVsX[0]) {
      effVsX[0]->SetPoint(ip,xp,stEffEdge[1]);
      effVsX[0]->SetPointError(ip,0.,0.,stEffEdge[1]-stEffEdge[0]+stEffEdgeErr[0][0],stEffEdgeErr[1][1]);
    }
    
  }
  
}


//---------------------------------------------------------------------------
void ComputeStation45Efficiency(TArrayD &chEff, TArrayD chEffErr[2], Double_t &st45Eff, Double_t st45EffErr[2])
{
  /// compute the station 4-5 efficiency and errors from chamber efficiencies and errors
  
  st45Eff = chEff[7]*chEff[8]*chEff[9] + chEff[7]*chEff[8]*chEff[10] + chEff[7]*chEff[9]*chEff[10] + chEff[8]*chEff[9]*chEff[10] - 3.*chEff[7]*chEff[8]*chEff[9]*chEff[10];
  
  Double_t d1 = chEff[8]*chEff[9] + chEff[8]*chEff[10] + chEff[9]*chEff[10] - 3.*chEff[8]*chEff[9]*chEff[10]; d1 *= d1;
  Double_t d2 = chEff[7]*chEff[9] + chEff[7]*chEff[10] + chEff[9]*chEff[10] - 3.*chEff[7]*chEff[9]*chEff[10]; d2 *= d2;
  Double_t d3 = chEff[7]*chEff[8] + chEff[7]*chEff[10] + chEff[8]*chEff[10] - 3.*chEff[7]*chEff[8]*chEff[10]; d3 *= d3;
  Double_t d4 = chEff[7]*chEff[8] + chEff[7]*chEff[9] + chEff[8]*chEff[9] - 3.*chEff[7]*chEff[8]*chEff[9]; d4 *= d4;
  Double_t d12 = chEff[9] + chEff[10] - 3.*chEff[9]*chEff[10]; d12 *= d12;
  Double_t d13 = chEff[8] + chEff[10] - 3.*chEff[8]*chEff[10]; d13 *= d13;
  Double_t d14 = chEff[8] + chEff[9] - 3.*chEff[8]*chEff[9]; d14 *= d14;
  Double_t d23 = chEff[7] + chEff[10] - 3.*chEff[7]*chEff[10]; d23 *= d23;
  Double_t d24 = chEff[7] + chEff[9] - 3.*chEff[7]*chEff[9]; d24 *= d24;
  Double_t d34 = chEff[7] + chEff[8] - 3.*chEff[7]*chEff[8]; d34 *= d34;
  Double_t d123 = 1. - 3.*chEff[10]; d123 *= d123;
  Double_t d124 = 1. - 3.*chEff[9]; d124 *= d124;
  Double_t d134 = 1. - 3.*chEff[8]; d134 *= d134;
  Double_t d234 = 1. - 3.*chEff[7]; d234 *= d234;
  Double_t d1234 = - 3.; d1234 *= d1234;
  
  for (Int_t i = 0; i < 2; ++i) {
    Double_t s1 = chEffErr[i][7] * chEffErr[i][7];
    Double_t s2 = chEffErr[i][8] * chEffErr[i][8];
    Double_t s3 = chEffErr[i][9] * chEffErr[i][9];
    Double_t s4 = chEffErr[i][10] * chEffErr[i][10];
    st45EffErr[i] = TMath::Sqrt(d1*s1 + d2*s2 + d3*s3 + d4*s4 + d12*s1*s2 + d13*s1*s3 + d14*s1*s4 + d23*s2*s3 + d24*s2*s4 + d34*s3*s4 + d123*s1*s2*s3 + d124*s1*s2*s4 + d134*s1*s3*s4 + d234*s2*s3*s4 + d1234*s1*s2*s3*s4);
  }
  
  st45EffErr[0] = TMath::Min(st45Eff, st45EffErr[0]);
  st45EffErr[1] = TMath::Min(1.-st45Eff, st45EffErr[1]);
  
}


//---------------------------------------------------------------------------
void GetStation45Efficiency(TArrayD &chEff, TArrayD chEffErr[2], TGraphAsymmErrors *effVsX[3], Int_t ip, Double_t xp)
{
  /// compute the station 4-5 efficiency and errors from chamber efficiencies and errors
  
  if (chEff[7] >= 0 && chEff[8] >= 0 && chEff[9] >= 0 && chEff[10] >= 0) {
    
    // compute station efficiency from known chamber efficiency
    Double_t stEff, stEffErr[2];
    ComputeStation45Efficiency(chEff, chEffErr, stEff, stEffErr);
    
    // fill graphs if required
    for (Int_t i = 0; i < 3; ++i) {
      if (effVsX[i]) {
        effVsX[i]->SetPoint(ip,xp,stEff);
        effVsX[i]->SetPointError(ip,0.,0.,stEffErr[0],stEffErr[1]);
      }
    }
    
  } else {
    
    Double_t edge[2] = {0., 1.};
    TArrayD chEffEdge[2];
    Double_t stEffEdge[2];
    Double_t stEffEdgeErr[2][2];
    
    for (Int_t i = 0; i < 2; ++i) {
      
      // set lower(upper) limit of chamber efficiency
      chEffEdge[i].Set(11);
      for (Int_t iCh = 7; iCh < 11; ++iCh) chEffEdge[i][iCh] = (chEff[iCh] < 0) ? edge[i] : chEff[iCh];
      
      // compute station efficiency
      ComputeStation45Efficiency(chEffEdge[i], chEffErr, stEffEdge[i], stEffEdgeErr[i]);
      
      // fill graph if required
      if (effVsX[i+1]) {
        effVsX[i+1]->SetPoint(ip,xp,stEffEdge[i]);
        effVsX[i+1]->SetPointError(ip,0.,0.,stEffEdgeErr[i][0],stEffEdgeErr[i][1]);
      }
      
    }
    
    // combine extreme cases to get station efficiency and fill graph if required
    if (effVsX[0]) {
      effVsX[0]->SetPoint(ip,xp,stEffEdge[1]);
      effVsX[0]->SetPointError(ip,0.,0.,stEffEdge[1]-stEffEdge[0]+stEffEdgeErr[0][0],stEffEdgeErr[1][1]);
    }
    
  }
  
}


//---------------------------------------------------------------------------
void ComputeTrackingEfficiency(Double_t stEff[6], Double_t stEffErr[6][2], Double_t &spectroEff, Double_t spectroEffErr[2])
{
  /// compute the spectrometer efficiency and errors from the station efficiencies and errors
  
  Double_t de[6][2];
  for (Int_t iSt = 0; iSt < 6; iSt++) de[iSt][0] = stEff[iSt]*stEff[iSt];
  
  if (moreTrackCandidates) {
    
    spectroEff = stEff[0] * stEff[1] * stEff[2] * stEff[3] * stEff[4];
    
    for (Int_t i = 0; i < 2; i++) {
      
      for (Int_t iSt = 0; iSt < 6; iSt++) de[iSt][1] = stEffErr[iSt][i]*stEffErr[iSt][i];
      
      spectroEffErr[i] = 0.;
      for (Int_t j = 1; j < 32; j++) {
	Double_t sigmaAdd = 1.;
	for (Int_t iSt = 0; iSt < 5; iSt++) sigmaAdd *= de[iSt][TESTBIT(j,iSt)];
	spectroEffErr[i] += sigmaAdd;
      }
      spectroEffErr[i] = TMath::Sqrt(spectroEffErr[i]);
      
    }
    
  } else {
    
    spectroEff = stEff[0] * stEff[1] * stEff[2] * stEff[5];
    
    for (Int_t i = 0; i < 2; i++) {
      
      for (Int_t iSt = 0; iSt < 6; iSt++) de[iSt][1] = stEffErr[iSt][i]*stEffErr[iSt][i];
      
      spectroEffErr[i] = 0.;
      for (Int_t j = 1; j < 16; j++) {
	Double_t sigmaAdd = de[5][TESTBIT(j,3)];
	for (Int_t iSt = 0; iSt < 3; iSt++) sigmaAdd *= de[iSt][TESTBIT(j,iSt)];
	spectroEffErr[i] += sigmaAdd;
      }
      spectroEffErr[i] = TMath::Sqrt(spectroEffErr[i]);
      
    }
    
  }
  
  spectroEffErr[0] = TMath::Min(spectroEff, spectroEffErr[0]);
  spectroEffErr[1] = TMath::Min(1.-spectroEff, spectroEffErr[1]);
  
}


//---------------------------------------------------------------------------
void GetTrackingEfficiency(TArrayD &chEff, TArrayD chEffErr[2], TGraphAsymmErrors *effVsSt[3],
                           TGraphAsymmErrors *effVsX[3], Int_t ip, Double_t xp, Bool_t print)
{
  /// compute the tracking efficiency from the individual chamber efficiencies
  
  Double_t stEff[6];
  Double_t stEffErr[6][2];
  
  // check if unknown chamber efficiency
  Bool_t allEffKnown = kTRUE;
  for (Int_t iCh = 1; iCh < 11 && allEffKnown; ++iCh) if (chEff[iCh] < 0) allEffKnown = kFALSE;
  
  if (allEffKnown) {
    
    // compute station efficiency
    for (Int_t iSt = 0; iSt < 5; ++iSt) ComputeStationEfficiency(chEff, chEffErr, iSt, stEff[iSt], stEffErr[iSt]);
    ComputeStation45Efficiency(chEff, chEffErr, stEff[5], stEffErr[5]);
    
    // compute spectrometer efficiency
    Double_t spectroEff, spectroEffErr[2];
    ComputeTrackingEfficiency(stEff, stEffErr, spectroEff, spectroEffErr);
    chEff[0] = spectroEff;
    chEffErr[0][0] = spectroEffErr[0];
    chEffErr[1][0] = spectroEffErr[1];
    
    // fill graphs if required
    for (Int_t i = 0; i < 3; ++i) {
      if (effVsX[i]) {
        effVsX[i]->SetPoint(ip,xp,chEff[0]);
        effVsX[i]->SetPointError(ip,0.,0.,chEffErr[0][0],chEffErr[1][0]);
      }
      if (effVsSt[i]) for (Int_t iSt = 0; iSt < 6; ++iSt) {
        effVsSt[i]->SetPoint(iSt,iSt+1,stEff[iSt]);
        effVsSt[i]->SetPointError(iSt,0.,0.,stEffErr[iSt][0],stEffErr[iSt][1]);
      }
    }
    
  } else {
    
    Double_t edge[2] = {0., 1.};
    TArrayD chEffEdge[2];
    Double_t stEffEdge[2][6];
    Double_t stEffEdgeErr[2][6][2];
    Double_t spectroEffEdge[2], spectroEffEdgeErr[2][2];
    
    for (Int_t i = 0; i < 2; ++i) {
      
      // set lower(upper) limit of chamber efficiency
      chEffEdge[i].Set(11);
      for (Int_t iCh = 1; iCh < 11; ++iCh) chEffEdge[i][iCh] = (chEff[iCh] < 0) ? edge[i] : chEff[iCh];
      
      // compute station efficiency
      for (Int_t iSt = 0; iSt < 5; ++iSt) ComputeStationEfficiency(chEffEdge[i], chEffErr, iSt, stEffEdge[i][iSt], stEffEdgeErr[i][iSt]);
      ComputeStation45Efficiency(chEffEdge[i], chEffErr, stEffEdge[i][5], stEffEdgeErr[i][5]);
      
      // compute spectrometer efficiency
      ComputeTrackingEfficiency(stEffEdge[i], stEffEdgeErr[i], spectroEffEdge[i], spectroEffEdgeErr[i]);
      
      // fill graph if required
      if (effVsX[i+1]) {
        effVsX[i+1]->SetPoint(ip,xp,spectroEffEdge[i]);
        effVsX[i+1]->SetPointError(ip,0.,0.,spectroEffEdgeErr[i][0],spectroEffEdgeErr[i][1]);
      }
      if (effVsSt[i+1]) for (Int_t iSt = 0; iSt < 6; ++iSt) {
        effVsSt[i+1]->SetPoint(iSt,iSt+1,stEffEdge[i][iSt]);
        effVsSt[i+1]->SetPointError(iSt,0.,0.,stEffEdgeErr[i][iSt][0],stEffEdgeErr[i][iSt][1]);
      }
      
    }
    
    // combine extreme cases to get station efficiency
    for (Int_t iSt = 0; iSt < 6; ++iSt) {
      stEff[iSt] = stEffEdge[1][iSt];
      stEffErr[iSt][0] = stEffEdge[1][iSt] - stEffEdge[0][iSt] + stEffEdgeErr[0][iSt][0];
      stEffErr[iSt][1] = stEffEdgeErr[1][iSt][1];
    }
    
    // combine extreme cases to get spectrometer efficiency
    chEff[0] = spectroEffEdge[1];
    chEffErr[0][0] = spectroEffEdge[1] - spectroEffEdge[0] + spectroEffEdgeErr[0][0];
    chEffErr[1][0] = spectroEffEdgeErr[1][1];
    
    // fill graph if required
    if (effVsX[0]) {
      effVsX[0]->SetPoint(ip,xp,chEff[0]);
      effVsX[0]->SetPointError(ip,0.,0.,chEffErr[0][0],chEffErr[1][0]);
    }
    if (effVsSt[0]) for (Int_t iSt = 0; iSt < 6; ++iSt) {
      effVsSt[0]->SetPoint(iSt,iSt+1,stEff[iSt]);
      effVsSt[0]->SetPointError(iSt,0.,0.,stEffErr[iSt][0],stEffErr[iSt][1]);
    }
    
  }
  
  // print results
  if (print) {
    for (Int_t iCh = 1; iCh < 11; ++iCh) {
      cout << "Efficiency chamber " << iCh << " : ";
      cout << chEff[iCh] << " + " << chEffErr[1][iCh] << " - " << chEffErr[0][iCh] << endl;
    }
    for (Int_t iSt = 0; iSt < 6; ++iSt) {
      if (iSt < 5) cout << "Station " << iSt+1 << " = ";
      else cout << "Station 45 = ";
      cout << stEff[iSt] << " + " << stEffErr[iSt][1] << " - " << stEffErr[iSt][0] << endl;
    }
    cout << "Total tracking efficiency : ";
    cout << chEff[0] << " + " << chEffErr[1][0] << " - " << chEffErr[0][0] << endl << endl;
  }
  
}


//---------------------------------------------------------------------------
void IntegrateMuonEfficiency(TGraphAsymmErrors &effVsRunLow, TGraphAsymmErrors &effVsRunUp,
                             TGraphAsymmErrors &effVsX, Int_t ip, Double_t xp)
{
  /// integrate efficiency over runs
  /// return kFALSE if efficiency unknown in all runs
  
  if (!runWeights) {
    effVsX.SetPoint(ip,xp,-1.);
    effVsX.SetPointError(ip,0.,0.,0.,0.);
    return;
  }
  
  // initialize
  Double_t rec[2] = {0., 0.};
  Double_t gen = 0.;
  Double_t intEffErr[2] = {0., 0.};
  Bool_t ok = kFALSE;
  
  // loop over runs
  Int_t nRuns = effVsRunLow.GetN();
  for (Int_t iRun = 0; iRun < nRuns; ++iRun) {
    
    // get run weight
    TString sRun = effVsRunLow.GetXaxis()->GetBinLabel(iRun+1);
    TParameter<Double_t> *weight = static_cast<TParameter<Double_t>*>(runWeights->FindObject(sRun.Data()));
    if (!weight) {
      printf("weight not found for run %s\n", sRun.Data());
      continue;
    }
    Double_t w = weight->GetVal();
    Double_t w2 = w*w;
    
    // get efficiency and error
    Double_t eff[2] = {effVsRunLow.GetY()[iRun], effVsRunUp.GetY()[iRun]};
    Double_t effErr[2] = {effVsRunLow.GetErrorYlow(iRun), effVsRunUp.GetErrorYhigh(iRun)};
    if (eff[0] < 0. || eff[1] < 0.) {
      printf("no efficiency measurement --> use 0(1) ± 0 as lower(upper) limit\n");
      eff[0] = 0.;
      eff[1] = 1.;
      effErr[0] = 0.;
      effErr[1] = 0.;
    } else ok = kTRUE;
    
    // integrate
    gen += w;
    for (Int_t i = 0; i < 2; ++i) {
      rec[i] += w*eff[i];
      intEffErr[i] += w2*effErr[i]*effErr[i];
    }
    
  }
  
  // compute efficiency
  if (gen > 0. && ok) {
    
    effVsX.SetPoint(ip,xp,rec[1]/gen);
    effVsX.SetPointError(ip,0.,0.,(rec[1]-rec[0]+TMath::Sqrt(intEffErr[0]))/gen,TMath::Sqrt(intEffErr[1])/gen);
    
  } else {
    
    if (gen <= 0.) printf("impossible to integrate, all weights = 0 or unknown ?!?\n");
    else printf("efficiency never measured --> return -1 ± 0\n");
    
    effVsX.SetPoint(ip,xp,-1.);
    effVsX.SetPointError(ip,0.,0.,0.,0.);
    
  }
  
}


//---------------------------------------------------------------------------
void LoadRunWeights(TString fileName)
{
  /// Set the number of interested events per run
  /// (used to weight the acc*eff correction integrated
  /// over run for any pt/y/centrality bins)
  
  if (runWeights) return;
  
  ifstream inFile(fileName.Data());
  if (!inFile.is_open()) {
    printf("cannot open file %s\n", fileName.Data());
    return;
  }
  
  runWeights = new THashList(1000);
  runWeights->SetOwner();
  
  TString line;
  while (! inFile.eof() ) {
    
    line.ReadLine(inFile,kTRUE);
    if(line.IsNull()) continue;
    
    TObjArray *param = line.Tokenize(" ");
    if (param->GetEntries() != 2) {
      printf("bad input line %s", line.Data());
      continue;
    }
    
    Int_t run = ((TObjString*)param->UncheckedAt(0))->String().Atoi();
    if (run < 0) {
      printf("invalid run number: %d", run);
      continue;
    }
    
    Float_t weight = ((TObjString*)param->UncheckedAt(1))->String().Atof();
    if (weight <= 0.) {
      printf("invalid weight: %g", weight);
      continue;
    }
    
    runWeights->Add(new TParameter<Double_t>(Form("%d",run), weight));
    
    delete param;
  }
  
  inFile.close();
  
}


//---------------------------------------------------------------------------
void SetCentPtCh(THnSparse& SparseData)
{
  /// Sets the centrality, pT and y range and the charge for integration
  /// If ptMax < 0, it sets ptMax as the maximum pT value on the THnSparse
  /// If yMin < -4, it sets yMin as the minimum y value on the THnSparse
  /// If yMax > -2.5, it sets yMax as the maximum y value on the THnSparse
  /// If phiMin < 0, it sets phiMin as the minimum phi value on the THnSparse
  /// If phiMax > 2π, it sets phiMax as the maximum phi value on the THnSparse
  
  if (ptMax >= 0 && ptMax < ptMin) {
    printf("Wrong pT range, ptMax must be higher than ptMin\n");
    exit(-1);
  }
  
  if (yMin > -2.5 || yMax < -4. || yMax < yMin) {
    printf("Wrong y range\n");
    exit(-1);
  }
  
  if (phiMin > TMath::TwoPi() || phiMax < 0. || phiMax < phiMin) {
    printf("Wrong phi range\n");
    exit(-1);
  }
  
  if (charge < -1 || charge > 1) {
    printf("Selected charge must be 0, 1 or -1\n");
    exit(-1);
  }
  
  // adjust centrality range
  centMax = TMath::Max(centMax-1.e-12, centMin);
  
  // set the centrality range for integration
  Int_t lowBin = SparseData.GetAxis(1)->FindBin(centMin);
  Int_t upBin = SparseData.GetAxis(1)->FindBin(centMax);
  SparseData.GetAxis(1)->SetRange(lowBin, upBin);
  
  // set the pt range for integration
  if (ptMin < 0.) lowBin = 0;
  else lowBin = SparseData.GetAxis(2)->FindBin(ptMin);
  if (ptMax < 0.) upBin = SparseData.GetAxis(2)->GetNbins()+1;
  else upBin = SparseData.GetAxis(2)->FindBin(ptMax);
  SparseData.GetAxis(2)->SetRange(lowBin, upBin);
  
  // set the y range for integration
  if (yMin < -4.) lowBin = 0;
  else lowBin = SparseData.GetAxis(3)->FindBin(yMin);
  if (yMax > -2.5) upBin = SparseData.GetAxis(3)->GetNbins()+1;
  else upBin = SparseData.GetAxis(3)->FindBin(yMax);
  SparseData.GetAxis(3)->SetRange(lowBin, upBin);
  
  // set the phi range for integration
  if (phiMin < 0.) lowBin = 0;
  else lowBin = SparseData.GetAxis(4)->FindBin(phiMin);
  if (phiMax > TMath::TwoPi()) upBin = SparseData.GetAxis(4)->GetNbins()+1;
  else upBin = SparseData.GetAxis(4)->FindBin(phiMax);
  SparseData.GetAxis(4)->SetRange(lowBin, upBin);
  
  // set the charge range
  if (charge != 0) {
    lowBin = SparseData.GetAxis(5)->FindBin(charge);
    SparseData.GetAxis(5)->SetRange(lowBin, lowBin);
  }
  
}


//---------------------------------------------------------------------------
TGraphAsymmErrors* CreateGraph(const char* name, const char* title, int value)
{
  /// create a graph of a given name and title
  
  TGraphAsymmErrors* g = new TGraphAsymmErrors();
  
  if ( value >= 0 ) {
    
    g->SetName(Form(name,value));
    g->SetTitle(Form(title,value));
    
  } else {
    
    g->SetName(name);
    g->SetTitle(title);
    
    
  }
  
  return g;
  
}


//---------------------------------------------------------------------------
void BeautifyGraph(TGraphAsymmErrors &g, const char* xAxisName, const char* yAxisName)
{
  /// beautify this graph
  
  g.SetLineStyle(1);
  g.SetLineColor(1);
  g.SetMarkerStyle(20);
  g.SetMarkerSize(0.7);
  g.SetMarkerColor(2);
  
  g.GetXaxis()->SetTitle(xAxisName);
  g.GetXaxis()->CenterTitle(kTRUE);
  g.GetXaxis()->SetLabelFont(22);
  g.GetXaxis()->SetTitleFont(22);
  
  g.GetYaxis()->SetTitle(yAxisName);
  g.GetYaxis()->CenterTitle(kTRUE);
  g.GetYaxis()->SetLabelFont(22);
  g.GetYaxis()->SetTitleFont(22);
  
}


//---------------------------------------------------------------------------
void BeautifyGraphs(TObjArray& array, const char* xAxisName, const char* yAxisName)
{
  /// beautify all the graphs in this array
  
  for ( Int_t i = 0; i <= array.GetLast(); ++i)
    BeautifyGraph(*static_cast<TGraphAsymmErrors*>(array.UncheckedAt(i)), xAxisName, yAxisName);
  
}


//---------------------------------------------------------------------------
void SetRunLabel(TGraphAsymmErrors &g, Int_t irun, const TList& runs)
{
  /// set the run label for all the graphs
  
  g.GetXaxis()->Set(irun+1, -0.5, irun+0.5);
  
  TIter nextRun(&runs);
  TObjString *srun = 0x0;
  Int_t iirun = 1;
  while ((srun = static_cast<TObjString*>(nextRun()))) {
    g.GetXaxis()->SetBinLabel(iirun, srun->GetName());
    iirun++;
  }
  
}


//---------------------------------------------------------------------------
void SetRunLabel(TObjArray& array, Int_t irun, const TList& runs)
{
  /// set the run label for all the graphs
  
  for (Int_t i = 0; i <= array.GetLast(); ++i)
    SetRunLabel(*static_cast<TGraphAsymmErrors*>(array.UncheckedAt(i)), irun, runs);
  
}


