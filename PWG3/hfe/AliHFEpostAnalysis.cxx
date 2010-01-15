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
// Post analysis code
// Drawing nice pictures containing
//  - Efficiency
//  - Signal/Background
//  - PID Performance
//  More Post Analysis code will be added in time
//
//  Autor: Markus Fasel
//
#include <TAxis.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TList.h>
#include "AliCFContainer.h"
#include "AliCFEffGrid.h"

#include "AliHFEcuts.h"
#include "AliHFEpostAnalysis.h"

ClassImp(AliHFEpostAnalysis)

//____________________________________________________________
AliHFEpostAnalysis::AliHFEpostAnalysis():
  TObject(),
  fResults(NULL),
  fAnalysisObjects(0),
  fEfficiencyContainer(NULL),
  fPIDperformance(NULL),
  fSignalToBackgroundMC(NULL)
{
  //
  // Default Constructor
  //
}

//____________________________________________________________
AliHFEpostAnalysis::AliHFEpostAnalysis(const AliHFEpostAnalysis &ref):
  TObject(ref),
  fResults(ref.fResults),
  fAnalysisObjects(ref.fAnalysisObjects),
  fEfficiencyContainer(ref.fEfficiencyContainer),
  fPIDperformance(ref.fPIDperformance),
  fSignalToBackgroundMC(ref.fSignalToBackgroundMC)
{
  //
  // Copy Constructor
  //
}

//____________________________________________________________
AliHFEpostAnalysis& AliHFEpostAnalysis::operator=(const AliHFEpostAnalysis &ref){
  //
  // Assignment Operator
  //
  TObject::operator=(ref);
  fResults = ref.fResults;
  fAnalysisObjects = ref.fAnalysisObjects;
  fPIDperformance = ref.fPIDperformance;
  fSignalToBackgroundMC = ref.fSignalToBackgroundMC;

  return *this;
}

//____________________________________________________________
AliHFEpostAnalysis::~AliHFEpostAnalysis(){
  //
  // Do not delete objects where we are not owner
  //
  if(fResults) delete fResults;
}

//____________________________________________________________
Int_t AliHFEpostAnalysis::SetResults(TList *input){
  //
  // Publish the results to the post analysis
  //
  Int_t nFound = 0;
  fEfficiencyContainer = dynamic_cast<AliCFContainer *>(input->FindObject("container"));
  if(!fEfficiencyContainer){
    AliError("Efficiency Correction Framework Container not found in the list of outputs");
  } else {
    SETBIT(fAnalysisObjects, kCFC);
    nFound++;
  }
  fPIDperformance = dynamic_cast<THnSparseF *>(input->FindObject("PIDperformance"));
  if(!fPIDperformance){
    AliError("Histogram fPIDperformance not found in the List of Outputs");
  } else {
    SETBIT(fAnalysisObjects, kPIDperf);
    nFound++;
  }
  fSignalToBackgroundMC = dynamic_cast<THnSparseF *>(input->FindObject("SignalToBackgroundMC"));
  if(!fSignalToBackgroundMC){
    AliError("Histogram fSignalToBackgroundMC not found in the list of outputs");
  } else {
    SETBIT(fAnalysisObjects, kSigBackg);
    nFound++;
  }
  AliInfo(Form("Found %d analysis objects", nFound));
  return nFound;
}

//____________________________________________________________
void AliHFEpostAnalysis::StoreOutput(const char *filename){
  //
  // Save the results produced in a rootfile
  //
  if(fResults){
    TFile *outfile = new TFile(filename, "RECREATE");
    outfile->cd();
    fResults->Write("HFEresults", TObject::kSingleKey);
    outfile->Close();
    delete outfile;
  }
}

//____________________________________________________________
void AliHFEpostAnalysis::DrawMCSignal2Background(){
  //
  // Draw the MC signal/background plots
  //
  if(!fSignalToBackgroundMC) return;

  // First Select everything within the first ITS Layer
  fSignalToBackgroundMC->GetAxis(4)->SetRange(2,2);
  // For Signal electrons we project axis 3 to everything > 0
  fSignalToBackgroundMC->GetAxis(3)->SetRange(2,3);
  TH1 *hSignal = fSignalToBackgroundMC->Projection(0);
  hSignal->SetName("hSignal");
  hSignal->SetTitle("Signal Electrons");
  // For Background studies project axis 3 to bin 0
  fSignalToBackgroundMC->GetAxis(3)->SetRange(1,1);
  TH1 *hBackground = fSignalToBackgroundMC->Projection(0); 
  hBackground->SetName("hBackground");
  hBackground->SetTitle("Background Electrons");
  // For All electrons we undo the projection of axis 3
  fSignalToBackgroundMC->GetAxis(3)->SetRange(0, fSignalToBackgroundMC->GetAxis(3)->GetNbins());
  TH1 *hAll = fSignalToBackgroundMC->Projection(0);
  hAll->SetName("hAll");
  hAll->SetTitle("All Electrons");
  // Prepare Canvas
  TCanvas *cEff = new TCanvas("cEff", "MC Sig/Backg studies", 800, 400);
  cEff->Divide(2);
  // Plot Efficiency
  TH1 *hEff = (TH1 *)hSignal->Clone();
  hEff->Divide(hAll);
  hEff->SetName("sigEff");
  hEff->SetTitle("Signal/(Signal+Background)");
  hEff->GetXaxis()->SetTitle("p_{T} / GeV/c");
  hEff->GetYaxis()->SetTitle("Efficiency");
  hEff->SetStats(kFALSE);
  hEff->SetMarkerStyle(22);
  hEff->SetMarkerColor(kBlue);
  cEff->cd(1);
  hEff->Draw("p");
  // Plot Signal/Background
  TH1 *hSB = (TH1 *)hSignal->Clone();
  hSB->Divide(hBackground);
  hSB->SetName("sigEff");
  hSB->SetTitle("Signal/Background");
  hSB->GetXaxis()->SetTitle("p_{T} / GeV/c");
  hSB->GetYaxis()->SetTitle("Signal/Background");
  hSB->SetStats(kFALSE);
  hSB->SetMarkerStyle(22);
  hSB->SetMarkerColor(kBlue);
  cEff->cd(2);
  hSB->Draw("p");

  // Undo projections
  fSignalToBackgroundMC->GetAxis(4)->SetRange(0, fSignalToBackgroundMC->GetAxis(4)->GetNbins());
}

//____________________________________________________________
void AliHFEpostAnalysis::DrawEfficiency(){
  //
  // Draw the Efficiency
  // We show: 
  // + InAcceptance / Generated
  // + Signal / Generated
  // + Selected / Generated
  // + Selected / InAcceptance (Reconstructible)
  //
  TCanvas *cEff = new TCanvas("cEff", "Efficiency", 800, 600);
  cEff->Divide(2,2);
  if(!fEfficiencyContainer) return;
  AliCFEffGrid *effCalc = new AliCFEffGrid("effCalc", "Efficiency Calculation Grid", *fEfficiencyContainer);
  effCalc->CalculateEfficiency(AliHFEcuts::kStepMCInAcceptance, AliHFEcuts::kStepMCGenerated);
  TH1 *effReconstructibleP = effCalc->Project(0);
  effReconstructibleP->SetName("effReconstructibleP");
  effReconstructibleP->SetTitle("Efficiency of reconstructible tracks");
  effReconstructibleP->GetXaxis()->SetTitle("p_{T} / GeV/c");
  effReconstructibleP->GetYaxis()->SetTitle("Efficiency");
  effReconstructibleP->GetYaxis()->SetRangeUser(0.,1.);
  effReconstructibleP->SetMarkerStyle(22);
  effReconstructibleP->SetMarkerColor(kBlue);
  effReconstructibleP->SetLineColor(kBlack);
  effReconstructibleP->SetStats(kFALSE);
  cEff->cd(1);
  effReconstructibleP->Draw("e");
  effCalc->CalculateEfficiency(AliHFEcuts::kStepMCsignal, AliHFEcuts::kStepMCGenerated);
  TH1 *effSignal = effCalc->Project(0);
  effSignal->SetName("effSignal");
  effSignal->SetTitle("Efficiency of Signal Electrons");
  effSignal->GetXaxis()->SetTitle("p_{T} / GeV/c");
  effSignal->GetYaxis()->SetTitle("Efficiency");
  effSignal->GetYaxis()->SetRangeUser(0., 1.);
  effSignal->SetMarkerStyle(22);
  effSignal->SetMarkerColor(kBlue);
  effSignal->SetLineColor(kBlack);
  effSignal->SetStats(kFALSE);
  cEff->cd(2);
  effSignal->Draw("e");
  effCalc->CalculateEfficiency(AliHFEcuts::kStepHFEcutsTRD + 1, AliHFEcuts::kStepMCGenerated);
  TH1 *effPIDP = effCalc->Project(0);
  effPIDP->SetName("effPIDP");
  effPIDP->SetTitle("Efficiency of selected tracks");
  effPIDP->GetXaxis()->SetTitle("p_{T} / GeV/c");
  effPIDP->GetYaxis()->SetTitle("Efficiency");
  effPIDP->GetYaxis()->SetRangeUser(0.,1.);
  effPIDP->SetMarkerStyle(22);
  effPIDP->SetMarkerColor(kBlue);
  effPIDP->SetLineColor(kBlack);
  effPIDP->SetStats(kFALSE);
  cEff->cd(3);
  effPIDP->Draw("e");
  effCalc->CalculateEfficiency(AliHFEcuts::kStepHFEcutsTRD + 1, AliHFEcuts::kStepMCInAcceptance);
  TH1 *effPIDAcc = effCalc->Project(0);
  effPIDAcc->SetName("effPIDAcc");
  effPIDAcc->SetTitle("Efficiency of selected tracks in acceptance");
  effPIDAcc->GetXaxis()->SetTitle("p_{T} / GeV/c");
  effPIDAcc->GetYaxis()->SetTitle("Efficiency");
  effPIDAcc->GetYaxis()->SetRangeUser(0.,1.);
  effPIDAcc->SetMarkerStyle(22);
  effPIDAcc->SetMarkerColor(kBlue);
  effPIDAcc->SetLineColor(kBlack);
  effPIDAcc->SetStats(kFALSE);
  cEff->cd(4);
  effPIDAcc->Draw("e");
  delete effCalc;
}

//____________________________________________________________
void AliHFEpostAnalysis::DrawPIDperformance(){
  //
  // Plotting Ratio histograms
  // + All electrons / all candidates (Purity for Electrons)
  // + All signal electrons / all electrons (Purity for signals)
  // For this the following pt-histograms have to be projected from the THnSparse
  // + All Electron candidates
  // + All Real electrons
  // + All Signal Electrons
  // + All misidentified electrons
  // Additionally we draw Efficiency histograms:
  // + Reconstructible Electrons
  // + Signal Electrons
  // + Selected Electrons
  // + Selected Electrons in Acceptance
  //

  if(!fPIDperformance) return;
  // Make projection
  // always project to pt dimension
  // get the histograms under our control
  TH1 *allCandidates = NULL, *allElectrons = NULL, *allSignals = NULL, *allFakes = NULL;
  allCandidates = fPIDperformance->Projection(0);
  allCandidates->SetName("hAllCandidates");
  allCandidates->SetTitle("All Candidates");
  allCandidates->Sumw2();
  Int_t firstDim3 = fPIDperformance->GetAxis(3)->GetFirst(), lastDim3 = fPIDperformance->GetAxis(3)->GetLast();
  fPIDperformance->GetAxis(3)->SetRange(firstDim3 + 1, lastDim3);
  allElectrons = fPIDperformance->Projection(0);
  allElectrons->Sumw2();
  allElectrons->SetName("hAllElectrons");
  allElectrons->SetTitle("All Electrons");
  Int_t firstDim4 = fPIDperformance->GetAxis(4)->GetFirst(), lastDim4 = fPIDperformance->GetAxis(4)->GetLast();
  fPIDperformance->GetAxis(4)->SetRange(firstDim4 + 1, lastDim4);
  allSignals = fPIDperformance->Projection(0);
  allSignals->Sumw2();
  allSignals->SetName("hAllSignals");
  allSignals->SetTitle("All Signal Electrons");
  fPIDperformance->GetAxis(4)->SetRange(firstDim4, lastDim4);       // Reset 4th axis
  fPIDperformance->GetAxis(3)->SetRange(firstDim3, firstDim3);  // Select fakes
  allFakes = fPIDperformance->Projection(0);
  allFakes->Sumw2();
  allFakes->SetName("hAllFakes");
  allFakes->SetTitle("All Fakes");
  fPIDperformance->GetAxis(3)->SetRange(firstDim3, lastDim3);       // Reset also 3rd axis

  // Make Ratios
  TH1D *electronPurity = dynamic_cast<TH1D *>(allElectrons->Clone());
  electronPurity->Divide(allCandidates);
  electronPurity->SetName("electronPurity");
  electronPurity->SetTitle("Electron Purity");
  electronPurity->GetXaxis()->SetTitle("p_{T} / GeV/c");
  electronPurity->GetYaxis()->SetTitle("Purity / %");
  TH1D *signalPurity = dynamic_cast<TH1D *>(allSignals->Clone());
  signalPurity->Divide(allElectrons);
  signalPurity->SetName("signalPurity");
  signalPurity->SetTitle("Purity of Electrons coming from Heavy flavours");
  signalPurity->GetXaxis()->SetTitle("p_{T} / GeV/c");
  signalPurity->GetYaxis()->SetTitle("Purity / %");
  TH1D *fakeContamination = dynamic_cast<TH1D *>(allFakes->Clone());
  fakeContamination->Divide(allCandidates);
  fakeContamination->SetName("fakeContamination");
  fakeContamination->SetTitle("Contamination of misidentified hadrons");
  fakeContamination->GetXaxis()->SetTitle("p_{T} / GeV/c");
  fakeContamination->GetYaxis()->SetTitle("Purity / %");

  // Graphics settings
  const Int_t nHistos = 7;
  TH1 *hptr[7] = {allCandidates, allElectrons, allSignals, allFakes, electronPurity, signalPurity, fakeContamination}, *histo;
  for(Int_t ihist = 0; ihist < nHistos; ihist++){
    (histo = hptr[ihist])->SetStats(kFALSE);
    histo->SetLineColor(kBlack);
    histo->SetMarkerColor(kBlue);
    histo->SetMarkerStyle(22);
  }

  // Make percent output
  electronPurity->Scale(100);
  signalPurity->Scale(100);
  fakeContamination->Scale(100);

  // Draw output
  TCanvas *cSpectra = new TCanvas("cSpectra","pt Spectra", 800, 600);
  cSpectra->Divide(2,2);
  TCanvas *cRatios = new TCanvas("cRatios", "Ratio Plots", 800, 600);
  cRatios->Divide(2,2);
  cSpectra->cd(1);
  gPad->SetLogy();
  allCandidates->GetXaxis()->SetTitle("p_{T} / GeV/c");
  allCandidates->GetYaxis()->SetTitle("dN/dp_{T} 1/GeV/c");
  allCandidates->Draw("e");
  cSpectra->cd(2);
  gPad->SetLogy();
  allElectrons->GetXaxis()->SetTitle("p_{T} / GeV/c");
  allElectrons->GetYaxis()->SetTitle("dN/dp_{T} 1/GeV/c");
  allElectrons->Draw("e");
  cSpectra->cd(3);
  gPad->SetLogy();
  allSignals->GetXaxis()->SetTitle("p_{T} / GeV/c");
  allSignals->GetYaxis()->SetTitle("dN/dp_{T} 1/GeV/c");
  allSignals->Draw("e");
  cSpectra->cd(4);
  gPad->SetLogy();
  allFakes->GetXaxis()->SetTitle("p_{T} / GeV/c");
  allFakes->GetYaxis()->SetTitle("dN/dp_{T} 1/GeV/c");
  allFakes->Draw("e");
  cRatios->cd(1);
  electronPurity->Draw("e");
  cRatios->cd(2);
  signalPurity->Draw("e");
  cRatios->cd(3);
  fakeContamination->Draw("e");
}

