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

#include "AliHFEcontainer.h"
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
Int_t AliHFEpostAnalysis::SetTaskQA(const TList *input){
  //
  // Publish the results to the post analysis
  //
  Int_t nFound = 0;
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
  fSignalToBackgroundMC->GetAxis(5)->SetRange(2,2);
  TH1 *hEff[3], *hSB[3];
  // Select for different charge
  hEff[0] = CreateHistoSignalToBackgroundMC(0, 0);
  hEff[1] = CreateHistoSignalToBackgroundMC(0, 1);
  hEff[2] = CreateHistoSignalToBackgroundMC(0, 2);

  hSB[0] = CreateHistoSignalToBackgroundMC(1, 0);
  hSB[1] = CreateHistoSignalToBackgroundMC(1, 1);
  hSB[2] = CreateHistoSignalToBackgroundMC(1, 2);
 
  // Undo projections
  fSignalToBackgroundMC->GetAxis(5)->SetRange(0, fSignalToBackgroundMC->GetAxis(4)->GetNbins());

  // Prepare Canvas
  TCanvas *cMCSB = new TCanvas("cMCSB", "MC Sig/Backg studies", 800, 400);
  cMCSB->Divide(2);
  TLegend *leg;
  TH1 **sample[2] = {&hEff[0], &hSB[0]};
  const char *chargename[3] = {"All Charges", "Negative Charge", "Positive Charge"};
  for(Int_t isample = 0; isample < 2; isample++){
    leg = new TLegend(0.7, 0.1, 0.89, 0.3);
    leg->SetBorderSize(1); 
    leg->SetFillColor(kWhite);
    cMCSB->cd(isample + 1);
    for(Int_t icharge = 0; icharge < 3; icharge++){
      sample[isample][icharge]->Draw(icharge > 0?  "psame" : "p");
      leg->AddEntry(sample[isample][icharge], chargename[icharge], "p");
    }
    leg->Draw();
    gPad->Update();
  }
}

//____________________________________________________________
TH1 *AliHFEpostAnalysis::CreateHistoSignalToBackgroundMC(Int_t mode, Int_t charge){ 
  //
  // Make Efficiency / SB histograms for different charges
  //
  TH1 *hNom = NULL, *hDenom = NULL;
  if(charge) fSignalToBackgroundMC->GetAxis(3)->SetRange(charge, charge);
  // For Signal electrons we project axis 4 to everything > 0
  fSignalToBackgroundMC->GetAxis(4)->SetRange(2,3);
  hNom = fSignalToBackgroundMC->Projection(0);
  hNom->SetName("nom");
  fSignalToBackgroundMC->GetAxis(4)->SetRange(0, fSignalToBackgroundMC->GetAxis(4)->GetLast() + 1);
  if(mode == 1) fSignalToBackgroundMC->GetAxis(4)->SetRange(1,1);
  hDenom =  fSignalToBackgroundMC->Projection(0);
  hDenom->SetName("denom");
  if(mode == 1)  fSignalToBackgroundMC->GetAxis(4)->SetRange(0, fSignalToBackgroundMC->GetAxis(4)->GetLast() + 1);
  if(charge) fSignalToBackgroundMC->GetAxis(3)->SetRange(0, fSignalToBackgroundMC->GetAxis(3)->GetLast() + 1);

  TH1 *hEff = dynamic_cast<TH1D *>(hNom->Clone());
  if(hEff){
    TString hname, cname;
    hname = mode ? "sigToBack" : "sigEff";
    Color_t mycolor = kBlack;
    switch(charge){
      case 0: mycolor = kBlue; cname = "All"; break;
      case 1: mycolor = kBlack; cname = "Neg"; break;
      case 2: mycolor = kRed; cname ="Pos"; break;
      default: break;
    }
    hname += cname;
    hEff->SetName(hname);
    hEff->SetTitle(mode ? "Signal/Background" : "Signal/(Signal+Background)");
    hEff->Divide(hDenom);

    // Make nice plots
    hEff->GetXaxis()->SetTitle("p_{T} / GeV/c");
    hEff->GetYaxis()->SetTitle("Efficiency");
    hEff->SetStats(kFALSE);
    hEff->SetLineColor(kBlack);
    hEff->SetLineWidth(1);
    hEff->SetMarkerStyle(22);
    hEff->SetMarkerColor(mycolor);
  }

  delete hNom; delete hDenom;
  return hEff;
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
  AliCFContainer *tracks = fEfficiencyContainer->MakeMergedCFContainer("trackContCombined", "MC + Rec(reco) Track Information", "MCTrackCont:recTrackContReco");
  AliCFEffGrid *effCalc = new AliCFEffGrid("effCalc", "Efficiency Calculation Grid", *tracks);
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
  effCalc->CalculateEfficiency(AliHFEcuts::kStepMCGeneratedZOutNoPileUpCentralityFine, AliHFEcuts::kStepMCGenerated);
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
  effCalc->CalculateEfficiency(tracks->GetNStep() - 1, AliHFEcuts::kStepMCGenerated);
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
  effCalc->CalculateEfficiency(tracks->GetNStep() - 1, AliHFEcuts::kStepMCInAcceptance);
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
  //

  if(!fPIDperformance) return;
  // Make projection
  TH1 *electronPurity[3], *signalPurity[3], *fakeContamination[3];
  electronPurity[0] = CreateHistoPIDperformance(0, 0);
  electronPurity[1] = CreateHistoPIDperformance(0, 1);
  electronPurity[2] = CreateHistoPIDperformance(0, 2);

  signalPurity[0] = CreateHistoPIDperformance(1, 0);
  signalPurity[1] = CreateHistoPIDperformance(1, 1);
  signalPurity[2] = CreateHistoPIDperformance(1, 2);

  fakeContamination[0] = CreateHistoPIDperformance(2, 0);
  fakeContamination[1] = CreateHistoPIDperformance(2, 1);
  fakeContamination[2] = CreateHistoPIDperformance(2, 2);

  // Draw output
  TCanvas *cRatios = new TCanvas("cRatios", "Ratio Plots", 800, 600);
  const char *chargename[3] = {"All Charges", "Negative Charge", "Positive Charge"};
  cRatios->Divide(2,2);
  TH1 **sample[3] = {&electronPurity[0], &signalPurity[0], &fakeContamination[0]};
  TLegend *leg;
  for(Int_t isample = 0; isample < 3; isample++){
    cRatios->cd(isample + 1);
    leg = new TLegend(0.7, 0.1, 0.89, 0.3);
    leg->SetBorderSize(1);
    leg->SetFillColor(kWhite);
    for(Int_t icharge = 0; icharge < 3; icharge++){
      leg->AddEntry(sample[isample][icharge], chargename[icharge], "p");
      sample[isample][icharge]->Draw(icharge > 0 ? "esame" : "e");
    }
    leg->Draw();
    gPad->Update();
  }
}

//____________________________________________________________
TH1 *AliHFEpostAnalysis::CreateHistoPIDperformance(Int_t mode, Int_t charge){
  //
  // Make Histograms for PID performance plots
  //
  fPIDperformance->GetAxis(4)->SetRange(0, fPIDperformance->GetAxis(4)->GetNbins()+1);
  fPIDperformance->GetAxis(3)->SetRange(0, fPIDperformance->GetAxis(3)->GetNbins() + 1);

  TH1 *hNom = NULL, *hDenom = NULL;
  TString hname, htitle, cname;
  Color_t mycolor = kBlack;
  if(charge) fPIDperformance->GetAxis(3)->SetRange(charge, charge);
  // Normalisation by all candidates - no restriction in axis 4 - only for mode == 1 
  if(mode == 1) fPIDperformance->GetAxis(4)->SetRange(2,3);
  hDenom = fPIDperformance->Projection(0);
  hDenom->Sumw2();
  hDenom->SetName("hDenom");
  if(mode == 1) fPIDperformance->GetAxis(4)->SetRange(0, fPIDperformance->GetAxis(4)->GetNbins() + 1);
  // Nominator need a restriction in the 4th axis
  switch(mode){
    case 0: // Electron purity
      fPIDperformance->GetAxis(4)->SetRange(2,3);
      hname = "electronPurity";
      htitle = "Electron Purity";
      break;
    case 1: // Signal purity
      fPIDperformance->GetAxis(4)->SetRange(3,3);   // here signal not divided into charm and beauty
      hname = "signalPurity";
      htitle = "Signal Purity";
      break;
    case 2: // Fake contamination
      fPIDperformance->GetAxis(4)->SetRange(1,1);
      hname = "fakeContamination";
      htitle = "Contamination of misidentified hadrons";
      break;
    default: break;
  }
  switch(charge){
    case 0: 
      cname = "All"; 
      mycolor = kBlue;
      break;
    case 1: 
      cname = "Neg"; 
      mycolor = kBlack;
      break;
    case 2: 
      cname = "Pos"; 
      mycolor = kRed;
      break;
  }
  hname += cname;
  hNom = fPIDperformance->Projection(0);
  hNom->Sumw2();
  hNom->SetName("hNom");
  // Reset axis
  fPIDperformance->GetAxis(4)->SetRange(0, fPIDperformance->GetAxis(4)->GetNbins()+1);
  if(charge) fPIDperformance->GetAxis(3)->SetRange(0, fPIDperformance->GetAxis(3)->GetNbins() + 1);

  // Create Efficiency histogram
  TH1 *hEff = dynamic_cast<TH1D *>(hNom->Clone());
  if(hEff){
    hEff->SetName(hname.Data());
    hEff->SetTitle(htitle.Data());
    hEff->Divide(hDenom);
    hEff->Scale(100.);
    hEff->GetXaxis()->SetTitle("p_{T} / GeV/c");
    hEff->GetYaxis()->SetTitle(mode < 2 ? "Purity / %" : "Contamination / %");
    hEff->GetYaxis()->SetRangeUser(0., 100.);
    hEff->SetStats(kFALSE);
    hEff->SetLineColor(kBlack);
    hEff->SetLineWidth(1);
    hEff->SetMarkerColor(mycolor);
    hEff->SetMarkerStyle(22);
    delete hNom; delete hDenom;
  }
  return hEff;
}

//____________________________________________________________
void AliHFEpostAnalysis::DrawCutEfficiency(Bool_t MC, Int_t source){
  //
  // Calculate efficiency for each cut step
  // Starting from MC steps 
  //
  TCanvas *output = new TCanvas("effCut", "Cut Step efficiency", 800, 600);
  output->cd();
  TLegend *leg = new TLegend(0.6, 0.7, 0.89, 0.89);
  leg->SetHeader("Cut Step:");
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetFillStyle(0);

  AliCFContainer *tracks = fEfficiencyContainer->MakeMergedCFContainer("mergedTracks", "Container for MC and reconstructed Track information", "MCTrackCont:recTrackContReco");
  Int_t nStepMC = fEfficiencyContainer->GetCFContainer("MCTrackCont")->GetNStep();
  AliDebug(1, Form("Number of MC Cut Steps: %d", nStepMC));
  const Int_t markerStart = 24;
  const Int_t colorStart = 1;
  TH1 *hTemp = NULL;

  if(MC){
    if(source > -1 && source < 4){
      AliInfo(Form("Setting source to %d", source))
      for(Int_t istep = 0; istep < tracks->GetNStep(); istep++) tracks->GetAxis(4, istep)->SetRange(source + 1, source + 1);
    }
  }
  AliCFEffGrid effcalc("cutEfficiency", "Cut step efficiency calculation", *tracks);
  Bool_t first = kTRUE;
  // Calculate efficiency for MC Steps
  Int_t histcounter = 0;
  if(MC){
    //
    // Draw plots related to MC values
    //
    effcalc.CalculateEfficiency(AliHFEcuts::kStepMCInAcceptance, AliHFEcuts::kStepMCGenerated);
    hTemp = effcalc.Project(0);
    hTemp->SetName(Form("hEff%d", histcounter));
    hTemp->SetTitle("Cut Step Efficiency");
    hTemp->SetMarkerColor(colorStart + histcounter);
    hTemp->SetMarkerStyle(markerStart + histcounter);
    hTemp->GetXaxis()->SetTitle("p / GeV/c");
    hTemp->GetYaxis()->SetTitle("Efficiency");
    hTemp->GetYaxis()->SetRangeUser(0., 2.);
    hTemp->SetStats(kFALSE);
    hTemp->Draw("ep");
    leg->AddEntry(hTemp, tracks->GetStepTitle(AliHFEcuts::kStepMCInAcceptance), "p");
    histcounter++;
    effcalc.CalculateEfficiency(nStepMC, AliHFEcuts::kStepMCGenerated);
    hTemp = effcalc.Project(0);
    hTemp->SetName("hEff2");
    hTemp->SetTitle("Cut Step Efficiency");
    hTemp->SetMarkerColor(colorStart + histcounter);
    hTemp->SetMarkerStyle(markerStart + histcounter);
    hTemp->GetXaxis()->SetTitle("p / GeV/c");
    hTemp->GetYaxis()->SetTitle("Efficiency");
    hTemp->GetYaxis()->SetRangeUser(0., 2.);
    hTemp->SetStats(kFALSE);
    hTemp->Draw("epsame");
    leg->AddEntry(hTemp, tracks->GetStepTitle(nStepMC), "p");
    first = kFALSE;
    histcounter++;
  }
  for(Int_t istep = nStepMC+1; istep < tracks->GetNStep(); istep++){
    effcalc.CalculateEfficiency(istep, istep - 1); 
    hTemp = effcalc.Project(0);
    hTemp->SetName(Form("hEff%d", istep));
    hTemp->SetTitle("Cut Step Efficiency");
    hTemp->SetMarkerColor(colorStart + histcounter);
    hTemp->SetMarkerStyle(markerStart + histcounter);
    hTemp->SetStats(kFALSE);
    hTemp->Draw(first ? "ep" : "epsame");
    hTemp->GetXaxis()->SetTitle("P / GeV/c");
    hTemp->GetYaxis()->SetTitle("Efficiency");
    hTemp->GetYaxis()->SetRangeUser(0., 2.);
    leg->AddEntry(hTemp, tracks->GetStepTitle(istep), "p");
    first = kFALSE;
    histcounter++;
  }
  leg->Draw();
  delete tracks;
}
