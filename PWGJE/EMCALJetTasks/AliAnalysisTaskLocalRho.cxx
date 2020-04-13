/**************************************************************************************
 * Copyright (C) 2013, Copyright Holders of the ALICE Collaboration                   *
 * All rights reserved.                                                               *
 *                                                                                    *
 * Redistribution and use in source and binary forms, with or without                 *
 * modification, are permitted provided that the following conditions are met:        *
 *     * Redistributions of source code must retain the above copyright               *
 *       notice, this list of conditions and the following disclaimer.                *
 *     * Redistributions in binary form must reproduce the above copyright            *
 *       notice, this list of conditions and the following disclaimer in the          *
 *       documentation and/or other materials provided with the distribution.         *
 *     * Neither the name of the <organization> nor the                               *
 *       names of its contributors may be used to endorse or promote products         *
 *       derived from this software without specific prior written permission.        *
 *                                                                                    *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND    *
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED      *
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE             *
 * DISCLAIMED. IN NO EVENT SHALL ALICE COLLABORATION BE LIABLE FOR ANY                *
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES         *
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;       *
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND        *
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT         *
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS      *
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                       *
 **************************************************************************************/
// root includes
#include <TStyle.h>
#include <TRandom3.h>
#include <TChain.h>
#include <TMath.h>
#include <TF1.h>
#include <TF2.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TProfile.h>
// aliroot includes
#include <AliAnalysisTask.h>
#include <AliAnalysisManager.h>
#include <AliCentrality.h>
#include <AliVVertex.h>
#include <AliESDEvent.h>
#include <AliAODEvent.h>
#include <AliAODTrack.h>
// emcal jet framework includes
#include <AliPicoTrack.h>
#include <AliEmcalJet.h>
#include <AliRhoParameter.h>
#include <AliLocalRhoParameter.h>
#include <AliAnalysisTaskLocalRho.h>
#include "AliMultSelection.h"

class AliAnalysisTaskLocalRho;
using namespace std;

ClassImp(AliAnalysisTaskLocalRho)

//_____________________________________________________________________________
AliAnalysisTaskLocalRho::AliAnalysisTaskLocalRho() : 
  AliAnalysisTaskEmcalJet("AliAnalysisTaskLocalRho", kTRUE), 
  fInitialized(0), fAttachToEvent(kTRUE), fFillHistograms(kFALSE), fNoEventWeightsForQC(kTRUE), 
  fUseScaledRho(0), fCentralityClasses(0), fUserSuppliedV2(0), fUserSuppliedV3(0), fUserSuppliedR2(0), 
  fUserSuppliedR3(0), fNAcceptedTracks(0), fNAcceptedTracksQCn(0), fInCentralitySelection(-1), 
  fFitModulationType(kNoFit), fQCRecovery(kTryFit), fUsePtWeight(kTRUE), fUsePtWeightErrorPropagation(kFALSE), fDetectorType(kTPC), 
  fFitModulationOptions("WLQI"), fRunModeType(kGrid), fFitModulation(0), fMinPvalue(0.01), fMaxPvalue(1), 
  fLocalJetMinEta(-10), fLocalJetMaxEta(-10), fLocalJetMinPhi(-10), fLocalJetMaxPhi(-10), fSoftTrackMinPt(0.15), 
  fSoftTrackMaxPt(5.), fHistPvalueCDF(0), fHistRhoStatusCent(0), fAbsVnHarmonics(kTRUE), fExcludeLeadingJetsFromFit(1.), 
  fRebinSwapHistoOnTheFly(kTRUE), fPercentageOfFits(10.), fUseV0EventPlaneFromHeader(kTRUE), fOutputList(0), 
  fOutputListGood(0), fOutputListBad(0), fHistSwap(0), fHistAnalysisSummary(0), fProfV2(0), fProfV2Cumulant(0), 
  fProfV3(0), fProfV3Cumulant(0) 
{
  // Default constructor

  for(Int_t i(0); i < 10; i++) {
    fHistPsi2[i] = 0; 
    fHistPsi3[i] = 0;
  }
}

//_____________________________________________________________________________
AliAnalysisTaskLocalRho::AliAnalysisTaskLocalRho(const char* name, runModeType type) : 
  AliAnalysisTaskEmcalJet(name, kTRUE),
  fInitialized(0), fAttachToEvent(kTRUE), fFillHistograms(kFALSE), fNoEventWeightsForQC(kTRUE), 
  fUseScaledRho(0), fCentralityClasses(0), fUserSuppliedV2(0), fUserSuppliedV3(0), fUserSuppliedR2(0), 
  fUserSuppliedR3(0), fNAcceptedTracks(0), fNAcceptedTracksQCn(0), fInCentralitySelection(-1), 
  fFitModulationType(kNoFit), fQCRecovery(kTryFit), fUsePtWeight(kTRUE), fUsePtWeightErrorPropagation(kFALSE), fDetectorType(kTPC), 
  fFitModulationOptions("WLQI"), fRunModeType(type), fFitModulation(0), fMinPvalue(0.01), fMaxPvalue(1), 
  fLocalJetMinEta(-10), fLocalJetMaxEta(-10), fLocalJetMinPhi(-10), fLocalJetMaxPhi(-10), fSoftTrackMinPt(0.15), 
  fSoftTrackMaxPt(5.), fHistPvalueCDF(0), fHistRhoStatusCent(0), fAbsVnHarmonics(kTRUE), fExcludeLeadingJetsFromFit(1.), 
  fRebinSwapHistoOnTheFly(kTRUE), fPercentageOfFits(10.), fUseV0EventPlaneFromHeader(kTRUE), fOutputList(0), 
  fOutputListGood(0), fOutputListBad(0), fHistSwap(0), fHistAnalysisSummary(0), fProfV2(0), fProfV2Cumulant(0), 
  fProfV3(0), fProfV3Cumulant(0) 
{
  // Constructor
  for(Int_t i(0); i < 10; i++) {
    fHistPsi2[i] = 0; 
    fHistPsi3[i] = 0;
  }

  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
  switch (fRunModeType) {
  case kLocal : {
    gStyle->SetOptFit(1);
    DefineOutput(2, TList::Class());
    DefineOutput(3, TList::Class());
  } break;
  default: break;   // suppress debug info explicitely when not running locally
  }
}

//_____________________________________________________________________________
AliAnalysisTaskLocalRho::~AliAnalysisTaskLocalRho()
{
  // destructor
  if(fOutputList)             delete fOutputList;
  if(fOutputListGood)         delete fOutputListGood;
  if(fOutputListBad)          delete fOutputListBad;
  if(fFitModulation)          delete fFitModulation;
  if(fHistSwap)               delete fHistSwap;
}

//_____________________________________________________________________________
void AliAnalysisTaskLocalRho::ExecOnce()
{
  #ifdef ALIANALYSISTASKLOCALRHO_DEBUG_FLAG_0
      printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
  #endif
  // Init the analysis
  if(fLocalRhoName=="") fLocalRhoName = Form("LocalRhoFrom_%s", GetName());
  fLocalRho = new AliLocalRhoParameter(fLocalRhoName.Data(), 0); 
  // add the local rho to the event if necessary
  if(fAttachToEvent) {
    if(!(InputEvent()->FindListObject(fLocalRho->GetName()))) {
      InputEvent()->AddObject(fLocalRho);
    } else {
      AliFatal(Form("%s: Container with same name %s already present. Aborting", GetName(), fLocalRho->GetName()));
    }
  }
  AliAnalysisTaskEmcalJet::ExecOnce();        // init the base clas
  if(fUseScaledRho) {
    // unscaled rho has been retrieved by the parent class, now we retrieve rho scaled
    fRho = dynamic_cast<AliRhoParameter*>(InputEvent()->FindListObject(Form("%s_Scaled", fRho->GetName())));
    if(!fRho) {
      AliFatal(Form("%s: Couldn't find container for scaled rho. Aborting !", GetName()));
    }
  }
  if(!GetJetContainer()) AliFatal(Form("%s: Couldn't get jet container. Aborting !", GetName()));
}

//_____________________________________________________________________________
Bool_t AliAnalysisTaskLocalRho::InitializeAnalysis() 
{
  // Initialize the anaysis
  #ifdef ALIANALYSISTASKLOCALRHO_DEBUG_FLAG_0
      printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
  #endif
  if(fLocalJetMinEta > -10 && fLocalJetMaxEta > -10) SetJetEtaLimits(fLocalJetMinEta, fLocalJetMaxEta);
  if(fLocalJetMinPhi > -10 && fLocalJetMaxPhi > -10) SetJetPhiLimits(fLocalJetMinPhi, fLocalJetMaxPhi);
  switch (fFitModulationType)  {
  case kNoFit : { SetModulationFit(new TF1("fit_kNoFit", "[0]", 0, TMath::TwoPi())); } break;
  case kV2 : {
    SetModulationFit(new TF1("fit_kV2", "[0]*([1]+[2]*[3]*TMath::Cos([2]*(x-[4])))", 0, TMath::TwoPi()));
    fFitModulation->SetParameter(0, 0.);        // normalization
    fFitModulation->SetParameter(3, 0.2);       // v2
    fFitModulation->FixParameter(1, 1.);        // constant
    fFitModulation->FixParameter(2, 2.);        // constant
  } break;
  case kV3: {
    SetModulationFit(new TF1("fit_kV3", "[0]*([1]+[2]*[3]*TMath::Cos([2]*(x-[4])))", 0, TMath::TwoPi()));
    fFitModulation->SetParameter(0, 0.);        // normalization
    fFitModulation->SetParameter(3, 0.2);       // v3
    fFitModulation->FixParameter(1, 1.);        // constant
    fFitModulation->FixParameter(2, 3.);        // constant
  } break;
  default : { // for the combined fit, the 'direct fourier series' or the user supplied vn values we use v2 and v3
    SetModulationFit(new TF1("fit_kCombined", "[0]*([1]+[2]*([3]*TMath::Cos([2]*(x-[4]))+[7]*TMath::Cos([5]*(x-[6]))))", 0, TMath::TwoPi()));
    fFitModulation->SetParameter(0, 0.);       // normalization
    fFitModulation->SetParameter(3, 0.2);      // v2
    fFitModulation->FixParameter(1, 1.);       // constant
    fFitModulation->FixParameter(2, 2.);       // constant
    fFitModulation->FixParameter(5, 3.);       // constant
    fFitModulation->SetParameter(7, 0.2);      // v3
  } break;
  }
  switch (fRunModeType) {
  case kGrid : { fFitModulationOptions += "N0"; } break;
  default : break;
  }
  FillAnalysisSummaryHistogram();
  return kTRUE;
}

//_____________________________________________________________________________
void AliAnalysisTaskLocalRho::UserCreateOutputObjects()
{
  // create output objects
  #ifdef ALIANALYSISTASKLOCALRHO_DEBUG_FLAG_0
      printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
  #endif
  fHistSwap = new TH1F("fHistSwap", "fHistSwap", 20, 0, TMath::TwoPi());
  if(!fCentralityClasses) {   // classes must be defined at this point
    Int_t c[] = {0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100};
    fCentralityClasses = new TArrayI(sizeof(c)/sizeof(c[0]), c);
  }
  fOutputList = new TList();
  fOutputList->SetOwner(kTRUE);
  // the analysis summary histo which stores all the analysis flags is always written to file
  fHistAnalysisSummary = BookTH1F("fHistAnalysisSummary", "flag", 51, -0.5, 51.5);
  if(!fFillHistograms) {
    PostData(1, fOutputList);
    return;
  }
  for(Int_t i(0); i < fCentralityClasses->GetSize()-1; i++) {     
    fHistPsi2[i] = BookTH1F("fHistPsi2", "#Psi_{2}", 100, -.5*TMath::Pi(), .5*TMath::Pi(), i);
    fHistPsi3[i] = BookTH1F("fHistPsi3", "#Psi_{3}", 100, -1.*TMath::Pi()/3., TMath::Pi()/3., i);
  }
  // cdf of chisquare distribution
  fHistPvalueCDF = BookTH1F("fHistPvalueCDF", "CDF #chi^{2}", 500, 0, 1);
  fHistRhoStatusCent = BookTH2F("fHistRhoStatusCent", "centrality", "status [0=ok, 1=failed]", 101, -1, 100, 2, -.5, 1.5);
  // vn profiles
  Float_t temp[fCentralityClasses->GetSize()];
  for(Int_t i(0); i < fCentralityClasses->GetSize(); i++) temp[i] = fCentralityClasses->At(i);
  fProfV2 = new TProfile("fProfV2", "fProfV2", fCentralityClasses->GetSize()-1, temp);
  fProfV3 = new TProfile("fProfV3", "fProfV3", fCentralityClasses->GetSize()-1, temp);
  fOutputList->Add(fProfV2);
  fOutputList->Add(fProfV3);
  switch (fFitModulationType) {
  case kQC2 : {
    fProfV2Cumulant = new TProfile("fProfV2Cumulant", "fProfV2Cumulant", fCentralityClasses->GetSize()-1, temp);
    fProfV3Cumulant = new TProfile("fProfV3Cumulant", "fProfV3Cumulant", fCentralityClasses->GetSize()-1, temp);
    fOutputList->Add(fProfV2Cumulant);
    fOutputList->Add(fProfV3Cumulant);
  } break;
  case kQC4 : {
    fProfV2Cumulant = new TProfile("fProfV2Cumulant", "fProfV2Cumulant", fCentralityClasses->GetSize()-1, temp);
    fProfV3Cumulant = new TProfile("fProfV3Cumulant", "fProfV3Cumulant", fCentralityClasses->GetSize()-1, temp);
    fOutputList->Add(fProfV2Cumulant);
    fOutputList->Add(fProfV3Cumulant);
  } break;
  default : break;
  }
  if(fUsePtWeight) fHistSwap->Sumw2();
  if(fUserSuppliedV2) fOutputList->Add(fUserSuppliedV2);
  if(fUserSuppliedV3) fOutputList->Add(fUserSuppliedV3);
  if(fUserSuppliedR2) fOutputList->Add(fUserSuppliedR2);
  if(fUserSuppliedR3) fOutputList->Add(fUserSuppliedR3);
  // increase readability of output list
  fOutputList->Sort();
  PostData(1, fOutputList);
  switch (fRunModeType) {
  case kLocal : {
    fOutputListGood = new TList();
    fOutputListGood->SetOwner(kTRUE);
    fOutputListBad = new TList();
    fOutputListBad->SetOwner(kTRUE);
    PostData(2, fOutputListGood);
    PostData(3, fOutputListBad);
  } break;
  default: break;
  }
}

//_____________________________________________________________________________
TH1F* AliAnalysisTaskLocalRho::BookTH1F(const char* name, const char* x, Int_t bins, Double_t min, Double_t max, Int_t c, Bool_t append)
{
  // Book a TH1F and connect it to the output container
  #ifdef ALIANALYSISTASKLOCALRHO_DEBUG_FLAG_0
      printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
  #endif
  if(!fOutputList) return 0x0;
  TString title(name);
  if(c!=-1) { // format centrality dependent histograms accordingly
    name = Form("%s_%i", name, c);
    title += Form("_%i-%i", fCentralityClasses->At(c), fCentralityClasses->At(1+c));
  }
  title += Form(";%s;[counts]", x);
  TH1F* histogram = new TH1F(name, title.Data(), bins, min, max);
  histogram->Sumw2();
  if(append) fOutputList->Add(histogram);
  return histogram;   
}

//_____________________________________________________________________________
TH2F* AliAnalysisTaskLocalRho::BookTH2F(const char* name, const char* x, const char*y, Int_t binsx, Double_t minx, Double_t maxx, 
					Int_t binsy, Double_t miny, Double_t maxy, Int_t c, Bool_t append)
{
  // Book a TH2F and connect it to the output container
  #ifdef ALIANALYSISTASKLOCALRHO_DEBUG_FLAG_0
      printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
  #endif
  if(!fOutputList) return 0x0;
  TString title(name);
  if(c!=-1) { // format centrality dependent histograms accordingly
    name = Form("%s_%i", name, c);
    title += Form("_%i-%i", fCentralityClasses->At(c), fCentralityClasses->At(1+c));
  }
  title += Form(";%s;%s", x, y);
  TH2F* histogram = new TH2F(name, title.Data(), binsx, minx, maxx, binsy, miny, maxy);
  histogram->Sumw2();
  if(append) fOutputList->Add(histogram);
  return histogram;   
}

//_____________________________________________________________________________
Bool_t AliAnalysisTaskLocalRho::Run()
{
  // Execute once for each event
  #ifdef ALIANALYSISTASKLOCALRHO_DEBUG_FLAG_0
      printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
  #endif
  if(!(InputEvent()||fTracks||fJets||fRho)) return kFALSE;
  if(!fInitialized) fInitialized = InitializeAnalysis();
  // get the centrality bin (necessary for some control histograms
  if(InputEvent()->GetRunNumber()>200000)
  {
      Float_t lPercentile = 300;
      AliMultSelection *MultSelection = 0x0;
      MultSelection = (AliMultSelection*)InputEvent()->FindListObject("MultSelection");
      if(!MultSelection) {
          //If you get this warning (and lPercentiles 300) please check that the AliMultSelectionTask actually ran (before your task)
          AliWarning("AliMultSelection object not found!");
      }else{
          lPercentile = MultSelection->GetMultiplicityPercentile("V0M");
      }
      fCent = lPercentile;
  }
  fInCentralitySelection = -1;

  for(Int_t i(0); i < fCentralityClasses->GetSize()-1; i++) {
    if(fCent >= fCentralityClasses->At(i) && fCent <= fCentralityClasses->At(1+i)) {
      fInCentralitySelection = i;
      break; }
  }
  if(fInCentralitySelection < 0) return kFALSE;
  // set the rho value 
  fLocalRho->SetVal(fRho->GetVal());
  // set the correct event plane accordign to the requested reference detector
  Double_t psi2(-1), psi3(-1);
  switch (fDetectorType) {    // determine the detector type for the rho fit
  case kTPC :     { 
    // [0] psi2         [1] psi3
    Double_t tpc[2];
    CalculateEventPlaneTPC(tpc);
    psi2 = tpc[0];         psi3 = tpc[1]; 
  } break;
  case kVZEROA :  { 
    // [0][0] psi2a     [1,0]   psi2c
    // [0][1] psi3a     [1,1]   psi3c
    Double_t vzero[2][2];
    CalculateEventPlaneVZERO(vzero);
    psi2 = vzero[0][0];    psi3 = vzero[0][1]; 
  }   break;  
  case kVZEROC :  { 
    // [0][0] psi2a     [1,0]   psi2c
    // [0][1] psi3a     [1,1]   psi3c
    Double_t vzero[2][2];
    CalculateEventPlaneVZERO(vzero);
    psi2 = vzero[1][0];    psi3 = vzero[1][1]; 
  }   break;
  case kVZEROComb : { 
    /* for the combined vzero event plane
     * [0] psi2         [1] psi3
     * not fully implmemented yet, use with caution ! */
    Double_t vzeroComb[2];
    CalculateEventPlaneCombinedVZERO(vzeroComb);
    psi2 = vzeroComb[0]; psi3 = vzeroComb[1];
  } break;
  default : break;
  }
  if(fFillHistograms) FillEventPlaneHistograms(psi2, psi3);
  switch (fFitModulationType) { // do the fits
  case kNoFit : { fFitModulation->FixParameter(0, fLocalRho->GetVal()); } break;
  case kV2 : {    // only v2
    if(CorrectRho(psi2, psi3)) {
      if(fFillHistograms) fProfV2->Fill(fCent, fFitModulation->GetParameter(3));
      if(fUserSuppliedR2) {
	Double_t r(fUserSuppliedR2->GetBinContent(fUserSuppliedR2->GetXaxis()->FindBin(fCent)));
	if(r > 0) fFitModulation->SetParameter(3, fFitModulation->GetParameter(3)/r);
      }
    }
  } break;
  case kV3 : {    // only v3
    if(CorrectRho(psi2, psi3)) {
      if(fUserSuppliedR3) {
	Double_t r(fUserSuppliedR3->GetBinContent(fUserSuppliedR3->GetXaxis()->FindBin(fCent)));
	if(r > 0) fFitModulation->SetParameter(3, fFitModulation->GetParameter(3)/r);
      }
      if(fFillHistograms) fProfV3->Fill(fCent, fFitModulation->GetParameter(3));
    }
  } break;
  case kQC2 : {   // qc2 analysis - NOTE: not a wise idea to use this !
    if(CorrectRho(psi2, psi3)) {
      if(fUserSuppliedR2 && fUserSuppliedR3) {
	// note for the qc method, resolution is REVERSED to go back to v2obs
	Double_t r2(fUserSuppliedR2->GetBinContent(fUserSuppliedR2->GetXaxis()->FindBin(fCent)));
	Double_t r3(fUserSuppliedR3->GetBinContent(fUserSuppliedR3->GetXaxis()->FindBin(fCent)));
	if(r2 > 0) fFitModulation->SetParameter(3, fFitModulation->GetParameter(3)*r2);
	if(r3 > 0) fFitModulation->SetParameter(7, fFitModulation->GetParameter(7)*r3);
      }
      if (fUsePtWeight) { // use weighted weights
	Double_t dQCnM11 = (fNoEventWeightsForQC) ? 1. : QCnM11();
	if(fFillHistograms) {
	  fProfV2->Fill(fCent, fFitModulation->GetParameter(3), dQCnM11);
	  fProfV3->Fill(fCent, fFitModulation->GetParameter(7), dQCnM11); 
	}
      } else {
	Double_t dQCnM = (fNoEventWeightsForQC) ? 2. : QCnM();
	if(fFillHistograms) {
	  fProfV2->Fill(fCent, fFitModulation->GetParameter(3), dQCnM*(dQCnM-1));
	  fProfV3->Fill(fCent, fFitModulation->GetParameter(7), dQCnM*(dQCnM-1));
	}
      }
    }
  } break;
  case kQC4 : {   // NOTE: see comment at kQC2
    if(CorrectRho(psi2, psi3)) {
      if(fUserSuppliedR2 && fUserSuppliedR3) {
	// note for the qc method, resolution is REVERSED to go back to v2obs   
	Double_t r2(fUserSuppliedR2->GetBinContent(fUserSuppliedR2->GetXaxis()->FindBin(fCent)));
	Double_t r3(fUserSuppliedR3->GetBinContent(fUserSuppliedR3->GetXaxis()->FindBin(fCent)));
	if(r2 > 0) fFitModulation->SetParameter(3, fFitModulation->GetParameter(3)*r2);
	if(r3 > 0) fFitModulation->SetParameter(7, fFitModulation->GetParameter(7)*r3);
      }
      if (fUsePtWeight) { // use weighted weights
	if(fFillHistograms) {
	  fProfV2->Fill(fCent, TMath::Power(fFitModulation->GetParameter(3),0.5)/*, QCnM1111()*/);
	  fProfV3->Fill(fCent, TMath::Power(fFitModulation->GetParameter(7),0.5)/*, QCnM1111()*/); 
	}
      } else {
	if(fFillHistograms) {
	  fProfV2->Fill(fCent, TMath::Power(fFitModulation->GetParameter(3),0.5)/*, QCnM()*(QCnM()-1)*(QCnM()-2)*(QCnM()-3)*/);
	  fProfV3->Fill(fCent, TMath::Power(fFitModulation->GetParameter(7),0.5)/*, QCnM()*(QCnM()-1)*(QCnM()-2)*(QCnM()-3)*/);
	}
      }
    }
  } break;
  default : {
    if(CorrectRho(psi2, psi3)) {
      if(fUserSuppliedR2 && fUserSuppliedR3) {
	Double_t r2(fUserSuppliedR2->GetBinContent(fUserSuppliedR2->GetXaxis()->FindBin(fCent)));
	Double_t r3(fUserSuppliedR3->GetBinContent(fUserSuppliedR3->GetXaxis()->FindBin(fCent)));
	if(r2 > 0) fFitModulation->SetParameter(3, fFitModulation->GetParameter(3)/r2);
	if(r3 > 0) fFitModulation->SetParameter(7, fFitModulation->GetParameter(7)/r3);
      }
      if(fFillHistograms) {
	fProfV2->Fill(fCent, fFitModulation->GetParameter(3));
	fProfV3->Fill(fCent, fFitModulation->GetParameter(7));
      }
    }
  } break;
  }
  // if all went well, add local rho
  fLocalRho->SetLocalRho(fFitModulation);
  PostData(1, fOutputList);
  return kTRUE;
}

//_____________________________________________________________________________
void AliAnalysisTaskLocalRho::CalculateEventPlaneVZERO(Double_t vzero[2][2]) const 
{
  #ifdef ALIANALYSISTASKLOCALRHO_DEBUG_FLAG_0
      printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
  #endif
  // Get the vzero event plane
  if(fUseV0EventPlaneFromHeader) {
    // use the vzero event plane from the event header
    // note: to use the calibrated vzero event plane, run 
    // $ALICE_ROOT/ANALYSIS/macros/AddTaskVZEROEPSelection.C
    // prior to this task (make sure the calibration is available for the dataset
    // you want to use)
    Double_t a(0), b(0), c(0), d(0), e(0), f(0), g(0), h(0);
    vzero[0][0] = InputEvent()->GetEventplane()->CalculateVZEROEventPlane(InputEvent(), 8, 2, a, b);
    vzero[1][0] = InputEvent()->GetEventplane()->CalculateVZEROEventPlane(InputEvent(), 9, 2, c, d);
    vzero[0][1] = InputEvent()->GetEventplane()->CalculateVZEROEventPlane(InputEvent(), 8, 3, e, f);
    vzero[1][1] = InputEvent()->GetEventplane()->CalculateVZEROEventPlane(InputEvent(), 9, 3, g, h);
    return;
  }
  // grab the vzero event plane without recentering
  if(fDebug > 0) printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
  Double_t qxa2(0), qya2(0), qxc2(0), qyc2(0);    // for psi2
  Double_t qxa3(0), qya3(0), qxc3(0), qyc3(0);    // for psi3
  for(Int_t iVZERO(0); iVZERO < 64; iVZERO++) {
    Double_t phi(TMath::PiOver4()*(.5+iVZERO%8)), /* eta(0), */ weight(InputEvent()->GetVZEROEqMultiplicity(iVZERO));
    //        (iVZERO<32) ? eta = -3.45+.5*(iVZERO/8) : eta = 4.8-.6*((iVZERO/8)-4);
    if(iVZERO<32) {
      qxa2 += weight*TMath::Cos(2.*phi);
      qya2 += weight*TMath::Sin(2.*phi);
      qxa3 += weight*TMath::Cos(3.*phi);
      qya3 += weight*TMath::Sin(3.*phi);
    }
    else {
      qxc2 += weight*TMath::Cos(2.*phi);
      qyc2 += weight*TMath::Sin(2.*phi);
      qxc3 += weight*TMath::Cos(3.*phi);
      qyc3 += weight*TMath::Sin(3.*phi);
    }
  }
  vzero[0][0] = .5*TMath::ATan2(qya2, qxa2);
  vzero[1][0] = .5*TMath::ATan2(qyc2, qxc2);
  vzero[0][1] = (1./3.)*TMath::ATan2(qya3, qxa3);
  vzero[1][1] = (1./3.)*TMath::ATan2(qyc3, qxc3);
}

//_____________________________________________________________________________
void AliAnalysisTaskLocalRho::CalculateEventPlaneTPC(Double_t* tpc)
{
  // Grab the TPC event plane. if parameter fExcludeLeadingJetsFromFit is larger than 0, 
  // strip in eta of width fExcludeLeadingJetsFromFit * GetJetContainer()->GetJetRadius() around the leading jet (before
  // subtraction of rho) will be exluded from the event plane estimate
  #ifdef ALIANALYSISTASKLOCALRHO_DEBUG_FLAG_0
      printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
  #endif
  fNAcceptedTracks = 0;                // reset the track counter
  Double_t qx2(0), qy2(0);     // for psi2
  Double_t qx3(0), qy3(0);     // for psi3
  if(fTracks) {
    Float_t excludeInEta = -999;
    if(fExcludeLeadingJetsFromFit > 0 ) {    // remove the leading jet from ep estimate
      AliEmcalJet* leadingJet(GetJetContainer()->GetLeadingJet());
      if(leadingJet) leadingJet->Eta();
    }
    Int_t iTracks(fTracks->GetEntriesFast());
    for(Int_t iTPC(0); iTPC < iTracks; iTPC++) {
      AliVTrack* track = static_cast<AliVTrack*>(fTracks->At(iTPC));
      if(!PassesCuts(track) || track->Pt() < fSoftTrackMinPt || track->Pt() > fSoftTrackMaxPt) continue;
      if(fExcludeLeadingJetsFromFit > 0 &&( (TMath::Abs(track->Eta() - excludeInEta) < GetJetContainer()->GetJetRadius()*fExcludeLeadingJetsFromFit ) || (TMath::Abs(track->Eta()) - GetJetContainer()->GetJetRadius() - GetJetContainer()->GetJetEtaMax() ) > 0 )) continue;
      fNAcceptedTracks++;
      qx2+= TMath::Cos(2.*track->Phi());
      qy2+= TMath::Sin(2.*track->Phi());
      qx3+= TMath::Cos(3.*track->Phi());
      qy3+= TMath::Sin(3.*track->Phi());
    }
  }
  tpc[0] = .5*TMath::ATan2(qy2, qx2);
  tpc[1] = (1./3.)*TMath::ATan2(qy3, qx3);
} 

//_____________________________________________________________________________
void AliAnalysisTaskLocalRho::CalculateEventPlaneCombinedVZERO(Double_t* comb) const
{
  // Grab the combined vzero event plane
  #ifdef ALIANALYSISTASKLOCALRHO_DEBUG_FLAG_0
      printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
  #endif
  Double_t a(0), b(0), c(0), d(0);
  comb[0] = InputEvent()->GetEventplane()->CalculateVZEROEventPlane(InputEvent(), 10, 2, a, b);
  comb[1] = InputEvent()->GetEventplane()->CalculateVZEROEventPlane(InputEvent(), 10, 3, c, d);
}

//_____________________________________________________________________________
Double_t AliAnalysisTaskLocalRho::CalculateQC2(Int_t harm) 
{
  // Get the second order q-cumulant, a -999 return will be caught in the qa routine of CorrectRho
  #ifdef ALIANALYSISTASKLOCALRHO_DEBUG_FLAG_0
      printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
  #endif
  Double_t reQ(0), imQ(0), modQ(0), M11(0), M(0);
  if(fUsePtWeight) {  // for the weighted 2-nd order q-cumulant
    QCnQnk(harm, 1, reQ, imQ);      // get the weighted 2-nd order q-vectors
    modQ = reQ*reQ+imQ*imQ;         // get abs Q-squared
    M11 = QCnM11();                 // equals S2,1 - S1,2
    return (M11 > 0) ? ((modQ - QCnS(1,2))/M11) : -999;
  } // else return the non-weighted 2-nd order q-cumulant
  QCnQnk(harm, 0, reQ, imQ);          // get the non-weighted 2-nd order q-vectors
  modQ = reQ*reQ+imQ*imQ;             // get abs Q-squared
  M = QCnM();
  return (M > 1) ? (modQ - M)/(M*(M-1)) : -999;
}

//_____________________________________________________________________________
Double_t AliAnalysisTaskLocalRho::CalculateQC4(Int_t harm) 
{
  // Get the fourth order q-cumulant, a -999 return will be caught in the qa routine of CorrectRho
  #ifdef ALIANALYSISTASKLOCALRHO_DEBUG_FLAG_0
      printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
  #endif
  Double_t reQn1(0), imQn1(0), reQ2n2(0), imQ2n2(0), reQn3(0), imQn3(0), M1111(0), M(0);
  Double_t a(0), b(0), c(0), d(0), e(0), f(0), g(0);  // terms of the calculation
  if(fUsePtWeight) {  // for the weighted 4-th order q-cumulant
    QCnQnk(harm, 1, reQn1, imQn1);
    QCnQnk(harm*2, 2, reQ2n2, imQ2n2);
    QCnQnk(harm, 3, reQn3, imQn3);
    // fill in the terms ...
    a = (reQn1*reQn1+imQn1*imQn1)*(reQn1*reQn1+imQn1*imQn1);
    b = reQ2n2*reQ2n2 + imQ2n2*imQ2n2;
    c = -2.*(reQ2n2*reQn1*reQn1-reQ2n2*imQn1*imQn1+2.*imQ2n2*reQn1*imQn1);
    d = 8.*(reQn3*reQn1+imQn3*imQn1);
    e = -4.*QCnS(1,2)*(reQn1*reQn1+imQn1*imQn1);
    f = -6.*QCnS(1,4);
    g = 2.*QCnS(2,2);
    M1111 = QCnM1111();
    return (M1111 > 0) ? (a+b+c+d+e+f+g)/M1111 : -999;
  }   // else return the unweighted case
  Double_t reQn(0), imQn(0), reQ2n(0), imQ2n(0);
  QCnQnk(harm, 0, reQn, imQn);
  QCnQnk(harm*2, 0, reQ2n, imQ2n);
  // fill in the terms ...
  M = QCnM();
  if(M < 4) return -999;
  a = (reQn*reQn+imQn*imQn)*(reQn*reQn+imQn*imQn);
  b = reQ2n*reQ2n + imQ2n*imQ2n;
  c = -2.*(reQ2n*reQn*reQn-reQ2n*imQn*imQn+2.*imQ2n*reQn*imQn);
  e = -4.*(M-2)*(reQn*reQn+imQn*imQn);
  f = 2.*M*(M-3);
  return (a+b+c+e+f)/(M*(M-1)*(M-2)*(M-3));
}

//_____________________________________________________________________________
void AliAnalysisTaskLocalRho::QCnQnk(Int_t n, Int_t k, Double_t &reQ, Double_t &imQ) 
{
  // Get the weighted n-th order q-vector, pass real and imaginary part as reference
  #ifdef ALIANALYSISTASKLOCALRHO_DEBUG_FLAG_0
      printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
  #endif
  if(!fTracks) return;
  fNAcceptedTracksQCn = 0;
  Int_t iTracks(fTracks->GetEntriesFast());
  for(Int_t iTPC(0); iTPC < iTracks; iTPC++) {
    AliVTrack* track = static_cast<AliVTrack*>(fTracks->At(iTPC));
    if(!PassesCuts(track) || track->Pt() < fSoftTrackMinPt || track->Pt() > fSoftTrackMaxPt) continue;
    fNAcceptedTracksQCn++;
    // for the unweighted case, k equals zero and the weight doesn't contribute to the equation below
    reQ += TMath::Power(track->Pt(), k) * TMath::Cos(((double)n)*track->Phi());
    imQ += TMath::Power(track->Pt(), k) * TMath::Sin(((double)n)*track->Phi());
  }
}

//_____________________________________________________________________________
Double_t AliAnalysisTaskLocalRho::QCnS(Int_t i, Int_t j) 
{
  // Get the weighted ij-th order autocorrelation correction
  #ifdef ALIANALYSISTASKLOCALRHO_DEBUG_FLAG_0
      printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
  #endif
  if(!fTracks || i <= 0 || j <= 0) return -999;
  Int_t iTracks(fTracks->GetEntriesFast());
  Double_t Sij(0);
  for(Int_t iTPC(0); iTPC < iTracks; iTPC++) {
    AliVTrack* track = static_cast<AliVTrack*>(fTracks->At(iTPC));
    if(!PassesCuts(track) || track->Pt() < fSoftTrackMinPt || track->Pt() > fSoftTrackMaxPt) continue;
    Sij+=TMath::Power(track->Pt(), j);
  }
  return TMath::Power(Sij, i);
}

//_____________________________________________________________________________
Double_t AliAnalysisTaskLocalRho::QCnM() 
{
  // Get multiplicity for unweighted q-cumulants. function QCnQnk should be called first
  #ifdef ALIANALYSISTASKLOCALRHO_DEBUG_FLAG_0
      printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
  #endif
  return (Double_t) fNAcceptedTracksQCn;
}

//_____________________________________________________________________________
Double_t AliAnalysisTaskLocalRho::QCnM11() 
{
  // Get multiplicity weights for the weighted two particle cumulant
  #ifdef ALIANALYSISTASKLOCALRHO_DEBUG_FLAG_0
      printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
  #endif
  return (QCnS(2,1) - QCnS(1,2));
}

//_____________________________________________________________________________
Double_t AliAnalysisTaskLocalRho::QCnM1111() 
{
  // Get multiplicity weights for the weighted four particle cumulant
  #ifdef ALIANALYSISTASKLOCALRHO_DEBUG_FLAG_0
      printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
  #endif
  return (QCnS(4,1)-6*QCnS(1,2)*QCnS(2,1)+8*QCnS(1,3)*QCnS(1,1)+3*QCnS(2,2)-6*QCnS(1,4));
}

//_____________________________________________________________________________
Bool_t AliAnalysisTaskLocalRho::QCnRecovery(Double_t psi2, Double_t psi3) 
{
  // Decides how to deal with the situation where c2 or c3 is negative 
  // Returns kTRUE depending on whether or not a modulated rho is used for the jet background
  #ifdef ALIANALYSISTASKLOCALRHO_DEBUG_FLAG_0
      printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
  #endif
  if(TMath::AreEqualAbs(fFitModulation->GetParameter(3), .0, 1e-10) && TMath::AreEqualAbs(fFitModulation->GetParameter(7), .0,1e-10)) {
    fFitModulation->SetParameter(7, 0);
    fFitModulation->SetParameter(3, 0);
    fFitModulation->SetParameter(0, fLocalRho->GetVal());
    return kTRUE;   // v2 and v3 have physical null values
  }
  switch (fQCRecovery) {
  case kFixedRho : {      // roll back to the original rho
    fFitModulation->SetParameter(7, 0);
    fFitModulation->SetParameter(3, 0);
    fFitModulation->SetParameter(0, fLocalRho->GetVal());
    return kFALSE;       // rho is forced to be fixed
  }
  case kNegativeVn : {
    Double_t c2(fFitModulation->GetParameter(3));
    Double_t c3(fFitModulation->GetParameter(7));
    if( c2 < 0 ) c2 = -1.*TMath::Sqrt(-1.*c2);
    if( c3 < 0 ) c3 = -1.*TMath::Sqrt(-1.*c3);
    fFitModulation->SetParameter(3, c2);
    fFitModulation->SetParameter(7, c3);
    return kTRUE;        // is this a physical quantity ?
  }
  case kTryFit : {
    fitModulationType tempType(fFitModulationType);  // store temporarily
    fFitModulationType = kCombined;
    fFitModulation->SetParameter(7, 0);
    fFitModulation->SetParameter(3, 0);
    Bool_t pass(CorrectRho(psi2, psi3));         // do the fit and all quality checks
    fFitModulationType = tempType;               // roll back for next event
    return pass;
  }
  default : return kFALSE;
  }
  return kFALSE;
}

//_____________________________________________________________________________
Bool_t AliAnalysisTaskLocalRho::CorrectRho(Double_t psi2, Double_t psi3) 
{
  // Get rho' -> rho(phi)
  // three routines are available, 1 and 2 can be used with or without pt weights
  //  [1] get vn from q-cumulants
  //      in case of cumulants, both cumulants and vn values are stored. in both cases, v2 and v3
  //      are expected. a check is performed to see if rho has no negative local minimum
  //      for full description, see Phys. Rev. C 83, 044913
  //      since the cn distribution has negative values, vn = sqrt(cn) can be imaginary sometimes
  //      in this case one can either roll back to the 'original' fixed rho, do a fit for vn or take use
  //      vn = - sqrt(|cn|) note that because of this, use of q-cumulants is not safe ! 
  //  [2] fitting a fourier expansion to the de/dphi distribution
  //      the fit can be done with either v2, v3 or a combination.
  //      in all cases, a cut can be made on the p-value of the chi-squared value of the fit
  //      and a check can be performed to see if rho has no negative local minimum
  //  [3] get v2 and v3 from user supplied histograms
  //      in this way, a fixed value of v2 and v3 is subtracted w.r.t. whichever event plane is requested
  #ifdef ALIANALYSISTASKLOCALRHO_DEBUG_FLAG_0
      printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
  #endif
  switch (fFitModulationType) {       // for approaches where no fitting is required
  case kQC2 : {
    fFitModulation->FixParameter(4, psi2); 
    fFitModulation->FixParameter(6, psi3);
    fFitModulation->FixParameter(3, CalculateQC2(2));   // set here with cn, vn = sqrt(cn)
    fFitModulation->FixParameter(7, CalculateQC2(3));
    // first fill the histos of the raw cumulant distribution
    if (fUsePtWeight) { // use weighted weights
      Double_t dQCnM11 = (fNoEventWeightsForQC) ? 1. : QCnM11();
      if(fFillHistograms) {
	fProfV2Cumulant->Fill(fCent, fFitModulation->GetParameter(3), dQCnM11);
	fProfV3Cumulant->Fill(fCent, fFitModulation->GetParameter(7), dQCnM11);
      }
    } else {
      Double_t dQCnM = (fNoEventWeightsForQC) ? 2. : QCnM();
      if(fFillHistograms) {
	fProfV2Cumulant->Fill(fCent, fFitModulation->GetParameter(3), dQCnM*(dQCnM-1));
	fProfV3Cumulant->Fill(fCent, fFitModulation->GetParameter(7), dQCnM*(dQCnM-1));
      }
    }
    // then see if one of the cn value is larger than zero and vn is readily available
    if(fFitModulation->GetParameter(3) > 0 && fFitModulation->GetParameter(7) > 0) {
      fFitModulation->FixParameter(3, TMath::Sqrt(fFitModulation->GetParameter(3)));
      fFitModulation->FixParameter(7, TMath::Sqrt(fFitModulation->GetParameter(7)));
    } else if (!QCnRecovery(psi2, psi3)) return kFALSE;  // try to recover the cumulant, this will set v2 and v3
    if(fAbsVnHarmonics && fFitModulation->GetMinimum(0, TMath::TwoPi()) < 0) {  // general check 
      fFitModulation->SetParameter(7, 0);
      fFitModulation->SetParameter(3, 0);
      fFitModulation->SetParameter(0, fLocalRho->GetVal());
      return kFALSE;
    }
    return kTRUE;
  } break;
  case kQC4 : {
    fFitModulation->FixParameter(4, psi2); 
    fFitModulation->FixParameter(6, psi3);
    fFitModulation->FixParameter(3, CalculateQC4(2));   // set here with cn, vn = sqrt(cn)
    fFitModulation->FixParameter(7, CalculateQC4(3));
    // first fill the histos of the raw cumulant distribution
    if (fUsePtWeight) { // use weighted weights
      if(fFillHistograms) {
	fProfV2Cumulant->Fill(fCent, fFitModulation->GetParameter(3)/*, QCnM1111()*/);
	fProfV3Cumulant->Fill(fCent, fFitModulation->GetParameter(7)/*, QCnM1111()*/);
      }
    } else {
      if(fFillHistograms) {
	fProfV2Cumulant->Fill(fCent, fFitModulation->GetParameter(3)/*, QCnM1111()*/);
	fProfV3Cumulant->Fill(fCent, fFitModulation->GetParameter(7)/*, QCnM1111()*/);
      }
    }
    // then see if one of the cn value is larger than zero and vn is readily available
    if(fFitModulation->GetParameter(3) > 0 && fFitModulation->GetParameter(7) > 0) {
      fFitModulation->FixParameter(3, TMath::Sqrt(fFitModulation->GetParameter(3)));
      fFitModulation->FixParameter(7, TMath::Sqrt(fFitModulation->GetParameter(7)));
    } else if (!QCnRecovery(psi2, psi3)) return kFALSE;  // try to recover the cumulant, this will set v2 and v3
    if(fAbsVnHarmonics && fFitModulation->GetMinimum(0, TMath::TwoPi()) < 0) {  // general check 
      fFitModulation->SetParameter(7, 0);
      fFitModulation->SetParameter(3, 0);
      fFitModulation->SetParameter(0, fLocalRho->GetVal());
      return kFALSE;
    }
  } break;
  case kIntegratedFlow : {
    // use v2 and v3 values from an earlier iteration over the data
    fFitModulation->FixParameter(3, fUserSuppliedV2->GetBinContent(fUserSuppliedV2->GetXaxis()->FindBin(fCent)));
    fFitModulation->FixParameter(4, psi2);
    fFitModulation->FixParameter(6, psi3);
    fFitModulation->FixParameter(7, fUserSuppliedV3->GetBinContent(fUserSuppliedV3->GetXaxis()->FindBin(fCent)));
    if(fAbsVnHarmonics && fFitModulation->GetMinimum(0, TMath::TwoPi()) < 0) { 
      fFitModulation->SetParameter(7, 0);
      fFitModulation->SetParameter(3, 0);
      fFitModulation->SetParameter(0, fLocalRho->GetVal());
      return kFALSE;
    }
    return kTRUE;
  }
  default : break;
  }
  TString detector("");
  switch (fDetectorType) {
  case kTPC : detector+="TPC";
    break;
  case kVZEROA : detector+="VZEROA";
    break;
  case kVZEROC : detector+="VZEROC";
    break;
  case kVZEROComb : detector+="VZEROComb";
    break; 
  default: break;
  }
  Int_t iTracks(fTracks->GetEntriesFast());
  Double_t excludeInEta = -999;
  Double_t excludeInPhi = -999;
  Double_t excludeInPt  = -999;
  if(iTracks <= 0 || fLocalRho->GetVal() <= 0 ) return kFALSE;   // no use fitting an empty event ...
  if(fExcludeLeadingJetsFromFit > 0 ) {
    AliEmcalJet* leadingJet = GetJetContainer()->GetLeadingJet();
    if(PassesCuts(leadingJet)) {
      excludeInEta = leadingJet->Eta();
      excludeInPhi = leadingJet->Phi();
      excludeInPt  = leadingJet->Pt();
    }
  }
  fHistSwap->Reset();                 // clear the histogram
  TH1F _tempSwap;     // on stack for quick access
  TH1F _tempSwapN;    // on stack for quick access, bookkeeping histogram
  if(fRebinSwapHistoOnTheFly) {
    if(fNAcceptedTracks < 49) fNAcceptedTracks = 49;       // avoid aliasing effects
    _tempSwap = TH1F("_tempSwap", "_tempSwap", TMath::CeilNint(TMath::Sqrt(fNAcceptedTracks)), 0, TMath::TwoPi());
    if(fUsePtWeightErrorPropagation) _tempSwapN = TH1F("_tempSwapN", "_tempSwapN", TMath::CeilNint(TMath::Sqrt(fNAcceptedTracks)), 0, TMath::TwoPi());
    if(fUsePtWeight) _tempSwap.Sumw2();
  }
  else _tempSwap = *fHistSwap;         // now _tempSwap holds the desired histo
  // non poissonian error when using pt weights
  Double_t totalpts(0.), totalptsquares(0.), totalns(0.);
  for(Int_t i(0); i < iTracks; i++) {
    AliVTrack* track = static_cast<AliVTrack*>(fTracks->At(i));
    if(fExcludeLeadingJetsFromFit > 0 &&( (TMath::Abs(track->Eta() - excludeInEta) < GetJetContainer()->GetJetRadius()*fExcludeLeadingJetsFromFit ) || (TMath::Abs(track->Eta()) - GetJetContainer()->GetJetRadius() - GetJetContainer()->GetJetEtaMax() ) > 0 )) continue;
    if(!PassesCuts(track) || track->Pt() > fSoftTrackMaxPt || track->Pt() < fSoftTrackMinPt) continue;
    if(fUsePtWeight) {
      _tempSwap.Fill(track->Phi(), track->Pt());
      if(fUsePtWeightErrorPropagation) {
        totalpts += track->Pt();
        totalptsquares += track->Pt()*track->Pt();
        totalns += 1;
        _tempSwapN.Fill(track->Phi());
      }
    }
  else _tempSwap.Fill(track->Phi());
  }
  if(fUsePtWeight && fUsePtWeightErrorPropagation) {
    // in the case of pt weights overwrite the poissonian error estimate which is assigned by root by a more sophisticated appraoch
    // the assumption here is that the bin error will be dominated by the uncertainty in the mean pt in a bin and in the uncertainty
    // of the number of tracks in a bin, the first of which will be estimated from the sample standard deviation of all tracks in the 
    // event, for the latter use a poissonian estimate. the two contrubitions are assumed to be uncorrelated
    if(totalns < 1) return kFALSE; // not one track passes the cuts
    for(Int_t l = 0; l < _tempSwap.GetNbinsX(); l++) {
      if(_tempSwapN.GetBinContent(l+1) == 0) {
        _tempSwap.SetBinContent(l+1,0);
        _tempSwap.SetBinError(l+1,0);
      }
      else {
        Double_t vartimesnsq = totalptsquares*totalns - totalpts*totalpts;
        Double_t variance = vartimesnsq/(totalns*(totalns-1.));
        Double_t SDOMSq = variance / _tempSwapN.GetBinContent(l+1);
        Double_t SDOMSqOverMeanSq = SDOMSq * _tempSwapN.GetBinContent(l+1) * _tempSwapN.GetBinContent(l+1) / (_tempSwapN.GetBinContent(l+1) * _tempSwapN.GetBinContent(l+1));
        Double_t poissonfrac = 1./_tempSwapN.GetBinContent(l+1);
        Double_t vartotalfrac = SDOMSqOverMeanSq + poissonfrac;
        Double_t vartotal = vartotalfrac * _tempSwap.GetBinContent(l+1) * _tempSwap.GetBinContent(l+1);
        if(vartotal > 0.0001) _tempSwap.SetBinError(l+1,TMath::Sqrt(vartotal));
        else {
          _tempSwap.SetBinContent(l+1,0);
          _tempSwap.SetBinError(l+1,0);
        }
      }
    }
  }

  fFitModulation->SetParameter(0, fLocalRho->GetVal());
  switch (fFitModulationType) {
  case kNoFit : { fFitModulation->FixParameter(0, fLocalRho->GetVal() ); 
  } break;
  case kV2 : { 
    fFitModulation->FixParameter(4, psi2); 
  } break;
  case kV3 : { 
    fFitModulation->FixParameter(4, psi3); 
  } break;
  case kCombined : {
    fFitModulation->FixParameter(4, psi2); 
    fFitModulation->FixParameter(6, psi3);
  } break;
  case kFourierSeries : {
    // in this approach, an explicit calculation will be made of vn = sqrt(xn^2+yn^2)
    // where x[y] = Integrate[r(phi)cos[sin](n phi)dphi, 0, 2pi]
    Double_t cos2(0), sin2(0), cos3(0), sin3(0), sumPt(0);
    for(Int_t i(0); i < iTracks; i++) {
      AliVTrack* track = static_cast<AliVTrack*>(fTracks->At(i));
      if(!PassesCuts(track) || track->Pt() > fSoftTrackMaxPt || track->Pt() < fSoftTrackMinPt) continue;
      sumPt += track->Pt();
      cos2 += track->Pt()*TMath::Cos(2*PhaseShift(track->Phi()-psi2)); 
      sin2 += track->Pt()*TMath::Sin(2*PhaseShift(track->Phi()-psi2));
      cos3 += track->Pt()*TMath::Cos(3*PhaseShift(track->Phi()-psi3)); 
      sin3 += track->Pt()*TMath::Sin(3*PhaseShift(track->Phi()-psi3));
    }
    fFitModulation->SetParameter(3, TMath::Sqrt(cos2*cos2+sin2*sin2)/fLocalRho->GetVal());
    fFitModulation->SetParameter(4, psi2);
    fFitModulation->SetParameter(6, psi3);
    fFitModulation->SetParameter(7, TMath::Sqrt(cos3*cos3+sin3*sin3)/fLocalRho->GetVal());
  } break;
  default : break;
  }
  _tempSwap.Fit(fFitModulation, fFitModulationOptions.Data(), "", 0, TMath::TwoPi());
  // the quality of the fit is evaluated from 1 - the cdf of the chi square distribution
  Double_t CDF(1.-ChiSquareCDF(fFitModulation->GetNDF(), fFitModulation->GetChisquare()));
  if(fFillHistograms) fHistPvalueCDF->Fill(CDF);
  if(CDF > fMinPvalue && CDF < fMaxPvalue && ( fAbsVnHarmonics && fFitModulation->GetMinimum(0, TMath::TwoPi()) > 0)) { // fit quality
    if(fFillHistograms) fHistRhoStatusCent->Fill(fCent, 0.);
    // for LOCAL didactic purposes, save the  best and the worst fits
    // this routine can produce a lot of output histograms (it's not memory 'safe') and will not work on GRID 
    // since the output will become unmergeable (i.e. different nodes may produce conflicting output)
    switch (fRunModeType) {
    case kLocal : {
      if(gRandom->Uniform(0, 100) > fPercentageOfFits) break;
      static Int_t didacticCounterBest(0);
      TProfile* didacticProfile = (TProfile*)_tempSwap.Clone(Form("Fit_%i_1-CDF_%.3f_cen_%i_%s", didacticCounterBest, CDF, fInCentralitySelection, detector.Data()));
      TF1* didactifFit = (TF1*)fFitModulation->Clone(Form("fit_%i_CDF_%.3f_cen_%i_%s", didacticCounterBest, CDF, fInCentralitySelection, detector.Data()));
      didacticProfile->GetListOfFunctions()->Add(didactifFit);
      fOutputListGood->Add(didacticProfile);
      didacticCounterBest++;
      TH2F* didacticSurface = BookTH2F(Form("surface_%s", didacticProfile->GetName()), "#phi", "#eta", 50, 0, TMath::TwoPi(), 50, -1, 1, -1, kFALSE);
      for(Int_t i(0); i < iTracks; i++) {
	AliVTrack* track = static_cast<AliVTrack*>(fTracks->At(i));
	if(PassesCuts(track)) {
	  if(fUsePtWeight) didacticSurface->Fill(track->Phi(), track->Eta(), track->Pt());
	  else didacticSurface->Fill(track->Phi(), track->Eta());
	}
      }
      if(fExcludeLeadingJetsFromFit) {       // visualize the excluded region
	TF2 *f2 = new TF2(Form("%s_LJ", didacticSurface->GetName()),"[0]*TMath::Gaus(x,[1],[2])*TMath::Gaus(y,[3],[4])", 0, TMath::TwoPi(), -1, 1);
	f2->SetParameters(excludeInPt/3.,excludeInPhi,.1,excludeInEta,.1);
	didacticSurface->GetListOfFunctions()->Add(f2);
      }
      fOutputListGood->Add(didacticSurface);
    } break;
    default : break;
    }
  } else {    // if the fit is of poor quality revert to the original rho estimate
    if(fFillHistograms) fHistRhoStatusCent->Fill(fCent, 1.);
    switch (fRunModeType) { // again see if we want to save the fit
    case kLocal : {
      static Int_t didacticCounterWorst(0);
      if(gRandom->Uniform(0, 100) > fPercentageOfFits) break;
      TProfile* didacticProfile = (TProfile*)_tempSwap.Clone(Form("Fit_%i_1-CDF_%.3f_cen_%i_%s", didacticCounterWorst, CDF, fInCentralitySelection, detector.Data() ));
      TF1* didactifFit = (TF1*)fFitModulation->Clone(Form("fit_%i_p_%.3f_cen_%i_%s", didacticCounterWorst, CDF, fInCentralitySelection, detector.Data()));
      didacticProfile->GetListOfFunctions()->Add(didactifFit);
      fOutputListBad->Add(didacticProfile);
      didacticCounterWorst++;
    } break;
    default : break;
    }
    switch (fFitModulationType) {
    case kNoFit : break;        // nothing to do
    case kCombined : fFitModulation->SetParameter(7, 0);        // no break
    case kFourierSeries : fFitModulation->SetParameter(7, 0);   // no break
    default : { // needs to be done if there was a poor fit
      fFitModulation->SetParameter(3, 0);
      fFitModulation->SetParameter(0, fLocalRho->GetVal());
    } break;
    }
    return kFALSE;  // return false if the fit is rejected
  }
  return kTRUE;
}

//_____________________________________________________________________________
void AliAnalysisTaskLocalRho::FillAnalysisSummaryHistogram() const
{
  // Fill the analysis summary histrogram, saves all relevant analysis settigns
  #ifdef ALIANALYSISTASKLOCALRHO_DEBUG_FLAG_0
      printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
  #endif
  fHistAnalysisSummary->GetXaxis()->SetBinLabel(2, "fJetRadius");
  fHistAnalysisSummary->SetBinContent(2, GetJetContainer()->GetJetRadius());
  fHistAnalysisSummary->GetXaxis()->SetBinLabel(3, "fJetEtaMin");
  fHistAnalysisSummary->SetBinContent(3, GetJetContainer()->GetJetEtaMin());
  fHistAnalysisSummary->GetXaxis()->SetBinLabel(4, "fJetEtaMax");
  fHistAnalysisSummary->SetBinContent(4, GetJetContainer()->GetJetEtaMax());
  fHistAnalysisSummary->GetXaxis()->SetBinLabel(5, "fJetPhiMin");
  fHistAnalysisSummary->SetBinContent(5, GetJetContainer()->GetJetPhiMin());
  fHistAnalysisSummary->GetXaxis()->SetBinLabel(6, "fJetPhiMax");
  fHistAnalysisSummary->SetBinContent(6, GetJetContainer()->GetJetPhiMin());
  fHistAnalysisSummary->GetXaxis()->SetBinLabel(34, "fitModulationType");
  fHistAnalysisSummary->SetBinContent(34, (int)fFitModulationType);
  fHistAnalysisSummary->GetXaxis()->SetBinLabel(35, "runModeType");
  fHistAnalysisSummary->SetBinContent(35, (int)fRunModeType);
  fHistAnalysisSummary->GetXaxis()->SetBinLabel(37, "iterator");
  fHistAnalysisSummary->SetBinContent(37, 1.);
  fHistAnalysisSummary->GetXaxis()->SetBinLabel(38, "fMinPvalue");
  fHistAnalysisSummary->SetBinContent(38, fMinPvalue);
  fHistAnalysisSummary->GetXaxis()->SetBinLabel(39, "fMaxPvalue");
  fHistAnalysisSummary->SetBinContent(39, fMaxPvalue);
  fHistAnalysisSummary->GetXaxis()->SetBinLabel(40, "fExcludeLeadingJetsFromFit");
  fHistAnalysisSummary->SetBinContent(40, fExcludeLeadingJetsFromFit);
  fHistAnalysisSummary->GetXaxis()->SetBinLabel(41, "fRebinSwapHistoOnTheFly");
  fHistAnalysisSummary->SetBinContent(41, (int)fRebinSwapHistoOnTheFly);
  fHistAnalysisSummary->GetXaxis()->SetBinLabel(42, "fUsePtWeight");
  fHistAnalysisSummary->SetBinContent(42, (int)fUsePtWeight);
  fHistAnalysisSummary->GetXaxis()->SetBinLabel(45, "fLocalJetMinEta");
  fHistAnalysisSummary->SetBinContent(45,fLocalJetMinEta );
  fHistAnalysisSummary->GetXaxis()->SetBinLabel(46, "fLocalJetMaxEta");
  fHistAnalysisSummary->SetBinContent(46, fLocalJetMaxEta);
  fHistAnalysisSummary->GetXaxis()->SetBinLabel(47, "fLocalJetMinPhi");
  fHistAnalysisSummary->SetBinContent(47, fLocalJetMinPhi);
  fHistAnalysisSummary->GetXaxis()->SetBinLabel(48, "fLocalJetMaxPhi");
  fHistAnalysisSummary->SetBinContent(48, fLocalJetMaxPhi);
  fHistAnalysisSummary->GetXaxis()->SetBinLabel(49, "fSoftTrackMinPt");
  fHistAnalysisSummary->SetBinContent(49, fSoftTrackMinPt);
  fHistAnalysisSummary->GetXaxis()->SetBinLabel(50, "fSoftTrackMaxPt");
  fHistAnalysisSummary->SetBinContent(50, fSoftTrackMaxPt);
  fHistAnalysisSummary->GetXaxis()->SetBinLabel(51, "fUseScaledRho");
  fHistAnalysisSummary->SetBinContent(51, fUseScaledRho);
}

//_____________________________________________________________________________
void AliAnalysisTaskLocalRho::FillEventPlaneHistograms(Double_t psi2, Double_t psi3) const
{
  // Fill event plane histograms
  #ifdef ALIANALYSISTASKLOCALRHO_DEBUG_FLAG_0
      printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
  #endif
  fHistPsi2[fInCentralitySelection]->Fill(psi2);
  fHistPsi3[fInCentralitySelection]->Fill(psi3);    
}

//_____________________________________________________________________________
void AliAnalysisTaskLocalRho::Terminate(Option_t *)
{
   #ifdef ALIANALYSISTASKLOCALRHO_DEBUG_FLAG_0
      printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
   #endif
 // Terminate
}

//_____________________________________________________________________________
void AliAnalysisTaskLocalRho::SetModulationFit(TF1* fit) 
{
  // Set function to fit modulation
  #ifdef ALIANALYSISTASKLOCALRHO_DEBUG_FLAG_0
      printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
  #endif
  if (fFitModulation) delete fFitModulation;
  fFitModulation = fit; 
}
