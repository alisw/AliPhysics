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

// AliAnalysisTaskCascadeTester:
// author: Redmer Alexander Bertens (rbertens@cern.ch)

#include "TChain.h"
#include "TH1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TObjArray.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliCentrality.h"
#include "AliAnalysisTaskCascadeTester.h"
#include "AliPIDCombined.h"
#include "AliMCEvent.h"
#include "TProfile.h"
#include "AliFlowCandidateTrack.h"
#include "AliFlowTrackCuts.h"
#include "AliFlowEventSimple.h"
#include "AliFlowTrackSimple.h"
#include "AliFlowCommonConstants.h"
#include "AliFlowEvent.h"
#include "TVector3.h"
#include "AliAODVZERO.h"
#include "AliPIDResponse.h"
#include "AliTPCPIDResponse.h"
#include "AliAODMCParticle.h"
#include "AliAnalysisTaskVnV0.h"
#include "AliEventPoolManager.h"
#include "AliMultSelection.h"
#include "TMatrixDSym.h"
#include "AliVVertex.h"
#include "AliAODcascade.h"

class AliFlowTrackCuts;

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskCascadeTester)

AliAnalysisTaskCascadeTester::AliAnalysisTaskCascadeTester() : AliAnalysisTaskSE(),
   fDebug(0), fIsMC(0), fEventMixing(0), fTypeMixing(0), fQA(0), fV0(0), fMassBins(1), fMinMass(-1.), fMaxMass(0.), fCutsRP(NULL), fNullCuts(0), fPIDResponse(0), fFlowEvent(0), fBayesianResponse(0), fCandidates(0), fCandidateEtaPtCut(0), fCandidateMinEta(0), fCandidateMaxEta(0), fCandidateMinPt(0), fCandidateMaxPt(0), fCandidateYCut(kFALSE), fCandidateMinY(-.5), fCandidateMaxY(.5), fNPtBins(18), fCentrality(999), fVertex(999), fAOD(0), fPoolManager(0), fOutputList(0), fEventStats(0), fCentralityPass(0), fCentralityNoPass(0), fNOPID(0), fPIDk(0),fNOPIDTOF(0), fPIDTOF(0), fPtP(0), fPtN(0), fPtKP(0), fPtKN(0), fMultCorAfterCuts(0), fMultvsCentr(0), fCentralityMin(0), fCentralityMax(100), fkCentralityMethodA(0), fkCentralityMethodB(0), fCentralityCut2010(0), fCentralityCut2011(0), fPOICuts(0), fVertexRange(0), fPhi(0), fPt(0), fEta(0), fVZEROA(0), fVZEROC(0), fTPCM(0)/*, fDeltaDipAngle(0), fDeltaDipPt(0), fApplyDeltaDipCut(0)*/, fDCAAll(0), fDCAXYQA(0), fDCAZQA(0), fDCAPrim(0), fDCASecondaryWeak(0), fDCAMaterial(0), fSubEventDPhiv2(0), fSkipEventSelection(0), fUsePidResponse(0), fPIDCombined(0), fPileUp(kTRUE), fOmega(kTRUE), fXi(kFALSE)
{
   // Default constructor
   for(Int_t i(0); i < 7; i++) fPIDConfig[i] = 1000.;
   for(Int_t i(0); i < 5; i++) fDCAConfig[i] = 0.;
   for(Int_t i(0); i < 20; i++) {
       fVertexMixingBins[i] = 0;
       fCentralityMixingBins[i] = 0;
   }
   fMixingParameters[0] = 1000; fMixingParameters[1] = 50000; fMixingParameters[2] = 5;
   for(Int_t i(0); i < 18; i++) {
       for(Int_t j(0); j < 2; j++) fV0Data[i][j] = 0;
       fInvMNP[i] = 0; fInvMNN[i] = 0; fInvMPP[i] = 0; fPtSpectra[i] = 0; fPtBins[i] = 0.;
   }
}
//_____________________________________________________________________________
AliAnalysisTaskCascadeTester::AliAnalysisTaskCascadeTester(const char *name) : AliAnalysisTaskSE(name),
   fDebug(0), fIsMC(0), fEventMixing(0), fTypeMixing(0), fQA(0), fV0(0), fMassBins(1), fMinMass(-1.), fMaxMass(0.), fCutsRP(NULL), fNullCuts(0), fPIDResponse(0), fFlowEvent(0), fBayesianResponse(0), fCandidates(0), fCandidateEtaPtCut(0), fCandidateMinEta(0), fCandidateMaxEta(0), fCandidateMinPt(0), fCandidateMaxPt(0), fCandidateYCut(kFALSE), fCandidateMinY(-.5), fCandidateMaxY(.5), fNPtBins(18), fCentrality(999), fVertex(999), fAOD(0), fPoolManager(0), fOutputList(0), fEventStats(0), fCentralityPass(0), fCentralityNoPass(0), fNOPID(0), fPIDk(0), fNOPIDTOF(0), fPIDTOF(0), fPtP(0), fPtN(0), fPtKP(0), fPtKN(0), fMultCorAfterCuts(0), fMultvsCentr(0), fCentralityMin(0), fCentralityMax(100), fkCentralityMethodA(0), fkCentralityMethodB(0), fCentralityCut2010(0), fCentralityCut2011(0), fPOICuts(0), fVertexRange(0), fPhi(0), fPt(0), fEta(0), fVZEROA(0), fVZEROC(0), fTPCM(0)/*, fDeltaDipAngle(0), fDeltaDipPt(0), fApplyDeltaDipCut(0)*/, fDCAAll(0), fDCAXYQA(0), fDCAZQA(0), fDCAPrim(0), fDCASecondaryWeak(0), fDCAMaterial(0), fSubEventDPhiv2(0), fSkipEventSelection(0), fUsePidResponse(0), fPIDCombined(0), fPileUp(kTRUE), fOmega(kTRUE), fXi(kFALSE)
{
   // Constructor
   for(Int_t i(0); i < 7; i++) fPIDConfig[i] = 1000.;
   for(Int_t i(0); i < 5; i++) fDCAConfig[i] = 0.;
   for(Int_t i(0); i < 20; i++) {
       fVertexMixingBins[i] = 0;
       fCentralityMixingBins[i] = 0;
   }
   fMixingParameters[0] = 1000; fMixingParameters[1] = 50000; fMixingParameters[2] = 5;
   for(Int_t i(0); i < 18; i++) {
       for(Int_t j(0); j < 2; j++) fV0Data[i][j] = 0;
       fInvMNP[i] = 0; fInvMNN[i] = 0; fInvMPP[i] = 0; fPtSpectra[i] = 0; fPtBins[i] = 0.;
   }
   DefineInput(0, TChain::Class());
   DefineOutput(1, TList::Class());
   DefineOutput(2, AliFlowEventSimple::Class());
   if(fDebug > 0) cout << " === Created instance of AliAnaysisTaskCascadeTester === " << endl;
}
//_____________________________________________________________________________
AliAnalysisTaskCascadeTester::~AliAnalysisTaskCascadeTester()
{
   // Destructor
   if (fNullCuts) delete fNullCuts;
   if (fOutputList) delete fOutputList;
   if (fCandidates) delete fCandidates;
   if (fFlowEvent) delete fFlowEvent;
   if (fEventMixing) delete fPoolManager;
   if (fPIDCombined) delete fPIDCombined;
   if (fDebug > 0) cout << " === Deleted instance of AliAnalysisTaskCascadeTester === " << endl;
}
//_____________________________________________________________________________
TH1F* AliAnalysisTaskCascadeTester::BookHistogram(const char* name)
{
   // Return a pointer to a TH1 with predefined binning
   if(fDebug > 0) cout << " *** BookHistogram() *** " << name << endl;
   TH1F *hist = 0x0;
   if(fOmega) hist = new TH1F(name, Form("M_{INV} (%s)", name), 60, 1.64, 1.71);
   else if(fXi) hist = new TH1F(name, Form("M_{INV} (%s)", name), 60, 1.29, 1.36);
   hist->GetXaxis()->SetTitle("M_{INV} (GeV / c^{2})");
   hist->GetYaxis()->SetTitle("No. of pairs");
   hist->SetMarkerStyle(kFullCircle);
   hist->Sumw2();
   fOutputList->Add(hist);
   return hist;
}
//_____________________________________________________________________________
TH2F* AliAnalysisTaskCascadeTester::BookPIDHistogram(const char* name, Bool_t TPC)
{
   // Return a pointer to a TH2 with predefined binning
   if(fDebug > 0) cout << " *** BookPIDHisotgram() *** " << endl;
   TH2F *hist = 0x0;
   if(TPC) {
       hist = new TH2F(name, Form("PID (%s)", name), 100, 0, 5, 100, 0, 1000);
       hist->GetYaxis()->SetTitle("dE/dx (a.u.)");
   }
   if(!TPC) {
       hist = new TH2F(name, Form("PID (%s)", name), 100, 0, 5, 100, 0, 1.5);
       hist->GetYaxis()->SetTitle("#beta");
   }
   hist->GetXaxis()->SetTitle("P (GeV / c)");
   fOutputList->Add(hist);
   return hist;
}
//_____________________________________________________________________________
TH1F* AliAnalysisTaskCascadeTester::InitPtSpectraHistograms(Float_t nmin, Float_t nmax)
{
   // intialize p_t histograms for each p_t bin
   if(fDebug > 0) cout << " *** InitPtSpectraHistograms() *** " << endl;
   TH1F* hist = new TH1F(Form("%4.2f p_{t} %4.2f", nmin, nmax), Form("%f p_{t} %f", nmin, nmax), 60, nmin, nmax);
   hist->GetXaxis()->SetTitle("p_{T} GeV / c");
   fOutputList->Add(hist);
   return hist;
}
//_____________________________________________________________________________
TH1F* AliAnalysisTaskCascadeTester::BookPtHistogram(const char* name)
{
   // Return a pointer to a p_T spectrum histogram
   if(fDebug > 0) cout << " *** BookPtHistogram() *** " << endl;
   TH1F* ratio = new TH1F(name, name, 100, 0, 7);
   ratio->GetXaxis()->SetTitle("p_{T} ( GeV / c^{2} )");
   ratio->GetYaxis()->SetTitle("No. of events");
   ratio->Sumw2();
   fOutputList->Add(ratio);
   return ratio;
}
//_____________________________________________________________________________
void AliAnalysisTaskCascadeTester::AddPhiIdentificationOutputObjects()
{
   // Add Phi Identification Output Objects
   if(fDebug > 0) cout << " ** AddPhiIdentificationOutputObjects() *** " << endl;
   if(fQA) {
       fEventStats = new TH1F("fHistStats", "Event Statistics", 18, -.25, 4.25);
       fEventStats->GetXaxis()->SetTitle("No. of events");
       fEventStats->GetYaxis()->SetTitle("Statistics");
       fOutputList->Add(fEventStats);
   }
   fCentralityPass = new TH1F("fCentralityPass", "Centrality Pass", 101, -1, 100);
   fOutputList->Add(fCentralityPass);
   if(fQA) {
       fCentralityNoPass = new TH1F("fCentralityNoPass", "Centrality No Pass", 101, -1, 100);
       fOutputList->Add(fCentralityNoPass);
       fNOPID = BookPIDHistogram("TPC signal, all particles", kTRUE);
       fPIDk = BookPIDHistogram("TPC signal, kaons", kTRUE);
       fNOPIDTOF = BookPIDHistogram("TOF signal, all particles", kFALSE);
       fPIDTOF = BookPIDHistogram("TOF signal, kaons", kFALSE);
   }
   for(Int_t ptbin(0); ptbin < fNPtBins; ptbin++) {
       fInvMNP[ptbin] = BookHistogram(Form("NP, %4.2f < p_{T} < %4.2f GeV", fPtBins[ptbin], fPtBins[ptbin+1]));
       fInvMNN[ptbin] = BookHistogram(Form("NN, %4.2f < p_{T} < %4.2f GeV", fPtBins[ptbin], fPtBins[ptbin+1]));
       fInvMPP[ptbin] = BookHistogram(Form("PP, %4.2f < p_{T} < %4.2f GeV", fPtBins[ptbin], fPtBins[ptbin+1]));
   }
   if(fQA) {
       for(Int_t ptbin(0); ptbin < fNPtBins; ptbin++) fPtSpectra[ptbin] = InitPtSpectraHistograms(fPtBins[ptbin], fPtBins[ptbin+1]);
       fPtP = BookPtHistogram("i^{+}");
       fPtN = BookPtHistogram("i^{-}");
       fPtKP = BookPtHistogram("K^{+}");
       fPtKN = BookPtHistogram("K^{-}");
       fPhi = new TH1F("fPhi", "#phi distribution", 100, -.5, 7);
       fOutputList->Add(fPhi);
       fPt = new TH1F("fPt", "p_{T}", 100, 0, 5.5);
       fOutputList->Add(fPt);
       fEta = new TH1F("fEta", "#eta distribution", 100, -1.1, 1.1);
       fOutputList->Add(fEta);
       fVZEROA = new TH1F("fVZEROA", "VZERO A Multiplicity", 1000, 0, 10000);
       fOutputList->Add(fVZEROA);
       fVZEROC = new TH1F("fVZEROC", "VZERO C Multiplicity", 1000, 0, 10000);
       fOutputList->Add(fVZEROC);
       fTPCM = new TH1F("fTPCM", "TPC multiplicity", 1000, 0, 10000);
       fOutputList->Add(fTPCM);
       fDCAXYQA = new TH1F("fDCAXYQA", "fDCAXYQA", 1000, -5, 5);
       fOutputList->Add(fDCAXYQA);
       fDCAZQA = new TH1F("fDCAZQA", "fDCAZQA", 1000, -5, 5);
       fOutputList->Add(fDCAZQA);
       if(fCentralityCut2010 || fCentralityCut2011) {
           fMultCorAfterCuts = new TH2F("fMultCorAfterCuts", "TPC vs Global multiplicity (After cuts); Global multiplicity; TPC multiplicity", 100, 0, 3000, 100, 0, 3000);
           fOutputList->Add(fMultCorAfterCuts);
           fMultvsCentr = new TH2F("fMultvsCentr", "Multiplicity vs centrality; centrality; Multiplicity", 9, -0.5, 100.5, 101, 0, 3000);
           fOutputList->Add(fMultvsCentr);
       }
   }
   if(fIsMC || fQA) {
       fDCAAll = new TH2F("fDCAAll", "fDCAAll", 1000, 0, 10, 1000, -5, 5);
       fOutputList->Add(fDCAAll);
       fDCAPrim = new TH2F("fDCAprim","fDCAprim", 1000, 0, 10, 1000, -5, 5);
       fOutputList->Add(fDCAPrim);
       fDCASecondaryWeak = new TH2F("fDCASecondaryWeak","fDCASecondaryWeak", 1000, 0, 10, 1000, -5, 5);
       fOutputList->Add(fDCASecondaryWeak);
       fDCAMaterial = new TH2F("fDCAMaterial","fDCAMaterial", 1000, 0, 10, 1000, -5, 5);
       fOutputList->Add(fDCAMaterial);
   }
   if(fV0) {
       fSubEventDPhiv2 = new TProfile("fSubEventDPhiv2", "fSubEventDPhiv2", 5, 0, 5);
       fSubEventDPhiv2->GetXaxis()->SetBinLabel(1, "<#Psi_{a} - #Psi_{b}>");
       fSubEventDPhiv2->GetXaxis()->SetBinLabel(2, "<#Psi_{a} - #Psi_{c}>");
       fSubEventDPhiv2->GetXaxis()->SetBinLabel(3, "<#Psi_{b} - #Psi_{c}>");
       fSubEventDPhiv2->GetXaxis()->SetBinLabel(4, "#sqrt{#frac{<#Psi_{a} - #Psi_{b}><#Psi_{a} - #Psi_{c}>}{<#Psi_{b} - #Psi_{c}>}}");
       fSubEventDPhiv2->GetXaxis()->SetBinLabel(5, "#sqrt{#frac{<#Psi_{a} - #Psi_{c}><#Psi_{b} - #Psi_{c}>}{<#Psi_{a} - #Psi_{b}>}}");
       fOutputList->Add(fSubEventDPhiv2);
       const char* V0[] = {"V0A", "V0C"};
       for(Int_t ptbin(0); ptbin < fNPtBins; ptbin++)
           for(Int_t iV0(0); iV0 < 2; iV0++) {
                   fV0Data[ptbin][iV0] = new TProfile(Form("%s v2 %4.2f < p_{T} < %4.2f GeV", V0[iV0], fPtBins[ptbin], fPtBins[ptbin+1]), Form("%s v2 %4.2f < p_{T} < %4.2f GeV", V0[iV0], fPtBins[ptbin], fPtBins[ptbin+1]), fMassBins, fMinMass, fMaxMass);
                   fOutputList->Add(fV0Data[ptbin][iV0]);
           }
   }
}
//_____________________________________________________________________________
void AliAnalysisTaskCascadeTester::UserCreateOutputObjects()
{
   // Create user defined output objects
   if(fDebug > 0) cout << " *** UserCreateOutputObjects() *** " << endl;
   fNullCuts = new AliFlowTrackCuts("null_cuts");
      // combined pid
   fPIDCombined = new AliPIDCombined;
   fPIDCombined->SetDefaultTPCPriors();
   fPIDCombined->SetDetectorMask(AliPIDResponse::kDetTPC|AliPIDResponse::kDetTOF);

   Double_t t(0);
   for(Int_t i = 0; i < 7; i++) t+=TMath::Abs(fPIDConfig[i]);
   if(t < 0.1) AliFatal("No valid PID procedure recognized -- terminating analysis !!!");
   if(fNPtBins > 18) AliFatal("Invalid number of pt bins initialied ( > 18 ) -- terminating analysis !!!");
   AliFlowCommonConstants* cc = AliFlowCommonConstants::GetMaster();
   cc->SetNbinsQ(500);           cc->SetNbinsPhi(180);           cc->SetNbinsMult(10000);
   cc->SetQMin(0.0);             cc->SetPhiMin(0.0);             cc->SetMultMin(0);
   cc->SetQMax(3.0);             cc->SetPhiMax(TMath::TwoPi());  cc->SetMultMax(10000);
   cc->SetNbinsMass(fMassBins);  cc->SetNbinsEta(200);           (fMassBins == 1) ? cc->SetNbinsPt(15) : cc->SetNbinsPt(100); // high pt
   cc->SetMassMin(fMinMass);     cc->SetEtaMin(-5.0);            cc->SetPtMin(0);
   cc->SetMassMax(fMaxMass);     cc->SetEtaMax(+5.0);            (fMassBins == 1) ? cc->SetPtMax(15) : cc->SetPtMax(10); // high pt
   AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
   if (man) {
      AliInputEventHandler* inputHandler = (AliInputEventHandler*)(man->GetInputEventHandler());
      if (inputHandler)   fPIDResponse = inputHandler->GetPIDResponse();
   }
   // Create all output objects and store them to a list
   fOutputList = new TList();
   fOutputList->SetOwner(kTRUE);
   // Create and post the Phi identification output objects
   AddPhiIdentificationOutputObjects();
   PostData(1, fOutputList);
   // create candidate array
   fCandidates = new TObjArray(1000);
   fCandidates->SetOwner(kTRUE);
   // create and post flowevent
   fFlowEvent = new AliFlowEvent(10000);
   PostData(2, fFlowEvent);
   if(fEventMixing) fPoolManager = InitializeEventMixing();
}
//_____________________________________________________________________________
AliEventPoolManager* AliAnalysisTaskCascadeTester::InitializeEventMixing()
{
   // initialize event mixing
  if(fDebug > 0) cout << " *** InitializeEventMixing() *** " << endl;
  Int_t _c(0), _v(0);
  for(Int_t i(0); i < 19; i++) {
      if (fCentralityMixingBins[i+1] < fCentralityMixingBins[i]) { _c = i; break; }
      else _c = 19;
  }
  for(Int_t i(0); i < 19; i++) {
      if (fVertexMixingBins[i+1] < fVertexMixingBins[i]) { _v = i; break; }
      else _v = 19;
  }
  if(fDebug > 0 ) cout << Form("   --> found %d centrality bins for mixing, %d vertex bins for mixing", _c, _v) << endl;
  Double_t centralityBins[_c];
  Double_t vertexBins[_v];
  for(Int_t i(0); i < _c + 1; i++) centralityBins[i] = fCentralityMixingBins[i];
  for(Int_t i(0); i < _v + 1; i++) vertexBins[i] = fVertexMixingBins[i];
  return new AliEventPoolManager(fMixingParameters[0], fMixingParameters[1], _c, (Double_t*)centralityBins, _v, (Double_t*)vertexBins);
}
//_____________________________________________________________________________
template <typename T> Double_t AliAnalysisTaskCascadeTester::InvariantMass(const T* track1, const T* track2) const
{
   // Return the invariant mass of two tracks, assuming both tracks are kaons
   if(fDebug > 1) cout << " *** InvariantMass() *** " << endl;
   if ((!track2) || (!track1)) return 0.;
   Double_t masss = TMath::Power(4.93676999999999977e-01, 2);
   Double_t pxs = TMath::Power((track1->Px() + track2->Px()), 2);
   Double_t pys = TMath::Power((track1->Py() + track2->Py()), 2);
   Double_t pzs = TMath::Power((track1->Pz() + track2->Pz()), 2);
   Double_t e1 = TMath::Sqrt(track1->P() * track1->P() + masss);
   Double_t e2 = TMath::Sqrt(track2->P() * track2->P() + masss);
   Double_t es = TMath::Power((e1 + e2), 2);
   if ((es - (pxs + pys + pzs)) < 0) return 0.;
   return TMath::Sqrt((es - (pxs + pys + pzs)));
}
//_____________________________________________________________________________
/*
template <typename T> Double_t AliAnalysisTaskCascadeTester::DeltaDipAngle(const T* track1, const T* track2) const
{
   // Calculate the delta dip angle between two particles
   if(fDebug > 1) cout << " *** DeltaDipAngle() *** " << endl;
   if (track1->P()*track2->P() == 0) return 999;
   return TMath::ACos(((track1->Pt() * track2->Pt()) + (track1->Pz() * track2->Pz())) / (track1->P() * track2->P()));
}
//_____________________________________________________________________________
template <typename T> Bool_t AliAnalysisTaskCascadeTester::CheckDeltaDipAngle(const T* track1, const T* track2) const
{
   // Check if pair passes delta dip angle cut within 0 < p_t < fDeltaDipPt
   if(fDebug > 1) cout << " *** CheckDeltaDipAngle() *** " << endl;
   if ((TMath::Abs(DeltaDipAngle(track1, track2)) < fDeltaDipAngle) && (PhiPt(track1, track2) < fDeltaDipPt)) return kFALSE;
   return kTRUE;
}
*/
//_____________________________________________________________________________
template <typename T> Bool_t AliAnalysisTaskCascadeTester::CheckCandidateEtaPtCut(const T* track1, const T* track2) const
{
   // Check if pair passes eta and pt cut
   if(fDebug > 1) cout << " *** CheckCandidateEtaPtCut() *** " << endl;
   if (fCandidateMinPt > PhiPt(track1, track2) || fCandidateMaxPt < PhiPt(track1, track2)) return kFALSE;
   TVector3 a(track1->Px(), track1->Py(), track1->Pz());
   TVector3 b(track2->Px(), track2->Py(), track2->Pz());
   TVector3 c = a + b;
   if (fCandidateMinEta > c.Eta() || fCandidateMaxEta < c.Eta()) return kFALSE;
   return kTRUE;
}
//_____________________________________________________________________________
template <typename T> Bool_t AliAnalysisTaskCascadeTester::EventCut(T* event)
{
   // Impose event cuts
   if(fDebug > 0) cout << " *** EventCut() *** " << endl;
   if (!event) return kFALSE;
//   if (fSkipEventSelection) return kTRUE;
   if (!CheckVertex(event)) return kFALSE;
   if (!CheckCentrality(event)) return kFALSE;
   if(fQA) PlotMultiplcities(event);
   return kTRUE;
}
//_____________________________________________________________________________
template <typename T> void AliAnalysisTaskCascadeTester::PlotMultiplcities(const T* event) const
{
   // QA multiplicity plots
   if(fDebug > 1) cout << " *** PlotMultiplcities() *** " << endl;
   fVZEROA->Fill(event->GetVZEROData()->GetMTotV0A());
   fVZEROC->Fill(event->GetVZEROData()->GetMTotV0C());
   fTPCM->Fill(event->GetNumberOfTracks());
}
//_____________________________________________________________________________
template <typename T> Bool_t AliAnalysisTaskCascadeTester::CheckVertex(const T* event)
{
   // Check if event vertex is within given range
   if(fDebug > 0) cout << " *** CheckVertex() *** " << endl;
   if (!event->GetPrimaryVertex()) return 0x0;
   fVertex = event->GetPrimaryVertex()->GetZ();
   if (TMath::Abs(fVertex) > fVertexRange) return 0x0;
   return kTRUE;
}
//_____________________________________________________________________________
template <typename T> Bool_t AliAnalysisTaskCascadeTester::CheckCentrality(T* event)
{
   // Check if event is within the set centrality range. Falls back to V0 centrality determination if no method is set
   if(fDebug > 0) cout << " *** CheckCentrality() *** " << endl;
   if (!fkCentralityMethodA) AliFatal("No centrality method set! FATAL ERROR!");

   // check if the AliMultSelection object is present. If so, we should invoke the
   // new centrality framework

   AliMultSelection *multSelection = 0x0; 
   multSelection = static_cast<AliMultSelection*>(event->FindListObject("MultSelection"));
   if(multSelection) {
       fCentrality = multSelection->GetMultiplicityPercentile("V0M");
       if(fCentrality > fCentralityMin && fCentrality < fCentralityMax) {
           fCentralityPass->Fill(fCentrality);
           return kTRUE;
       } else {
           fCentralityNoPass->Fill(fCentrality);
           return kFALSE;
       }
   }

   else  fCentrality = event->GetCentrality()->GetCentralityPercentile(fkCentralityMethodA);
   Double_t cenB(-999);
   // if a second centrality estimator is requited, set it
   (fkCentralityMethodB) ? cenB = event->GetCentrality()->GetCentralityPercentile(fkCentralityMethodB) : cenB = fCentrality;
   if (TMath::Abs(fCentrality-cenB) > 5 || cenB >= 80 || cenB < 0 || fCentrality <= fCentralityMin || fCentrality > fCentralityMax) {
      if(fQA) fCentralityNoPass->Fill(fCentrality) ;
      return kFALSE;
   }
   const Int_t nGoodTracks = event->GetNumberOfTracks();
   if(fCentralityCut2010) { // cut on outliers
      Float_t multTPC(0.); // tpc mult estimate
      Float_t multGlob(0.); // global multiplicity
      for(Int_t iTracks = 0; iTracks < nGoodTracks; iTracks++) { // fill tpc mult
          AliAODTrack* trackAOD = dynamic_cast<AliAODTrack*>(event->GetTrack(iTracks));
          if(!trackAOD) AliFatal("Not a standard AOD");
          if (!trackAOD) continue;
          if (!(trackAOD->TestFilterBit(1))) continue;
          if ((trackAOD->Pt() < .2) || (trackAOD->Pt() > 5.0) || (TMath::Abs(trackAOD->Eta()) > .8) || (trackAOD->GetTPCNcls() < 70)  || (trackAOD->GetDetPid()->GetTPCsignal() < 10.0) || (trackAOD->Chi2perNDF() < 0.2)) continue;
          multTPC++;
      }
      for(Int_t iTracks = 0; iTracks < nGoodTracks; iTracks++) { // fill global mult
          AliAODTrack* trackAOD = dynamic_cast<AliAODTrack*>(event->GetTrack(iTracks));
          if(!trackAOD) AliFatal("Not a standard AOD");
          if (!trackAOD) continue;
          if (!(trackAOD->TestFilterBit(16))) continue;
          if ((trackAOD->Pt() < .2) || (trackAOD->Pt() > 5.0) || (TMath::Abs(trackAOD->Eta()) > .8) || (trackAOD->GetTPCNcls() < 70) || (trackAOD->GetDetPid()->GetTPCsignal() < 10.0) || (trackAOD->Chi2perNDF() < 0.1)) continue;
          Double_t b[2] = {-99., -99.};
          Double_t bCov[3] = {-99., -99., -99.};
          AliAODTrack copy(*trackAOD);
          if (!(copy.PropagateToDCA(event->GetPrimaryVertex(), event->GetMagneticField(), 100., b, bCov))) continue;
          if ((TMath::Abs(b[0]) > 0.3) || (TMath::Abs(b[1]) > 0.3)) continue;
          multGlob++;
      } //track loop
 //     printf(" mult TPC %.2f, mult Glob %.2f \n", multTPC, multGlob);
      if(! (multTPC > (-40.3+1.22*multGlob) && multTPC < (32.1+1.59*multGlob))) return kFALSE;
      if(fQA) {  
          fMultCorAfterCuts->Fill(multGlob, multTPC);  
          fMultvsCentr->Fill(fCentrality, multTPC);
      }
   }

   if(fCentralityCut2011) { // cut on outliers
      Float_t multTPC(0.); // tpc mult estimate
      Float_t multGlob(0.); // global multiplicity
      for(Int_t iTracks = 0; iTracks < nGoodTracks; iTracks++) { // fill tpc mult
          AliAODTrack* trackAOD = dynamic_cast<AliAODTrack*>(event->GetTrack(iTracks));
          if(!trackAOD) AliFatal("Not a standard AOD");
          if (!trackAOD) continue;
          if (!(trackAOD->TestFilterBit(1))) continue;
          if ((trackAOD->Pt() < .2) || (trackAOD->Pt() > 5.0) || (TMath::Abs(trackAOD->Eta()) > .8) || (trackAOD->GetTPCNcls() < 70)  || (trackAOD->GetDetPid()->GetTPCsignal() < 10.0) || (trackAOD->Chi2perNDF() < 0.2)) continue;
          multTPC++;
      }
      for(Int_t iTracks = 0; iTracks < nGoodTracks; iTracks++) { // fill global mult
          AliAODTrack* trackAOD = dynamic_cast<AliAODTrack*>(event->GetTrack(iTracks));
          if(!trackAOD) AliFatal("Not a standard AOD");
          if (!trackAOD) continue;
          if (!(trackAOD->TestFilterBit(16))) continue;
          if ((trackAOD->Pt() < .2) || (trackAOD->Pt() > 5.0) || (TMath::Abs(trackAOD->Eta()) > .8) || (trackAOD->GetTPCNcls() < 70) || (trackAOD->GetDetPid()->GetTPCsignal() < 10.0) || (trackAOD->Chi2perNDF() < 0.1)) continue;
          Double_t b[2] = {-99., -99.};
          Double_t bCov[3] = {-99., -99., -99.};
          AliAODTrack copy(*trackAOD);
          if (!(copy.PropagateToDCA(event->GetPrimaryVertex(), event->GetMagneticField(), 100., b, bCov))) continue;
          if ((TMath::Abs(b[0]) > 0.3) || (TMath::Abs(b[1]) > 0.3)) continue;
          multGlob++;
      } //track loop
      //printf(" mult TPC %.2f, mult Glob %.2f \n", multTPC, multGlob);
      if(! (multTPC > (-36.73 + 1.48*multGlob) && multTPC < (62.87 + 1.78*multGlob))) return kFALSE;
      if(fQA) {  
          fMultCorAfterCuts->Fill(multGlob, multTPC);  
          fMultvsCentr->Fill(fCentrality, multTPC);
      }
   }

   fCentralityPass->Fill(fCentrality);
   return kTRUE;
}
//_____________________________________________________________________________
template <typename T> Double_t AliAnalysisTaskCascadeTester::PhiPt(const T* track1, const T* track2) const
{
   // return p_t of track pair
   TVector3 a(track1->Px(), track1->Py(), track1->Pz());
   TVector3 b(track2->Px(), track2->Py(), track2->Pz());
   TVector3 c = a + b;
   return c.Pt();
}
//_____________________________________________________________________________
template <typename T> void AliAnalysisTaskCascadeTester::PtSelector(Int_t tracktype, const T* track1, const T* track2) const
{
   // plot m_inv spectra of like- and unlike-sign kaon pairs for each pt bin
   Double_t pt = PhiPt(track1, track2);
   if (tracktype == 0) {
       for(Int_t ptbin(0); ptbin < fNPtBins; ptbin++) {
           if ((fPtBins[ptbin] <= pt) && (fPtBins[ptbin+1] > pt )) {
               fInvMNP[ptbin]->Fill(InvariantMass(track1, track2));
               if(fQA) fPtSpectra[ptbin]->Fill(pt);
           }
       }
   }
   if (tracktype == 1) {
       for(Int_t ptbin(0); ptbin < fNPtBins; ptbin++) {
           if ((fPtBins[ptbin] <= pt) && (fPtBins[ptbin+1] > pt )) {
               fInvMPP[ptbin]->Fill(InvariantMass(track1, track2));
           }
       }
   }
   if (tracktype == 2) {
       for(Int_t ptbin(0); ptbin < fNPtBins; ptbin++) {
           if ((fPtBins[ptbin] <= pt) && (fPtBins[ptbin+1] > pt )) {
               fInvMNN[ptbin]->Fill(InvariantMass(track1, track2));
           }
       }
   }
}
//_____________________________________________________________________________
template <typename T> Bool_t AliAnalysisTaskCascadeTester::PhiTrack(T* track) const
{
   // check if track meets quality cuts
   if(!track) return kFALSE;
   return fPOICuts->IsSelected(track);
}
//_____________________________________________________________________________
template <typename T> void AliAnalysisTaskCascadeTester::SetNullCuts(T* event)
{
   // Set null cuts
   if (fDebug > 0) cout << " *** SetNullCuts() *** " << fCutsRP << endl;
   fCutsRP->SetEvent(event, MCEvent());
   fNullCuts->SetParamType(AliFlowTrackCuts::kGlobal);
   fNullCuts->SetPtRange(+1, -1); // select nothing QUICK
   fNullCuts->SetEtaRange(+1, -1); // select nothing VZERO
   fNullCuts->SetEvent(event, MCEvent());
}
//_____________________________________________________________________________
void AliAnalysisTaskCascadeTester::PrepareFlowEvent(Int_t iMulti)
{
   // Prepare flow events
   if (fDebug > 0 ) cout << " *** PrepareFlowEvent() *** " << endl;
   fFlowEvent->ClearFast();
   fFlowEvent->Fill(fCutsRP, fNullCuts);
   fFlowEvent->SetReferenceMultiplicity(iMulti);
   fFlowEvent->DefineDeadZone(0, 0, 0, 0);
}
//_____________________________________________________________________________
void AliAnalysisTaskCascadeTester::UserExec(Option_t *)
{
   // UserExec: called for each event. Commented where necessary
   if(fDebug > 0 ) cout << " *** UserExec() *** " << endl;
   TObjArray* MixingCandidates = 0x0; // init null pointer for event mixing
   if(fEventMixing) {
       MixingCandidates = new TObjArray();
       MixingCandidates->SetOwner(kTRUE); // mixing candidates has ownership of objects in array
   }
   if (!fPIDResponse) {
      if(fDebug > 0 ) cout << " Could not get PID response " << endl;
      return;
   }
   fAOD = dynamic_cast<AliAODEvent*>(InputEvent()); // check for aod data type
   if (fAOD) {
      if (!EventCut(fAOD)) return; // check for event cuts


      // add extra pileup cuts for high intensity runs
      // courtesy of alex dobrin
      if (fPileUp){

          if (plpMV(fAOD))
              return;
	
          Short_t isPileup = fAOD->IsPileupFromSPD(3);
          if (isPileup != 0)
              return;
	
          if (((AliAODHeader*)fAOD->GetHeader())->GetRefMultiplicityComb08() < 0)
              return;
	
          /*
          Int_t bc2 = ((AliAODHeader*)fAOD->GetHeader())->GetIRInt2ClosestInteractionMap();
          if (bc2 != 0)
              return;
	
          Int_t bc1 = ((AliAODHeader*)fAOD->GetHeader())->GetIRInt1ClosestInteractionMap();
          if (bc1 != 0)
              return;
           */
	
      }
       
      //new function for 2015 to remove incomplete events
      if (fAOD->IsIncompleteDAQ())
          return;


      SetNullCuts(fAOD);
      PrepareFlowEvent(fAOD->GetNumberOfTracks());
      fCandidates->SetLast(-1);
      if(fIsMC) IsMC(); // launch mc stuff
      if(fQA) fEventStats->Fill(0);


      const Int_t ncasc = fAOD->GetNumberOfCascades();
      Double_t  lPrimaryVtxPosition[3];
      fAOD->GetPrimaryVertex()->GetXYZ(lPrimaryVtxPosition);

       //cascades vn
      for (Int_t jCasc = 0; jCasc < ncasc; jCasc++) {
          
          AliAODcascade* casc = (AliAODcascade*) fAOD->GetCascade(jCasc);
          
          if (!casc){
              delete casc;
              continue;
          }
                      
          Double_t lDcaXiDaughters          = casc->DcaXiDaughters();
          Double_t lDcaBachToPrimVertexXi   = casc->DcaBachToPrimVertex();
          Double_t lXiCosineOfPointingAngle = casc->CosPointingAngleXi(lPrimaryVtxPosition[0], lPrimaryVtxPosition[1], lPrimaryVtxPosition[2]);
          
          Double_t lPosXi[3] = { -1000.0, -1000.0, -1000.0 };
          lPosXi[0] = casc->DecayVertexXiX();
          lPosXi[1] = casc->DecayVertexXiY();
          lPosXi[2] = casc->DecayVertexXiZ();
          
          Double_t lDcaV0DaughtersXi = casc->DcaV0Daughters();
          Double_t lV0toXiCosineOfPointingAngle = casc->CosPointingAngle(lPosXi);
          
          Double_t lPosV0Xi[3] = { -1000.0, -1000.0, -1000.0 };
          lPosV0Xi[0] = casc->DecayVertexV0X();
          lPosV0Xi[1] = casc->DecayVertexV0Y();
          lPosV0Xi[2] = casc->DecayVertexV0Z();
          
          Double_t lDcaV0ToPrimVertexXi = casc->DcaV0ToPrimVertex();
          Double_t lDcaPosToPrimVertexXi = casc->DcaPosToPrimVertex();
          Double_t lDcaNegToPrimVertexXi = casc->DcaNegToPrimVertex();
          
          Double_t lInvMassLambdaAsCascDghter = casc->MassLambda();
          
          Double_t lXiRadius = TMath::Sqrt(lPosXi[0]*lPosXi[0] + lPosXi[1]*lPosXi[1]);
          Double_t lV0RadiusXi = TMath::Sqrt(lPosV0Xi[0]*lPosV0Xi[0] + lPosV0Xi[1]*lPosV0Xi[1]);

                      
          if (lDcaXiDaughters > 0.3)
              continue; // in AliCascadeVertexer: max allowed DCA between the V0 and the bachelor
          
          if (lXiCosineOfPointingAngle < 0.999 )
              continue; // in AliCascadeVertexer: min allowed cosine of the cascade pointing angle
          
          if (lDcaV0ToPrimVertexXi < 0.05)
              continue; // in AliCascadeVertexer: min allowed V0 impact parameter
          
          if (lDcaBachToPrimVertexXi < 0.03)
              continue; // in AliCascadeVertexer: min allowed bachelor's impact parameter
          
          if (TMath::Abs(lInvMassLambdaAsCascDghter-1.11568) > 0.008 )
              continue;  // in AliCascadeVertexer: "window" around the Lambda mass
          
          if(lXiRadius < 0.9)
              continue; // in AliCascadeVertexer: min radius of the fiducial volume
          
          if(lXiRadius > 100.)
              continue; // in AliCascadeVertexer: max radius of the fiducial volume
     
          
          if (lDcaV0DaughtersXi > 1.)
              continue; // in AliV0vertexer: max allowed DCA between the daughter tracks
          
          if (lV0toXiCosineOfPointingAngle < 0.998)
              continue; // in AliV0vertexer: min allowed cosine of V0's pointing angle
          
          if (lDcaPosToPrimVertexXi < 0.1)
              continue; // in AliV0vertexer: min allowed impact parameter for the 1st daughter
          
          if (lDcaNegToPrimVertexXi < 0.1)
              continue; // in AliV0vertexer: min allowed impact parameter for the 2nd daughter
          
          if(lV0RadiusXi < 0.9)
              continue; // in AliV0vertexer: min radius of the fiducial volume
          
          if(lV0RadiusXi > 100.)
              continue; // in AliV0vertexer: max radius of the fiducial volume

                      
          if ((casc->Pt() < 0.2) || (casc->Pt() >= 8.) || (TMath::Abs(casc->Eta()) >= 0.8))
              continue;
          
          
          
          AliAODTrack* pTrackXi    = (AliAODTrack*)casc->GetDaughter(0);
          AliAODTrack* nTrackXi    = (AliAODTrack*)casc->GetDaughter(1);
          AliAODTrack* bachTrackXi = (AliAODTrack*)casc->GetDecayVertexXi()->GetDaughter(0);
          
          if (!pTrackXi || !nTrackXi || !bachTrackXi)
              continue;
          
          
          UInt_t lIdxPosXi  = (UInt_t) TMath::Abs(pTrackXi->GetID());
          UInt_t lIdxNegXi  = (UInt_t) TMath::Abs(nTrackXi->GetID());
          UInt_t lBachIdx   = (UInt_t) TMath::Abs(bachTrackXi->GetID());
          if(lBachIdx == lIdxNegXi)
              continue;
          
          if(lBachIdx == lIdxPosXi)
              continue;
          
          
          if (!(pTrackXi->IsOn(AliAODTrack::kTPCrefit)))
              continue;
          
          if (!(nTrackXi->IsOn(AliAODTrack::kTPCrefit)))
              continue;
      
          if (!(bachTrackXi->IsOn(AliAODTrack::kTPCrefit)))
              continue;
      
          
          if (pTrackXi->GetTPCNcls() < 70 || nTrackXi->GetTPCNcls() < 70 || bachTrackXi->GetTPCNcls() < 70)
              continue;
          
          if (TMath::Abs(pTrackXi->Eta()) > 0.8 || (TMath::Abs(nTrackXi->Eta()) > 0.8 || TMath::Abs(bachTrackXi->Eta()) > 0.8))
              continue;
              
          if (pTrackXi->Pt() < 0.15 || nTrackXi->Pt() < 0.15 || bachTrackXi->Pt() < 0.15)
              continue;
          
          
          //PID cuts
          if (bachTrackXi->GetTPCsignalN() < 70 || nTrackXi->GetTPCsignalN() < 70 || pTrackXi->GetTPCsignalN() < 70)
              continue;
          
          // Bachelor
          Double_t nSigmaBachK = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(bachTrackXi, AliPID::kKaon));
          Double_t nSigmaBachPi = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(bachTrackXi, AliPID::kPion));
          
          // Negative V0 daughter
          Double_t nSigmaPiNeg = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(nTrackXi, AliPID::kPion));
          Double_t nSigmaPNeg = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(nTrackXi,AliPID::kProton));
          
          // Positive V0 daughter
          Double_t nSigmaPiPos = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(pTrackXi,AliPID::kPion));
          Double_t nSigmaPPos = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(pTrackXi,AliPID::kProton));

          
          Double_t lInvMassOmega = casc->MassOmega();
          Double_t lInvMassXi = casc->MassXi();
          
          if (lInvMassOmega >  1.64 && lInvMassOmega < 1.71 && nSigmaBachK < 3.){
              if ((nSigmaPiNeg < 3. && nSigmaPPos < 3.) || (nSigmaPNeg < 3. && nSigmaPiPos < 3.)){
              if(fOmega) MakeTrack(lInvMassOmega, casc->Pt(), casc->Phi(), casc->Eta());
              }
              
          }
          
          if (lInvMassXi > 1.29 && lInvMassXi < 1.36 && nSigmaBachPi < 3.){
              
              if ((nSigmaPiNeg < 3. && nSigmaPPos < 3.) || (nSigmaPNeg < 3. && nSigmaPiPos < 3.)){
              
              if(fXi) MakeTrack(lInvMassXi,casc->Pt(), casc->Phi(), casc->Eta());
              //fill Xi invariant mass and v2 vs invariant mass hists
              }
              
          }
      }



      if (fDebug > 0)  printf("I received %d candidates\n", fCandidates->GetEntriesFast()); // insert candidates into flow events
      for (int iCand = 0; iCand != fCandidates->GetEntriesFast(); ++iCand) {
         AliFlowCandidateTrack *cand = dynamic_cast<AliFlowCandidateTrack*>(fCandidates->At(iCand));
         if (!cand) continue;
         if (fDebug > 1) printf(" --> Checking at candidate %d with %d daughters: mass %f\n", iCand, cand->GetNDaughters(), cand->Mass());
         for (int iDau = 0; iDau != cand->GetNDaughters(); ++iDau) {
            if (fDebug>1) printf("     *** Daughter %d with fID %d ***", iDau, cand->GetIDDaughter(iDau));
            for (int iRPs = 0; iRPs != fFlowEvent->NumberOfTracks(); ++iRPs) {
               AliFlowTrack *iRP = dynamic_cast<AliFlowTrack*>(fFlowEvent->GetTrack(iRPs));
               if (!iRP) continue;
               if (!iRP->InRPSelection()) continue;
               if (cand->GetIDDaughter(iDau) == iRP->GetID()) {
                  if (fDebug > 1) printf("      was in RP set");
                  iRP->SetForRPSelection(kFALSE);
                  fFlowEvent->SetNumberOfRPs(fFlowEvent->GetNumberOfRPs() - 1);
               }
            }
            if (fDebug > 1) printf("\n");
         }
         cand->SetForPOISelection(kTRUE);
         fFlowEvent->InsertTrack(((AliFlowTrack*) cand));
      }
      if (fDebug > 0) printf("TPCevent %d\n", fFlowEvent->NumberOfTracks());
  
  /*    if(!fEventMixing) { // combinatorial background
          for (Int_t pTracks = 0; pTracks < unp ; pTracks++) {
             for (Int_t nTracks = pTracks + 1; nTracks < unp ; nTracks++) {
//                if (fApplyDeltaDipCut && (!CheckDeltaDipAngle(up[pTracks], up[nTracks]))) continue;
                if (fCandidateEtaPtCut && (!CheckCandidateEtaPtCut(up[pTracks], up[nTracks]))) continue;
                PtSelector(1, up[pTracks], up[nTracks]);
             }
          }
          for (Int_t nTracks = 0; nTracks < unn ; nTracks++) {
             for (Int_t pTracks = nTracks + 1; pTracks < unn ; pTracks++) {
//                if (fApplyDeltaDipCut && (!CheckDeltaDipAngle(un[nTracks], un[pTracks]))) continue;
                if (fCandidateEtaPtCut && (!CheckCandidateEtaPtCut(un[nTracks], un[pTracks]))) continue;
                PtSelector(2, un[nTracks], un[pTracks]);
             }
          }
      }
      if(fEventMixing) ReconstructionWithEventMixing(MixingCandidates);
      */
      PostData(1, fOutputList);
      PostData(2, fFlowEvent);
   }
}
//_____________________________________________________________________________
void AliAnalysisTaskCascadeTester::Exec(Option_t* c)
{
    // skip the event selection for SE task (e.g. for MC productions)
    if(fSkipEventSelection) AliAnalysisTaskCascadeTester::UserExec(c);
    // exec of task se will do event selection and call UserExec 
    else AliAnalysisTaskSE::Exec(c);
}
//_____________________________________________________________________________
void AliAnalysisTaskCascadeTester::ReconstructionWithEventMixing(TObjArray* MixingCandidates) const
{
    // perform phi reconstruction with event mixing
    if(fDebug > 0) cout << " *** ReconstructionWithEventMixing() *** " << endl;
    AliEventPool* pool = fPoolManager->GetEventPool(fCentrality, fVertex);
    if(!pool) AliFatal(Form("No pool found for centrality = %f, zVtx = %f", fCentrality, fVertex));
    if(pool->IsReady() || pool->NTracksInPool() > fMixingParameters[1] / 10 || pool->GetCurrentNEvents() >= fMixingParameters[2]) {
        Int_t nEvents = pool->GetCurrentNEvents();
        if(fDebug > 0) cout << " --> " << nEvents << " events in mixing buffer ... " << endl;
        for (Int_t iEvent(0); iEvent < nEvents; iEvent++) {
            TObjArray* mixed_candidates = pool->GetEvent(iEvent);
            if(!mixed_candidates) continue; // this should NEVER happen
            Int_t bufferTracks = mixed_candidates->GetEntriesFast(); // buffered candidates
            Int_t candidates = MixingCandidates->GetEntriesFast(); // mixing candidates
            if(fDebug > 0) cout << Form("   - mixing %d tracks with %d buffered tracks from event %d ... ", candidates, bufferTracks, iEvent) << endl;
            AliPhiMesonHelperTrack* buffer_un[bufferTracks]; // set values for buffered candidates
            AliPhiMesonHelperTrack* buffer_up[bufferTracks];
            Int_t buffer_unp(0);
            Int_t buffer_unn(0);
            AliPhiMesonHelperTrack* mix_un[candidates];// set values for mixing candidates
            AliPhiMesonHelperTrack* mix_up[candidates];
            Int_t mix_unp(0);
            Int_t mix_unn(0);
            for (Int_t iTracks = 0; iTracks < candidates; iTracks++) { // distinguish between k+ and k- for mixing candidates
                AliPhiMesonHelperTrack* track = (AliPhiMesonHelperTrack*)MixingCandidates->At(iTracks);
                if(!track) continue;
                if (track->Charge() > 0) {
                    mix_up[mix_unp] = track;
                    mix_unp++;
                }
                else if (track->Charge() < 0 ) {
                   mix_un[mix_unn] = track;
                   mix_unn++;
                }
            }
            for (Int_t iTracks = 0; iTracks < bufferTracks; iTracks++) { // distinguish between k+ and k- for buffered candidates
                AliPhiMesonHelperTrack* track = (AliPhiMesonHelperTrack*)mixed_candidates->At(iTracks);
                if(!track) continue;
                if (track->Charge() > 0) {
                    buffer_up[buffer_unp] = track;
                    buffer_unp++;
                }
                else if (track->Charge() < 0 ) {
                   buffer_un[buffer_unn] = track;
                   buffer_unn++;
                }
            }
            for (Int_t pMix = 0; pMix < mix_unp; pMix++) { // mix k+ (current event) k+ (buffer)
                if(fDebug > 1 ) cout << Form("mixing current k+ track %d with", pMix);
                if(!fTypeMixing) { // mix with like-sign kaons
                    for(Int_t pBuf = 0; pBuf < buffer_unp; pBuf++) {
                        if(fDebug > 1 ) cout << Form(" buffered k+ track %d", pBuf) << endl;
//                        if (fApplyDeltaDipCut && (!CheckDeltaDipAngle(mix_up[pMix], buffer_up[pBuf]))) continue;
                        if (fCandidateMinEta && (!CheckCandidateEtaPtCut(mix_up[pMix], buffer_up[pBuf]))) continue;
                        PtSelector(1, mix_up[pMix], buffer_up[pBuf]);
                    }
                }
                else if(fTypeMixing) { // mix with unlike-sign kaons
                    for(Int_t nBuf = 0; nBuf < buffer_unn; nBuf++) {
                        if(fDebug > 1 ) cout << Form(" buffered k- track %d", nBuf) << endl;
//                        if (fApplyDeltaDipCut && (!CheckDeltaDipAngle(mix_up[pMix], buffer_un[nBuf]))) continue;
                        if (fCandidateMinEta && (!CheckCandidateEtaPtCut(mix_up[pMix], buffer_un[nBuf]))) continue;
                        PtSelector(2, mix_up[pMix], buffer_un[nBuf]);
                    }
                }
            }
            for (Int_t nMix = 0; nMix < mix_unn; nMix++) { // mix k- (current event) k- (buffer)
                if(fDebug > 1 ) cout << Form("mixing current k- track %d with", nMix);
                if(!fTypeMixing) { // mix with like-sign kaons
                    for(Int_t nBuf = 0; nBuf < buffer_unn; nBuf++) {
                        if(fDebug > 1 ) cout << Form(" buffered k- track %d", nBuf) << endl;
//                        if (fApplyDeltaDipCut && (!CheckDeltaDipAngle(mix_un[nMix], buffer_un[nBuf]))) continue;
                        if (fCandidateMinEta && (!CheckCandidateEtaPtCut(mix_un[nMix], buffer_un[nBuf]))) continue;
                        PtSelector(2, mix_un[nMix], buffer_un[nBuf]);
                    }
                }
                else if(fTypeMixing) { // mix with unlike-sign kaons
                    for(Int_t pBuf = 0; pBuf < buffer_unp; pBuf++) {
                        if(fDebug > 1 ) cout << Form(" buffered k+ track %d", pBuf) << endl;
//                        if (fApplyDeltaDipCut && (!CheckDeltaDipAngle(mix_un[nMix], buffer_up[pBuf]))) continue;
                        if (fCandidateMinEta && (!CheckCandidateEtaPtCut(mix_un[nMix], buffer_up[pBuf]))) continue;
                        PtSelector(1, mix_un[nMix], buffer_up[pBuf]);
                    }
                }
            }
        } // end of mixed events loop
    } // end of checking to see whether pool is filled correctly
    pool->UpdatePool(MixingCandidates); // update pool with current mixing candidates (for next event)
    if(fDebug > 0 ) pool->PrintInfo();
}
//_____________________________________________________________________________
void AliAnalysisTaskCascadeTester::Terminate(Option_t *)
{
    // terminate
    if(fDebug > 0) cout << " *** Terminate() *** " << endl;
}
//______________________________________________________________________________
void  AliAnalysisTaskCascadeTester::MakeTrack(Double_t mass, Double_t pt, Double_t phi, Double_t eta) const
{
   // Construct Flow Candidate Track from two selected candidates
   if(fDebug > 1 ) cout << " *** MakeTrack() *** " << endl;
   // if requested, check rapidity at this point
   Bool_t overwrite = kTRUE;
   AliFlowCandidateTrack *sTrack = static_cast<AliFlowCandidateTrack*>(fCandidates->At(fCandidates->GetLast() + 1));
   if (!sTrack) {
      sTrack = new AliFlowCandidateTrack(); //deleted by fCandidates
      overwrite = kFALSE;
   }
   else sTrack->ClearMe();
   sTrack->SetMass(mass);
   sTrack->SetPt(pt);
   sTrack->SetPhi(phi);
   sTrack->SetEta(eta);
   sTrack->SetForPOISelection(kTRUE);
   sTrack->SetForRPSelection(kFALSE);
   if (overwrite) fCandidates->SetLast(fCandidates->GetLast() + 1);
   else fCandidates->AddLast(sTrack);
   return;
}
//_____________________________________________________________________________
void AliAnalysisTaskCascadeTester::IsMC()
{
    // Fill QA histos for MC analysis
   TClonesArray *arrayMC = 0;
   if(fDebug > 0) cout << " -> Switching to MC mode <- " << endl;
   // fill array with mc tracks 
   arrayMC = (TClonesArray*) fAOD->GetList()->FindObject(AliAODMCParticle::StdBranchName());
   if (!arrayMC) AliFatal("Error: MC particles branch not found!\n");
   for (Int_t iTracks = 0; iTracks < fAOD->GetNumberOfTracks(); iTracks++) {
     AliAODTrack* track = dynamic_cast<AliAODTrack*>(fAOD->GetTrack(iTracks));
     if(!track) AliFatal("Not a standard AOD");
     if (fDebug>1) cout << " Received MC kaon " << endl;
     Double_t b[2] = { -99., -99.};
     Double_t bCov[3] = { -99., -99., -99.};
     AliAODTrack copy(*track);
     if(!copy.PropagateToDCA(fAOD->GetPrimaryVertex(), fAOD->GetMagneticField(), 100., b, bCov)) return;
     // find corresponding mc particle
     AliAODMCParticle *partMC = (AliAODMCParticle*) arrayMC->At(TMath::Abs(track->GetLabel()));
     if (!partMC) {
         if(fDebug > 1) cout << "Cannot get MC particle" << endl;
         continue;
     }
     // Check if it is primary, secondary from material or secondary from weak decay
     Bool_t isPrimary           = partMC->IsPhysicalPrimary();
     Bool_t isSecondaryMaterial = kFALSE;
     Bool_t isSecondaryWeak     = kFALSE;
     if (!isPrimary) {
         Int_t mfl = -999, codemoth = -999;
         Int_t indexMoth = partMC->GetMother();
         if (indexMoth >= 0) { // is not fake
            AliAODMCParticle* moth = (AliAODMCParticle*) arrayMC->At(indexMoth);
            codemoth = TMath::Abs(moth->GetPdgCode());
            mfl = Int_t(codemoth / TMath::Power(10, Int_t(TMath::Log10(codemoth))));
         }
         if (mfl == 3) isSecondaryWeak     = kTRUE;
         else       isSecondaryMaterial = kTRUE;
      }
      if (isPrimary) {
          fDCAPrim->Fill(track->Pt(), b[0]);
          fDCAXYQA->Fill(b[0]);
          fDCAZQA->Fill(b[1]);
      }
      if (isSecondaryWeak)  fDCASecondaryWeak->Fill(track->Pt(), b[0]);
      if (isSecondaryMaterial) fDCAMaterial->Fill(track->Pt(), b[0]);
   }
}



//_____________________________________________________________________________
Double_t AliAnalysisTaskCascadeTester::GetWDist(const AliVVertex* v0, const AliVVertex* v1)
{
    
    // calculate sqrt of weighted distance to other vertex
    if (!v0 || !v1) {
        printf("One of vertices is not valid\n");
        return 0;
    }
    static TMatrixDSym vVb(3);
    double dist = -1;
    double dx = v0->GetX()-v1->GetX();
    double dy = v0->GetY()-v1->GetY();
    double dz = v0->GetZ()-v1->GetZ();
    double cov0[6],cov1[6];
    v0->GetCovarianceMatrix(cov0);
    v1->GetCovarianceMatrix(cov1);
 //
 // fQxavsV0[0]          fQxnmV0A
 // fQyavsV0[0]          fQynmV0A
 // fQxavsV0[1]          fQxnsV0A
 // fQyavsV0[1]          fQynsV0A
 // fQxcvsV0[0]          fQxnmV0C
 // fQycvsV0[0]          fQynmV0C
 // fQxcvsV0[1]          fQxnsV0C
 // fQycvsV0[1]          fQynsV0C    vVb(0,0) = cov0[0]+cov1[0];
    vVb(1,1) = cov0[2]+cov1[2];
    vVb(2,2) = cov0[5]+cov1[5];
    vVb(1,0) = vVb(0,1) = cov0[1]+cov1[1];
    vVb(0,2) = vVb(1,2) = vVb(2,0) = vVb(2,1) = 0.;
    vVb.InvertFast();
    if (!vVb.IsValid()) {printf("Singular Matrix\n"); return dist;}
    dist = vVb(0,0)*dx*dx + vVb(1,1)*dy*dy + vVb(2,2)*dz*dz
    +    2*vVb(0,1)*dx*dy + 2*vVb(0,2)*dx*dz + 2*vVb(1,2)*dy*dz;
    return dist>0 ? TMath::Sqrt(dist) : -1;
    
}




//_____________________________________________________________________________
Bool_t AliAnalysisTaskCascadeTester::plpMV(const AliAODEvent* aod)
{
    // check for multi-vertexer pile-up
    //
    const int    kMinPlpContrib = 5;
    const double kMaxPlpChi2 = 5.0;
    const double kMinWDist = 15;
    //
    const AliVVertex* vtPrm = 0;
    const AliVVertex* vtPlp = 0;
    int nPlp = 0;
    //
    if ( !(nPlp=aod->GetNumberOfPileupVerticesTracks()) ) return kFALSE;
    vtPrm = aod->GetPrimaryVertex();
    if (vtPrm == aod->GetPrimaryVertexSPD()) return kTRUE; // there are pile-up vertices but no primary
    
    //int bcPrim = vtPrm->GetBC();
    //
    for (int ipl=0;ipl<nPlp;ipl++) {
        vtPlp = (const AliVVertex*)aod->GetPileupVertexTracks(ipl);
        //
        if (vtPlp->GetNContributors() < kMinPlpContrib) continue;
        if (vtPlp->GetChi2perNDF() > kMaxPlpChi2) continue;
        //  int bcPlp = vtPlp->GetBC();
        //  if (bcPlp!=AliVTrack::kTOFBCNA && TMath::Abs(bcPlp-bcPrim)>2) return kTRUE; // pile-up from other BC
        //
        double wDst = GetWDist(vtPrm,vtPlp);
        if (wDst<kMinWDist) continue;
        //
        return kTRUE; // pile-up: well separated vertices
    }
    //
    return kFALSE;
    //
}

//_____________________________________________________________________________
//---------------------------------------------------
Short_t AliAnalysisTaskCascadeTester::FindMinNSigma(Double_t nSpi, Double_t nSk, Double_t nSp) const
{
    
    Short_t kPID = 0;
    
    if((nSk == nSpi) && (nSk == nSp))
        return kPID;
    if((nSk < nSpi) && (nSk < nSp) && (nSk < fPIDConfig[0]))
        kPID = 2;
    
    if((nSpi < nSk) && (nSpi < nSp) && (nSpi < fPIDConfig[0]))
        kPID = 1;
    
    if((nSp < nSk) && (nSp < nSpi) && (nSp < fPIDConfig[0]))
        kPID = 3;
    return kPID;
    
}


//_____________________________________________________________________________
Bool_t AliAnalysisTaskCascadeTester::GetDoubleCountingK(Double_t nSk, Short_t minNSigma) const
{
    
    Bool_t kDC = kFALSE;
    
    if (nSk < fPIDConfig[0] && minNSigma != 2)
        kDC = kTRUE;
    
    return kDC;
    
}

