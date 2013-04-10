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
// AliAnalysisTwoParticleResonanceFlowTask:
// author: Redmer Alexander Bertens (rbertens@cern.ch)
//         You Zhou                 (you.zhou@cern.ch) 
// analyis task for 
// AliResonanceFlowHelperTrack provides a lightweight helper track for reconstruction
// This analysis DOES NOT support ESDs
// Version: v2012.9.7


#include "TChain.h"
#include "TTree.h"
#include "TH1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TObject.h"
#include "TObjArray.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliCentrality.h"
#include "AliVEvent.h"
#include "AliAnalysisTwoParticleResonanceFlowTask.h"
#include "AliFlowBayesianPID.h"
#include "AliStack.h"
#include "AliMCEvent.h"
#include "TProfile.h"
#include "TDirectoryFile.h"
#include "AliFlowCandidateTrack.h"
#include "AliFlowTrackCuts.h"
#include "AliFlowEventSimple.h"
#include "AliFlowTrackSimple.h"
#include "AliFlowCommonConstants.h"
#include "AliFlowEvent.h"
#include "TVector3.h"
#include "AliAODVZERO.h"
#include "AliPID.h"
#include "AliPIDResponse.h"
#include "AliAODMCParticle.h"
#include "AliAnalysisTaskVnV0.h"
#include "AliEventPoolManager.h"

class AliFlowTrackCuts;

using std::cout;
using std::endl;

ClassImp(AliAnalysisTwoParticleResonanceFlowTask)

AliAnalysisTwoParticleResonanceFlowTask::AliAnalysisTwoParticleResonanceFlowTask() : AliAnalysisTaskSE(),
  fSpeciesA(0), fSpeciesB(0), fChargeA(0), fChargeB(0), fMassA(0), fMassB(0), fMinPtA(0), fMaxPtA(0), fMinPtB(0), fMaxPtB(0), fIsMC(0), fEventMixing(0), fPhiMinusPsiMethod(0), fQA(0), fV0(0), fMassBins(1), fMinMass(-1.), fMaxMass(0.), fCutsRP(NULL), fNullCuts(0), fPIDResponse(0), fFlowEvent(0), fBayesianResponse(0), fCandidates(0),  fCandidateEtaPtCut(0), fCandidateMinEta(0), fCandidateMaxEta(0), fCandidateMinPt(0), fCandidateMaxPt(0), fNPtBins(18), fNdPhiBins(18), fCentrality(999), fVertex(999), fAOD(0), fPoolManager(0), fOutputList(0), fEventStats(0), fCentralityPass(0), fCentralityNoPass(0), fNOPID(0), fPIDk(0), fPIDp(0), fPtP(0), fPtN(0), fPtSpeciesA(0), fPtSpeciesB(0), fMultCorAfterCuts(0), fMultvsCentr(0), fCentralityMin(0), fCentralityMax(100), fkCentralityMethodA(0), fkCentralityMethodB(0), fCentralityCut2010(0), fCentralityCut2011(0), fPOICuts(0), fVertexRange(0), fPhi(0), fEta(0), fVZEROA(0), fVZEROC(0), fTPCM(0), fDeltaDipAngle(0), fDeltaDipPt(0), fApplyDeltaDipCut(0), fDCAAll(0), fDCAXYQA(0), fDCAZQA(0), fDCAPrim(0), fDCASecondaryWeak(0), fDCAMaterial(0), fSubEventDPhiv2(0), fAnalysisSummary(0)
{
   // Default constructor
   for(Int_t i(0); i < 7; i++) fPIDConfig[i] = 0.;
   for(Int_t i(0); i < 5; i++) fDCAConfig[i] = 0.;
   for(Int_t i(0); i < 20; i++) {
       fVertexMixingBins[i] = 0;
       fCentralityMixingBins[i] = 0;
   }
   fMixingParameters[0] = 1000; fMixingParameters[1] = 50000; fMixingParameters[2] = 5;
   for(Int_t i(0); i < 18; i++) {
       for(Int_t j(0); j < 2; j++) {
           fV0Data[i][j] = 0;
           for(Int_t k(0); k < 15; k++) fPhiMinusPsiDataContainer[k][i][j] = 0x0;
           for(Int_t k(0); k < 15; k++) fPhiMinusPsiBackgroundContainer[k][i][j] = 0x0;
       }
       fResonanceSignal[i] = 0; fResonanceBackground[i] = 0; fPtSpectra[i] = 0; fPtBins[i] = 0.; fdPhiBins[i] = 0.;
   }
   for(Int_t i(0); i < 12; i++) fAddTaskMacroSummary[i] = 0.;
}
//_____________________________________________________________________________
AliAnalysisTwoParticleResonanceFlowTask::AliAnalysisTwoParticleResonanceFlowTask(const char *name) : AliAnalysisTaskSE(name), fSpeciesA(0), fSpeciesB(0), fChargeA(0), fChargeB(0), fMassA(0), fMassB(0), fMinPtA(0), fMaxPtA(0), fMinPtB(0), fMaxPtB(0), fIsMC(0), fEventMixing(0), fPhiMinusPsiMethod(0), fQA(0), fV0(0), fMassBins(1), fMinMass(-1.), fMaxMass(0.), fCutsRP(NULL), fNullCuts(0), fPIDResponse(0), fFlowEvent(0),fBayesianResponse(0), fCandidates(0),  fCandidateEtaPtCut(0), fCandidateMinEta(0), fCandidateMaxEta(0), fCandidateMinPt(0), fCandidateMaxPt(0), fNPtBins(18), fNdPhiBins(18), fCentrality(999), fVertex(999), fAOD(0), fPoolManager(0), fOutputList(0), fEventStats(0), fCentralityPass(0), fCentralityNoPass(0), fNOPID(0), fPIDk(0), fPIDp(0), fPtP(0), fPtN(0), fPtSpeciesA(0), fPtSpeciesB(0), fMultCorAfterCuts(0), fMultvsCentr(0), fCentralityMin(0), fCentralityMax(100), fkCentralityMethodA(0), fkCentralityMethodB(0), fCentralityCut2010(0), fCentralityCut2011(0), fPOICuts(0), fVertexRange(0), fPhi(0), fEta(0), fVZEROA(0), fVZEROC(0), fTPCM(0), fDeltaDipAngle(0), fDeltaDipPt(0), fApplyDeltaDipCut(0), fDCAAll(0), fDCAXYQA(0), fDCAZQA(0), fDCAPrim(0), fDCASecondaryWeak(0), fDCAMaterial(0), fSubEventDPhiv2(0), fAnalysisSummary(0)
{
   // Constructor
  for(Int_t i(0); i < 7; i++) fPIDConfig[i] = 0.;
  for(Int_t i(0); i < 5; i++) fDCAConfig[i] = 0.;
  for(Int_t i(0); i < 20; i++) {
      fVertexMixingBins[i] = 0;
      fCentralityMixingBins[i] = 0;
  }
  fMixingParameters[0] = 1000; fMixingParameters[1] = 50000; fMixingParameters[2] = 5;
  for(Int_t i(0); i < 18; i++) {
      for(Int_t j(0); j < 2; j++) {
          fV0Data[i][j] = 0;
          for(Int_t k(0); k < 15; k++) fPhiMinusPsiDataContainer[k][i][j] = 0x0;
          for(Int_t k(0); k < 15; k++) fPhiMinusPsiBackgroundContainer[k][i][j] = 0x0;
      }
      fResonanceSignal[i] = 0; fResonanceBackground[i] = 0; fPtSpectra[i] = 0; fPtBins[i] = 0.; fdPhiBins[i] = 0.;
  }
  for(Int_t i(0); i < 12; i++) fAddTaskMacroSummary[i] = 0.;
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
  DefineOutput(2, AliFlowEventSimple::Class());
}
//_____________________________________________________________________________
AliAnalysisTwoParticleResonanceFlowTask::~AliAnalysisTwoParticleResonanceFlowTask()
{
   // Destructor
   if (fNullCuts) delete fNullCuts;
   if (fOutputList) delete fOutputList;
   if (fCandidates) delete fCandidates;
   if (fFlowEvent) delete fFlowEvent;
   if (fBayesianResponse) delete fBayesianResponse;
   if (fEventMixing) delete fPoolManager;
}
//_____________________________________________________________________________
void AliAnalysisTwoParticleResonanceFlowTask::ForceExit(Int_t type, const char* message)
{
    // force exit
    if(type==0) cout << " --> Illegal options in AddTask <-- " << endl;
    if(type==1) cout << " --> Unknown error <-- " << endl;
    AliFatal(message);
}
//_____________________________________________________________________________
TH1F* AliAnalysisTwoParticleResonanceFlowTask::BookHistogram(const char* name)
{
   // Return a pointer to a TH1 with predefined binning
   TH1F *hist = new TH1F(name, Form("M_{INV} (%s)", name), 2*fMassBins, fMinMass, fMaxMass);
   hist->GetXaxis()->SetTitle("M_{INV} (GeV / c^{2})");
   hist->GetYaxis()->SetTitle("No. of pairs");
   hist->SetMarkerStyle(kFullCircle);
   hist->Sumw2();
   fOutputList->Add(hist);
   return hist;
}
//_____________________________________________________________________________
TH2F* AliAnalysisTwoParticleResonanceFlowTask::BookPIDHistogram(const char* name, Bool_t TPC)
{
   // Return a pointer to a TH2 with predefined binning
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
TH1F* AliAnalysisTwoParticleResonanceFlowTask::InitPtSpectraHistograms(Float_t nmin, Float_t nmax)
{
   // intialize p_t histograms for each p_t bin
   TH1F* hist = new TH1F(Form("%4.2f p_{t} %4.2f", nmin, nmax), Form("%f p_{t} %f", nmin, nmax), 2*fMassBins, nmin, nmax);
   hist->GetXaxis()->SetTitle("p_{T} GeV / c");
   fOutputList->Add(hist);
   return hist;
}
//_____________________________________________________________________________
TH1F* AliAnalysisTwoParticleResonanceFlowTask::BookPtHistogram(const char* name)
{
   // Return a pointer to a p_T spectrum histogram
   TH1F* ratio = new TH1F(name, name, 100, 0, 7);
   ratio->GetXaxis()->SetTitle("p_{T} ( GeV / c^{2} )");
   ratio->GetYaxis()->SetTitle("No. of events");
   ratio->Sumw2();
   fOutputList->Add(ratio);
   return ratio;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTwoParticleResonanceFlowTask::InitializeAnalysis()
{
    // set up the analysis
   fNullCuts = new AliFlowTrackCuts("null_cuts");
   fBayesianResponse = new AliFlowBayesianPID();
   Float_t t(0);
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
   fMassA = fMassA*fMassA;  //we always use mass squared in calculations, so doing it once in the init reduces cpu load
   fMassB = fMassB*fMassB;
   return kTRUE;
}
//_____________________________________________________________________________
void AliAnalysisTwoParticleResonanceFlowTask::UserCreateOutputObjects()
{
   // Create user defined output objects
   if(!InitializeAnalysis()) ForceExit(0, " > Couldn't initialize the analysis !!! <");
     // Create all output objects and store them to a list
   fOutputList = new TList();
   fOutputList->SetOwner(kTRUE);
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
       fPIDk = BookPIDHistogram("TPC signal, Species A", kTRUE);
       fPIDp = BookPIDHistogram("TPC signal, Species B", kTRUE);
   }
   for(Int_t ptbin(0); ptbin < fNPtBins; ptbin++) {
       if(!fPhiMinusPsiMethod) {
           fResonanceSignal[ptbin] = BookHistogram(Form("NP, %4.2f < p_{T} < %4.2f GeV", fPtBins[ptbin], fPtBins[ptbin+1]));
           fResonanceBackground[ptbin] = BookHistogram(Form("PP, %4.2f < p_{T} < %4.2f GeV", fPtBins[ptbin], fPtBins[ptbin+1]));
       }
       else if(fPhiMinusPsiMethod) {
           for(Int_t i(0); i<fNdPhiBins; i++) {// for all delta phi bins
                   fPhiMinusPsiDataContainer[i][ptbin][0] = BookHistogram(Form("%4.2f < p_{T} < %4.2f GeV, %4.2f < (#phi - #Psi_{A}) < %4.2f", fPtBins[ptbin], fPtBins[ptbin+1], fdPhiBins[i], fdPhiBins[i+1]));
                   fPhiMinusPsiDataContainer[i][ptbin][1] = BookHistogram(Form("%4.2f < p_{T} < %4.2f GeV, %4.2f < (#phi - #Psi_{C}) < %4.2f", fPtBins[ptbin], fPtBins[ptbin+1], fdPhiBins[i], fdPhiBins[i+1]));
                   fPhiMinusPsiBackgroundContainer[i][ptbin][0] = BookHistogram(Form("BG, %4.2f < p_{T} < %4.2f GeV, %4.2f < (#phi - #Psi_{A}) < %4.2f", fPtBins[ptbin], fPtBins[ptbin+1], fdPhiBins[i], fdPhiBins[i+1]));
                   fPhiMinusPsiBackgroundContainer[i][ptbin][1] = BookHistogram(Form("BG, %4.2f < p_{T} < %4.2f GeV, %4.2f < (#phi - #Psi_{C}) < %4.2f", fPtBins[ptbin], fPtBins[ptbin+1], fdPhiBins[i], fdPhiBins[i+1]));
           }
       }
   }
   if(fQA) {
       for(Int_t ptbin(0); ptbin < fNPtBins; ptbin++) fPtSpectra[ptbin] = InitPtSpectraHistograms(fPtBins[ptbin], fPtBins[ptbin+1]);
       fPtP = BookPtHistogram("i^{+}");
       fPtN = BookPtHistogram("i^{-}");
       fPtSpeciesA = BookPtHistogram("p_{T} spectrum Species A");
       fPtSpeciesB = BookPtHistogram("p_{T} spectrum Species B");
       fPhi = new TH1F("fPhi", "#phi distribution", 100, -.5, 7);
       fOutputList->Add(fPhi);
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
   if(fV0||fPhiMinusPsiMethod) {
       fSubEventDPhiv2 = new TProfile("fSubEventDPhiv2", "fSubEventDPhiv2", 3, 0, 3);
       fSubEventDPhiv2->GetXaxis()->SetBinLabel(1, "<cos(2(#Psi_{a} - #Psi_{b}))>");
       fSubEventDPhiv2->GetXaxis()->SetBinLabel(2, "<cos(2(#Psi_{a} - #Psi_{c}>))");
       fSubEventDPhiv2->GetXaxis()->SetBinLabel(3, "<cos(2(#Psi_{b} - #Psi_{c}>))");
       fOutputList->Add(fSubEventDPhiv2);
       const char* V0[] = {"V0A", "V0C"};
       for(Int_t ptbin(0); ptbin < fNPtBins; ptbin++)
           for(Int_t iV0(0); iV0 < 2; iV0++) {
                   fV0Data[ptbin][iV0] = new TProfile(Form("%s v2 %4.2f < p_{T} < %4.2f GeV", V0[iV0], fPtBins[ptbin], fPtBins[ptbin+1]), Form("%s v2 %4.2f < p_{T} < %4.2f GeV", V0[iV0], fPtBins[ptbin], fPtBins[ptbin+1]), fMassBins, fMinMass, fMaxMass);
                   fOutputList->Add(fV0Data[ptbin][iV0]);
           }
   }
   TString name = "fAnalysisSummary_";
   name+=this->GetName();
   fAnalysisSummary = new TH1F(name.Data(), name.Data(), 34, 0, 34);
   fAnalysisSummary->GetXaxis()->SetBinLabel(1, "fIsMC");
   fAnalysisSummary->GetXaxis()->SetBinLabel(2, "fEventMixing");
   fAnalysisSummary->GetXaxis()->SetBinLabel(3, "fAODAnalysis");
   if(fPIDConfig[6] > 0) fAnalysisSummary->GetXaxis()->SetBinLabel(4, "bayesian purity");
   else fAnalysisSummary->GetXaxis()->SetBinLabel(4, "TPC TOF PID");
   fAnalysisSummary->GetXaxis()->SetBinLabel(5, "dca mode");
   fAnalysisSummary->GetXaxis()->SetBinLabel(6, "fVertex");
   fAnalysisSummary->GetXaxis()->SetBinLabel(7, "fCentralityMin");
   fAnalysisSummary->GetXaxis()->SetBinLabel(8, "fCentralityMax");
   fAnalysisSummary->GetXaxis()->SetBinLabel(9, "centrality");
   fAnalysisSummary->GetXaxis()->SetBinLabel(10, fkCentralityMethodA);
   fAnalysisSummary->GetXaxis()->SetBinLabel(11, "fApplyDeltaDipCut");
   fAnalysisSummary->GetXaxis()->SetBinLabel(12, "POI filterbit");
   fAnalysisSummary->GetXaxis()->SetBinLabel(13, "RP filterbit");
   fAnalysisSummary->GetXaxis()->SetBinLabel(14, "PhiMinusPsiMethod");
   fAnalysisSummary->GetXaxis()->SetBinLabel(15, "SP");
   fAnalysisSummary->GetXaxis()->SetBinLabel(16, "SPSUB");
   fAnalysisSummary->GetXaxis()->SetBinLabel(17, "QC");
   fAnalysisSummary->GetXaxis()->SetBinLabel(18, "EP");
   fAnalysisSummary->GetXaxis()->SetBinLabel(19, "EP3sub");
   fAnalysisSummary->GetXaxis()->SetBinLabel(20, "V0_SP");
   fAnalysisSummary->GetXaxis()->SetBinLabel(21, "eta gap");
   fAnalysisSummary->GetXaxis()->SetBinLabel(22, "harm");
   fAnalysisSummary->GetXaxis()->SetBinLabel(23, "highPtMode");
   fAnalysisSummary->GetXaxis()->SetBinLabel(24, "deltaDipAngle");
   fAnalysisSummary->GetXaxis()->SetBinLabel(26, "SpeciesA");
   fAnalysisSummary->GetXaxis()->SetBinLabel(27, "chargeA");
   fAnalysisSummary->GetXaxis()->SetBinLabel(28, "fMinPtA");
   fAnalysisSummary->GetXaxis()->SetBinLabel(29, "fMaxPtA");
   fAnalysisSummary->GetXaxis()->SetBinLabel(30, "SpeciesB");
   fAnalysisSummary->GetXaxis()->SetBinLabel(31, "chargeB");
   fAnalysisSummary->GetXaxis()->SetBinLabel(32, "fMinPtB");
   fAnalysisSummary->GetXaxis()->SetBinLabel(33, "fMaxPtB");
   fAnalysisSummary->GetXaxis()->SetBinLabel(34, "iterator");
   fOutputList->Add(fAnalysisSummary);
   fAnalysisSummary->SetBinContent(1, (Float_t)fIsMC);
   fAnalysisSummary->SetBinContent(2, (Float_t)fEventMixing);
   fAnalysisSummary->SetBinContent(3, (Float_t)1);
   if(fPIDConfig[6] > 0) fAnalysisSummary->SetBinContent(4, (Float_t)fPIDConfig[6]);
   fAnalysisSummary->SetBinContent(5, (Float_t)fDCAConfig[0]);
   fAnalysisSummary->SetBinContent(6, (Float_t)fVertexRange);
   fAnalysisSummary->SetBinContent(7, (Float_t)fCentralityMin);
   fAnalysisSummary->SetBinContent(8, (Float_t)fCentralityMax);
   fAnalysisSummary->SetBinContent(11, (Float_t)fApplyDeltaDipCut);
   fAnalysisSummary->SetBinContent(12, (Float_t)fPOICuts->GetAODFilterBit());
   fAnalysisSummary->SetBinContent(13, (Float_t)fCutsRP->GetAODFilterBit());
   for(Int_t i(0); i<12;i++)  fAnalysisSummary->SetBinContent(14+i, fAddTaskMacroSummary[i]);
   fAnalysisSummary->SetBinContent(26, (Float_t)fSpeciesA);
   fAnalysisSummary->SetBinContent(27, (Float_t)fChargeA);
   fAnalysisSummary->SetBinContent(28, (Float_t)fMinPtA);
   fAnalysisSummary->SetBinContent(29, (Float_t)fMaxPtA);
   fAnalysisSummary->SetBinContent(30, (Float_t)fSpeciesB);
   fAnalysisSummary->SetBinContent(31, (Float_t)fChargeB);
   fAnalysisSummary->SetBinContent(32, (Float_t)fMinPtB);
   fAnalysisSummary->SetBinContent(33, (Float_t)fMaxPtB);
   fAnalysisSummary->SetBinContent(34, (Float_t)1);
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
AliEventPoolManager* AliAnalysisTwoParticleResonanceFlowTask::InitializeEventMixing()
{
   // initialize event mixing
  Int_t _c(0), _v(0);
  for(Int_t i(0); i < 19; i++) {
      if (fCentralityMixingBins[i+1] < fCentralityMixingBins[i]) { _c = i; break; }
      else _c = 19;
  }
  for(Int_t i(0); i < 19; i++) {
      if (fVertexMixingBins[i+1] < fVertexMixingBins[i]) { _v = i; break; }
      else _v = 19;
  }
  Double_t centralityBins[_c];
  Double_t vertexBins[_v];
  for(Int_t i(0); i < _c + 1; i++) centralityBins[i] = fCentralityMixingBins[i];
  for(Int_t i(0); i < _v + 1; i++) vertexBins[i] = fVertexMixingBins[i];
  return new AliEventPoolManager(fMixingParameters[0], fMixingParameters[1], _c, (Double_t*)centralityBins, _v, (Double_t*)vertexBins);
}
//_____________________________________________________________________________
template <typename T> Float_t AliAnalysisTwoParticleResonanceFlowTask::InvariantMass(const T* track1, const T* track2) const
{
   // Return the invariant mass of two tracks
   if ((!track2) || (!track1)) return 0.;
   Float_t pxs = TMath::Power((track1->Px() + track2->Px()), 2);
   Float_t pys = TMath::Power((track1->Py() + track2->Py()), 2);
   Float_t pzs = TMath::Power((track1->Pz() + track2->Pz()), 2);
   Float_t e1 = TMath::Sqrt(track1->P() * track1->P() + track1->Mass());
   Float_t e2 = TMath::Sqrt(track2->P() * track2->P() + track2->Mass());
   Float_t es = TMath::Power((e1 + e2), 2);
   if ((es - (pxs + pys + pzs)) < 0) return 0.;
   return TMath::Sqrt((es - (pxs + pys + pzs)));
}
//_____________________________________________________________________________
template <typename T> Float_t AliAnalysisTwoParticleResonanceFlowTask::DeltaDipAngle(const T* track1, const T* track2) const
{
   // Calculate the delta dip angle between two particles
   if (track1->P()*track2->P() == 0) return 999;
   return TMath::ACos(((track1->Pt() * track2->Pt()) + (track1->Pz() * track2->Pz())) / (track1->P() * track2->P()));
}
//_____________________________________________________________________________
template <typename T> Bool_t AliAnalysisTwoParticleResonanceFlowTask::CheckDeltaDipAngle(const T* track1, const T* track2) const
{
   // Check if pair passes delta dip angle cut within 0 < p_t < fDeltaDipPt
   if ((TMath::Abs(DeltaDipAngle(track1, track2)) < fDeltaDipAngle) && (PairPt(track1, track2) < fDeltaDipPt)) return kFALSE;
   return kTRUE;
}
//_____________________________________________________________________________
template <typename T> Bool_t AliAnalysisTwoParticleResonanceFlowTask::CheckCandidateEtaPtCut(const T* track1, const T* track2) const
{
   // Check if pair passes eta and pt cut
   if (fCandidateMinPt > PairPt(track1, track2) || fCandidateMaxPt < PairPt(track1, track2)) return kFALSE;
   TVector3 a(track1->Px(), track1->Py(), track1->Pz());
   TVector3 b(track2->Px(), track2->Py(), track2->Pz());
   TVector3 c = a + b;
   if (fCandidateMinEta > c.Eta() || fCandidateMaxEta < c.Eta()) return kFALSE;
   return kTRUE;
}
//_____________________________________________________________________________
template <typename T> Bool_t AliAnalysisTwoParticleResonanceFlowTask::EventCut(T* event)
{
   // Impose event cuts
   if (!event) return kFALSE;
   if (!CheckVertex(event)) return kFALSE;
   if (!CheckCentrality(event)) return kFALSE;
   if(fQA) PlotMultiplcities(event);
   return kTRUE;
}
//_____________________________________________________________________________
template <typename T> void AliAnalysisTwoParticleResonanceFlowTask::PlotMultiplcities(const T* event) const
{
   // QA multiplicity plots
   fVZEROA->Fill(event->GetVZEROData()->GetMTotV0A());
   fVZEROC->Fill(event->GetVZEROData()->GetMTotV0C());
   fTPCM->Fill(event->GetNumberOfTracks());
}
//_____________________________________________________________________________
template <typename T> Bool_t AliAnalysisTwoParticleResonanceFlowTask::CheckVertex(const T* event)
{
   // Check if event vertex is within given range
   if (!event->GetPrimaryVertex()) return 0x0;
   fVertex = event->GetPrimaryVertex()->GetZ();
   if (TMath::Abs(fVertex) > fVertexRange) return 0x0;
   return kTRUE;
}
//_____________________________________________________________________________
template <typename T> Bool_t AliAnalysisTwoParticleResonanceFlowTask::CheckCentrality(T* event)
{
   // Check if event is within the set centrality range. Falls back to V0 centrality determination if no method is set
   if(fDebug > 0) cout << " *** CheckCentrality() *** " << endl;
   if (!fkCentralityMethodA) AliFatal("No centrality method set! FATAL ERROR!");
   fCentrality = event->GetCentrality()->GetCentralityPercentile(fkCentralityMethodA);
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
       AliAODTrack* trackAOD = event->GetTrack(iTracks);
       if (!trackAOD) continue;
       if (!(trackAOD->TestFilterBit(1))) continue;
       if ((trackAOD->Pt() < .2) || (trackAOD->Pt() > 5.0) || (TMath::Abs(trackAOD->Eta()) > .8) || (trackAOD->GetTPCNcls() < 70)  || (trackAOD->GetDetPid()->GetTPCsignal() < 10.0) || (trackAOD->Chi2perNDF() < 0.2)) continue;
       multTPC++;
     }
     for(Int_t iTracks = 0; iTracks < nGoodTracks; iTracks++) { // fill global mult
       AliAODTrack* trackAOD = event->GetTrack(iTracks);
       if (!trackAOD) continue;
       if (!(trackAOD->TestFilterBit(16))) continue;
       if ((trackAOD->Pt() < .2) || (trackAOD->Pt() > 5.0) || (TMath::Abs(trackAOD->Eta()) > .8) || (trackAOD->GetTPCNcls() < 70) || (trackAOD->GetDetPid()->GetTPCsignal() < 10.0) || (trackAOD->Chi2perNDF() < 0.1)) continue;
       Double_t b[2] = {-99., -99.};
       Double_t bCov[3] = {-99., -99., -99.};
       if (!(trackAOD->PropagateToDCA(event->GetPrimaryVertex(), event->GetMagneticField(), 100., b, bCov))) continue;
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

   if(fCentralityCut2011)
     { // cut on outliers
     Float_t multTPC(0.); // tpc mult estimate
     Float_t multGlob(0.); // global multiplicity
     for(Int_t iTracks = 0; iTracks < nGoodTracks; iTracks++) { // fill tpc mult
       AliAODTrack* trackAOD = event->GetTrack(iTracks);
       if (!trackAOD) continue;
       if (!(trackAOD->TestFilterBit(1))) continue;
       if ((trackAOD->Pt() < .2) || (trackAOD->Pt() > 5.0) || (TMath::Abs(trackAOD->Eta()) > .8) || (trackAOD->GetTPCNcls() < 70)  || (trackAOD->GetDetPid()->GetTPCsignal() < 10.0) || (trackAOD->Chi2perNDF() < 0.2)) continue;
       multTPC++;
     }
     for(Int_t iTracks = 0; iTracks < nGoodTracks; iTracks++) { // fill global mult
       AliAODTrack* trackAOD = event->GetTrack(iTracks);
       if (!trackAOD) continue;
       if (!(trackAOD->TestFilterBit(16))) continue;
       if ((trackAOD->Pt() < .2) || (trackAOD->Pt() > 5.0) || (TMath::Abs(trackAOD->Eta()) > .8) || (trackAOD->GetTPCNcls() < 70) || (trackAOD->GetDetPid()->GetTPCsignal() < 10.0) || (trackAOD->Chi2perNDF() < 0.1)) continue;
       Double_t b[2] = {-99., -99.};
       Double_t bCov[3] = {-99., -99., -99.};
       if (!(trackAOD->PropagateToDCA(event->GetPrimaryVertex(), event->GetMagneticField(), 100., b, bCov))) continue;
       if ((TMath::Abs(b[0]) > 0.3) || (TMath::Abs(b[1]) > 0.3)) continue;
       multGlob++;
     } //track loop
     //     printf(" mult TPC %.2f, mult Glob %.2f \n", multTPC, multGlob);
     if(! (multTPC > (-36.73+1.48*multGlob) && multTPC < (62.87+1.78*multGlob))) return kFALSE;
     if(fQA) {  
       fMultCorAfterCuts->Fill(multGlob, multTPC);  
       fMultvsCentr->Fill(fCentrality, multTPC);
     }
   }
   
   fCentralityPass->Fill(fCentrality);
   return kTRUE;
}
//_____________________________________________________________________________
void AliAnalysisTwoParticleResonanceFlowTask::InitializeBayesianPID(AliAODEvent* event)
{
   // Initialize the Bayesian PID object for AOD
   fBayesianResponse->SetDetResponse(event, fCentrality);
}
//_____________________________________________________________________________
template <typename T> Bool_t AliAnalysisTwoParticleResonanceFlowTask::PassesTPCbayesianCut(T* track, Int_t species) const
{
   // Check if the particle passes the TPC TOF bayesian cut.
   fBayesianResponse->ComputeProb(track);
   if (!fBayesianResponse->GetCurrentMask(0)) return kFALSE; // return false if TPC has no response
   Int_t contamination(0); //type of particle we don't want in the sample
   if(species==3)  contamination = 2;  //we want to reject pions when selection kaons
   if(species==2)  contamination = 3; // we want to reject kaons when selecting pions
   Float_t *probabilities = fBayesianResponse->GetProb();
   return (probabilities[species] > fPIDConfig[6]  && (probabilities[contamination] < probabilities[species]));// 3 is kaon, 2 is pion
}
//_____________________________________________________________________________
Bool_t AliAnalysisTwoParticleResonanceFlowTask::PassesDCACut(AliAODTrack* track) const
{
    // check if track passes dca cut according to dca routine
    // setup the routine as follows:
    // fDCAConfig[0] < -1 no pt dependence
    // fDCAConfig[0] =  0 do nothing
    // fDCAConfig[0] >  1 pt dependent dca cut
    if(fIsMC) return kTRUE;
    if( (fDCAConfig[0] < 0.1) && (fDCAConfig[0] > -0.1) ) {
        if(fQA) {
            Double_t b[2] = { -99., -99.};
            Double_t bCov[3] = { -99., -99., -99.};
            track->PropagateToDCA(fAOD->GetPrimaryVertex(), fAOD->GetMagneticField(), 100., b, bCov);
            fDCAXYQA->Fill(b[0]);
            fDCAZQA->Fill(b[1]);
            fDCAPrim->Fill(track->Pt(), b[0]);
        }
        return kTRUE;
    }
    Double_t b[2] = { -99., -99.};
    Double_t bCov[3] = { -99., -99., -99.};
    track->PropagateToDCA(fAOD->GetPrimaryVertex(), fAOD->GetMagneticField(), 100., b, bCov);
    if((!fIsMC)&&fQA) fDCAMaterial->Fill(track->Pt(), b[0]);
    if( (fDCAConfig[0] < -.9) && ( (TMath::Abs(b[0]) > fDCAConfig[1]) || (TMath::Abs(b[1]) > fDCAConfig[2])) ) return kFALSE;
    if(fDCAConfig[0] > .9) {
        if(fDCAConfig[4] < TMath::Abs(b[1])) return kFALSE;
        Float_t denom = TMath::Power(track->Pt(), fDCAConfig[3]);
        if( denom < 0.0000001 ) return kFALSE; // avoid division by zero
        if( (fDCAConfig[1] + fDCAConfig[2] / denom) < TMath::Abs(b[0]) ) return kFALSE;
    }
    if(fQA) {
        fDCAXYQA->Fill(b[0]);
        fDCAZQA->Fill(b[1]);
        fDCAPrim->Fill(track->Pt(), b[0]);
    }
    return kTRUE;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTwoParticleResonanceFlowTask::DoOwnPID(AliAODTrack* track, Int_t species) const
{
    // do custom pid based on tpc and tof response
   Float_t nSigmasTPC        = (species==fSpeciesA) ? fPIDConfig[0] : fPIDConfig[3]; // number that is used when only tpc is available (5)
   Float_t nSigmasTPCwithTOF = (species==fSpeciesA) ? fPIDConfig[1] : fPIDConfig[4]; // number that is used for tpc when tof is available (5)
   Float_t nSigmasTOF        = (species==fSpeciesA) ? fPIDConfig[2] : fPIDConfig[5]; // number that is used for tof (if available) (3)
   if ((track->GetStatus() & AliESDtrack::kTOFout)&&(track->GetStatus() & AliESDtrack::kTIME)) { // check if the track is matched by tof signal
       if (TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track, (AliPID::EParticleType)species)) > nSigmasTPCwithTOF) return kFALSE; // reject with tpc
       if (track->P() < 1.5) return (TMath::Abs(fPIDResponse->NumberOfSigmasTOF(track, (AliPID::EParticleType)species)) <= nSigmasTOF);
       return (TMath::Abs(fPIDResponse->NumberOfSigmasTOF(track, (AliPID::EParticleType)species)) <= (nSigmasTOF-.333*nSigmasTOF));      
   }
   else { // switch to tpc pid if no tof hit is found
       if(track->Pt() < .35) return (TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track, (AliPID::EParticleType)species)) <= nSigmasTPC);
       else if(track->Pt() < .5) return (TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track, (AliPID::EParticleType)species)) < (nSigmasTPC-.4*nSigmasTPC));
       return (TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track, (AliPID::EParticleType)species)) <= (nSigmasTPC-.6*nSigmasTPC));
   }
   return kFALSE;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTwoParticleResonanceFlowTask::AcceptTrack(AliAODTrack* track, Int_t species) const
{
   if(fQA) TrackQA(track, species, kTRUE); // do qa before cuts
   if(species==fSpeciesA&&((track->Pt()<fMinPtA)||(track->Pt()>fMaxPtA))) return kFALSE; // least expensive check first
   if(species==fSpeciesB&&((track->Pt()<fMinPtB)||(track->Pt()>fMaxPtB))) return kFALSE;
   if(!PassesDCACut(track)) return kFALSE;
   if(fPIDConfig[6] < 0) {
       if(DoOwnPID(track, species)) {
           if(fQA) TrackQA(track, species, kFALSE);
           return kTRUE;
       }
       return kFALSE;
   }
   else {
       if (PassesTPCbayesianCut(track, species)) {
           if(fQA) TrackQA(track, species, kFALSE);
           return kTRUE;
       }
   }
   return kFALSE;
}
//_____________________________________________________________________________
template <typename T> Float_t AliAnalysisTwoParticleResonanceFlowTask::PairPt(const T* track1, const T* track2, Bool_t phi) const
{
   // return p_t of track pair, select phi to return the azimuthal angle instead
   TVector3 a(track1->Px(), track1->Py(), track1->Pz());
   TVector3 b(track2->Px(), track2->Py(), track2->Pz());
   TVector3 c = a + b;
   if(phi) return c.Phi();
   return c.Pt();
}
//_____________________________________________________________________________
template <typename T> Float_t AliAnalysisTwoParticleResonanceFlowTask::PtSelector(Int_t tracktype, const T* track1, const T* track2, Float_t mass) const
{
   // plot m_inv spectra of like- and unlike sign pairs
   Float_t pt = PairPt(track1, track2);
   if (tracktype == 0) { // signal histograms
       for(Int_t ptbin(0); ptbin < fNPtBins; ptbin++) {
           if ((fPtBins[ptbin] <= pt) && (fPtBins[ptbin+1] > pt )) {
               fResonanceSignal[ptbin]->Fill(mass);
               if(fQA) fPtSpectra[ptbin]->Fill(pt);
           }
       }
   }
   if (tracktype == 1) for(Int_t ptbin(0); ptbin < fNPtBins; ptbin++) if ((fPtBins[ptbin] <= pt) && (fPtBins[ptbin+1] > pt )) fResonanceBackground[ptbin]->Fill(mass);
   return pt;
}
//_____________________________________________________________________________
template <typename T> Bool_t AliAnalysisTwoParticleResonanceFlowTask::QualityCheck(T* track) const
{
   // check if track meets quality cuts
   if(!track) return kFALSE;
   return fPOICuts->IsSelected(track);
}
//_____________________________________________________________________________
void AliAnalysisTwoParticleResonanceFlowTask::TrackQA(AliAODTrack* track, Int_t species, Bool_t allChargedParticles) const
{
    // fill qa histo's
    if(allChargedParticles) { // track didn't pass identification yet
        fNOPID->Fill(track->P(), track->GetTPCsignal());
        if(track->Charge() > 0) {fEventStats->Fill(1); fPtP->Fill(track->Pt());}
        if(track->Charge() < 0) {fEventStats->Fill(2); fPtN->Fill(track->Pt());}
    }
    else if(!allChargedParticles) { // track has been identified
        if(species==fSpeciesA) {
            fPtSpeciesA->Fill(track->Pt());
            fPIDk->Fill(track->P(), track->GetTPCsignal());
        }
        else if(species==fSpeciesB) {
            fPtSpeciesB->Fill(track->Pt());
            fPIDp->Fill(track->P(), track->GetTPCsignal());
        }
        fPhi->Fill(track->Phi()); 
        fEta->Fill(track->Eta());
    }
}
//_____________________________________________________________________________
template <typename T> void AliAnalysisTwoParticleResonanceFlowTask::SetNullCuts(T* event)
{
   // Set null cuts
   fCutsRP->SetEvent(event, MCEvent());
   if(event) fNullCuts->SetParamType(AliFlowTrackCuts::kGlobal);
   fNullCuts->SetPtRange(+1, -1); // select nothing QUICK
   fNullCuts->SetEtaRange(+1, -1); // select nothing VZERO
   fNullCuts->SetEvent(event, MCEvent());
}
//_____________________________________________________________________________
void AliAnalysisTwoParticleResonanceFlowTask::PrepareFlowEvent(Int_t iMulti)
{
   // Prepare flow events
   fFlowEvent->ClearFast();
   fFlowEvent->Fill(fCutsRP, fNullCuts);
   fFlowEvent->SetReferenceMultiplicity(iMulti);
   fFlowEvent->DefineDeadZone(0, 0, 0, 0);
}
//_____________________________________________________________________________
void AliAnalysisTwoParticleResonanceFlowTask::PhiMinusPsiMethod(TObjArray* MixingCandidates)
{
    // extract v2 using the phi minus psi method
    if(!AliAnalysisTaskVnV0::IsPsiComputed()) { // AliAnalysisTaskVnV0 must be added to analysis que before this task !!!
        return;
    }
    // retrieve data
    Float_t abcPsi2[] = {AliAnalysisTaskVnV0::GetPsi2V0A(), AliAnalysisTaskVnV0::GetPsi2TPC(), AliAnalysisTaskVnV0::GetPsi2V0C()};
    fSubEventDPhiv2->Fill(0.5, TMath::Cos(2.*(abcPsi2[0]-abcPsi2[1]))); // vzeroa - tpc
    fSubEventDPhiv2->Fill(1.5, TMath::Cos(2.*(abcPsi2[0]-abcPsi2[2]))); // vzeroa - vzeroc
    fSubEventDPhiv2->Fill(2.5, TMath::Cos(2.*(abcPsi2[1]-abcPsi2[2]))); // tpc - vzeroc
    // if event plane stuff went ok, fill the histograms necessary for flow analysis with phi - psi method
    AliEventPool* pool = fPoolManager->GetEventPool(fCentrality, fVertex);
    if(!pool) AliFatal(Form("No pool found for centrality = %f, zVtx = %f, fatal error! ", fCentrality, fVertex));
    else if(pool->IsReady() || pool->NTracksInPool() > fMixingParameters[1] / 10 || pool->GetCurrentNEvents() >= fMixingParameters[2]) {
        Int_t nEvents = pool->GetCurrentNEvents();
        for (Int_t iEvent(0); iEvent < nEvents; iEvent++) {
            TObjArray* mixed_candidates = pool->GetEvent(iEvent);
            if(!mixed_candidates) continue; // this should NEVER happen
            Int_t bufferTracks = mixed_candidates->GetEntriesFast(); // buffered candidates
            Int_t candidates = MixingCandidates->GetEntriesFast(); // mixing candidates
            TObjArray* SpeciesA = new TObjArray();
            SpeciesA->SetOwner(kTRUE);
            TObjArray* SpeciesB = new TObjArray();
            SpeciesB->SetOwner(kTRUE);
            TObjArray* SpeciesAFromBuffer = new TObjArray();
            SpeciesAFromBuffer->SetOwner(kTRUE);
            TObjArray* SpeciesBFromBuffer = new TObjArray();
            SpeciesBFromBuffer->SetOwner(kTRUE);
            for (Int_t iTracks = 0; iTracks < candidates; iTracks++) { 
                AliResonanceFlowHelperTrack* track = (AliResonanceFlowHelperTrack*)MixingCandidates->At(iTracks);
                if (!track) continue;
                if (track->Charge() == fChargeA && track->Species() == fSpeciesA) SpeciesA->Add(track);
                else if (track->Charge() == fChargeB && track->Species() == fSpeciesB) SpeciesB->Add(track);
            }
            for (Int_t iTracks = 0; iTracks < bufferTracks; iTracks++) { 
                AliResonanceFlowHelperTrack* track = (AliResonanceFlowHelperTrack*)mixed_candidates->At(iTracks);
                if (!track) continue;
                if (track->Charge() == fChargeA && track->Species() == fSpeciesA ) SpeciesAFromBuffer->Add(track);
                else if (track->Charge() == fChargeB && track->Species() == fSpeciesB) SpeciesBFromBuffer->Add(track);
            }
            PhiMinusPsiMethodWriteData(kTRUE, SpeciesA, SpeciesB, abcPsi2); // write signal information to histograms
            PhiMinusPsiMethodWriteData(kFALSE, SpeciesA, SpeciesBFromBuffer, abcPsi2);
            PhiMinusPsiMethodWriteData(kFALSE, SpeciesAFromBuffer, SpeciesB, abcPsi2);
        } // end of mixed events loop
    } // end of checking to see whether pool is filled correctly
    pool->UpdatePool(MixingCandidates); // update pool with current mixing candidates (for next event)
}
//_____________________________________________________________________________
void AliAnalysisTwoParticleResonanceFlowTask::PhiMinusPsiMethodWriteData(Bool_t signal, TObjArray* SpeciesA, TObjArray* SpeciesB, Float_t* abcPsi2)
{
    // write data for phi minus psi method into correct  histograms
    for (Int_t i = 0; i < SpeciesA->GetEntriesFast(); i++) { // create pairs
        for(Int_t j = 0; j < SpeciesB->GetEntriesFast(); j++) { 
            AliResonanceFlowHelperTrack* trackA = (AliResonanceFlowHelperTrack*)(SpeciesA->At(i));
            AliResonanceFlowHelperTrack* trackB = (AliResonanceFlowHelperTrack*)(SpeciesB->At(j));
            if(!(trackA&&trackB)) continue; // shouldn't happen
            if(trackA->ID() == trackB->ID()) continue; // remove autocorrelations 
            if(fApplyDeltaDipCut && (!CheckDeltaDipAngle(trackA, trackB))) continue;
    	    if(fCandidateEtaPtCut && (!CheckCandidateEtaPtCut(trackA, trackB))) continue;
            Float_t minv = InvariantMass(trackA, trackB);
            Float_t pt = PtSelector(2, trackA, trackB, minv);
            TVector3 a(trackA->Px(), trackA->Py(), trackA->Pz());
	    TVector3 b(trackB->Px(), trackB->Py(), trackB->Pz());
	    TVector3 c = a + b;
	    Float_t phi = c.Phi();
//            if(phi < 0) phi += TMath::Pi();
            Float_t dPhi_a = phi - abcPsi2[0];
            Float_t dPhi_c = phi - abcPsi2[2];
//            if (dPhi_a < 0) dPhi_a += TMath::Pi();
//            if (dPhi_c < 0) dPhi_c += TMath::Pi();
            // now all necessary numbers are here, so we can put them into histograms
            for(Int_t k(0); k < fNdPhiBins; k++) { // dPhi bin loop
                if(dPhi_a >= fdPhiBins[k] && dPhi_a < fdPhiBins[k+1]) // see if track is in desired dPhi_a range
                for(Int_t ptbin(0); ptbin < fNPtBins; ptbin++) { // loop over pt bins
                    if(pt <= fPtBins[ptbin] || pt > fPtBins[ptbin+1]) continue;
                    signal ? fPhiMinusPsiDataContainer[k][ptbin][0]->Fill(minv) : fPhiMinusPsiBackgroundContainer[k][ptbin][0]->Fill(minv);
                }
                if(dPhi_c >= fdPhiBins[k] && dPhi_c < fdPhiBins[k+1]) // see if track is in desired dPhi_c range
                for(Int_t ptbin(0); ptbin < fNPtBins; ptbin++) { // loop over pt bins
                    if(pt <= fPtBins[ptbin] || pt > fPtBins[ptbin+1]) continue;
                    signal ? fPhiMinusPsiDataContainer[k][ptbin][1]->Fill(minv) : fPhiMinusPsiBackgroundContainer[k][ptbin][1]->Fill(minv);
                }
            }
        } // end of loop on i- tracks
    } // end of loop on i+ tracks
}
//_____________________________________________________________________________
void AliAnalysisTwoParticleResonanceFlowTask::VZEROSubEventAnalysis()
{
    // vzero event plane analysis using three subevents
    if(!AliAnalysisTaskVnV0::IsPsiComputed()) { // AliAnalysisTaskVnV0 must be added to analysis que before this task !!!
        return;
    }
    // retrieve data
    Float_t abcPsi2[] = {AliAnalysisTaskVnV0::GetPsi2V0A(), AliAnalysisTaskVnV0::GetPsi2TPC(), AliAnalysisTaskVnV0::GetPsi2V0C()};
    for(Int_t i(0); i < 3; i++) if(abcPsi2[i] == 0)  {
            return;
    }
    // save info for resolution calculations
    fSubEventDPhiv2->Fill(0.5, TMath::Cos(2.*(abcPsi2[0]-abcPsi2[1]))); // vzeroa - tpc
    fSubEventDPhiv2->Fill(1.5, TMath::Cos(2.*(abcPsi2[0]-abcPsi2[2]))); // vzeroa - vzeroc
    fSubEventDPhiv2->Fill(2.5, TMath::Cos(2.*(abcPsi2[1]-abcPsi2[2]))); // tpc - vzeroc
    Float_t minv[fMassBins+1];
    Float_t _dphi[fMassBins][fNPtBins][2]; //minv, pt, v0a-c
    Int_t occurence[fMassBins][fNPtBins]; //minv, pt
    Float_t _inc = (fMaxMass - fMinMass) / (Float_t)fMassBins;
    if(_inc==0) return; // avoid division by zero later on
    for(Int_t mb(0); mb < fMassBins+1; mb++) minv[mb] = fMinMass + mb * _inc;
    for(Int_t i(0); i < fMassBins; i++) for (Int_t j(0); j < fNPtBins; j++) for(Int_t k(0); k < 2; k ++) {
        _dphi[i][j][k] = 0;
        if(k==0) occurence[i][j] = 0;
    }
    for(Int_t iCand(0); iCand < fCandidates->GetEntriesFast(); iCand++) { // loop over unlike sign tracks
        AliFlowTrackSimple *track = dynamic_cast<AliFlowTrackSimple*>(fCandidates->At(iCand));
        if(!track) {
            continue;
        }
        for(Int_t mb(0); mb < fMassBins; mb++) { // loop over massbands
            if(track->Mass() <= minv[mb] || track->Mass() > minv[mb+1]) continue;
            for(Int_t ptbin(0); ptbin < fNPtBins; ptbin++) { // loop over pt bins
                if(track->Pt() <= fPtBins[ptbin] || track->Pt() > fPtBins[ptbin+1]) continue;
                _dphi[mb][ptbin][0]+=TMath::Cos(2.*(track->Phi() - abcPsi2[0]));
                _dphi[mb][ptbin][1]+=TMath::Cos(2.*(track->Phi() - abcPsi2[2]));
                occurence[mb][ptbin]+=1;
            }
        }
        for(Int_t mb(0); mb < fMassBins; mb++) // write vn values to tprofiles
            for(Int_t ptbin(0); ptbin < fNPtBins; ptbin++) {
                if(occurence[mb][ptbin]==0) continue;
                fV0Data[ptbin][0]->Fill(mb*_inc+fMinMass+0.5*_inc, _dphi[mb][ptbin][0]/(Float_t)occurence[mb][ptbin]);
                fV0Data[ptbin][1]->Fill(mb*_inc+fMinMass+0.5*_inc, _dphi[mb][ptbin][1]/(Float_t)occurence[mb][ptbin]);
            }
    }
}
//_____________________________________________________________________________
void AliAnalysisTwoParticleResonanceFlowTask::DoAnalysisOnTheFly(AliFlowEventSimple* event)
{
    // initialize the task for on the fly analysis and call the user exec
    if(fFlowEvent) delete fFlowEvent;           // clear out the old flow event 
    fFlowEvent = (AliFlowEvent*)event;          // cast the input event to its derived type
    UserExec("fly");                            // call the UserExec with flag 'on the fly'
}
//_____________________________________________________________________________
void AliAnalysisTwoParticleResonanceFlowTask::DoAnalysisOnTheFly(TDirectoryFile* outputFile)
{
    // write the anlaysis to an output file
    outputFile->Add(fOutputList);
    outputFile->Write(outputFile->GetName(), TObject::kSingleKey);
}
//_____________________________________________________________________________
void AliAnalysisTwoParticleResonanceFlowTask::DoAnalysisOnTheFly(TObjArray* MixingCandidates, TObjArray* SpeciesA, TObjArray* ocSpeciesA, TObjArray* SpeciesB, TObjArray* ocSpeciesB)
{
    // do the flow analysis on the fly. 
     Int_t iTracks = fFlowEvent->NumberOfTracks();
     fCandidates->SetLast(-1);                  // clean out candidate array
     Int_t  charge(-1), tID(0);                 // we'll use this guy to store charge
     if(fQA) fEventStats->Fill(0);              // event counter
     for (Int_t i = 0; i < iTracks; i++) {      // track loop. iterator i is used as unique track id (necessary later on to avoid auto-correlations)
        AliFlowTrackSimple* track = fFlowEvent->GetTrack(i);
        if(!track) continue;
        tID = track->GetID();                    // store ID
        (tID > 0) ? charge = 1 : charge = -1 ;          // get the charge of the track
        Double_t pt(track->Pt()), phi(track->Phi()), px(pt*TMath::Cos(phi)), py(pt*TMath::Sin(phi)), pz(track->Weight()), p(TMath::Sqrt(px*px+py*py+pz*pz)); // TODO ugly, but pz is stored as weight ...
        if (charge == fChargeA && TMath::Abs(tID)==TMath::Abs(fSpeciesA)) {     // store species a
            SpeciesA->Add(new AliResonanceFlowHelperTrack(track->Eta(), phi, p, px, py, pz, pt, charge, fMassA, i, fSpeciesA));
            if(fEventMixing) MixingCandidates->Add(new AliResonanceFlowHelperTrack(track->Eta(), phi, p, px, py, pz, pt, charge, fMassA, i, fSpeciesA));
        }
        if (charge == -1*fChargeA && TMath::Abs(tID)==TMath::Abs(fSpeciesA)) { // store opposite charge species a
           ocSpeciesA->Add(new AliResonanceFlowHelperTrack(track->Eta(), phi, p, px, py, pz, pt, charge, fMassA, i,fSpeciesA));                                             
           if(fEventMixing) MixingCandidates->Add(new AliResonanceFlowHelperTrack(track->Eta(), phi, p, px, py, pz, pt, charge, fMassA, i, fSpeciesA));
        }
        if (charge == fChargeB && TMath::Abs(tID)==TMath::Abs(fSpeciesB)) { // store species b
           SpeciesB->Add(new AliResonanceFlowHelperTrack(track->Eta(), phi, p, px, py, pz, pt, charge, fMassB, i, fSpeciesB));
           if(fEventMixing) MixingCandidates->Add(new AliResonanceFlowHelperTrack(track->Eta(), phi, p, px, py, pz, pt, charge, fMassB, i, fSpeciesB));
        }
        if (charge == -1*fChargeB && TMath::Abs(tID)==TMath::Abs(fSpeciesB)) { // store opposite charge species b
           ocSpeciesB->Add(new AliResonanceFlowHelperTrack(track->Eta(), phi, p, px, py, pz, pt, charge, fMassB, i, fSpeciesB));
           if(fEventMixing) MixingCandidates->Add(new AliResonanceFlowHelperTrack(track->Eta(), phi, p, px, py, pz, pt, charge, fMassB, i, fSpeciesB));
        }
        // at the end: convert the flow track to an 'actual' flow track with id and charge
        track->SetID(i);                                // set the unique id
        track->SetCharge(charge);                       // set charge
     }
}
//_____________________________________________________________________________
void AliAnalysisTwoParticleResonanceFlowTask::UserExec(Option_t * option)
{
  // UserExec: called for each event. Commented where necessary
   TObjArray* MixingCandidates = 0x0; // init null pointer for event mixing
   if(fEventMixing || fPhiMinusPsiMethod) {
       MixingCandidates = new TObjArray();
       MixingCandidates->SetOwner(kTRUE); // mixing candidates has ownership of objects in array
   }
   TObjArray* SpeciesA = new TObjArray(); // create arrays for the helper tracks
   SpeciesA->SetOwner(kTRUE);
   TObjArray* ocSpeciesA = new TObjArray(); // opposite charge particles
   ocSpeciesA->SetOwner(kTRUE);
   TObjArray* SpeciesB = new TObjArray();
   SpeciesB->SetOwner(kTRUE);
   TObjArray* ocSpeciesB = new TObjArray();
   ocSpeciesB->SetOwner(kTRUE);
   // check the options
   char chopt[128];
   strlcpy(chopt,option,128); 
   char *onTheFly;
   onTheFly = strstr(chopt,"fly");
   if(onTheFly) { // do the on the fly analysis
       printf(" > we're ready to fly ... ! \n"); 
       MixingCandidates = new TObjArray();
       MixingCandidates->SetOwner(kTRUE);
       DoAnalysisOnTheFly(MixingCandidates, SpeciesA, ocSpeciesA, SpeciesB, ocSpeciesB); 
   }
   if(!onTheFly) { // do analysis in the regular way (with the manager, etc)
      if (!fPIDResponse) { // kill the event if there isn't a pid response
         return;
      }
   
      fAOD = dynamic_cast<AliAODEvent*>(InputEvent()); // check for aod data type
      if (fAOD) {
         if (!EventCut(fAOD)) return; // check for event cuts
         InitializeBayesianPID(fAOD); // init event objects
         SetNullCuts(fAOD);
         Int_t iTracks = fAOD->GetNumberOfTracks();
         PrepareFlowEvent(iTracks); // does number of tracks correspond to the set filterbit ??? FIXME !!!
         fCandidates->SetLast(-1);
         if(fIsMC) IsMC(); // launch mc stuff FIXME
         if(fQA) fEventStats->Fill(0);
         for (Int_t i = 0; i < iTracks; i++) { // select analysis candidates
            AliAODTrack* track = fAOD->GetTrack(i);
            if (!QualityCheck(track)) continue; // reject poor quality tracks
            if (track->Charge() == fChargeA && AcceptTrack(track, fSpeciesA)) { // store species a
                SpeciesA->Add(new AliResonanceFlowHelperTrack(track->Eta(), track->Phi(), track->P(), track->Px(), track->Py(), track->Pz(), 
                                                              track->Pt(), track->Charge(), fMassA, track->GetID(), fSpeciesA));
   	     if(fEventMixing) MixingCandidates->Add(new AliResonanceFlowHelperTrack(track->Eta(), track->Phi(), track->P(),  track->Px(), 
                                                                                       track->Py(), track->Pz(), track->Pt(), track->Charge(), fMassA, 
                                                                                       track->GetID(), fSpeciesA));
   	 }
            if (track->Charge() == -1*fChargeA && AcceptTrack(track, fSpeciesA)) { // store opposite charge species a
                ocSpeciesA->Add(new AliResonanceFlowHelperTrack(track->Eta(), track->Phi(), track->P(), track->Px(), track->Py(), track->Pz(), 
                                                                track->Pt(), track->Charge(), fMassA, track->GetID(), fSpeciesA));
   	     if(fEventMixing) MixingCandidates->Add(new AliResonanceFlowHelperTrack(track->Eta(), track->Phi(), track->P(), track->Px(), 
                                                                                       track->Py(), track->Pz(), track->Pt(), track->Charge(), fMassA, 
                                                                                       track->GetID(), fSpeciesA));
    	 }
   	 if (track->Charge() == fChargeB && AcceptTrack(track, fSpeciesB)) { // store species b
                SpeciesB->Add(new AliResonanceFlowHelperTrack(track->Eta(), track->Phi(), track->P(), track->Px(), track->Py(), track->Pz(), 
                                                              track->Pt(), track->Charge(), fMassB, track->GetID(), fSpeciesB));
   	     if(fEventMixing) MixingCandidates->Add(new AliResonanceFlowHelperTrack(track->Eta(), track->Phi(), track->P(), track->Px(), 
                                                                                       track->Py(), track->Pz(), track->Pt(), track->Charge(), fMassB, 
                                                                                       track->GetID(), fSpeciesB));
   	 }
            if (track->Charge() == -1*fChargeB && AcceptTrack(track, fSpeciesB)) { // store opposite charge species b
                ocSpeciesB->Add(new AliResonanceFlowHelperTrack(track->Eta(), track->Phi(), track->P(), track->Px(), track->Py(), track->Pz(), 
                                                                track->Pt(), track->Charge(), fMassB, track->GetID(), fSpeciesB));
   	     if(fEventMixing) MixingCandidates->Add(new AliResonanceFlowHelperTrack(track->Eta(), track->Phi(), track->P(), track->Px(), 
                                                                                       track->Py(), track->Pz(), track->Pt(), track->Charge(), fMassB, 
                                                                                       track->GetID(), fSpeciesB));
   	     
            }
         }
      } // end the loop for event manager's aod's
   } // end of data loop
   // do the phi minus psi method if that's specified FIXME currenlty only supports event mixing
   if (fPhiMinusPsiMethod) { // for the phi minus psi method, no need for flow package etc 
       PhiMinusPsiMethod(MixingCandidates); // call the method
       PostData(1, fOutputList); // save the output data
       return; // return to retrieve next event
   }
   // start the resonance reconstruction, this fills the fCandidates array with flow tracks 
   ResonanceSignal(SpeciesA, SpeciesB); 
   if(fV0) VZEROSubEventAnalysis();
   for (int iCand = 0; iCand != fCandidates->GetEntriesFast(); ++iCand) {
      AliFlowCandidateTrack *cand = dynamic_cast<AliFlowCandidateTrack*>(fCandidates->At(iCand));
      if (!cand) continue;
      for (int iDau = 0; iDau != cand->GetNDaughters(); ++iDau) {
         for (int iRPs = 0; iRPs != fFlowEvent->NumberOfTracks(); ++iRPs) {
            AliFlowTrack *iRP = dynamic_cast<AliFlowTrack*>(fFlowEvent->GetTrack(iRPs));
            if (!iRP) continue;
            if (!iRP->InRPSelection()) continue;
            if (cand->GetIDDaughter(iDau) == iRP->GetID()) {
               iRP->SetForRPSelection(kFALSE);
               fFlowEvent->SetNumberOfRPs(fFlowEvent->GetNumberOfRPs() - 1);
            }
         }
      }
      cand->SetForPOISelection(kTRUE);
      fFlowEvent->InsertTrack(((AliFlowTrack*) cand));
   }
   if(!fEventMixing) { // do the combinatorial background
       ResonanceBackground(ocSpeciesA, SpeciesB); // mix opposite charge species A with species B
       if(fSpeciesA!=fSpeciesB) ResonanceBackground(SpeciesA, ocSpeciesB); // mix species A with opposite charge species B
   }
   else ReconstructionWithEventMixing(MixingCandidates);
   if(!onTheFly) {
       PostData(1, fOutputList);
       PostData(2, fFlowEvent);
   }
   delete SpeciesA;
   delete SpeciesB;
   delete ocSpeciesA;
   delete ocSpeciesB; // clean heap memory
}
//_____________________________________________________________________________
void AliAnalysisTwoParticleResonanceFlowTask::ResonanceSignal(TObjArray* SpeciesA, TObjArray* SpeciesB) const
{
    // fill signal histograms
    Int_t spA(SpeciesA->GetEntries()), spB(SpeciesB->GetEntries());
      for (Int_t i(0); i < spA; i++) { //track loop over species A
	for (Int_t j(0); j < spB; j++) { //track loop over species B
          AliResonanceFlowHelperTrack* trackA = (AliResonanceFlowHelperTrack*)SpeciesA->At(i);
          AliResonanceFlowHelperTrack* trackB = (AliResonanceFlowHelperTrack*)SpeciesB->At(j);
          if(!(trackA&&trackB)) continue; // shouldn't happen
          if(trackA->ID()==trackB->ID()) continue; // beware of autocorrelations
	  if (fApplyDeltaDipCut && (!CheckDeltaDipAngle(trackA, trackB))) continue;
	  if (fCandidateEtaPtCut && (!CheckCandidateEtaPtCut(trackA, trackB))) continue;
	  Float_t mass = InvariantMass(trackA, trackB);
	  Float_t pt = PtSelector(0, trackA, trackB, mass);
	  TVector3 a(trackA->Px(), trackA->Py(), trackA->Pz());
	  TVector3 b(trackB->Px(), trackB->Py(), trackB->Pz());
	  TVector3 c = a + b;
	  Float_t phi = c.Phi();
	  Float_t eta = c.Eta();
	  Int_t nIDs[2];
	  nIDs[0] = trackA->ID();
	  nIDs[1] = trackB->ID();
	  MakeTrack(mass, pt, phi, eta, 2, nIDs);
        }
    }
}
//_____________________________________________________________________________
void AliAnalysisTwoParticleResonanceFlowTask::ResonanceBackground(TObjArray* SpeciesA, TObjArray* SpeciesB, Bool_t checkAutoCorrelations) const
{
    // fill background histograms
     for (Int_t i(0); i < SpeciesA->GetEntriesFast(); i++) { //track loop over species A
	for (Int_t j(0); j < SpeciesB->GetEntriesFast(); j++) { //track loop over species B
          AliResonanceFlowHelperTrack* trackA = (AliResonanceFlowHelperTrack*)SpeciesA->At(i);
          AliResonanceFlowHelperTrack* trackB = (AliResonanceFlowHelperTrack*)SpeciesB->At(j);
          if(!(trackA&&trackB)) continue; // shouldn't happen
          if((trackA->ID() == trackB->ID()) && checkAutoCorrelations) continue; // remove autocorrelations 
	  if (fApplyDeltaDipCut && (!CheckDeltaDipAngle(trackA, trackB))) continue;
	  if (fCandidateEtaPtCut && (!CheckCandidateEtaPtCut(trackA, trackB))) continue;
          PtSelector(1, trackA, trackB, InvariantMass(trackA, trackB));
        }
    }
}
//_____________________________________________________________________________
void AliAnalysisTwoParticleResonanceFlowTask::ReconstructionWithEventMixing(TObjArray* MixingCandidates) const
{
    // perform reconstruction with event mixing
    AliEventPool* pool = fPoolManager->GetEventPool(fCentrality, fVertex);
    if(!pool) AliFatal(Form("No pool found for centrality = %f, zVtx = %f, fatal error!", fCentrality, fVertex));
    else if(pool->IsReady() || pool->NTracksInPool() > fMixingParameters[1] / 10 || pool->GetCurrentNEvents() >= fMixingParameters[2]) {
        Int_t nEvents = pool->GetCurrentNEvents();
        for (Int_t iEvent(0); iEvent < nEvents; iEvent++) {
            TObjArray* mixed_candidates = pool->GetEvent(iEvent);
            if(!mixed_candidates) continue; // this should NEVER happen
            Int_t bufferTracks = mixed_candidates->GetEntriesFast(); // buffered candidates
            Int_t candidates = MixingCandidates->GetEntriesFast(); // mixing candidates
            TObjArray* SpeciesA = new TObjArray();
            SpeciesA->SetOwner(kTRUE);
            TObjArray* SpeciesB = new TObjArray();
            SpeciesB->SetOwner(kTRUE);
            TObjArray* SpeciesAFromBuffer = new TObjArray();
            SpeciesAFromBuffer->SetOwner(kTRUE);
            TObjArray* SpeciesBFromBuffer = new TObjArray();
            SpeciesBFromBuffer->SetOwner(kTRUE);
            for (Int_t iTracks = 0; iTracks < candidates; iTracks++) { 
                AliResonanceFlowHelperTrack* track = (AliResonanceFlowHelperTrack*)MixingCandidates->At(iTracks);
                if (!track) continue;
                if (track->Charge() == fChargeA && track->Species() == fSpeciesA) SpeciesA->Add(track);
                else if (track->Charge() == fChargeB && track->Species() == fSpeciesB) SpeciesB->Add(track);
            }
            for (Int_t iTracks = 0; iTracks < bufferTracks; iTracks++) { 
                AliResonanceFlowHelperTrack* track = (AliResonanceFlowHelperTrack*)mixed_candidates->At(iTracks);
                if (!track) continue;
                if (track->Charge() == fChargeA && track->Species() == fSpeciesA ) SpeciesAFromBuffer->Add(track);
                else if (track->Charge() == fChargeB && track->Species() == fSpeciesB) SpeciesBFromBuffer->Add(track);
            }
            ResonanceBackground(SpeciesA, SpeciesBFromBuffer, kFALSE);
            ResonanceBackground(SpeciesB, SpeciesAFromBuffer, kFALSE); 
        } // end of mixed events loop
    } // end of checking to see whether pool is filled correctly
    pool->UpdatePool(MixingCandidates); // update pool with current mixing candidates (for next event)
}
//_____________________________________________________________________________
void AliAnalysisTwoParticleResonanceFlowTask::Terminate(Option_t *)
{
    // terminate
}
//______________________________________________________________________________
void  AliAnalysisTwoParticleResonanceFlowTask::MakeTrack(Float_t mass, Float_t pt, Float_t phi, Float_t eta, Int_t nDau, Int_t iID[]) const
{
   // Construct Flow Candidate Track from two selected candidates
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
   for (Int_t iDau = 0; iDau != nDau; ++iDau) sTrack->AddDaughter(iID[iDau]);
   sTrack->SetForPOISelection(kTRUE);
   sTrack->SetForRPSelection(kFALSE);
   if (overwrite) fCandidates->SetLast(fCandidates->GetLast() + 1);
   else fCandidates->AddLast(sTrack);
   return;
}
//_____________________________________________________________________________
void AliAnalysisTwoParticleResonanceFlowTask::IsMC()
{
    ForceExit(1,"No MC mode available at this stage of developement ... ");
}
//_____________________________________________________________________________
