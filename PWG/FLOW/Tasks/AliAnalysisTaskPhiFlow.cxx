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

// AliAnalysisTaskPhiFlow:
// author: Redmer Alexander Bertens (rbertens@nikhef.nl)
// analyis task for phi-meson reconstruction and determination of V2
// handles aod's and esd's transparently

#include "TChain.h"
#include "TTree.h"
#include "TH1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TObjArray.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliESDInputHandler.h"
#include "AliAODInputHandler.h"
#include "AliCentrality.h"
#include "AliVEvent.h"
#include "AliAnalysisTaskPhiFlow.h"
#include "AliFlowBayesianPID.h"
#include "AliStack.h"
#include "AliESDtrackCuts.h"
#include "AliMCEvent.h"
#include "TProfile.h"
#include "AliFlowCandidateTrack.h"
#include "AliFlowTrackCuts.h"
#include "AliFlowEventSimple.h"
#include "AliFlowCommonConstants.h"
#include "AliFlowEvent.h"
#include "TVector3.h"
#include "TRandom2.h"
#include "AliESDVZERO.h"
#include "AliAODVZERO.h"
#include "AliPID.h"
#include "AliPIDResponse.h"

class AliFlowTrackCuts;

ClassImp(AliAnalysisTaskPhiFlow)

AliAnalysisTaskPhiFlow::AliAnalysisTaskPhiFlow() : AliAnalysisTaskSE(),
   fDebug(0), fAODAnalysis(0), fMassBins(0), fMinMass(0), fMaxMass(0), fCutsRP(NULL), fNullCuts(0), fPIDResponse(0), fFlowEvent(0), fBayesianResponse(0), fCandidates(0), fOldTrackParam(0), fRequireTPCStandAlone(0), fStrictKaonCuts(0), fCandidateEtaPtCut(0), fCandidateMinEta(0), fCandidateMaxEta(0), fCandidateMinPt(0), fCandidateMaxPt(0), fCentrality(0), fESD(0), fAOD(0), fOutputList(0), fEventStats(0), fCentralityPass(0), fCentralityNoPass(0), fNOPID(0), fPIDk(0), fInvMNP03(0), fInvMPP03(0), fInvMNN03(0), fInvMNP36(0), fInvMPP36(0), fInvMNN36(0), fInvMNP69(0), fInvMPP69(0), fInvMNN69(0), fInvMNP912(0), fInvMPP912(0), fInvMNN912(0), fInvMNP1215(0), fInvMPP1215(0), fInvMNN1215(0), fInvMNP1518(0), fInvMPP1518(0), fInvMNN1518(0), fInvMNP1821(0), fInvMPP1821(0), fInvMNN1821(0), fInvMNP2124(0), fInvMPP2124(0), fInvMNN2124(0), fInvMNP2427(0), fInvMPP2427(0), fInvMNN2427(0), fInvMNP2730(0), fInvMPP2730(0), fInvMNN2730(0), fInvMNP3035(0), fInvMPP3035(0), fInvMNN3035(0), fInvMNP3540(0), fInvMPP3540(0), fInvMNN3540(0), fInvMNP4045(0), fInvMPP4045(0), fInvMNN4045(0), fInvMNP4550(0), fInvMPP4550(0), fInvMNN4550(0), fInvMNP5055(0), fInvMPP5055(0), fInvMNN5055(0), fInvMNP5560(0), fInvMPP5560(0), fInvMNN5560(0), fInvMNP6065(0), fInvMPP6065(0), fInvMNN6065(0), fInvMNP6570(0), fInvMPP6570(0), fInvMNN6570(0), fPtSpectra03(0), fPtSpectra36(0), fPtSpectra69(0), fPtSpectra912(0), fPtSpectra1215(0), fPtSpectra1518(0), fPtSpectra1821(0), fPtSpectra2124(0), fPtSpectra2427(0), fPtSpectra2730(0), fPtSpectra3035(0), fPtSpectra3540(0), fPtSpectra4045(0), fPtSpectra4550(0), fPtSpectra5055(0), fPtSpectra5560(0), fPtSpectra6065(0), fPtSpectra6570(0), fPtP(0), fPtN(0), fPtKP(0), fPtKN(0), fCentralityMin(0), fCentralityMax(100), fkCentralityMethod(0), fPOICuts(0), fVertexRange(0), fPhi(0), fPt(0), fEta(0), fVZEROA(0), fVZEROC(0), fTPCM(0), fDeltaDipAngle(0), fDeltaDipPt(0), fApplyDeltaDipCut(0), fDCAXY(0), fDCAZ(0), fDCA(0), fDCAXYQA(0), fDCAZQA(0)
{
   // Default constructor
   for(Int_t i = 0; i < 7; i++) fPIDConfig[i] = 0.;
}
//_____________________________________________________________________________
AliAnalysisTaskPhiFlow::AliAnalysisTaskPhiFlow(const char *name) : AliAnalysisTaskSE(name),
   fDebug(0), fAODAnalysis(0), fMassBins(0), fMinMass(0), fMaxMass(0), fCutsRP(NULL), fNullCuts(0), fPIDResponse(0), fFlowEvent(0), fBayesianResponse(0), fCandidates(0), fOldTrackParam(0), fRequireTPCStandAlone(0), fStrictKaonCuts(0), fCandidateEtaPtCut(0), fCandidateMinEta(0), fCandidateMaxEta(0), fCandidateMinPt(0), fCandidateMaxPt(0), fCentrality(0), fESD(0), fAOD(0), fOutputList(0), fEventStats(0), fCentralityPass(0), fCentralityNoPass(0), fNOPID(0), fPIDk(0), fInvMNP03(0), fInvMPP03(0), fInvMNN03(0), fInvMNP36(0), fInvMPP36(0), fInvMNN36(0), fInvMNP69(0), fInvMPP69(0), fInvMNN69(0), fInvMNP912(0), fInvMPP912(0), fInvMNN912(0), fInvMNP1215(0), fInvMPP1215(0), fInvMNN1215(0), fInvMNP1518(0), fInvMPP1518(0), fInvMNN1518(0), fInvMNP1821(0), fInvMPP1821(0), fInvMNN1821(0), fInvMNP2124(0), fInvMPP2124(0), fInvMNN2124(0), fInvMNP2427(0), fInvMPP2427(0), fInvMNN2427(0), fInvMNP2730(0), fInvMPP2730(0), fInvMNN2730(0), fInvMNP3035(0), fInvMPP3035(0), fInvMNN3035(0), fInvMNP3540(0), fInvMPP3540(0), fInvMNN3540(0), fInvMNP4045(0), fInvMPP4045(0), fInvMNN4045(0), fInvMNP4550(0), fInvMPP4550(0), fInvMNN4550(0), fInvMNP5055(0), fInvMPP5055(0), fInvMNN5055(0), fInvMNP5560(0), fInvMPP5560(0), fInvMNN5560(0), fInvMNP6065(0), fInvMPP6065(0), fInvMNN6065(0), fInvMNP6570(0), fInvMPP6570(0), fInvMNN6570(0), fPtSpectra03(0), fPtSpectra36(0), fPtSpectra69(0), fPtSpectra912(0), fPtSpectra1215(0), fPtSpectra1518(0), fPtSpectra1821(0), fPtSpectra2124(0), fPtSpectra2427(0), fPtSpectra2730(0), fPtSpectra3035(0), fPtSpectra3540(0), fPtSpectra4045(0), fPtSpectra4550(0), fPtSpectra5055(0), fPtSpectra5560(0), fPtSpectra6065(0), fPtSpectra6570(0), fPtP(0), fPtN(0), fPtKP(0), fPtKN(0), fCentralityMin(0), fCentralityMax(100), fkCentralityMethod(0), fPOICuts(0), fVertexRange(0), fPhi(0), fPt(0), fEta(0), fVZEROA(0), fVZEROC(0), fTPCM(0), fDeltaDipAngle(0), fDeltaDipPt(0), fApplyDeltaDipCut(0), fDCAXY(0), fDCAZ(0), fDCA(0), fDCAXYQA(0), fDCAZQA(0)
{
   // Constructor
   for(Int_t i = 0; i < 7; i++) fPIDConfig[i] = 0.;

   DefineInput(0, TChain::Class());
   DefineOutput(1, TList::Class());
   DefineOutput(2, AliFlowEventSimple::Class());
}
//_____________________________________________________________________________
AliAnalysisTaskPhiFlow::~AliAnalysisTaskPhiFlow()
{
   // Destructor
   if (fNullCuts) delete fNullCuts;
   if (fOutputList) delete fOutputList;
   if (fCandidates) delete fCandidates;
   if (fFlowEvent) delete fFlowEvent;
   if (fBayesianResponse) delete fBayesianResponse;
}
//_____________________________________________________________________________
TH1F* AliAnalysisTaskPhiFlow::BookHistogram(const char* name)
{
   // Return a pointer to a TH1 with predefined binning
   TH1F *hist = new TH1F(name, Form("M_{INV} (%s)", name), 60, .99, 1.092);
   hist->GetXaxis()->SetTitle("M_{INV} (GeV / c^{2})");
   hist->GetYaxis()->SetTitle("No. of pairs");
   hist->SetMarkerStyle(kFullCircle);
   hist->Sumw2();
   fOutputList->Add(hist);
   return hist;
}
//_____________________________________________________________________________
TH2F* AliAnalysisTaskPhiFlow::BookPIDHistogram(const char* name)
{
   // Return a pointer to a TH2 with predefined binning
   TH2F *hist = new TH2F(name, Form("PID (%s)", name), 100, 0, 5, 100, 0, 1000);
   hist->GetXaxis()->SetTitle("P (GeV / c)");
   hist->GetYaxis()->SetTitle("dE/dx (a.u.)");
   fOutputList->Add(hist);
   return hist;
}
//_____________________________________________________________________________
TH1F* AliAnalysisTaskPhiFlow::InitPtSpectraHistograms(Int_t i)
{
   // intialize p_t histograms for each p_t bin
   Double_t nmin(0);
   Double_t nmax(0);
   (i < 11) ? nmin = 0.3 * (i - 1) : nmin = 0.5 * (i - 11) + 3;
   (i < 11) ? nmax = 0.3 * (i - 1) + 0.3 : nmax = 0.5 * (i - 11) + 3.5;
   TH1F* hist = new TH1F(Form("%f p_{t} %f", nmin, nmax), Form("%f p_{t} %f", nmin, nmax), 60, nmin, nmax);
   hist->GetXaxis()->SetTitle("p_{T} GeV / c");
   fOutputList->Add(hist);
   return hist;
}
//_____________________________________________________________________________
TH1F* AliAnalysisTaskPhiFlow::BookPtHistogram(const char* name)
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
void AliAnalysisTaskPhiFlow::AddPhiIdentificationOutputObjects()
{
   // Add Phi Identification Output Objects
   fEventStats = new TH1F("fHistStats", "Event Statistics", 18, -.25, 4.25);
   fEventStats->GetXaxis()->SetTitle("No. of events");
   fEventStats->GetYaxis()->SetTitle("Statistics");
   fOutputList->Add(fEventStats);
   fCentralityPass = new TH1F("fCentralityPass", "Centrality Pass", 101, -1, 100);
   fCentralityNoPass = new TH1F("fCentralityNoPass", "Centrality No Pass", 101, -1, 100);
   fOutputList->Add(fCentralityPass);
   fOutputList->Add(fCentralityNoPass);

   fNOPID = BookPIDHistogram("TPC signal, all particles");
   fPIDk = BookPIDHistogram("TPC signal, kaons");
   fInvMNP03 = BookHistogram("NP, 0 < p_{T} < 0.3 GeV");
   fInvMPP03 = BookHistogram("PP, 0 < p_{T} < 0.3 GeV");
   fInvMNN03 = BookHistogram("NN, 0 < p_{T} < 0.3 GeV");
   fInvMNP36 = BookHistogram("NP, 0.3 < p_{T} < 0.6 GeV");
   fInvMPP36 = BookHistogram("PP, 0.3 < p_{T} < 0.6 GeV");
   fInvMNN36 = BookHistogram("NN, 0.3 < p_{T} < 0.6 GeV");
   fInvMNP69 = BookHistogram("NP, 0.6 < p_{T} < 0.9 GeV");
   fInvMPP69 = BookHistogram("PP, 0.6 < p_{T} < 0.9 GeV");
   fInvMNN69 = BookHistogram("NN, 0.6 < p_{T} < 0.9 GeV");
   fInvMNP912 = BookHistogram("NP, 0.9 < p_{T} < 1.2 GeV");
   fInvMPP912 = BookHistogram("PP, 0.9 < p_{T} < 1.2 GeV");
   fInvMNN912 = BookHistogram("NN, 0.9 < p_{T} < 1.2 GeV");
   fInvMNP1215 = BookHistogram("NP, 1.2 < p_{T} < 1.5 GeV");
   fInvMPP1215 = BookHistogram("PP, 1.2 < p_{T} < 1.5 GeV");
   fInvMNN1215 = BookHistogram("NN, 1.2 < p_{T} < 1.5 GeV");
   fInvMNP1518 = BookHistogram("NP, 1.5 < p_{T} < 1.8 GeV");
   fInvMPP1518 = BookHistogram("PP, 1.5 < p_{T} < 1.8 GeV");
   fInvMNN1518 = BookHistogram("NN, 1.5 < p_{T} < 1.8 GeV");
   fInvMNP1821 = BookHistogram("NP, 1.8 < p_{T} < 2.1 GeV");
   fInvMPP1821 = BookHistogram("PP, 1.8 < p_{T} < 2.1 GeV");
   fInvMNN1821 = BookHistogram("NN, 1.8 < p_{T} < 2.1 GeV");
   fInvMNP2124 = BookHistogram("NP, 2.1 < p_{T} < 2.4 GeV");
   fInvMPP2124 = BookHistogram("PP, 2.1 < p_{T} < 2.4 GeV");
   fInvMNN2124 = BookHistogram("NN, 2.1 < p_{T} < 2.4 GeV");
   fInvMNP2427 = BookHistogram("NP, 2.4 < p_{T} < 2.7 GeV");
   fInvMPP2427 = BookHistogram("PP, 2.4 < p_{T} < 2.7 GeV");
   fInvMNN2427 = BookHistogram("NN, 2.4 < p_{T} < 2.7 GeV");
   fInvMNP2730 = BookHistogram("NP, 2.7 < p_{T} < 3.0 GeV");
   fInvMPP2730 = BookHistogram("PP, 2.7 < p_{T} < 3.0 GeV");
   fInvMNN2730 = BookHistogram("NN, 2.7 < p_{T} < 3.0 GeV");
   fInvMNP3035 = BookHistogram("NP, 3.0 < p_{T} < 3.5 GeV");
   fInvMPP3035 = BookHistogram("PP, 3.0 < p_{T} < 3.5 GeV");
   fInvMNN3035 = BookHistogram("NN, 3.0 < p_{T} < 3.5 GeV");
   fInvMNP3540 = BookHistogram("NP, 3.5 < p_{T} < 4.0 GeV");
   fInvMPP3540 = BookHistogram("PP, 3.5 < p_{T} < 4.0 GeV");
   fInvMNN3540 = BookHistogram("NN, 3.5 < p_{T} < 4.0 GeV");
   fInvMNP4045 = BookHistogram("NP, 4.0 < p_{T} < 4.5 GeV");
   fInvMPP4045 = BookHistogram("PP, 4.0 < p_{T} < 4.5 GeV");
   fInvMNN4045 = BookHistogram("NN, 4.0 < p_{T} < 4.5 GeV");
   fInvMNP4550 = BookHistogram("NP, 4.5 < p_{T} < 5.0 GeV");
   fInvMPP4550 = BookHistogram("PP, 4.5 < p_{T} < 5.0 GeV");
   fInvMNN4550 = BookHistogram("NN, 4.5 < p_{T} < 5.0 GeV");
   fInvMNP5055 = BookHistogram("NP, 5.0 < p_{T} < 5.5 GeV");
   fInvMPP5055 = BookHistogram("PP, 5.0 < p_{T} < 5.5 GeV");
   fInvMNN5055 = BookHistogram("NN, 5.0 < p_{T} < 5.5 GeV");
   fInvMNP5560 = BookHistogram("NP, 5.5 < p_{T} < 6.0 GeV");
   fInvMPP5560 = BookHistogram("PP, 5.5 < p_{T} < 6.0 GeV");
   fInvMNN5560 = BookHistogram("NN, 5.5 < p_{T} < 6.0 GeV");
   fInvMNP6065 = BookHistogram("NP, 6.0 < p_{T} < 6.5 GeV");
   fInvMPP6065 = BookHistogram("PP, 6.0 < p_{T} < 6.5 GeV");
   fInvMNN6065 = BookHistogram("NN, 6.0 < p_{T} < 6.5 GeV");
   fInvMNP6570 = BookHistogram("NP, 6.5 < p_{T} < 7.0 GeV");
   fInvMPP6570 = BookHistogram("PP, 6.5 < p_{T} < 7.0 GeV");
   fInvMNN6570 = BookHistogram("NN, 6.5 < p_{T} < 7.0 GeV");

   fPtSpectra03 = InitPtSpectraHistograms(1);
   fPtSpectra36 = InitPtSpectraHistograms(2);
   fPtSpectra69 = InitPtSpectraHistograms(3);
   fPtSpectra912 = InitPtSpectraHistograms(4);
   fPtSpectra1215 = InitPtSpectraHistograms(5);
   fPtSpectra1518 = InitPtSpectraHistograms(6);
   fPtSpectra1821 = InitPtSpectraHistograms(7);
   fPtSpectra2124 = InitPtSpectraHistograms(8);
   fPtSpectra2427 = InitPtSpectraHistograms(9);
   fPtSpectra2730 = InitPtSpectraHistograms(10);
   fPtSpectra3035 = InitPtSpectraHistograms(11);
   fPtSpectra3540 = InitPtSpectraHistograms(12);
   fPtSpectra4045 = InitPtSpectraHistograms(13);
   fPtSpectra4550 = InitPtSpectraHistograms(14);
   fPtSpectra5055 = InitPtSpectraHistograms(15);
   fPtSpectra5560 = InitPtSpectraHistograms(16);
   fPtSpectra6065 = InitPtSpectraHistograms(17);
   fPtSpectra6570 = InitPtSpectraHistograms(18);

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
}
//_____________________________________________________________________________
void AliAnalysisTaskPhiFlow::UserCreateOutputObjects()
{
   // Create user defined output objects
   fNullCuts = new AliFlowTrackCuts("null_cuts");
   fBayesianResponse = new AliFlowBayesianPID();
   Double_t t(0);
   for(Int_t i = 0; i < 7; i++) t+=TMath::Abs(fPIDConfig[i]);
   if(t < 0.1) AliFatal("No valid PID procedure recognized -- terminating analysis!!!");

   AliFlowCommonConstants* cc = AliFlowCommonConstants::GetMaster();
   cc->SetNbinsMult(10000);
   cc->SetMultMin(0);
   cc->SetMultMax(10000);

   cc->SetNbinsPt(100);
   cc->SetPtMin(0);
   cc->SetPtMax(10);

   cc->SetNbinsPhi(180);
   cc->SetPhiMin(0.0);
   cc->SetPhiMax(TMath::TwoPi());

   cc->SetNbinsEta(200);
   cc->SetEtaMin(-5.0);
   cc->SetEtaMax(+5.0);

   cc->SetNbinsQ(500);
   cc->SetQMin(0.0);
   cc->SetQMax(3.0);

   cc->SetNbinsMass(fMassBins);
   cc->SetMassMin(fMinMass);
   cc->SetMassMax(fMaxMass);

   // setup initial state of PID response object
   if (!fOldTrackParam) fBayesianResponse->SetNewTrackParam();

   AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
   if (man)
   {
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
}
//_____________________________________________________________________________
template <typename T> Double_t AliAnalysisTaskPhiFlow::InvariantMass(const T* track1, const T* track2) const
{
   // Return the invariant mass of two tracks, assuming both tracks are kaons
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
template <typename T> Double_t AliAnalysisTaskPhiFlow::DeltaDipAngle(const T* track1, const T* track2) const
{
   // Calculate the delta dip angle between two particles
   if (track1->P()*track2->P() == 0) return 999;
   return TMath::ACos(((track1->Pt() * track2->Pt()) + (track1->Pz() * track2->Pz())) / (track1->P() * track2->P()));
}
//_____________________________________________________________________________
template <typename T> Bool_t AliAnalysisTaskPhiFlow::CheckDeltaDipAngle(const T* track1, const T* track2) const
{
   // Check if pair passes delta dip angle cut within 0 < p_t < fDeltaDipPt
   if ((TMath::Abs(DeltaDipAngle(track1, track2)) < fDeltaDipAngle) && (PhiPt(track1, track2) < fDeltaDipPt)) return kFALSE;
   return kTRUE;
}
//_____________________________________________________________________________
template <typename T> Bool_t AliAnalysisTaskPhiFlow::CheckCandidateEtaPtCut(const T* track1, const T* track2) const
{
   // Check if pair passes eta and pt cut
   if (fCandidateMinPt > PhiPt(track1, track2) || fCandidateMaxPt < PhiPt(track1, track2)) return kFALSE;
   TVector3 a(track1->Px(), track1->Py(), track1->Pz());
   TVector3 b(track2->Px(), track2->Py(), track2->Pz());
   TVector3 c = a + b;
   if (fCandidateMinEta > c.Eta() || fCandidateMaxEta < c.Eta()) return kFALSE;
   return kTRUE;
}
//_____________________________________________________________________________
void AliAnalysisTaskPhiFlow::SetCentralityParameters(Double_t CentralityMin, Double_t CentralityMax, const char* CentralityMethod)
{
   // Set a centrality range ]min, max] and define the method to use for centrality selection
   fCentralityMin = CentralityMin;
   fCentralityMax = CentralityMax;
   fkCentralityMethod = CentralityMethod;
}
//_____________________________________________________________________________
template <typename T> Bool_t AliAnalysisTaskPhiFlow::EventCut(T* event)
{
   // Impose event cuts
   if (!event) return kFALSE;
   if (!CheckVertex(event)) return kFALSE;
   if (!CheckCentrality(event)) return kFALSE;
   PlotVZeroMultiplcities(event);
   return kTRUE;
}
//_____________________________________________________________________________
template <typename T> void AliAnalysisTaskPhiFlow::PlotVZeroMultiplcities(const T* event) const
{
   // QA multiplicity plots
   fVZEROA->Fill(event->GetVZEROData()->GetMTotV0A());
   fVZEROC->Fill(event->GetVZEROData()->GetMTotV0C());
}
//_____________________________________________________________________________
template <typename T> Bool_t AliAnalysisTaskPhiFlow::CheckVertex(const T* event) const
{
   // Check if event vertex is within given range
   if (!event->GetPrimaryVertex()) return 0x0;
   if (TMath::Abs((event->GetPrimaryVertex())->GetZ()) > fVertexRange) return 0x0;
   return kTRUE;
}
//_____________________________________________________________________________
template <typename T> Bool_t AliAnalysisTaskPhiFlow::CheckCentrality(T* event)
{
   // Check if event is within the set centrality range. Falls back to V0 centrality determination if no method is set
   if (!fkCentralityMethod) AliFatal("No centrality method set! FATAL ERROR!");
   fCentrality = event->GetCentrality()->GetCentralityPercentile(fkCentralityMethod);
   if ((fCentrality <= fCentralityMin) || (fCentrality > fCentralityMax))
   {
      fCentralityNoPass->Fill(fCentrality) ;
      return kFALSE;
   }
   fCentralityPass->Fill(fCentrality);
   return kTRUE;
}
//_____________________________________________________________________________
void AliAnalysisTaskPhiFlow::InitializeBayesianPID(AliESDEvent* event)
{
   // Initialize the Bayesian PID object for ESD analysis
   fBayesianResponse->SetDetResponse(event, fCentrality, AliESDpid::kTOF_T0, kTRUE);
   if (fDebug) cout << " --> Initialized Bayesian ESD PID object " << endl;
}
//_____________________________________________________________________________
void AliAnalysisTaskPhiFlow::InitializeBayesianPID(AliAODEvent* event)
{
   // Initialize the Bayesian PID object for AOD
   fBayesianResponse->SetDetResponse(event, fCentrality);
   if (fDebug) cout << " --> Initialized Bayesian AOD PID object " << endl;
}
//_____________________________________________________________________________
template <typename T> Bool_t AliAnalysisTaskPhiFlow::PassesTPCbayesianCut(T* track) const
{
   // Check if the particle passes the TPC TOF bayesian cut.
   if ((!fStrictKaonCuts) && (!PassesStrictKaonCuts(track))) return kFALSE;
   fBayesianResponse->ComputeProb(track);
   if (!fBayesianResponse->GetCurrentMask(0)) return kFALSE; // return false if TPC has no response
   Float_t *probabilities = fBayesianResponse->GetProb();
   if (probabilities[3] > fPIDConfig[6])
   {
      fPhi->Fill(track->Phi());
      fPt->Fill(track->Pt());
      fEta->Fill(track->Eta());
      return kTRUE;
   }
   return kFALSE;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskPhiFlow::PassesStrictKaonCuts(AliESDtrack* track) const
{
   // propagate dca from TPC
   Double_t b[2] = { -99., -99.};
   Double_t bCov[3] = { -99., -99., -99.};
   if (!track->PropagateToDCA(fESD->GetPrimaryVertex(), fESD->GetMagneticField(), 100., b, bCov)) return kFALSE;
   if ((TMath::Abs(b[0]) > 3.0) || (TMath::Abs(b[1]) > 2.4)) return kFALSE;
   return kTRUE;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskPhiFlow::PassesStrictKaonCuts(AliAODTrack* track) const
{
   // propagate dca from TPC
   Double_t b[2] = { -99., -99.};
   Double_t bCov[3] = { -99., -99., -99.};
   if (!track->PropagateToDCA(fAOD->GetPrimaryVertex(), fAOD->GetMagneticField(), 100., b, bCov)) return kFALSE;
   if ((TMath::Abs(b[0]) > 3.0) || (TMath::Abs(b[1]) > 2.4)) return kFALSE;
   return kTRUE;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskPhiFlow::IsKaon(AliESDtrack* track) const
{
   // Check if particle is a kaon according to method set in steering macro
   if (fRequireTPCStandAlone && (track->GetStatus()&AliESDtrack::kTPCin) == 0) return kFALSE;
   if (fPIDConfig[1]!=0 || fPIDConfig[4]!=0) AliFatal("TPC || ITS PID not available in ESD anlaysis -- terminating analysis !!!");
   if (PassesTPCbayesianCut(track))
   {
     fPIDk->Fill(track->P(), track->GetTPCsignal());
     return kTRUE;
   }
   return kFALSE;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskPhiFlow::IsKaon(AliAODTrack* track) const
{
   // Kaon identification routine, based on multiple detectors and approaches
   if (fRequireTPCStandAlone && (!track->TestFilterBit(1))) return kFALSE;
   fNOPID->Fill(track->P(), track->GetTPCsignal());
   if(track->Pt() < fPIDConfig[1])
   {
       if(fDebug) cout << " ITS received track with p_t " << track->Pt() << endl;
       // if tpc control is disabled, pure its pid
       if(fPIDConfig[2] < 0.)
       {
           if (TMath::Abs(fPIDResponse->NumberOfSigmasITS(track, AliPID::kKaon)) < fPIDConfig[0]) return kTRUE;
           return kFALSE;
       }
       // else, switch to ITS pid with TPC rejection of protons and pions
       if (TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton)) < fPIDConfig[3]) return kFALSE;
       else if (TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion)) < fPIDConfig[3]) return kFALSE;
       else if (TMath::Abs(fPIDResponse->NumberOfSigmasITS(track, AliPID::kKaon)) < fPIDConfig[0])
       {
           fPIDk->Fill(track->P(), track->GetTPCsignal());
           return kTRUE;
       }
       return kFALSE;
   }
   if ((track->Pt() > fPIDConfig[1]) && (track->Pt() < fPIDConfig[4]))
   {
       if(fDebug) cout << " TPC received track with p_t " << track->Pt() << endl;
       // if its control is disabled, pure tpc pid
       if(fPIDConfig[5] < 0.)
       {
           if (TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track, AliPID::kKaon)) < fPIDConfig[3]) return kTRUE;
           return kFALSE;
       }
       // else, switch to TPC pid with ITS rejection of protons and pions
       if (TMath::Abs(fPIDResponse->NumberOfSigmasITS(track, AliPID::kProton)) < fPIDConfig[0]) return kFALSE;
       else if (TMath::Abs(fPIDResponse->NumberOfSigmasITS(track, AliPID::kPion)) < fPIDConfig[0]) return kFALSE;
       else if (TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track, AliPID::kKaon)) < fPIDConfig[3])
       {
           fPIDk->Fill(track->P(), track->GetTPCsignal());
           return kTRUE;
       }
       return kFALSE;
   }
   if(fDebug) cout << " Bayesian method received track with p_t " << track->Pt() << endl;
   // switch to bayesian PID
   if (PassesTPCbayesianCut(track))
   {
       fPIDk->Fill(track->P(), track->GetTPCsignal());
       return kTRUE;
   }
   return kFALSE;
}
//_____________________________________________________________________________
template <typename T> Double_t AliAnalysisTaskPhiFlow::PhiPt(const T* track1, const T* track2) const
{
   // Calculate transverse momentum (p_T) of two tracks
   TVector3 a(track1->Px(), track1->Py(), track1->Pz());
   TVector3 b(track2->Px(), track2->Py(), track2->Pz());
   TVector3 c = a + b;
   return c.Pt();
}
//_____________________________________________________________________________
template <typename T> void AliAnalysisTaskPhiFlow::PtSelector(Int_t tracktype, const T* track1, const T* track2) const
{
   // Request transverse momentum (p_T), then fill invariant mass histograms as a function of p_T
   Double_t pt = PhiPt(track1, track2);
   if (tracktype == 0)
   {
      if ((0.0 <= pt) && (0.3 > pt))
      {
         fInvMNP03->Fill(InvariantMass(track1, track2));
         fPtSpectra03->Fill(pt);
      }
      if ((0.3 <= pt) && (0.6 > pt))
      {
         fInvMNP36->Fill(InvariantMass(track1, track2));
         fPtSpectra36->Fill(pt);
      }
      if ((0.6 <= pt) && (0.9 > pt))
      {
         fInvMNP69->Fill(InvariantMass(track1, track2));
         fPtSpectra69->Fill(pt);
      }
      if ((0.9 <= pt) && (1.2 > pt))
      {
         fInvMNP912->Fill(InvariantMass(track1, track2));
         fPtSpectra912->Fill(pt);
      }
      if ((1.2 <= pt) && (1.5 > pt))
      {
         fInvMNP1215->Fill(InvariantMass(track1, track2));
         fPtSpectra1215->Fill(pt);
      }
      if ((1.5 <= pt) && (1.8 > pt))
      {
         fInvMNP1518->Fill(InvariantMass(track1, track2));
         fPtSpectra1518->Fill(pt);
      }
      if ((1.8 <= pt) && (2.1 > pt))
      {
         fInvMNP1821->Fill(InvariantMass(track1, track2));
         fPtSpectra1821->Fill(pt);
      }
      if ((2.1 <= pt) && (2.4 > pt))
      {
         fInvMNP2124->Fill(InvariantMass(track1, track2));
         fPtSpectra2124->Fill(pt);
      }
      if ((2.4 <= pt) && (2.7 > pt))
      {
         fInvMNP2427->Fill(InvariantMass(track1, track2));
         fPtSpectra2427->Fill(pt);
      }
      if ((2.7 <= pt) && (3.0 > pt))
      {
         fInvMNP2730->Fill(InvariantMass(track1, track2));
         fPtSpectra2730->Fill(pt);
      }
      if ((3.0 <= pt) && (3.5 > pt))
      {
         fInvMNP3035->Fill(InvariantMass(track1, track2));
         fPtSpectra3035->Fill(pt);
      }
      if ((3.5 <= pt) && (4.0 > pt))
      {
         fInvMNP3540->Fill(InvariantMass(track1, track2));
         fPtSpectra3540->Fill(pt);
      }
      if ((4.0 <= pt) && (4.5 > pt))
      {
         fInvMNP4045->Fill(InvariantMass(track1, track2));
         fPtSpectra4045->Fill(pt);
      }
      if ((4.5 <= pt) && (5.0 > pt))
      {
         fInvMNP4550->Fill(InvariantMass(track1, track2));
         fPtSpectra4550->Fill(pt);
      }
      if ((5.0 <= pt) && (5.5 > pt))
      {
         fInvMNP5055->Fill(InvariantMass(track1, track2));
         fPtSpectra5055->Fill(pt);
      }
      if ((5.5 <= pt) && (6.0 > pt))
      {
         fInvMNP5560->Fill(InvariantMass(track1, track2));
         fPtSpectra5560->Fill(pt);
      }
      if ((6.0 <= pt) && (6.5 > pt))
      {
         fInvMNP6065->Fill(InvariantMass(track1, track2));
         fPtSpectra6065->Fill(pt);
      }
      if ((6.5 <= pt) && (7.0 > pt))
      {
         fInvMNP6570->Fill(InvariantMass(track1, track2));
         fPtSpectra6570->Fill(pt);
      }
   }
   if (tracktype == 1)
   {
      if ((0.0 <= pt) && (0.3 > pt)) fInvMPP03->Fill(InvariantMass(track1, track2));
      if ((0.3 <= pt) && (0.6 > pt)) fInvMPP36->Fill(InvariantMass(track1, track2));
      if ((0.6 <= pt) && (0.9 > pt)) fInvMPP69->Fill(InvariantMass(track1, track2));
      if ((0.9 <= pt) && (1.2 > pt)) fInvMPP912->Fill(InvariantMass(track1, track2));
      if ((1.2 <= pt) && (1.5 > pt)) fInvMPP1215->Fill(InvariantMass(track1, track2));
      if ((1.5 <= pt) && (1.8 > pt)) fInvMPP1518->Fill(InvariantMass(track1, track2));
      if ((1.8 <= pt) && (2.1 > pt)) fInvMPP1821->Fill(InvariantMass(track1, track2));
      if ((2.1 <= pt) && (2.4 > pt)) fInvMPP2124->Fill(InvariantMass(track1, track2));
      if ((2.4 <= pt) && (2.7 > pt)) fInvMPP2427->Fill(InvariantMass(track1, track2));
      if ((2.7 <= pt) && (3.0 > pt)) fInvMPP2730->Fill(InvariantMass(track1, track2));
      if ((3.0 <= pt) && (3.5 > pt)) fInvMPP3035->Fill(InvariantMass(track1, track2));
      if ((3.5 <= pt) && (4.0 > pt)) fInvMPP3540->Fill(InvariantMass(track1, track2));
      if ((4.0 <= pt) && (4.5 > pt)) fInvMPP4045->Fill(InvariantMass(track1, track2));
      if ((4.5 <= pt) && (5.0 > pt)) fInvMPP4550->Fill(InvariantMass(track1, track2));
      if ((5.0 <= pt) && (5.5 > pt)) fInvMPP5055->Fill(InvariantMass(track1, track2));
      if ((5.5 <= pt) && (6.0 > pt)) fInvMPP5560->Fill(InvariantMass(track1, track2));
      if ((6.0 <= pt) && (6.5 > pt)) fInvMPP6065->Fill(InvariantMass(track1, track2));
      if ((6.5 <= pt) && (7.0 > pt)) fInvMPP6570->Fill(InvariantMass(track1, track2));
   }
   if (tracktype == 2)
   {
      if ((0.0 <= pt) && (0.3 > pt)) fInvMNN03->Fill(InvariantMass(track1, track2));
      if ((0.3 <= pt) && (0.6 > pt)) fInvMNN36->Fill(InvariantMass(track1, track2));
      if ((0.6 <= pt) && (0.9 > pt)) fInvMNN69->Fill(InvariantMass(track1, track2));
      if ((0.9 <= pt) && (1.2 > pt)) fInvMNN912->Fill(InvariantMass(track1, track2));
      if ((1.2 <= pt) && (1.5 > pt)) fInvMNN1215->Fill(InvariantMass(track1, track2));
      if ((1.5 <= pt) && (1.8 > pt)) fInvMNN1518->Fill(InvariantMass(track1, track2));
      if ((1.8 <= pt) && (2.1 > pt)) fInvMNN1821->Fill(InvariantMass(track1, track2));
      if ((2.1 <= pt) && (2.4 > pt)) fInvMNN2124->Fill(InvariantMass(track1, track2));
      if ((2.4 <= pt) && (2.7 > pt)) fInvMNN2427->Fill(InvariantMass(track1, track2));
      if ((2.7 <= pt) && (3.0 > pt)) fInvMNN2730->Fill(InvariantMass(track1, track2));
      if ((3.0 <= pt) && (3.5 > pt)) fInvMNN3035->Fill(InvariantMass(track1, track2));
      if ((3.5 <= pt) && (4.0 > pt)) fInvMNN3540->Fill(InvariantMass(track1, track2));
      if ((4.0 <= pt) && (4.5 > pt)) fInvMNN4045->Fill(InvariantMass(track1, track2));
      if ((4.5 <= pt) && (5.0 > pt)) fInvMNN4550->Fill(InvariantMass(track1, track2));
      if ((5.0 <= pt) && (5.5 > pt)) fInvMNN5055->Fill(InvariantMass(track1, track2));
      if ((5.5 <= pt) && (6.0 > pt)) fInvMNN5560->Fill(InvariantMass(track1, track2));
      if ((6.0 <= pt) && (6.5 > pt)) fInvMNN6065->Fill(InvariantMass(track1, track2));
      if ((6.5 <= pt) && (7.0 > pt)) fInvMNN6570->Fill(InvariantMass(track1, track2));
   }
}
//_____________________________________________________________________________
template <typename T> Bool_t AliAnalysisTaskPhiFlow::PhiTrack(T* track) const
{
   // Check if track is suitable for phi flow analysis
   if(!track) return kFALSE;
   if(!fPOICuts->IsSelected(track)) return kFALSE;
   if(fAODAnalysis && fDCA) {
       // ONLY FOR AOD ANALYSIS
       // if flagged force propagation of DCA to primary vertex
       Double_t b[] = { -99., -99.};
       Double_t bCov[] = { -99., -99., -99.};
       track->PropagateToDCA(fAOD->GetPrimaryVertex(), fAOD->GetMagneticField(), 100., b, bCov);
       if( (TMath::Abs(b[0]) > fDCAXY) || (TMath::Abs(b[1]) > fDCAZ) ) return kFALSE;
       fDCAXYQA->Fill(b[0]);
       fDCAZQA->Fill(b[1]);
   }
   return kTRUE;
}
//_____________________________________________________________________________
template <typename T> void AliAnalysisTaskPhiFlow::SetNullCuts(T* event)
{
   // Set null cuts
   if (fDebug) cout << " fCutsRP " << fCutsRP << endl;
   fCutsRP->SetEvent(event, MCEvent());
   fNullCuts->SetParamType(AliFlowTrackCuts::kGlobal);
   fNullCuts->SetPtRange(+1, -1); // select nothing QUICK
   fNullCuts->SetEtaRange(+1, -1); // select nothing VZERO
   fNullCuts->SetEvent(event, MCEvent());
}
//_____________________________________________________________________________
void AliAnalysisTaskPhiFlow::PrepareFlowEvent(Int_t iMulti)
{
   // Prepare flow events
   fFlowEvent->ClearFast();
   fFlowEvent->Fill(fCutsRP, fNullCuts);
   fFlowEvent->SetReferenceMultiplicity(iMulti);
   fFlowEvent->DefineDeadZone(0, 0, 0, 0);
}
//_____________________________________________________________________________
void AliAnalysisTaskPhiFlow::UserExec(Option_t *)
{
   // UserExec: called for each event. Commented where necessary
   // check for AOD data type
   if (!fPIDResponse)
   {
      AliError("Cannot get pid response");
      return;
   }

   fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
   if (fAOD)
   {
      fAODAnalysis = kTRUE;
      // Check whether event passes event cuts
      if (!EventCut(fAOD)) return;
      // initialize event objects
      InitializeBayesianPID(fAOD);
      SetNullCuts(fAOD);
      PrepareFlowEvent(fAOD->GetNumberOfTracks());
      fCandidates->SetLast(-1);
      // Calculate event plane Q vectors and event plane resolution
      fEventStats->Fill(0);
      Int_t unTracks = fAOD->GetNumberOfTracks();
      AliAODTrack* un[unTracks];
      AliAODTrack* up[unTracks];
      Int_t unp(0);
      Int_t unn(0);
      // Loop through tracks, check for species (Kaons), fill arrays according to charge
      for (Int_t iTracks = 0; iTracks < unTracks; iTracks++)
      {
         AliAODTrack* track = fAOD->GetTrack(iTracks);
         if (!PhiTrack(track)) continue;
         if (fStrictKaonCuts && (!PassesStrictKaonCuts(track))) continue;
         Bool_t charge = kFALSE;
         if (track->Charge() > 0)
         {
            charge = kTRUE;
            fEventStats->Fill(1);
            fPtP->Fill(track->Pt());
         }
         if (track->Charge() < 0)
         {
            fEventStats->Fill(2);
            fPtN->Fill(track->Pt());
         }
         if (IsKaon(track))
         {
            if (charge)
            {
               up[unp] = track;
               unp++;
               fEventStats->Fill(3);
               fPtKP->Fill(track->Pt());
            }
            if (!charge)
            {
               un[unn] = track;
               unn++;
               fEventStats->Fill(4);
               fPtKN->Fill(track->Pt());
            }
         }
      }
      // Calculate invariant mass of like- and unlike sign pairs as a function of p_T, store data in histograms
      for (Int_t pTracks = 0; pTracks < unp ; pTracks++)
      {
         for (Int_t nTracks = 0; nTracks < unn ; nTracks++)
         {
            if (fApplyDeltaDipCut && (!CheckDeltaDipAngle(up[pTracks], un[nTracks]))) continue;
            if (fCandidateEtaPtCut && (!CheckCandidateEtaPtCut(up[pTracks], un[nTracks]))) continue;
            PtSelector(0, up[pTracks], un[nTracks]);
            Double_t pt = PhiPt(up[pTracks], un[nTracks]);
            Double_t mass = InvariantMass(up[pTracks], un[nTracks]);
            TVector3 a(up[pTracks]->Px(), up[pTracks]->Py(), up[pTracks]->Pz());
            TVector3 b(un[nTracks]->Px(), un[nTracks]->Py(), up[pTracks]->Pz());
            TVector3 c = a + b;
            Double_t phi = c.Phi();
            Double_t eta = c.Eta();
            Int_t nIDs[2];
            nIDs[0] = up[pTracks]->GetID();
            nIDs[1] = un[nTracks]->GetID();
            MakeTrack(mass, pt, phi, eta, 2, nIDs);
         }
      }

      if (fDebug)  printf("I received %d candidates\n", fCandidates->GetEntriesFast());
      for (int iCand = 0; iCand != fCandidates->GetEntriesFast(); ++iCand)
      {
         AliFlowCandidateTrack *cand = dynamic_cast<AliFlowCandidateTrack*>(fCandidates->At(iCand));
         if (!cand) continue;
         if (fDebug) printf(" >Checking at candidate %d with %d daughters: mass %f\n", iCand, cand->GetNDaughters(), cand->Mass());
         for (int iDau = 0; iDau != cand->GetNDaughters(); ++iDau)
         {
            if (fDebug) printf("  >Daughter %d with fID %d", iDau, cand->GetIDDaughter(iDau));
            for (int iRPs = 0; iRPs != fFlowEvent->NumberOfTracks(); ++iRPs)
            {
               AliFlowTrack *iRP = dynamic_cast<AliFlowTrack*>(fFlowEvent->GetTrack(iRPs));
               if (!iRP) continue;
               if (!iRP->InRPSelection()) continue;
               if (cand->GetIDDaughter(iDau) == iRP->GetID())
               {
                  if (fDebug) printf(" was in RP set");
                  iRP->SetForRPSelection(kFALSE);
                  fFlowEvent->SetNumberOfRPs(fFlowEvent->GetNumberOfRPs() - 1);
               }
            }
            if (fDebug) printf("\n");
         }
         cand->SetForPOISelection(kTRUE);
         fFlowEvent->InsertTrack(((AliFlowTrack*) cand));
      }
      if (fDebug) printf("TPCevent %d\n", fFlowEvent->NumberOfTracks());


      for (Int_t pTracks = 0; pTracks < unp ; pTracks++)
      {
         for (Int_t nTracks = pTracks + 1; nTracks < unp ; nTracks++)
         {
            if (fApplyDeltaDipCut && (!CheckDeltaDipAngle(up[pTracks], up[nTracks]))) continue;
            if (fCandidateEtaPtCut && (!CheckCandidateEtaPtCut(up[pTracks], up[nTracks]))) continue;
            PtSelector(1, up[pTracks], up[nTracks]);
         }
      }
      for (Int_t nTracks = 0; nTracks < unn ; nTracks++)
      {
         for (Int_t pTracks = nTracks + 1; pTracks < unn ; pTracks++)
         {
            if (fApplyDeltaDipCut && (!CheckDeltaDipAngle(un[nTracks], un[pTracks]))) continue;
            if (fCandidateEtaPtCut && (!CheckCandidateEtaPtCut(un[nTracks], un[pTracks]))) continue;
            PtSelector(2, un[nTracks], un[pTracks]);
         }
      }
      PostData(1, fOutputList);
      PostData(2, fFlowEvent);
   }

   fESD = dynamic_cast<AliESDEvent*>(InputEvent());
   if (fESD)
   {
      fAODAnalysis = kFALSE;
      // Check whether event passes event cuts
      if (!EventCut(fESD)) return;
      InitializeBayesianPID(fESD);
      SetNullCuts(fESD);
      PrepareFlowEvent(fESD->GetNumberOfTracks());
      // Calculate event plane Q vectors and event plane resolution
      fEventStats->Fill(0);
      Int_t unTracks = fESD->GetNumberOfTracks();
      AliESDtrack* un[unTracks];
      AliESDtrack* up[unTracks];
      Int_t unp(0);
      Int_t unn(0);
      // Loop through tracks, check for species (Kaons), fill arrays according to charge
      for (Int_t iTracks = 0; iTracks < unTracks; iTracks++)
      {
         AliESDtrack* track = fESD->GetTrack(iTracks);
         if (!PhiTrack(track)) continue;
         Bool_t charge = kFALSE;
         if (track->Charge() > 0)
         {
            charge = kTRUE;
            fEventStats->Fill(1);
            fPtP->Fill(track->Pt());
         }
         if (track->Charge() < 0)
         {
            fEventStats->Fill(2);
            fPtN->Fill(track->Pt());
         }
         if (IsKaon(track))
         {
            if (charge)
            {
               up[unp] = track;
               unp++;
               fEventStats->Fill(3);
               fPtKP->Fill(track->Pt());
            }
            if (!charge)
            {
               un[unn] = track;
               unn++;
               fEventStats->Fill(4);
               fPtKN->Fill(track->Pt());
            }
         }
      }
      // Calculate invariant mass of like- and unlike sign pairs as a function of p_T, store data in histograms
        for (Int_t pTracks = 0; pTracks < unp ; pTracks++)
      {
         for (Int_t nTracks = 0; nTracks < unn ; nTracks++)
         {
            if (fApplyDeltaDipCut && (!CheckDeltaDipAngle(up[pTracks], un[nTracks]))) continue;
            if (fCandidateEtaPtCut && (!CheckCandidateEtaPtCut(up[pTracks], un[nTracks]))) continue;
            PtSelector(0, up[pTracks], un[nTracks]);
            Double_t pt = PhiPt(up[pTracks], un[nTracks]);
            Double_t mass = InvariantMass(up[pTracks], un[nTracks]);
            TVector3 a(up[pTracks]->Px(), up[pTracks]->Py(), up[pTracks]->Pz());
            TVector3 b(un[nTracks]->Px(), un[nTracks]->Py(), up[pTracks]->Pz());
            TVector3 c = a + b;
            Double_t phi = c.Phi();
            Double_t eta = c.Eta();
            Int_t nIDs[2];
            nIDs[0] = up[pTracks]->GetID();
            nIDs[1] = un[nTracks]->GetID();
            MakeTrack(mass, pt, phi, eta, 2, nIDs);
         }
      }

      if (fDebug)  printf("I received %d candidates\n", fCandidates->GetEntriesFast());
      for (int iCand = 0; iCand != fCandidates->GetEntriesFast(); ++iCand)
      {
         AliFlowCandidateTrack *cand = dynamic_cast<AliFlowCandidateTrack*>(fCandidates->At(iCand));
         if (!cand) continue;
         if (fDebug) printf(" >Checking at candidate %d with %d daughters: mass %f\n", iCand, cand->GetNDaughters(), cand->Mass());
         for (int iDau = 0; iDau != cand->GetNDaughters(); ++iDau)
         {
            if (fDebug) printf("  >Daughter %d with fID %d", iDau, cand->GetIDDaughter(iDau));
            for (int iRPs = 0; iRPs != fFlowEvent->NumberOfTracks(); ++iRPs)
            {
               AliFlowTrack *iRP = dynamic_cast<AliFlowTrack*>(fFlowEvent->GetTrack(iRPs));
               if (!iRP) continue;
               if (!iRP->InRPSelection()) continue;
               if (cand->GetIDDaughter(iDau) == iRP->GetID())
               {
                  if (fDebug) printf(" was in RP set");
                  iRP->SetForRPSelection(kFALSE);
                  fFlowEvent->SetNumberOfRPs(fFlowEvent->GetNumberOfRPs() - 1);
               }
            }
            if (fDebug) printf("\n");
         }
         cand->SetForPOISelection(kTRUE);
         fFlowEvent->InsertTrack(((AliFlowTrack*) cand));
      }
      if (fDebug) printf("TPCevent %d\n", fFlowEvent->NumberOfTracks());

      for (Int_t pTracks = 0; pTracks < unp ; pTracks++)
      {
         for (Int_t nTracks = pTracks + 1; nTracks < unp ; nTracks++)
         {
            if (fApplyDeltaDipCut && (!CheckDeltaDipAngle(up[pTracks], up[nTracks]))) continue;
            if (fCandidateMinEta && (!CheckCandidateEtaPtCut(up[pTracks], up[nTracks]))) continue;
            PtSelector(1, up[pTracks], up[nTracks]);
         }
      }
      for (Int_t nTracks = 0; nTracks < unn ; nTracks++)
      {
         for (Int_t pTracks = nTracks + 1; pTracks < unn ; pTracks++)
         {
            if (fApplyDeltaDipCut && (!CheckDeltaDipAngle(un[nTracks], un[pTracks]))) continue;
            if (fCandidateEtaPtCut && (!CheckCandidateEtaPtCut(un[nTracks], un[pTracks]))) continue;
            PtSelector(2, un[nTracks], un[pTracks]);
         }
      }
      PostData(1, fOutputList);
      PostData(2, fFlowEvent);
   }
}
//_____________________________________________________________________________
void AliAnalysisTaskPhiFlow::Terminate(Option_t *)
{
   // Terminate
}
//______________________________________________________________________________
void  AliAnalysisTaskPhiFlow::MakeTrack(Double_t mass,
                                              Double_t pt, Double_t phi, Double_t eta,
                                              Int_t nDau, Int_t iID[]) const
{
   // Consruct Flow Candidate Track from two selected candidates
   Bool_t overwrite = kTRUE;
   AliFlowCandidateTrack *sTrack = static_cast<AliFlowCandidateTrack*>(fCandidates->At(fCandidates->GetLast() + 1));
   if (!sTrack)
   {
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
void AliAnalysisTaskPhiFlow::SetCommonConstants(Int_t massBins, Double_t minMass, Double_t maxMass) 
{
  // set common constants
  fMassBins = massBins;
  fMinMass = minMass;
  fMaxMass = maxMass;
}
//_____________________________________________________________________________

