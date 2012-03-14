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
#include "AliESDpid.h"
#include "AliAODpidUtil.h"

class AliFlowTrackCuts;

ClassImp(AliAnalysisTaskPhiFlow)

AliAnalysisTaskPhiFlow::AliAnalysisTaskPhiFlow() : AliAnalysisTaskSE(),
   fAODAnalysis(0), fEtaMinA(0), fEtaMaxA(0), fEtaMinB(0), fEtaMaxB(0), fCutsRP(NULL), fNullCuts(0), fBayesianResponse(0), fOldTrackParam(0), fRequireTPCStandAlone(0), fStrictKaonCuts(0), fCandidateEtaPtCut(0), fCandidateMinEta(0), fCandidateMaxEta(0), fCandidateMinPt(0), fCandidateMaxPt(0), fParticleProbability(0.5), fCentrality(0), fESD(0), fAOD(0), fOutputList(0), fEventStats(0), fCentralityPass(0), fCentralityNoPass(0), fNOPID(0), fPIDk(0), fPIDDeltaDip(0), fInvMNP03(0), fInvMPP03(0), fInvMNN03(0), fInvMNP36(0), fInvMPP36(0), fInvMNN36(0), fInvMNP69(0), fInvMPP69(0), fInvMNN69(0), fInvMNP912(0), fInvMPP912(0), fInvMNN912(0), fInvMNP1215(0), fInvMPP1215(0), fInvMNN1215(0), fInvMNP1518(0), fInvMPP1518(0), fInvMNN1518(0), fInvMNP1821(0), fInvMPP1821(0), fInvMNN1821(0), fInvMNP2124(0), fInvMPP2124(0), fInvMNN2124(0), fInvMNP2427(0), fInvMPP2427(0), fInvMNN2427(0), fInvMNP2730(0), fInvMPP2730(0), fInvMNN2730(0), fInvMNP3035(0), fInvMPP3035(0), fInvMNN3035(0), fInvMNP3540(0), fInvMPP3540(0), fInvMNN3540(0), fInvMNP4045(0), fInvMPP4045(0), fInvMNN4045(0), fInvMNP4550(0), fInvMPP4550(0), fInvMNN4550(0), fInvMNP5055(0), fInvMPP5055(0), fInvMNN5055(0), fInvMNP5560(0), fInvMPP5560(0), fInvMNN5560(0), fInvMNP6065(0), fInvMPP6065(0), fInvMNN6065(0), fInvMNP6570(0), fInvMPP6570(0), fInvMNN6570(0), fDeltaPhiPsiNP03(0), fDeltaPhiPsiNP36(0), fDeltaPhiPsiNP69(0), fDeltaPhiPsiNP912(0), fDeltaPhiPsiNP1215(0), fDeltaPhiPsiNP1518(0), fDeltaPhiPsiNP1821(0), fDeltaPhiPsiNP2124(0), fDeltaPhiPsiNP2427(0), fDeltaPhiPsiNP2730(0), fDeltaPhiPsiNP3035(0), fDeltaPhiPsiNP3540(0), fDeltaPhiPsiNP4045(0), fDeltaPhiPsiNP4550(0), fDeltaPhiPsiNP5055(0), fDeltaPhiPsiNP5560(0), fDeltaPhiPsiNP6065(0), fDeltaPhiPsiNP6570(0), fProfV2(0), fProfV2Sin(0), fProfV2InvM03(0), fProfV2InvM36(0), fProfV2InvM69(0), fProfV2InvM912(0), fProfV2InvM1215(0), fProfV2InvM1518(0), fProfV2InvM1821(0), fProfV2InvM2124(0), fProfV2InvM2427(0), fProfV2InvM2730(0), fProfV2InvM3035(0), fProfV2InvM3540(0), fProfV2InvM4045(0), fProfV2InvM4550(0), fProfV2InvM5055(0), fProfV2InvM5560(0), fProfV2InvM6065(0), fProfV2InvM6570(0), fProfV2SinInvM03(0), fProfV2SinInvM36(0), fProfV2SinInvM69(0), fProfV2SinInvM912(0), fProfV2SinInvM1215(0), fProfV2SinInvM1518(0), fProfV2SinInvM1821(0), fProfV2SinInvM2124(0), fProfV2SinInvM2427(0), fProfV2SinInvM2730(0), fProfV2SinInvM3035(0), fProfV2SinInvM3540(0), fProfV2SinInvM4045(0), fProfV2SinInvM4550(0), fProfV2SinInvM5055(0), fProfV2SinInvM5560(0), fProfV2SinInvM6065(0), fProfV2SinInvM6570(0), fPtSpectra03(0), fPtSpectra36(0), fPtSpectra69(0), fPtSpectra912(0), fPtSpectra1215(0), fPtSpectra1518(0), fPtSpectra1821(0), fPtSpectra2124(0), fPtSpectra2427(0), fPtSpectra2730(0), fPtSpectra3035(0), fPtSpectra3540(0), fPtSpectra4045(0), fPtSpectra4550(0), fPtSpectra5055(0), fPtSpectra5560(0), fPtSpectra6065(0), fPtSpectra6570(0), fEventPlaneSTAR(0), fEventPlaneResolutionRandom(0), fEventPlaneResolutionEta(0), fPtP(0), fPtN(0), fPtKP(0), fPtKN(0), fCentralityMin(0), fCentralityMax(100), fkCentralityMethod(0), fRPCuts(0), fPOICuts(0), fVertexRange(0), fQx(0), fQy(0), fEventPlane(0), fPhi(0), fPt(0), fEta(0), fVZEROA(0), fVZEROC(0), fTPCM(0), fPIDtype(kCombined), fDeltaDipAngle(0), fDeltaDipPt(0), fApplyDeltaDipCut(0), fPairLoss(0), fEventPlanePtCut(0), fSetEventPlanePtCut(0)
{
   // Dummy constructor
   fNullCuts = new AliFlowTrackCuts("null_cuts");
   for (Int_t i = 0; i < 2; i++)
   {
      for (Int_t j = 0; j < 30; j++)
      {
         fFlowBands[i][j] = 0;
         if (i == 0)  fFlowEvent[j] = new AliFlowEvent(10000);
      }
   }
   if (fPIDtype == kCombined) fBayesianResponse = new AliFlowBayesianPID();
}
//_____________________________________________________________________________
AliAnalysisTaskPhiFlow::AliAnalysisTaskPhiFlow(const char *name) : AliAnalysisTaskSE(name),
   fAODAnalysis(0), fEtaMinA(0), fEtaMaxA(0), fEtaMinB(0), fEtaMaxB(0), fCutsRP(NULL), fNullCuts(0), fBayesianResponse(0), fOldTrackParam(0), fRequireTPCStandAlone(0), fStrictKaonCuts(0), fCandidateEtaPtCut(0), fCandidateMinEta(0), fCandidateMaxEta(0), fCandidateMinPt(0), fCandidateMaxPt(0), fParticleProbability(0.5), fCentrality(0), fESD(0), fAOD(0),fOutputList(0), fEventStats(0), fCentralityPass(0), fCentralityNoPass(0), fNOPID(0), fPIDk(0), fPIDDeltaDip(0), fInvMNP03(0), fInvMPP03(0), fInvMNN03(0), fInvMNP36(0), fInvMPP36(0), fInvMNN36(0), fInvMNP69(0), fInvMPP69(0), fInvMNN69(0), fInvMNP912(0), fInvMPP912(0), fInvMNN912(0), fInvMNP1215(0), fInvMPP1215(0), fInvMNN1215(0), fInvMNP1518(0), fInvMPP1518(0), fInvMNN1518(0), fInvMNP1821(0), fInvMPP1821(0), fInvMNN1821(0), fInvMNP2124(0), fInvMPP2124(0), fInvMNN2124(0), fInvMNP2427(0), fInvMPP2427(0), fInvMNN2427(0), fInvMNP2730(0), fInvMPP2730(0), fInvMNN2730(0), fInvMNP3035(0), fInvMPP3035(0), fInvMNN3035(0), fInvMNP3540(0), fInvMPP3540(0), fInvMNN3540(0), fInvMNP4045(0), fInvMPP4045(0), fInvMNN4045(0), fInvMNP4550(0), fInvMPP4550(0), fInvMNN4550(0), fInvMNP5055(0), fInvMPP5055(0), fInvMNN5055(0), fInvMNP5560(0), fInvMPP5560(0), fInvMNN5560(0), fInvMNP6065(0), fInvMPP6065(0), fInvMNN6065(0), fInvMNP6570(0), fInvMPP6570(0), fInvMNN6570(0), fDeltaPhiPsiNP03(0), fDeltaPhiPsiNP36(0), fDeltaPhiPsiNP69(0), fDeltaPhiPsiNP912(0), fDeltaPhiPsiNP1215(0), fDeltaPhiPsiNP1518(0), fDeltaPhiPsiNP1821(0), fDeltaPhiPsiNP2124(0), fDeltaPhiPsiNP2427(0), fDeltaPhiPsiNP2730(0), fDeltaPhiPsiNP3035(0), fDeltaPhiPsiNP3540(0), fDeltaPhiPsiNP4045(0), fDeltaPhiPsiNP4550(0), fDeltaPhiPsiNP5055(0), fDeltaPhiPsiNP5560(0), fDeltaPhiPsiNP6065(0), fDeltaPhiPsiNP6570(0), fProfV2(0), fProfV2Sin(0), fProfV2InvM03(0), fProfV2InvM36(0), fProfV2InvM69(0), fProfV2InvM912(0), fProfV2InvM1215(0), fProfV2InvM1518(0), fProfV2InvM1821(0), fProfV2InvM2124(0), fProfV2InvM2427(0), fProfV2InvM2730(0), fProfV2InvM3035(0), fProfV2InvM3540(0), fProfV2InvM4045(0), fProfV2InvM4550(0), fProfV2InvM5055(0), fProfV2InvM5560(0), fProfV2InvM6065(0), fProfV2InvM6570(0), fProfV2SinInvM03(0), fProfV2SinInvM36(0), fProfV2SinInvM69(0), fProfV2SinInvM912(0), fProfV2SinInvM1215(0), fProfV2SinInvM1518(0), fProfV2SinInvM1821(0), fProfV2SinInvM2124(0), fProfV2SinInvM2427(0), fProfV2SinInvM2730(0), fProfV2SinInvM3035(0), fProfV2SinInvM3540(0), fProfV2SinInvM4045(0), fProfV2SinInvM4550(0), fProfV2SinInvM5055(0), fProfV2SinInvM5560(0), fProfV2SinInvM6065(0), fProfV2SinInvM6570(0), fPtSpectra03(0), fPtSpectra36(0), fPtSpectra69(0), fPtSpectra912(0), fPtSpectra1215(0), fPtSpectra1518(0), fPtSpectra1821(0), fPtSpectra2124(0), fPtSpectra2427(0), fPtSpectra2730(0), fPtSpectra3035(0), fPtSpectra3540(0), fPtSpectra4045(0), fPtSpectra4550(0), fPtSpectra5055(0), fPtSpectra5560(0), fPtSpectra6065(0), fPtSpectra6570(0), fEventPlaneSTAR(0), fEventPlaneResolutionRandom(0), fEventPlaneResolutionEta(0), fPtP(0), fPtN(0), fPtKP(0), fPtKN(0), fCentralityMin(0), fCentralityMax(100), fkCentralityMethod(0), fRPCuts(0), fPOICuts(0), fVertexRange(0), fQx(0), fQy(0), fEventPlane(0), fPhi(0), fPt(0), fEta(0), fVZEROA(0), fVZEROC(0), fTPCM(0), fPIDtype(kCombined), fDeltaDipAngle(0), fDeltaDipPt(0), fApplyDeltaDipCut(0), fPairLoss(0), fEventPlanePtCut(0), fSetEventPlanePtCut(0)
{
   // Constructor
   fNullCuts = new AliFlowTrackCuts("null_cuts");
   for (int i = 0; i < 2; i++)
   {
      for (int j = 0; j < 30; j++)
      {
         fFlowBands[i][j] = 0;
         if (i == 0)  fFlowEvent[j] = new AliFlowEvent(10000);
      }
   }
   if (fPIDtype == kCombined) fBayesianResponse = new AliFlowBayesianPID();
   DefineInput(0, TChain::Class());
   DefineOutput(1, TList::Class());
   for (Int_t i = (0 + 2); i < (30 + 2); i++)
   {
      DefineOutput(i, AliFlowEventSimple::Class());
   }
}
//_____________________________________________________________________________
AliAnalysisTaskPhiFlow::~AliAnalysisTaskPhiFlow()
{
   // Destructor
   for (Int_t i = 0; i < 30; i++)
   {
      delete fFlowEvent[i];
   }
   delete fNullCuts;
   if (fPIDtype == kCombined) delete fBayesianResponse;
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
TH1F* AliAnalysisTaskPhiFlow::BookDPhiPsiHistogram(const char* name)
{
   // Return a pointer to a TH1 with predefined binning
   TH1F *hist = new TH1F(Form("#phi - #Psi (%s)", name), Form("#phi - #Psi (%s)", name), 60, -7, 7);
   hist->GetXaxis()->SetTitle("#phi - #Psi (rad)");
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
   (i < 11) ? nmin = 0.3 * (i-1) : nmin = 0.5 * (i-11) + 3;
   (i < 11) ? nmax = 0.3 * (i-1) + 0.3 : nmax = 0.5 * (i-11) + 3.5;
   TH1F* hist = new TH1F(Form("%f p_{t} %f", nmin, nmax), Form("%f p_{t} %f", nmin, nmax), 60, nmin, nmax);
   hist->GetXaxis()->SetTitle("p_{T} GeV / c");
   fOutputList->Add(hist);
   return hist;
}
//_____________________________________________________________________________
TProfile* AliAnalysisTaskPhiFlow::BookV2Profile(const char* name, Bool_t pt, Bool_t cos)
{
   // Return a pointer to a v_2 TProfile with predefined binning
   TProfile* profile = 0x0;
   if (pt)
   {
      profile = new TProfile(name, Form("v_{c, 2}^{pair} (%s)", name), 100, 0., 4);
      profile->GetXaxis()->SetTitle("p_{T} GeV / c");
   }
   if (!pt)
   {
      if (cos)
      {
         profile = new TProfile(Form("v_{c, 2}^{pair} (%s)", name), Form("v_{c, 2}^{pair} (%s)", name), 30, 0.99, 1.092);
         profile->GetXaxis()->SetTitle("M_{INV} GeV / c^{2}");
      }
      if (!cos)
      {
         profile = new TProfile(Form("v_{s, 2}^{pair} (%s)", name), Form("v_{s, 2}^{pair} (%s)", name), 30, 0.99, 1.092);
         profile->GetXaxis()->SetTitle("M_{INV} GeV / c^{2}");
      }
   }
   profile->SetMarkerStyle(kOpenCross);
   if (cos) profile->GetYaxis()->SetTitle("v_{2} <cos n[#phi_{pair} - #Psi_{RP}]>");
   if (!cos) profile->GetYaxis()->SetTitle("v_{2} <sin n[#phi_{pair} - #Psi_{RP}]>");
   profile->Sumw2();
   fOutputList->Add(profile);
   return profile;
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

   fEventPlaneSTAR = new TH1F("fEventPlaneSTAR", "Event plane orientation", 200, -.2, 4);
   fEventPlaneSTAR->GetXaxis()->SetTitle("#Psi_{2}");
   fEventPlaneSTAR->GetYaxis()->SetTitle("Events");
   fEventPlaneResolutionRandom = new TProfile("fEventPlaneResolutionRandom", "RP kaon pairs, Random", 1, 0.5, 1.5);
   fEventPlaneResolutionEta = new TProfile("fEventPlaneResolutionEta", "RP kaon pairs, Eta", 1, 0.5, 1.5);

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
   fProfV2 = BookV2Profile("v_{2} cos terms", kTRUE, kTRUE);
   fProfV2Sin = BookV2Profile("v_{2} sin terms", kTRUE, kFALSE);
   fProfV2InvM03 = BookV2Profile("0.0 < p_{T} < 0.3 GeV", kFALSE, kTRUE);
   fProfV2InvM36 = BookV2Profile("0.3 < p_{T} < 0.6 GeV", kFALSE, kTRUE);
   fProfV2InvM69 = BookV2Profile("0.6 < p_{T} < 0.9 GeV", kFALSE, kTRUE);
   fProfV2InvM912 = BookV2Profile("0.9 < p_{T} < 1.2 GeV", kFALSE, kTRUE);
   fProfV2InvM1215 = BookV2Profile("1.1 < p_{T} < 1.5 GeV", kFALSE, kTRUE);
   fProfV2InvM1518 = BookV2Profile("1.5 < p_{T} < 1.8 GeV", kFALSE, kTRUE);
   fProfV2InvM1821 = BookV2Profile("1.8 < p_{T} < 2.1 GeV", kFALSE, kTRUE);
   fProfV2InvM2124 = BookV2Profile("2.1 < p_{T} < 2.4 GeV", kFALSE, kTRUE);
   fProfV2InvM2427 = BookV2Profile("2.4 < p_{T} < 2.7 GeV", kFALSE, kTRUE);
   fProfV2InvM2730 = BookV2Profile("2.7 < p_{T} < 3.0 GeV", kFALSE, kTRUE);
   fProfV2InvM3035 = BookV2Profile("3.0 < p_{T} < 3.5 GeV", kFALSE, kTRUE);
   fProfV2InvM3540 = BookV2Profile("3.5 < p_{T} < 4.0 GeV", kFALSE, kTRUE);
   fProfV2InvM4045 = BookV2Profile("4.0 < p_{T} < 4.5 GeV", kFALSE, kTRUE);
   fProfV2InvM4550 = BookV2Profile("4.5 < p_{T} < 5.0 GeV", kFALSE, kTRUE);
   fProfV2InvM5055 = BookV2Profile("5.0 < p_{T} < 5.5 GeV", kFALSE, kTRUE);
   fProfV2InvM5560 = BookV2Profile("5.5 < p_{T} < 6.0 GeV", kFALSE, kTRUE);
   fProfV2InvM6065 = BookV2Profile("6.0 < p_{T} < 6.5 GeV", kFALSE, kTRUE);
   fProfV2InvM6570 = BookV2Profile("6.5 < p_{T} < 7.0 GeV", kFALSE, kTRUE);
   fProfV2SinInvM03 = BookV2Profile("0.0 < p_{T} < 0.3 GeV", kFALSE, kFALSE);
   fProfV2SinInvM36 = BookV2Profile("0.3 < p_{T} < 0.6 GeV", kFALSE, kFALSE);
   fProfV2SinInvM69 = BookV2Profile("0.6 < p_{T} < 0.9 GeV", kFALSE, kFALSE);
   fProfV2SinInvM912 = BookV2Profile("0.9 < p_{T} < 1.2 GeV", kFALSE, kFALSE);
   fProfV2SinInvM1215 = BookV2Profile("1.2 < p_{T} < 1.5 GeV", kFALSE, kFALSE);
   fProfV2SinInvM1518 = BookV2Profile("1.5 < p_{T} < 1.8 GeV", kFALSE, kFALSE);
   fProfV2SinInvM1821 = BookV2Profile("1.8 < p_{T} < 2.1 GeV", kFALSE, kFALSE);
   fProfV2SinInvM2124 = BookV2Profile("2.1 < p_{T} < 2.4 GeV", kFALSE, kFALSE);
   fProfV2SinInvM2427 = BookV2Profile("2.4 < p_{T} < 2.7 GeV", kFALSE, kFALSE);
   fProfV2SinInvM2730 = BookV2Profile("2.5 < p_{T} < 3.0 GeV", kFALSE, kFALSE);
   fProfV2SinInvM3035 = BookV2Profile("3.0 < p_{T} < 3.5 GeV", kFALSE, kFALSE);
   fProfV2SinInvM3540 = BookV2Profile("3.5 < p_{T} < 4.0 GeV", kFALSE, kFALSE);
   fProfV2SinInvM4045 = BookV2Profile("4.0 < p_{T} < 4.5 GeV", kFALSE, kFALSE);
   fProfV2SinInvM4550 = BookV2Profile("4.5 < p_{T} < 5.0 GeV", kFALSE, kFALSE);
   fProfV2SinInvM5055 = BookV2Profile("5.0 < p_{T} < 5.5 GeV", kFALSE, kFALSE);
   fProfV2SinInvM5560 = BookV2Profile("5.5 < p_{T} < 6.0 GeV", kFALSE, kFALSE);
   fProfV2SinInvM6065 = BookV2Profile("6.0 < p_{T} < 6.5 GeV", kFALSE, kFALSE);
   fProfV2SinInvM6570 = BookV2Profile("6.5 < p_{T} < 7.0 GeV", kFALSE, kFALSE);

   fDeltaPhiPsiNP03 = BookDPhiPsiHistogram("NP, 0 < p_{T} < 0.3 GeV");
   fDeltaPhiPsiNP36 = BookDPhiPsiHistogram("NP, 0.3 < p_{T} < 0.6 GeV");
   fDeltaPhiPsiNP69 = BookDPhiPsiHistogram("NP, 0.6 < p_{T} < 0.9 GeV");
   fDeltaPhiPsiNP912 = BookDPhiPsiHistogram("NP, 0.9 < p_{T} < 1.2 GeV");
   fDeltaPhiPsiNP1215 = BookDPhiPsiHistogram("NP, 1.2 < p_{T} < 1.5 GeV");
   fDeltaPhiPsiNP1518 = BookDPhiPsiHistogram("NP, 1.5 < p_{T} < 1.8 GeV");
   fDeltaPhiPsiNP1821 = BookDPhiPsiHistogram("NP, 1.8 < p_{T} < 2.1 GeV");
   fDeltaPhiPsiNP2124 = BookDPhiPsiHistogram("NP, 2.1 < p_{T} < 2.4 GeV");
   fDeltaPhiPsiNP2427 = BookDPhiPsiHistogram("NP, 2.4 < p_{T} < 2.7 GeV");
   fDeltaPhiPsiNP2730 = BookDPhiPsiHistogram("NP, 2.7 < p_{T} < 3.0 GeV");
   fDeltaPhiPsiNP3035 = BookDPhiPsiHistogram("NP, 3.0 < p_{T} < 3.5 GeV");
   fDeltaPhiPsiNP3540 = BookDPhiPsiHistogram("NP, 3.5 < p_{T} < 4.0 GeV");
   fDeltaPhiPsiNP4045 = BookDPhiPsiHistogram("NP, 4.0 < p_{T} < 4.5 GeV");
   fDeltaPhiPsiNP4550 = BookDPhiPsiHistogram("NP, 4.5 < p_{T} < 5.0 GeV");
   fDeltaPhiPsiNP5055 = BookDPhiPsiHistogram("NP, 5.0 < p_{T} < 5.5 GeV");
   fDeltaPhiPsiNP5560 = BookDPhiPsiHistogram("NP, 5.5 < p_{T} < 6.0 GeV");
   fDeltaPhiPsiNP6065 = BookDPhiPsiHistogram("NP, 6.0 < p_{T} < 6.5 GeV");
   fDeltaPhiPsiNP6570 = BookDPhiPsiHistogram("NP, 6.5 < p_{T} < 7.0 GeV");


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

   fOutputList->Add(fEventPlaneSTAR);
   fPtP = BookPtHistogram("i^{+}");
   fPtN = BookPtHistogram("i^{-}");
   fPtKP = BookPtHistogram("K^{+}");
   fPtKN = BookPtHistogram("K^{-}");

   fEventPlane = new TH1F("fEventPlane", "Event plane", 200, -0.2, 4.);
   fEventPlane->GetXaxis()->SetTitle("#Psi_{2}");
   fEventPlane->GetYaxis()->SetTitle("Events");
   fOutputList->Add(fEventPlane);

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

   fPairLoss = new TH2F("fPairLoss", "DphiDeta", 100, -0.2, 0.2, 100, -0.2, 0.2);
   fPairLoss->GetXaxis()->SetTitle("#Delta #eta");
   fPairLoss->GetYaxis()->SetTitle("#Delta #phi");

   fOutputList->Add(fEventPlaneResolutionRandom);
   fOutputList->Add(fEventPlaneResolutionEta);
   fOutputList->Add(fPairLoss);
   fPIDDeltaDip = BookPIDHistogram("TPC signal after delta dip cut");
}
//_____________________________________________________________________________
void AliAnalysisTaskPhiFlow::UserCreateOutputObjects()
{
   // Create user defined output objects

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

   // setup initial state of PID response object
   if (!fOldTrackParam) fBayesianResponse->SetNewTrackParam();

   // Create all output objects and store them to a list
   fOutputList = new TList();
   fOutputList->SetOwner();

   // Create and post the Phi identification output objects
   AddPhiIdentificationOutputObjects();
   PostData(1, fOutputList);

   // post flow events (see ctor)
   for (Int_t i = 0; i < 30; i++) PostData(2 + i, fFlowEvent[i]);
}
//_____________________________________________________________________________
template <typename T> void AliAnalysisTaskPhiFlow::PairLoss(const T* track1, const T* track2) const
{
   if (track1->Pt() > track2->Pt()) fPairLoss->Fill((track1->Eta() - track2->Eta()), (track1->Phi() - track2->Phi() + TMath::ASin(0.075 * (1.2 / track2->Pt())) - TMath::ASin(0.075 * (1.2 / track1->Pt()))));
   if (track2->Pt() > track1->Pt()) fPairLoss->Fill((track2->Eta() - track1->Eta()), (track2->Phi() - track1->Phi() + TMath::ASin(0.075 * (1.2 / track1->Pt())) - TMath::ASin(0.075 * (1.2 / track2->Pt()))));
}
//_____________________________________________________________________________
template <typename T> Double_t AliAnalysisTaskPhiFlow::InvariantMass(const T* track1, const T* track2) const
{
   // Return the invariant mass of two tracks, assuming both tracks are kaons
   if ((!track2) || (!track1)) return 0.;
   PairLoss(track1, track2);
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
   // Calculate the delta dip angle between two particles (the opening angle of a pair
   // in the p_t p_z plane)
   // A cut can be used to remove electron / positron contamination from photon conversion
   // 999 (no cut) is returned when angle cannot be calculated
   if (track1->P()*track2->P() == 0) return 999;
   return TMath::ACos(((track1->Pt() * track2->Pt()) + (track1->Pz() * track2->Pz())) / (track1->P() * track2->P()));
}
//_____________________________________________________________________________
template <typename T> Bool_t AliAnalysisTaskPhiFlow::CheckDeltaDipAngle(const T* track1, const T* track2) const
{
   // Check if pair passes delta dip angle cut within 0 < p_t < fDeltaDipPt
   if ((TMath::Abs(DeltaDipAngle(track1, track2)) < fDeltaDipAngle) && (PhiPt(track1, track2) < fDeltaDipPt)) return kFALSE;
   fPIDDeltaDip->Fill(track1->P(), track1->GetTPCsignal());
   fPIDDeltaDip->Fill(track2->P(), track2->GetTPCsignal());
   return kTRUE;
}
//_____________________________________________________________________________
template <typename T> Bool_t AliAnalysisTaskPhiFlow::CheckCandidateEtaPtCut(const T* track1, const T* track2) const
{
    // Check if pair passes eta and pt cut
    if (fCandidateMinPt > PhiPt(track1, track2) || fCandidateMaxPt < PhiPt(track1, track2) ) return kFALSE;
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
   // miminum bias trigger for AOD's
   if (fAODAnalysis) 
   {
       AliInputEventHandler* eh = dynamic_cast<AliInputEventHandler *>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
       if(!eh) return kFALSE;
       if(!(eh->IsEventSelected() & AliVEvent::kMB)) return kFALSE;
   }
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
}
//_____________________________________________________________________________
void AliAnalysisTaskPhiFlow::InitializeBayesianPID(AliAODEvent* event)
{
   // Initialize the Bayesian PID object for AOD
   fBayesianResponse->SetDetResponse(event, fCentrality);
}
//_____________________________________________________________________________
template <typename T> Bool_t AliAnalysisTaskPhiFlow::PassesTPCbayesianCut(T* track) const
{
   // Check if the particle passes the TPC TOF bayesian cut.
   if ((!fStrictKaonCuts) && (!PassesStrictKaonCuts(track))) return kFALSE;
   fBayesianResponse->ComputeProb(track);
   if (!fBayesianResponse->GetCurrentMask(0)) return kFALSE; // return false if TPC has no response
   Float_t *probabilities = fBayesianResponse->GetProb();
   if (probabilities[3] > fParticleProbability)
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
   // see if track passes additional kaon selection cuts
   // with many thanks to francesco noferini
   Double_t b[2] = { -99., -99.};
   Double_t bCov[3] = { -99., -99., -99.};
   if (!track->PropagateToDCA(fESD->GetPrimaryVertex(), fESD->GetMagneticField(), 100., b, bCov)) return kFALSE;
   if ((TMath::Abs(b[0]) > 3.0) || (TMath::Abs(b[1]) > 2.4)) return kFALSE;
   return kTRUE;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskPhiFlow::PassesStrictKaonCuts(AliAODTrack* track) const
{
   // see if track passes additional kaon selection cuts
   // with many thanks to francesco noferini
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
   switch (fPIDtype)
   {
      case kTPC:
         // general tpc pid
         fNOPID->Fill(track->P(), track->GetTPCsignal());
         Double_t p[AliPID::kSPECIES];
         track->GetTPCpid(p);
         if(p[3] > fParticleProbability)
         {
            fPIDk->Fill(track->P(), track->GetTPCsignal());
            return kTRUE;
         }
         break;
      case kCombined:
         // Use null_cuts dummy to check if track passes bayesian TPC TOF pid cut
         fNOPID->Fill(track->P(), track->GetTPCsignal());
         if (PassesTPCbayesianCut(track))
         {
            fPIDk->Fill(track->P(), track->GetTPCsignal());
            return kTRUE;
         }
         break;
      default:
         AliFatal("No PID procedure set or available. Analysis will terminate!");
         return kFALSE;
         break;
   }
   return kFALSE;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskPhiFlow::IsKaon(AliAODTrack* track) const
{
   // Check if particle is a kaon according to method set in steering macro

   if (fRequireTPCStandAlone && (!track->TestFilterBit(1))) return kFALSE;
   switch (fPIDtype)
   {
      case kTPC:
         fNOPID->Fill(track->P(), track->GetTPCsignal());
//        Double_t p[AliPID::kSPECIES];
//         AliAODpidUtil* pid;
//         pid->MakeTPCPID(track, p);
//         if(p[3] > fParticleProbability)
//         {
//            fPIDk->Fill(track->P(), track->GetTPCsignal());
//            return kTRUE;
//         }
         return kFALSE;
         break;
      case kCombined:
         fNOPID->Fill(track->P(), track->GetTPCsignal());
         if (PassesTPCbayesianCut(track))
         {
            fPIDk->Fill(track->P(), track->GetTPCsignal());
            return kTRUE;
         }
         break;
      default:
         AliFatal("No PID procedure set or available. Analysis will terminate!");
         return kFALSE;
         break;
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
      if ((0.0 <= pt) && (0.3 > pt)) { fInvMNP03->Fill(InvariantMass(track1, track2)); fPtSpectra03->Fill(pt);}
      if ((0.3 <= pt) && (0.6 > pt)) { fInvMNP36->Fill(InvariantMass(track1, track2)); fPtSpectra36->Fill(pt); }
      if ((0.6 <= pt) && (0.9 > pt)) { fInvMNP69->Fill(InvariantMass(track1, track2)); fPtSpectra69->Fill(pt); }
      if ((0.9 <= pt) && (1.2 > pt)) { fInvMNP912->Fill(InvariantMass(track1, track2)); fPtSpectra912->Fill(pt); }
      if ((1.2 <= pt) && (1.5 > pt)) { fInvMNP1215->Fill(InvariantMass(track1, track2)); fPtSpectra1215->Fill(pt); }
      if ((1.5 <= pt) && (1.8 > pt)) { fInvMNP1518->Fill(InvariantMass(track1, track2)); fPtSpectra1518->Fill(pt); }
      if ((1.8 <= pt) && (2.1 > pt)) { fInvMNP1821->Fill(InvariantMass(track1, track2)); fPtSpectra1821->Fill(pt); }
      if ((2.1 <= pt) && (2.4 > pt)) { fInvMNP2124->Fill(InvariantMass(track1, track2)); fPtSpectra2124->Fill(pt); }
      if ((2.4 <= pt) && (2.7 > pt)) { fInvMNP2427->Fill(InvariantMass(track1, track2)); fPtSpectra2427->Fill(pt); }
      if ((2.7 <= pt) && (3.0 > pt)) { fInvMNP2730->Fill(InvariantMass(track1, track2)); fPtSpectra2730->Fill(pt); }
      if ((3.0 <= pt) && (3.5 > pt)) { fInvMNP3035->Fill(InvariantMass(track1, track2)); fPtSpectra3035->Fill(pt); }
      if ((3.5 <= pt) && (4.0 > pt)) { fInvMNP3540->Fill(InvariantMass(track1, track2)); fPtSpectra3540->Fill(pt); }
      if ((4.0 <= pt) && (4.5 > pt)) { fInvMNP4045->Fill(InvariantMass(track1, track2)); fPtSpectra4045->Fill(pt); }
      if ((4.5 <= pt) && (5.0 > pt)) { fInvMNP4550->Fill(InvariantMass(track1, track2)); fPtSpectra4550->Fill(pt); }
      if ((5.0 <= pt) && (5.5 > pt)) { fInvMNP5055->Fill(InvariantMass(track1, track2)); fPtSpectra5055->Fill(pt); }
      if ((5.5 <= pt) && (6.0 > pt)) { fInvMNP5560->Fill(InvariantMass(track1, track2)); fPtSpectra5560->Fill(pt); }
      if ((6.0 <= pt) && (6.5 > pt)) { fInvMNP6065->Fill(InvariantMass(track1, track2)); fPtSpectra6065->Fill(pt); }
      if ((6.5 <= pt) && (7.0 > pt)) { fInvMNP6570->Fill(InvariantMass(track1, track2)); fPtSpectra6570->Fill(pt); }
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
void AliAnalysisTaskPhiFlow::EventPlane(const AliESDEvent* event)
{
   // Calculate Q-vectors and event plane resolution (two different methods)
   fQy = 0;
   fQx = 0;
   Double_t ceventsin(0);
   Double_t ceventcos(0);
   Double_t deventsin(0);
   Double_t deventcos(0);
   Double_t eeventsin(0);
   Double_t eeventcos(0);
   Double_t feventsin(0);
   Double_t feventcos(0);
   Int_t nTracks = event->GetNumberOfTracks();
   Int_t tpcMultiplicity(0);
   for (Int_t i = 0; i < nTracks; i++)
   {
      AliESDtrack* track = event->GetTrack(i);
      if (!EventPlaneTrack(track)) continue;
      if (fSetEventPlanePtCut && (track->Pt() > fEventPlanePtCut)) continue;
      tpcMultiplicity++;
      Double_t r = gRandom->Uniform(0, 1);
      if (r < 0.5)
      {
         ceventsin += TMath::Sin(2 * track->Phi());
         ceventcos += TMath::Cos(2 * track->Phi());
      }
      if (r > 0.5)
      {
         deventsin += TMath::Sin(2 * track->Phi());
         deventcos += TMath::Cos(2 * track->Phi());
      }
      if (track->Eta() > 0)
      {
         eeventsin += TMath::Sin(2 * track->Phi());
         eeventcos += TMath::Cos(2 * track->Phi());
      }
      if (track->Eta() < 0)
      {
         feventsin += TMath::Sin(2 * track->Phi());
         feventcos += TMath::Cos(2 * track->Phi());
      }
      fQy += TMath::Sin(2 * track->Phi());
      fQx += TMath::Cos(2 * track->Phi());
   }
   fTPCM->Fill(tpcMultiplicity);
   EventPlaneResolution(kTRUE, ceventcos, ceventsin, deventcos, deventsin, fQx, fQy, eeventcos, eeventsin, feventcos, feventsin);
   EventPlaneResolution(kFALSE, ceventcos, ceventsin, deventcos, deventsin, fQx, fQy, eeventcos, eeventsin, feventcos, feventsin);
   fEventPlane->Fill((TMath::ATan2(fQy, fQx) / 2) + TMath::Pi() / 2.);
}
//_____________________________________________________________________________
void AliAnalysisTaskPhiFlow::EventPlane(const AliAODEvent* event)
{
   // Overloaded function - see comment at EventPlane(const AliESDEvent* event);
   fQy = 0;
   fQx = 0;
   Double_t ceventsin(0);
   Double_t ceventcos(0);
   Double_t deventsin(0);
   Double_t deventcos(0);
   Double_t eeventsin(0);
   Double_t eeventcos(0);
   Double_t feventsin(0);
   Double_t feventcos(0);
   Int_t nTracks = event->GetNumberOfTracks();
   Int_t tpcMultiplicity(0);
   for (Int_t i = 0; i < nTracks; i++)
   {
      AliAODTrack* track = event->GetTrack(i);
      if (!EventPlaneTrack(track)) continue;
      if (fSetEventPlanePtCut && (track->Pt() > fEventPlanePtCut)) continue;
      tpcMultiplicity++;
      Double_t r = gRandom->Uniform(0, 1);
      if (r < 0.5)
      {
         ceventsin += TMath::Sin(2 * track->Phi());
         ceventcos += TMath::Cos(2 * track->Phi());
      }
      if (r > 0.5)
      {
         deventsin += TMath::Sin(2 * track->Phi());
         deventcos += TMath::Cos(2 * track->Phi());
      }
      if (track->Eta() > 0)
      {
         eeventsin += TMath::Sin(2 * track->Phi());
         eeventcos += TMath::Cos(2 * track->Phi());
      }
      if (track->Eta() < 0)
      {
         feventsin += TMath::Sin(2 * track->Phi());
         feventcos += TMath::Cos(2 * track->Phi());
      }
      fQy += TMath::Sin(2 * track->Phi());
      fQx += TMath::Cos(2 * track->Phi());
   }
   fTPCM->Fill(tpcMultiplicity);
   EventPlaneResolution(kTRUE, ceventcos, ceventsin, deventcos, deventsin, fQx, fQy, eeventcos, eeventsin, feventcos, feventsin);
   EventPlaneResolution(kFALSE, ceventcos, ceventsin, deventcos, deventsin, fQx, fQy, eeventcos, eeventsin, feventcos, feventsin);
   fEventPlane->Fill((TMath::ATan2(fQy, fQx) / 2 + TMath::Pi() / 2.));
}
//_____________________________________________________________________________
void AliAnalysisTaskPhiFlow::EventPlaneResolution(Bool_t random,
                                                     Double_t cossumQ1,
                                                     Double_t sinsumQ1,
                                                     Double_t cossumQ2,
                                                     Double_t sinsumQ2,
                                                     Double_t cossum,
                                                     Double_t sinsum,
                                                     Double_t cossumQwest,
                                                     Double_t sinsumQwest,
                                                     Double_t cossumQeast,
                                                     Double_t sinsumQeast) const
{
   // Calculate the event plane resolution (called twice, for different methods)
   if (random)
   {
      Double_t res2Sub = EventPlaneStar(kTRUE, cossumQ1, sinsumQ1, cossumQ2, sinsumQ2, cossum, sinsum, cossumQwest, sinsumQwest, cossumQeast, sinsumQeast);
      if (res2Sub < -900) return;
      Double_t chiSub2 = Chi(res2Sub);
      Double_t res2 = ResolutionAsFunctionOfChi(TMath::Sqrt(2.) * chiSub2); // full event plane res.
      fEventPlaneResolutionRandom->Fill(1, res2);
   }
   if (!random)
   {
      Double_t res2Sub = EventPlaneStar(kFALSE, cossumQ1, sinsumQ1, cossumQ2, sinsumQ2, cossum, sinsum, cossumQwest, sinsumQwest, cossumQeast, sinsumQeast);
      if (res2Sub < -900) return;
      Double_t chiSub2 = Chi(res2Sub);
      Double_t res2 = ResolutionAsFunctionOfChi(TMath::Sqrt(2.) * chiSub2); // full event plane res.
      fEventPlaneResolutionEta->Fill(1, res2);
   }
}
//_____________________________________________________________________________
Double_t AliAnalysisTaskPhiFlow::EventPlaneStar(Bool_t random,
                                                   Double_t cossumQ1,
                                                   Double_t sinsumQ1,
                                                   Double_t cossumQ2,
                                                   Double_t sinsumQ2,
                                                   Double_t cossum,
                                                   Double_t sinsum,
                                                   Double_t cossumQwest,
                                                   Double_t sinsumQwest,
                                                   Double_t cossumQeast,
                                                   Double_t sinsumQeast) const
{
   // Calculate subevent resolution
   // This method also calculates the event plane
   // This implementation is copied from STAR analysis code
   TVector2 *q1 = new TVector2(cossumQ1, sinsumQ1);
   TVector2 *q2 = new TVector2(cossumQ2, sinsumQ2);
   TVector2 *qt = new TVector2(cossum, sinsum);
   TVector2 *qwest  = new TVector2(cossumQwest, sinsumQwest);
   TVector2 *qeast = new TVector2(cossumQeast, sinsumQeast);
   if (qt->Mod() == 0)
   {
      delete qt;
      delete q2;
      delete q1;
      delete qwest;
      delete qeast;
      return -999.; // make sure length > 0
   }
   Double_t twopi = TMath::TwoPi();

   Double_t rxnPlane, rxnPlane1, rxnPlane2;
   rxnPlane = qt->Phi() / 2.; //reaction plane angle ranges from 0 to pi
   rxnPlane1 = q1->Phi() / 2.;
   rxnPlane2 = q2->Phi() / 2.;
   if ((qt->Mod()) != 0.0)
   {
      Double_t  eventPlane = 0.5 * qt->Phi();
      if (random) fEventPlaneSTAR->Fill(eventPlane);
   }
   if (random)
   {
      if ((q1->Mod2()) != 0.0 && (q2->Mod2()) != 0.0)
      {
         Double_t eventDiff = 0.5 * q1->Phi() - 0.5 * q2->Phi();
         if (eventDiff < 0.0) eventDiff += twopi / 2.0;
         return TMath::Cos(2.*eventDiff);
      }
   }
   if (!random)
   {
      if ((qwest->Mod2()) != 0.0 && (qeast->Mod2()) != 0.0)
      {
         Double_t eventDiff = 0.5 * qwest->Phi() - 0.5 * qeast->Phi();
         if (eventDiff < 0.0) eventDiff += twopi / 2.0;
         return TMath::Cos(2.*eventDiff);
      }
   }
   return -999.;
}
//_____________________________________________________________________________
Double_t AliAnalysisTaskPhiFlow::Chi(Double_t res) const
{
   // Iterative algorithm to solve R(chi) for chi
   Double_t chi   = 2.0;
   Double_t delta = 1.0;
   for (Int_t i = 0; i < 15; i++)
   {
      chi   = (ResolutionAsFunctionOfChi(chi) < res) ? chi + delta : chi - delta;
      delta = delta / 2.;
   }
   return chi;
}
//_____________________________________________________________________________
Double_t AliAnalysisTaskPhiFlow::ResolutionAsFunctionOfChi(Double_t chi) const
{
   // Definition of R(chi)

   Double_t con = 0.626657;                   // sqrt(pi/2)/2
   Double_t arg = chi * chi / 4.;
   Double_t res = con * chi * exp(-arg) * (TMath::BesselI0(arg) + TMath::BesselI1(arg));
   return res;
}
//_____________________________________________________________________________
template <typename T> void AliAnalysisTaskPhiFlow::EllipticFlow(Double_t* const v2, const T* track1, const T* track2) const
{
   // 1) Calculate elliptic flow coefficient v_2 - cos terms
   // 2) Store v_2 in array for further analysis
   // 3) Fill phi - Psi histograms per p_t bin
   Double_t y = fQy - (TMath::Sin(2 * track1->Phi()) + TMath::Sin(2 * track2->Phi()));
   Double_t x = fQx - (TMath::Cos(2 * track1->Phi()) + TMath::Cos(2 * track2->Phi()));
   Double_t psi2 = (TMath::ATan2(y, x)) / 2.;
   Double_t psi2_sub(0);
   psi2 > 0 ? psi2_sub = psi2: psi2_sub = psi2 + TMath::Pi();
   Double_t pt = PhiPt(track1, track2);
   // to calculate the angle between the particles, firstly, construct the p_t vector of the phi
   TVector3 a(track1->Px(), track1->Py(), track1->Pz());
   TVector3 b(track2->Px(), track2->Py(), track2->Pz());
   TVector3 c = a + b;
   Double_t phi = c.Phi();
   Double_t invm  = InvariantMass(track1, track2);
   Double_t flow = TMath::Cos(2 * (phi - psi2));
   if ((0 <= pt) && (0.3 > pt))
   {
      v2[0] += flow;
      v2[18]++;
      fProfV2InvM03->Fill(invm, flow);
      fDeltaPhiPsiNP03->Fill(phi - psi2_sub);
   }
   if ((0.3 <= pt) && (0.6 > pt))
   {
      v2[1] += flow;
      v2[19]++;
      fProfV2InvM36->Fill(invm, flow);
      fDeltaPhiPsiNP36->Fill(phi - psi2_sub);
   }
   if ((0.6 <= pt) && (0.9 > pt))
   {
      v2[2] += flow;
      v2[20]++;
      fProfV2InvM69->Fill(invm, flow);
      fDeltaPhiPsiNP69->Fill(phi - psi2_sub);
   }
   if ((0.9 <= pt) && (1.2 > pt))
   {
      v2[3] += flow;
      v2[21]++;
      fProfV2InvM912->Fill(invm, flow);
      fDeltaPhiPsiNP912->Fill(phi - psi2_sub);
   }
   if ((1.2 <= pt) && (1.5 > pt))
   {
      v2[4] += flow;
      v2[22]++;
      fProfV2InvM1215->Fill(invm, flow);
      fDeltaPhiPsiNP1215->Fill(phi - psi2_sub);
   }
   if ((1.5 <= pt) && (1.8 > pt))
   {
      v2[5] += flow;
      v2[23]++;
      fProfV2InvM1518->Fill(invm, flow);
      fDeltaPhiPsiNP1518->Fill(phi - psi2_sub);
   }
   if ((1.8 <= pt) && (2.1 > pt))
   {
      v2[6] += flow;
      v2[24]++;
      fProfV2InvM1821->Fill(invm, flow);
      fDeltaPhiPsiNP1821->Fill(phi - psi2_sub);
   }
   if ((2.1 <= pt) && (2.4 > pt))
   {
      v2[7] += flow;
      v2[25]++;
      fProfV2InvM2124->Fill(invm, flow);
      fDeltaPhiPsiNP2124->Fill(phi - psi2_sub);
   }
   if ((2.4 <= pt) && (2.7 > pt))
   {
      v2[8] += flow;
      v2[26]++;
      fProfV2InvM2427->Fill(invm, flow);
      fDeltaPhiPsiNP2427->Fill(phi - psi2_sub);
   }
   if ((2.7 <= pt) && (3.0 > pt))
   {
      v2[9] += flow;
      v2[27]++;
      fProfV2InvM2730->Fill(invm, flow);
      fDeltaPhiPsiNP2730->Fill(phi - psi2_sub);
   }
   if ((3.0 <= pt) && (3.5 > pt))
   {
      v2[10] += flow;
      v2[28]++;
      fProfV2InvM3035->Fill(invm, flow);
      fDeltaPhiPsiNP3035->Fill(phi - psi2_sub);
   }
   if ((3.5 <= pt) && (4.0 > pt))
   {
      v2[11] += flow;
      v2[29]++;
      fProfV2InvM3540->Fill(invm, flow);
      fDeltaPhiPsiNP3540->Fill(phi - psi2_sub);
   }
   if ((4.0 <= pt) && (4.5 > pt))
   {
      v2[12] += flow;
      v2[30]++;
      fProfV2InvM4045->Fill(invm, flow);
      fDeltaPhiPsiNP4045->Fill(phi - psi2_sub);
   }
   if ((4.5 <= pt) && (5.0 > pt))
   {
      v2[13] += flow;
      v2[31]++;
      fProfV2InvM4550->Fill(invm, flow);
      fDeltaPhiPsiNP4550->Fill(phi - psi2_sub);
   }
   if ((5.0 <= pt) && (5.5 > pt))
   {
      v2[14] += flow;
      v2[32]++;
      fProfV2InvM5055->Fill(invm, flow);
      fDeltaPhiPsiNP5055->Fill(phi - psi2_sub);
   }
   if ((5.5 <= pt) && (6.0 > pt))
   {
      v2[15] += flow;
      v2[33]++;
      fProfV2InvM5560->Fill(invm, flow);
      fDeltaPhiPsiNP5560->Fill(phi - psi2_sub);
   }
   if ((6.0 <= pt) && (6.5 > pt))
   {
      v2[16] += flow;
      v2[34]++;
      fProfV2InvM6065->Fill(invm, flow);
      fDeltaPhiPsiNP6065->Fill(phi - psi2_sub);
   }
   if ((6.5 <= pt) && (7.0 > pt))
   {
      v2[17] += flow;
      v2[35]++;
      fProfV2InvM6570->Fill(invm, flow);
      fDeltaPhiPsiNP6570->Fill(phi - psi2_sub);
   }
}
//_____________________________________________________________________________
template <typename T> void AliAnalysisTaskPhiFlow::EllipticFlowSin(Double_t* const v2Sin, const T* track1, const T* track2) const
{
   // Calculate elliptic flow coefficient v_2 - sin terms
   // store v_2 in array for further analysis
   Double_t y = fQy - (TMath::Sin(2 * track1->Phi()) + TMath::Sin(2 * track2->Phi()));
   Double_t x = fQx - (TMath::Cos(2 * track1->Phi()) + TMath::Cos(2 * track2->Phi()));
   Double_t psi2 = (TMath::ATan2(y, x)) / 2.;
   Double_t pt = PhiPt(track1, track2);
   TVector3 a(track1->Px(), track1->Py(), track1->Pz());
   TVector3 b(track2->Px(), track2->Py(), track2->Pz());
   TVector3 c = a + b;
   Double_t phi = c.Phi();
   Double_t invm = InvariantMass(track1, track2);
   Double_t flow = TMath::Sin(2 * (phi - psi2));
   if ((0 <= pt) && (0.3 > pt))
   {
      v2Sin[0] += flow;
      v2Sin[18]++;
      fProfV2SinInvM03->Fill(invm, flow);
   }
   if ((0.3 <= pt) && (0.6 > pt))
   {
      v2Sin[1] += flow;
      v2Sin[19]++;
      fProfV2SinInvM36->Fill(invm, flow);
   }
   if ((0.6 <= pt) && (0.9 > pt))
   {
      v2Sin[2] += flow;
      v2Sin[20]++;
      fProfV2SinInvM69->Fill(invm, flow);
   }
   if ((0.9 <= pt) && (1.2 > pt))
   {
      v2Sin[3] += flow;
      v2Sin[21]++;
      fProfV2SinInvM912->Fill(invm, flow);
   }
   if ((1.2 <= pt) && (1.5 > pt))
   {
      v2Sin[4] += flow;
      v2Sin[22]++;
      fProfV2SinInvM1215->Fill(invm, flow);
   }
   if ((1.5 <= pt) && (1.8 > pt))
   {
      v2Sin[5] += flow;
      v2Sin[23]++;
      fProfV2SinInvM1518->Fill(invm, flow);
   }
   if ((1.8 <= pt) && (2.1 > pt))
   {
      v2Sin[6] += flow;
      v2Sin[24]++;
      fProfV2SinInvM1821->Fill(invm, flow);
   }
   if ((2.1 <= pt) && (2.4 > pt))
   {
      v2Sin[7] += flow;
      v2Sin[25]++;
      fProfV2SinInvM2124->Fill(invm, flow);
   }
   if ((2.4 <= pt) && (2.7 > pt))
   {
      v2Sin[8] += flow;
      v2Sin[26]++;
      fProfV2SinInvM2427->Fill(invm, flow);
   }
   if ((2.7 <= pt) && (3.0 > pt))
   {
      v2Sin[9] += flow;
      v2Sin[27]++;
      fProfV2SinInvM2730->Fill(invm, flow);
   }

   if ((3.0 <= pt) && (3.5 > pt))
   {
      v2Sin[10] += flow;
      v2Sin[28]++;
      fProfV2SinInvM3035->Fill(invm, flow);
   }
   if ((3.5 <= pt) && (4.0 > pt))
   {
      v2Sin[11] += flow;
      v2Sin[29]++;
      fProfV2SinInvM3540->Fill(invm, flow);
   }
   if ((4.0 <= pt) && (4.5 > pt))
   {
      v2Sin[12] += flow;
      v2Sin[30]++;
      fProfV2SinInvM4045->Fill(invm, flow);
   }
   if ((4.5 <= pt) && (5.0 > pt))
   {
      v2Sin[13] += flow;
      v2Sin[31]++;
      fProfV2SinInvM4550->Fill(invm, flow);
   }
   if ((5.0 <= pt) && (5.5 > pt))
   {
      v2Sin[14] += flow;
      v2Sin[32]++;
      fProfV2SinInvM5055->Fill(invm, flow);
   }
   if ((5.5 <= pt) && (6.0 > pt))
   {
      v2Sin[15] += flow;
      v2Sin[33]++;
      fProfV2SinInvM5560->Fill(invm, flow);
   }
   if ((6.0 <= pt) && (6.5 > pt))
   {
      v2Sin[16] += flow;
      v2Sin[34]++;
      fProfV2SinInvM6065->Fill(invm, flow);
   }
   if ((6.5 <= pt) && (7.0 > pt))
   {
      v2Sin[17] += flow;
      v2Sin[35]++;
      fProfV2SinInvM6570->Fill(invm, flow);
   }
}
//_____________________________________________________________________________
template <typename T> Bool_t AliAnalysisTaskPhiFlow::EventPlaneTrack(T* track) const
{
   // Check if track is suitable for event plane estimation
   if (!track) return kFALSE;
   return fCutsRP->IsSelected(track);
}
//_____________________________________________________________________________
template <typename T> Bool_t AliAnalysisTaskPhiFlow::PhiTrack(T* track) const
{
   // Check if track is suitable for phi flow analysis
   if (!track) return kFALSE;
   return fPOICuts->IsSelected(track);
}
//_____________________________________________________________________________
void AliAnalysisTaskPhiFlow::FlowFinish(Double_t* const v2) const
{
   // Push track averaged (event) flow into flow array - cos terms
   Double_t pt[] = {150, 450, 750, 1050, 1350, 1650, 1950, 2250, 2550, 2850, 3250, 3750, 4250, 4750, 5250, 5750, 6250, 6750};
   for (Int_t i = 0; i < 18; i++)
   {
      if (v2[i + 18] == 0) continue;
      fProfV2->Fill(pt[i] / 1000., (v2[i] / v2[i + 18]));
   }
}
//_____________________________________________________________________________
void AliAnalysisTaskPhiFlow::FlowFinishSin(Double_t* const v2Sin) const
{
   // Push track averaged (event) flow into flow array - sin terms
   Double_t pt[] = {150, 450, 750, 1050, 1350, 1650, 1950, 2250, 2550, 2850, 3250, 3750, 4250, 4750, 5250, 5750, 6250, 6750};
   for (Int_t i = 0; i < 18; i++)
   {
      if (v2Sin[i + 18] == 0) continue;
      fProfV2Sin->Fill(pt[i] / 1000., (v2Sin[i] / v2Sin[i + 18]));
   }
}
//_____________________________________________________________________________
template <typename T> void AliAnalysisTaskPhiFlow::SetNullCuts(T* event)
{
   // Set null cuts
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
   for (int r = 0; r != 30; ++r)
   {
      fFlowEvent[r]->ClearFast();
      fFlowEvent[r]->Fill(fCutsRP, fNullCuts);
      fFlowEvent[r]->SetReferenceMultiplicity(iMulti);
      fFlowEvent[r]->DefineDeadZone(0, 0, 0, 0);
      fFlowEvent[r]->TagSubeventsInEta(fEtaMinA, fEtaMaxA,
                                       fEtaMinB, fEtaMaxB);
   }
}
//_____________________________________________________________________________
void AliAnalysisTaskPhiFlow::UserExec(Option_t *)
{
   // UserExec: execute for each event. Commented where necessary
   // check for AOD data type

   fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
   if (fAOD)
   {
      fAODAnalysis = kTRUE;
      // Check whether event passes event cuts
      if (!EventCut(fAOD)) return;
      InitializeBayesianPID(fAOD);
      SetNullCuts(fAOD);
      PrepareFlowEvent(fAOD->GetNumberOfTracks());
      // Calculate event plane Q vectors and event plane resolution
      EventPlane(fAOD);
      fEventStats->Fill(0);
      Int_t unTracks = fAOD->GetNumberOfTracks();
      AliAODTrack* un[unTracks];
      AliAODTrack* up[unTracks];
      Int_t unp(0);
      Int_t unn(0);
      Double_t v2[36];
      Double_t v2Sin[36];
      // Flush flow arrays
      for (Int_t i = 0; i < 36; i++)
      {
         v2[i] = 0.;
         v2Sin[i] = 0.;
      }
      // Loop through tracks, check for species (Kaons), fill arrays according to charge
      for (Int_t iTracks = 0; iTracks < unTracks; iTracks++)
      {
         AliAODTrack* track = fAOD->GetTrack(iTracks);
         if (!PhiTrack(track)) continue;
         if (fStrictKaonCuts&&(!PassesStrictKaonCuts(track))) continue;
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
            EllipticFlow(v2, up[pTracks], un[nTracks]);
            EllipticFlowSin(v2Sin, up[pTracks], un[nTracks]);
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
            for (Int_t r = 0; r != 30; ++r)
               if ((mass >= fFlowBands[0][r]) && (mass < fFlowBands[1][r]))
               {
                  AliFlowCandidateTrack *sTrack = (AliFlowCandidateTrack*)
                                                  MakeTrack(mass, pt, phi, eta, 2, nIDs);
                  for (Int_t iDau = 0; iDau != 2; ++iDau)
                     for (Int_t iRPs = 0; iRPs != fFlowEvent[r]->NumberOfTracks(); ++iRPs)
                     {
                        AliFlowTrack *iRP = (AliFlowTrack*)(fFlowEvent[r]->GetTrack(iRPs));
                        if (!iRP->InRPSelection()) continue;
                        if (fabs(sTrack->GetIDDaughter(iDau)) == fabs(iRP->GetID()))
                        {
                           sTrack->SetDaughter(iDau, iRP);
                           iRP->SetForRPSelection(kFALSE);
                        }
                     }
                  fFlowEvent[r]->AddTrack(sTrack);
               }
         }
      }
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
      //Push event averaged track values into flow members
      FlowFinish(v2);
      FlowFinishSin(v2Sin);

      PostData(1, fOutputList);
      for (int m = 0; m != 30; ++m) PostData(2 + m, fFlowEvent[m]);
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
      EventPlane(fESD);
      fEventStats->Fill(0);
      Int_t unTracks = fESD->GetNumberOfTracks();
      AliESDtrack* un[unTracks];
      AliESDtrack* up[unTracks];
      Int_t unp(0);
      Int_t unn(0);
      Double_t v2[36];
      Double_t v2Sin[36];
      // Flush flow arrays
      for (Int_t i = 0; i < 36; i++)
      {
         v2[i] = 0.;
         v2Sin[i] = 0.;
      }
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
            EllipticFlow(v2, up[pTracks], un[nTracks]);
            EllipticFlowSin(v2Sin, up[pTracks], un[nTracks]);
            Double_t pt = PhiPt(up[pTracks], un[nTracks]);
            Double_t mass = InvariantMass(up[pTracks], un[nTracks]);
            TVector3 a(up[pTracks]->Px(), up[pTracks]->Py(), up[pTracks]->Pz());
            TVector3 b(un[nTracks]->Px(), un[nTracks]->Py(), up[pTracks]->Pz());
            TVector3 c = a + b;
            Double_t phi = c.Phi();
            Double_t eta = c.Eta();
            int nIDs[2];
            nIDs[0] = up[pTracks]->GetID();
            nIDs[1] = un[nTracks]->GetID();
            for (Int_t r = 0; r != 30; ++r)
               if ((mass >= fFlowBands[0][r]) && (mass < fFlowBands[1][r]))
               {
                  AliFlowCandidateTrack *sTrack = (AliFlowCandidateTrack*)
                                                  MakeTrack(mass, pt, phi, eta, 2, nIDs);
                  if (0) printf("   ᶫInjecting phi candidate on band %d \n", r);
                  for (Int_t iDau = 0; iDau != 2; ++iDau)
                     for (Int_t iRPs = 0; iRPs != fFlowEvent[r]->NumberOfTracks(); ++iRPs)
                     {
                        AliFlowTrack *iRP = (AliFlowTrack*)(fFlowEvent[r]->GetTrack(iRPs));
                        if (!iRP->InRPSelection()) continue;
                        if (fabs(sTrack->GetIDDaughter(iDau)) == fabs(iRP->GetID()))
                        {
                           sTrack->SetDaughter(iDau, iRP);
                           iRP->SetForRPSelection(kFALSE);
                           if (0) printf("    ᶫdaughter%d with fID %d was removed from this RP set\n", iDau, sTrack->GetIDDaughter(iDau));
                        }
                     }
                  fFlowEvent[r]->AddTrack(sTrack);
               }
         }
      }
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
      //Push event averaged track values into flow members
      FlowFinish(v2);
      FlowFinishSin(v2Sin);

      PostData(1, fOutputList);
      for (int m = 0; m != 30; ++m) PostData(2 + m, fFlowEvent[m]);
   }
}
//_____________________________________________________________________________
void AliAnalysisTaskPhiFlow::Terminate(Option_t *)
{
   // Terminate
}
//______________________________________________________________________________
AliFlowCandidateTrack*  AliAnalysisTaskPhiFlow::MakeTrack(Double_t mass,
      Double_t pt, Double_t phi, Double_t eta,
      Int_t nDau, Int_t iID[]) const
{
   // Consruct Flow Candidate Track from two selected candidates
   AliFlowCandidateTrack *sTrack = new AliFlowCandidateTrack();
   sTrack->SetMass(mass);
   sTrack->SetPt(pt);
   sTrack->SetPhi(phi);
   sTrack->SetEta(eta);
   for (Int_t iDau = 0; iDau != nDau; ++iDau)
      sTrack->AddDaughter(iID[iDau]);
   sTrack->SetForPOISelection(kTRUE);
   sTrack->SetForRPSelection(kFALSE);
   return sTrack;
}
//_____________________________________________________________________________

