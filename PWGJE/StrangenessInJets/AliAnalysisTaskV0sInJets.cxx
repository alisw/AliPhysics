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

// task for analysis of V0s (K0S, (anti-)Lambda) in charged jets
// Author: Vit Kucera (vit.kucera@cern.ch)

#include "TChain.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "THnSparse.h"
#include "TCanvas.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include <TDatabasePDG.h>
#include <TPDGCode.h>
#include "AliPIDResponse.h"
#include "AliInputEventHandler.h"
#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"
#include "TClonesArray.h"
//#include "AliEventInfoObject.cxx"
//#include "AliV0Object.cxx"
//#include "AliJetObject.cxx"
#include "TRandom3.h"

#include "AliAnalysisTaskV0sInJets.h"

ClassImp(AliAnalysisTaskV0sInJets)

// upper edges of centrality bins
//const Int_t AliAnalysisTaskV0sInJets::fgkiCentBinRanges[AliAnalysisTaskV0sInJets::fgkiNBinsCent] = {10, 30, 50, 80}; // Alice Zimmermann
//const Int_t AliAnalysisTaskV0sInJets::fgkiCentBinRanges[AliAnalysisTaskV0sInJets::fgkiNBinsCent] = {10, 20, 40, 60, 80}; // Vit Kucera, initial binning
//const Int_t AliAnalysisTaskV0sInJets::fgkiCentBinRanges[AliAnalysisTaskV0sInJets::fgkiNBinsCent] = {5, 10, 20, 40, 60, 80}; // Iouri Belikov, LF analysis
const Int_t AliAnalysisTaskV0sInJets::fgkiCentBinRanges[AliAnalysisTaskV0sInJets::fgkiNBinsCent] = {10}; // only central

// axis: pT of V0
const Double_t AliAnalysisTaskV0sInJets::fgkdBinsPtV0[2] = {0, 12};
const Int_t AliAnalysisTaskV0sInJets::fgkiNBinsPtV0 = sizeof(AliAnalysisTaskV0sInJets::fgkdBinsPtV0) / sizeof((AliAnalysisTaskV0sInJets::fgkdBinsPtV0)[0]) - 1;
const Int_t AliAnalysisTaskV0sInJets::fgkiNBinsPtV0Init = int(((AliAnalysisTaskV0sInJets::fgkdBinsPtV0)[AliAnalysisTaskV0sInJets::fgkiNBinsPtV0] - (AliAnalysisTaskV0sInJets::fgkdBinsPtV0)[0]) / 0.1); // bin width 0.1 GeV/c
// axis: pT of jets
const Double_t AliAnalysisTaskV0sInJets::fgkdBinsPtJet[2] = {0, 100};
const Int_t AliAnalysisTaskV0sInJets::fgkiNBinsPtJet = sizeof(AliAnalysisTaskV0sInJets::fgkdBinsPtJet) / sizeof(AliAnalysisTaskV0sInJets::fgkdBinsPtJet[0]) - 1;
const Int_t AliAnalysisTaskV0sInJets::fgkiNBinsPtJetInit = int(((AliAnalysisTaskV0sInJets::fgkdBinsPtJet)[AliAnalysisTaskV0sInJets::fgkiNBinsPtJet] - (AliAnalysisTaskV0sInJets::fgkdBinsPtJet)[0]) / 5.); // bin width 5 GeV/c
// axis: K0S invariant mass
const Int_t AliAnalysisTaskV0sInJets::fgkiNBinsMassK0s = 300;
const Double_t AliAnalysisTaskV0sInJets::fgkdMassK0sMin = 0.35;
const Double_t AliAnalysisTaskV0sInJets::fgkdMassK0sMax = 0.65;
// axis: Lambda invariant mass
const Int_t AliAnalysisTaskV0sInJets::fgkiNBinsMassLambda = 200;
const Double_t AliAnalysisTaskV0sInJets::fgkdMassLambdaMin = 1.05;
const Double_t AliAnalysisTaskV0sInJets::fgkdMassLambdaMax = 1.25;


// Default constructor
AliAnalysisTaskV0sInJets::AliAnalysisTaskV0sInJets():
  AliAnalysisTaskSE(),
  fAODIn(0),
  fAODOut(0),
  fOutputListStd(0),
  fOutputListQA(0),
  fOutputListCuts(0),
  fOutputListMC(0),
//  ftreeOut(0),

  fiAODAnalysis(1),
  fbIsPbPb(1),

  fdCutDCAToPrimVtxMin(0.1),
  fdCutDCADaughtersMax(1.),
  fdCutNSigmadEdxMax(3),
  fdCutCPAMin(0.998),
  fdCutNTauMax(5),

  fsJetBranchName(0),
  fsJetBgBranchName(0),
  fdCutPtJetMin(0),
  fdCutPtTrackMin(5),
  fdRadiusJet(0.4),
  fdRadiusJetBg(0.4),
  fbJetSelection(0),
  fbMCAnalysis(0),
//  fbTreeOutput(0),
  fRandom(0),

  fdCutVertexZ(10),
  fdCutVertexR2(1),
  fdCutCentLow(0),
  fdCutCentHigh(80),

/*
  fBranchV0Rec(0),
  fBranchV0Gen(0),
  fBranchJet(0),
  fEventInfo(0),
*/
  fdCentrality(0),
  fh1EventCounterCut(0),
  fh1EventCent(0),
  fh1EventCent2(0),
  fh1EventCent2Jets(0),
  fh1EventCent2NoJets(0),
  fh2EventCentTracks(0),
  fh1V0CandPerEvent(0),
  fh1NRndConeCent(0),
  fh1NMedConeCent(0),
  fh1AreaExcluded(0),

  fh2CCK0s(0),
  fh2CCLambda(0),
  fh3CCMassCorrelBoth(0),
  fh3CCMassCorrelKNotL(0),
  fh3CCMassCorrelLNotK(0)
{
  for(Int_t i = 0; i < fgkiNQAIndeces; i++)
  {
    fh1QAV0Status[i] = 0;
    fh1QAV0TPCRefit[i] = 0;
    fh1QAV0TPCRows[i] = 0;
    fh1QAV0TPCFindable[i] = 0;
    fh1QAV0TPCRowsFind[i] = 0;
    fh1QAV0Eta[i] = 0;
    fh2QAV0EtaRows[i] = 0;
    fh2QAV0PtRows[i] = 0;
    fh2QAV0PhiRows[i] = 0;
    fh2QAV0NClRows[i] = 0;
    fh2QAV0EtaNCl[i] = 0;

    fh2QAV0EtaPtK0sPeak[i] = 0;
    fh2QAV0EtaEtaK0s[i] = 0;
    fh2QAV0PhiPhiK0s[i] = 0;
    fh1QAV0RapK0s[i] = 0;
    fh2QAV0PtPtK0sPeak[i] = 0;
    fh2ArmPodK0s[i] = 0;

    fh2QAV0EtaPtLambdaPeak[i] = 0;
    fh2QAV0EtaEtaLambda[i] = 0;
    fh2QAV0PhiPhiLambda[i] = 0;
    fh1QAV0RapLambda[i] = 0;
    fh2QAV0PtPtLambdaPeak[i] = 0;
    fh2ArmPodLambda[i] = 0;

    fh2QAV0EtaPtALambdaPeak[i] = 0;
    fh2QAV0EtaEtaALambda[i] = 0;
    fh2QAV0PhiPhiALambda[i] = 0;
    fh1QAV0RapALambda[i] = 0;
    fh2QAV0PtPtALambdaPeak[i] = 0;
    fh2ArmPodALambda[i] = 0;

    fh1QAV0Pt[i] = 0;
    fh1QAV0Charge[i] = 0;
    fh1QAV0DCAVtx[i] = 0;
    fh1QAV0DCAV0[i] = 0;
    fh1QAV0Cos[i] = 0;
    fh1QAV0R[i] = 0;
    fh1QACTau2D[i] = 0;
    fh1QACTau3D[i] = 0;

    fh2ArmPod[i] = 0;

    /*
    fh2CutTPCRowsK0s[i] = 0;
    fh2CutTPCRowsLambda[i] = 0;
    fh2CutPtPosK0s[i] = 0;
    fh2CutPtNegK0s[i] = 0;
    fh2CutPtPosLambda[i] = 0;
    fh2CutPtNegLambda[i] = 0;
    fh2CutDCAVtx[i] = 0;
    fh2CutDCAV0[i] = 0;
    fh2CutCos[i] = 0;
    fh2CutR[i] = 0;
    fh2CutEtaK0s[i] = 0;
    fh2CutEtaLambda[i] = 0;
    fh2CutRapK0s[i] = 0;
    fh2CutRapLambda[i] = 0;
    fh2CutCTauK0s[i] = 0;
    fh2CutCTauLambda[i] = 0;
    fh2CutPIDPosK0s[i] = 0;
    fh2CutPIDNegK0s[i] = 0;
    fh2CutPIDPosLambda[i] = 0;
    fh2CutPIDNegLambda[i] = 0;

    fh2Tau3DVs2D[i] = 0;
    */
  }
  for(Int_t i = 0; i < fgkiNCategV0; i++)
  {
    fh1V0InvMassK0sAll[i] = 0;
    fh1V0InvMassLambdaAll[i] = 0;
    fh1V0InvMassALambdaAll[i] = 0;
  }
  for(Int_t i = 0; i < fgkiNBinsCent; i++)
  {
    fh1EventCounterCutCent[i] = 0;
    fh1V0CounterCentK0s[i] = 0;
    fh1V0CounterCentLambda[i] = 0;
    fh1V0CounterCentALambda[i] = 0;
    fh1V0CandPerEventCentK0s[i] = 0;
    fh1V0CandPerEventCentLambda[i] = 0;
    fh1V0CandPerEventCentALambda[i] = 0;
    fh1V0InvMassK0sCent[i] = 0;
    fh1V0InvMassLambdaCent[i] = 0;
    fh1V0InvMassALambdaCent[i] = 0;
    fh1V0K0sPtMCGen[i] = 0;
    fh2V0K0sPtMassMCRec[i] = 0;
    fh1V0K0sPtMCRecFalse[i] = 0;
    fh2V0K0sEtaPtMCGen[i] = 0;
    fh3V0K0sEtaPtMassMCRec[i] = 0;
    fh2V0K0sInJetPtMCGen[i] = 0;
    fh3V0K0sInJetPtMassMCRec[i] = 0;
    fh3V0K0sInJetEtaPtMCGen[i] = 0;
    fh4V0K0sInJetEtaPtMassMCRec[i] = 0;
    fh2V0K0sMCResolMPt[i] = 0;
    fh2V0K0sMCPtGenPtRec[i] = 0;
    fh1V0LambdaPtMCGen[i] = 0;
    fh2V0LambdaPtMassMCRec[i] = 0;
    fh1V0LambdaPtMCRecFalse[i] = 0;
    fh2V0LambdaEtaPtMCGen[i] = 0;
    fh3V0LambdaEtaPtMassMCRec[i] = 0;
    fh2V0LambdaInJetPtMCGen[i] = 0;
    fh3V0LambdaInJetPtMassMCRec[i] = 0;
    fh3V0LambdaInJetEtaPtMCGen[i] = 0;
    fh4V0LambdaInJetEtaPtMassMCRec[i] = 0;
    fh2V0LambdaMCResolMPt[i] = 0;
    fh2V0LambdaMCPtGenPtRec[i] = 0;
    fhnV0LambdaInclMCFD[i] = 0;
    fhnV0LambdaInJetsMCFD[i] = 0;
    fhnV0LambdaBulkMCFD[i] = 0;
    fh1V0XiPtMCGen[i] = 0;
    fh1V0ALambdaPt[i] = 0;
    fh1V0ALambdaPtMCGen[i] = 0;
    fh1V0ALambdaPtMCRec[i] = 0;
    fh2V0ALambdaPtMassMCRec[i] = 0;
    fh1V0ALambdaPtMCRecFalse[i] = 0;
    fh2V0ALambdaEtaPtMCGen[i] = 0;
    fh3V0ALambdaEtaPtMassMCRec[i] = 0;
    fh2V0ALambdaInJetPtMCGen[i] = 0;
    fh2V0ALambdaInJetPtMCRec[i] = 0;
    fh3V0ALambdaInJetPtMassMCRec[i] = 0;
    fh3V0ALambdaInJetEtaPtMCGen[i] = 0;
    fh4V0ALambdaInJetEtaPtMassMCRec[i] = 0;
    fh2V0ALambdaMCResolMPt[i] = 0;
    fh2V0ALambdaMCPtGenPtRec[i] = 0;
    fhnV0ALambdaInclMCFD[i] = 0;
    fhnV0ALambdaInJetsMCFD[i] = 0;
    fhnV0ALambdaBulkMCFD[i] = 0;
    fh1V0AXiPtMCGen[i] = 0;

    // eta daughters
//      fhnV0K0sInclDaughterEtaPtPtMCGen[i] = 0;
    fhnV0K0sInclDaughterEtaPtPtMCRec[i] = 0;
//      fhnV0K0sInJetsDaughterEtaPtPtMCGen[i] = 0;
    fhnV0K0sInJetsDaughterEtaPtPtMCRec[i] = 0;
//      fhnV0LambdaInclDaughterEtaPtPtMCGen[i] = 0;
    fhnV0LambdaInclDaughterEtaPtPtMCRec[i] = 0;
//      fhnV0LambdaInJetsDaughterEtaPtPtMCGen[i] = 0;
    fhnV0LambdaInJetsDaughterEtaPtPtMCRec[i] = 0;
//      fhnV0ALambdaInclDaughterEtaPtPtMCGen[i] = 0;
    fhnV0ALambdaInclDaughterEtaPtPtMCRec[i] = 0;
//      fhnV0ALambdaInJetsDaughterEtaPtPtMCGen[i] = 0;
    fhnV0ALambdaInJetsDaughterEtaPtPtMCRec[i] = 0;

    // Inclusive
    fhnV0InclusiveK0s[i] = 0;
    fhnV0InclusiveLambda[i] = 0;
    fhnV0InclusiveALambda[i] = 0;
    // Cones
    fhnV0InJetK0s[i] = 0;
    fhnV0InPerpK0s[i] = 0;
    fhnV0InRndK0s[i] = 0;
    fhnV0InMedK0s[i] = 0;
    fhnV0OutJetK0s[i] = 0;
    fhnV0NoJetK0s[i] = 0;
    fhnV0InJetLambda[i] = 0;
    fhnV0InPerpLambda[i] = 0;
    fhnV0InRndLambda[i] = 0;
    fhnV0InMedLambda[i] = 0;
    fhnV0OutJetLambda[i] = 0;
    fhnV0NoJetLambda[i] = 0;
    fhnV0InJetALambda[i] = 0;
    fhnV0InPerpALambda[i] = 0;
    fhnV0InRndALambda[i] = 0;
    fhnV0InMedALambda[i] = 0;
    fhnV0OutJetALambda[i] = 0;
    fhnV0NoJetALambda[i] = 0;

    fh2V0PtJetAngleK0s[i] = 0;
    fh2V0PtJetAngleLambda[i] = 0;
    fh2V0PtJetAngleALambda[i] = 0;
    fh1DCAInK0s[i] = 0;
    fh1DCAInLambda[i] = 0;
    fh1DCAInALambda[i] = 0;
    fh1DCAOutK0s[i] = 0;
    fh1DCAOutLambda[i] = 0;
    fh1DCAOutALambda[i] = 0;
//    fh1DeltaZK0s[i] = 0;
//    fh1DeltaZLambda[i] = 0;
//    fh1DeltaZALambda[i] = 0;

    fh1PtJet[i] = 0;
    fh1EtaJet[i] = 0;
    fh2EtaPtJet[i] = 0;
    fh1PhiJet[i] = 0;
    fh1NJetPerEvent[i] = 0;
    fh2EtaPhiRndCone[i] = 0;
    fh2EtaPhiMedCone[i] = 0;

    fh1VtxZ[i] = 0;
    fh2VtxXY[i] = 0;
  }
}

// Constructor
AliAnalysisTaskV0sInJets::AliAnalysisTaskV0sInJets(const char* name):
  AliAnalysisTaskSE(name),
  fAODIn(0),
  fAODOut(0),
  fOutputListStd(0),
  fOutputListQA(0),
  fOutputListCuts(0),
  fOutputListMC(0),
//  ftreeOut(0),

  fiAODAnalysis(1),
  fbIsPbPb(1),

  fdCutDCAToPrimVtxMin(0.1),
  fdCutDCADaughtersMax(1.),
  fdCutNSigmadEdxMax(3),
  fdCutCPAMin(0.998),
  fdCutNTauMax(5),

  fsJetBranchName(0),
  fsJetBgBranchName(0),
  fdCutPtJetMin(0),
  fdCutPtTrackMin(5),
  fdRadiusJet(0.4),
  fdRadiusJetBg(0.4),
  fbJetSelection(0),
  fbMCAnalysis(0),
//  fbTreeOutput(0),
  fRandom(0),

  fdCutVertexZ(10),
  fdCutVertexR2(1),
  fdCutCentLow(0),
  fdCutCentHigh(80),
/*
  fBranchV0Rec(0),
  fBranchV0Gen(0),
  fBranchJet(0),
  fEventInfo(0),
*/
  fdCentrality(0),
  fh1EventCounterCut(0),
  fh1EventCent(0),
  fh1EventCent2(0),
  fh1EventCent2Jets(0),
  fh1EventCent2NoJets(0),
  fh2EventCentTracks(0),
  fh1V0CandPerEvent(0),
  fh1NRndConeCent(0),
  fh1NMedConeCent(0),
  fh1AreaExcluded(0),

  fh2CCK0s(0),
  fh2CCLambda(0),
  fh3CCMassCorrelBoth(0),
  fh3CCMassCorrelKNotL(0),
  fh3CCMassCorrelLNotK(0)
{
  for(Int_t i = 0; i < fgkiNQAIndeces; i++)
  {
    fh1QAV0Status[i] = 0;
    fh1QAV0TPCRefit[i] = 0;
    fh1QAV0TPCRows[i] = 0;
    fh1QAV0TPCFindable[i] = 0;
    fh1QAV0TPCRowsFind[i] = 0;
    fh1QAV0Eta[i] = 0;
    fh2QAV0EtaRows[i] = 0;
    fh2QAV0PtRows[i] = 0;
    fh2QAV0PhiRows[i] = 0;
    fh2QAV0NClRows[i] = 0;
    fh2QAV0EtaNCl[i] = 0;

    fh2QAV0EtaPtK0sPeak[i] = 0;
    fh2QAV0EtaEtaK0s[i] = 0;
    fh2QAV0PhiPhiK0s[i] = 0;
    fh1QAV0RapK0s[i] = 0;
    fh2QAV0PtPtK0sPeak[i] = 0;
    fh2ArmPodK0s[i] = 0;

    fh2QAV0EtaPtLambdaPeak[i] = 0;
    fh2QAV0EtaEtaLambda[i] = 0;
    fh2QAV0PhiPhiLambda[i] = 0;
    fh1QAV0RapLambda[i] = 0;
    fh2QAV0PtPtLambdaPeak[i] = 0;
    fh2ArmPodLambda[i] = 0;

    fh2QAV0EtaPtALambdaPeak[i] = 0;
    fh2QAV0EtaEtaALambda[i] = 0;
    fh2QAV0PhiPhiALambda[i] = 0;
    fh1QAV0RapALambda[i] = 0;
    fh2QAV0PtPtALambdaPeak[i] = 0;
    fh2ArmPodALambda[i] = 0;

    fh1QAV0Pt[i] = 0;
    fh1QAV0Charge[i] = 0;
    fh1QAV0DCAVtx[i] = 0;
    fh1QAV0DCAV0[i] = 0;
    fh1QAV0Cos[i] = 0;
    fh1QAV0R[i] = 0;
    fh1QACTau2D[i] = 0;
    fh1QACTau3D[i] = 0;

    fh2ArmPod[i] = 0;

    /*
    fh2CutTPCRowsK0s[i] = 0;
    fh2CutTPCRowsLambda[i] = 0;
    fh2CutPtPosK0s[i] = 0;
    fh2CutPtNegK0s[i] = 0;
    fh2CutPtPosLambda[i] = 0;
    fh2CutPtNegLambda[i] = 0;
    fh2CutDCAVtx[i] = 0;
    fh2CutDCAV0[i] = 0;
    fh2CutCos[i] = 0;
    fh2CutR[i] = 0;
    fh2CutEtaK0s[i] = 0;
    fh2CutEtaLambda[i] = 0;
    fh2CutRapK0s[i] = 0;
    fh2CutRapLambda[i] = 0;
    fh2CutCTauK0s[i] = 0;
    fh2CutCTauLambda[i] = 0;
    fh2CutPIDPosK0s[i] = 0;
    fh2CutPIDNegK0s[i] = 0;
    fh2CutPIDPosLambda[i] = 0;
    fh2CutPIDNegLambda[i] = 0;

    fh2Tau3DVs2D[i] = 0;
    */
  }
  for(Int_t i = 0; i < fgkiNCategV0; i++)
  {
    fh1V0InvMassK0sAll[i] = 0;
    fh1V0InvMassLambdaAll[i] = 0;
    fh1V0InvMassALambdaAll[i] = 0;
  }
  for(Int_t i = 0; i < fgkiNBinsCent; i++)
  {
    fh1EventCounterCutCent[i] = 0;
    fh1V0CounterCentK0s[i] = 0;
    fh1V0CounterCentLambda[i] = 0;
    fh1V0CounterCentALambda[i] = 0;
    fh1V0CandPerEventCentK0s[i] = 0;
    fh1V0CandPerEventCentLambda[i] = 0;
    fh1V0CandPerEventCentALambda[i] = 0;
    fh1V0InvMassK0sCent[i] = 0;
    fh1V0InvMassLambdaCent[i] = 0;
    fh1V0InvMassALambdaCent[i] = 0;
    fh1V0K0sPtMCGen[i] = 0;
    fh2V0K0sPtMassMCRec[i] = 0;
    fh1V0K0sPtMCRecFalse[i] = 0;
    fh2V0K0sEtaPtMCGen[i] = 0;
    fh3V0K0sEtaPtMassMCRec[i] = 0;
    fh2V0K0sInJetPtMCGen[i] = 0;
    fh3V0K0sInJetPtMassMCRec[i] = 0;
    fh3V0K0sInJetEtaPtMCGen[i] = 0;
    fh4V0K0sInJetEtaPtMassMCRec[i] = 0;
    fh2V0K0sMCResolMPt[i] = 0;
    fh2V0K0sMCPtGenPtRec[i] = 0;
    fh1V0LambdaPtMCGen[i] = 0;
    fh2V0LambdaPtMassMCRec[i] = 0;
    fh1V0LambdaPtMCRecFalse[i] = 0;
    fh2V0LambdaEtaPtMCGen[i] = 0;
    fh3V0LambdaEtaPtMassMCRec[i] = 0;
    fh2V0LambdaInJetPtMCGen[i] = 0;
    fh3V0LambdaInJetPtMassMCRec[i] = 0;
    fh3V0LambdaInJetEtaPtMCGen[i] = 0;
    fh4V0LambdaInJetEtaPtMassMCRec[i] = 0;
    fh2V0LambdaMCResolMPt[i] = 0;
    fh2V0LambdaMCPtGenPtRec[i] = 0;
    fhnV0LambdaInclMCFD[i] = 0;
    fhnV0LambdaInJetsMCFD[i] = 0;
    fhnV0LambdaBulkMCFD[i] = 0;
    fh1V0XiPtMCGen[i] = 0;
    fh1V0ALambdaPt[i] = 0;
    fh1V0ALambdaPtMCGen[i] = 0;
    fh1V0ALambdaPtMCRec[i] = 0;
    fh2V0ALambdaPtMassMCRec[i] = 0;
    fh1V0ALambdaPtMCRecFalse[i] = 0;
    fh2V0ALambdaEtaPtMCGen[i] = 0;
    fh3V0ALambdaEtaPtMassMCRec[i] = 0;
    fh2V0ALambdaInJetPtMCGen[i] = 0;
    fh2V0ALambdaInJetPtMCRec[i] = 0;
    fh3V0ALambdaInJetPtMassMCRec[i] = 0;
    fh3V0ALambdaInJetEtaPtMCGen[i] = 0;
    fh4V0ALambdaInJetEtaPtMassMCRec[i] = 0;
    fh2V0ALambdaMCResolMPt[i] = 0;
    fh2V0ALambdaMCPtGenPtRec[i] = 0;
    fhnV0ALambdaInclMCFD[i] = 0;
    fhnV0ALambdaInJetsMCFD[i] = 0;
    fhnV0ALambdaBulkMCFD[i] = 0;
    fh1V0AXiPtMCGen[i] = 0;

    // eta daughters
//      fhnV0K0sInclDaughterEtaPtPtMCGen[i] = 0;
    fhnV0K0sInclDaughterEtaPtPtMCRec[i] = 0;
//      fhnV0K0sInJetsDaughterEtaPtPtMCGen[i] = 0;
    fhnV0K0sInJetsDaughterEtaPtPtMCRec[i] = 0;
//      fhnV0LambdaInclDaughterEtaPtPtMCGen[i] = 0;
    fhnV0LambdaInclDaughterEtaPtPtMCRec[i] = 0;
//      fhnV0LambdaInJetsDaughterEtaPtPtMCGen[i] = 0;
    fhnV0LambdaInJetsDaughterEtaPtPtMCRec[i] = 0;
//      fhnV0ALambdaInclDaughterEtaPtPtMCGen[i] = 0;
    fhnV0ALambdaInclDaughterEtaPtPtMCRec[i] = 0;
//      fhnV0ALambdaInJetsDaughterEtaPtPtMCGen[i] = 0;
    fhnV0ALambdaInJetsDaughterEtaPtPtMCRec[i] = 0;

    // Inclusive
    fhnV0InclusiveK0s[i] = 0;
    fhnV0InclusiveLambda[i] = 0;
    fhnV0InclusiveALambda[i] = 0;
    // Cones
    fhnV0InJetK0s[i] = 0;
    fhnV0InPerpK0s[i] = 0;
    fhnV0InRndK0s[i] = 0;
    fhnV0InMedK0s[i] = 0;
    fhnV0OutJetK0s[i] = 0;
    fhnV0NoJetK0s[i] = 0;
    fhnV0InJetLambda[i] = 0;
    fhnV0InPerpLambda[i] = 0;
    fhnV0InRndLambda[i] = 0;
    fhnV0InMedLambda[i] = 0;
    fhnV0OutJetLambda[i] = 0;
    fhnV0NoJetLambda[i] = 0;
    fhnV0InJetALambda[i] = 0;
    fhnV0InPerpALambda[i] = 0;
    fhnV0InRndALambda[i] = 0;
    fhnV0InMedALambda[i] = 0;
    fhnV0OutJetALambda[i] = 0;
    fhnV0NoJetALambda[i] = 0;

    fh2V0PtJetAngleK0s[i] = 0;
    fh2V0PtJetAngleLambda[i] = 0;
    fh2V0PtJetAngleALambda[i] = 0;
    fh1DCAInK0s[i] = 0;
    fh1DCAInLambda[i] = 0;
    fh1DCAInALambda[i] = 0;
    fh1DCAOutK0s[i] = 0;
    fh1DCAOutLambda[i] = 0;
    fh1DCAOutALambda[i] = 0;
//    fh1DeltaZK0s[i] = 0;
//    fh1DeltaZLambda[i] = 0;
//    fh1DeltaZALambda[i] = 0;

    fh1PtJet[i] = 0;
    fh1EtaJet[i] = 0;
    fh2EtaPtJet[i] = 0;
    fh1PhiJet[i] = 0;
    fh1NJetPerEvent[i] = 0;
    fh2EtaPhiRndCone[i] = 0;
    fh2EtaPhiMedCone[i] = 0;

    fh1VtxZ[i] = 0;
    fh2VtxXY[i] = 0;
  }
  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 id reserved by the base class for AOD
  // Output slot #1 writes into a TList container
  DefineOutput(1, TList::Class());
  DefineOutput(2, TList::Class());
  DefineOutput(3, TList::Class());
  DefineOutput(4, TList::Class());
  DefineOutput(5, TTree::Class());
}

AliAnalysisTaskV0sInJets::~AliAnalysisTaskV0sInJets()
{
/*
  if (fBranchV0Rec)
    fBranchV0Rec->Delete();
  delete fBranchV0Rec;
  fBranchV0Rec = 0;
  if (fBranchV0Gen)
    fBranchV0Gen->Delete();
  delete fBranchV0Gen;
  fBranchV0Gen = 0;
  if (fBranchJet)
    fBranchJet->Delete();
  delete fBranchJet;
  fBranchJet = 0;
  if (fEventInfo)
    fEventInfo->Delete();
  delete fEventInfo;
  fEventInfo = 0;
*/
  delete fRandom;
  fRandom = 0;
}

void AliAnalysisTaskV0sInJets::UserCreateOutputObjects()
{
  // Create histograms
  // Called once

  fRandom = new TRandom3(0);

/*
  if (!fBranchV0Rec && fbTreeOutput)
    {
//      fBranchV0Rec = new TClonesArray("AliAODv0",0);
      fBranchV0Rec = new TClonesArray("AliV0Object",0);
      fBranchV0Rec->SetName("branch_V0Rec");
    }
  if (!fBranchV0Gen && fbTreeOutput)
    {
      fBranchV0Gen = new TClonesArray("AliAODMCParticle",0);
      fBranchV0Gen->SetName("branch_V0Gen");
    }
  if (!fBranchJet && fbTreeOutput)
    {
//      fBranchJet = new TClonesArray("AliAODJet",0);
      fBranchJet = new TClonesArray("AliJetObject",0);
      fBranchJet->SetName("branch_Jet");
    }
  if (!fEventInfo && fbTreeOutput)
    {
      fEventInfo = new AliEventInfoObject();
      fEventInfo->SetName("eventInfo");
    }
  Int_t dSizeBuffer = 32000; // default 32000
  if (fbTreeOutput)
  {
  ftreeOut = new TTree("treeV0","Tree V0");
  ftreeOut->Branch("branch_V0Rec",&fBranchV0Rec,dSizeBuffer,2);
  ftreeOut->Branch("branch_V0Gen",&fBranchV0Gen,dSizeBuffer,2);
  ftreeOut->Branch("branch_Jet",&fBranchJet,dSizeBuffer,2);
  ftreeOut->Branch("eventInfo",&fEventInfo,dSizeBuffer,2);
  }
*/

  fOutputListStd = new TList();
  fOutputListStd->SetOwner();
  fOutputListQA = new TList();
  fOutputListQA->SetOwner();
  fOutputListCuts = new TList();
  fOutputListCuts->SetOwner();
  fOutputListMC = new TList();
  fOutputListMC->SetOwner();

  // event categories
  const Int_t iNCategEvent = 6;
  TString categEvent[iNCategEvent] = {"coll. candid.", "AOD OK", "vtx & cent", "with V0", "with jets", "jet selection"};
  // labels for stages of V0 selection
  TString categV0[fgkiNCategV0] = {"all"/*0*/, "mass range"/*1*/, "rec. method"/*2*/, "tracks TPC"/*3*/, "track pt"/*4*/, "DCA prim v"/*5*/, "DCA daughters"/*6*/, "CPA"/*7*/, "volume"/*8*/, "track #it{#eta}"/*9*/, "V0 #it{y} & #it{#eta}"/*10*/, "lifetime"/*11*/, "PID"/*12*/, "Arm.-Pod."/*13*/, "inclusive"/*14*/, "in jet event"/*15*/, "in jet"/*16*/};

  fh1EventCounterCut = new TH1D("fh1EventCounterCut", "Number of events after filtering;selection filter;counts", iNCategEvent, 0, iNCategEvent);
  for(Int_t i = 0; i < iNCategEvent; i++)
    fh1EventCounterCut->GetXaxis()->SetBinLabel(i + 1, categEvent[i].Data());
  fh1EventCent2 = new TH1D("fh1EventCent2", "Number of events vs centrality;centrality;counts", 100, 0, 100);
  fh1EventCent2Jets = new TH1D("fh1EventCent2Jets", "Number of sel.-jet events vs centrality;centrality;counts", 100, 0, 100);
  fh1EventCent2NoJets = new TH1D("fh1EventCent2NoJets", "Number of no-jet events vs centrality;centrality;counts", 100, 0, 100);
  fh2EventCentTracks = new TH2D("fh2EventCentTracks", "Number of tracks vs centrality;centrality;tracks;counts", 100, 0, 100, 150, 0, 15e3);
  fh1EventCent = new TH1D("fh1EventCent", "Number of events in centrality bins;centrality;counts", fgkiNBinsCent, 0, fgkiNBinsCent);
  for(Int_t i = 0; i < fgkiNBinsCent; i++)
    fh1EventCent->GetXaxis()->SetBinLabel(i + 1, GetCentBinLabel(i).Data());
  fh1NRndConeCent = new TH1D("fh1NRndConeCent", "Number of rnd. cones in centrality bins;centrality;counts", fgkiNBinsCent, 0, fgkiNBinsCent);
  for(Int_t i = 0; i < fgkiNBinsCent; i++)
    fh1NRndConeCent->GetXaxis()->SetBinLabel(i + 1, GetCentBinLabel(i).Data());
  fh1NMedConeCent = new TH1D("fh1NMedConeCent", "Number of med.-cl. cones in centrality bins;centrality;counts", fgkiNBinsCent, 0, fgkiNBinsCent);
  for(Int_t i = 0; i < fgkiNBinsCent; i++)
    fh1NMedConeCent->GetXaxis()->SetBinLabel(i + 1, GetCentBinLabel(i).Data());
  fh1AreaExcluded = new TH1D("fh1AreaExcluded", "Area of excluded cones in centrality bins;centrality;area", fgkiNBinsCent, 0, fgkiNBinsCent);
  for(Int_t i = 0; i < fgkiNBinsCent; i++)
    fh1AreaExcluded->GetXaxis()->SetBinLabel(i + 1, GetCentBinLabel(i).Data());
  fOutputListStd->Add(fh1EventCounterCut);
  fOutputListStd->Add(fh1EventCent);
  fOutputListStd->Add(fh1EventCent2);
  fOutputListStd->Add(fh1EventCent2Jets);
  fOutputListStd->Add(fh1EventCent2NoJets);
  fOutputListStd->Add(fh1NRndConeCent);
  fOutputListStd->Add(fh1NMedConeCent);
  fOutputListStd->Add(fh1AreaExcluded);
  fOutputListStd->Add(fh2EventCentTracks);

  fh1V0CandPerEvent = new TH1D("fh1V0CandPerEvent", "Number of all V0 candidates per event;candidates;events", 1000, 0, 1000);
  fOutputListStd->Add(fh1V0CandPerEvent);

  for(Int_t i = 0; i < fgkiNBinsCent; i++)
  {
    fh1EventCounterCutCent[i] = new TH1D(Form("fh1EventCounterCutCent_%d", i), Form("Number of events after filtering, cent %s;selection filter;counts", GetCentBinLabel(i).Data()), iNCategEvent, 0, iNCategEvent);
    for(Int_t j = 0; j < iNCategEvent; j++)
      fh1EventCounterCutCent[i]->GetXaxis()->SetBinLabel(j + 1, categEvent[j].Data());
    fh1V0CandPerEventCentK0s[i] = new TH1D(Form("fh1V0CandPerEventCentK0s_%d", i), Form("Number of selected K0s candidates per event, cent %s;candidates;events", GetCentBinLabel(i).Data()), 100, 0, 100);
    fh1V0CandPerEventCentLambda[i] = new TH1D(Form("fh1V0CandPerEventCentLambda_%d", i), Form("Number of selected Lambda candidates per event, cent %s;candidates;events", GetCentBinLabel(i).Data()), 100, 0, 100);
    fh1V0CandPerEventCentALambda[i] = new TH1D(Form("fh1V0CandPerEventCentALambda_%d", i), Form("Number of selected ALambda candidates per event, cent %s;candidates;events", GetCentBinLabel(i).Data()), 100, 0, 100);
    fh1V0CounterCentK0s[i] = new TH1D(Form("fh1V0CounterCentK0s_%d", i), Form("Number of K0s candidates after cuts, cent %s;cut;counts", GetCentBinLabel(i).Data()), fgkiNCategV0, 0, fgkiNCategV0);
    fh1V0CounterCentLambda[i] = new TH1D(Form("fh1V0CounterCentLambda_%d", i), Form("Number of Lambda candidates after cuts, cent %s;cut;counts", GetCentBinLabel(i).Data()), fgkiNCategV0, 0, fgkiNCategV0);
    fh1V0CounterCentALambda[i] = new TH1D(Form("fh1V0CounterCentALambda_%d", i), Form("Number of ALambda candidates after cuts, cent %s;cut;counts", GetCentBinLabel(i).Data()), fgkiNCategV0, 0, fgkiNCategV0);
    for(Int_t j = 0; j < fgkiNCategV0; j++)
    {
      fh1V0CounterCentK0s[i]->GetXaxis()->SetBinLabel(j + 1, categV0[j].Data());
      fh1V0CounterCentLambda[i]->GetXaxis()->SetBinLabel(j + 1, categV0[j].Data());
      fh1V0CounterCentALambda[i]->GetXaxis()->SetBinLabel(j + 1, categV0[j].Data());
    }
    fOutputListStd->Add(fh1EventCounterCutCent[i]);
    fOutputListStd->Add(fh1V0CandPerEventCentK0s[i]);
    fOutputListStd->Add(fh1V0CandPerEventCentLambda[i]);
    fOutputListStd->Add(fh1V0CandPerEventCentALambda[i]);
    fOutputListStd->Add(fh1V0CounterCentK0s[i]);
    fOutputListStd->Add(fh1V0CounterCentLambda[i]);
    fOutputListStd->Add(fh1V0CounterCentALambda[i]);
  }
  // pt binning for V0 and jets
  Int_t iNBinsPtV0 = fgkiNBinsPtV0Init;
  Double_t dPtV0Min = fgkdBinsPtV0[0];
  Double_t dPtV0Max = fgkdBinsPtV0[fgkiNBinsPtV0];
  Int_t iNJetPtBins = fgkiNBinsPtJetInit;
  Double_t dJetPtMin = fgkdBinsPtJet[0];
  Double_t dJetPtMax = fgkdBinsPtJet[fgkiNBinsPtJet];

  fh2CCK0s = new TH2D("fh2CCK0s", "K0s candidates in Lambda peak", fgkiNBinsMassK0s, fgkdMassK0sMin, fgkdMassK0sMax, iNBinsPtV0, dPtV0Min, dPtV0Max);
  fh2CCLambda = new TH2D("fh2CCLambda", "Lambda candidates in K0s peak", fgkiNBinsMassLambda, fgkdMassLambdaMin, fgkdMassLambdaMax, iNBinsPtV0, dPtV0Min, dPtV0Max);
  Int_t binsCorrel[3] = {fgkiNBinsMassK0s, fgkiNBinsMassLambda, iNBinsPtV0};
  Double_t xminCorrel[3] = {fgkdMassK0sMin, fgkdMassLambdaMin, dPtV0Min};
  Double_t xmaxCorrel[3] = {fgkdMassK0sMax, fgkdMassLambdaMax, dPtV0Max};
//  Int_t binsCorrel[3] = {200, 200, iNBinsPtV0};
//  Double_t xminCorrel[3] = {0, 0, dPtV0Min};
//  Double_t xmaxCorrel[3] = {2, 2, dPtV0Max};
  fh3CCMassCorrelBoth = new THnSparseD("fh3CCMassCorrelBoth", "Mass correlation: K0S && Lambda;m K0S;m Lambda;pT", 3, binsCorrel, xminCorrel, xmaxCorrel);
  fh3CCMassCorrelKNotL = new THnSparseD("fh3CCMassCorrelKNotL", "Mass correlation: K0S, not Lambda;m K0S;m Lambda;pT", 3, binsCorrel, xminCorrel, xmaxCorrel);
  fh3CCMassCorrelLNotK = new THnSparseD("fh3CCMassCorrelLNotK", "Mass correlation: Lambda, not K0S;m K0S;m Lambda;pT", 3, binsCorrel, xminCorrel, xmaxCorrel);
  fOutputListQA->Add(fh2CCK0s);
  fOutputListQA->Add(fh2CCLambda);
  fOutputListQA->Add(fh3CCMassCorrelBoth);
  fOutputListQA->Add(fh3CCMassCorrelKNotL);
  fOutputListQA->Add(fh3CCMassCorrelLNotK);

  Double_t dStepEtaV0 = 0.025;
  Double_t dRangeEtaV0Max = 0.8;
  const Int_t iNBinsEtaV0 = 2 * Int_t(dRangeEtaV0Max / dStepEtaV0);
  // inclusive
  const Int_t iNDimIncl = 3;
  Int_t binsKIncl[iNDimIncl] = {fgkiNBinsMassK0s, iNBinsPtV0, iNBinsEtaV0};
  Double_t xminKIncl[iNDimIncl] = {fgkdMassK0sMin, dPtV0Min, -dRangeEtaV0Max};
  Double_t xmaxKIncl[iNDimIncl] = {fgkdMassK0sMax, dPtV0Max, dRangeEtaV0Max};
  Int_t binsLIncl[iNDimIncl] = {fgkiNBinsMassLambda, iNBinsPtV0, iNBinsEtaV0};
  Double_t xminLIncl[iNDimIncl] = {fgkdMassLambdaMin, dPtV0Min, -dRangeEtaV0Max};
  Double_t xmaxLIncl[iNDimIncl] = {fgkdMassLambdaMax, dPtV0Max, dRangeEtaV0Max};
  // binning in jets
  const Int_t iNDimInJC = 4;
  Int_t binsKInJC[iNDimInJC] = {fgkiNBinsMassK0s, iNBinsPtV0, iNBinsEtaV0, iNJetPtBins};
  Double_t xminKInJC[iNDimInJC] = {fgkdMassK0sMin, dPtV0Min, -dRangeEtaV0Max, dJetPtMin};
  Double_t xmaxKInJC[iNDimInJC] = {fgkdMassK0sMax, dPtV0Max, dRangeEtaV0Max, dJetPtMax};
  Int_t binsLInJC[iNDimInJC] = {fgkiNBinsMassLambda, iNBinsPtV0, iNBinsEtaV0, iNJetPtBins};
  Double_t xminLInJC[iNDimInJC] = {fgkdMassLambdaMin, dPtV0Min, -dRangeEtaV0Max, dJetPtMin};
  Double_t xmaxLInJC[iNDimInJC] = {fgkdMassLambdaMax, dPtV0Max, dRangeEtaV0Max, dJetPtMax};

  // binning eff inclusive vs eta-pT
  Double_t dStepDeltaEta = 0.1;
  Double_t dRangeDeltaEtaMax = 0.5;
  const Int_t iNBinsDeltaEta = 2 * Int_t(dRangeDeltaEtaMax / dStepDeltaEta);
  Int_t binsEtaK[3] = {fgkiNBinsMassK0s, iNBinsPtV0, iNBinsEtaV0};
  Double_t xminEtaK[3] = {fgkdMassK0sMin, dPtV0Min, -dRangeEtaV0Max};
  Double_t xmaxEtaK[3] = {fgkdMassK0sMax, dPtV0Max, dRangeEtaV0Max};
  Int_t binsEtaL[3] = {fgkiNBinsMassLambda, iNBinsPtV0, iNBinsEtaV0};
  Double_t xminEtaL[3] = {fgkdMassLambdaMin, dPtV0Min, -dRangeEtaV0Max};
  Double_t xmaxEtaL[3] = {fgkdMassLambdaMax, dPtV0Max, dRangeEtaV0Max};
  // binning eff in jets vs eta-pT
  // associated
  Int_t binsEtaKInRec[5] = {fgkiNBinsMassK0s, iNBinsPtV0, iNBinsEtaV0, iNJetPtBins, iNBinsDeltaEta};
  Double_t xminEtaKInRec[5] = {fgkdMassK0sMin, dPtV0Min, -dRangeEtaV0Max, dJetPtMin, -dRangeDeltaEtaMax};
  Double_t xmaxEtaKInRec[5] = {fgkdMassK0sMax, dPtV0Max, dRangeEtaV0Max, dJetPtMax, dRangeDeltaEtaMax};
  Int_t binsEtaLInRec[5] = {fgkiNBinsMassLambda, iNBinsPtV0, iNBinsEtaV0, iNJetPtBins, iNBinsDeltaEta};
  Double_t xminEtaLInRec[5] = {fgkdMassLambdaMin, dPtV0Min, -dRangeEtaV0Max, dJetPtMin, -dRangeDeltaEtaMax};
  Double_t xmaxEtaLInRec[5] = {fgkdMassLambdaMax, dPtV0Max, dRangeEtaV0Max, dJetPtMax, dRangeDeltaEtaMax};
  // generated
  Int_t binsEtaInGen[4] = {iNBinsPtV0, iNBinsEtaV0, iNJetPtBins, iNBinsDeltaEta};
  Double_t xminEtaInGen[4] = {dPtV0Min, -dRangeEtaV0Max, dJetPtMin, -dRangeDeltaEtaMax};
  Double_t xmaxEtaInGen[4] = {dPtV0Max, dRangeEtaV0Max, dJetPtMax, dRangeDeltaEtaMax};
  // daughter eta: charge-etaD-ptD-etaV0-ptV0-ptJet
  const Int_t iNDimEtaD = 6;
  Int_t binsEtaDaughter[iNDimEtaD] = {2, 20, iNBinsPtV0, iNBinsEtaV0, iNBinsPtV0, iNJetPtBins};
  Double_t xminEtaDaughter[iNDimEtaD] = {0, -1, dPtV0Min, -dRangeEtaV0Max, dPtV0Min, dJetPtMin};
  Double_t xmaxEtaDaughter[iNDimEtaD] = {2, 1, dPtV0Max, dRangeEtaV0Max, dPtV0Max, dJetPtMax};

  for(Int_t i = 0; i < fgkiNBinsCent; i++)
  {
    fh1V0InvMassK0sCent[i] = new TH1D(Form("fh1V0InvMassK0sCent_%d", i), Form("K0s: V0 invariant mass, cent %s;#it{m}_{inv} (GeV/#it{c}^{2});counts", GetCentBinLabel(i).Data()), fgkiNBinsMassK0s, fgkdMassK0sMin, fgkdMassK0sMax);
    fh1V0InvMassLambdaCent[i] = new TH1D(Form("fh1V0InvMassLambdaCent_%d", i), Form("Lambda: V0 invariant mass, cent %s;#it{m}_{inv} (GeV/#it{c}^{2});counts", GetCentBinLabel(i).Data()), fgkiNBinsMassLambda, fgkdMassLambdaMin, fgkdMassLambdaMax);
    fh1V0InvMassALambdaCent[i] = new TH1D(Form("fh1V0InvMassALambdaCent_%d", i), Form("ALambda: V0 invariant mass, cent %s;#it{m}_{inv} (GeV/#it{c}^{2});counts", GetCentBinLabel(i).Data()), fgkiNBinsMassLambda, fgkdMassLambdaMin, fgkdMassLambdaMax);
    fOutputListStd->Add(fh1V0InvMassK0sCent[i]);
    fOutputListStd->Add(fh1V0InvMassLambdaCent[i]);
    fOutputListStd->Add(fh1V0InvMassALambdaCent[i]);
    // Inclusive
    fhnV0InclusiveK0s[i] = new THnSparseD(Form("fhnV0InclusiveK0s_C%d", i), "K0s: V0 invariant mass vs pt;#it{m}_{inv} (GeV/#it{c}^{2});pt (GeV/#it{c});counts", iNDimIncl, binsKIncl, xminKIncl, xmaxKIncl);
    fhnV0InclusiveLambda[i] = new THnSparseD(Form("fhnV0InclusiveLambda_C%d", i), "Lambda: V0 invariant mass vs pt;#it{m}_{inv} (GeV/#it{c}^{2});pt (GeV/#it{c});counts", iNDimIncl, binsLIncl, xminLIncl, xmaxLIncl);
    fhnV0InclusiveALambda[i] = new THnSparseD(Form("fhnV0InclusiveALambda_C%d", i), "ALambda: V0 invariant mass vs pt;#it{m}_{inv} (GeV/#it{c}^{2});pt (GeV/#it{c});counts", iNDimIncl, binsLIncl, xminLIncl, xmaxLIncl);
    fOutputListStd->Add(fhnV0InclusiveK0s[i]);
    fOutputListStd->Add(fhnV0InclusiveLambda[i]);
    fOutputListStd->Add(fhnV0InclusiveALambda[i]);
    // In cones
    fhnV0InJetK0s[i] = new THnSparseD(Form("fhnV0InJetK0s_%d", i), Form("K0s: Mass vs Pt in jets, cent: %s;#it{m}_{inv} (GeV/#it{c}^{2});#it{p}_{T}^{V0} (GeV/#it{c});#it{p}_{T}^{jet} (GeV/#it{c})", GetCentBinLabel(i).Data()), iNDimInJC, binsKInJC, xminKInJC, xmaxKInJC);
    fOutputListStd->Add(fhnV0InJetK0s[i]);
    fhnV0InPerpK0s[i] = new THnSparseD(Form("fhnV0InPerpK0s_%d", i), Form("K0s: Mass vs Pt in perp. cones, cent: %s;#it{m}_{inv} (GeV/#it{c}^{2});#it{p}_{T}^{V0} (GeV/#it{c});#it{p}_{T}^{jet} (GeV/#it{c})", GetCentBinLabel(i).Data()), iNDimInJC, binsKInJC, xminKInJC, xmaxKInJC);
    fOutputListStd->Add(fhnV0InPerpK0s[i]);
    fhnV0InRndK0s[i] = new THnSparseD(Form("fhnV0InRndK0s_%d", i), Form("K0s: Mass vs Pt in rnd. cones, cent: %s;#it{m}_{inv} (GeV/#it{c}^{2});#it{p}_{T}^{V0} (GeV/#it{c})", GetCentBinLabel(i).Data()), iNDimIncl, binsKIncl, xminKIncl, xmaxKIncl);
    fOutputListStd->Add(fhnV0InRndK0s[i]);
    fhnV0InMedK0s[i] = new THnSparseD(Form("fhnV0InMedK0s_%d", i), Form("K0s: Mass vs Pt in med.-cl. cones, cent: %s;#it{m}_{inv} (GeV/#it{c}^{2});#it{p}_{T}^{V0} (GeV/#it{c})", GetCentBinLabel(i).Data()), iNDimIncl, binsKIncl, xminKIncl, xmaxKIncl);
    fOutputListStd->Add(fhnV0InMedK0s[i]);
    fhnV0OutJetK0s[i] = new THnSparseD(Form("fhnV0OutJetK0s_%d", i), Form("K0s: Pt outside jets, cent: %s;#it{m}_{inv} (GeV/#it{c}^{2});#it{p}_{T}^{V0} (GeV/#it{c})", GetCentBinLabel(i).Data()), iNDimIncl, binsKIncl, xminKIncl, xmaxKIncl);
    fOutputListStd->Add(fhnV0OutJetK0s[i]);
    fhnV0NoJetK0s[i] = new THnSparseD(Form("fhnV0NoJetK0s_%d", i), Form("K0s: Pt in jet-less events, cent: %s;#it{m}_{inv} (GeV/#it{c}^{2});#it{p}_{T}^{V0} (GeV/#it{c})", GetCentBinLabel(i).Data()), iNDimIncl, binsKIncl, xminKIncl, xmaxKIncl);
    fOutputListStd->Add(fhnV0NoJetK0s[i]);
    fhnV0InJetLambda[i] = new THnSparseD(Form("fhnV0InJetLambda_%d", i), Form("Lambda: Mass vs Pt in jets, cent: %s;#it{m}_{inv} (GeV/#it{c}^{2});#it{p}_{T}^{V0} (GeV/#it{c});#it{p}_{T}^{jet} (GeV/#it{c})", GetCentBinLabel(i).Data()), iNDimInJC, binsLInJC, xminLInJC, xmaxLInJC);
    fOutputListStd->Add(fhnV0InJetLambda[i]);
    fhnV0InPerpLambda[i] = new THnSparseD(Form("fhnV0InPerpLambda_%d", i), Form("Lambda: Mass vs Pt in perp. cones, cent: %s;#it{m}_{inv} (GeV/#it{c}^{2});#it{p}_{T}^{V0} (GeV/#it{c});#it{p}_{T}^{jet} (GeV/#it{c})", GetCentBinLabel(i).Data()), iNDimInJC, binsLInJC, xminLInJC, xmaxLInJC);
    fOutputListStd->Add(fhnV0InPerpLambda[i]);
    fhnV0InRndLambda[i] = new THnSparseD(Form("fhnV0InRndLambda_%d", i), Form("Lambda: Mass vs Pt in rnd. cones, cent: %s;#it{m}_{inv} (GeV/#it{c}^{2});#it{p}_{T}^{V0} (GeV/#it{c})", GetCentBinLabel(i).Data()), iNDimIncl, binsLIncl, xminLIncl, xmaxLIncl);
    fOutputListStd->Add(fhnV0InRndLambda[i]);
    fhnV0InMedLambda[i] = new THnSparseD(Form("fhnV0InMedLambda_%d", i), Form("Lambda: Mass vs Pt in med.-cl. cones, cent: %s;#it{m}_{inv} (GeV/#it{c}^{2});#it{p}_{T}^{V0} (GeV/#it{c})", GetCentBinLabel(i).Data()), iNDimIncl, binsLIncl, xminLIncl, xmaxLIncl);
    fOutputListStd->Add(fhnV0InMedLambda[i]);
    fhnV0OutJetLambda[i] = new THnSparseD(Form("fhnV0OutJetLambda_%d", i), Form("Lambda: Pt outside jets, cent: %s;#it{m}_{inv} (GeV/#it{c}^{2});#it{p}_{T}^{V0} (GeV/#it{c})", GetCentBinLabel(i).Data()), iNDimIncl, binsLIncl, xminLIncl, xmaxLIncl);
    fOutputListStd->Add(fhnV0OutJetLambda[i]);
    fhnV0NoJetLambda[i] = new THnSparseD(Form("fhnV0NoJetLambda_%d", i), Form("Lambda: Pt in jet-less events, cent: %s;#it{m}_{inv} (GeV/#it{c}^{2});#it{p}_{T}^{V0} (GeV/#it{c})", GetCentBinLabel(i).Data()), iNDimIncl, binsLIncl, xminLIncl, xmaxLIncl);
    fOutputListStd->Add(fhnV0NoJetLambda[i]);
    fhnV0InJetALambda[i] = new THnSparseD(Form("fhnV0InJetALambda_%d", i), Form("ALambda: Mass vs Pt in jets, cent: %s;#it{m}_{inv} (GeV/#it{c}^{2});#it{p}_{T}^{V0} (GeV/#it{c});#it{p}_{T}^{jet} (GeV/#it{c})", GetCentBinLabel(i).Data()), iNDimInJC, binsLInJC, xminLInJC, xmaxLInJC);
    fOutputListStd->Add(fhnV0InJetALambda[i]);
    fhnV0InPerpALambda[i] = new THnSparseD(Form("fhnV0InPerpALambda_%d", i), Form("ALambda: Mass vs Pt in perp. cones, cent: %s;#it{m}_{inv} (GeV/#it{c}^{2});#it{p}_{T}^{V0} (GeV/#it{c});#it{p}_{T}^{jet} (GeV/#it{c})", GetCentBinLabel(i).Data()), iNDimInJC, binsLInJC, xminLInJC, xmaxLInJC);
    fOutputListStd->Add(fhnV0InPerpALambda[i]);
    fhnV0InRndALambda[i] = new THnSparseD(Form("fhnV0InRndALambda_%d", i), Form("ALambda: Mass vs Pt in rnd. cones, cent: %s;#it{m}_{inv} (GeV/#it{c}^{2});#it{p}_{T}^{V0} (GeV/#it{c})", GetCentBinLabel(i).Data()), iNDimIncl, binsLIncl, xminLIncl, xmaxLIncl);
    fOutputListStd->Add(fhnV0InRndALambda[i]);
    fhnV0InMedALambda[i] = new THnSparseD(Form("fhnV0InMedALambda_%d", i), Form("ALambda: Mass vs Pt in med.-cl. cones, cent: %s;#it{m}_{inv} (GeV/#it{c}^{2});#it{p}_{T}^{V0} (GeV/#it{c})", GetCentBinLabel(i).Data()), iNDimIncl, binsLIncl, xminLIncl, xmaxLIncl);
    fOutputListStd->Add(fhnV0InMedALambda[i]);
    fhnV0OutJetALambda[i] = new THnSparseD(Form("fhnV0OutJetALambda_%d", i), Form("ALambda: Pt outside jets, cent: %s;#it{m}_{inv} (GeV/#it{c}^{2});#it{p}_{T}^{V0} (GeV/#it{c})", GetCentBinLabel(i).Data()), iNDimIncl, binsLIncl, xminLIncl, xmaxLIncl);
    fOutputListStd->Add(fhnV0OutJetALambda[i]);
    fhnV0NoJetALambda[i] = new THnSparseD(Form("fhnV0NoJetALambda_%d", i), Form("ALambda: Pt in jet-less events, cent: %s;#it{m}_{inv} (GeV/#it{c}^{2});#it{p}_{T}^{V0} (GeV/#it{c})", GetCentBinLabel(i).Data()), iNDimIncl, binsLIncl, xminLIncl, xmaxLIncl);
    fOutputListStd->Add(fhnV0NoJetALambda[i]);

    fh2V0PtJetAngleK0s[i] = new TH2D(Form("fh2V0PtJetAngleK0s_%d", i), Form("K0s: #it{p}_{T}^{jet} vs angle V0-jet, cent: %s;#it{p}_{T}^{jet};#it{#alpha}", GetCentBinLabel(i).Data()), iNJetPtBins, dJetPtMin, dJetPtMax, 100, 0, fdRadiusJet + 0.1);
    fOutputListStd->Add(fh2V0PtJetAngleK0s[i]);
    fh2V0PtJetAngleLambda[i] = new TH2D(Form("fh2V0PtJetAngleLambda_%d", i), Form("Lambda: #it{p}_{T}^{jet} vs angle V0-jet, cent: %s;#it{p}_{T}^{jet};#it{#alpha}", GetCentBinLabel(i).Data()), iNJetPtBins, dJetPtMin, dJetPtMax, 100, 0, fdRadiusJet + 0.1);
    fOutputListStd->Add(fh2V0PtJetAngleLambda[i]);
    fh2V0PtJetAngleALambda[i] = new TH2D(Form("fh2V0PtJetAngleALambda_%d", i), Form("ALambda: #it{p}_{T}^{jet} vs angle V0-jet, cent: %s;#it{p}_{T}^{jet};#it{#alpha}", GetCentBinLabel(i).Data()), iNJetPtBins, dJetPtMin, dJetPtMax, 100, 0, fdRadiusJet + 0.1);
    fOutputListStd->Add(fh2V0PtJetAngleALambda[i]);

    fh1DCAInK0s[i] = new TH1D(Form("fh1DCAInK0s_%d", i), Form("K0s in jets: DCA daughters, cent %s;DCA (#sigma)", GetCentBinLabel(i).Data()), 50, 0, 1);
    fOutputListQA->Add(fh1DCAInK0s[i]);
    fh1DCAInLambda[i] = new TH1D(Form("fh1DCAInLambda_%d", i), Form("Lambda in jets: DCA daughters, cent %s;DCA (#sigma)", GetCentBinLabel(i).Data()), 50, 0, 1);
    fOutputListQA->Add(fh1DCAInLambda[i]);
    fh1DCAInALambda[i] = new TH1D(Form("fh1DCAInALambda_%d", i), Form("ALambda in jets: DCA daughters, cent %s;DCA (#sigma)", GetCentBinLabel(i).Data()), 50, 0, 1);
    fOutputListQA->Add(fh1DCAInALambda[i]);

    fh1DCAOutK0s[i] = new TH1D(Form("fh1DCAOutK0s_%d", i), Form("K0s outside jets: DCA daughters, cent %s;DCA (#sigma)", GetCentBinLabel(i).Data()), 50, 0, 1);
    fOutputListQA->Add(fh1DCAOutK0s[i]);
    fh1DCAOutLambda[i] = new TH1D(Form("fh1DCAOutLambda_%d", i), Form("Lambda outside jets: DCA daughters, cent %s;DCA (#sigma)", GetCentBinLabel(i).Data()), 50, 0, 1);
    fOutputListQA->Add(fh1DCAOutLambda[i]);
    fh1DCAOutALambda[i] = new TH1D(Form("fh1DCAOutALambda_%d", i), Form("ALambda outside jets: DCA daughters, cent %s;DCA (#sigma)", GetCentBinLabel(i).Data()), 50, 0, 1);
    fOutputListQA->Add(fh1DCAOutALambda[i]);

    /*
    fh1DeltaZK0s[i] = new TH1D(Form("fh1DeltaZK0s_%d", i), Form("K0s: #Delta#it{z} vertices, cent %s;#it{z} (cm)", GetCentBinLabel(i).Data()), 50, -10, 10);
    fOutputListQA->Add(fh1DeltaZK0s[i]);
    fh1DeltaZLambda[i] = new TH1D(Form("fh1DeltaZLambda_%d", i), Form("Lambda: #Delta#it{z} vertices, cent %s;#it{z} (cm)", GetCentBinLabel(i).Data()), 50, -10, 10);
    fOutputListQA->Add(fh1DeltaZLambda[i]);
    fh1DeltaZALambda[i] = new TH1D(Form("fh1DeltaZALambda_%d", i), Form("ALambda: #Delta#it{z} vertices, cent %s;#it{z} (cm)", GetCentBinLabel(i).Data()), 50, -10, 10);
    fOutputListQA->Add(fh1DeltaZALambda[i]);
    */

    // jet histograms
    fh1PtJet[i] = new TH1D(Form("fh1PtJet_%d", i), Form("Jet pt spectrum, cent: %s;#it{p}_{T} jet (GeV/#it{c})", GetCentBinLabel(i).Data()), iNJetPtBins, dJetPtMin, dJetPtMax);
    fOutputListStd->Add(fh1PtJet[i]);
    fh1EtaJet[i] = new TH1D(Form("fh1EtaJet_%d", i), Form("Jet eta spectrum, cent: %s;#it{#eta} jet", GetCentBinLabel(i).Data()), 80, -1., 1.);
    fOutputListStd->Add(fh1EtaJet[i]);
    fh2EtaPtJet[i] = new TH2D(Form("fh2EtaPtJet_%d", i), Form("Jet eta vs pT spectrum, cent: %s;#it{#eta} jet;#it{p}_{T} jet (GeV/#it{c})", GetCentBinLabel(i).Data()), 80, -1., 1., iNJetPtBins, dJetPtMin, dJetPtMax);
    fOutputListStd->Add(fh2EtaPtJet[i]);
    fh2EtaPhiRndCone[i] = new TH2D(Form("fh2EtaPhiRndCone_%d", i), Form("Rnd. cones: eta vs phi, cent: %s;#it{#eta} cone;#it{#phi} cone", GetCentBinLabel(i).Data()), 80, -1., 1., 100, 0., TMath::TwoPi());
    fOutputListStd->Add(fh2EtaPhiRndCone[i]);
    fh2EtaPhiMedCone[i] = new TH2D(Form("fh2EtaPhiMedCone_%d", i), Form("Med.-cl. cones: eta vs phi, cent: %s;#it{#eta} cone;#it{#phi} cone", GetCentBinLabel(i).Data()), 80, -1., 1., 100, 0., TMath::TwoPi());
    fOutputListStd->Add(fh2EtaPhiMedCone[i]);
    fh1PhiJet[i] = new TH1D(Form("fh1PhiJet_%d", i), Form("Jet phi spectrum, cent: %s;#it{#phi} jet", GetCentBinLabel(i).Data()), 100, 0., TMath::TwoPi());
    fOutputListStd->Add(fh1PhiJet[i]);
    fh1NJetPerEvent[i] = new TH1D(Form("fh1NJetPerEvent_%d", i), Form("Number of selected jets per event, cent: %s;# jets;# events", GetCentBinLabel(i).Data()), 100, 0., 100.);
    fOutputListStd->Add(fh1NJetPerEvent[i]);
    // event histograms
    fh1VtxZ[i] = new TH1D(Form("fh1VtxZ_%d", i), Form("#it{z} coordinate of the primary vertex, cent: %s;#it{z} (cm)", GetCentBinLabel(i).Data()), 150, -1.5 * fdCutVertexZ, 1.5 * fdCutVertexZ);
    fOutputListQA->Add(fh1VtxZ[i]);
    fh2VtxXY[i] = new TH2D(Form("fh2VtxXY_%d", i), Form("#it{xy} coordinate of the primary vertex, cent: %s;#it{x} (cm);#it{y} (cm)", GetCentBinLabel(i).Data()), 200, -0.2, 0.2, 500, -0.5, 0.5);
    fOutputListQA->Add(fh2VtxXY[i]);
    // fOutputListStd->Add([i]);
    if(fbMCAnalysis)
    {
      // inclusive pt
      fh1V0K0sPtMCGen[i] = new TH1D(Form("fh1V0K0sPtMCGen_%d", i), Form("MC K0s generated: pt spectrum, cent: %s;MC #it{p}_{T} (GeV/#it{c})", GetCentBinLabel(i).Data()), iNBinsPtV0, dPtV0Min, dPtV0Max);
      fOutputListMC->Add(fh1V0K0sPtMCGen[i]);
      fh2V0K0sPtMassMCRec[i] = new TH2D(Form("fh2V0K0sPtMassMCRec_%d", i), Form("MC K0s associated: pt-m spectrum, cent: %s;MC #it{p}_{T} (GeV/#it{c});#it{m}_{inv} (GeV/#it{c}^{2})", GetCentBinLabel(i).Data()), iNBinsPtV0, dPtV0Min, dPtV0Max, fgkiNBinsMassK0s, fgkdMassK0sMin, fgkdMassK0sMax);
      fOutputListMC->Add(fh2V0K0sPtMassMCRec[i]);
      fh1V0K0sPtMCRecFalse[i] = new TH1D(Form("fh1V0K0sPtMCRecFalse_%d", i), Form("MC K0s false: pt spectrum, cent: %s;MC #it{p}_{T} (GeV/#it{c})", GetCentBinLabel(i).Data()), iNBinsPtV0, dPtV0Min, dPtV0Max);
      fOutputListMC->Add(fh1V0K0sPtMCRecFalse[i]);
      // inclusive pt-eta
      fh2V0K0sEtaPtMCGen[i] = new TH2D(Form("fh2V0K0sEtaPtMCGen_%d", i), Form("MC K0s generated: pt-eta spectrum, cent: %s;MC #it{p}_{T} (GeV/#it{c});#eta", GetCentBinLabel(i).Data()), iNBinsPtV0, dPtV0Min, dPtV0Max, iNBinsEtaV0, -dRangeEtaV0Max, dRangeEtaV0Max);
      fOutputListMC->Add(fh2V0K0sEtaPtMCGen[i]);
      fh3V0K0sEtaPtMassMCRec[i] = new THnSparseD(Form("fh3V0K0sEtaPtMassMCRec_%d", i), Form("MC K0s associated: m-pt-eta spectrum, cent: %s;#it{m}_{inv} (GeV/#it{c}^{2});MC #it{p}_{T} (GeV/#it{c});#eta", GetCentBinLabel(i).Data()), 3, binsEtaK, xminEtaK, xmaxEtaK);
      fOutputListMC->Add(fh3V0K0sEtaPtMassMCRec[i]);
      // in jet pt
      fh2V0K0sInJetPtMCGen[i] = new TH2D(Form("fh2V0K0sInJetPtMCGen_%d", i), Form("MC K0s in jet generated: pt-ptJet spectrum, cent: %s;MC #it{p}_{T} (GeV/#it{c});#it{p}_{T}^{jet} (GeV/#it{c})", GetCentBinLabel(i).Data()), iNBinsPtV0, dPtV0Min, dPtV0Max, iNJetPtBins, dJetPtMin, dJetPtMax);
      fOutputListMC->Add(fh2V0K0sInJetPtMCGen[i]);
      fh3V0K0sInJetPtMassMCRec[i] = new THnSparseD(Form("fh3V0K0sInJetPtMassMCRec_%d", i), Form("MC K0s in jet associated: m-pt-ptJet spectrum, cent: %s;#it{m}_{inv} (GeV/#it{c}^{2});MC #it{p}_{T} (GeV/#it{c});#it{p}_{T}^{jet} (GeV/#it{c})", GetCentBinLabel(i).Data()), iNDimInJC, binsKInJC, xminKInJC, xmaxKInJC);
      fOutputListMC->Add(fh3V0K0sInJetPtMassMCRec[i]);
      // in jet pt-eta
      fh3V0K0sInJetEtaPtMCGen[i] = new THnSparseD(Form("fh3V0K0sInJetEtaPtMCGen_%d", i), Form("MC K0s generated: pt-eta-ptJet spectrum, cent: %s;MC #it{p}_{T} (GeV/#it{c});#eta;#it{p}_{T}^{jet} (GeV/#it{c})", GetCentBinLabel(i).Data()), 4, binsEtaInGen, xminEtaInGen, xmaxEtaInGen);
      fOutputListMC->Add(fh3V0K0sInJetEtaPtMCGen[i]);
      fh4V0K0sInJetEtaPtMassMCRec[i] = new THnSparseD(Form("fh4V0K0sInJetEtaPtMassMCRec_%d", i), Form("MC K0s associated: m-pt-eta-ptJet spectrum, cent: %s;#it{m}_{inv} (GeV/#it{c}^{2});MC #it{p}_{T} (GeV/#it{c});#eta;#it{p}_{T}^{jet} (GeV/#it{c})", GetCentBinLabel(i).Data()), 5, binsEtaKInRec, xminEtaKInRec, xmaxEtaKInRec);
      fOutputListMC->Add(fh4V0K0sInJetEtaPtMassMCRec[i]);

      fh2V0K0sMCResolMPt[i] = new TH2D(Form("fh2V0K0sMCResolMPt_%d", i), Form("MC K0s associated: #Delta#it{m} vs pt, cent %s;#Delta#it{m} = #it{m}_{inv} - #it{m}_{true} (GeV/#it{c}^{2});#it{p}_{T}^{rec} (GeV/#it{c})", GetCentBinLabel(i).Data()), 100, -0.02, 0.02, iNBinsPtV0, dPtV0Min, dPtV0Max);
      fOutputListMC->Add(fh2V0K0sMCResolMPt[i]);
      fh2V0K0sMCPtGenPtRec[i] = new TH2D(Form("fh2V0K0sMCPtGenPtRec_%d", i), Form("MC K0s associated: pt gen vs pt rec, cent %s;#it{p}_{T}^{gen} (GeV/#it{c});#it{p}_{T}^{rec} (GeV/#it{c})", GetCentBinLabel(i).Data()), iNBinsPtV0, dPtV0Min, dPtV0Max, iNBinsPtV0, dPtV0Min, dPtV0Max);
      fOutputListMC->Add(fh2V0K0sMCPtGenPtRec[i]);

      // inclusive pt
      fh1V0LambdaPtMCGen[i] = new TH1D(Form("fh1V0LambdaPtMCGen_%d", i), Form("MC Lambda generated: pt spectrum, cent: %s;MC #it{p}_{T} (GeV/#it{c})", GetCentBinLabel(i).Data()), iNBinsPtV0, dPtV0Min, dPtV0Max);
      fOutputListMC->Add(fh1V0LambdaPtMCGen[i]);
      fh2V0LambdaPtMassMCRec[i] = new TH2D(Form("fh2V0LambdaPtMassMCRec_%d", i), Form("MC Lambda associated: pt-m spectrum, cent: %s;MC #it{p}_{T} (GeV/#it{c});#it{m}_{inv} (GeV/#it{c}^{2})", GetCentBinLabel(i).Data()), iNBinsPtV0, dPtV0Min, dPtV0Max, fgkiNBinsMassLambda, fgkdMassLambdaMin, fgkdMassLambdaMax);
      fOutputListMC->Add(fh2V0LambdaPtMassMCRec[i]);
      fh1V0LambdaPtMCRecFalse[i] = new TH1D(Form("fh1V0LambdaPtMCRecFalse_%d", i), Form("MC Lambda false: pt spectrum, cent: %s;MC #it{p}_{T} (GeV/#it{c})", GetCentBinLabel(i).Data()), iNBinsPtV0, dPtV0Min, dPtV0Max);
      fOutputListMC->Add(fh1V0LambdaPtMCRecFalse[i]);
      // inclusive pt-eta
      fh2V0LambdaEtaPtMCGen[i] = new TH2D(Form("fh2V0LambdaEtaPtMCGen_%d", i), Form("MC Lambda generated: pt-eta spectrum, cent: %s;MC #it{p}_{T} (GeV/#it{c});#eta", GetCentBinLabel(i).Data()), iNBinsPtV0, dPtV0Min, dPtV0Max, iNBinsEtaV0, -dRangeEtaV0Max, dRangeEtaV0Max);
      fOutputListMC->Add(fh2V0LambdaEtaPtMCGen[i]);
      fh3V0LambdaEtaPtMassMCRec[i] = new THnSparseD(Form("fh3V0LambdaEtaPtMassMCRec_%d", i), Form("MC Lambda associated: m-pt-eta spectrum, cent: %s;#it{m}_{inv} (GeV/#it{c}^{2});MC #it{p}_{T} (GeV/#it{c});#eta", GetCentBinLabel(i).Data()), 3, binsEtaL, xminEtaL, xmaxEtaL);
      fOutputListMC->Add(fh3V0LambdaEtaPtMassMCRec[i]);
      // in jet pt
      fh2V0LambdaInJetPtMCGen[i] = new TH2D(Form("fh2V0LambdaInJetPtMCGen_%d", i), Form("MC Lambda in jet generated: pt-ptJet spectrum, cent: %s;MC #it{p}_{T} (GeV/#it{c});#it{p}_{T}^{jet} (GeV/#it{c})", GetCentBinLabel(i).Data()), iNBinsPtV0, dPtV0Min, dPtV0Max, iNJetPtBins, dJetPtMin, dJetPtMax);
      fOutputListMC->Add(fh2V0LambdaInJetPtMCGen[i]);
      fh3V0LambdaInJetPtMassMCRec[i] = new THnSparseD(Form("fh3V0LambdaInJetPtMassMCRec_%d", i), Form("MC Lambda in jet associated: m-pt-ptJet spectrum, cent: %s;#it{m}_{inv} (GeV/#it{c}^{2});MC #it{p}_{T} (GeV/#it{c});#it{p}_{T}^{jet} (GeV/#it{c})", GetCentBinLabel(i).Data()), iNDimInJC, binsLInJC, xminLInJC, xmaxLInJC);
      fOutputListMC->Add(fh3V0LambdaInJetPtMassMCRec[i]);
      // in jet pt-eta
      fh3V0LambdaInJetEtaPtMCGen[i] = new THnSparseD(Form("fh3V0LambdaInJetEtaPtMCGen_%d", i), Form("MC Lambda generated: pt-eta-ptJet spectrum, cent: %s;MC #it{p}_{T} (GeV/#it{c});#eta;#it{p}_{T}^{jet} (GeV/#it{c})", GetCentBinLabel(i).Data()), 4, binsEtaInGen, xminEtaInGen, xmaxEtaInGen);
      fOutputListMC->Add(fh3V0LambdaInJetEtaPtMCGen[i]);
      fh4V0LambdaInJetEtaPtMassMCRec[i] = new THnSparseD(Form("fh4V0LambdaInJetEtaPtMassMCRec_%d", i), Form("MC Lambda associated: m-pt-eta-ptJet spectrum, cent: %s;#it{m}_{inv} (GeV/#it{c}^{2});MC #it{p}_{T} (GeV/#it{c});#eta;#it{p}_{T}^{jet} (GeV/#it{c})", GetCentBinLabel(i).Data()), 5, binsEtaLInRec, xminEtaLInRec, xmaxEtaLInRec);
      fOutputListMC->Add(fh4V0LambdaInJetEtaPtMassMCRec[i]);

      fh2V0LambdaMCResolMPt[i] = new TH2D(Form("fh2V0LambdaMCResolMPt_%d", i), Form("MC Lambda associated: #Delta#it{m} vs pt, cent %s;#Delta#it{m} = #it{m}_{inv} - #it{m}_{true} (GeV/#it{c}^{2});#it{p}_{T}^{rec} (GeV/#it{c})", GetCentBinLabel(i).Data()), 100, -0.02, 0.02, iNBinsPtV0, dPtV0Min, dPtV0Max);
      fOutputListMC->Add(fh2V0LambdaMCResolMPt[i]);
      fh2V0LambdaMCPtGenPtRec[i] = new TH2D(Form("fh2V0LambdaMCPtGenPtRec_%d", i), Form("MC Lambda associated: pt gen vs pt rec, cent %s;#it{p}_{T}^{gen} (GeV/#it{c});#it{p}_{T}^{rec} (GeV/#it{c})", GetCentBinLabel(i).Data()), iNBinsPtV0, dPtV0Min, dPtV0Max, iNBinsPtV0, dPtV0Min, dPtV0Max);
      fOutputListMC->Add(fh2V0LambdaMCPtGenPtRec[i]);

      // inclusive pt
      fh1V0ALambdaPtMCGen[i] = new TH1D(Form("fh1V0ALambdaPtMCGen_%d", i), Form("MC ALambda generated: pt spectrum, cent: %s;MC #it{p}_{T} (GeV/#it{c})", GetCentBinLabel(i).Data()), iNBinsPtV0, dPtV0Min, dPtV0Max);
      fOutputListMC->Add(fh1V0ALambdaPtMCGen[i]);
      fh2V0ALambdaPtMassMCRec[i] = new TH2D(Form("fh2V0ALambdaPtMassMCRec_%d", i), Form("MC ALambda associated: pt-m spectrum, cent: %s;MC #it{p}_{T} (GeV/#it{c});#it{m}_{inv} (GeV/#it{c}^{2})", GetCentBinLabel(i).Data()), iNBinsPtV0, dPtV0Min, dPtV0Max, fgkiNBinsMassLambda, fgkdMassLambdaMin, fgkdMassLambdaMax);
      fOutputListMC->Add(fh2V0ALambdaPtMassMCRec[i]);
      fh1V0ALambdaPtMCRecFalse[i] = new TH1D(Form("fh1V0ALambdaPtMCRecFalse_%d", i), Form("MC ALambda false: pt spectrum, cent: %s;MC #it{p}_{T} (GeV/#it{c})", GetCentBinLabel(i).Data()), iNBinsPtV0, dPtV0Min, dPtV0Max);
      fOutputListMC->Add(fh1V0ALambdaPtMCRecFalse[i]);
      // inclusive pt-eta
      fh2V0ALambdaEtaPtMCGen[i] = new TH2D(Form("fh2V0ALambdaEtaPtMCGen_%d", i), Form("MC ALambda generated: pt-eta spectrum, cent: %s;MC #it{p}_{T} (GeV/#it{c});#eta", GetCentBinLabel(i).Data()), iNBinsPtV0, dPtV0Min, dPtV0Max, iNBinsEtaV0, -dRangeEtaV0Max, dRangeEtaV0Max);
      fOutputListMC->Add(fh2V0ALambdaEtaPtMCGen[i]);
      fh3V0ALambdaEtaPtMassMCRec[i] = new THnSparseD(Form("fh3V0ALambdaEtaPtMassMCRec_%d", i), Form("MC ALambda associated: m-pt-eta spectrum, cent: %s;#it{m}_{inv} (GeV/#it{c}^{2});MC #it{p}_{T} (GeV/#it{c});#eta", GetCentBinLabel(i).Data()), 3, binsEtaL, xminEtaL, xmaxEtaL);
      fOutputListMC->Add(fh3V0ALambdaEtaPtMassMCRec[i]);
      // in jet pt
      fh2V0ALambdaInJetPtMCGen[i] = new TH2D(Form("fh2V0ALambdaInJetPtMCGen_%d", i), Form("MC ALambda in jet generated: pt-ptJet spectrum, cent: %s;MC #it{p}_{T} (GeV/#it{c});#it{p}_{T}^{jet} (GeV/#it{c})", GetCentBinLabel(i).Data()), iNBinsPtV0, dPtV0Min, dPtV0Max, iNJetPtBins, dJetPtMin, dJetPtMax);
      fOutputListMC->Add(fh2V0ALambdaInJetPtMCGen[i]);
      fh3V0ALambdaInJetPtMassMCRec[i] = new THnSparseD(Form("fh3V0ALambdaInJetPtMassMCRec_%d", i), Form("MC ALambda in jet associated: m-pt-ptJet spectrum, cent: %s;#it{m}_{inv} (GeV/#it{c}^{2});MC #it{p}_{T} (GeV/#it{c});#it{p}_{T}^{jet} (GeV/#it{c})", GetCentBinLabel(i).Data()), iNDimInJC, binsLInJC, xminLInJC, xmaxLInJC);
      fOutputListMC->Add(fh3V0ALambdaInJetPtMassMCRec[i]);
      // in jet pt-eta
      fh3V0ALambdaInJetEtaPtMCGen[i] = new THnSparseD(Form("fh3V0ALambdaInJetEtaPtMCGen_%d", i), Form("MC ALambda generated: pt-eta-ptJet spectrum, cent: %s;MC #it{p}_{T} (GeV/#it{c});#eta;#it{p}_{T}^{jet} (GeV/#it{c})", GetCentBinLabel(i).Data()), 4, binsEtaInGen, xminEtaInGen, xmaxEtaInGen);
      fOutputListMC->Add(fh3V0ALambdaInJetEtaPtMCGen[i]);
      fh4V0ALambdaInJetEtaPtMassMCRec[i] = new THnSparseD(Form("fh4V0ALambdaInJetEtaPtMassMCRec_%d", i), Form("MC ALambda associated: m-pt-eta-ptJet spectrum, cent: %s;#it{m}_{inv} (GeV/#it{c}^{2});MC #it{p}_{T} (GeV/#it{c});#eta;#it{p}_{T}^{jet} (GeV/#it{c})", GetCentBinLabel(i).Data()), 5, binsEtaLInRec, xminEtaLInRec, xmaxEtaLInRec);
      fOutputListMC->Add(fh4V0ALambdaInJetEtaPtMassMCRec[i]);

      fh2V0ALambdaMCResolMPt[i] = new TH2D(Form("fh2V0ALambdaMCResolMPt_%d", i), Form("MC ALambda associated: #Delta#it{m} vs pt, cent %s;#Delta#it{m} = #it{m}_{inv} - #it{m}_{true} (GeV/#it{c}^{2});#it{p}_{T}^{rec} (GeV/#it{c})", GetCentBinLabel(i).Data()), 100, -0.02, 0.02, iNBinsPtV0, dPtV0Min, dPtV0Max);
      fOutputListMC->Add(fh2V0ALambdaMCResolMPt[i]);
      fh2V0ALambdaMCPtGenPtRec[i] = new TH2D(Form("fh2V0ALambdaMCPtGenPtRec_%d", i), Form("MC ALambda associated: pt gen vs pt rec, cent %s;#it{p}_{T}^{gen} (GeV/#it{c});#it{p}_{T}^{rec} (GeV/#it{c})", GetCentBinLabel(i).Data()), iNBinsPtV0, dPtV0Min, dPtV0Max, iNBinsPtV0, dPtV0Min, dPtV0Max);
      fOutputListMC->Add(fh2V0ALambdaMCPtGenPtRec[i]);

      Int_t iNBinsPtXi = 80;
      Double_t dPtXiMin = 0;
      Double_t dPtXiMax = 8;
      const Int_t iNDimFD = 3;
      Int_t binsFD[iNDimFD] = {iNBinsPtV0, iNBinsPtXi, iNJetPtBins};
      Double_t xminFD[iNDimFD] = {dPtV0Min, dPtXiMin, dJetPtMin};
      Double_t xmaxFD[iNDimFD] = {dPtV0Max, dPtXiMax, dJetPtMax};
      fhnV0LambdaInclMCFD[i] = new THnSparseD(Form("fhnV0LambdaInclMCFD_%d", i), Form("MC Lambda associated, inclusive, from Xi: pt-pt, cent %s;#it{p}_{T}^{#Lambda,gen.} (GeV/#it{c});#it{p}_{T}^{#Xi,gen.} (GeV/#it{c});#it{p}_{T}^{jet} (GeV/#it{c})", GetCentBinLabel(i).Data()), iNDimFD, binsFD, xminFD, xmaxFD);
      fOutputListMC->Add(fhnV0LambdaInclMCFD[i]);
      fhnV0LambdaInJetsMCFD[i] = new THnSparseD(Form("fhnV0LambdaInJetsMCFD_%d", i), Form("MC Lambda associated, in JC, from Xi: pt-pt-ptJet, cent %s;#it{p}_{T}^{#Lambda,gen.} (GeV/#it{c});#it{p}_{T}^{#Xi,gen.} (GeV/#it{c});#it{p}_{T}^{jet} (GeV/#it{c})", GetCentBinLabel(i).Data()), iNDimFD, binsFD, xminFD, xmaxFD);
      fOutputListMC->Add(fhnV0LambdaInJetsMCFD[i]);
      fhnV0LambdaBulkMCFD[i] = new THnSparseD(Form("fhnV0LambdaBulkMCFD_%d", i), Form("MC Lambda associated, in no jet events, from Xi: pt-pt, cent %s;#it{p}_{T}^{#Lambda,gen.} (GeV/#it{c});#it{p}_{T}^{#Xi,gen.} (GeV/#it{c});#it{p}_{T}^{jet} (GeV/#it{c})", GetCentBinLabel(i).Data()), iNDimFD, binsFD, xminFD, xmaxFD);
      fOutputListMC->Add(fhnV0LambdaBulkMCFD[i]);
      fh1V0XiPtMCGen[i] = new TH1D(Form("fh1V0XiPtMCGen_%d", i), Form("MC Xi^{-} generated: Pt spectrum, cent %s;#it{p}_{T}^{#Xi^{-},gen.} (GeV/#it{c})", GetCentBinLabel(i).Data()), iNBinsPtXi, dPtXiMin, dPtXiMax);
      fOutputListMC->Add(fh1V0XiPtMCGen[i]);
      fhnV0ALambdaInclMCFD[i] = new THnSparseD(Form("fhnV0ALambdaInclMCFD_%d", i), Form("MC ALambda associated, from AXi: pt-pt, cent %s;#it{p}_{T}^{A#Lambda,gen.} (GeV/#it{c});#it{p}_{T}^{A#Xi,gen.} (GeV/#it{c});#it{p}_{T}^{jet} (GeV/#it{c})", GetCentBinLabel(i).Data()), iNDimFD, binsFD, xminFD, xmaxFD);
      fOutputListMC->Add(fhnV0ALambdaInclMCFD[i]);
      fhnV0ALambdaInJetsMCFD[i] = new THnSparseD(Form("fhnV0ALambdaInJetsMCFD_%d", i), Form("MC ALambda associated, in JC, from AXi: pt-pt-ptJet, cent %s;#it{p}_{T}^{A#Lambda,gen.} (GeV/#it{c});#it{p}_{T}^{A#Xi,gen.} (GeV/#it{c});#it{p}_{T}^{jet} (GeV/#it{c})", GetCentBinLabel(i).Data()), iNDimFD, binsFD, xminFD, xmaxFD);
      fOutputListMC->Add(fhnV0ALambdaInJetsMCFD[i]);
      fhnV0ALambdaBulkMCFD[i] = new THnSparseD(Form("fhnV0ALambdaBulkMCFD_%d", i), Form("MC ALambda associated, in no jet events, from AXi: pt-pt-ptJet, cent %s;#it{p}_{T}^{A#Lambda,gen.} (GeV/#it{c});#it{p}_{T}^{A#Xi,gen.} (GeV/#it{c});#it{p}_{T}^{jet} (GeV/#it{c})", GetCentBinLabel(i).Data()), iNDimFD, binsFD, xminFD, xmaxFD);
      fOutputListMC->Add(fhnV0ALambdaBulkMCFD[i]);
      fh1V0AXiPtMCGen[i] = new TH1D(Form("fh1V0AXiPtMCGen_%d", i), Form("MC AXi^{-} generated: Pt spectrum, cent %s;#it{p}_{T}^{A#Xi^{-},gen.} (GeV/#it{c})", GetCentBinLabel(i).Data()), iNBinsPtXi, dPtXiMin, dPtXiMax);
      fOutputListMC->Add(fh1V0AXiPtMCGen[i]);

      // daughter eta
//          fhnV0K0sInclDaughterEtaPtPtMCGen[i] = new THnSparseD(Form("fhnV0K0sInclDaughterEtaPtPtMCGen_%d",i),Form("MC K0S, inclusive, gen., daughters: charge-etaD-ptD-etaV0-ptV0-ptJet, cent: %s;charge;eta daughter;pT daughter;eta V0;pT V0;pT jet",GetCentBinLabel(i).Data()),iNDimEtaD,binsEtaDaughter,xminEtaDaughter,xmaxEtaDaughter);
      fhnV0K0sInclDaughterEtaPtPtMCRec[i] = new THnSparseD(Form("fhnV0K0sInclDaughterEtaPtPtMCRec_%d", i), Form("MC K0S, inclusive, assoc., daughters: charge-etaD-ptD-etaV0-ptV0-ptJet, cent: %s;charge;eta daughter;pT daughter;eta V0;pT V0;pT jet", GetCentBinLabel(i).Data()), iNDimEtaD, binsEtaDaughter, xminEtaDaughter, xmaxEtaDaughter);
//          fhnV0K0sInJetsDaughterEtaPtPtMCGen[i] = new THnSparseD(Form("fhnV0K0sInJetsDaughterEtaPtPtMCGen_%d",i),Form("MC K0S, in JC, gen., daughters: charge-etaD-ptD-etaV0-ptV0-ptJet, cent: %s;charge;eta daughter;pT daughter;eta V0;pT V0;pT jet",GetCentBinLabel(i).Data()),iNDimEtaD,binsEtaDaughter,xminEtaDaughter,xmaxEtaDaughter);
      fhnV0K0sInJetsDaughterEtaPtPtMCRec[i] = new THnSparseD(Form("fhnV0K0sInJetsDaughterEtaPtPtMCRec_%d", i), Form("MC K0S, in JC, assoc., daughters: charge-etaD-ptD-etaV0-ptV0-ptJet, cent: %s;charge;eta daughter;pT daughter;eta V0;pT V0;pT jet", GetCentBinLabel(i).Data()), iNDimEtaD, binsEtaDaughter, xminEtaDaughter, xmaxEtaDaughter);
//          fhnV0LambdaInclDaughterEtaPtPtMCGen[i] = new THnSparseD(Form("fhnV0LambdaInclDaughterEtaPtPtMCGen_%d",i),Form("MC Lambda, inclusive, gen., daughters: charge-etaD-ptD-etaV0-ptV0-ptJet, cent: %s;charge;eta daughter;pT daughter;eta V0;pT V0;pT jet",GetCentBinLabel(i).Data()),iNDimEtaD,binsEtaDaughter,xminEtaDaughter,xmaxEtaDaughter);
      fhnV0LambdaInclDaughterEtaPtPtMCRec[i] = new THnSparseD(Form("fhnV0LambdaInclDaughterEtaPtPtMCRec_%d", i), Form("MC Lambda, inclusive, assoc., daughters: charge-etaD-ptD-etaV0-ptV0-ptJet, cent: %s;charge;eta daughter;pT daughter;eta V0;pT V0;pT jet", GetCentBinLabel(i).Data()), iNDimEtaD, binsEtaDaughter, xminEtaDaughter, xmaxEtaDaughter);
//          fhnV0LambdaInJetsDaughterEtaPtPtMCGen[i] = new THnSparseD(Form("fhnV0LambdaInJetsDaughterEtaPtPtMCGen_%d",i),Form("MC Lambda, in JC, gen., daughters: charge-etaD-ptD-etaV0-ptV0-ptJet, cent: %s;charge;eta daughter;pT daughter;eta V0;pT V0;pT jet",GetCentBinLabel(i).Data()),iNDimEtaD,binsEtaDaughter,xminEtaDaughter,xmaxEtaDaughter);
      fhnV0LambdaInJetsDaughterEtaPtPtMCRec[i] = new THnSparseD(Form("fhnV0LambdaInJetsDaughterEtaPtPtMCRec_%d", i), Form("MC Lambda, in JC, assoc., daughters: charge-etaD-ptD-etaV0-ptV0-ptJet, cent: %s;charge;eta daughter;pT daughter;eta V0;pT V0;pT jet", GetCentBinLabel(i).Data()), iNDimEtaD, binsEtaDaughter, xminEtaDaughter, xmaxEtaDaughter);
//          fhnV0ALambdaInclDaughterEtaPtPtMCGen[i] = new THnSparseD(Form("fhnV0ALambdaInclDaughterEtaPtPtMCGen_%d",i),Form("MC ALambda, inclusive, gen., daughters: charge-etaD-ptD-etaV0-ptV0-ptJet, cent: %s;charge;eta daughter;pT daughter;eta V0;pT V0;pT jet",GetCentBinLabel(i).Data()),iNDimEtaD,binsEtaDaughter,xminEtaDaughter,xmaxEtaDaughter);
      fhnV0ALambdaInclDaughterEtaPtPtMCRec[i] = new THnSparseD(Form("fhnV0ALambdaInclDaughterEtaPtPtMCRec_%d", i), Form("MC ALambda, inclusive, assoc., daughters: charge-etaD-ptD-etaV0-ptV0-ptJet, cent: %s;charge;eta daughter;pT daughter;eta V0;pT V0;pT jet", GetCentBinLabel(i).Data()), iNDimEtaD, binsEtaDaughter, xminEtaDaughter, xmaxEtaDaughter);
//          fhnV0ALambdaInJetsDaughterEtaPtPtMCGen[i] = new THnSparseD(Form("fhnV0ALambdaInJetsDaughterEtaPtPtMCGen_%d",i),Form("MC ALambda, in JC, gen., daughters: charge-etaD-ptD-etaV0-ptV0-ptJet, cent: %s;charge;eta daughter;pT daughter;eta V0;pT V0;pT jet",GetCentBinLabel(i).Data()),iNDimEtaD,binsEtaDaughter,xminEtaDaughter,xmaxEtaDaughter);
      fhnV0ALambdaInJetsDaughterEtaPtPtMCRec[i] = new THnSparseD(Form("fhnV0ALambdaInJetsDaughterEtaPtPtMCRec_%d", i), Form("MC ALambda, in JC, assoc., daughters: charge-etaD-ptD-etaV0-ptV0-ptJet, cent: %s;charge;eta daughter;pT daughter;eta V0;pT V0;pT jet", GetCentBinLabel(i).Data()), iNDimEtaD, binsEtaDaughter, xminEtaDaughter, xmaxEtaDaughter);

//          fOutputListMC->Add(fhnV0K0sInclDaughterEtaPtPtMCGen[i]);
      fOutputListMC->Add(fhnV0K0sInclDaughterEtaPtPtMCRec[i]);
//          fOutputListMC->Add(fhnV0K0sInJetsDaughterEtaPtPtMCGen[i]);
      fOutputListMC->Add(fhnV0K0sInJetsDaughterEtaPtPtMCRec[i]);
//          fOutputListMC->Add(fhnV0LambdaInclDaughterEtaPtPtMCGen[i]);
      fOutputListMC->Add(fhnV0LambdaInclDaughterEtaPtPtMCRec[i]);
//          fOutputListMC->Add(fhnV0LambdaInJetsDaughterEtaPtPtMCGen[i]);
      fOutputListMC->Add(fhnV0LambdaInJetsDaughterEtaPtPtMCRec[i]);
//          fOutputListMC->Add(fhnV0ALambdaInclDaughterEtaPtPtMCGen[i]);
      fOutputListMC->Add(fhnV0ALambdaInclDaughterEtaPtPtMCRec[i]);
//          fOutputListMC->Add(fhnV0ALambdaInJetsDaughterEtaPtPtMCGen[i]);
      fOutputListMC->Add(fhnV0ALambdaInJetsDaughterEtaPtPtMCRec[i]);
    }
  }

  // QA Histograms
  for(Int_t i = 0; i < fgkiNQAIndeces; i++)
  {
//    [i] = new TH1D(Form("%d",i),";;Counts",,,);
    fh1QAV0Status[i] = new TH1D(Form("fh1QAV0Status_%d", i), "QA: V0 status", 2, 0, 2);
    fh1QAV0TPCRefit[i] = new TH1D(Form("fh1QAV0TPCRefit_%d", i), "QA: TPC refit", 2, 0, 2);
    fh1QAV0TPCRows[i] = new TH1D(Form("fh1QAV0TPCRows_%d", i), "QA: TPC Rows", 160, 0, 160);
    fh1QAV0TPCFindable[i] = new TH1D(Form("fh1QAV0TPCFindable_%d", i), "QA: TPC Findable", 160, 0, 160);
    fh1QAV0TPCRowsFind[i] = new TH1D(Form("fh1QAV0TPCRowsFind_%d", i), "QA: TPC Rows/Findable", 100, 0, 2);
    fh1QAV0Eta[i] = new TH1D(Form("fh1QAV0Eta_%d", i), "QA: Daughter Eta", 200, -2, 2);
    fh2QAV0EtaRows[i] = new TH2D(Form("fh2QAV0EtaRows_%d", i), "QA: Daughter Eta vs TPC rows;#eta;TPC rows", 200, -2, 2, 160, 0, 160);
    fh2QAV0PtRows[i] = new TH2D(Form("fh2QAV0PtRows_%d", i), "QA: Daughter Pt vs TPC rows;pt;TPC rows", 100, 0, 10, 160, 0, 160);
    fh2QAV0PhiRows[i] = new TH2D(Form("fh2QAV0PhiRows_%d", i), "QA: Daughter Phi vs TPC rows;#phi;TPC rows", 100, 0, TMath::TwoPi(), 160, 0, 160);
    fh2QAV0NClRows[i] = new TH2D(Form("fh2QAV0NClRows_%d", i), "QA: Daughter NCl vs TPC rows;findable clusters;TPC rows", 100, 0, 160, 160, 0, 160);
    fh2QAV0EtaNCl[i] = new TH2D(Form("fh2QAV0EtaNCl_%d", i), "QA: Daughter Eta vs NCl;#eta;findable clusters", 200, -2, 2, 160, 0, 160);

    fh2QAV0EtaPtK0sPeak[i] = new TH2D(Form("fh2QAV0EtaPtK0sPeak_%d", i), "QA: K0s: Daughter Eta vs V0 pt, peak;track eta;V0 pt", 200, -2, 2, iNBinsPtV0, dPtV0Min, dPtV0Max);
    fh2QAV0EtaEtaK0s[i] = new TH2D(Form("fh2QAV0EtaEtaK0s_%d", i), "QA: K0s: Eta vs Eta Daughter", 200, -2, 2, 200, -2, 2);
    fh2QAV0PhiPhiK0s[i] = new TH2D(Form("fh2QAV0PhiPhiK0s_%d", i), "QA: K0s: Phi vs Phi Daughter", 200, 0, TMath::TwoPi(), 200, 0, TMath::TwoPi());
    fh1QAV0RapK0s[i] = new TH1D(Form("fh1QAV0RapK0s_%d", i), "QA: K0s: V0 Rapidity", 200, -2, 2);
    fh2QAV0PtPtK0sPeak[i] = new TH2D(Form("fh2QAV0PtPtK0sPeak_%d", i), "QA: K0s: Daughter Pt vs Pt;neg pt;pos pt", 100, 0, 5, 100, 0, 5);

    fh2QAV0EtaPtLambdaPeak[i] = new TH2D(Form("fh2QAV0EtaPtLambdaPeak_%d", i), "QA: Lambda: Daughter Eta vs V0 pt, peak;track eta;V0 pt", 200, -2, 2, iNBinsPtV0, dPtV0Min, dPtV0Max);
    fh2QAV0EtaEtaLambda[i] = new TH2D(Form("fh2QAV0EtaEtaLambda_%d", i), "QA: Lambda: Eta vs Eta Daughter", 200, -2, 2, 200, -2, 2);
    fh2QAV0PhiPhiLambda[i] = new TH2D(Form("fh2QAV0PhiPhiLambda_%d", i), "QA: Lambda: Phi vs Phi Daughter", 200, 0, TMath::TwoPi(), 200, 0, TMath::TwoPi());
    fh1QAV0RapLambda[i] = new TH1D(Form("fh1QAV0RapLambda_%d", i), "QA: Lambda: V0 Rapidity", 200, -2, 2);
    fh2QAV0PtPtLambdaPeak[i] = new TH2D(Form("fh2QAV0PtPtLambdaPeak_%d", i), "QA: Lambda: Daughter Pt vs Pt;neg pt;pos pt", 100, 0, 5, 100, 0, 5);

    fh1QAV0Pt[i] = new TH1D(Form("fh1QAV0Pt_%d", i), "QA: Daughter Pt", 100, 0, 5);
    fh1QAV0Charge[i] = new TH1D(Form("fh1QAV0Charge_%d", i), "QA: V0 Charge", 3, -1, 2);
    fh1QAV0DCAVtx[i] = new TH1D(Form("fh1QAV0DCAVtx_%d", i), "QA: DCA daughters to primary vertex", 100, 0, 10);
    fh1QAV0DCAV0[i] = new TH1D(Form("fh1QAV0DCAV0_%d", i), "QA: DCA daughters", 100, 0, 2);
    fh1QAV0Cos[i] = new TH1D(Form("fh1QAV0Cos_%d", i), "QA: CPA", 10000, 0.9, 1);
    fh1QAV0R[i] = new TH1D(Form("fh1QAV0R_%d", i), "QA: R", 1500, 0, 150);
    fh1QACTau2D[i] = new TH1D(Form("fh1QACTau2D_%d", i), "QA: K0s: c#tau 2D;mR/pt#tau", 100, 0, 10);
    fh1QACTau3D[i] = new TH1D(Form("fh1QACTau3D_%d", i), "QA: K0s: c#tau 3D;mL/p#tau", 100, 0, 10);

    fh2ArmPod[i] = new TH2D(Form("fh2ArmPod_%d", i), "Armenteros-Podolanski;#alpha;#it{p}_{T}^{Arm}", 100, -1., 1., 50, 0., 0.25);
    fh2ArmPodK0s[i] = new TH2D(Form("fh2ArmPodK0s_%d", i), "K0s: Armenteros-Podolanski;#alpha;#it{p}_{T}^{Arm}", 100, -1., 1., 50, 0., 0.25);
    fh2ArmPodLambda[i] = new TH2D(Form("fh2ArmPodLambda_%d", i), "Lambda: Armenteros-Podolanski;#alpha;#it{p}_{T}^{Arm}", 100, -1., 1., 50, 0., 0.25);
    fh2ArmPodALambda[i] = new TH2D(Form("fh2ArmPodALambda_%d", i), "ALambda: Armenteros-Podolanski;#alpha;#it{p}_{T}^{Arm}", 100, -1., 1., 50, 0., 0.25);

    fOutputListQA->Add(fh1QAV0Status[i]);
    fOutputListQA->Add(fh1QAV0TPCRefit[i]);
    fOutputListQA->Add(fh1QAV0TPCRows[i]);
    fOutputListQA->Add(fh1QAV0TPCFindable[i]);
    fOutputListQA->Add(fh1QAV0TPCRowsFind[i]);
    fOutputListQA->Add(fh1QAV0Eta[i]);
    fOutputListQA->Add(fh2QAV0EtaRows[i]);
    fOutputListQA->Add(fh2QAV0PtRows[i]);
    fOutputListQA->Add(fh2QAV0PhiRows[i]);
    fOutputListQA->Add(fh2QAV0NClRows[i]);
    fOutputListQA->Add(fh2QAV0EtaNCl[i]);

    fOutputListQA->Add(fh2QAV0EtaPtK0sPeak[i]);
    fOutputListQA->Add(fh2QAV0EtaEtaK0s[i]);
    fOutputListQA->Add(fh2QAV0PhiPhiK0s[i]);
    fOutputListQA->Add(fh1QAV0RapK0s[i]);
    fOutputListQA->Add(fh2QAV0PtPtK0sPeak[i]);

    fOutputListQA->Add(fh2QAV0EtaPtLambdaPeak[i]);
    fOutputListQA->Add(fh2QAV0EtaEtaLambda[i]);
    fOutputListQA->Add(fh2QAV0PhiPhiLambda[i]);
    fOutputListQA->Add(fh1QAV0RapLambda[i]);
    fOutputListQA->Add(fh2QAV0PtPtLambdaPeak[i]);

    fOutputListQA->Add(fh1QAV0Pt[i]);
    fOutputListQA->Add(fh1QAV0Charge[i]);
    fOutputListQA->Add(fh1QAV0DCAVtx[i]);
    fOutputListQA->Add(fh1QAV0DCAV0[i]);
    fOutputListQA->Add(fh1QAV0Cos[i]);
    fOutputListQA->Add(fh1QAV0R[i]);
    fOutputListQA->Add(fh1QACTau2D[i]);
    fOutputListQA->Add(fh1QACTau3D[i]);

    fOutputListQA->Add(fh2ArmPod[i]);
    fOutputListQA->Add(fh2ArmPodK0s[i]);
    fOutputListQA->Add(fh2ArmPodLambda[i]);
    fOutputListQA->Add(fh2ArmPodALambda[i]);

    /*
    fh2CutTPCRowsK0s[i] = new TH2D(Form("fh2CutTPCRowsK0s_%d", i), "Cuts: K0s: TPC Rows vs mass;#it{m}_{inv} (GeV/#it{c}^{2});TPC rows", fgkiNBinsMassK0s, fgkdMassK0sMin, fgkdMassK0sMax, 160, 0, 160);
    fh2CutTPCRowsLambda[i] = new TH2D(Form("fh2CutTPCRowsLambda_%d", i), "Cuts: Lambda: TPC Rows vs mass;#it{m}_{inv} (GeV/#it{c}^{2});TPC rows", fgkiNBinsMassLambda, fgkdMassLambdaMin, fgkdMassLambdaMax, 160, 0, 160);
    fh2CutPtPosK0s[i] = new TH2D(Form("fh2CutPtPosK0s_%d", i), "Cuts: K0s: Pt pos;#it{m}_{inv} (GeV/#it{c}^{2});pt pos", fgkiNBinsMassK0s, fgkdMassK0sMin, fgkdMassK0sMax, 100, 0, 5);
    fh2CutPtNegK0s[i] = new TH2D(Form("fh2CutPtNegK0s_%d", i), "Cuts: K0s: Pt neg;#it{m}_{inv} (GeV/#it{c}^{2});pt neg", fgkiNBinsMassK0s, fgkdMassK0sMin, fgkdMassK0sMax, 100, 0, 5);
    fh2CutPtPosLambda[i] = new TH2D(Form("fh2CutPtPosLambda_%d", i), "Cuts: Lambda: Pt pos;#it{m}_{inv} (GeV/#it{c}^{2});pt pos", fgkiNBinsMassLambda, fgkdMassLambdaMin, fgkdMassLambdaMax, 100, 0, 5);
    fh2CutPtNegLambda[i] = new TH2D(Form("fh2CutPtNegLambda_%d", i), "Cuts: Lambda: Pt neg;#it{m}_{inv} (GeV/#it{c}^{2});pt neg", fgkiNBinsMassLambda, fgkdMassLambdaMin, fgkdMassLambdaMax, 100, 0, 5);
    fh2CutDCAVtx[i] = new TH2D(Form("fh2CutDCAVtx_%d", i), "Cuts: DCA daughters to prim. vtx.;#it{m}_{inv} (GeV/#it{c}^{2});DCA daughter to prim. vtx. (cm)", fgkiNBinsMassK0s, fgkdMassK0sMin, fgkdMassK0sMax, 100, 0, 10);
    fh2CutDCAV0[i] = new TH2D(Form("fh2CutDCAV0_%d", i), "Cuts: DCA daughters;#it{m}_{inv} (GeV/#it{c}^{2});DCA daughters / #sigma_{TPC}", fgkiNBinsMassK0s, fgkdMassK0sMin, fgkdMassK0sMax, 100, 0, 2);
    fh2CutCos[i] = new TH2D(Form("fh2CutCos_%d", i), "Cuts: CPA;#it{m}_{inv} (GeV/#it{c}^{2});CPA", fgkiNBinsMassK0s, fgkdMassK0sMin, fgkdMassK0sMax, 10000, 0.9, 1);
    fh2CutR[i] = new TH2D(Form("fh2CutR_%d", i), "Cuts: R;#it{m}_{inv} (GeV/#it{c}^{2});R (cm)", fgkiNBinsMassK0s, fgkdMassK0sMin, fgkdMassK0sMax, 1500, 0, 150);
    fh2CutEtaK0s[i] = new TH2D(Form("fh2CutEtaK0s_%d", i), "Cuts: K0s: Eta;#it{m}_{inv} (GeV/#it{c}^{2});#eta", fgkiNBinsMassK0s, fgkdMassK0sMin, fgkdMassK0sMax, 200, -2, 2);
    fh2CutEtaLambda[i] = new TH2D(Form("fh2CutEtaLambda_%d", i), "Cuts: Lambda: Eta;#it{m}_{inv} (GeV/#it{c}^{2});#eta", fgkiNBinsMassLambda, fgkdMassLambdaMin, fgkdMassLambdaMax, 200, -2, 2);
    fh2CutRapK0s[i] = new TH2D(Form("fh2CutRapK0s_%d", i), "Cuts: K0s: Rapidity;#it{m}_{inv} (GeV/#it{c}^{2});y", fgkiNBinsMassK0s, fgkdMassK0sMin, fgkdMassK0sMax, 200, -2, 2);
    fh2CutRapLambda[i] = new TH2D(Form("fh2CutRapLambda_%d", i), "Cuts: Lambda: Rapidity;#it{m}_{inv} (GeV/#it{c}^{2});y", fgkiNBinsMassLambda, fgkdMassLambdaMin, fgkdMassLambdaMax, 200, -2, 2);
    fh2CutCTauK0s[i] = new TH2D(Form("fh2CutCTauK0s_%d", i), "Cuts: K0s: #it{c#tau};#it{m}_{inv} (GeV/#it{c}^{2});#it{mL/p#tau}", fgkiNBinsMassK0s, fgkdMassK0sMin, fgkdMassK0sMax, 100, 0, 10);
    fh2CutCTauLambda[i] = new TH2D(Form("fh2CutCTauLambda_%d", i), "Cuts: Lambda: #it{c#tau};#it{m}_{inv} (GeV/#it{c}^{2});#it{mL/p#tau}", fgkiNBinsMassLambda, fgkdMassLambdaMin, fgkdMassLambdaMax, 100, 0, 10);
    fh2CutPIDPosK0s[i] = new TH2D(Form("fh2CutPIDPosK0s_%d", i), "Cuts: K0s: PID pos;#it{m}_{inv} (GeV/#it{c}^{2});##sigma_{d#it{E}/d#it{x}}", fgkiNBinsMassK0s, fgkdMassK0sMin, fgkdMassK0sMax, 100, 0, 10);
    fh2CutPIDNegK0s[i] = new TH2D(Form("fh2CutPIDNegK0s_%d", i), "Cuts: K0s: PID neg;#it{m}_{inv} (GeV/#it{c}^{2});##sigma_{d#it{E}/d#it{x}}", fgkiNBinsMassK0s, fgkdMassK0sMin, fgkdMassK0sMax, 100, 0, 10);
    fh2CutPIDPosLambda[i] = new TH2D(Form("fh2CutPIDPosLambda_%d", i), "Cuts: Lambda: PID pos;#it{m}_{inv} (GeV/#it{c}^{2});##sigma_{d#it{E}/d#it{x}}", fgkiNBinsMassLambda, fgkdMassLambdaMin, fgkdMassLambdaMax, 100, 0, 10);
    fh2CutPIDNegLambda[i] = new TH2D(Form("fh2CutPIDNegLambda_%d", i), "Cuts: Lambda: PID neg;#it{m}_{inv} (GeV/#it{c}^{2});##sigma_{d#it{E}/d#it{x}}", fgkiNBinsMassLambda, fgkdMassLambdaMin, fgkdMassLambdaMax, 100, 0, 10);
    fh2Tau3DVs2D[i] = new TH2D(Form("fh2Tau3DVs2D_%d", i), "Decay 3D vs 2D;pt;3D/2D", 100, 0, 10, 200, 0.5, 1.5);

    fOutputListCuts->Add(fh2CutTPCRowsK0s[i]);
    fOutputListCuts->Add(fh2CutTPCRowsLambda[i]);
    fOutputListCuts->Add(fh2CutPtPosK0s[i]);
    fOutputListCuts->Add(fh2CutPtNegK0s[i]);
    fOutputListCuts->Add(fh2CutPtPosLambda[i]);
    fOutputListCuts->Add(fh2CutPtNegLambda[i]);
    fOutputListCuts->Add(fh2CutDCAVtx[i]);
    fOutputListCuts->Add(fh2CutDCAV0[i]);
    fOutputListCuts->Add(fh2CutCos[i]);
    fOutputListCuts->Add(fh2CutR[i]);
    fOutputListCuts->Add(fh2CutEtaK0s[i]);
    fOutputListCuts->Add(fh2CutEtaLambda[i]);
    fOutputListCuts->Add(fh2CutRapK0s[i]);
    fOutputListCuts->Add(fh2CutRapLambda[i]);
    fOutputListCuts->Add(fh2CutCTauK0s[i]);
    fOutputListCuts->Add(fh2CutCTauLambda[i]);
    fOutputListCuts->Add(fh2CutPIDPosK0s[i]);
    fOutputListCuts->Add(fh2CutPIDNegK0s[i]);
    fOutputListCuts->Add(fh2CutPIDPosLambda[i]);
    fOutputListCuts->Add(fh2CutPIDNegLambda[i]);
    fOutputListCuts->Add(fh2Tau3DVs2D[i]);
    */
  }

  for(Int_t i = 0; i < fgkiNCategV0; i++)
  {
    fh1V0InvMassK0sAll[i] = new TH1D(Form("fh1V0InvMassK0sAll_%d", i), Form("K0s: V0 invariant mass, %s;#it{m}_{inv} (GeV/#it{c}^{2});counts", categV0[i].Data()), fgkiNBinsMassK0s, fgkdMassK0sMin, fgkdMassK0sMax);
    fh1V0InvMassLambdaAll[i] = new TH1D(Form("fh1V0InvMassLambdaAll_%d", i), Form("Lambda: V0 invariant mass, %s;#it{m}_{inv} (GeV/#it{c}^{2});counts", categV0[i].Data()), fgkiNBinsMassLambda, fgkdMassLambdaMin, fgkdMassLambdaMax);
    fh1V0InvMassALambdaAll[i] = new TH1D(Form("fh1V0InvMassALambdaAll_%d", i), Form("ALambda: V0 invariant mass, %s;#it{m}_{inv} (GeV/#it{c}^{2});counts", categV0[i].Data()), fgkiNBinsMassLambda, fgkdMassLambdaMin, fgkdMassLambdaMax);
    fOutputListStd->Add(fh1V0InvMassK0sAll[i]);
    fOutputListStd->Add(fh1V0InvMassLambdaAll[i]);
    fOutputListStd->Add(fh1V0InvMassALambdaAll[i]);
  }

  for(Int_t i = 0; i < fOutputListStd->GetEntries(); ++i)
  {
    TH1* h1 = dynamic_cast<TH1*>(fOutputListStd->At(i));
    if(h1)
    {
      h1->Sumw2();
      continue;
    }
    THnSparse* hn = dynamic_cast<THnSparse*>(fOutputListStd->At(i));
    if(hn) hn->Sumw2();
  }
  for(Int_t i = 0; i < fOutputListQA->GetEntries(); ++i)
  {
    TH1* h1 = dynamic_cast<TH1*>(fOutputListQA->At(i));
    if(h1)
    {
      h1->Sumw2();
      continue;
    }
    THnSparse* hn = dynamic_cast<THnSparse*>(fOutputListQA->At(i));
    if(hn) hn->Sumw2();
  }
  for(Int_t i = 0; i < fOutputListCuts->GetEntries(); ++i)
  {
    TH1* h1 = dynamic_cast<TH1*>(fOutputListCuts->At(i));
    if(h1)
    {
      h1->Sumw2();
      continue;
    }
    THnSparse* hn = dynamic_cast<THnSparse*>(fOutputListCuts->At(i));
    if(hn) hn->Sumw2();
  }
  for(Int_t i = 0; i < fOutputListMC->GetEntries(); ++i)
  {
    TH1* h1 = dynamic_cast<TH1*>(fOutputListMC->At(i));
    if(h1)
    {
      h1->Sumw2();
      continue;
    }
    THnSparse* hn = dynamic_cast<THnSparse*>(fOutputListMC->At(i));
    if(hn) hn->Sumw2();
  }

  PostData(1, fOutputListStd);
  PostData(2, fOutputListQA);
  PostData(3, fOutputListCuts);
  PostData(4, fOutputListMC);
//  if (fbTreeOutput)
//    PostData(5,ftreeOut);
}

void AliAnalysisTaskV0sInJets::UserExec(Option_t*)
{
  // Main loop, called for each event
  if(fDebug > 5) printf("TaskV0sInJets: UserExec: Start\n");
/*
  // reset branches for each event
  if (fBranchV0Rec)
    fBranchV0Rec->Clear();
  if (fBranchV0Gen)
    fBranchV0Gen->Clear();
  if (fBranchJet)
    fBranchJet->Clear();
  if (fEventInfo)
    fEventInfo->Reset();
*/
  if(!fiAODAnalysis)
    return;

  if(fDebug > 2) printf("TaskV0sInJets: AOD analysis\n");
  fh1EventCounterCut->Fill(0); // all available selected events (collision candidates)

  if(fDebug > 5) printf("TaskV0sInJets: UserExec: Loading AOD\n");
  fAODIn = dynamic_cast<AliAODEvent*>(InputEvent()); // input AOD
  fAODOut = AODEvent(); // output AOD
  if(!fAODOut)
  {
    if(fDebug > 0) printf("TaskV0sInJets: No output AOD found\n");
    return;
  }
  if(!fAODIn)
  {
    if(fDebug > 0) printf("TaskV0sInJets: No input AOD found\n");
    return;
  }
  if(fDebug > 5) printf("TaskV0sInJets: UserExec: Loading AOD OK\n");

  TClonesArray* arrayMC = 0; // array particles in the MC event
  AliAODMCHeader* headerMC = 0; // MC header
  Int_t iNTracksMC = 0; // number of MC tracks
  Double_t dPrimVtxMCX = 0., dPrimVtxMCY = 0., dPrimVtxMCZ = 0.; // position of the MC primary vertex

  // Simulation info
  if(fbMCAnalysis)
  {
    arrayMC = (TClonesArray*)fAODIn->FindListObject(AliAODMCParticle::StdBranchName());
    if(!arrayMC)
    {
      if(fDebug > 0) printf("TaskV0sInJets: No MC array found\n");
      return;
    }
    if(fDebug > 5) printf("TaskV0sInJets: MC array found\n");
    iNTracksMC = arrayMC->GetEntriesFast();
    if(fDebug > 5) printf("TaskV0sInJets: There are %d MC tracks in this event\n", iNTracksMC);
//      if (!iNTracksMC)
//        return;
    headerMC = (AliAODMCHeader*)fAODIn->FindListObject(AliAODMCHeader::StdBranchName());
    if(!headerMC)
    {
      if(fDebug > 0) printf("TaskV0sInJets: No MC header found\n");
      return;
    }
    // get position of the MC primary vertex
    dPrimVtxMCX = headerMC->GetVtxX();
    dPrimVtxMCY = headerMC->GetVtxY();
    dPrimVtxMCZ = headerMC->GetVtxZ();
  }

  // PID Response Task object
  AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = (AliInputEventHandler*)mgr->GetInputEventHandler();
  AliPIDResponse* fPIDResponse = inputHandler->GetPIDResponse();
  if(!fPIDResponse)
  {
    if(fDebug > 0) printf("TaskV0sInJets: No PID response object found\n");
    return;
  }

  // AOD files are OK
  fh1EventCounterCut->Fill(1);

  // Event selection
  if(!IsSelectedForJets(fAODIn, fdCutVertexZ, fdCutVertexR2, fdCutCentLow, fdCutCentHigh, 1, 0.1)) // cut on |delta z| in  2011 data between SPD vertex and nominal primary vertex
//  if (!IsSelectedForJets(fAODIn,fdCutVertexZ,fdCutVertexR2,fdCutCentLow,fdCutCentHigh)) // no need for cutting in 2010 data
  {
    if(fDebug > 5) printf("TaskV0sInJets: Event rejected\n");
    return;
  }

//  fdCentrality = ((AliVAODHeader*)fAODIn->GetHeader())->GetCentrality(); // event centrality
  fdCentrality = ((AliVAODHeader*)fAODIn->GetHeader())->GetCentralityP()->GetCentralityPercentile("V0M"); // event centrality
  if(!fbIsPbPb)
    fdCentrality = 0.;
  Int_t iCentIndex = GetCentralityBinIndex(fdCentrality); // get index of centrality bin
  if(iCentIndex < 0)
  {
    if(fDebug > 5) printf("TaskV0sInJets: Event is out of histogram range\n");
    return;
  }
  fh1EventCounterCut->Fill(2); // selected events (vertex, centrality)
  fh1EventCounterCutCent[iCentIndex]->Fill(2);

  UInt_t iNTracks = fAODIn->GetNumberOfTracks(); // get number of tracks in event
  if(fDebug > 5) printf("TaskV0sInJets: There are %d tracks in this event\n", iNTracks);
//  if (!iNTracks)
//    return;

  Int_t iNV0s = fAODIn->GetNumberOfV0s(); // get the number of V0 candidates
  if(!iNV0s)
  {
    if(fDebug > 2) printf("TaskV0sInJets: No V0s found in event\n");
//      return;
  }

  //===== Event is OK for the analysis =====
  fh1EventCent->Fill(iCentIndex);
  fh1EventCent2->Fill(fdCentrality);
  fh2EventCentTracks->Fill(fdCentrality, iNTracks);

//  if (fbTreeOutput)
//    fEventInfo->SetAODEvent(fAODIn);
//  printf("V0sInJets: EventInfo: Centrality: %f\n",fEventInfo->GetCentrality());

  if(iNV0s)
  {
    fh1EventCounterCut->Fill(3); // events with V0s
    fh1EventCounterCutCent[iCentIndex]->Fill(3);
  }

//  Int_t iNV0SelV0Rec = 0;
//  Int_t iNV0SelV0Gen = 0;

  AliAODv0* v0 = 0; // pointer to V0 candidates
//  AliV0Object* objectV0 = 0;
  TVector3 vecV0Momentum; // 3D vector of V0 momentum
  Double_t dMassV0K0s = 0; // invariant mass of the K0s candidate
  Double_t dMassV0Lambda = 0; // invariant mass of the Lambda candidate
  Double_t dMassV0ALambda = 0; // invariant mass of the Lambda candidate
  Int_t iNV0CandTot = 0; // counter of all V0 candidates at the beginning
  Int_t iNV0CandK0s = 0; // counter of K0s candidates at the end
  Int_t iNV0CandLambda = 0; // counter of Lambda candidates at the end
  Int_t iNV0CandALambda = 0; // counter of Lambda candidates at the end

  Bool_t bUseOldCuts = 0; // old reconstruction cuts
  Bool_t bUseAliceCuts = 0; // cuts used by Alice Zimmermann
  Bool_t bUseIouriCuts = 0; // cuts used by Iouri
  Bool_t bPrintCuts = 0; // print out which cuts are applied
  Bool_t bPrintJetSelection = 0; // print out which jets are selected

  // Values of V0 reconstruction cuts:
  // Daughter tracks
  Int_t iRefit = AliAODTrack::kTPCrefit; // TPC refit for daughter tracks
  Double_t dDCAToPrimVtxMin = fdCutDCAToPrimVtxMin; // 0.1; // [cm] min DCA of daughters to the prim vtx
  Double_t dDCADaughtersMax = fdCutDCADaughtersMax; // 1.; // [sigma of TPC tracking] max DCA between daughters
  Double_t dEtaDaughterMax = 0.8; // max |pseudorapidity| of daughter tracks
  Double_t dNSigmadEdxMax = fdCutNSigmadEdxMax;// 3.; // [sigma dE/dx] max difference between measured and expected signal of dE/dx in the TPC
  Double_t dPtProtonPIDMax = 1.; // [GeV/c] maxium pT of proton for applying PID cut
  // V0 candidate
  Bool_t bOnFly = 0; // on-the-fly (yes) or offline (no) reconstructed
  Double_t dCPAMin = fdCutCPAMin;// 0.998; // min cosine of the pointing angle
  Double_t dRadiusDecayMin = 5.; // [cm] min radial distance of the decay vertex
  Double_t dRadiusDecayMax = 100.; // [cm] max radial distance of the decay vertex
  Double_t dEtaMax = 0.7; // max |pseudorapidity| of V0
  Double_t dNTauMax = fdCutNTauMax; // 5.0; // [tau] max proper lifetime in multiples of the mean lifetime

  // Old cuts Start
  Double_t dNCrossedRowsTPCMin = 70.; // min number of crossed TPC rows (turned off)
//  Double_t dCrossedRowsOverFindMin = 0.8; // min ratio crossed rows / findable clusters (turned off)
//  Double_t dCrossedRowsOverFindMax = 1e3; // max ratio crossed rows / findable clusters (turned off)
  Double_t dPtDaughterMin = 0.150; // [GeV/c] min transverse momentum of daughter tracks (turned off)
  Double_t dRapMax = 0.75; // max |rapidity| of V0 (turned off)
  // Old cuts End

  // Other cuts
  Double_t dNSigmaMassMax = 3.; // [sigma m] max difference between candidate mass and real particle mass (used only for mass peak method of signal extraction)
  Double_t dDistPrimaryMax = 0.01; // [cm] max distance of production point to the primary vertex (criterion for choice of MC particles considered as primary)

  // Selection of active cuts
  Bool_t bCutEtaDaughter = 1; // daughter pseudorapidity
  Bool_t bCutRapV0 = 0; // V0 rapidity
  Bool_t bCutEtaV0 = 1; // V0 pseudorapidity
  Bool_t bCutTau = 1; // V0 lifetime
  Bool_t bCutPid = 1; // PID (TPC dE/dx)
  Bool_t bCutArmPod = 1; // Armenteros-Podolanski for K0S
//  Bool_t bCutCross = 0; // cross contamination

  if(bUseOldCuts)
  {
    bCutRapV0 = 1;
    dEtaMax = 0.75;
    dNTauMax = 3.0;
  }
  else if(bUseAliceCuts)
  {
//      bOnFly = 1;
    dEtaMax = 0.75;
    dNTauMax = 5.0;
  }
  else if(bUseIouriCuts)
  {
    bCutRapV0 = 1;
    bCutEtaV0 = 0;
    dNTauMax = 3.0;
    dRapMax = 0.5;
  }

  Double_t dCTauK0s = 2.6844; // [cm] c tau of K0S
  Double_t dCTauLambda = 7.89; // [cm] c tau of Lambda

  // Load PDG values of particle masses
  Double_t dMassPDGK0s = TDatabasePDG::Instance()->GetParticle(kK0Short)->Mass();
  Double_t dMassPDGLambda = TDatabasePDG::Instance()->GetParticle(kLambda0)->Mass();

  // PDG codes of used particles
  Int_t iPdgCodePion = 211;
  Int_t iPdgCodeProton = 2212;
  Int_t iPdgCodeK0s = 310;
  Int_t iPdgCodeLambda = 3122;

  // Jet selection: fdCutPtJetMin, fdCutPtTrackMin
  Double_t dJetEtaWindow = dEtaMax - fdRadiusJet; // max jet |pseudorapidity|, to make sure that V0s can appear in the entire jet area
  Double_t dCutJetAreaMin = 0.6 * TMath::Pi() * fdRadiusJet * fdRadiusJet; // minimum jet area
  Double_t dRadiusExcludeCone = 2 * fdRadiusJet; // radius of cones around jets excluded for V0 outside jets
  Bool_t bLeadingJetOnly = 0;

  if(bUseAliceCuts)
  {
    fdCutPtJetMin = 5;
    fdCutPtTrackMin = 5;
    dCutJetAreaMin = 0;
    bLeadingJetOnly = 0;
  }

//  Int_t iNJetAll = 0; // number of reconstructed jets in fBranchJet
//  iNJetAll = fBranchJet->GetEntriesFast(); // number of reconstructed jets
  TClonesArray* jetArray = 0; // object where the input jets are stored
  TClonesArray* jetArrayBg = 0; // object where the kt clusters are stored
  Int_t iNJet = 0; // number of reconstructed jets in the input
  TClonesArray* jetArraySel = new TClonesArray("AliAODJet", 0); // object where the selected jets are copied
  Int_t iNJetSel = 0; // number of selected reconstructed jets
//  iNJetSel = jetArraySel->GetEntriesFast(); // number of selected reconstructed jets
  TClonesArray* jetArrayPerp = new TClonesArray("AliAODJet", 0); // object where the perp. cones are stored
  Int_t iNJetPerp = 0; // number of perpendicular cones

  AliAODJet* jet = 0; // pointer to a jet
//  AliJetObject* objectJet = 0;
  AliAODJet* jetPerp = 0; // pointer to a perp. cone
  AliAODJet* jetRnd = 0; // pointer to a rand. cone
  AliAODJet* jetMed = 0; // pointer to a median cluster
  TVector3 vecJetMomentum; // 3D vector of jet momentum
//  TVector3 vecPerpConeMomentum; // 3D vector of perpendicular cone momentum
//  TVector3 vecRndConeMomentum; // 3D vector of random cone momentum
  Bool_t bJetEventGood = kTRUE; // indicator of good jet events

//  printf("iNJetAll, iNJetSel: %d %d\n",iNJetAll,iNJetSel);

  if(fbJetSelection)  // analysis of V0s in jets is switched on
  {
    jetArray = (TClonesArray*)(fAODOut->FindListObject(fsJetBranchName.Data())); // find object with jets in the output AOD
    if(!jetArray)
    {
      if(fDebug > 0) printf("TaskV0sInJets: No array of name: %s\n", fsJetBranchName.Data());
      bJetEventGood = kFALSE;
    }
    if(bJetEventGood)
      iNJet = jetArray->GetEntriesFast();
    if(bJetEventGood && !iNJet)  // check whether there are some jets
    {
      if(fDebug > 2) printf("TaskV0sInJets: No jets in array\n");
      bJetEventGood = kFALSE;
    }
    if(bJetEventGood)
    {
//          printf("TaskV0sInJets: Loading bg array of name: %s\n",fsJetBgBranchName.Data());
      jetArrayBg = (TClonesArray*)(fAODOut->FindListObject(fsJetBgBranchName.Data())); // find object with jets in the output AOD
      if(!jetArrayBg)
      {
        if(fDebug > 0) printf("TaskV0sInJets: No bg array of name: %s\n", fsJetBgBranchName.Data());
//              bJetEventGood = kFALSE;
      }
    }
  }
  else // no in-jet analysis
    bJetEventGood = kFALSE;

  // select good jets and copy them to another array
  if(bJetEventGood)
  {
    if(bLeadingJetOnly)
      iNJet = 1; // only leading jets
    if(fDebug > 5) printf("TaskV0sInJets: Jet selection for %d jets\n", iNJet);
    for(Int_t iJet = 0; iJet < iNJet; iJet++)
    {
      AliAODJet* jetSel = (AliAODJet*)jetArray->At(iJet); // load a jet in the list
      if(!jetSel)
      {
        if(fDebug > 0) printf("TaskV0sInJets: Cannot load jet %d\n", iJet);
        continue;
      }
      if(bPrintJetSelection)
        if(fDebug > 7) printf("jet: i = %d, pT = %f, eta = %f, phi = %f, pt lead tr = %f ", iJet, jetSel->Pt(), jetSel->Eta(), jetSel->Phi(), jetSel->GetPtLeading());
//          printf("TaskV0sInJets: Checking pt > %.2f for jet %d with pt %.2f\n",fdCutPtJetMin,iJet,jetSel->Pt());
      if(jetSel->Pt() < fdCutPtJetMin)  // selection of high-pt jets
      {
        if(bPrintJetSelection)
          if(fDebug > 7) printf("rejected (pt)\n");
        continue;
      }
//          printf("TaskV0sInJets: Checking |eta| < %.2f for jet %d with |eta| %.2f\n",dJetEtaWindow,iJet,TMath::Abs(jetSel->Eta()));
      if(TMath::Abs(jetSel->Eta()) > dJetEtaWindow)  // selection of jets in the chosen pseudorapidity range
      {
        if(bPrintJetSelection)
          if(fDebug > 7) printf("rejected (eta)\n");
        continue;
      }
      if(!bUseOldCuts)
      {
        if(jetSel->EffectiveAreaCharged() < dCutJetAreaMin)
          continue;
      }
      Int_t iNTracksInJet = 0;
      Double_t dPtLeadTrack = 0; // pt of the leading track
//          Int_t iLeadTrack = 0;
      iNTracksInJet = jetSel->GetRefTracks()->GetEntriesFast(); // number od tracks that constitute the jet
//          printf("TaskV0sInJets: Searching for leading track from %d tracks in jet %d\n",iNTracksInJet,iJet);
      if(fdCutPtTrackMin > 0)  // a positive min leading track pt is set
      {
        for(Int_t j = 0; j < iNTracksInJet; j++)  // find the track with the highest pt
        {
          AliAODTrack* track = (AliAODTrack*)jetSel->GetTrack(j); // is this the leading track?
          if(!track)
            continue;
//                  printf("TaskV0sInJets: %d: %.2f\n",j,track->Pt());
          if(track->Pt() > dPtLeadTrack)
          {
            dPtLeadTrack = track->Pt();
//                      iLeadTrack = j;
          }
        }
//        printf("Leading track pT: my: %f, ali: %f\n",dPtLeadTrack,jetSel->GetPtLeading());
//        printf("TaskV0sInJets: Checking leading track pt > %.2f for pt %.2f of track %d in jet %d\n",fdCutPtTrackMin,dPtLeadTrack,iLeadTrack,iJet);
        if(dPtLeadTrack < fdCutPtTrackMin)  // selection of high-pt jet-track events
        {
          if(bPrintJetSelection)
            if(fDebug > 7) printf("rejected (track pt)\n");
          continue;
        }
      }
      if(bPrintJetSelection)
        if(fDebug > 7) printf("accepted\n");
      if(fDebug > 5) printf("TaskV0sInJets: Jet %d with pt %.2f passed selection\n", iJet, jetSel->Pt());
/*
      if (fbTreeOutput)
      {
//      new ((*fBranchJet)[iNJetAll++]) AliAODJet(*((AliAODJet*)jetSel));
        objectJet = new ((*fBranchJet)[iNJetAll++]) AliJetObject(jetSel); // copy selected jet to the array
//      objectJet->SetPtLeadingTrack(dPtLeadTrack);
        objectJet->SetRadius(fdRadiusJet);
        objectJet = 0;
      }
*/
      TLorentzVector vecPerpPlus(*(jetSel->MomentumVector()));
      vecPerpPlus.RotateZ(TMath::Pi() / 2.); // rotate vector by 90 deg around z
      TLorentzVector vecPerpMinus(*(jetSel->MomentumVector()));
      vecPerpMinus.RotateZ(-TMath::Pi() / 2.); // rotate vector by -90 deg around z
//      AliAODJet jetTmp = AliAODJet(vecPerp);
      if(fDebug > 5) printf("TaskV0sInJets: Adding perp. cones number %d, %d\n", iNJetPerp, iNJetPerp + 1);
//      printf("TaskV0sInJets: Adding perp. cone number %d: pT %f, phi %f, eta %f, pT %f, phi %f, eta %f\n",iNJetSel,vecPerp.Pt(),vecPerp.Phi(),vecPerp.Eta(),jetTmp.Pt(),jetTmp.Phi(),jetTmp.Eta());
      new((*jetArrayPerp)[iNJetPerp++]) AliAODJet(vecPerpPlus);  // write perp. cone to the array
      new((*jetArrayPerp)[iNJetPerp++]) AliAODJet(vecPerpMinus);  // write perp. cone to the array
      if(fDebug > 5) printf("TaskV0sInJets: Adding jet number %d\n", iNJetSel);
//      printf("TaskV0sInJets: Adding jet number %d: pT %f, phi %f, eta %f\n",iNJetSel,jetSel->Pt(),jetSel->Phi(),jetSel->Eta());
      new((*jetArraySel)[iNJetSel++]) AliAODJet(*((AliAODJet*)jetSel));  // copy selected jet to the array
    }
    if(fDebug > 5) printf("TaskV0sInJets: Added jets: %d\n", iNJetSel);
    iNJetSel = jetArraySel->GetEntriesFast();
    if(fDebug > 2) printf("TaskV0sInJets: Selected jets in array: %d\n", iNJetSel);
    fh1NJetPerEvent[iCentIndex]->Fill(iNJetSel);
    // fill jet spectra
    for(Int_t iJet = 0; iJet < iNJetSel; iJet++)
    {
      jet = (AliAODJet*)jetArraySel->At(iJet); // load a jet in the list
      fh1PtJet[iCentIndex]->Fill(jet->Pt()); // pt spectrum of selected jets
      fh1EtaJet[iCentIndex]->Fill(jet->Eta()); // eta spectrum of selected jets
      fh2EtaPtJet[iCentIndex]->Fill(jet->Eta(), jet->Pt()); // eta-pT spectrum of selected jets
      fh1PhiJet[iCentIndex]->Fill(jet->Phi()); // phi spectrum of selected jets
      Double_t dAreaExcluded = TMath::Pi() * dRadiusExcludeCone * dRadiusExcludeCone; // area of the cone
      dAreaExcluded -= AreaCircSegment(dRadiusExcludeCone, dEtaMax - jet->Eta()); // positive eta overhang
      dAreaExcluded -= AreaCircSegment(dRadiusExcludeCone, dEtaMax + jet->Eta()); // negative eta overhang
      fh1AreaExcluded->Fill(iCentIndex, dAreaExcluded);
    }
    jet = 0;
  }

  if(bJetEventGood)  // there should be some reconstructed jets
  {
    fh1EventCounterCut->Fill(4); // events with jet(s)
    fh1EventCounterCutCent[iCentIndex]->Fill(4); // events with jet(s)
    if(iNJetSel)
    {
      fh1EventCounterCut->Fill(5); // events with selected jets
      fh1EventCounterCutCent[iCentIndex]->Fill(5);
    }
  }
  if(iNJetSel)
    fh1EventCent2Jets->Fill(fdCentrality);
  else
    fh1EventCent2NoJets->Fill(fdCentrality);

  if(iNJetSel)
  {
    jetRnd = GetRandomCone(jetArraySel, dJetEtaWindow, 2 * fdRadiusJet);
    if(jetRnd)
    {
      fh1NRndConeCent->Fill(iCentIndex);
      fh2EtaPhiRndCone[iCentIndex]->Fill(jetRnd->Eta(), jetRnd->Phi());
    }
    jetMed = GetMedianCluster(jetArrayBg, dJetEtaWindow);
    if(jetMed)
    {
      fh1NMedConeCent->Fill(iCentIndex);
      fh2EtaPhiMedCone[iCentIndex]->Fill(jetMed->Eta(), jetMed->Phi());
    }
  }

  // Loading primary vertex info
  AliAODVertex* primVtx = fAODIn->GetPrimaryVertex(); // get the primary vertex
  Double_t dPrimVtxPos[3]; // primary vertex position {x,y,z}
  primVtx->GetXYZ(dPrimVtxPos);
  fh1VtxZ[iCentIndex]->Fill(dPrimVtxPos[2]);
  fh2VtxXY[iCentIndex]->Fill(dPrimVtxPos[0], dPrimVtxPos[1]);

  //===== Start of loop over V0 candidates =====
  if(fDebug > 2) printf("TaskV0sInJets: Start of V0 loop\n");
  for(Int_t iV0 = 0; iV0 < iNV0s; iV0++)
  {
    v0 = fAODIn->GetV0(iV0); // get next candidate from the list in AOD
    if(!v0)
      continue;

    iNV0CandTot++;

    // Initialization of status indicators
    Bool_t bIsCandidateK0s = kTRUE; // candidate for K0s
    Bool_t bIsCandidateLambda = kTRUE; // candidate for Lambda
    Bool_t bIsCandidateALambda = kTRUE; // candidate for Lambda
    Bool_t bIsInPeakK0s = kFALSE; // candidate within the K0s mass peak
    Bool_t bIsInPeakLambda = kFALSE; // candidate within the Lambda mass peak
    Bool_t bIsInConeJet = kFALSE; // candidate within the jet cones
    Bool_t bIsInConePerp = kFALSE; // candidate within the perpendicular cone
    Bool_t bIsInConeRnd = kFALSE; // candidate within the random cone
    Bool_t bIsInConeMed = kFALSE; // candidate within the median-cluster cone
    Bool_t bIsOutsideCones = kFALSE; // candidate outside excluded cones

    // Invariant mass calculation
    dMassV0K0s = v0->MassK0Short();
    dMassV0Lambda = v0->MassLambda();
    dMassV0ALambda = v0->MassAntiLambda();

    Int_t iCutIndex = 0; // indicator of current selection step
    // 0
    // All V0 candidates
    FillCandidates(dMassV0K0s, dMassV0Lambda, dMassV0ALambda, bIsCandidateK0s, bIsCandidateLambda, bIsCandidateALambda, iCutIndex, iCentIndex);
    iCutIndex++;

    // Skip candidates outside the histogram range
    if((dMassV0K0s < fgkdMassK0sMin) || (dMassV0K0s >= fgkdMassK0sMax))
      bIsCandidateK0s = kFALSE;
    if((dMassV0Lambda < fgkdMassLambdaMin) || (dMassV0Lambda >= fgkdMassLambdaMax))
      bIsCandidateLambda = kFALSE;
    if((dMassV0ALambda < fgkdMassLambdaMin) || (dMassV0ALambda >= fgkdMassLambdaMax))
      bIsCandidateALambda = kFALSE;
    if(!bIsCandidateK0s && !bIsCandidateLambda && !bIsCandidateALambda)
      continue;

    Double_t dPtV0 = TMath::Sqrt(v0->Pt2V0()); // transverse momentum of V0
    vecV0Momentum = TVector3(v0->Px(), v0->Py(), v0->Pz()); // set the vector of V0 momentum

    // Sigma of the mass peak window
    Double_t dMassPeakWindowK0s = dNSigmaMassMax * MassPeakSigmaOld(dPtV0, 0);
    Double_t dMassPeakWindowLambda = dNSigmaMassMax * MassPeakSigmaOld(dPtV0, 1);
//      Double_t dMassPeakWindowK0s = dNSigmaMassMax*MassPeakSigma(iCentIndex,dPtV0,0);
//      Double_t dMassPeakWindowLambda = dNSigmaMassMax*MassPeakSigma(iCentIndex,dPtV0,1);

    // Invariant mass peak selection
    if(TMath::Abs(dMassV0K0s - dMassPDGK0s) < dMassPeakWindowK0s)
      bIsInPeakK0s = kTRUE;
    if(TMath::Abs(dMassV0Lambda - dMassPDGLambda) < dMassPeakWindowLambda)
      bIsInPeakLambda = kTRUE;

    // Retrieving all relevant properties of the V0 candidate
    Bool_t bOnFlyStatus = v0->GetOnFlyStatus(); // online (on fly) reconstructed vs offline reconstructed
    const AliAODTrack* trackPos = (AliAODTrack*)v0->GetDaughter(0); // positive daughter track
    const AliAODTrack* trackNeg = (AliAODTrack*)v0->GetDaughter(1); // negative daughter track
    Double_t dPtDaughterPos = trackPos->Pt(); // transverse momentum of a daughter track
    Double_t dPtDaughterNeg = trackNeg->Pt();
    Double_t dNRowsPos = trackPos->GetTPCClusterInfo(2, 1); // crossed TPC pad rows of a daughter track
    Double_t dNRowsNeg = trackNeg->GetTPCClusterInfo(2, 1);
    Double_t dDCAToPrimVtxPos = TMath::Abs(v0->DcaPosToPrimVertex()); // dca of a daughter to the primary vertex
    Double_t dDCAToPrimVtxNeg = TMath::Abs(v0->DcaNegToPrimVertex());
    Double_t dDCADaughters = v0->DcaV0Daughters(); // dca between daughters
    Double_t dCPA = v0->CosPointingAngle(primVtx); // cosine of the pointing angle
    Double_t dSecVtxPos[3]; // V0 vertex position {x,y,z}
//      Double_t dSecVtxPos[3] = {v0->DecayVertexV0X(),v0->DecayVertexV0Y(),v0->DecayVertexV0Z()}; // V0 vertex position
    v0->GetSecondaryVtx(dSecVtxPos);
    Double_t dRadiusDecay = TMath::Sqrt(dSecVtxPos[0] * dSecVtxPos[0] + dSecVtxPos[1] * dSecVtxPos[1]); // distance of the V0 vertex from the z-axis
    Double_t dEtaDaughterNeg = trackNeg->Eta(); // = v0->EtaProng(1), pseudorapidity of a daughter track
    Double_t dEtaDaughterPos = trackPos->Eta(); // = v0->EtaProng(0)
    Double_t dRapK0s = v0->RapK0Short(); // rapidity calculated for K0s assumption
    Double_t dRapLambda = v0->RapLambda(); // rapidity calculated for Lambda assumption
    Double_t dEtaV0 = v0->Eta(); // V0 pseudorapidity
//      Double_t dPhiV0 = v0->Phi(); // V0 pseudorapidity
    Double_t dDecayPath[3];
    for(Int_t iPos = 0; iPos < 3; iPos++)
      dDecayPath[iPos] = dSecVtxPos[iPos] - dPrimVtxPos[iPos]; // vector of the V0 path
    Double_t dDecLen = TMath::Sqrt(dDecayPath[0] * dDecayPath[0] + dDecayPath[1] * dDecayPath[1] + dDecayPath[2] * dDecayPath[2]); // path length L
    Double_t dDecLen2D = TMath::Sqrt(dDecayPath[0] * dDecayPath[0] + dDecayPath[1] * dDecayPath[1]); // transverse path length R
    Double_t dLOverP = dDecLen / v0->P(); // L/p
    Double_t dROverPt = dDecLen2D / dPtV0; // R/pT
    Double_t dMLOverPK0s = dMassPDGK0s * dLOverP; // m*L/p = c*(proper lifetime)
//      Double_t dMLOverPLambda = dMassPDGLambda*dLOverP; // m*L/p
    Double_t dMROverPtK0s = dMassPDGK0s * dROverPt; // m*R/pT
    Double_t dMROverPtLambda = dMassPDGLambda * dROverPt; // m*R/pT
    Double_t dNSigmaPosPion   = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(trackPos, AliPID::kPion)); // difference between measured and expected signal of the dE/dx in the TPC
    Double_t dNSigmaPosProton = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(trackPos, AliPID::kProton));
    Double_t dNSigmaNegPion   = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(trackNeg, AliPID::kPion));
    Double_t dNSigmaNegProton = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(trackNeg, AliPID::kProton));
    Double_t dAlpha = v0->AlphaV0(); // Armenteros-Podolanski alpha
    Double_t dPtArm = v0->PtArmV0(); // Armenteros-Podolanski pT
    AliAODVertex* prodVtxDaughterPos = (AliAODVertex*)(trackPos->GetProdVertex()); // production vertex of the positive daughter track
    Char_t cTypeVtxProdPos = prodVtxDaughterPos->GetType(); // type of the production vertex
    AliAODVertex* prodVtxDaughterNeg = (AliAODVertex*)(trackNeg->GetProdVertex()); // production vertex of the negative daughter track
    Char_t cTypeVtxProdNeg = prodVtxDaughterNeg->GetType(); // type of the production vertex

//    fh2Tau3DVs2D[0]->Fill(dPtV0, dLOverP / dROverPt);

    // QA histograms before cuts
    FillQAHistogramV0(primVtx, v0, 0, bIsCandidateK0s, bIsCandidateLambda, bIsInPeakK0s, bIsInPeakLambda);
    // Cut vs mass histograms before cuts
    /*
    if(bIsCandidateK0s)
    {
      fh2CutTPCRowsK0s[0]->Fill(dMassV0K0s, dNRowsPos);
      fh2CutTPCRowsK0s[0]->Fill(dMassV0K0s, dNRowsNeg);
      fh2CutPtPosK0s[0]->Fill(dMassV0K0s, dPtDaughterPos);
      fh2CutPtNegK0s[0]->Fill(dMassV0K0s, dPtDaughterNeg);
      fh2CutDCAVtx[0]->Fill(dMassV0K0s, dDCAToPrimVtxPos);
      fh2CutDCAVtx[0]->Fill(dMassV0K0s, dDCAToPrimVtxNeg);
      fh2CutDCAV0[0]->Fill(dMassV0K0s, dDCADaughters);
      fh2CutCos[0]->Fill(dMassV0K0s, dCPA);
      fh2CutR[0]->Fill(dMassV0K0s, dRadiusDecay);
      fh2CutEtaK0s[0]->Fill(dMassV0K0s, dEtaDaughterPos);
      fh2CutEtaK0s[0]->Fill(dMassV0K0s, dEtaDaughterNeg);
      fh2CutRapK0s[0]->Fill(dMassV0K0s, dRapK0s);
      fh2CutCTauK0s[0]->Fill(dMassV0K0s, dMROverPtK0s / dCTauK0s);
      fh2CutPIDPosK0s[0]->Fill(dMassV0K0s, dNSigmaPosPion);
      fh2CutPIDNegK0s[0]->Fill(dMassV0K0s, dNSigmaNegPion);
    }
    if(bIsCandidateLambda)
    {
      fh2CutTPCRowsLambda[0]->Fill(dMassV0Lambda, dNRowsPos);
      fh2CutTPCRowsLambda[0]->Fill(dMassV0Lambda, dNRowsNeg);
      fh2CutPtPosLambda[0]->Fill(dMassV0Lambda, dPtDaughterPos);
      fh2CutPtNegLambda[0]->Fill(dMassV0Lambda, dPtDaughterNeg);
      fh2CutEtaLambda[0]->Fill(dMassV0Lambda, dEtaDaughterPos);
      fh2CutEtaLambda[0]->Fill(dMassV0Lambda, dEtaDaughterNeg);
      fh2CutRapLambda[0]->Fill(dMassV0Lambda, dRapLambda);
      fh2CutCTauLambda[0]->Fill(dMassV0Lambda, dMROverPtLambda / dCTauLambda);
      fh2CutPIDPosLambda[0]->Fill(dMassV0Lambda, dNSigmaPosProton);
      fh2CutPIDNegLambda[0]->Fill(dMassV0Lambda, dNSigmaNegPion);
    }
    */

    //===== Start of reconstruction cutting =====

    // 1
    // All V0 candidates
    FillCandidates(dMassV0K0s, dMassV0Lambda, dMassV0ALambda, bIsCandidateK0s, bIsCandidateLambda, bIsCandidateALambda, iCutIndex, iCentIndex);
    iCutIndex++;

    // Start of global cuts
    // 2
    // Reconstruction method
    if(bPrintCuts) printf("Rec: Applying cut: Reconstruction method: on-the-fly? %s\n", (bOnFly ? "yes" : "no"));
    if(bOnFlyStatus != bOnFly)
      continue;
    FillCandidates(dMassV0K0s, dMassV0Lambda, dMassV0ALambda, bIsCandidateK0s, bIsCandidateLambda, bIsCandidateALambda, iCutIndex, iCentIndex);
    iCutIndex++;

    // 3
    // Tracks TPC OK
    if(bPrintCuts) printf("Rec: Applying cut: Correct charge of daughters\n");
    if(!trackNeg || !trackPos)
      continue;
    if(trackNeg->Charge() == trackPos->Charge())  // daughters have different charge?
      continue;
    if(trackNeg->Charge() != -1)  // daughters have expected charge?
      continue;
    if(trackPos->Charge() != 1)  // daughters have expected charge?
      continue;

    if(bPrintCuts) printf("Rec: Applying cut: TPC refit: %d\n", iRefit);
    if(!trackNeg->IsOn(iRefit))  // TPC refit is ON?
      continue;
    if(bPrintCuts) printf("Rec: Applying cut: Type of production vertex of daughter: Not %d\n", AliAODVertex::kKink);
    if(cTypeVtxProdNeg == AliAODVertex::kKink) // kink daughter rejection
      continue;
    // Old cuts Start
    if(bUseOldCuts)
    {
      if(bPrintCuts) printf("Rec: Applying cut: Number of TPC rows: > %f\n", dNCrossedRowsTPCMin);
      if(dNRowsNeg < dNCrossedRowsTPCMin)  // Crossed TPC padrows
        continue;
//      Int_t findable = trackNeg->GetTPCNclsF(); // Findable clusters
//      if (findable <= 0)
//        continue;
//      if (dNRowsNeg/findable < dCrossedRowsOverFindMin)
//        continue;
//      if (dNRowsNeg/findable > dCrossedRowsOverFindMax)
//        continue;
    }
    // Old cuts End

    if(!trackPos->IsOn(iRefit))
      continue;
    if(cTypeVtxProdPos == AliAODVertex::kKink) // kink daughter rejection
      continue;
    // Old cuts Start
    if(bUseOldCuts)
    {
      if(dNRowsPos < dNCrossedRowsTPCMin)
        continue;
//      findable = trackPos->GetTPCNclsF();
//      if (findable <= 0)
//        continue;
//      if (dNRowsPos/findable < dCrossedRowsOverFindMin)
//        continue;
//      if (dNRowsPos/findable > dCrossedRowsOverFindMax)
//        continue;
    }
    // Old cuts End

    FillCandidates(dMassV0K0s, dMassV0Lambda, dMassV0ALambda, bIsCandidateK0s, bIsCandidateLambda, bIsCandidateALambda, iCutIndex, iCentIndex);
    iCutIndex++;

    // 4
    // Daughters: transverse momentum cut
    if(bUseOldCuts)
    {
      if(bPrintCuts) printf("Rec: Applying cut: Daughter pT: > %f\n", dPtDaughterMin);
      if((dPtDaughterNeg < dPtDaughterMin) || (dPtDaughterPos < dPtDaughterMin))
        continue;
      FillCandidates(dMassV0K0s, dMassV0Lambda, dMassV0ALambda, bIsCandidateK0s, bIsCandidateLambda, bIsCandidateALambda, iCutIndex, iCentIndex);
    }
    iCutIndex++;

    // 5
    // Daughters: Impact parameter of daughters to prim vtx
    if(bPrintCuts) printf("Rec: Applying cut: Daughter DCA to prim vtx: > %f\n", dDCAToPrimVtxMin);
    if((dDCAToPrimVtxNeg < dDCAToPrimVtxMin) || (dDCAToPrimVtxPos < dDCAToPrimVtxMin))
      continue;
    FillCandidates(dMassV0K0s, dMassV0Lambda, dMassV0ALambda, bIsCandidateK0s, bIsCandidateLambda, bIsCandidateALambda, iCutIndex, iCentIndex);
    iCutIndex++;

    // 6
    // Daughters: DCA
    if(bPrintCuts) printf("Rec: Applying cut: DCA between daughters: < %f\n", dDCADaughtersMax);
    if(dDCADaughters > dDCADaughtersMax)
      continue;
    FillCandidates(dMassV0K0s, dMassV0Lambda, dMassV0ALambda, bIsCandidateK0s, bIsCandidateLambda, bIsCandidateALambda, iCutIndex, iCentIndex);
    iCutIndex++;

    // 7
    // V0: Cosine of the pointing angle
    if(bPrintCuts) printf("Rec: Applying cut: CPA: > %f\n", dCPAMin);
    if(dCPA < dCPAMin)
      continue;
    FillCandidates(dMassV0K0s, dMassV0Lambda, dMassV0ALambda, bIsCandidateK0s, bIsCandidateLambda, bIsCandidateALambda, iCutIndex, iCentIndex);
    iCutIndex++;

    // 8
    // V0: Fiducial volume
    if(bPrintCuts) printf("Rec: Applying cut: Decay radius: > %f, < %f\n", dRadiusDecayMin, dRadiusDecayMax);
    if((dRadiusDecay < dRadiusDecayMin) || (dRadiusDecay > dRadiusDecayMax))
      continue;
    FillCandidates(dMassV0K0s, dMassV0Lambda, dMassV0ALambda, bIsCandidateK0s, bIsCandidateLambda, bIsCandidateALambda, iCutIndex, iCentIndex);
    iCutIndex++;

    // 9
    // Daughters: pseudorapidity cut
    if(bCutEtaDaughter)
    {
      if(bPrintCuts) printf("Rec: Applying cut: Daughter |eta|: < %f\n", dEtaDaughterMax);
      if((TMath::Abs(dEtaDaughterNeg) > dEtaDaughterMax) || (TMath::Abs(dEtaDaughterPos) > dEtaDaughterMax))
        continue;
      FillCandidates(dMassV0K0s, dMassV0Lambda, dMassV0ALambda, bIsCandidateK0s, bIsCandidateLambda, bIsCandidateALambda, iCutIndex, iCentIndex);
    }
    iCutIndex++;
    // End of global cuts

    // Start of particle-dependent cuts
    // 10
    // V0: rapidity cut & pseudorapidity cut
    if(bCutRapV0)
    {
      if(bPrintCuts) printf("Rec: Applying cut: V0 |y|: < %f\n", dRapMax);
      if(TMath::Abs(dRapK0s) > dRapMax)
        bIsCandidateK0s = kFALSE;
      if(TMath::Abs(dRapLambda) > dRapMax)
      {
        bIsCandidateLambda = kFALSE;
        bIsCandidateALambda = kFALSE;
      }
    }
    if(bCutEtaV0)
    {
      if(bPrintCuts) printf("Rec: Applying cut: V0 |eta|: < %f\n", dEtaMax);
      if(TMath::Abs(dEtaV0) > dEtaMax)
      {
        bIsCandidateK0s = kFALSE;
        bIsCandidateLambda = kFALSE;
        bIsCandidateALambda = kFALSE;
      }
      FillCandidates(dMassV0K0s, dMassV0Lambda, dMassV0ALambda, bIsCandidateK0s, bIsCandidateLambda, bIsCandidateALambda, iCutIndex, iCentIndex);
    }
    iCutIndex++;

    // 11
    // Lifetime cut
    if(bCutTau)
    {
      if(bPrintCuts) printf("Rec: Applying cut: Proper lifetime: < %f\n", dNTauMax);
      if(dMROverPtK0s > dNTauMax * dCTauK0s)
        bIsCandidateK0s = kFALSE;
      if(dMROverPtLambda > dNTauMax * dCTauLambda)
      {
        bIsCandidateLambda = kFALSE;
        bIsCandidateALambda = kFALSE;
      }
      FillCandidates(dMassV0K0s, dMassV0Lambda, dMassV0ALambda, bIsCandidateK0s, bIsCandidateLambda, bIsCandidateALambda, iCutIndex, iCentIndex);
    }
    iCutIndex++;

    // 12
    // Daughter PID
    if(bCutPid)
    {
      if(bUseOldCuts)
      {
        if(bPrintCuts) printf("Rec: Applying cut: Delta dE/dx (both daughters): < %f\n", dNSigmadEdxMax);
        if(dNSigmaPosPion > dNSigmadEdxMax || dNSigmaNegPion > dNSigmadEdxMax)  // pi+, pi-
          bIsCandidateK0s = kFALSE;
        if(dNSigmaPosProton > dNSigmadEdxMax || dNSigmaNegPion > dNSigmadEdxMax)  // p+, pi-
          bIsCandidateLambda = kFALSE;
        if(dNSigmaNegProton > dNSigmadEdxMax || dNSigmaPosPion > dNSigmadEdxMax)  // p-, pi+
          bIsCandidateALambda = kFALSE;
      }
      else
      {
        if(bPrintCuts) printf("Rec: Applying cut: Delta dE/dx (proton below %f GeV/c): < %f\n", dPtProtonPIDMax, dNSigmadEdxMax);
        if((dPtDaughterPos < dPtProtonPIDMax) && (dNSigmaPosProton > dNSigmadEdxMax))    // p+
          bIsCandidateLambda = kFALSE;
        if((dPtDaughterNeg < dPtProtonPIDMax) && (dNSigmaNegProton > dNSigmadEdxMax))    // p-
          bIsCandidateALambda = kFALSE;
      }
      FillCandidates(dMassV0K0s, dMassV0Lambda, dMassV0ALambda, bIsCandidateK0s, bIsCandidateLambda, bIsCandidateALambda, iCutIndex, iCentIndex);
    }
    iCutIndex++;

    Double_t valueCorrel[3] = {dMassV0K0s, dMassV0Lambda, dPtV0};
    if(bIsCandidateK0s && bIsCandidateLambda)
      fh3CCMassCorrelBoth->Fill(valueCorrel); // correlation of mass distribution of candidates selected as both K0s and Lambda
    if(bIsCandidateK0s && !bIsCandidateLambda)
      fh3CCMassCorrelKNotL->Fill(valueCorrel); // correlation of mass distribution of candidates selected as K0s and not Lambda
    if(!bIsCandidateK0s && bIsCandidateLambda)
      fh3CCMassCorrelLNotK->Fill(valueCorrel); // correlation of mass distribution of candidates selected as not K0s and Lambda

    // 13
    // Armenteros-Podolanski cut
    if(bCutArmPod)
    {
      if(bPrintCuts) printf("Rec: Applying cut: Armenteros-Podolanski (K0S): pT > %f * |alpha|\n", 0.2);
      if(dPtArm < TMath::Abs(0.2 * dAlpha))
        bIsCandidateK0s = kFALSE;
//      if(dPtArm < 0.025)
//      {
//        bIsCandidateLambda = kFALSE;
//        bIsCandidateALambda = kFALSE;
//      }
      FillCandidates(dMassV0K0s, dMassV0Lambda, dMassV0ALambda, bIsCandidateK0s, bIsCandidateLambda, bIsCandidateALambda, iCutIndex, iCentIndex);
    }
    iCutIndex++;

    // Cross contamination
    if(bIsInPeakK0s)
    {
      if(bIsCandidateLambda)  // Lambda candidates in K0s peak, excluded from Lambda candidates by CC cut
        fh2CCLambda->Fill(dMassV0Lambda, dPtV0);
    }
    if(bIsInPeakLambda)
    {
      if(bIsCandidateK0s)  // K0s candidates in Lambda peak, excluded from K0s candidates by CC cut
        fh2CCK0s->Fill(dMassV0K0s, dPtV0);
    }
//      if (bCutCross)
//        {
//          if (bIsInPeakK0s)
//            bIsCandidateLambda = kFALSE;
//          if (bIsInPeakLambda)
//            bIsCandidateK0s = kFALSE;
//          FillCandidates(dMassV0K0s, dMassV0Lambda, dMassV0ALambda, bIsCandidateK0s, bIsCandidateLambda, bIsCandidateALambda, iCutIndex, iCentIndex);
//        }
//      iCutIndex++;

    // End of particle-dependent cuts

    //===== End of reconstruction cutting =====

    if(!bIsCandidateK0s && !bIsCandidateLambda && !bIsCandidateALambda)
      continue;

/*
    if(fDebug>5) printf("TaskV0sInJets: Adding selected V0 to branch\n");
    // Add selected candidates to the output tree branch
    if ((bIsCandidateK0s || bIsCandidateLambda || bIsCandidateALambda) && fbTreeOutput)
      {
        objectV0 = new ((*fBranchV0Rec)[iNV0SelV0Rec++]) AliV0Object(v0,primVtx);
//        new ((*fBranchV0Rec)[iNV0SelV0Rec++]) AliAODv0(*((AliAODv0*)v0));
        objectV0->SetIsCandidateK0S(bIsCandidateK0s);
        objectV0->SetIsCandidateLambda(bIsCandidateLambda);
        objectV0->SetIsCandidateALambda(bIsCandidateALambda);
        objectV0->SetNSigmaPosProton(dNSigmaPosProton);
        objectV0->SetNSigmaNegProton(dNSigmaNegProton);
      }
*/

    // Selection of V0s in jet cones, perpendicular cones, random cones, outside cones
    if(bJetEventGood && iNJetSel && (bIsCandidateK0s || bIsCandidateLambda || bIsCandidateALambda))
    {
      // Selection of V0s in jet cones
      if(fDebug > 5) printf("TaskV0sInJets: Searching for V0 %d %d in %d jet cones\n", bIsCandidateK0s, bIsCandidateLambda, iNJetSel);
      for(Int_t iJet = 0; iJet < iNJetSel; iJet++)
      {
        jet = (AliAODJet*)jetArraySel->At(iJet); // load a jet in the list
        vecJetMomentum = TVector3(jet->Px(), jet->Py(), jet->Pz()); // set the vector of jet momentum
        if(fDebug > 5) printf("TaskV0sInJets: Checking if V0 %d %d in jet cone %d\n", bIsCandidateK0s, bIsCandidateLambda, iJet);
        if(IsParticleInCone(v0, jet, fdRadiusJet)) // If good jet in event, find out whether V0 is in that jet
        {
          if(fDebug > 5) printf("TaskV0sInJets: V0 %d %d found in jet cone %d\n", bIsCandidateK0s, bIsCandidateLambda, iJet);
          bIsInConeJet = kTRUE;
          break;
        }
      }
      // Selection of V0s in perp. cones
      if(fDebug > 5) printf("TaskV0sInJets: Searching for V0 %d %d in %d perp. cones\n", bIsCandidateK0s, bIsCandidateLambda, iNJetSel);
      for(Int_t iJet = 0; iJet < iNJetPerp; iJet++)
      {
        jetPerp = (AliAODJet*)jetArrayPerp->At(iJet); // load a jet in the list
        if(fDebug > 5) printf("TaskV0sInJets: Checking if V0 %d %d in perp. cone %d\n", bIsCandidateK0s, bIsCandidateLambda, iJet);
        if(IsParticleInCone(v0, jetPerp, fdRadiusJet)) // V0 in perp. cone
        {
          if(fDebug > 5) printf("TaskV0sInJets: V0 %d %d found in perp. cone %d\n", bIsCandidateK0s, bIsCandidateLambda, iJet);
          bIsInConePerp = kTRUE;
          break;
        }
      }
      // Selection of V0s in random cones
      if(jetRnd)
      {
        if(fDebug > 5) printf("TaskV0sInJets: Searching for V0 %d %d in the rnd. cone\n", bIsCandidateK0s, bIsCandidateLambda);
        if(IsParticleInCone(v0, jetRnd, fdRadiusJet)) // V0 in rnd. cone?
        {
          if(fDebug > 5) printf("TaskV0sInJets: V0 %d %d found in the rnd. cone\n", bIsCandidateK0s, bIsCandidateLambda);
          bIsInConeRnd = kTRUE;
        }
      }
      // Selection of V0s in median-cluster cones
      if(jetMed)
      {
        if(fDebug > 5) printf("TaskV0sInJets: Searching for V0 %d %d in the med. cone\n", bIsCandidateK0s, bIsCandidateLambda);
        if(IsParticleInCone(v0, jetMed, fdRadiusJet)) // V0 in med. cone?
        {
          if(fDebug > 5) printf("TaskV0sInJets: V0 %d %d found in the med. cone\n", bIsCandidateK0s, bIsCandidateLambda);
          bIsInConeMed = kTRUE;
        }
      }
      // Selection of V0s outside jet cones
      if(fDebug > 5) printf("TaskV0sInJets: Searching for V0 %d %d outside jet cones\n", bIsCandidateK0s, bIsCandidateLambda);
      if(!OverlapWithJets(jetArraySel, v0, dRadiusExcludeCone)) // V0 oustide jet cones
      {
        if(fDebug > 5) printf("TaskV0sInJets: V0 %d %d found outside jet cones\n", bIsCandidateK0s, bIsCandidateLambda);
        bIsOutsideCones = kTRUE;
      }
    }

    // QA histograms after cuts
    FillQAHistogramV0(primVtx, v0, 1, bIsCandidateK0s, bIsCandidateLambda, bIsInPeakK0s, bIsInPeakLambda);
    // Cut vs mass histograms after cuts
    /*
    if(bIsCandidateK0s)
    {
      fh2CutTPCRowsK0s[1]->Fill(dMassV0K0s, dNRowsPos);
      fh2CutTPCRowsK0s[1]->Fill(dMassV0K0s, dNRowsNeg);
      fh2CutPtPosK0s[1]->Fill(dMassV0K0s, dPtDaughterPos);
      fh2CutPtNegK0s[1]->Fill(dMassV0K0s, dPtDaughterNeg);
      fh2CutDCAVtx[1]->Fill(dMassV0K0s, dDCAToPrimVtxPos);
      fh2CutDCAVtx[1]->Fill(dMassV0K0s, dDCAToPrimVtxNeg);
      fh2CutDCAV0[1]->Fill(dMassV0K0s, dDCADaughters);
      fh2CutCos[1]->Fill(dMassV0K0s, dCPA);
      fh2CutR[1]->Fill(dMassV0K0s, dRadiusDecay);
      fh2CutEtaK0s[1]->Fill(dMassV0K0s, dEtaDaughterPos);
      fh2CutEtaK0s[1]->Fill(dMassV0K0s, dEtaDaughterNeg);
      fh2CutRapK0s[1]->Fill(dMassV0K0s, dRapK0s);
      fh2CutCTauK0s[1]->Fill(dMassV0K0s, dMROverPtK0s / dCTauK0s);
      fh2CutPIDPosK0s[1]->Fill(dMassV0K0s, dNSigmaPosPion);
      fh2CutPIDNegK0s[1]->Fill(dMassV0K0s, dNSigmaNegPion);
      fh1DeltaZK0s[iCentIndex]->Fill(dDecayPath[2]);
    }
    if(bIsCandidateLambda)
    {
      fh2CutTPCRowsLambda[1]->Fill(dMassV0Lambda, dNRowsPos);
      fh2CutTPCRowsLambda[1]->Fill(dMassV0Lambda, dNRowsNeg);
      fh2CutPtPosLambda[1]->Fill(dMassV0Lambda, dPtDaughterPos);
      fh2CutPtNegLambda[1]->Fill(dMassV0Lambda, dPtDaughterNeg);
      fh2CutEtaLambda[1]->Fill(dMassV0Lambda, dEtaDaughterPos);
      fh2CutEtaLambda[1]->Fill(dMassV0Lambda, dEtaDaughterNeg);
      fh2CutRapLambda[1]->Fill(dMassV0Lambda, dRapLambda);
      fh2CutCTauLambda[1]->Fill(dMassV0Lambda, dMROverPtLambda / dCTauLambda);
      fh2CutPIDPosLambda[1]->Fill(dMassV0Lambda, dNSigmaPosProton);
      fh2CutPIDNegLambda[1]->Fill(dMassV0Lambda, dNSigmaNegPion);
      fh1DeltaZLambda[iCentIndex]->Fill(dDecayPath[2]);
    }
    */

    //===== Start of filling V0 spectra =====

    Double_t dAngle = TMath::Pi(); // angle between V0 momentum and jet momentum
    if(bIsInConeJet)
    {
      dAngle = vecV0Momentum.Angle(vecJetMomentum);
    }

    // iCutIndex = 14
    if(bIsCandidateK0s)
    {
      // 14 K0s candidates after cuts
//          printf("K0S: i = %d, m = %f, pT = %f, eta = %f, phi = %f\n",iV0,dMassV0K0s,dPtV0,dEtaV0,dPhiV0);
      FillCandidates(dMassV0K0s, dMassV0Lambda, dMassV0ALambda, bIsCandidateK0s, kFALSE, kFALSE, iCutIndex, iCentIndex);
      Double_t valueKIncl[3] = {dMassV0K0s, dPtV0, dEtaV0};
      fhnV0InclusiveK0s[iCentIndex]->Fill(valueKIncl);
      fh1V0InvMassK0sCent[iCentIndex]->Fill(dMassV0K0s);

      fh1QACTau2D[1]->Fill(dMROverPtK0s / dCTauK0s);
      fh1QACTau3D[1]->Fill(dMLOverPK0s / dCTauK0s);
//      fh2Tau3DVs2D[1]->Fill(dPtV0, dLOverP / dROverPt);

      if(iNJetSel)
      {
        // 15 K0s in jet events
        FillCandidates(dMassV0K0s, dMassV0Lambda, dMassV0ALambda, bIsCandidateK0s, kFALSE, kFALSE, iCutIndex + 1, iCentIndex);
      }
      if(bIsInConeJet)
      {
        // 16 K0s in jets
        FillCandidates(dMassV0K0s, dMassV0Lambda, dMassV0ALambda, bIsCandidateK0s, kFALSE, kFALSE, iCutIndex + 2, iCentIndex);
        Double_t valueKInJC[4] = {dMassV0K0s, dPtV0, dEtaV0, jet->Pt()};
        fhnV0InJetK0s[iCentIndex]->Fill(valueKInJC);
        fh2V0PtJetAngleK0s[iCentIndex]->Fill(jet->Pt(), dAngle);
      }
      if(bIsOutsideCones)
      {
        Double_t valueKOutJC[3] = {dMassV0K0s, dPtV0, dEtaV0};
        fhnV0OutJetK0s[iCentIndex]->Fill(valueKOutJC);
      }
      if(bIsInConePerp)
      {
        Double_t valueKInPC[4] = {dMassV0K0s, dPtV0, dEtaV0, jetPerp->Pt()};
        fhnV0InPerpK0s[iCentIndex]->Fill(valueKInPC);
      }
      if(bIsInConeRnd)
      {
        Double_t valueKInRnd[3] = {dMassV0K0s, dPtV0, dEtaV0};
        fhnV0InRndK0s[iCentIndex]->Fill(valueKInRnd);
      }
      if(bIsInConeMed)
      {
        Double_t valueKInMed[3] = {dMassV0K0s, dPtV0, dEtaV0};
        fhnV0InMedK0s[iCentIndex]->Fill(valueKInMed);
      }
      if(!iNJetSel)
      {
        Double_t valueKNoJet[3] = {dMassV0K0s, dPtV0, dEtaV0};
        fhnV0NoJetK0s[iCentIndex]->Fill(valueKNoJet);
      }
      iNV0CandK0s++;
    }
    if(bIsCandidateLambda)
    {
      // 14 Lambda candidates after cuts
//          printf("La: i = %d, m = %f, pT = %f, eta = %f, phi = %f\n",iV0,dMassV0Lambda,dPtV0,dEtaV0,dPhiV0);
      FillCandidates(dMassV0K0s, dMassV0Lambda, dMassV0ALambda, kFALSE, bIsCandidateLambda, kFALSE, iCutIndex, iCentIndex);
      Double_t valueLIncl[3] = {dMassV0Lambda, dPtV0, dEtaV0};
      fhnV0InclusiveLambda[iCentIndex]->Fill(valueLIncl);
      fh1V0InvMassLambdaCent[iCentIndex]->Fill(dMassV0Lambda);
      if(iNJetSel)
      {
        // 15 Lambda in jet events
        FillCandidates(dMassV0K0s, dMassV0Lambda, dMassV0ALambda, kFALSE, bIsCandidateLambda, kFALSE, iCutIndex + 1, iCentIndex);
      }
      if(bIsInConeJet)
      {
        // 16 Lambda in jets
        FillCandidates(dMassV0K0s, dMassV0Lambda, dMassV0ALambda, kFALSE, bIsCandidateLambda, kFALSE, iCutIndex + 2, iCentIndex);
        Double_t valueLInJC[4] = {dMassV0Lambda, dPtV0, dEtaV0, jet->Pt()};
        fhnV0InJetLambda[iCentIndex]->Fill(valueLInJC);
        fh2V0PtJetAngleLambda[iCentIndex]->Fill(jet->Pt(), dAngle);
      }
      if(bIsOutsideCones)
      {
        Double_t valueLOutJet[3] = {dMassV0Lambda, dPtV0, dEtaV0};
        fhnV0OutJetLambda[iCentIndex]->Fill(valueLOutJet);
      }
      if(bIsInConePerp)
      {
        Double_t valueLInPC[4] = {dMassV0Lambda, dPtV0, dEtaV0, jetPerp->Pt()};
        fhnV0InPerpLambda[iCentIndex]->Fill(valueLInPC);
      }
      if(bIsInConeRnd)
      {
        Double_t valueLInRnd[3] = {dMassV0Lambda, dPtV0, dEtaV0};
        fhnV0InRndLambda[iCentIndex]->Fill(valueLInRnd);
      }
      if(bIsInConeMed)
      {
        Double_t valueLInMed[3] = {dMassV0Lambda, dPtV0, dEtaV0};
        fhnV0InMedLambda[iCentIndex]->Fill(valueLInMed);
      }
      if(!iNJetSel)
      {
        Double_t valueLNoJet[3] = {dMassV0Lambda, dPtV0, dEtaV0};
        fhnV0NoJetLambda[iCentIndex]->Fill(valueLNoJet);
      }
      iNV0CandLambda++;
    }
    if(bIsCandidateALambda)
    {
      // 14 ALambda candidates after cuts
//          printf("AL: i = %d, m = %f, pT = %f, eta = %f, phi = %f\n",iV0,dMassV0ALambda,dPtV0,dEtaV0,dPhiV0);
      FillCandidates(dMassV0K0s, dMassV0Lambda, dMassV0ALambda, kFALSE, kFALSE, bIsCandidateALambda, iCutIndex, iCentIndex);
      Double_t valueALIncl[3] = {dMassV0ALambda, dPtV0, dEtaV0};
      fhnV0InclusiveALambda[iCentIndex]->Fill(valueALIncl);
      fh1V0InvMassALambdaCent[iCentIndex]->Fill(dMassV0ALambda);
      if(iNJetSel)
      {
        // 15 ALambda in jet events
        FillCandidates(dMassV0K0s, dMassV0Lambda, dMassV0ALambda, kFALSE, kFALSE, bIsCandidateALambda, iCutIndex + 1, iCentIndex);
      }
      if(bIsInConeJet)
      {
        // 16 ALambda in jets
        FillCandidates(dMassV0K0s, dMassV0Lambda, dMassV0ALambda, kFALSE, kFALSE, bIsCandidateALambda, iCutIndex + 2, iCentIndex);
        Double_t valueLInJC[4] = {dMassV0ALambda, dPtV0, dEtaV0, jet->Pt()};
        fhnV0InJetALambda[iCentIndex]->Fill(valueLInJC);
        fh2V0PtJetAngleALambda[iCentIndex]->Fill(jet->Pt(), dAngle);
      }
      if(bIsOutsideCones)
      {
        Double_t valueALOutJet[3] = {dMassV0ALambda, dPtV0, dEtaV0};
        fhnV0OutJetALambda[iCentIndex]->Fill(valueALOutJet);
      }
      if(bIsInConePerp)
      {
        Double_t valueLInPC[4] = {dMassV0ALambda, dPtV0, dEtaV0, jetPerp->Pt()};
        fhnV0InPerpALambda[iCentIndex]->Fill(valueLInPC);
      }
      if(bIsInConeRnd)
      {
        Double_t valueALInRnd[3] = {dMassV0ALambda, dPtV0, dEtaV0};
        fhnV0InRndALambda[iCentIndex]->Fill(valueALInRnd);
      }
      if(bIsInConeMed)
      {
        Double_t valueALInMed[3] = {dMassV0ALambda, dPtV0, dEtaV0};
        fhnV0InMedALambda[iCentIndex]->Fill(valueALInMed);
      }
      if(!iNJetSel)
      {
        Double_t valueALNoJet[3] = {dMassV0ALambda, dPtV0, dEtaV0};
        fhnV0NoJetALambda[iCentIndex]->Fill(valueALNoJet);
      }
      iNV0CandALambda++;
    }
    //===== End of filling V0 spectra =====


    //===== Association of reconstructed V0 candidates with MC particles =====
    if(fbMCAnalysis)
    {
      // Associate selected candidates only
//          if ( !(bIsCandidateK0s && bIsInPeakK0s) && !(bIsCandidateLambda && bIsInPeakLambda) ) // signal candidates
      if(!(bIsCandidateK0s) && !(bIsCandidateLambda)  && !(bIsCandidateALambda))    // chosen candidates with any mass
        continue;

      // Get MC labels of reconstructed daughter tracks
      Int_t iLabelPos = TMath::Abs(trackPos->GetLabel());
      Int_t iLabelNeg = TMath::Abs(trackNeg->GetLabel());

      // Make sure MC daughters are in the array range
      if((iLabelNeg < 0) || (iLabelNeg >= iNTracksMC) || (iLabelPos < 0) || (iLabelPos >= iNTracksMC))
        continue;

      // Get MC particles corresponding to reconstructed daughter tracks
      AliAODMCParticle* particleMCDaughterNeg = (AliAODMCParticle*)arrayMC->At(iLabelNeg);
      AliAODMCParticle* particleMCDaughterPos = (AliAODMCParticle*)arrayMC->At(iLabelPos);
      if(!particleMCDaughterNeg || !particleMCDaughterPos)
        continue;

      // Make sure MC daughter particles are not physical primary
      if((particleMCDaughterNeg->IsPhysicalPrimary()) || (particleMCDaughterPos->IsPhysicalPrimary()))
        continue;

      // Get identities of MC daughter particles
      Int_t iPdgCodeDaughterPos = particleMCDaughterPos->GetPdgCode();
      Int_t iPdgCodeDaughterNeg = particleMCDaughterNeg->GetPdgCode();

      // Get index of the mother particle for each MC daughter particle
      Int_t iIndexMotherPos = particleMCDaughterPos->GetMother();
      Int_t iIndexMotherNeg = particleMCDaughterNeg->GetMother();

      if((iIndexMotherNeg < 0) || (iIndexMotherNeg >= iNTracksMC) || (iIndexMotherPos < 0) || (iIndexMotherPos >= iNTracksMC))
        continue;

      // Check whether MC daughter particles have the same mother
      if(iIndexMotherNeg != iIndexMotherPos)
        continue;

      // Get the MC mother particle of both MC daughter particles
      AliAODMCParticle* particleMCMother = (AliAODMCParticle*)arrayMC->At(iIndexMotherPos);
      if(!particleMCMother)
        continue;

      // Get identity of the MC mother particle
      Int_t iPdgCodeMother = particleMCMother->GetPdgCode();

      // Skip not interesting particles
      if((iPdgCodeMother != iPdgCodeK0s) && (TMath::Abs(iPdgCodeMother) != iPdgCodeLambda))
        continue;

      // Check identity of the MC mother particle and the decay channel
      // Is MC mother particle K0S?
      Bool_t bV0MCIsK0s = ((iPdgCodeMother == iPdgCodeK0s) && (iPdgCodeDaughterPos == +iPdgCodePion) && (iPdgCodeDaughterNeg == -iPdgCodePion));
      // Is MC mother particle Lambda?
      Bool_t bV0MCIsLambda = ((iPdgCodeMother == +iPdgCodeLambda) && (iPdgCodeDaughterPos == +iPdgCodeProton) && (iPdgCodeDaughterNeg == -iPdgCodePion));
      // Is MC mother particle anti-Lambda?
      Bool_t bV0MCIsALambda = ((iPdgCodeMother == -iPdgCodeLambda) && (iPdgCodeDaughterPos == +iPdgCodePion) && (iPdgCodeDaughterNeg == -iPdgCodeProton));

      Double_t dPtV0Gen = particleMCMother->Pt();
//          Double_t dRapV0MC = particleMCMother->Y();
      Double_t dEtaV0Gen = particleMCMother->Eta();
//          Double_t dPhiV0Gen = particleMCMother->Phi();

      // Is MC mother particle physical primary? Attention!! Definition of IsPhysicalPrimary may change!!
//          Bool_t bV0MCIsPrimary = particleMCMother->IsPhysicalPrimary();
      // Get the MC mother particle of the MC mother particle
      Int_t iIndexMotherOfMother = particleMCMother->GetMother();
      AliAODMCParticle* particleMCMotherOfMother = 0;
      if(iIndexMotherOfMother >= 0)
        particleMCMotherOfMother = (AliAODMCParticle*)arrayMC->At(iIndexMotherOfMother);
      // Get identity of the MC mother particle of the MC mother particle if it exists
      Int_t iPdgCodeMotherOfMother = 0;
      if(particleMCMotherOfMother)
        iPdgCodeMotherOfMother = particleMCMotherOfMother->GetPdgCode();
      // Check if the MC mother particle of the MC mother particle is a physical primary Sigma (3212 - Sigma0, 3224 - Sigma*+, 3214 - Sigma*0, 3114 - Sigma*-)
//          Bool_t bV0MCComesFromSigma = kFALSE; // Is MC mother particle daughter of a Sigma?
//          if ( (particleMCMotherOfMother && particleMCMotherOfMother->IsPhysicalPrimary()) && ( (TMath::Abs(iPdgCodeMotherOfMother)==3212) || (TMath::Abs(iPdgCodeMotherOfMother)==3224) || (TMath::Abs(iPdgCodeMotherOfMother)==3214) || (TMath::Abs(iPdgCodeMotherOfMother)==3114) ) )
//            bV0MCComesFromSigma = kTRUE;
      // Should MC mother particle be considered as primary when it is Lambda?
//          Bool_t bV0MCIsPrimaryLambda = (bV0MCIsPrimary || bV0MCComesFromSigma);
      // Check if the MC mother particle of the MC mother particle is a Xi (3322 - Xi0, 3312 - Xi-)
      Bool_t bV0MCComesFromXi = ((particleMCMotherOfMother) && ((iPdgCodeMotherOfMother == 3322) || (iPdgCodeMotherOfMother == 3312))); // Is MC mother particle daughter of a Xi?
      Bool_t bV0MCComesFromAXi = ((particleMCMotherOfMother) && ((iPdgCodeMotherOfMother == -3322) || (iPdgCodeMotherOfMother == -3312))); // Is MC mother particle daughter of a anti-Xi?

      // Get the distance between production point of the MC mother particle and the primary vertex
      Double_t dx = dPrimVtxMCX - particleMCMother->Xv();
      Double_t dy = dPrimVtxMCY - particleMCMother->Yv();
      Double_t dz = dPrimVtxMCZ - particleMCMother->Zv();
      Double_t dDistPrimary = TMath::Sqrt(dx * dx + dy * dy + dz * dz);
      Bool_t bV0MCIsPrimaryDist = (dDistPrimary < dDistPrimaryMax); // Is close enough to be considered primary-like?

      // phi, eta resolution for V0-reconstruction
//          Double_t dResolutionV0Eta = particleMCMother->Eta()-v0->Eta();
//          Double_t dResolutionV0Phi = particleMCMother->Phi()-v0->Phi();

/*
      if (fbTreeOutput)
      {
      objectV0->SetPtTrue(dPtV0Gen);
      objectV0->SetEtaTrue(dEtaV0Gen);
      objectV0->SetPhiTrue(dPhiV0Gen);
      objectV0->SetPDGCode(iPdgCodeMother);
      objectV0->SetPDGCodeMother(iPdgCodeMotherOfMother);
      }
*/

      // K0s
//          if (bIsCandidateK0s && bIsInPeakK0s) // selected candidates in peak
      if(bIsCandidateK0s)  // selected candidates with any mass
      {
//              if (bV0MCIsK0s && bV0MCIsPrimary) // well reconstructed candidates
        if(bV0MCIsK0s && bV0MCIsPrimaryDist)  // well reconstructed candidates
        {
//                  if (fbTreeOutput)
//                    objectV0->SetOrigin(1);
          fh2V0K0sPtMassMCRec[iCentIndex]->Fill(dPtV0Gen, dMassV0K0s);
          Double_t valueEtaK[3] = {dMassV0K0s, dPtV0Gen, dEtaV0Gen};
          fh3V0K0sEtaPtMassMCRec[iCentIndex]->Fill(valueEtaK);

          Double_t valueEtaDKNeg[6] = {0, particleMCDaughterNeg->Eta(), particleMCDaughterNeg->Pt(), dEtaV0Gen, dPtV0Gen, 0};
          fhnV0K0sInclDaughterEtaPtPtMCRec[iCentIndex]->Fill(valueEtaDKNeg);
          Double_t valueEtaDKPos[6] = {1, particleMCDaughterPos->Eta(), particleMCDaughterPos->Pt(), dEtaV0Gen, dPtV0Gen, 0};
          fhnV0K0sInclDaughterEtaPtPtMCRec[iCentIndex]->Fill(valueEtaDKPos);

          fh2V0K0sMCResolMPt[iCentIndex]->Fill(dMassV0K0s - dMassPDGK0s, dPtV0);
          fh2V0K0sMCPtGenPtRec[iCentIndex]->Fill(dPtV0Gen, dPtV0);
          if(bIsInConeJet)  // true V0 associated to a candidate in jet
          {
            Double_t valueKInJCMC[4] = {dMassV0K0s, dPtV0Gen, dEtaV0Gen, jet->Pt()};
            fh3V0K0sInJetPtMassMCRec[iCentIndex]->Fill(valueKInJCMC);
            Double_t valueEtaKIn[5] = {dMassV0K0s, dPtV0Gen, dEtaV0Gen, jet->Pt(), dEtaV0Gen - jet->Eta()};
            fh4V0K0sInJetEtaPtMassMCRec[iCentIndex]->Fill(valueEtaKIn);

            Double_t valueEtaDKJCNeg[6] = {0, particleMCDaughterNeg->Eta(), particleMCDaughterNeg->Pt(), dEtaV0Gen, dPtV0Gen, jet->Pt()};
            fhnV0K0sInJetsDaughterEtaPtPtMCRec[iCentIndex]->Fill(valueEtaDKJCNeg);
            Double_t valueEtaDKJCPos[6] = {1, particleMCDaughterPos->Eta(), particleMCDaughterPos->Pt(), dEtaV0Gen, dPtV0Gen, jet->Pt()};
            fhnV0K0sInJetsDaughterEtaPtPtMCRec[iCentIndex]->Fill(valueEtaDKJCPos);
          }
        }
        if(bV0MCIsK0s && !bV0MCIsPrimaryDist)  // not primary K0s
        {
//                  if (fbTreeOutput)
//                    objectV0->SetOrigin(-1);
          fh1V0K0sPtMCRecFalse[iCentIndex]->Fill(dPtV0Gen);
        }
      }
      // Lambda
//          if (bIsCandidateLambda && bIsInPeakLambda) // selected candidates in peak
      if(bIsCandidateLambda)  // selected candidates with any mass
      {
//              if (bV0MCIsLambda && bV0MCIsPrimaryLambda) // well reconstructed candidates
        if(bV0MCIsLambda && bV0MCIsPrimaryDist)  // well reconstructed candidates
        {
//                  if (fbTreeOutput)
//                    objectV0->SetOrigin(1);
          fh2V0LambdaPtMassMCRec[iCentIndex]->Fill(dPtV0Gen, dMassV0Lambda);
          Double_t valueEtaL[3] = {dMassV0Lambda, dPtV0Gen, dEtaV0Gen};
          fh3V0LambdaEtaPtMassMCRec[iCentIndex]->Fill(valueEtaL);

          Double_t valueEtaDLNeg[6] = {0, particleMCDaughterNeg->Eta(), particleMCDaughterNeg->Pt(), dEtaV0Gen, dPtV0Gen, 0};
          fhnV0LambdaInclDaughterEtaPtPtMCRec[iCentIndex]->Fill(valueEtaDLNeg);
          Double_t valueEtaDLPos[6] = {1, particleMCDaughterPos->Eta(), particleMCDaughterPos->Pt(), dEtaV0Gen, dPtV0Gen, 0};
          fhnV0LambdaInclDaughterEtaPtPtMCRec[iCentIndex]->Fill(valueEtaDLPos);

          fh2V0LambdaMCResolMPt[iCentIndex]->Fill(dMassV0Lambda - dMassPDGLambda, dPtV0);
          fh2V0LambdaMCPtGenPtRec[iCentIndex]->Fill(dPtV0Gen, dPtV0);
          if(bIsInConeJet)  // true V0 associated to a reconstructed candidate in jet
          {
            Double_t valueLInJCMC[4] = {dMassV0Lambda, dPtV0Gen, dEtaV0Gen, jet->Pt()};
            fh3V0LambdaInJetPtMassMCRec[iCentIndex]->Fill(valueLInJCMC);
            Double_t valueEtaLIn[5] = {dMassV0Lambda, dPtV0Gen, dEtaV0Gen, jet->Pt(), dEtaV0Gen - jet->Eta()};
            fh4V0LambdaInJetEtaPtMassMCRec[iCentIndex]->Fill(valueEtaLIn);

            Double_t valueEtaDLJCNeg[6] = {0, particleMCDaughterNeg->Eta(), particleMCDaughterNeg->Pt(), dEtaV0Gen, dPtV0Gen, jet->Pt()};
            fhnV0LambdaInJetsDaughterEtaPtPtMCRec[iCentIndex]->Fill(valueEtaDLJCNeg);
            Double_t valueEtaDLJCPos[6] = {1, particleMCDaughterPos->Eta(), particleMCDaughterPos->Pt(), dEtaV0Gen, dPtV0Gen, jet->Pt()};
            fhnV0LambdaInJetsDaughterEtaPtPtMCRec[iCentIndex]->Fill(valueEtaDLJCPos);
          }
        }
        // Fill the feed-down histograms
        if(bV0MCIsLambda && bV0MCComesFromXi)
        {
//                  if (fbTreeOutput)
//                    objectV0->SetOrigin(2);
          Double_t valueFDLIncl[3] = {dPtV0Gen, particleMCMotherOfMother->Pt(), 0.};
          fhnV0LambdaInclMCFD[iCentIndex]->Fill(valueFDLIncl);
          if(bIsInConeRnd)
          {
            fhnV0LambdaBulkMCFD[iCentIndex]->Fill(valueFDLIncl);
          }
          if(bIsInConeJet)
          {
            Double_t valueFDLInJets[3] = {dPtV0Gen, particleMCMotherOfMother->Pt(), jet->Pt()};
            fhnV0LambdaInJetsMCFD[iCentIndex]->Fill(valueFDLInJets);
          }
        }
        if(bV0MCIsLambda && !bV0MCIsPrimaryDist && !bV0MCComesFromXi)  // not primary Lambda
        {
//                  if (fbTreeOutput)
//                    objectV0->SetOrigin(-1);
          fh1V0LambdaPtMCRecFalse[iCentIndex]->Fill(dPtV0Gen);
        }
      }
      // anti-Lambda
//          if (bIsCandidateALambda && bIsInPeakALambda) // selected candidates in peak
      if(bIsCandidateALambda)  // selected candidates with any mass
      {
//              if (bV0MCIsALambda && bV0MCIsPrimaryALambda) // well reconstructed candidates
        if(bV0MCIsALambda && bV0MCIsPrimaryDist)  // well reconstructed candidates
        {
//                  if (fbTreeOutput)
//                    objectV0->SetOrigin(1);
          fh2V0ALambdaPtMassMCRec[iCentIndex]->Fill(dPtV0Gen, dMassV0ALambda);
          Double_t valueEtaAL[3] = {dMassV0ALambda, dPtV0Gen, dEtaV0Gen};
          fh3V0ALambdaEtaPtMassMCRec[iCentIndex]->Fill(valueEtaAL);

          Double_t valueEtaDALNeg[6] = {0, particleMCDaughterNeg->Eta(), particleMCDaughterNeg->Pt(), dEtaV0Gen, dPtV0Gen, 0};
          fhnV0ALambdaInclDaughterEtaPtPtMCRec[iCentIndex]->Fill(valueEtaDALNeg);
          Double_t valueEtaDALPos[6] = {1, particleMCDaughterPos->Eta(), particleMCDaughterPos->Pt(), dEtaV0Gen, dPtV0Gen, 0};
          fhnV0ALambdaInclDaughterEtaPtPtMCRec[iCentIndex]->Fill(valueEtaDALPos);

          fh2V0ALambdaMCResolMPt[iCentIndex]->Fill(dMassV0ALambda - dMassPDGLambda, dPtV0);
          fh2V0ALambdaMCPtGenPtRec[iCentIndex]->Fill(dPtV0Gen, dPtV0);
          if(bIsInConeJet)  // true V0 associated to a reconstructed candidate in jet
          {
            Double_t valueALInJCMC[4] = {dMassV0ALambda, dPtV0Gen, dEtaV0Gen, jet->Pt()};
            fh3V0ALambdaInJetPtMassMCRec[iCentIndex]->Fill(valueALInJCMC);
            Double_t valueEtaALIn[5] = {dMassV0ALambda, dPtV0Gen, dEtaV0Gen, jet->Pt(), dEtaV0Gen - jet->Eta()};
            fh4V0ALambdaInJetEtaPtMassMCRec[iCentIndex]->Fill(valueEtaALIn);

            Double_t valueEtaDALJCNeg[6] = {0, particleMCDaughterNeg->Eta(), particleMCDaughterNeg->Pt(), dEtaV0Gen, dPtV0Gen, jet->Pt()};
            fhnV0ALambdaInJetsDaughterEtaPtPtMCRec[iCentIndex]->Fill(valueEtaDALJCNeg);
            Double_t valueEtaDALJCPos[6] = {1, particleMCDaughterPos->Eta(), particleMCDaughterPos->Pt(), dEtaV0Gen, dPtV0Gen, jet->Pt()};
            fhnV0ALambdaInJetsDaughterEtaPtPtMCRec[iCentIndex]->Fill(valueEtaDALJCPos);
          }
        }
        // Fill the feed-down histograms
        if(bV0MCIsALambda && bV0MCComesFromAXi)
        {
//                  if (fbTreeOutput)
//                    objectV0->SetOrigin(2);
          Double_t valueFDALIncl[3] = {dPtV0Gen, particleMCMotherOfMother->Pt(), 0.};
          fhnV0ALambdaInclMCFD[iCentIndex]->Fill(valueFDALIncl);
          if(bIsInConeRnd)
          {
            fhnV0ALambdaBulkMCFD[iCentIndex]->Fill(valueFDALIncl);
          }
          if(bIsInConeJet)
          {
            Double_t valueFDALInJets[3] = {dPtV0Gen, particleMCMotherOfMother->Pt(), jet->Pt()};
            fhnV0ALambdaInJetsMCFD[iCentIndex]->Fill(valueFDALInJets);
          }
        }
        if(bV0MCIsALambda && !bV0MCIsPrimaryDist && !bV0MCComesFromAXi)  // not primary anti-Lambda
        {
//                  if (fbTreeOutput)
//                    objectV0->SetOrigin(-1);
          fh1V0ALambdaPtMCRecFalse[iCentIndex]->Fill(dPtV0Gen);
        }
      }
    }
    //===== End Association of reconstructed V0 candidates with MC particles =====
  }
  //===== End of V0 loop =====

  fh1V0CandPerEvent->Fill(iNV0CandTot);
  fh1V0CandPerEventCentK0s[iCentIndex]->Fill(iNV0CandK0s);
  fh1V0CandPerEventCentLambda[iCentIndex]->Fill(iNV0CandLambda);
  fh1V0CandPerEventCentALambda[iCentIndex]->Fill(iNV0CandALambda);

  if(fDebug > 2) printf("TaskV0sInJets: End of V0 loop\n");

  // Spectra of generated particles
  if(fbMCAnalysis)
  {
    for(Int_t iPartMC = 0; iPartMC < iNTracksMC; iPartMC++)
    {
      // Get MC particle
      AliAODMCParticle* particleMC = (AliAODMCParticle*)arrayMC->At(iPartMC);
      if(!particleMC)
        continue;

      // Get identity of MC particle
      Int_t iPdgCodeParticleMC = particleMC->GetPdgCode();
      // Fill Xi spectrum (3322 - Xi0, 3312 - Xi-)
//          if ( (iPdgCodeParticleMC==3322) || (iPdgCodeParticleMC==3312) )
      if((iPdgCodeParticleMC == 3312) && (TMath::Abs(particleMC->Y()) < 0.5))
      {
//              if (fbTreeOutput)
//                new ((*fBranchV0Gen)[iNV0SelV0Gen++]) AliAODMCParticle(*((AliAODMCParticle*)particleMC));
        fh1V0XiPtMCGen[iCentIndex]->Fill(particleMC->Pt());
      }
      if((iPdgCodeParticleMC == -3312) && (TMath::Abs(particleMC->Y()) < 0.5))
      {
//              if (fbTreeOutput)
//                new ((*fBranchV0Gen)[iNV0SelV0Gen++]) AliAODMCParticle(*((AliAODMCParticle*)particleMC));
        fh1V0AXiPtMCGen[iCentIndex]->Fill(particleMC->Pt());
      }
      // Skip not interesting particles
      if((iPdgCodeParticleMC != iPdgCodeK0s) && (TMath::Abs(iPdgCodeParticleMC) != iPdgCodeLambda))
        continue;

      // Check identity of the MC V0 particle
      // Is MC V0 particle K0S?
      Bool_t bV0MCIsK0s = (iPdgCodeParticleMC == iPdgCodeK0s);
      // Is MC V0 particle Lambda?
      Bool_t bV0MCIsLambda = (iPdgCodeParticleMC == +iPdgCodeLambda);
      // Is MC V0 particle anti-Lambda?
      Bool_t bV0MCIsALambda = (iPdgCodeParticleMC == -iPdgCodeLambda);

      Double_t dPtV0Gen = particleMC->Pt();
      Double_t dRapV0Gen = particleMC->Y();
      Double_t dEtaV0Gen = particleMC->Eta();

      // V0 rapidity cut
      if(bCutRapV0)
      {
        if(bPrintCuts) printf("Gen: Applying cut: V0 |y|: < %f\n", dRapMax);
        if((TMath::Abs(dRapV0Gen) > dRapMax))
          continue;
      }
      // V0 pseudorapidity cut
      if(bCutEtaV0)
      {
        if(bPrintCuts) printf("Gen: Applying cut: V0 |eta|: < %f\n", dEtaMax);
        if((TMath::Abs(dEtaV0Gen) > dEtaMax))
          continue;
      }
/*
      // Is MC V0 particle physical primary? Attention!! Definition of IsPhysicalPrimary may change!!
      Bool_t bV0MCIsPrimary = particleMC->IsPhysicalPrimary();

      // Get the MC mother particle of the MC V0 particle
      Int_t iIndexMotherOfMother = particleMC->GetMother();
      AliAODMCParticle* particleMCMotherOfMother = 0;
      if (iIndexMotherOfMother >= 0)
        particleMCMotherOfMother = (AliAODMCParticle*)arrayMC->At(iIndexMotherOfMother);
      // Get identity of the MC mother particle of the MC V0 particle if it exists
      Int_t iPdgCodeMotherOfMother = 0;
      if (particleMCMotherOfMother)
        iPdgCodeMotherOfMother = particleMCMotherOfMother->GetPdgCode();
      // Check if the MC mother particle is a physical primary Sigma
      Bool_t bV0MCComesFromSigma = kFALSE;
      if ((particleMCMotherOfMother && particleMCMotherOfMother->IsPhysicalPrimary()) && (TMath::Abs(iPdgCodeMotherOfMother)==3212) || (TMath::Abs(iPdgCodeMotherOfMother)==3224) || (TMath::Abs(iPdgCodeMotherOfMother)==3214) || (TMath::Abs(iPdgCodeMotherOfMother)==3114) )
        bV0MCComesFromSigma = kTRUE;
      // Should the MC V0 particle be considered as primary when it is Lambda?
      Bool_t bV0MCIsPrimaryLambda = (bV0MCIsPrimary || bV0MCComesFromSigma);
*/
      // Reject non primary particles
//          if (!bV0MCIsPrimaryLambda)
//            continue;

      // Get the distance between the production point of the MC V0 particle and the primary vertex
      Double_t dx = dPrimVtxMCX - particleMC->Xv();
      Double_t dy = dPrimVtxMCY - particleMC->Yv();
      Double_t dz = dPrimVtxMCZ - particleMC->Zv();
      Double_t dDistPrimary = TMath::Sqrt(dx * dx + dy * dy + dz * dz);
      Bool_t bV0MCIsPrimaryDist = (dDistPrimary < dDistPrimaryMax); // Is close enough to be considered primary-like?

      // Check whether the MC V0 particle is in a MC jet
      AliAODJet* jetMC = 0;
      Bool_t bIsMCV0InJet = kFALSE;
      if(iNJetSel)
      {
        if(fDebug > 5) printf("TaskV0sInJets: Searching for gen V0 in %d MC jets\n", iNJetSel);
        for(Int_t iJet = 0; iJet < iNJetSel; iJet++)
        {
          jetMC = (AliAODJet*)jetArraySel->At(iJet); // load a jet in the list
          if(fDebug > 5) printf("TaskV0sInJets: Checking if gen V0 in MC jet %d\n", iJet);
          if(IsParticleInCone(particleMC, jetMC, fdRadiusJet)) // If good jet in event, find out whether V0 is in that jet
          {
            if(fDebug > 5) printf("TaskV0sInJets: gen V0 found in MC jet %d\n", iJet);
            bIsMCV0InJet = kTRUE;
            break;
          }
        }
      }

      // Select only primary-like MC V0 particles
      // K0s
//          if (bV0MCIsK0s && bV0MCIsPrimary) // well reconstructed candidates
      if(bV0MCIsK0s && bV0MCIsPrimaryDist)  // well reconstructed candidates
      {
//              if (fbTreeOutput)
//                new ((*fBranchV0Gen)[iNV0SelV0Gen++]) AliAODMCParticle(*((AliAODMCParticle*)particleMC));
        fh1V0K0sPtMCGen[iCentIndex]->Fill(dPtV0Gen);
        fh2V0K0sEtaPtMCGen[iCentIndex]->Fill(dPtV0Gen, dEtaV0Gen);
        if(bIsMCV0InJet)
        {
          fh2V0K0sInJetPtMCGen[iCentIndex]->Fill(dPtV0Gen, jetMC->Pt());
          Double_t valueEtaKInGen[4] = {dPtV0Gen, dEtaV0Gen, jetMC->Pt(), dEtaV0Gen - jetMC->Eta()};
          fh3V0K0sInJetEtaPtMCGen[iCentIndex]->Fill(valueEtaKInGen);
        }
      }
      // Lambda
//          if (bV0MCIsLambda && bV0MCIsPrimaryLambda) // well reconstructed candidates
      if(bV0MCIsLambda && bV0MCIsPrimaryDist)  // well reconstructed candidates
      {
//              if (fbTreeOutput)
//                new ((*fBranchV0Gen)[iNV0SelV0Gen++]) AliAODMCParticle(*((AliAODMCParticle*)particleMC));
        fh1V0LambdaPtMCGen[iCentIndex]->Fill(dPtV0Gen);
        fh2V0LambdaEtaPtMCGen[iCentIndex]->Fill(dPtV0Gen, dEtaV0Gen);
        if(bIsMCV0InJet)
        {
          fh2V0LambdaInJetPtMCGen[iCentIndex]->Fill(dPtV0Gen, jetMC->Pt());
          Double_t valueEtaLInGen[4] = {dPtV0Gen, dEtaV0Gen, jetMC->Pt(), dEtaV0Gen - jetMC->Eta()};
          fh3V0LambdaInJetEtaPtMCGen[iCentIndex]->Fill(valueEtaLInGen);
        }
      }
      // anti-Lambda
//          if (bV0MCIsALambda && bV0MCIsPrimaryALambda) // well reconstructed candidates
      if(bV0MCIsALambda && bV0MCIsPrimaryDist)  // well reconstructed candidates
      {
//              if (fbTreeOutput)
//                new ((*fBranchV0Gen)[iNV0SelV0Gen++]) AliAODMCParticle(*((AliAODMCParticle*)particleMC));
        fh1V0ALambdaPtMCGen[iCentIndex]->Fill(dPtV0Gen);
        fh2V0ALambdaEtaPtMCGen[iCentIndex]->Fill(dPtV0Gen, dEtaV0Gen);
        if(bIsMCV0InJet)
        {
          fh2V0ALambdaInJetPtMCGen[iCentIndex]->Fill(dPtV0Gen, jetMC->Pt());
          Double_t valueEtaALInGen[4] = {dPtV0Gen, dEtaV0Gen, jetMC->Pt(), dEtaV0Gen - jetMC->Eta()};
          fh3V0ALambdaInJetEtaPtMCGen[iCentIndex]->Fill(valueEtaALInGen);
        }
      }
    }
  }

//  if (fbTreeOutput)
//    ftreeOut->Fill();

  jetArraySel->Delete();
  delete jetArraySel;
  jetArrayPerp->Delete();
  delete jetArrayPerp;
  if(jetRnd)
    delete jetRnd;
  jetRnd = 0;

  PostData(1, fOutputListStd);
  PostData(2, fOutputListQA);
  PostData(3, fOutputListCuts);
  PostData(4, fOutputListMC);
//  if (fbTreeOutput)
//    PostData(5,ftreeOut);
//  if(fDebug>5) printf("TaskV0sInJets: UserExec: End\n");
}

void AliAnalysisTaskV0sInJets::FillQAHistogramV0(AliAODVertex* vtx, const AliAODv0* vZero, Int_t iIndexHisto, Bool_t IsCandK0s, Bool_t IsCandLambda, Bool_t IsInPeakK0s, Bool_t IsInPeakLambda)
{
  if(!IsCandK0s && !IsCandLambda)
    return;

//  Double_t dMassK0s = vZero->MassK0Short();
//  Double_t dMassLambda = vZero->MassLambda();

  fh1QAV0Status[iIndexHisto]->Fill(vZero->GetOnFlyStatus());

  AliAODTrack* trackNeg = (AliAODTrack*)vZero->GetDaughter(1); // negative track
  AliAODTrack* trackPos = (AliAODTrack*)vZero->GetDaughter(0); // positive track

  Short_t fTotalCharge = 0;
  for(Int_t i = 0; i < 2; i++)
  {
    AliAODTrack* track = (AliAODTrack*)vZero->GetDaughter(i); // track
    // Tracks TPC OK
    fh1QAV0TPCRefit[iIndexHisto]->Fill(track->IsOn(AliAODTrack::kTPCrefit));
    Double_t nCrossedRowsTPC = track->GetTPCClusterInfo(2, 1);
    fh1QAV0TPCRows[iIndexHisto]->Fill(nCrossedRowsTPC);
    Int_t findable = track->GetTPCNclsF();
    fh1QAV0TPCFindable[iIndexHisto]->Fill(findable);
    if(findable != 0)
    {
      fh1QAV0TPCRowsFind[iIndexHisto]->Fill(nCrossedRowsTPC / findable);
    }
    // Daughters: pseudo-rapidity cut
    fh1QAV0Eta[iIndexHisto]->Fill(track->Eta());
    if((nCrossedRowsTPC > (160. / (250. - 85.) * (255.*TMath::Abs(tan(track->Theta())) - 85.)) + 20.) && (track->Eta() < 0) && (track->Pt() > 0.15))
//      if (IsCandK0s)
    {
      fh2QAV0EtaRows[iIndexHisto]->Fill(track->Eta(), nCrossedRowsTPC);
      fh2QAV0PtRows[iIndexHisto]->Fill(track->Pt(), nCrossedRowsTPC);
      fh2QAV0PhiRows[iIndexHisto]->Fill(track->Phi(), nCrossedRowsTPC);
      fh2QAV0NClRows[iIndexHisto]->Fill(findable, nCrossedRowsTPC);
      fh2QAV0EtaNCl[iIndexHisto]->Fill(track->Eta(), findable);
    }

    // Daughters: transverse momentum cut
    fh1QAV0Pt[iIndexHisto]->Fill(track->Pt());
    fTotalCharge += track->Charge();
  }
  fh1QAV0Charge[iIndexHisto]->Fill(fTotalCharge);

  // Daughters: Impact parameter of daughters to prim vtx
  fh1QAV0DCAVtx[iIndexHisto]->Fill(TMath::Abs(vZero->DcaNegToPrimVertex()));
  fh1QAV0DCAVtx[iIndexHisto]->Fill(TMath::Abs(vZero->DcaPosToPrimVertex()));
//  fh2CutDCAVtx[iIndexHisto]->Fill(dMassK0s,TMath::Abs(vZero->DcaNegToPrimVertex()));

  // Daughters: DCA
  fh1QAV0DCAV0[iIndexHisto]->Fill(vZero->DcaV0Daughters());
//  fh2CutDCAV0[iIndexHisto]->Fill(dMassK0s,vZero->DcaV0Daughters());

  // V0: Cosine of the pointing angle
  fh1QAV0Cos[iIndexHisto]->Fill(vZero->CosPointingAngle(vtx));
//  fh2CutCos[iIndexHisto]->Fill(dMassK0s,vZero->CosPointingAngle(vtx));

  // V0: Fiducial volume
  Double_t xyz[3];
  vZero->GetSecondaryVtx(xyz);
  Double_t r2 = xyz[0] * xyz[0] + xyz[1] * xyz[1];
  fh1QAV0R[iIndexHisto]->Fill(TMath::Sqrt(r2));

  Double_t dAlpha = vZero->AlphaV0();
  Double_t dPtArm = vZero->PtArmV0();

  if(IsCandK0s)
  {
    if(IsInPeakK0s)
    {
//          fh2QAV0EtaPtK0sPeak[iIndexHisto]->Fill(trackNeg->Eta(),vZero->Pt());
//          fh2QAV0EtaPtK0sPeak[iIndexHisto]->Fill(trackPos->Eta(),vZero->Pt());
      fh2QAV0EtaPtK0sPeak[iIndexHisto]->Fill(vZero->Eta(), vZero->Pt());
      fh2QAV0PtPtK0sPeak[iIndexHisto]->Fill(trackNeg->Pt(), trackPos->Pt());
      fh2ArmPodK0s[iIndexHisto]->Fill(dAlpha, dPtArm);
    }
    fh2QAV0EtaEtaK0s[iIndexHisto]->Fill(trackNeg->Eta(), trackPos->Eta());
    fh2QAV0PhiPhiK0s[iIndexHisto]->Fill(trackNeg->Phi(), trackPos->Phi());
    fh1QAV0RapK0s[iIndexHisto]->Fill(vZero->RapK0Short());
  }

  if(IsCandLambda)
  {
    if(IsInPeakLambda)
    {
//          fh2QAV0EtaPtLambdaPeak[iIndexHisto]->Fill(trackNeg->Eta(),vZero->Pt());
//          fh2QAV0EtaPtLambdaPeak[iIndexHisto]->Fill(trackPos->Eta(),vZero->Pt());
      fh2QAV0EtaPtLambdaPeak[iIndexHisto]->Fill(vZero->Eta(), vZero->Pt());
      fh2QAV0PtPtLambdaPeak[iIndexHisto]->Fill(trackNeg->Pt(), trackPos->Pt());
      fh2ArmPodLambda[iIndexHisto]->Fill(dAlpha, dPtArm);
    }
    fh2QAV0EtaEtaLambda[iIndexHisto]->Fill(trackNeg->Eta(), trackPos->Eta());
    fh2QAV0PhiPhiLambda[iIndexHisto]->Fill(trackNeg->Phi(), trackPos->Phi());
    fh1QAV0RapLambda[iIndexHisto]->Fill(vZero->RapLambda());
  }

  fh2ArmPod[iIndexHisto]->Fill(dAlpha, dPtArm);

}

void AliAnalysisTaskV0sInJets::FillCandidates(Double_t mK, Double_t mL, Double_t mAL, Bool_t isK, Bool_t isL, Bool_t isAL, Int_t iCut/*cut index*/, Int_t iCent/*cent index*/)
{
  if(isK)
  {
    fh1V0CounterCentK0s[iCent]->Fill(iCut);
    fh1V0InvMassK0sAll[iCut]->Fill(mK);
  }
  if(isL)
  {
    fh1V0CounterCentLambda[iCent]->Fill(iCut);
    fh1V0InvMassLambdaAll[iCut]->Fill(mL);
  }
  if(isAL)
  {
    fh1V0CounterCentALambda[iCent]->Fill(iCut);
    fh1V0InvMassALambdaAll[iCut]->Fill(mAL);
  }
}

Bool_t AliAnalysisTaskV0sInJets::IsParticleInCone(const AliVParticle* part1, const AliVParticle* part2, Double_t dRMax) const
{
// decides whether a particle is inside a jet cone
  if(!part1 || !part2)
    return kFALSE;

  TVector3 vecMom2(part2->Px(), part2->Py(), part2->Pz());
  TVector3 vecMom1(part1->Px(), part1->Py(), part1->Pz());
  Double_t dR = vecMom2.DeltaR(vecMom1); // = sqrt(dEta*dEta+dPhi*dPhi)
  if(dR < dRMax) // momentum vectors of part1 and part2 are closer than dRMax
    return kTRUE;
  return kFALSE;
}

Bool_t AliAnalysisTaskV0sInJets::OverlapWithJets(const TClonesArray* array, const AliVParticle* part, Double_t dDistance) const
{
// decides whether a cone overlaps with other jets
  if(!part)
  {
    if(fDebug > 0) printf("AliAnalysisTaskV0sInJets::OverlapWithJets: Error: No part\n");
    return kFALSE;
  }
  if(!array)
  {
    if(fDebug > 0) printf("AliAnalysisTaskV0sInJets::OverlapWithJets: Error: No array\n");
    return kFALSE;
  }
  Int_t iNJets = array->GetEntriesFast();
  if(iNJets <= 0)
  {
    if(fDebug > 2) printf("AliAnalysisTaskV0sInJets::OverlapWithJets: Warning: No jets\n");
    return kFALSE;
  }
  AliVParticle* jet = 0;
  for(Int_t iJet = 0; iJet < iNJets; iJet++)
  {
    jet = (AliVParticle*)array->At(iJet);
    if(!jet)
    {
      if(fDebug > 0) printf("AliAnalysisTaskV0sInJets::OverlapWithJets: Error: Failed to load jet %d/%d\n", iJet, iNJets);
      continue;
    }
    if(IsParticleInCone(part, jet, dDistance))
      return kTRUE;
  }
  return kFALSE;
}

AliAODJet* AliAnalysisTaskV0sInJets::GetRandomCone(const TClonesArray* array, Double_t dEtaConeMax, Double_t dDistance) const
{
// generate a random cone which does not overlap with selected jets
//  printf("Generating random cone...\n");
  TLorentzVector vecCone;
  AliAODJet* part = 0;
  Double_t dEta, dPhi;
  Int_t iNTrialsMax = 10;
  Bool_t bStatus = kFALSE;
  for(Int_t iTry = 0; iTry < iNTrialsMax; iTry++)
  {
//      printf("Try %d\n",iTry);
    dEta = dEtaConeMax * (2 * fRandom->Rndm() - 1.); // random eta in [-dEtaConeMax,+dEtaConeMax]
    dPhi = TMath::TwoPi() * fRandom->Rndm(); // random phi in [0,2*Pi]
    vecCone.SetPtEtaPhiM(1., dEta, dPhi, 0.);
    part = new AliAODJet(vecCone);
    if(!OverlapWithJets(array, part, dDistance))
    {
      bStatus = kTRUE;
//          printf("Success\n");
      break;
    }
    else
      delete part;
  }
  if(!bStatus)
    part = 0;
  return part;
}

AliAODJet* AliAnalysisTaskV0sInJets::GetMedianCluster(const TClonesArray* array, Double_t dEtaConeMax) const
{
// sort kt clusters by pT/area and return the middle one, based on code in AliAnalysisTaskJetChem
  if(!array)
  {
    if(fDebug > 0) printf("AliAnalysisTaskV0sInJets::GetMedianCluster: Error: No array\n");
    return NULL;
  }
  Int_t iNClTot = array->GetEntriesFast(); // number of all clusters in the array
  Int_t iNCl = 0; // number of accepted clusters

  // get list of densities
  std::vector<std::vector<Double_t> > vecListClusters; // vector that contains pairs [ index, density ]
//  printf("AliAnalysisTaskV0sInJets::GetMedianCluster: Loop over %d clusters.\n", iNClTot);
  for(Int_t ij = 0; ij < iNClTot; ij++)
  {
    AliAODJet* clusterBg = (AliAODJet*)(array->At(ij));
    if(!clusterBg)
    {
//      printf("AliAnalysisTaskV0sInJets::GetMedianCluster: cluster %d/%d not ok\n", ij, iNClTot);
      return NULL;
    }
    if(TMath::Abs(clusterBg->Eta()) > 0.9 - fdRadiusJetBg)
      continue;
//    fh2EtaPhiMedCone[0]->Fill(clusterBg->Eta(), clusterBg->Phi());
//    printf("AliAnalysisTaskV0sInJets::GetMedianCluster: Cluster %d/%d used as accepted cluster %d.\n", ij, iNClTot, int(vecListClusters.size()));
    Double_t dPtBg = clusterBg->Pt();
    Double_t dAreaBg = clusterBg->EffectiveAreaCharged();
    Double_t dDensityBg = 0;
    if(dAreaBg > 0)
      dDensityBg = dPtBg / dAreaBg;
    std::vector<Double_t> vecCluster;
    vecCluster.push_back(ij);
    vecCluster.push_back(dDensityBg);
    vecListClusters.push_back(vecCluster);
  }
  iNCl = vecListClusters.size();
  if(iNCl < 3) // need at least 3 clusters (skipping 2 highest)
  {
//    if(fDebug > 2) printf("AliAnalysisTaskV0sInJets::GetMedianCluster: Warning: Too little clusters\n");
    return NULL;
  }

//  printf("AliAnalysisTaskV0sInJets::GetMedianCluster: Original lists:\n");
//  for(Int_t i = 0; i < iNCl; i++)
//    printf("%g %g\n", (vecListClusters[i])[0], (vecListClusters[i])[1]);

  // sort list of indeces by density in descending order
  std::sort(vecListClusters.begin(), vecListClusters.end(), CompareClusters);

//  printf("AliAnalysisTaskV0sInJets::GetMedianCluster: Sorted lists:\n");
//  for(Int_t i = 0; i < iNCl; i++)
//    printf("%g %g\n", (vecListClusters[i])[0], (vecListClusters[i])[1]);

  // get median cluster with median density
  AliAODJet* clusterMed = 0;
  Int_t iIndex = 0; // index of the median cluster in the sorted list
  Int_t iIndexMed = 0; // index of the median cluster in the original array
  if(TMath::Odd(iNCl))  // odd number of clusters
  {
    iIndex = (Int_t)(0.5 * (iNCl + 1)); // = (n - skip + 1)/2 + 1, skip = 2
//    printf("AliAnalysisTaskV0sInJets::GetMedianCluster: Odd, median index = %d/%d\n", iIndex, iNCl);
  }
  else // even number: picking randomly one of the two closest to median
  {
    Int_t iIndex1 = (Int_t)(0.5 * iNCl); // = (n - skip)/2 + 1, skip = 2
    Int_t iIndex2 = (Int_t)(0.5 * iNCl + 1); // = (n - skip)/2 + 1 + 1, skip = 2
    iIndex = ((fRandom->Rndm() > 0.5) ? iIndex1 : iIndex2);
//    printf("AliAnalysisTaskV0sInJets::GetMedianCluster: Even, median index = %d or %d -> %d/%d\n", iIndex1, iIndex2, iIndex, iNCl);
  }
  iIndexMed = Int_t((vecListClusters[iIndex])[0]);

//  printf("AliAnalysisTaskV0sInJets::GetMedianCluster: getting median cluster %d/%d ok, rho = %g\n", iIndexMed, iNClTot, (vecListClusters[iIndex])[1]);
  clusterMed = (AliAODJet*)(array->At(iIndexMed));

  if(TMath::Abs(clusterMed->Eta()) > dEtaConeMax)
    return NULL;

  return clusterMed;
}

Double_t AliAnalysisTaskV0sInJets::AreaCircSegment(Double_t dRadius, Double_t dDistance) const
{
// calculate area of a circular segment defined by the circle radius and the (oriented) distance between the secant line and the circle centre
  Double_t dEpsilon = 1e-2;
  Double_t dR = dRadius;
  Double_t dD = dDistance;
  if(TMath::Abs(dR) < dEpsilon)
  {
    if(fDebug > 0) printf("AliAnalysisTaskV0sInJets::AreaCircSegment: Error: Too small radius: %f < %f\n", dR, dEpsilon);
    return 0.;
  }
  if(dD >= dR)
    return 0.;
  if(dD <= -dR)
    return TMath::Pi() * dR * dR;
  return dR * dR * TMath::ACos(dD / dR) - dD * TMath::Sqrt(dR * dR - dD * dD);
}

Bool_t AliAnalysisTaskV0sInJets::IsSelectedForJets(AliAODEvent* fAOD, Double_t dVtxZCut, Double_t dVtxR2Cut, Double_t dCentCutLo, Double_t dCentCutUp, Bool_t bCutDeltaZ, Double_t dDeltaZMax)
{
// event selection
  AliAODVertex* vertex = fAOD->GetPrimaryVertex();
  if(!vertex)
    return kFALSE;
  Int_t iNContribMin = 3;
  if(!fbIsPbPb)
    iNContribMin = 2;
  if(vertex->GetNContributors() < iNContribMin)
    return kFALSE;
  TString vtxTitle(vertex->GetTitle());
  if(vtxTitle.Contains("TPCVertex"))
    return kFALSE;
  Double_t zVertex = vertex->GetZ();
  if(TMath::Abs(zVertex) > dVtxZCut)
    return kFALSE;
  if(bCutDeltaZ)
  {
    AliAODVertex* vertexSPD = fAOD->GetPrimaryVertexSPD();
    if(!vertexSPD)
    {
//          printf("IsSelectedForJets: Error: No SPD vertex\n");
      return kFALSE;
    }
    Double_t zVertexSPD = vertexSPD->GetZ();
    if(TMath::Abs(zVertex - zVertexSPD) > dDeltaZMax)
    {
//          printf("IsSelectedForJets: Rejecting event due to delta z = %f - %f = %f\n",zVertex,zVertexSPD,zVertex-zVertexSPD);
      return kFALSE;
    }
//      printf("IsSelectedForJets: Event OK: %f - %f = %f\n",zVertex,zVertexSPD,zVertex-zVertexSPD);
  }
  Double_t xVertex = vertex->GetX();
  Double_t yVertex = vertex->GetY();
  Double_t radiusSq = yVertex * yVertex + xVertex * xVertex;
  if(radiusSq > dVtxR2Cut)
    return kFALSE;
  Double_t centrality;
//  centrality = fAOD->GetHeader()->GetCentrality();
  centrality = ((AliVAODHeader*)fAOD->GetHeader())->GetCentralityP()->GetCentralityPercentile("V0M");
  if(fbIsPbPb)
  {
    if(centrality < 0)
      return kFALSE;
    if((dCentCutUp < 0) || (dCentCutLo < 0) || (dCentCutUp > 100) || (dCentCutLo > 100) || (dCentCutLo > dCentCutUp))
      return kFALSE;
    if((centrality < dCentCutLo) || (centrality > dCentCutUp))
      return kFALSE;
  }
  else
  {
    if(centrality != -1)
      return kFALSE;
  }
  return kTRUE;
}

Int_t AliAnalysisTaskV0sInJets::GetCentralityBinIndex(Double_t centrality)
{
// returns index of the centrality bin corresponding to the provided value of centrality
  if(centrality < 0 || centrality > fgkiCentBinRanges[fgkiNBinsCent - 1])
    return -1;
  for(Int_t i = 0; i < fgkiNBinsCent; i++)
  {
    if(centrality <= fgkiCentBinRanges[i])
      return i;
  }
  return -1;
}

Int_t AliAnalysisTaskV0sInJets::GetCentralityBinEdge(Int_t index)
{
// returns the upper edge of the centrality bin corresponding to the provided value of index
  if(index < 0 || index >= fgkiNBinsCent)
    return -1;
  return fgkiCentBinRanges[index];
}

TString AliAnalysisTaskV0sInJets::GetCentBinLabel(Int_t index)
{
// get string with centrality range for given bin
  TString lowerEdge = ((index == 0) ? "0" : Form("%d", GetCentralityBinEdge(index - 1)));
  TString upperEdge = Form("%d", GetCentralityBinEdge(index));
  return Form("%s-%s %%", lowerEdge.Data(), upperEdge.Data());
}

Double_t AliAnalysisTaskV0sInJets::MassPeakSigmaOld(Double_t pt, Int_t particle)
{
// estimation of the sigma of the invariant-mass peak as a function of pT and particle type
  switch(particle)
  {
    case 0: // K0S
      return 0.0044 + 0.0004 * (pt - 1.);
      break;
    case 1: // Lambda
      return 0.0023 + 0.00034 * (pt - 1.);
      break;
    default:
      return 0;
      break;
  }
}

bool AliAnalysisTaskV0sInJets::CompareClusters(const std::vector<Double_t> cluster1, const std::vector<Double_t> cluster2)
{
  return (cluster1[1] > cluster2[1]);
}
