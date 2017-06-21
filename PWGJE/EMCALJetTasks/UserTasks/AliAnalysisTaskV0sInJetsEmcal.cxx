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

//-------------------------------------------------------------------------
//     task for analysis of V0s (K0S, (anti-)Lambda) in charged jets
//     fork of AliAnalysisTaskV0sInJets for the EMCal framework
//     Author: Vit Kucera (vit.kucera@cern.ch)
//-------------------------------------------------------------------------

#include <TChain.h>
#include <TH1D.h>
#include <TH2D.h>
#include <THnSparse.h>
#include <TClonesArray.h>
#include <TRandom3.h>
#include <TDatabasePDG.h>
#include <TPDGCode.h>

#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliAODEvent.h"
#include "AliAODMCHeader.h"
#include "AliLog.h"
#include "AliEventPoolManager.h"
#include "AliPIDResponse.h"
#include "AliAODTrack.h"
#include "AliAODMCParticle.h"
#include "AliMCEvent.h"

#include "AliEmcalJet.h"
#include "AliJetContainer.h"
#include "AliParticleContainer.h"
//#include "AliClusterContainer.h"
//#include "AliPicoTrack.h"

#include "AliAnalysisTaskV0sInJetsEmcal.h"

ClassImp(AliAnalysisTaskV0sInJetsEmcal)

// upper edges of centrality bins
//const Int_t AliAnalysisTaskV0sInJetsEmcal::fgkiCentBinRanges[] = {10, 30, 50, 80}; // Alice Zimmermann
//const Int_t AliAnalysisTaskV0sInJetsEmcal::fgkiCentBinRanges[] = {10, 20, 40, 60, 80}; // Vit Kucera, initial binning
//const Int_t AliAnalysisTaskV0sInJetsEmcal::fgkiCentBinRanges[] = {5, 10, 20, 40, 60, 80}; // Iouri Belikov, LF analysis
const Int_t AliAnalysisTaskV0sInJetsEmcal::fgkiCentBinRanges[] = {10}; // only central
// edges of centrality bins for event mixing
Double_t AliAnalysisTaskV0sInJetsEmcal::fgkiCentMixBinRanges[] = {0, 5, 10}; // only central
// edges of z_vtx bins for event mixing
Double_t AliAnalysisTaskV0sInJetsEmcal::fgkiZVtxMixBinRanges[] = { -10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10};
// axis: pT of V0
const Double_t AliAnalysisTaskV0sInJetsEmcal::fgkdBinsPtV0[] = {0, 12};
const Int_t AliAnalysisTaskV0sInJetsEmcal::fgkiNBinsPtV0 = sizeof(AliAnalysisTaskV0sInJetsEmcal::fgkdBinsPtV0) / sizeof((AliAnalysisTaskV0sInJetsEmcal::fgkdBinsPtV0)[0]) - 1;
const Int_t AliAnalysisTaskV0sInJetsEmcal::fgkiNBinsPtV0Init = int(((AliAnalysisTaskV0sInJetsEmcal::fgkdBinsPtV0)[AliAnalysisTaskV0sInJetsEmcal::fgkiNBinsPtV0] - (AliAnalysisTaskV0sInJetsEmcal::fgkdBinsPtV0)[0]) / 0.1); // bin width 0.1 GeV/c
const Int_t AliAnalysisTaskV0sInJetsEmcal::fgkiNBinsPtV0InitInJet = int(((AliAnalysisTaskV0sInJetsEmcal::fgkdBinsPtV0)[AliAnalysisTaskV0sInJetsEmcal::fgkiNBinsPtV0] - (AliAnalysisTaskV0sInJetsEmcal::fgkdBinsPtV0)[0]) / 0.5); // bin width 0.5 GeV/c
// axis: pT of jets
const Double_t AliAnalysisTaskV0sInJetsEmcal::fgkdBinsPtJet[] = {0, 100};
const Int_t AliAnalysisTaskV0sInJetsEmcal::fgkiNBinsPtJet = sizeof(AliAnalysisTaskV0sInJetsEmcal::fgkdBinsPtJet) / sizeof(AliAnalysisTaskV0sInJetsEmcal::fgkdBinsPtJet[0]) - 1;
//const Int_t AliAnalysisTaskV0sInJetsEmcal::fgkiNBinsPtJetInit = int(((AliAnalysisTaskV0sInJetsEmcal::fgkdBinsPtJet)[AliAnalysisTaskV0sInJetsEmcal::fgkiNBinsPtJet] - (AliAnalysisTaskV0sInJetsEmcal::fgkdBinsPtJet)[0]) / 5.); // bin width 5 GeV/c
const Int_t AliAnalysisTaskV0sInJetsEmcal::fgkiNBinsPtJetInit = int(((AliAnalysisTaskV0sInJetsEmcal::fgkdBinsPtJet)[AliAnalysisTaskV0sInJetsEmcal::fgkiNBinsPtJet] - (AliAnalysisTaskV0sInJetsEmcal::fgkdBinsPtJet)[0]) / 10.); // bin width 10 GeV/c
// axis: K0S invariant mass
const Int_t AliAnalysisTaskV0sInJetsEmcal::fgkiNBinsMassK0s = 300;
const Double_t AliAnalysisTaskV0sInJetsEmcal::fgkdMassK0sMin = 0.35; // [GeV/c^2]
const Double_t AliAnalysisTaskV0sInJetsEmcal::fgkdMassK0sMax = 0.65; // [GeV/c^2]
// axis: Lambda invariant mass
const Int_t AliAnalysisTaskV0sInJetsEmcal::fgkiNBinsMassLambda = 200;
const Double_t AliAnalysisTaskV0sInJetsEmcal::fgkdMassLambdaMin = 1.05; // [GeV/c^2]
const Double_t AliAnalysisTaskV0sInJetsEmcal::fgkdMassLambdaMax = 1.25; // [GeV/c^2]
// delta phi range
const Double_t AliAnalysisTaskV0sInJetsEmcal::fgkdDeltaPhiMin = -TMath::PiOver2(); // minimum delta-phi_V0-jet
const Double_t AliAnalysisTaskV0sInJetsEmcal::fgkdDeltaPhiMax = 3.*TMath::PiOver2(); // maximum delta-phi_V0-jet


// Default constructor
AliAnalysisTaskV0sInJetsEmcal::AliAnalysisTaskV0sInJetsEmcal():
  AliAnalysisTaskEmcalJet(),
  fAODIn(0),
  fAODOut(0),
  fEventMC(0),
  fRandom(0),
  fPoolMgr(0),
  fOutputListStd(0),
  fOutputListQA(0),
  fOutputListCuts(0),
  fOutputListMC(0),

  fbIsPbPb(0),
  fbMCAnalysis(0),
  fsGeneratorName(""),

  fdCutVertexZ(0),
  fdCutVertexR2(0),
  fdCutCentLow(0),
  fdCutCentHigh(0),
  fdCutDeltaZMax(0),
  fdCentrality(0),

  fbCorrelations(0),
  fiSizePool(0),
  fiNJetsPerPool(0),
  ffFractionMin(0),
  fiNEventsMin(0),
  fdDeltaEtaMax(0),

  fbTPCRefit(0),
  fbRejectKinks(0),
  fbFindableClusters(0),
  fdCutNCrossedRowsTPCMin(0),
  fdCutCrossedRowsOverFindMin(0),
  fdCutCrossedRowsOverFindMax(0),
  fdCutPtDaughterMin(0),
  fdCutDCAToPrimVtxMin(0),
  fdCutDCADaughtersMax(0),
  fdCutEtaDaughterMax(0),
  fdCutNSigmadEdxMax(0),
  fdPtProtonPIDMax(0),
  fbOnFly(0),
  fdCutCPAKMin(0),
  fdCutCPALMin(0),
  fdCutRadiusDecayMin(0),
  fdCutRadiusDecayMax(0),
  fdCutEtaV0Max(0),
  fdCutRapV0Max(0),
  fdCutNTauKMax(0),
  fdCutNTauLMax(0),
  fbCutArmPod(0),
  fbCutCross(0),

  fbJetSelection(0),
  fdCutPtJetMin(0),
  fdCutPtTrackJetMin(0),
  fdCutAreaPercJetMin(0),
  fdDistanceV0JetMax(0.4),
  fiBgSubtraction(0),

  fbCompareTriggers(0),

  fJetsCont(0),
  fJetsBgCont(0),
  fTracksCont(0),

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
    fh2QAV0PtNCls[i] = 0;
    fh2QAV0PtChi[i] = 0;
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

    fh2QAV0PhiPtK0sPeak[i] = 0;
    fh2QAV0PhiPtLambdaPeak[i] = 0;
    fh2QAV0PhiPtALambdaPeak[i] = 0;
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
    fhnPtDaughterK0s[i] = 0;
    fhnPtDaughterLambda[i] = 0;
    fhnPtDaughterALambda[i] = 0;
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
    fh1V0ALambdaPtMCGen[i] = 0;
    fh2V0ALambdaPtMassMCRec[i] = 0;
    fh1V0ALambdaPtMCRecFalse[i] = 0;
    fh2V0ALambdaEtaPtMCGen[i] = 0;
    fh3V0ALambdaEtaPtMassMCRec[i] = 0;
    fh2V0ALambdaInJetPtMCGen[i] = 0;
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
    fhnV0K0sInclDaughterEtaPtPtMCRec[i] = 0;
    fhnV0K0sInJetsDaughterEtaPtPtMCRec[i] = 0;
    fhnV0LambdaInclDaughterEtaPtPtMCRec[i] = 0;
    fhnV0LambdaInJetsDaughterEtaPtPtMCRec[i] = 0;
    fhnV0ALambdaInclDaughterEtaPtPtMCRec[i] = 0;
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
    // correlations
    fhnV0CorrelSEK0s[i] = 0;
    fhnV0CorrelMEK0s[i] = 0;
    fhnV0CorrelSELambda[i] = 0;
    fhnV0CorrelMELambda[i] = 0;

    fh2V0PtJetAngleK0s[i] = 0;
    fh2V0PtJetAngleLambda[i] = 0;
    fh2V0PtJetAngleALambda[i] = 0;

    fh1PtJet[i] = 0;
    fh1EtaJet[i] = 0;
    fh2EtaPtJet[i] = 0;
    fh1PhiJet[i] = 0;
    fh2PtJetPtTrackLeading[i] = 0;
    fh2PtJetPtTrigger[i] = 0;
    fh1PtTrigger[i] = 0;
    fh1NJetPerEvent[i] = 0;
    fh2EtaPhiRndCone[i] = 0;
    fh2EtaPhiMedCone[i] = 0;
    fh1DistanceJets[i] = 0;
    fh1DistanceV0JetsK0s[i] = 0;
    fh1DistanceV0JetsLambda[i] = 0;
    fh1DistanceV0JetsALambda[i] = 0;

    fh1VtxZ[i] = 0;
    fh1VtxZME[i] = 0;
    fh2VtxXY[i] = 0;
  }
}

// Constructor
AliAnalysisTaskV0sInJetsEmcal::AliAnalysisTaskV0sInJetsEmcal(const char* name):
  AliAnalysisTaskEmcalJet(name, kTRUE),
  fAODIn(0),
  fAODOut(0),
  fEventMC(0),
  fRandom(0),
  fPoolMgr(0),
  fOutputListStd(0),
  fOutputListQA(0),
  fOutputListCuts(0),
  fOutputListMC(0),

  fbIsPbPb(0),
  fbMCAnalysis(0),
  fsGeneratorName(""),

  fdCutVertexZ(0),
  fdCutVertexR2(0),
  fdCutCentLow(0),
  fdCutCentHigh(0),
  fdCutDeltaZMax(0),
  fdCentrality(0),

  fbCorrelations(0),
  fiSizePool(0),
  fiNJetsPerPool(0),
  ffFractionMin(0),
  fiNEventsMin(0),
  fdDeltaEtaMax(0),

  fbTPCRefit(0),
  fbRejectKinks(0),
  fbFindableClusters(0),
  fdCutNCrossedRowsTPCMin(0),
  fdCutCrossedRowsOverFindMin(0),
  fdCutCrossedRowsOverFindMax(0),
  fdCutPtDaughterMin(0),
  fdCutDCAToPrimVtxMin(0),
  fdCutDCADaughtersMax(0),
  fdCutEtaDaughterMax(0),
  fdCutNSigmadEdxMax(0),
  fdPtProtonPIDMax(0),
  fbOnFly(0),
  fdCutCPAKMin(0),
  fdCutCPALMin(0),
  fdCutRadiusDecayMin(0),
  fdCutRadiusDecayMax(0),
  fdCutEtaV0Max(0),
  fdCutRapV0Max(0),
  fdCutNTauKMax(0),
  fdCutNTauLMax(0),
  fbCutArmPod(0),
  fbCutCross(0),

  fbJetSelection(0),
  fdCutPtJetMin(0),
  fdCutPtTrackJetMin(0),
  fdCutAreaPercJetMin(0),
  fdDistanceV0JetMax(0.4),
  fiBgSubtraction(0),

  fbCompareTriggers(0),

  fJetsCont(0),
  fJetsBgCont(0),
  fTracksCont(0),

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
    fh2QAV0PtNCls[i] = 0;
    fh2QAV0PtChi[i] = 0;
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

    fh2QAV0PhiPtK0sPeak[i] = 0;
    fh2QAV0PhiPtLambdaPeak[i] = 0;
    fh2QAV0PhiPtALambdaPeak[i] = 0;
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
    fhnPtDaughterK0s[i] = 0;
    fhnPtDaughterLambda[i] = 0;
    fhnPtDaughterALambda[i] = 0;
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
    fh1V0ALambdaPtMCGen[i] = 0;
    fh2V0ALambdaPtMassMCRec[i] = 0;
    fh1V0ALambdaPtMCRecFalse[i] = 0;
    fh2V0ALambdaEtaPtMCGen[i] = 0;
    fh3V0ALambdaEtaPtMassMCRec[i] = 0;
    fh2V0ALambdaInJetPtMCGen[i] = 0;
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
    fhnV0K0sInclDaughterEtaPtPtMCRec[i] = 0;
    fhnV0K0sInJetsDaughterEtaPtPtMCRec[i] = 0;
    fhnV0LambdaInclDaughterEtaPtPtMCRec[i] = 0;
    fhnV0LambdaInJetsDaughterEtaPtPtMCRec[i] = 0;
    fhnV0ALambdaInclDaughterEtaPtPtMCRec[i] = 0;
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
    // correlations
    fhnV0CorrelSEK0s[i] = 0;
    fhnV0CorrelMEK0s[i] = 0;
    fhnV0CorrelSELambda[i] = 0;
    fhnV0CorrelMELambda[i] = 0;

    fh2V0PtJetAngleK0s[i] = 0;
    fh2V0PtJetAngleLambda[i] = 0;
    fh2V0PtJetAngleALambda[i] = 0;

    fh1PtJet[i] = 0;
    fh1EtaJet[i] = 0;
    fh2EtaPtJet[i] = 0;
    fh1PhiJet[i] = 0;
    fh2PtJetPtTrackLeading[i] = 0;
    fh2PtJetPtTrigger[i] = 0;
    fh1PtTrigger[i] = 0;
    fh1NJetPerEvent[i] = 0;
    fh2EtaPhiRndCone[i] = 0;
    fh2EtaPhiMedCone[i] = 0;
    fh1DistanceJets[i] = 0;
    fh1DistanceV0JetsK0s[i] = 0;
    fh1DistanceV0JetsLambda[i] = 0;
    fh1DistanceV0JetsALambda[i] = 0;

    fh1VtxZ[i] = 0;
    fh1VtxZME[i] = 0;
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
}

AliAnalysisTaskV0sInJetsEmcal::~AliAnalysisTaskV0sInJetsEmcal()
{
  delete fRandom;
  fRandom = 0;
}

void AliAnalysisTaskV0sInJetsEmcal::UserCreateOutputObjects()
{
  // Called once

  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

  if(fbJetSelection)
  {
    fJetsCont = GetJetContainer(0);
    fJetsBgCont = GetJetContainer(1);
    if(fbCompareTriggers)
    {
      fTracksCont = dynamic_cast<AliTrackContainer*>(GetParticleContainer(0));
    }
  }

  // Initialise random-number generator
  fRandom = new TRandom3(0);

  // Create event pool manager
  if(fbCorrelations)
  {
    fPoolMgr = new AliEventPoolManager(fiSizePool, fiNJetsPerPool, fgkiNBinsCentMix, fgkiCentMixBinRanges, fgkiNBinsZVtxMix, fgkiZVtxMixBinRanges);
    fPoolMgr->SetTargetValues(fiNJetsPerPool, ffFractionMin, fiNEventsMin);
  }

  // Create histograms

  fOutputListStd = new TList();
  fOutputListStd->SetOwner();
  fOutputListQA = new TList();
  fOutputListQA->SetOwner();
  fOutputListCuts = new TList();
  fOutputListCuts->SetOwner();
  fOutputListMC = new TList();
  fOutputListMC->SetOwner();

  // event categories
  const Int_t iNCategEvent = 12;
  TString categEvent[iNCategEvent] = {
    "coll. candid.",
    "AOD OK",
    "no SPD pile-up",
    "PV contrib.",
    "not TPC PV",
    "|z| PV",
    "|#Deltaz| SPD PV",
    "r PV",
    "cent.", //2
    "with V0",
    "with jets",
    "jet selection"
  };

  // labels for stages of V0 selection
  TString categV0[fgkiNCategV0] = {
    "all"/*0*/,
    "mass range"/*1*/,
    "rec. method"/*2*/,
    "tracks TPC"/*3*/,
    "track pt"/*4*/,
    "DCA prim v"/*5*/,
    "DCA daughters"/*6*/,
    "CPA"/*7*/,
    "volume"/*8*/,
    "track #it{#eta}"/*9*/,
    "V0 #it{y} & #it{#eta}"/*10*/,
    "lifetime"/*11*/,
    "PID"/*12*/,
    "Arm.-Pod."/*13*/,
    "cross-cont."/*14*/,
    "inclusive"/*15*/,
    "in jet event"/*16*/,
    "in jet"/*17*/
  };

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

  fh1V0CandPerEvent = new TH1D("fh1V0CandPerEvent", "Number of all V0 candidates per event;candidates;events", 300, 0, 30000);
  fOutputListStd->Add(fh1V0CandPerEvent);

  for(Int_t i = 0; i < fgkiNBinsCent; i++)
  {
    fh1EventCounterCutCent[i] = new TH1D(Form("fh1EventCounterCutCent_%d", i), Form("Number of events after filtering, cent %s;selection filter;counts", GetCentBinLabel(i).Data()), iNCategEvent, 0, iNCategEvent);
    for(Int_t j = 0; j < iNCategEvent; j++)
      fh1EventCounterCutCent[i]->GetXaxis()->SetBinLabel(j + 1, categEvent[j].Data());
    fh1V0CandPerEventCentK0s[i] = new TH1D(Form("fh1V0CandPerEventCentK0s_%d", i), Form("Number of selected K0s candidates per event, cent %s;candidates;events", GetCentBinLabel(i).Data()), 100, 0, 200);
    fh1V0CandPerEventCentLambda[i] = new TH1D(Form("fh1V0CandPerEventCentLambda_%d", i), Form("Number of selected Lambda candidates per event, cent %s;candidates;events", GetCentBinLabel(i).Data()), 100, 0, 200);
    fh1V0CandPerEventCentALambda[i] = new TH1D(Form("fh1V0CandPerEventCentALambda_%d", i), Form("Number of selected ALambda candidates per event, cent %s;candidates;events", GetCentBinLabel(i).Data()), 100, 0, 200);
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
  Int_t iNBinsPtV0InJet = fgkiNBinsPtV0InitInJet;
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
  fh3CCMassCorrelBoth = new THnSparseD("fh3CCMassCorrelBoth", "Mass correlation: K0S && Lambda;m K0S;m Lambda;pT", 3, binsCorrel, xminCorrel, xmaxCorrel);
  fh3CCMassCorrelKNotL = new THnSparseD("fh3CCMassCorrelKNotL", "Mass correlation: K0S, not Lambda;m K0S;m Lambda;pT", 3, binsCorrel, xminCorrel, xmaxCorrel);
  fh3CCMassCorrelLNotK = new THnSparseD("fh3CCMassCorrelLNotK", "Mass correlation: Lambda, not K0S;m K0S;m Lambda;pT", 3, binsCorrel, xminCorrel, xmaxCorrel);
  fOutputListQA->Add(fh2CCK0s);
  fOutputListQA->Add(fh2CCLambda);
  fOutputListQA->Add(fh3CCMassCorrelBoth);
  fOutputListQA->Add(fh3CCMassCorrelKNotL);
  fOutputListQA->Add(fh3CCMassCorrelLNotK);

  Double_t dStepEtaV0 = 0.05;
  Double_t dRangeEtaV0Max = 0.8;
  const Int_t iNBinsEtaV0 = 2 * Int_t(dRangeEtaV0Max / dStepEtaV0);
//  printf("%s: %d\n", "iNBinsEtaV0", iNBinsEtaV0);
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
  Int_t binsKInJC[iNDimInJC] = {fgkiNBinsMassK0s, iNBinsPtV0InJet, iNBinsEtaV0, iNJetPtBins};
  Double_t xminKInJC[iNDimInJC] = {fgkdMassK0sMin, dPtV0Min, -dRangeEtaV0Max, dJetPtMin};
  Double_t xmaxKInJC[iNDimInJC] = {fgkdMassK0sMax, dPtV0Max, dRangeEtaV0Max, dJetPtMax};
  Int_t binsLInJC[iNDimInJC] = {fgkiNBinsMassLambda, iNBinsPtV0InJet, iNBinsEtaV0, iNJetPtBins};
  Double_t xminLInJC[iNDimInJC] = {fgkdMassLambdaMin, dPtV0Min, -dRangeEtaV0Max, dJetPtMin};
  Double_t xmaxLInJC[iNDimInJC] = {fgkdMassLambdaMax, dPtV0Max, dRangeEtaV0Max, dJetPtMax};
  // binning in V0-jet correlations
  Int_t iNBinsDeltaPhi = 60;
  const Int_t iNDimCorrel = 6;
//  Double_t fdDeltaEtaMax = 2. * 0.7; // >= 2*eta_V0^max
//  Int_t iNBinsDeltaEtaCorrel = 2 * Int_t(fdDeltaEtaMax / 0.1); // does not work if fdDeltaEtaMax is Float_t
//  Int_t iNBinsDeltaEtaCorrel = 2 * Int_t(fdDeltaEtaMax / 0.7); // does not work if fdDeltaEtaMax is Float_t
//  Int_t iNBinsDeltaEtaCorrel = 2 * Int_t(fdDeltaEtaMax / fdDeltaEtaMax); // does not work if fdDeltaEtaMax is Float_t
  Int_t iNBinsDeltaEtaCorrel = 1;
//  printf("%s: %d\n", "iNBinsDeltaEtaCorrel", iNBinsDeltaEtaCorrel);
  Int_t binsKCorrel[iNDimCorrel] = {fgkiNBinsMassK0s / 2, iNBinsPtV0InJet, iNBinsEtaV0, iNJetPtBins, iNBinsDeltaPhi, iNBinsDeltaEtaCorrel};
  Double_t xminKCorrel[iNDimCorrel] = {fgkdMassK0sMin, dPtV0Min, -dRangeEtaV0Max, dJetPtMin, fgkdDeltaPhiMin, -fdDeltaEtaMax};
  Double_t xmaxKCorrel[iNDimCorrel] = {fgkdMassK0sMax, dPtV0Max, dRangeEtaV0Max, dJetPtMax, fgkdDeltaPhiMax, fdDeltaEtaMax};
  Int_t binsLCorrel[iNDimCorrel] = {fgkiNBinsMassLambda / 2, iNBinsPtV0InJet, iNBinsEtaV0, iNJetPtBins, iNBinsDeltaPhi, iNBinsDeltaEtaCorrel};
  Double_t xminLCorrel[iNDimCorrel] = {fgkdMassLambdaMin, dPtV0Min, -dRangeEtaV0Max, dJetPtMin, fgkdDeltaPhiMin, -fdDeltaEtaMax};
  Double_t xmaxLCorrel[iNDimCorrel] = {fgkdMassLambdaMax, dPtV0Max, dRangeEtaV0Max, dJetPtMax, fgkdDeltaPhiMax, fdDeltaEtaMax};
//  printf("binsKCorrel: %d %d %d %d %d %d\n", fgkiNBinsMassK0s / 2, iNBinsPtV0InJet, iNBinsEtaV0, iNJetPtBins, iNBinsDeltaPhi, iNBinsDeltaEtaCorrel);

  // binning eff inclusive vs eta-pT
  Double_t dStepDeltaEta = 0.1;
  Double_t dRangeDeltaEtaMax = 0.5;
  const Int_t iNBinsDeltaEta = 2 * Int_t(dRangeDeltaEtaMax / dStepDeltaEta);
//  printf("%s: %d\n", "iNBinsDeltaEta", iNBinsDeltaEta);
  Int_t binsEtaK[3] = {fgkiNBinsMassK0s, iNBinsPtV0, iNBinsEtaV0};
  Double_t xminEtaK[3] = {fgkdMassK0sMin, dPtV0Min, -dRangeEtaV0Max};
  Double_t xmaxEtaK[3] = {fgkdMassK0sMax, dPtV0Max, dRangeEtaV0Max};
  Int_t binsEtaL[3] = {fgkiNBinsMassLambda, iNBinsPtV0, iNBinsEtaV0};
  Double_t xminEtaL[3] = {fgkdMassLambdaMin, dPtV0Min, -dRangeEtaV0Max};
  Double_t xmaxEtaL[3] = {fgkdMassLambdaMax, dPtV0Max, dRangeEtaV0Max};
  // binning eff in jets vs eta-pT
  // associated
  Int_t binsEtaKInRec[5] = {fgkiNBinsMassK0s, iNBinsPtV0InJet, iNBinsEtaV0, iNJetPtBins, iNBinsDeltaEta};
  Double_t xminEtaKInRec[5] = {fgkdMassK0sMin, dPtV0Min, -dRangeEtaV0Max, dJetPtMin, -dRangeDeltaEtaMax};
  Double_t xmaxEtaKInRec[5] = {fgkdMassK0sMax, dPtV0Max, dRangeEtaV0Max, dJetPtMax, dRangeDeltaEtaMax};
  Int_t binsEtaLInRec[5] = {fgkiNBinsMassLambda, iNBinsPtV0InJet, iNBinsEtaV0, iNJetPtBins, iNBinsDeltaEta};
  Double_t xminEtaLInRec[5] = {fgkdMassLambdaMin, dPtV0Min, -dRangeEtaV0Max, dJetPtMin, -dRangeDeltaEtaMax};
  Double_t xmaxEtaLInRec[5] = {fgkdMassLambdaMax, dPtV0Max, dRangeEtaV0Max, dJetPtMax, dRangeDeltaEtaMax};
  // generated
  Int_t binsEtaInGen[4] = {iNBinsPtV0InJet, iNBinsEtaV0, iNJetPtBins, iNBinsDeltaEta};
  Double_t xminEtaInGen[4] = {dPtV0Min, -dRangeEtaV0Max, dJetPtMin, -dRangeDeltaEtaMax};
  Double_t xmaxEtaInGen[4] = {dPtV0Max, dRangeEtaV0Max, dJetPtMax, dRangeDeltaEtaMax};
  // daughter eta: charge-etaD-ptD-etaV0-ptV0-ptJet
  const Int_t iNDimEtaD = 6;
  Int_t binsEtaDaughter[iNDimEtaD] = {2, 20, iNBinsPtV0, iNBinsEtaV0, iNBinsPtV0, iNJetPtBins};
  Double_t xminEtaDaughter[iNDimEtaD] = {0, -1, dPtV0Min, -dRangeEtaV0Max, dPtV0Min, dJetPtMin};
  Double_t xmaxEtaDaughter[iNDimEtaD] = {2, 1, dPtV0Max, dRangeEtaV0Max, dPtV0Max, dJetPtMax};
  // daughter track pt: pos-neg-ptV0-ptJet-ptLead
  const Int_t iNDimDaughter = 5;
  Int_t binsDaughter[iNDimDaughter] = {80, 80, iNBinsPtV0, iNJetPtBins, 80};
  Double_t xminDaughter[iNDimDaughter] = {0, 0, dPtV0Min, dJetPtMin, 0};
  Double_t xmaxDaughter[iNDimDaughter] = {20, 20, dPtV0Max, dJetPtMax, 20};

  for(Int_t i = 0; i < fgkiNBinsCent; i++)
  {
    // Daughter track Pt vs pt
    fhnPtDaughterK0s[i] = new THnSparseD(Form("fhnPtDaughterK0s_%d", i), Form("K0s: Daughter Pt vs Pt, cent: %s;neg #it{p}_{T} (GeV/#it{c});pos #it{p}_{T} (GeV/#it{c});p_{T}^{V0} (GeV/#it{c});#it{p}_{T}^{jet} (GeV/#it{c});#it{p}_{T} leading track (GeV/#it{c});counts", GetCentBinLabel(i).Data()), iNDimDaughter, binsDaughter, xminDaughter, xmaxDaughter);
    fhnPtDaughterLambda[i] = new THnSparseD(Form("fhnPtDaughterLambda_%d", i), Form("Lambda: Daughter Pt vs Pt, cent: %s;neg #it{p}_{T} (GeV/#it{c});pos #it{p}_{T} (GeV/#it{c});p_{T}^{V0} (GeV/#it{c});#it{p}_{T}^{jet} (GeV/#it{c});#it{p}_{T} leading track (GeV/#it{c});counts", GetCentBinLabel(i).Data()), iNDimDaughter, binsDaughter, xminDaughter, xmaxDaughter);
    fhnPtDaughterALambda[i] = new THnSparseD(Form("fhnPtDaughterALambda_%d", i), Form("ALambda: Daughter Pt vs Pt, cent: %s;neg #it{p}_{T} (GeV/#it{c});pos #it{p}_{T} (GeV/#it{c});p_{T}^{V0} (GeV/#it{c});#it{p}_{T}^{jet} (GeV/#it{c});#it{p}_{T} leading track (GeV/#it{c});counts", GetCentBinLabel(i).Data()), iNDimDaughter, binsDaughter, xminDaughter, xmaxDaughter);
    fOutputListStd->Add(fhnPtDaughterK0s[i]);
    fOutputListStd->Add(fhnPtDaughterLambda[i]);
    fOutputListStd->Add(fhnPtDaughterALambda[i]);

    fh1V0InvMassK0sCent[i] = new TH1D(Form("fh1V0InvMassK0sCent_%d", i), Form("K0s: V0 invariant mass, cent %s;#it{m}_{inv} (GeV/#it{c}^{2});counts", GetCentBinLabel(i).Data()), fgkiNBinsMassK0s, fgkdMassK0sMin, fgkdMassK0sMax);
    fh1V0InvMassLambdaCent[i] = new TH1D(Form("fh1V0InvMassLambdaCent_%d", i), Form("Lambda: V0 invariant mass, cent %s;#it{m}_{inv} (GeV/#it{c}^{2});counts", GetCentBinLabel(i).Data()), fgkiNBinsMassLambda, fgkdMassLambdaMin, fgkdMassLambdaMax);
    fh1V0InvMassALambdaCent[i] = new TH1D(Form("fh1V0InvMassALambdaCent_%d", i), Form("ALambda: V0 invariant mass, cent %s;#it{m}_{inv} (GeV/#it{c}^{2});counts", GetCentBinLabel(i).Data()), fgkiNBinsMassLambda, fgkdMassLambdaMin, fgkdMassLambdaMax);
    fOutputListStd->Add(fh1V0InvMassK0sCent[i]);
    fOutputListStd->Add(fh1V0InvMassLambdaCent[i]);
    fOutputListStd->Add(fh1V0InvMassALambdaCent[i]);
    // Inclusive
    fhnV0InclusiveK0s[i] = new THnSparseD(Form("fhnV0InclusiveK0s_C%d", i), "K0s: V0 invariant mass vs pt;#it{m}_{inv} (GeV/#it{c}^{2});#it{p}_{T}^{V0} (GeV/#it{c});#it{#eta}_{V0};counts", iNDimIncl, binsKIncl, xminKIncl, xmaxKIncl);
    fhnV0InclusiveLambda[i] = new THnSparseD(Form("fhnV0InclusiveLambda_C%d", i), "Lambda: V0 invariant mass vs pt;#it{m}_{inv} (GeV/#it{c}^{2});#it{p}_{T}^{V0} (GeV/#it{c});#it{#eta}_{V0};counts", iNDimIncl, binsLIncl, xminLIncl, xmaxLIncl);
    fhnV0InclusiveALambda[i] = new THnSparseD(Form("fhnV0InclusiveALambda_C%d", i), "ALambda: V0 invariant mass vs pt;#it{m}_{inv} (GeV/#it{c}^{2});#it{p}_{T}^{V0} (GeV/#it{c});#it{#eta}_{V0};counts", iNDimIncl, binsLIncl, xminLIncl, xmaxLIncl);
    fOutputListStd->Add(fhnV0InclusiveK0s[i]);
    fOutputListStd->Add(fhnV0InclusiveLambda[i]);
    fOutputListStd->Add(fhnV0InclusiveALambda[i]);
    // In cones
    fhnV0InJetK0s[i] = new THnSparseD(Form("fhnV0InJetK0s_%d", i), Form("K0s: Mass vs Pt in jets, cent: %s;#it{m}_{inv} (GeV/#it{c}^{2});#it{p}_{T}^{V0} (GeV/#it{c});#it{#eta}_{V0};#it{p}_{T}^{jet} (GeV/#it{c})", GetCentBinLabel(i).Data()), iNDimInJC, binsKInJC, xminKInJC, xmaxKInJC);
    fOutputListStd->Add(fhnV0InJetK0s[i]);
    fhnV0InPerpK0s[i] = new THnSparseD(Form("fhnV0InPerpK0s_%d", i), Form("K0s: Mass vs Pt in perp. cones, cent: %s;#it{m}_{inv} (GeV/#it{c}^{2});#it{p}_{T}^{V0} (GeV/#it{c});#it{#eta}_{V0};#it{p}_{T}^{jet} (GeV/#it{c})", GetCentBinLabel(i).Data()), iNDimInJC, binsKInJC, xminKInJC, xmaxKInJC);
    fOutputListStd->Add(fhnV0InPerpK0s[i]);
    fhnV0InRndK0s[i] = new THnSparseD(Form("fhnV0InRndK0s_%d", i), Form("K0s: Mass vs Pt in rnd. cones, cent: %s;#it{m}_{inv} (GeV/#it{c}^{2});#it{p}_{T}^{V0} (GeV/#it{c});#it{#eta}_{V0}", GetCentBinLabel(i).Data()), iNDimIncl, binsKIncl, xminKIncl, xmaxKIncl);
    fOutputListStd->Add(fhnV0InRndK0s[i]);
    fhnV0InMedK0s[i] = new THnSparseD(Form("fhnV0InMedK0s_%d", i), Form("K0s: Mass vs Pt in med.-cl. cones, cent: %s;#it{m}_{inv} (GeV/#it{c}^{2});#it{p}_{T}^{V0} (GeV/#it{c});#it{#eta}_{V0}", GetCentBinLabel(i).Data()), iNDimIncl, binsKIncl, xminKIncl, xmaxKIncl);
    fOutputListStd->Add(fhnV0InMedK0s[i]);
    fhnV0OutJetK0s[i] = new THnSparseD(Form("fhnV0OutJetK0s_%d", i), Form("K0s: Pt outside jets, cent: %s;#it{m}_{inv} (GeV/#it{c}^{2});#it{p}_{T}^{V0} (GeV/#it{c});#it{#eta}_{V0}", GetCentBinLabel(i).Data()), iNDimIncl, binsKIncl, xminKIncl, xmaxKIncl);
    fOutputListStd->Add(fhnV0OutJetK0s[i]);
    fhnV0NoJetK0s[i] = new THnSparseD(Form("fhnV0NoJetK0s_%d", i), Form("K0s: Pt in jet-less events, cent: %s;#it{m}_{inv} (GeV/#it{c}^{2});#it{p}_{T}^{V0} (GeV/#it{c});#it{#eta}_{V0}", GetCentBinLabel(i).Data()), iNDimIncl, binsKIncl, xminKIncl, xmaxKIncl);
    fOutputListStd->Add(fhnV0NoJetK0s[i]);
    fhnV0InJetLambda[i] = new THnSparseD(Form("fhnV0InJetLambda_%d", i), Form("Lambda: Mass vs Pt in jets, cent: %s;#it{m}_{inv} (GeV/#it{c}^{2});#it{p}_{T}^{V0} (GeV/#it{c});#it{#eta}_{V0};#it{p}_{T}^{jet} (GeV/#it{c})", GetCentBinLabel(i).Data()), iNDimInJC, binsLInJC, xminLInJC, xmaxLInJC);
    fOutputListStd->Add(fhnV0InJetLambda[i]);
    fhnV0InPerpLambda[i] = new THnSparseD(Form("fhnV0InPerpLambda_%d", i), Form("Lambda: Mass vs Pt in perp. cones, cent: %s;#it{m}_{inv} (GeV/#it{c}^{2});#it{p}_{T}^{V0} (GeV/#it{c});#it{#eta}_{V0};#it{p}_{T}^{jet} (GeV/#it{c})", GetCentBinLabel(i).Data()), iNDimInJC, binsLInJC, xminLInJC, xmaxLInJC);
    fOutputListStd->Add(fhnV0InPerpLambda[i]);
    fhnV0InRndLambda[i] = new THnSparseD(Form("fhnV0InRndLambda_%d", i), Form("Lambda: Mass vs Pt in rnd. cones, cent: %s;#it{m}_{inv} (GeV/#it{c}^{2});#it{p}_{T}^{V0} (GeV/#it{c});#it{#eta}_{V0}", GetCentBinLabel(i).Data()), iNDimIncl, binsLIncl, xminLIncl, xmaxLIncl);
    fOutputListStd->Add(fhnV0InRndLambda[i]);
    fhnV0InMedLambda[i] = new THnSparseD(Form("fhnV0InMedLambda_%d", i), Form("Lambda: Mass vs Pt in med.-cl. cones, cent: %s;#it{m}_{inv} (GeV/#it{c}^{2});#it{p}_{T}^{V0} (GeV/#it{c});#it{#eta}_{V0}", GetCentBinLabel(i).Data()), iNDimIncl, binsLIncl, xminLIncl, xmaxLIncl);
    fOutputListStd->Add(fhnV0InMedLambda[i]);
    fhnV0OutJetLambda[i] = new THnSparseD(Form("fhnV0OutJetLambda_%d", i), Form("Lambda: Pt outside jets, cent: %s;#it{m}_{inv} (GeV/#it{c}^{2});#it{p}_{T}^{V0} (GeV/#it{c});#it{#eta}_{V0}", GetCentBinLabel(i).Data()), iNDimIncl, binsLIncl, xminLIncl, xmaxLIncl);
    fOutputListStd->Add(fhnV0OutJetLambda[i]);
    fhnV0NoJetLambda[i] = new THnSparseD(Form("fhnV0NoJetLambda_%d", i), Form("Lambda: Pt in jet-less events, cent: %s;#it{m}_{inv} (GeV/#it{c}^{2});#it{p}_{T}^{V0} (GeV/#it{c});#it{#eta}_{V0}", GetCentBinLabel(i).Data()), iNDimIncl, binsLIncl, xminLIncl, xmaxLIncl);
    fOutputListStd->Add(fhnV0NoJetLambda[i]);
    fhnV0InJetALambda[i] = new THnSparseD(Form("fhnV0InJetALambda_%d", i), Form("ALambda: Mass vs Pt in jets, cent: %s;#it{m}_{inv} (GeV/#it{c}^{2});#it{p}_{T}^{V0} (GeV/#it{c});#it{#eta}_{V0};#it{p}_{T}^{jet} (GeV/#it{c})", GetCentBinLabel(i).Data()), iNDimInJC, binsLInJC, xminLInJC, xmaxLInJC);
    fOutputListStd->Add(fhnV0InJetALambda[i]);
    fhnV0InPerpALambda[i] = new THnSparseD(Form("fhnV0InPerpALambda_%d", i), Form("ALambda: Mass vs Pt in perp. cones, cent: %s;#it{m}_{inv} (GeV/#it{c}^{2});#it{p}_{T}^{V0} (GeV/#it{c});#it{#eta}_{V0};#it{p}_{T}^{jet} (GeV/#it{c})", GetCentBinLabel(i).Data()), iNDimInJC, binsLInJC, xminLInJC, xmaxLInJC);
    fOutputListStd->Add(fhnV0InPerpALambda[i]);
    fhnV0InRndALambda[i] = new THnSparseD(Form("fhnV0InRndALambda_%d", i), Form("ALambda: Mass vs Pt in rnd. cones, cent: %s;#it{m}_{inv} (GeV/#it{c}^{2});#it{p}_{T}^{V0} (GeV/#it{c});#it{#eta}_{V0}", GetCentBinLabel(i).Data()), iNDimIncl, binsLIncl, xminLIncl, xmaxLIncl);
    fOutputListStd->Add(fhnV0InRndALambda[i]);
    fhnV0InMedALambda[i] = new THnSparseD(Form("fhnV0InMedALambda_%d", i), Form("ALambda: Mass vs Pt in med.-cl. cones, cent: %s;#it{m}_{inv} (GeV/#it{c}^{2});#it{p}_{T}^{V0} (GeV/#it{c});#it{#eta}_{V0}", GetCentBinLabel(i).Data()), iNDimIncl, binsLIncl, xminLIncl, xmaxLIncl);
    fOutputListStd->Add(fhnV0InMedALambda[i]);
    fhnV0OutJetALambda[i] = new THnSparseD(Form("fhnV0OutJetALambda_%d", i), Form("ALambda: Pt outside jets, cent: %s;#it{m}_{inv} (GeV/#it{c}^{2});#it{p}_{T}^{V0} (GeV/#it{c});#it{#eta}_{V0}", GetCentBinLabel(i).Data()), iNDimIncl, binsLIncl, xminLIncl, xmaxLIncl);
    fOutputListStd->Add(fhnV0OutJetALambda[i]);
    fhnV0NoJetALambda[i] = new THnSparseD(Form("fhnV0NoJetALambda_%d", i), Form("ALambda: Pt in jet-less events, cent: %s;#it{m}_{inv} (GeV/#it{c}^{2});#it{p}_{T}^{V0} (GeV/#it{c});#it{#eta}_{V0}", GetCentBinLabel(i).Data()), iNDimIncl, binsLIncl, xminLIncl, xmaxLIncl);
    fOutputListStd->Add(fhnV0NoJetALambda[i]);
    // correlations
    fhnV0CorrelSEK0s[i] = new THnSparseD(Form("fhnV0CorrelSEK0s_%d", i), Form("K0s: Correlations with jets in same events, cent: %s;#it{m}_{inv} (GeV/#it{c}^{2});#it{p}_{T}^{V0} (GeV/#it{c});#it{#eta}_{V0};#it{p}_{T}^{jet} (GeV/#it{c});#Delta#it{#phi}_{V0-jet};#Delta#it{#eta}_{V0-jet}", GetCentBinLabel(i).Data()), iNDimCorrel, binsKCorrel, xminKCorrel, xmaxKCorrel);
    fOutputListStd->Add(fhnV0CorrelSEK0s[i]);
    fhnV0CorrelMEK0s[i] = new THnSparseD(Form("fhnV0CorrelMEK0s_%d", i), Form("K0s: Correlations with jets in mixed events, cent: %s;#it{m}_{inv} (GeV/#it{c}^{2});#it{p}_{T}^{V0} (GeV/#it{c});#it{#eta}_{V0};#it{p}_{T}^{jet} (GeV/#it{c});#Delta#it{#phi}_{V0-jet};#Delta#it{#eta}_{V0-jet}", GetCentBinLabel(i).Data()), iNDimCorrel, binsKCorrel, xminKCorrel, xmaxKCorrel);
    fOutputListStd->Add(fhnV0CorrelMEK0s[i]);
    fhnV0CorrelSELambda[i] = new THnSparseD(Form("fhnV0CorrelSELambda_%d", i), Form("Lambda: Correlations with jets in same events, cent: %s;#it{m}_{inv} (GeV/#it{c}^{2});#it{p}_{T}^{V0} (GeV/#it{c});#it{#eta}_{V0};#it{p}_{T}^{jet} (GeV/#it{c});#Delta#it{#phi}_{V0-jet};#Delta#it{#eta}_{V0-jet}", GetCentBinLabel(i).Data()), iNDimCorrel, binsLCorrel, xminLCorrel, xmaxLCorrel);
    fOutputListStd->Add(fhnV0CorrelSELambda[i]);
    fhnV0CorrelMELambda[i] = new THnSparseD(Form("fhnV0CorrelMELambda_%d", i), Form("Lambda: Correlations with jets in mixed events, cent: %s;#it{m}_{inv} (GeV/#it{c}^{2});#it{p}_{T}^{V0} (GeV/#it{c});#it{#eta}_{V0};#it{p}_{T}^{jet} (GeV/#it{c});#Delta#it{#phi}_{V0-jet};#Delta#it{#eta}_{V0-jet}", GetCentBinLabel(i).Data()), iNDimCorrel, binsLCorrel, xminLCorrel, xmaxLCorrel);
    fOutputListStd->Add(fhnV0CorrelMELambda[i]);

    fh2V0PtJetAngleK0s[i] = new TH2D(Form("fh2V0PtJetAngleK0s_%d", i), Form("K0s: #it{p}_{T}^{jet} vs angle V0-jet, cent: %s;#it{p}_{T}^{jet};#it{#alpha}", GetCentBinLabel(i).Data()), 2 * iNJetPtBins, dJetPtMin, dJetPtMax, 100, 0, fdDistanceV0JetMax + 0.1);
    fOutputListStd->Add(fh2V0PtJetAngleK0s[i]);
    fh2V0PtJetAngleLambda[i] = new TH2D(Form("fh2V0PtJetAngleLambda_%d", i), Form("Lambda: #it{p}_{T}^{jet} vs angle V0-jet, cent: %s;#it{p}_{T}^{jet};#it{#alpha}", GetCentBinLabel(i).Data()), 2 * iNJetPtBins, dJetPtMin, dJetPtMax, 100, 0, fdDistanceV0JetMax + 0.1);
    fOutputListStd->Add(fh2V0PtJetAngleLambda[i]);
    fh2V0PtJetAngleALambda[i] = new TH2D(Form("fh2V0PtJetAngleALambda_%d", i), Form("ALambda: #it{p}_{T}^{jet} vs angle V0-jet, cent: %s;#it{p}_{T}^{jet};#it{#alpha}", GetCentBinLabel(i).Data()), 2 * iNJetPtBins, dJetPtMin, dJetPtMax, 100, 0, fdDistanceV0JetMax + 0.1);
    fOutputListStd->Add(fh2V0PtJetAngleALambda[i]);

    // jet histograms
    fh1PtJet[i] = new TH1D(Form("fh1PtJet_%d", i), Form("Jet pt spectrum, cent: %s;#it{p}_{T} jet (GeV/#it{c})", GetCentBinLabel(i).Data()), 2 * iNJetPtBins, dJetPtMin, dJetPtMax);
    fOutputListStd->Add(fh1PtJet[i]);
    fh1EtaJet[i] = new TH1D(Form("fh1EtaJet_%d", i), Form("Jet eta spectrum, cent: %s;#it{#eta} jet", GetCentBinLabel(i).Data()), 80, -1., 1.);
    fOutputListStd->Add(fh1EtaJet[i]);
    fh2EtaPtJet[i] = new TH2D(Form("fh2EtaPtJet_%d", i), Form("Jet eta vs pT spectrum, cent: %s;#it{#eta} jet;#it{p}_{T} jet (GeV/#it{c})", GetCentBinLabel(i).Data()), 80, -1., 1., 2 * iNJetPtBins, dJetPtMin, dJetPtMax);
    fOutputListStd->Add(fh2EtaPtJet[i]);
    fh1PhiJet[i] = new TH1D(Form("fh1PhiJet_%d", i), Form("Jet phi spectrum, cent: %s;#it{#phi} jet", GetCentBinLabel(i).Data()), 90, 0., TMath::TwoPi());
    fOutputListStd->Add(fh1PhiJet[i]);
    fh2PtJetPtTrackLeading[i] = new TH2D(Form("fh2PtJetPtTrackLeading_%d", i), Form("jet pt vs leading track pt, cent: %s;#it{p}_{T}^{jet} (GeV/#it{c});#it{p}_{T} leading track (GeV/#it{c})", GetCentBinLabel(i).Data()), 4 * iNJetPtBins, dJetPtMin, dJetPtMax, 200, 0., 20);
    fOutputListStd->Add(fh2PtJetPtTrackLeading[i]);
    fh2PtJetPtTrigger[i] = new TH2D(Form("fh2PtJetPtTrigger_%d", i), Form("jet pt vs trigger track pt, cent: %s;#it{p}_{T}^{jet} (GeV/#it{c});#it{p}_{T} trigger track (GeV/#it{c})", GetCentBinLabel(i).Data()), 4 * iNJetPtBins, dJetPtMin, dJetPtMax, 200, 0., 20);
    fOutputListStd->Add(fh2PtJetPtTrigger[i]);
    fh1PtTrigger[i] = new TH1D(Form("fh1PtTrigger_%d", i), Form("trigger track pt, cent: %s;#it{p}_{T} trigger track (GeV/#it{c})", GetCentBinLabel(i).Data()), 200, 0., 20);
    fOutputListStd->Add(fh1PtTrigger[i]);
    fh1NJetPerEvent[i] = new TH1D(Form("fh1NJetPerEvent_%d", i), Form("Number of selected jets per event, cent: %s;# jets;# events", GetCentBinLabel(i).Data()), 100, 0., 100.);
    fOutputListStd->Add(fh1NJetPerEvent[i]);
    fh2EtaPhiRndCone[i] = new TH2D(Form("fh2EtaPhiRndCone_%d", i), Form("Rnd. cones: eta vs phi, cent: %s;#it{#eta} cone;#it{#phi} cone", GetCentBinLabel(i).Data()), 80, -1., 1., 90, 0., TMath::TwoPi());
    fOutputListStd->Add(fh2EtaPhiRndCone[i]);
    fh2EtaPhiMedCone[i] = new TH2D(Form("fh2EtaPhiMedCone_%d", i), Form("Med.-cl. cones: eta vs phi, cent: %s;#it{#eta} cone;#it{#phi} cone", GetCentBinLabel(i).Data()), 80, -1., 1., 90, 0., TMath::TwoPi());
    fOutputListStd->Add(fh2EtaPhiMedCone[i]);
    fh1DistanceJets[i] = new TH1D(Form("fh1DistanceJets_%d", i), Form("Distance between jets in #eta-#phi, cent: %s;#it{D}", GetCentBinLabel(i).Data()), 40, 0., 4.);
    fOutputListStd->Add(fh1DistanceJets[i]);
    fh1DistanceV0JetsK0s[i] = new TH1D(Form("fh1DistanceV0JetsK0s_%d", i), Form("Distance between V0 and the closest jet in #eta-#phi, cent: %s;#it{D}", GetCentBinLabel(i).Data()), 80, 0., 2.);
    fh1DistanceV0JetsLambda[i] = new TH1D(Form("fh1DistanceV0JetsLambda_%d", i), Form("Distance between V0 and the closest jet in #eta-#phi, cent: %s;#it{D}", GetCentBinLabel(i).Data()), 80, 0., 2.);
    fh1DistanceV0JetsALambda[i] = new TH1D(Form("fh1DistanceV0JetsALambda_%d", i), Form("Distance between V0 and the closest jet in #eta-#phi, cent: %s;#it{D}", GetCentBinLabel(i).Data()), 80, 0., 2.);
    fOutputListStd->Add(fh1DistanceV0JetsK0s[i]);
    fOutputListStd->Add(fh1DistanceV0JetsLambda[i]);
    fOutputListStd->Add(fh1DistanceV0JetsALambda[i]);
    // event histograms
    fh1VtxZ[i] = new TH1D(Form("fh1VtxZ_%d", i), Form("#it{z} coordinate of the primary vertex, cent: %s;#it{z} (cm)", GetCentBinLabel(i).Data()), 150, -1.5 * fdCutVertexZ, 1.5 * fdCutVertexZ);
    fOutputListQA->Add(fh1VtxZ[i]);
    fh1VtxZME[i] = new TH1D(Form("fh1VtxZME_%d", i), Form("#it{z} coordinate of the primary vertex, cent: %s;#it{z} (cm)", GetCentBinLabel(i).Data()), 150, -1.5 * fdCutVertexZ, 1.5 * fdCutVertexZ);
    fOutputListQA->Add(fh1VtxZME[i]);
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
      fh2V0K0sInJetPtMCGen[i] = new TH2D(Form("fh2V0K0sInJetPtMCGen_%d", i), Form("MC K0s in jet generated: pt-ptJet spectrum, cent: %s;MC #it{p}_{T} (GeV/#it{c});#it{p}_{T}^{jet} (GeV/#it{c})", GetCentBinLabel(i).Data()), iNBinsPtV0InJet, dPtV0Min, dPtV0Max, iNJetPtBins, dJetPtMin, dJetPtMax);
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
      fh2V0LambdaInJetPtMCGen[i] = new TH2D(Form("fh2V0LambdaInJetPtMCGen_%d", i), Form("MC Lambda in jet generated: pt-ptJet spectrum, cent: %s;MC #it{p}_{T} (GeV/#it{c});#it{p}_{T}^{jet} (GeV/#it{c})", GetCentBinLabel(i).Data()), iNBinsPtV0InJet, dPtV0Min, dPtV0Max, iNJetPtBins, dJetPtMin, dJetPtMax);
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
      fh2V0ALambdaInJetPtMCGen[i] = new TH2D(Form("fh2V0ALambdaInJetPtMCGen_%d", i), Form("MC ALambda in jet generated: pt-ptJet spectrum, cent: %s;MC #it{p}_{T} (GeV/#it{c});#it{p}_{T}^{jet} (GeV/#it{c})", GetCentBinLabel(i).Data()), iNBinsPtV0InJet, dPtV0Min, dPtV0Max, iNJetPtBins, dJetPtMin, dJetPtMax);
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
      fhnV0K0sInclDaughterEtaPtPtMCRec[i] = new THnSparseD(Form("fhnV0K0sInclDaughterEtaPtPtMCRec_%d", i), Form("MC K0S, inclusive, assoc., daughters: charge-etaD-ptD-etaV0-ptV0-ptJet, cent: %s;charge;eta gen daughter;pT gen daughter;eta gen V0;pT gen V0;pT rec jet", GetCentBinLabel(i).Data()), iNDimEtaD, binsEtaDaughter, xminEtaDaughter, xmaxEtaDaughter);
      fhnV0K0sInJetsDaughterEtaPtPtMCRec[i] = new THnSparseD(Form("fhnV0K0sInJetsDaughterEtaPtPtMCRec_%d", i), Form("MC K0S, in JC, assoc., daughters: charge-etaD-ptD-etaV0-ptV0-ptJet, cent: %s;charge;eta gen daughter;pT gen daughter;eta gen V0;pT gen V0;pT rec jet", GetCentBinLabel(i).Data()), iNDimEtaD, binsEtaDaughter, xminEtaDaughter, xmaxEtaDaughter);
      fhnV0LambdaInclDaughterEtaPtPtMCRec[i] = new THnSparseD(Form("fhnV0LambdaInclDaughterEtaPtPtMCRec_%d", i), Form("MC Lambda, inclusive, assoc., daughters: charge-etaD-ptD-etaV0-ptV0-ptJet, cent: %s;charge;eta gen daughter;pT gen daughter;eta gen V0;pT gen V0;pT rec jet", GetCentBinLabel(i).Data()), iNDimEtaD, binsEtaDaughter, xminEtaDaughter, xmaxEtaDaughter);
      fhnV0LambdaInJetsDaughterEtaPtPtMCRec[i] = new THnSparseD(Form("fhnV0LambdaInJetsDaughterEtaPtPtMCRec_%d", i), Form("MC Lambda, in JC, assoc., daughters: charge-etaD-ptD-etaV0-ptV0-ptJet, cent: %s;charge;eta gen daughter;pT gen daughter;eta gen V0;pT gen V0;pT rec jet", GetCentBinLabel(i).Data()), iNDimEtaD, binsEtaDaughter, xminEtaDaughter, xmaxEtaDaughter);
      fhnV0ALambdaInclDaughterEtaPtPtMCRec[i] = new THnSparseD(Form("fhnV0ALambdaInclDaughterEtaPtPtMCRec_%d", i), Form("MC ALambda, inclusive, assoc., daughters: charge-etaD-ptD-etaV0-ptV0-ptJet, cent: %s;charge;eta gen daughter;pT gen daughter;eta gen V0;pT gen V0;pT rec jet", GetCentBinLabel(i).Data()), iNDimEtaD, binsEtaDaughter, xminEtaDaughter, xmaxEtaDaughter);
      fhnV0ALambdaInJetsDaughterEtaPtPtMCRec[i] = new THnSparseD(Form("fhnV0ALambdaInJetsDaughterEtaPtPtMCRec_%d", i), Form("MC ALambda, in JC, assoc., daughters: charge-etaD-ptD-etaV0-ptV0-ptJet, cent: %s;charge;eta gen daughter;pT gen daughter;eta gen V0;pT gen V0;pT rec jet", GetCentBinLabel(i).Data()), iNDimEtaD, binsEtaDaughter, xminEtaDaughter, xmaxEtaDaughter);

      fOutputListMC->Add(fhnV0K0sInclDaughterEtaPtPtMCRec[i]);
      fOutputListMC->Add(fhnV0K0sInJetsDaughterEtaPtPtMCRec[i]);
      fOutputListMC->Add(fhnV0LambdaInclDaughterEtaPtPtMCRec[i]);
      fOutputListMC->Add(fhnV0LambdaInJetsDaughterEtaPtPtMCRec[i]);
      fOutputListMC->Add(fhnV0ALambdaInclDaughterEtaPtPtMCRec[i]);
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
    fh2QAV0PtNCls[i] = new TH2D(Form("fh2QAV0PtNCls_%d", i), "QA: Daughter Pt vs TPC clusters;pt;# TPC clusters", 100, 0, 10, 160, 0, 160);
    fh2QAV0PtChi[i] = new TH2D(Form("fh2QAV0PtChi_%d", i), "QA: Daughter Pt vs Chi2/ndf;pt;Chi2/ndf", 100, 0, 10, 100, 0, 100);
    fh1QAV0TPCRowsFind[i] = new TH1D(Form("fh1QAV0TPCRowsFind_%d", i), "QA: TPC Rows/Findable", 100, 0, 2);
    fh1QAV0Eta[i] = new TH1D(Form("fh1QAV0Eta_%d", i), "QA: Daughter Eta", 200, -2, 2);
    fh2QAV0EtaRows[i] = new TH2D(Form("fh2QAV0EtaRows_%d", i), "QA: Daughter Eta vs TPC rows;#eta;TPC rows", 200, -2, 2, 160, 0, 160);
    fh2QAV0PtRows[i] = new TH2D(Form("fh2QAV0PtRows_%d", i), "QA: Daughter Pt vs TPC rows;pt;TPC rows", 100, 0, 10, 160, 0, 160);
    fh2QAV0PhiRows[i] = new TH2D(Form("fh2QAV0PhiRows_%d", i), "QA: Daughter Phi vs TPC rows;#phi;TPC rows", 90, 0, TMath::TwoPi(), 160, 0, 160);
    fh2QAV0NClRows[i] = new TH2D(Form("fh2QAV0NClRows_%d", i), "QA: Daughter NCl vs TPC rows;findable clusters;TPC rows", 100, 0, 160, 160, 0, 160);
    fh2QAV0EtaNCl[i] = new TH2D(Form("fh2QAV0EtaNCl_%d", i), "QA: Daughter Eta vs NCl;#eta;findable clusters", 200, -2, 2, 160, 0, 160);

    fh2QAV0EtaPtK0sPeak[i] = new TH2D(Form("fh2QAV0EtaPtK0sPeak_%d", i), "QA: K0s: Daughter Eta vs V0 pt, peak;track eta;V0 pt", 200, -2, 2, iNBinsPtV0, dPtV0Min, dPtV0Max);
    fh2QAV0EtaEtaK0s[i] = new TH2D(Form("fh2QAV0EtaEtaK0s_%d", i), "QA: K0s: Eta vs Eta Daughter", 200, -2, 2, 200, -2, 2);
    fh2QAV0PhiPhiK0s[i] = new TH2D(Form("fh2QAV0PhiPhiK0s_%d", i), "QA: K0s: Phi vs Phi Daughter", 90, 0, TMath::TwoPi(), 90, 0, TMath::TwoPi());
    fh1QAV0RapK0s[i] = new TH1D(Form("fh1QAV0RapK0s_%d", i), "QA: K0s: V0 Rapidity", 200, -2, 2);
    fh2QAV0PtPtK0sPeak[i] = new TH2D(Form("fh2QAV0PtPtK0sPeak_%d", i), "QA: K0s: Daughter Pt vs Pt;neg pt;pos pt", 100, 0, 5, 100, 0, 5);

    fh2QAV0EtaPtLambdaPeak[i] = new TH2D(Form("fh2QAV0EtaPtLambdaPeak_%d", i), "QA: Lambda: Daughter Eta vs V0 pt, peak;track eta;V0 pt", 200, -2, 2, iNBinsPtV0, dPtV0Min, dPtV0Max);
    fh2QAV0EtaEtaLambda[i] = new TH2D(Form("fh2QAV0EtaEtaLambda_%d", i), "QA: Lambda: Eta vs Eta Daughter", 200, -2, 2, 200, -2, 2);
    fh2QAV0PhiPhiLambda[i] = new TH2D(Form("fh2QAV0PhiPhiLambda_%d", i), "QA: Lambda: Phi vs Phi Daughter", 90, 0, TMath::TwoPi(), 90, 0, TMath::TwoPi());
    fh1QAV0RapLambda[i] = new TH1D(Form("fh1QAV0RapLambda_%d", i), "QA: Lambda: V0 Rapidity", 200, -2, 2);
    fh2QAV0PtPtLambdaPeak[i] = new TH2D(Form("fh2QAV0PtPtLambdaPeak_%d", i), "QA: Lambda: Daughter Pt vs Pt;neg pt;pos pt", 100, 0, 5, 100, 0, 5);

    fh2QAV0EtaPtALambdaPeak[i] = new TH2D(Form("fh2QAV0EtaPtALambdaPeak_%d", i), "QA: anti-Lambda: Daughter Eta vs V0 pt, peak;track eta;V0 pt", 200, -2, 2, iNBinsPtV0, dPtV0Min, dPtV0Max);
    fh2QAV0EtaEtaALambda[i] = new TH2D(Form("fh2QAV0EtaEtaALambda_%d", i), "QA: anti-Lambda: Eta vs Eta Daughter", 200, -2, 2, 200, -2, 2);
    fh2QAV0PhiPhiALambda[i] = new TH2D(Form("fh2QAV0PhiPhiALambda_%d", i), "QA: anti-Lambda: Phi vs Phi Daughter", 90, 0, TMath::TwoPi(), 90, 0, TMath::TwoPi());
    fh1QAV0RapALambda[i] = new TH1D(Form("fh1QAV0RapALambda_%d", i), "QA: anti-Lambda: V0 Rapidity", 200, -2, 2);
    fh2QAV0PtPtALambdaPeak[i] = new TH2D(Form("fh2QAV0PtPtALambdaPeak_%d", i), "QA: anti-Lambda: Daughter Pt vs Pt;neg pt;pos pt", 100, 0, 5, 100, 0, 5);

    fh2QAV0PhiPtK0sPeak[i] = new TH2D(Form("fh2QAV0PhiPtK0sPeak_%d", i), "QA: K0S: #phi-pt;#phi;pt", 90, 0, TMath::TwoPi(), iNBinsPtV0, dPtV0Min, dPtV0Max);
    fh2QAV0PhiPtLambdaPeak[i] = new TH2D(Form("fh2QAV0PhiPtLambdaPeak_%d", i), "QA: Lambda: #phi-pt;#phi;pt", 90, 0, TMath::TwoPi(), iNBinsPtV0, dPtV0Min, dPtV0Max);
    fh2QAV0PhiPtALambdaPeak[i] = new TH2D(Form("fh2QAV0PhiPtALambdaPeak_%d", i), "QA: anti-Lambda: #phi-pt;#phi;pt", 90, 0, TMath::TwoPi(), iNBinsPtV0, dPtV0Min, dPtV0Max);
    fh1QAV0Pt[i] = new TH1D(Form("fh1QAV0Pt_%d", i), "QA: Daughter Pt", 100, 0, 5);
    fh1QAV0Charge[i] = new TH1D(Form("fh1QAV0Charge_%d", i), "QA: V0 Charge", 3, -1, 2);
    fh1QAV0DCAVtx[i] = new TH1D(Form("fh1QAV0DCAVtx_%d", i), "QA: DCA daughters to primary vertex", 1000, 0, 10);
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
    fOutputListQA->Add(fh2QAV0PtNCls[i]);
    fOutputListQA->Add(fh2QAV0PtChi[i]);
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

    fOutputListQA->Add(fh2QAV0EtaPtALambdaPeak[i]);
    fOutputListQA->Add(fh2QAV0EtaEtaALambda[i]);
    fOutputListQA->Add(fh2QAV0PhiPhiALambda[i]);
    fOutputListQA->Add(fh1QAV0RapALambda[i]);
    fOutputListQA->Add(fh2QAV0PtPtALambdaPeak[i]);

    fOutputListQA->Add(fh2QAV0PhiPtK0sPeak[i]);
    fOutputListQA->Add(fh2QAV0PhiPtLambdaPeak[i]);
    fOutputListQA->Add(fh2QAV0PhiPtALambdaPeak[i]);
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
}

void AliAnalysisTaskV0sInJetsEmcal::ExecOnce()
{
  AliAnalysisTaskEmcalJet::ExecOnce();

  if(fbJetSelection)
  {
    if(fJetsCont && fJetsCont->GetArray() == 0)
      fJetsCont = 0;
    if(fJetsBgCont && fJetsBgCont->GetArray() == 0)
      fJetsBgCont = 0;
    if(fbCompareTriggers)
    {
      if(fTracksCont && fTracksCont->GetArray() == 0)
        fTracksCont = 0;
    }
  }

  // Jet selection
  Double_t dCutEtaJetMax = fdCutEtaV0Max - fdDistanceV0JetMax; // max jet |pseudorapidity|, to make sure that V0s can appear in the entire jet area
  if(fbJetSelection && fJetsCont)
  {
    if(fdCutPtTrackJetMin > 0.) fJetsCont->SetPtBiasJetTrack(fdCutPtTrackJetMin);
    if(fdCutAreaPercJetMin > 0.) fJetsCont->SetPercAreaCut(fdCutAreaPercJetMin);
    if(fdCutEtaV0Max > 0. && fdDistanceV0JetMax > 0.) fJetsCont->SetJetEtaLimits(-dCutEtaJetMax, dCutEtaJetMax);
    // Trigger selection
    if(fbCompareTriggers && fTracksCont)
    {
      fTracksCont->SetFilterHybridTracks(kTRUE);
      fTracksCont->SetParticlePtCut(5.);
      fTracksCont->SetParticleEtaLimits(-0.8, 0.8);
    }
  }

  printf("=======================================================\n");
  printf("%s: Configuration summary:\n", ClassName());
  printf("task name: %s\n", GetName());
  printf("-------------------------------------------------------\n");
  printf("collision system: %s\n", fbIsPbPb ? "Pb+Pb" : "p+p");
  printf("data type: %s\n", fbMCAnalysis ? "MC" : "real");
  if(fbMCAnalysis)
    printf("MC generator: %s\n", fsGeneratorName.Length() ? fsGeneratorName.Data() : "any");
  if(fbIsPbPb)
    printf("centrality range: %g-%g %%\n", fdCutCentLow, fdCutCentHigh);
  if(fdCutVertexZ > 0.) printf("max |z| of the prim vtx [cm]: %g\n", fdCutVertexZ);
  if(fdCutVertexR2 > 0.) printf("max r^2 of the prim vtx [cm^2]: %g\n", fdCutVertexR2);
  if(fdCutDeltaZMax > 0.) printf("max |Delta z| between nominal prim vtx and SPD vtx [cm]: %g\n", fdCutDeltaZMax);
  printf("-------------------------------------------------------\n");
  if(fbTPCRefit) printf("TPC refit for daughter tracks\n");
  if(fbRejectKinks) printf("reject kink-like production vertices of daughter tracks\n");
  if(fbFindableClusters) printf("require positive number of findable clusters\n");
  if(fdCutNCrossedRowsTPCMin > 0.) printf("min number of crossed TPC rows: %g\n", fdCutNCrossedRowsTPCMin);
  if(fdCutCrossedRowsOverFindMin > 0.) printf("min ratio crossed rows / findable clusters: %g\n", fdCutCrossedRowsOverFindMin);
  if(fdCutCrossedRowsOverFindMax > 0.) printf("max ratio crossed rows / findable clusters: %g\n", fdCutCrossedRowsOverFindMax);
  if(fdCutPtDaughterMin > 0.) printf("min pt of daughter tracks [GeV/c]: %g\n", fdCutPtDaughterMin);
  if(fdCutDCAToPrimVtxMin > 0.) printf("min DCA of daughters to the prim vtx [cm]: %g\n", fdCutDCAToPrimVtxMin);
  if(fdCutDCADaughtersMax > 0.) printf("max DCA between daughters [sigma of TPC tracking]: %g\n", fdCutDCADaughtersMax);
  if(fdCutEtaDaughterMax > 0.) printf("max |eta| of daughter tracks: %g\n", fdCutEtaDaughterMax);
  if(fdCutNSigmadEdxMax > 0. && (!fbIsPbPb || (fbIsPbPb && fdPtProtonPIDMax > 0.))) printf("max |Delta(dE/dx)| in the TPC [sigma dE/dx]: %g\n", fdCutNSigmadEdxMax);
  if(fdCutNSigmadEdxMax > 0. && fbIsPbPb && fdPtProtonPIDMax > 0.) printf("max pt of proton for applying PID cut [GeV/c]: %g\n", fdPtProtonPIDMax);
  printf("V0 reconstruction method: %s\n", fbOnFly ? "on-the-fly" : "offline");
  if(fdCutCPAKMin > 0.) printf("min CPA, K0S: %g\n", fdCutCPAKMin);
  if(fdCutCPALMin > 0.) printf("min CPA, (A)Lambda: %g\n", fdCutCPALMin);
  if(fdCutRadiusDecayMin > 0. && fdCutRadiusDecayMax > 0.) printf("R of the decay vertex [cm]: %g-%g\n", fdCutRadiusDecayMin, fdCutRadiusDecayMax);
  if(fdCutEtaV0Max > 0.) printf("max |eta| of V0: %g\n", fdCutEtaV0Max);
  if(fdCutRapV0Max > 0.) printf("max |y| of V0: %g\n", fdCutRapV0Max);
  if(fdCutNTauKMax > 0.) printf("max proper lifetime, K0S [tau]: %g\n", fdCutNTauKMax);
  if(fdCutNTauLMax > 0.) printf("max proper lifetime, (A)Lambda [tau]: %g\n", fdCutNTauLMax);
  if(fbCutArmPod) printf("Armenteros-Podolanski cut for K0S\n");
  if(fbCutCross) printf("cross-contamination cut\n");
  printf("-------------------------------------------------------\n");
  printf("analysis of V0s in jets: %s\n", fbJetSelection ? "yes" : "no");
  if(fbJetSelection)
  {
    printf("rho subtraction: ");
    if(fiBgSubtraction == 0) printf("none\n");
    else if(fiBgSubtraction == 1) printf("scalar\n");
    else if(fiBgSubtraction == 2) printf("vector\n");
    if(fdCutPtJetMin > 0.) printf("min jet pt [GeV/c]: %g\n", fdCutPtJetMin);
    if(fdCutPtTrackJetMin > 0.) printf("min pt of leading jet-track [GeV/c]: %g\n", fdCutPtTrackJetMin);
    if(fdCutAreaPercJetMin > 0.) printf("min area of jet [pi*R^2]: %g\n", fdCutAreaPercJetMin);
    if(fdCutEtaV0Max > 0. && fdDistanceV0JetMax > 0.) printf("max |eta| of jet: %g\n", dCutEtaJetMax);
    printf("max distance between V0 and jet axis: %g\n", fdDistanceV0JetMax);
    printf("angular correlations of V0s with jets: %s\n", fbCorrelations ? "yes" : "no");
    if(fbCorrelations) printf("max |delta-eta_V0-jet|: %g\n", fdDeltaEtaMax);
    printf("pt correlations of jets with trigger tracks: %s\n", fbCompareTriggers ? "yes" : "no");
    printf("-------------------------------------------------------\n");
    if(fJetsCont)
    {
      printf("Signal jet container parameters\n");
      printf("Jet R = %g\n", fJetsCont->GetJetRadius());
      fJetsCont->PrintCuts();
    }
    else
      printf("No signal jet container!\n");
    if(fJetsBgCont)
    {
      printf("Background jet container parameters\n");
      printf("Jet R = %g\n", fJetsBgCont->GetJetRadius());
      fJetsBgCont->PrintCuts();
    }
    else
      printf("No background jet container\n");
  }
  printf("=======================================================\n");
}

Bool_t AliAnalysisTaskV0sInJetsEmcal::Run()
{
// Run analysis code here, if needed. It will be executed before FillHistograms().

  return kTRUE; // If return kFALSE FillHistogram() will NOT be executed.
}

Bool_t AliAnalysisTaskV0sInJetsEmcal::FillHistograms()
{
  // Main loop, called for each event
  if(fDebug > 0) printf("%s %s::%s: %s\n", GetName(), ClassName(), __func__, "Start");

  fh1EventCounterCut->Fill(0); // all available selected events (collision candidates)

  fAODIn = dynamic_cast<AliAODEvent*>(InputEvent()); // input AOD
  if(!fAODIn)
  {
    AliError("No input AOD found!");
    return kFALSE;
  }
  if(fDebug > 1) printf("%s %s::%s: %s\n", GetName(), ClassName(), __func__, "Loading AOD OK");

  TClonesArray* arrayMC = 0; // array particles in the MC event
  AliAODMCHeader* headerMC = 0; // MC header
  Int_t iNTracksMC = 0; // number of MC tracks
  Double_t dPrimVtxMCX = 0., dPrimVtxMCY = 0., dPrimVtxMCZ = 0.; // position of the MC primary vertex

  // Simulation info
  if(fbMCAnalysis)
  {
    fEventMC = MCEvent();
    arrayMC = (TClonesArray*)fAODIn->FindListObject(AliAODMCParticle::StdBranchName());
    if(!arrayMC || !fEventMC)
    {
      AliError("No MC array/event found!");
      return kFALSE;
    }
    if(fDebug > 1) printf("%s %s::%s: %s\n", GetName(), ClassName(), __func__, "MC array found");
    iNTracksMC = arrayMC->GetEntriesFast();
    if(fDebug > 2) printf("%s %s::%s: %s\n", GetName(), ClassName(), __func__, Form("There are %d MC tracks in this event", iNTracksMC));
    headerMC = (AliAODMCHeader*)fAODIn->FindListObject(AliAODMCHeader::StdBranchName());
    if(!headerMC)
    {
      AliError("No MC header found!");
      return kFALSE;
    }
    // get position of the MC primary vertex
    dPrimVtxMCX = headerMC->GetVtxX();
    dPrimVtxMCY = headerMC->GetVtxY();
    dPrimVtxMCZ = headerMC->GetVtxZ();
  }

  // PID Response Task object
  AliPIDResponse* fPIDResponse = 0;
  if(fdCutNSigmadEdxMax > 0. && (!fbIsPbPb || (fbIsPbPb && fdPtProtonPIDMax > 0.)))
  {
    AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
    AliInputEventHandler* inputHandler = (AliInputEventHandler*)mgr->GetInputEventHandler();
    fPIDResponse = inputHandler->GetPIDResponse();
    if(!fPIDResponse)
    {
      AliError("No PID response object found!");
      return kFALSE;
    }
  }

  // AOD files are OK
  fh1EventCounterCut->Fill(1);

  // Event selection
  if(!IsSelectedForJets(fAODIn, fdCutVertexZ, fdCutVertexR2, fdCutCentLow, fdCutCentHigh, fdCutDeltaZMax))
  {
    if(fDebug > 0) printf("%s %s::%s: %s\n", GetName(), ClassName(), __func__, "Event rejected");
    return kFALSE;
  }
  if(fDebug > 0) printf("%s %s::%s: %s\n", GetName(), ClassName(), __func__, Form("Event accepted: cent. %g", fdCentrality));
  if(!fbIsPbPb)
    fdCentrality = 0.; // select the first bin for p+p data
  Int_t iCentIndex = GetCentralityBinIndex(fdCentrality); // get index of centrality bin
  if(iCentIndex < 0)
  {
    AliError(Form("Event is out of histogram range: cent: %g!", fdCentrality));
    return kFALSE;
  }
  fh1EventCounterCut->Fill(8); // selected events (centrality OK)
  fh1EventCounterCutCent[iCentIndex]->Fill(8);

  UInt_t iNTracks = fAODIn->GetNumberOfTracks(); // get number of tracks in event
  if(fDebug > 2) printf("%s %s::%s: %s\n", GetName(), ClassName(), __func__, Form("There are %d tracks in this event", iNTracks));

//  Double_t dMagField = fAODIn->GetMagneticField();
//  if(fDebug > 2) printf("%s %s::%s: %s\n", GetName(), ClassName(), __func__, Form("Magnetic field: %g", dMagField));

  Int_t iNV0s = fAODIn->GetNumberOfV0s(); // get the number of V0 candidates
  if(!iNV0s)
  {
    if(fDebug > 0) printf("%s %s::%s: %s\n", GetName(), ClassName(), __func__, "No V0s found in event");
  }

  //===== Event is OK for the analysis =====
  fh1EventCent->Fill(iCentIndex);
  fh1EventCent2->Fill(fdCentrality);
  fh2EventCentTracks->Fill(fdCentrality, iNTracks);

  AliEventPool* pool = 0;
  Bool_t bPoolReady = kFALSE; // status of pool
  Double_t dZVtxME = -100; // z_vtx for events used in mixed events
  if(fbCorrelations)
  {
    AliAODVertex* vertex = fAODIn->GetPrimaryVertex();
    dZVtxME = vertex->GetZ();
    pool = fPoolMgr->GetEventPool(fdCentrality, dZVtxME);
    if(!pool)
    {
      AliError(Form("No pool found for centrality = %g, z_vtx = %g!", fdCentrality, dZVtxME));
      return kFALSE;
    }
    bPoolReady = pool->IsReady();
    if(fDebug > 2)
    {
      //      pool->SetDebug(1);
      //      pool->PrintInfo();
      printf("%s %s::%s: %s\n", GetName(), ClassName(), __func__, Form("Pool %d-%d: events in pool = %d, jets in pool = %d", pool->MultBinIndex(), pool->ZvtxBinIndex(), pool->GetCurrentNEvents(), pool->NTracksInPool()));
    }
  }

  if(iNV0s)
  {
    fh1EventCounterCut->Fill(9); // events with V0s
    fh1EventCounterCutCent[iCentIndex]->Fill(9);
  }

  AliAODv0* v0 = 0; // pointer to V0 candidates
  TVector3 vecV0Momentum; // 3D vector of V0 momentum
  Double_t dMassV0K0s = 0; // invariant mass of the K0s candidate
  Double_t dMassV0Lambda = 0; // invariant mass of the Lambda candidate
  Double_t dMassV0ALambda = 0; // invariant mass of the Lambda candidate
  Int_t iNV0CandTot = 0; // counter of all V0 candidates at the beginning
  Int_t iNV0CandK0s = 0; // counter of K0s candidates at the end
  Int_t iNV0CandLambda = 0; // counter of Lambda candidates at the end
  Int_t iNV0CandALambda = 0; // counter of Lambda candidates at the end

  Bool_t bPrintCuts = 0; // print out which cuts are applied
  Bool_t bPrintJetSelection = 0; // print out which jets are selected

  // Other cuts
  Double_t dNSigmaMassMax = 3.; // [sigma m] max difference between candidate mass and real particle mass (used only for mass peak method of signal extraction)
  Double_t dDistPrimaryMax = 0.01; // [cm] max distance of production point to the primary vertex (criterion for choice of MC particles considered as primary)

  // Mean lifetime
  Double_t dCTauK0s = 2.6844; // [cm] c*tau of K0S
  Double_t dCTauLambda = 7.89; // [cm] c*tau of Lambda

  // particle masses from PDG
  Double_t dMassPDGK0s = TDatabasePDG::Instance()->GetParticle(kK0Short)->Mass();
  Double_t dMassPDGLambda = TDatabasePDG::Instance()->GetParticle(kLambda0)->Mass();

  // PDG codes of used particles
  Int_t iPdgCodePion = 211;
  Int_t iPdgCodeProton = 2212;
  Int_t iPdgCodeK0s = 310;
  Int_t iPdgCodeLambda = 3122;

  Double_t dCutEtaJetMax = fdCutEtaV0Max - fdDistanceV0JetMax; // max jet |pseudorapidity|, to make sure that V0s can appear in the entire jet area
  Double_t dRadiusExcludeCone = 2 * fdDistanceV0JetMax; // radius of cones around jets excluded for V0 outside jets
  Bool_t bLeadingJetOnly = 0; // consider only leading jets

  Int_t iNJet = 0; // number of reconstructed jets in the input
  TClonesArray* arrayJetSel = new TClonesArray("AliAODJet", 0); // object where the selected jets are copied as AliAODJet
  Int_t iNJetSel = 0; // number of selected reconstructed jets
  TClonesArray* arrayJetPerp = new TClonesArray("AliAODJet", 0); // object where the perp. cones are stored
  Int_t iNJetPerp = 0; // number of perpendicular cones
  TObjArray* arrayMixedEventAdd = new TObjArray; // array of jets from this event to be added to the pool
  arrayMixedEventAdd->SetOwner(kTRUE);

  AliAODJet* jet = 0; // pointer to a jet
  AliAODJet* jetPerp = 0; // pointer to a perp. cone
  AliAODJet* jetRnd = 0; // pointer to a rand. cone
  AliEmcalJet* jetMed = 0; // pointer to a median cluster
  TVector3 vecJetMomentum; // 3D vector of jet momentum
  TVector3 vecJetMomentumPair; // 3D vector of another jet momentum for calculating distance between jets
  Bool_t bJetEventGood = kTRUE; // indicator of good jet events
  Double_t dRho = 0; // average bg pt density
  TLorentzVector vecJetSel; // 4-momentum of selected jet
  TLorentzVector vecPerpPlus; // 4-momentum of perpendicular cone plus
  TLorentzVector vecPerpMinus; // 4-momentum of perpendicular cone minus
  Double_t dPtJet = 0; // pt of jet associated with V0
  Double_t dPtJetTrackLeading = 0; // pt of leading track of jet associated with V0

  if(fbJetSelection) // analysis of V0s in jets is switched on
  {
    if(!fJetsCont)
    {
      AliError("No signal jet container!");
      bJetEventGood = kFALSE;
    }
    if(bJetEventGood)
      iNJet = fJetsCont->GetNJets();
    if(bJetEventGood && !iNJet) // check whether there are some jets
    {
      if(fDebug > 0) printf("%s %s::%s: %s\n", GetName(), ClassName(), __func__, "No jets in array");
      bJetEventGood = kFALSE;
    }
    if(fbIsPbPb && bJetEventGood && !fJetsBgCont)
    {
      AliWarning("No bg jet container!");
//      bJetEventGood = kFALSE;
    }
  }
  else // no in-jet analysis
    bJetEventGood = kFALSE;

  // select good jets and copy them to another array
  if(bJetEventGood)
  {
    if(fiBgSubtraction)
    {
      dRho = fJetsCont->GetRhoVal();
      if(fDebug > 4) printf("%s::%s: %s\n", ClassName(), __func__, Form("Loaded rho value: %g", dRho));
    }
    if(bLeadingJetOnly)
      iNJet = 1; // only leading jets
    if(fDebug > 2) printf("%s %s::%s: %s\n", GetName(), ClassName(), __func__, Form("Jet selection for %d jets", iNJet));
    for(Int_t iJet = 0; iJet < iNJet; iJet++)
    {
      AliEmcalJet* jetSel = (AliEmcalJet*)(fJetsCont->GetAcceptJet(iJet)); // load a jet in the list
      if(bLeadingJetOnly)
        jetSel = fJetsCont->GetLeadingJet();
      if(!jetSel)
      {
        if(fDebug > 4) printf("%s %s::%s: %s\n", GetName(), ClassName(), __func__, Form("Jet %d not accepted in container", iJet));
        continue;
      }

      Double_t dPtJetCorr = jetSel->Pt(); // raw pt
      Double_t dEtaJetCorr = jetSel->Eta(); // raw eta
      Double_t dPhiJetCorr = jetSel->Phi(); // raw phi
      Double_t dPtTrackJet = fJetsCont->GetLeadingHadronPt(jetSel); // pt of leading track

      if(fiBgSubtraction == 1) // scalar subtraction: correct pt only
      {
        dPtJetCorr = jetSel->PtSub(dRho);
      }
      else if(fiBgSubtraction == 2) // vector subtraction: correct momentum vector
      {
        vecJetSel = jetSel->SubtractRhoVect(dRho);
        dPtJetCorr = vecJetSel.Pt();
        dEtaJetCorr = vecJetSel.Eta();
        dPhiJetCorr = TVector2::Phi_0_2pi(vecJetSel.Phi());
      }

      if(bPrintJetSelection)
        if(fDebug > 4) printf("jet: i = %d, pT = %g, eta = %g, phi = %g, pt lead tr = %g ", iJet, dPtJetCorr, dEtaJetCorr, dPhiJetCorr, dPtTrackJet);
      if(fdCutPtJetMin > 0. && dPtJetCorr < fdCutPtJetMin) // selection of high-pt jets, needs to be applied on the pt after bg subtraction
      {
        if(bPrintJetSelection)
          if(fDebug > 4) printf("rejected (pt)\n");
        continue;
      }
      if(bPrintJetSelection)
        if(fDebug > 4) printf("accepted\n");
      if(fDebug > 4) printf("%s %s::%s: %s\n", GetName(), ClassName(), __func__, Form("Jet %d with pt %g passed selection", iJet, dPtJetCorr));

      vecJetSel.SetPtEtaPhiM(dPtJetCorr, dEtaJetCorr, dPhiJetCorr, 0.);
      vecPerpPlus = vecJetSel;
      vecPerpMinus = vecJetSel;
      vecPerpPlus.RotateZ(TMath::Pi() / 2.); // rotate vector by +90 deg around z
      vecPerpMinus.RotateZ(-TMath::Pi() / 2.); // rotate vector by -90 deg around z
      if(fDebug > 4) printf("%s %s::%s: %s\n", GetName(), ClassName(), __func__, Form("Adding perp. cones number %d, %d, pt %g", iNJetPerp, iNJetPerp + 1, vecPerpPlus.Pt()));
      new((*arrayJetPerp)[iNJetPerp++]) AliAODJet(vecPerpPlus); // write perp. cone to the array
      new((*arrayJetPerp)[iNJetPerp++]) AliAODJet(vecPerpMinus); // write perp. cone to the array
      if(fDebug > 4) printf("%s %s::%s: %s\n", GetName(), ClassName(), __func__, Form("Adding jet number %d", iNJetSel));
      new((*arrayJetSel)[iNJetSel++]) AliAODJet(vecJetSel); // copy selected jet to the array
      ((AliAODJet*)arrayJetSel->At(iNJetSel - 1))->SetPtLeading(dPtTrackJet);
      fh2PtJetPtTrackLeading[iCentIndex]->Fill(dPtJetCorr, dPtTrackJet); // pt_jet vs pt of leading jet track
      if(fbCorrelations)
        arrayMixedEventAdd->Add(new TLorentzVector(vecJetSel)); // copy selected jet to the list of new jets for event mixing
    }
    if(fDebug > 3) printf("%s %s::%s: %s\n", GetName(), ClassName(), __func__, Form("Added jets: %d", iNJetSel));
    iNJetSel = arrayJetSel->GetEntriesFast();
    if(fDebug > 2) printf("%s %s::%s: %s\n", GetName(), ClassName(), __func__, Form("Selected jets in array: %d", iNJetSel));
    // fill jet spectra
    for(Int_t iJet = 0; iJet < iNJetSel; iJet++)
    {
      jet = (AliAODJet*)arrayJetSel->At(iJet); // load a jet in the list
      fh1PtJet[iCentIndex]->Fill(jet->Pt()); // pt spectrum of selected jets
      fh1EtaJet[iCentIndex]->Fill(jet->Eta()); // eta spectrum of selected jets
      fh2EtaPtJet[iCentIndex]->Fill(jet->Eta(), jet->Pt()); // eta-pT spectrum of selected jets
      fh1PhiJet[iCentIndex]->Fill(jet->Phi()); // phi spectrum of selected jets
      Double_t dAreaExcluded = TMath::Pi() * dRadiusExcludeCone * dRadiusExcludeCone; // area of the cone
      dAreaExcluded -= AreaCircSegment(dRadiusExcludeCone, fdCutEtaV0Max - jet->Eta()); // positive eta overhang
      dAreaExcluded -= AreaCircSegment(dRadiusExcludeCone, fdCutEtaV0Max + jet->Eta()); // negative eta overhang
      fh1AreaExcluded->Fill(iCentIndex, dAreaExcluded);
      // calculate distances between all pairs of selected jets in the event
      vecJetMomentum.SetXYZ(jet->Px(), jet->Py(), jet->Pz()); // set the vector of jet momentum
      for(Int_t iJetPair = iJet + 1; iJetPair < iNJetSel; iJetPair++)
      {
        if(fDebug > 4) printf("%s %s::%s: %s\n", GetName(), ClassName(), __func__, Form("Calculating distance for jet pair: %d-%d/%d\n", iJet, iJetPair, iNJetSel));
        jet = (AliAODJet*)arrayJetSel->At(iJetPair); // load another jet in the list
        vecJetMomentumPair.SetXYZ(jet->Px(), jet->Py(), jet->Pz()); // set the vector of the second jet momentum
        fh1DistanceJets[iCentIndex]->Fill(vecJetMomentum.DeltaR(vecJetMomentumPair));
      }
    }
    jet = 0; // just to be sure
  }
  fh1NJetPerEvent[iCentIndex]->Fill(iNJetSel);

  if(bJetEventGood) // there should be some reconstructed jets
  {
    fh1EventCounterCut->Fill(10); // events with jet(s)
    fh1EventCounterCutCent[iCentIndex]->Fill(10); // events with jet(s)
    if(iNJetSel)
    {
      fh1EventCounterCut->Fill(11); // events with selected jets
      fh1EventCounterCutCent[iCentIndex]->Fill(11);
    }
  }
  if(iNJetSel)
    fh1EventCent2Jets->Fill(fdCentrality);
  else
    fh1EventCent2NoJets->Fill(fdCentrality);

  if(iNJetSel)
  {
    jetRnd = GetRandomCone(arrayJetSel, dCutEtaJetMax, 2 * fdDistanceV0JetMax);
    if(jetRnd)
    {
      fh1NRndConeCent->Fill(iCentIndex);
      fh2EtaPhiRndCone[iCentIndex]->Fill(jetRnd->Eta(), jetRnd->Phi());
    }
    if(fJetsBgCont)
    {
      jetMed = GetMedianCluster(fJetsBgCont, dCutEtaJetMax);
      if(jetMed)
      {
        fh1NMedConeCent->Fill(iCentIndex);
        fh2EtaPhiMedCone[iCentIndex]->Fill(jetMed->Eta(), jetMed->Phi());
      }
    }
  }

  if(fbJetSelection && fbCompareTriggers) // Correlations of pt_jet with pt_trigger-track
  {
    if(!fTracksCont)
      AliError("No track container!");
    else
    {
      AliAODJet* jetTrig = 0;
      fTracksCont->ResetCurrentID();
      AliVTrack* track = static_cast<AliVTrack*>(fTracksCont->GetNextAcceptParticle());
      if(track) // there are some accepted trigger tracks
      {
        while(track) // loop over selected tracks
        {
          fh1PtTrigger[iCentIndex]->Fill(track->Pt());
          if(iNJetSel) // there are some accepted jets
          {
            for(Int_t iJet = 0; iJet < iNJetSel; iJet++)
            {
              jetTrig = (AliAODJet*)arrayJetSel->At(iJet);
              fh2PtJetPtTrigger[iCentIndex]->Fill(jetTrig->Pt(), track->Pt());
            }
          }
          else // there are no accepted jets
          {
            fh2PtJetPtTrigger[iCentIndex]->Fill(0., track->Pt());
          }
          track = static_cast<AliVTrack*>(fTracksCont->GetNextAcceptParticle()); // load next accepted trigger track
        }
      }
      else // there are no accepted trigger tracks
      {
        if(iNJetSel) // there are some accepted jets
        {
          for(Int_t iJet = 0; iJet < iNJetSel; iJet++)
          {
            jetTrig = (AliAODJet*)arrayJetSel->At(iJet);
            fh2PtJetPtTrigger[iCentIndex]->Fill(jetTrig->Pt(), 0.);
          }
        }
        else // there are no accepted jets
          fh2PtJetPtTrigger[iCentIndex]->Fill(0., 0.);
      }
    }
  }

  // Loading primary vertex info
  AliAODVertex* primVtx = fAODIn->GetPrimaryVertex(); // get the primary vertex
  Double_t dPrimVtxPos[3]; // primary vertex position {x,y,z}
  primVtx->GetXYZ(dPrimVtxPos);
  fh1VtxZ[iCentIndex]->Fill(dPrimVtxPos[2]);
  fh2VtxXY[iCentIndex]->Fill(dPrimVtxPos[0], dPrimVtxPos[1]);

  //===== Start of loop over V0 candidates =====
  if(fDebug > 0) printf("%s %s::%s: %s\n", GetName(), ClassName(), __func__, "Start of V0 loop");
  for(Int_t iV0 = 0; iV0 < iNV0s; iV0++)
  {
    v0 = fAODIn->GetV0(iV0); // get next candidate from the list in AOD
    if(!v0)
      continue;

    iNV0CandTot++;

    // Initialization of status indicators
    Bool_t bIsCandidateK0s = kTRUE; // candidate for K0s
    Bool_t bIsCandidateLambda = kTRUE; // candidate for Lambda
    Bool_t bIsCandidateALambda = kTRUE; // candidate for anti-Lambda
    Bool_t bIsInPeakK0s = kFALSE; // candidate within the K0s mass peak
    Bool_t bIsInPeakLambda = kFALSE; // candidate within the Lambda mass peak
    Bool_t bIsInPeakALambda = kFALSE; // candidate within the anti-Lambda mass peak
    Bool_t bIsInConeJet = kFALSE; // candidate within the jet cones
    Bool_t bIsInConePerp = kFALSE; // candidate within a perpendicular cone
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

    Double_t dPtV0 = TMath::Sqrt(v0->Pt2V0()); // transverse momentum of V0
    vecV0Momentum = TVector3(v0->Px(), v0->Py(), v0->Pz()); // set the vector of V0 momentum

    // Sigma of the mass peak window
    Double_t dMassPeakWindowK0s = dNSigmaMassMax * MassPeakSigmaOld(dPtV0, 0);
    Double_t dMassPeakWindowLambda = dNSigmaMassMax * MassPeakSigmaOld(dPtV0, 1);
    if(!fbIsPbPb) // p-p
    {
      dMassPeakWindowK0s = 0.010; // LF p-p
      dMassPeakWindowLambda = 0.005; // LF p-p
    }

    // Invariant mass peak selection
    if(TMath::Abs(dMassV0K0s - dMassPDGK0s) < dMassPeakWindowK0s)
      bIsInPeakK0s = kTRUE;
    if(TMath::Abs(dMassV0Lambda - dMassPDGLambda) < dMassPeakWindowLambda)
      bIsInPeakLambda = kTRUE;
    if(TMath::Abs(dMassV0ALambda - dMassPDGLambda) < dMassPeakWindowLambda)
      bIsInPeakALambda = kTRUE;

    // QA histograms before cuts
    FillQAHistogramV0(primVtx, v0, 0, bIsCandidateK0s, bIsCandidateLambda, bIsCandidateALambda, bIsInPeakK0s, bIsInPeakLambda, bIsInPeakALambda);

    // Skip candidates outside the histogram range
    if((dMassV0K0s < fgkdMassK0sMin) || (dMassV0K0s >= fgkdMassK0sMax))
      bIsCandidateK0s = kFALSE;
    if((dMassV0Lambda < fgkdMassLambdaMin) || (dMassV0Lambda >= fgkdMassLambdaMax))
      bIsCandidateLambda = kFALSE;
    if((dMassV0ALambda < fgkdMassLambdaMin) || (dMassV0ALambda >= fgkdMassLambdaMax))
      bIsCandidateALambda = kFALSE;
    if(!bIsCandidateK0s && !bIsCandidateLambda && !bIsCandidateALambda)
      continue;

    // Retrieving all relevant properties of the V0 candidate
    Bool_t bOnFlyStatus = v0->GetOnFlyStatus(); // online (on fly) reconstructed vs offline reconstructed
    const AliAODTrack* trackPos = (AliAODTrack*)v0->GetDaughter(0); // positive daughter track
    const AliAODTrack* trackNeg = (AliAODTrack*)v0->GetDaughter(1); // negative daughter track
    Double_t dPtDaughterPos = trackPos->Pt(); // transverse momentum of a daughter track calculated as if primary, != v0->PtProng(0)
    Double_t dPtDaughterNeg = trackNeg->Pt(); // != v0->PtProng(1)
    Double_t dNRowsPos = trackPos->GetTPCClusterInfo(2, 1); // crossed TPC pad rows of a daughter track
    Double_t dNRowsNeg = trackNeg->GetTPCClusterInfo(2, 1);
    Double_t dFindablePos = Double_t(trackPos->GetTPCNclsF()); // Findable clusters
    Double_t dFindableNeg = Double_t(trackNeg->GetTPCNclsF());
    Double_t dDCAToPrimVtxPos = TMath::Abs(v0->DcaPosToPrimVertex()); // dca of a daughter to the primary vertex
    Double_t dDCAToPrimVtxNeg = TMath::Abs(v0->DcaNegToPrimVertex());
    Double_t dDCADaughters = v0->DcaV0Daughters(); // dca between daughters
    Double_t dCPA = v0->CosPointingAngle(primVtx); // cosine of the pointing angle
    Double_t dSecVtxPos[3]; // V0 vertex position {x,y,z}
//      Double_t dSecVtxPos[3] = {v0->DecayVertexV0X(),v0->DecayVertexV0Y(),v0->DecayVertexV0Z()}; // V0 vertex position
    v0->GetSecondaryVtx(dSecVtxPos);
    Double_t dRadiusDecay = TMath::Sqrt(dSecVtxPos[0] * dSecVtxPos[0] + dSecVtxPos[1] * dSecVtxPos[1]); // distance of the V0 vertex from the z-axis
    Double_t dEtaDaughterPos = trackPos->Eta(); // pseudorapidity of a daughter track calculated as if primary, != v0->EtaProng(0)
    Double_t dEtaDaughterNeg = trackNeg->Eta(); // != v0->EtaProng(1);
    Double_t dRapK0s = v0->RapK0Short(); // rapidity calculated for K0s assumption
    Double_t dRapLambda = v0->RapLambda(); // rapidity calculated for Lambda assumption
    Double_t dEtaV0 = v0->Eta(); // V0 pseudorapidity
    Double_t dPhiV0 = v0->Phi(); // V0 azimuth
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
    Double_t dNSigmaPosPion   = (fPIDResponse ? TMath::Abs(fPIDResponse->NumberOfSigmasTPC(trackPos, AliPID::kPion)) : 0.); // difference between measured and expected signal of the dE/dx in the TPC
    Double_t dNSigmaPosProton = (fPIDResponse ? TMath::Abs(fPIDResponse->NumberOfSigmasTPC(trackPos, AliPID::kProton)) : 0.);
    Double_t dNSigmaNegPion   = (fPIDResponse ? TMath::Abs(fPIDResponse->NumberOfSigmasTPC(trackNeg, AliPID::kPion)) : 0.);
    Double_t dNSigmaNegProton = (fPIDResponse ? TMath::Abs(fPIDResponse->NumberOfSigmasTPC(trackNeg, AliPID::kProton)) : 0.);
    Double_t dAlpha = v0->AlphaV0(); // Armenteros-Podolanski alpha
    Double_t dPtArm = v0->PtArmV0(); // Armenteros-Podolanski pT
    AliAODVertex* prodVtxDaughterPos = (AliAODVertex*)(trackPos->GetProdVertex()); // production vertex of the positive daughter track
    Char_t cTypeVtxProdPos = prodVtxDaughterPos->GetType(); // type of the production vertex
    AliAODVertex* prodVtxDaughterNeg = (AliAODVertex*)(trackNeg->GetProdVertex()); // production vertex of the negative daughter track
    Char_t cTypeVtxProdNeg = prodVtxDaughterNeg->GetType(); // type of the production vertex

//    fh2Tau3DVs2D[0]->Fill(dPtV0, dLOverP / dROverPt);

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
    if(bPrintCuts) printf("Rec: Applying cut: Reconstruction method: %s\n", (fbOnFly ? "on-the-fly" : "offline"));
    if(bOnFlyStatus != fbOnFly)
      continue;
    FillCandidates(dMassV0K0s, dMassV0Lambda, dMassV0ALambda, bIsCandidateK0s, bIsCandidateLambda, bIsCandidateALambda, iCutIndex, iCentIndex);
    iCutIndex++;

    // 3
    // Tracks TPC OK
    if(bPrintCuts) printf("Rec: Applying cut: Correct charge of daughters\n");
    if(!trackNeg || !trackPos)
      continue;
    if(trackNeg->Charge() == trackPos->Charge()) // daughters have different charge?
      continue;
    if(trackNeg->Charge() != -1) // daughters have expected charge?
      continue;
    if(trackPos->Charge() != 1) // daughters have expected charge?
      continue;

    if(fbTPCRefit)
    {
      if(bPrintCuts) printf("Rec: Applying cut: TPC refit\n");
      if(!trackNeg->IsOn(AliAODTrack::kTPCrefit)) // TPC refit is ON?
        continue;
      if(!trackPos->IsOn(AliAODTrack::kTPCrefit))
        continue;
    }

    if(fbRejectKinks)
    {
      if(bPrintCuts) printf("Rec: Applying cut: Type of production vertex of daughter: No kinks\n");
      if(cTypeVtxProdNeg == AliAODVertex::kKink) // kink daughter rejection
        continue;
      if(cTypeVtxProdPos == AliAODVertex::kKink)
        continue;
    }

    if(fbFindableClusters)
    {
      if(bPrintCuts) printf("Rec: Applying cut: Positive number of findable clusters\n");
      if(dFindableNeg <= 0.)
        continue;
      if(dFindablePos <= 0.)
        continue;
    }

    if(fdCutNCrossedRowsTPCMin > 0.)
    {
      if(bPrintCuts) printf("Rec: Applying cut: Number of TPC rows >= %g\n", fdCutNCrossedRowsTPCMin);
      if(dNRowsNeg < fdCutNCrossedRowsTPCMin) // Crossed TPC padrows
        continue;
      if(dNRowsPos < fdCutNCrossedRowsTPCMin)
        continue;
    }

    if(fdCutCrossedRowsOverFindMin > 0.)
    {
      if(bPrintCuts) printf("Rec: Applying cut: rows/findable >= %g\n", fdCutCrossedRowsOverFindMin);
      if(dNRowsNeg / dFindableNeg < fdCutCrossedRowsOverFindMin)
        continue;
      if(dNRowsPos / dFindablePos < fdCutCrossedRowsOverFindMin)
        continue;
    }

    if(fdCutCrossedRowsOverFindMax > 0.)
    {
      if(bPrintCuts) printf("Rec: Applying cut: rows/findable <= %g\n", fdCutCrossedRowsOverFindMax);
      if(dNRowsNeg / dFindableNeg > fdCutCrossedRowsOverFindMax)
        continue;
      if(dNRowsPos / dFindablePos > fdCutCrossedRowsOverFindMax)
        continue;
    }

    FillCandidates(dMassV0K0s, dMassV0Lambda, dMassV0ALambda, bIsCandidateK0s, bIsCandidateLambda, bIsCandidateALambda, iCutIndex, iCentIndex);
    iCutIndex++;

    // 4
    // Daughters: transverse momentum cut
    if(fdCutPtDaughterMin > 0.)
    {
      if(bPrintCuts) printf("Rec: Applying cut: Daughter pt >= %g\n", fdCutPtDaughterMin);
      if((dPtDaughterNeg < fdCutPtDaughterMin) || (dPtDaughterPos < fdCutPtDaughterMin))
        continue;
      FillCandidates(dMassV0K0s, dMassV0Lambda, dMassV0ALambda, bIsCandidateK0s, bIsCandidateLambda, bIsCandidateALambda, iCutIndex, iCentIndex);
    }
    iCutIndex++;

    // 5
    // Daughters: Impact parameter of daughters to prim vtx
    if(fdCutDCAToPrimVtxMin > 0.)
    {
      if(bPrintCuts) printf("Rec: Applying cut: Daughter DCA to prim vtx >= %g\n", fdCutDCAToPrimVtxMin);
      if((dDCAToPrimVtxNeg < fdCutDCAToPrimVtxMin) || (dDCAToPrimVtxPos < fdCutDCAToPrimVtxMin))
        continue;
      FillCandidates(dMassV0K0s, dMassV0Lambda, dMassV0ALambda, bIsCandidateK0s, bIsCandidateLambda, bIsCandidateALambda, iCutIndex, iCentIndex);
    }
    iCutIndex++;

    // 6
    // Daughters: DCA
    if(fdCutDCADaughtersMax > 0.)
    {
      if(bPrintCuts) printf("Rec: Applying cut: DCA between daughters <= %g\n", fdCutDCADaughtersMax);
      if(dDCADaughters > fdCutDCADaughtersMax)
        continue;
      FillCandidates(dMassV0K0s, dMassV0Lambda, dMassV0ALambda, bIsCandidateK0s, bIsCandidateLambda, bIsCandidateALambda, iCutIndex, iCentIndex);
    }
    iCutIndex++;

    // 7
    // V0: Cosine of the pointing angle
    if(fdCutCPAKMin > 0.)
    {
      if(bPrintCuts) printf("Rec: Applying cut: CPA >= %g (K)\n", fdCutCPAKMin);
      if(dCPA < fdCutCPAKMin)
        bIsCandidateK0s = kFALSE;
    }
    if(fdCutCPALMin > 0.)
    {
      if(bPrintCuts) printf("Rec: Applying cut: CPA >= %g (L, AL)\n", fdCutCPALMin);
      if(dCPA < fdCutCPALMin)
      {
        bIsCandidateLambda = kFALSE;
        bIsCandidateALambda = kFALSE;
      }
    }
    if(fdCutCPAKMin > 0. || fdCutCPALMin > 0.)
      FillCandidates(dMassV0K0s, dMassV0Lambda, dMassV0ALambda, bIsCandidateK0s, bIsCandidateLambda, bIsCandidateALambda, iCutIndex, iCentIndex);
    iCutIndex++;

    // 8
    // V0: Fiducial volume
    if(fdCutRadiusDecayMin > 0. && fdCutRadiusDecayMax > 0.)
    {
      if(bPrintCuts) printf("Rec: Applying cut: Decay radius >= %g, <= %g\n", fdCutRadiusDecayMin, fdCutRadiusDecayMax);
      if((dRadiusDecay < fdCutRadiusDecayMin) || (dRadiusDecay > fdCutRadiusDecayMax))
        continue;
      FillCandidates(dMassV0K0s, dMassV0Lambda, dMassV0ALambda, bIsCandidateK0s, bIsCandidateLambda, bIsCandidateALambda, iCutIndex, iCentIndex);
    }
    iCutIndex++;

    // 9
    // Daughters: pseudorapidity cut
    if(fdCutEtaDaughterMax > 0.)
    {
      if(bPrintCuts) printf("Rec: Applying cut: Daughter |eta| < %g\n", fdCutEtaDaughterMax);
      if((TMath::Abs(dEtaDaughterNeg) > fdCutEtaDaughterMax) || (TMath::Abs(dEtaDaughterPos) > fdCutEtaDaughterMax))
        continue;
      FillCandidates(dMassV0K0s, dMassV0Lambda, dMassV0ALambda, bIsCandidateK0s, bIsCandidateLambda, bIsCandidateALambda, iCutIndex, iCentIndex);
    }
    iCutIndex++;
    // End of global cuts

    // Start of particle-dependent cuts
    // 10
    // V0: pseudorapidity cut & rapidity cut
    if(fdCutEtaV0Max > 0.)
    {
      if(bPrintCuts) printf("Rec: Applying cut: V0 |eta| < %g\n", fdCutEtaV0Max);
      if(TMath::Abs(dEtaV0) > fdCutEtaV0Max)
      {
        bIsCandidateK0s = kFALSE;
        bIsCandidateLambda = kFALSE;
        bIsCandidateALambda = kFALSE;
      }
    }
    if(fdCutRapV0Max > 0.)
    {
      if(bPrintCuts) printf("Rec: Applying cut: V0 |y| < %g\n", fdCutRapV0Max);
      if(TMath::Abs(dRapK0s) > fdCutRapV0Max)
        bIsCandidateK0s = kFALSE;
      if(TMath::Abs(dRapLambda) > fdCutRapV0Max)
      {
        bIsCandidateLambda = kFALSE;
        bIsCandidateALambda = kFALSE;
      }
    }
    if(fdCutEtaV0Max > 0. || fdCutRapV0Max > 0.)
      FillCandidates(dMassV0K0s, dMassV0Lambda, dMassV0ALambda, bIsCandidateK0s, bIsCandidateLambda, bIsCandidateALambda, iCutIndex, iCentIndex);
    iCutIndex++;

    // 11
    // Lifetime cut
    if(fdCutNTauKMax > 0.)
    {
      if(bPrintCuts) printf("Rec: Applying cut: Proper lifetime < %g (K)\n", fdCutNTauKMax);
      if(dMROverPtK0s > fdCutNTauKMax * dCTauK0s)
        bIsCandidateK0s = kFALSE;
    }
    if(fdCutNTauLMax > 0.)
    {
      if(bPrintCuts) printf("Rec: Applying cut: Proper lifetime < %g (L, AL)\n", fdCutNTauLMax);
      if(dMROverPtLambda > fdCutNTauLMax * dCTauLambda)
      {
        bIsCandidateLambda = kFALSE;
        bIsCandidateALambda = kFALSE;
      }
    }
    if(fdCutNTauKMax > 0. || fdCutNTauLMax > 0.)
      FillCandidates(dMassV0K0s, dMassV0Lambda, dMassV0ALambda, bIsCandidateK0s, bIsCandidateLambda, bIsCandidateALambda, iCutIndex, iCentIndex);
    iCutIndex++;

    // 12
    // Daughter PID
    if(fdCutNSigmadEdxMax > 0.)
    {
      if(fbIsPbPb && fdPtProtonPIDMax > 0.) // Pb-Pb
      {
        if(bPrintCuts) printf("Rec: Applying cut: Delta dE/dx (proton below %g GeV/c) < %g\n", fdPtProtonPIDMax, fdCutNSigmadEdxMax);
        if((dPtDaughterPos < fdPtProtonPIDMax) && (dNSigmaPosProton > fdCutNSigmadEdxMax)) // p+
          bIsCandidateLambda = kFALSE;
        if((dPtDaughterNeg < fdPtProtonPIDMax) && (dNSigmaNegProton > fdCutNSigmadEdxMax)) // p-
          bIsCandidateALambda = kFALSE;
      }
      else // p-p
      {
        if(bPrintCuts) printf("Rec: Applying cut: Delta dE/dx (both daughters): < %g\n", fdCutNSigmadEdxMax);
        if(dNSigmaPosPion > fdCutNSigmadEdxMax || dNSigmaNegPion > fdCutNSigmadEdxMax) // pi+, pi-
          bIsCandidateK0s = kFALSE;
        if(dNSigmaPosProton > fdCutNSigmadEdxMax || dNSigmaNegPion > fdCutNSigmadEdxMax) // p+, pi-
          bIsCandidateLambda = kFALSE;
        if(dNSigmaNegProton > fdCutNSigmadEdxMax || dNSigmaPosPion > fdCutNSigmadEdxMax) // p-, pi+
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
    if(fbCutArmPod)
    {
      if(bPrintCuts) printf("Rec: Applying cut: Armenteros-Podolanski (K0S) pT > %g * |alpha|\n", 0.2);
      if(dPtArm < TMath::Abs(0.2 * dAlpha))
        bIsCandidateK0s = kFALSE;
      FillCandidates(dMassV0K0s, dMassV0Lambda, dMassV0ALambda, bIsCandidateK0s, bIsCandidateLambda, bIsCandidateALambda, iCutIndex, iCentIndex);
    }
    iCutIndex++;

    // 14
    // Cross-contamination
    /*
    if(bIsInPeakK0s)
    {
      if(bIsCandidateLambda) // Lambda candidates in K0s peak, excluded from Lambda candidates by CC cut
        fh2CCLambda->Fill(dMassV0Lambda, dPtV0);
    }
    if(bIsInPeakLambda)
    {
      if(bIsCandidateK0s) // K0s candidates in Lambda peak, excluded from K0s candidates by CC cut
        fh2CCK0s->Fill(dMassV0K0s, dPtV0);
    }
    */
    if(fbCutCross)
    {
      if(bIsInPeakK0s)
      {
        bIsCandidateLambda = kFALSE;
        bIsCandidateALambda = kFALSE;
      }
      if(bIsInPeakLambda)
      {
        bIsCandidateK0s = kFALSE;
      }
      if(bIsInPeakALambda)
      {
        bIsCandidateK0s = kFALSE;
      }
      FillCandidates(dMassV0K0s, dMassV0Lambda, dMassV0ALambda, bIsCandidateK0s, bIsCandidateLambda, bIsCandidateALambda, iCutIndex, iCentIndex);
    }
    iCutIndex++;
    // End of particle-dependent cuts

    //===== End of reconstruction cutting =====

    if(!bIsCandidateK0s && !bIsCandidateLambda && !bIsCandidateALambda)
      continue;

    // Selection of V0s in jet cones, perpendicular cones, random cones, outside cones
    if(iNJetSel && (bIsCandidateK0s || bIsCandidateLambda || bIsCandidateALambda))
    {
      Double_t dDMin, dD = 0; // minimal / current value of V0-jet Distance (used for estimation of closest jet to V0)
      dDMin = 20.;
      for(Int_t iJet = 0; iJet < iNJetSel; iJet++) // could be included in loop beneath
      {
        // finding the closest jet to the v0
        jet = (AliAODJet*)arrayJetSel->At(iJet); // load a jet in the list
        if(!jet)
          continue;
        dD = GetD(v0, jet);
        if(dD < dDMin)
          dDMin = dD;
      }
      if(bIsCandidateK0s)
        fh1DistanceV0JetsK0s[iCentIndex]->Fill(dDMin);
      if(bIsCandidateLambda)
        fh1DistanceV0JetsLambda[iCentIndex]->Fill(dDMin);
      if(bIsCandidateALambda)
        fh1DistanceV0JetsALambda[iCentIndex]->Fill(dDMin);

      // Selection of V0s in jet cones
      if(fDebug > 4) printf("%s %s::%s: %s\n", GetName(), ClassName(), __func__, Form("Searching for V0 %d %d in %d jet cones", bIsCandidateK0s, bIsCandidateLambda, iNJetSel));
      for(Int_t iJet = 0; iJet < iNJetSel; iJet++)
      {
        jet = (AliAODJet*)arrayJetSel->At(iJet); // load a jet in the list
        if(!jet)
          continue;
        vecJetMomentum.SetXYZ(jet->Px(), jet->Py(), jet->Pz()); // set the vector of jet momentum
        if(fDebug > 4) printf("%s %s::%s: %s\n", GetName(), ClassName(), __func__, Form("Checking if V0 %d %d in jet cone %d", bIsCandidateK0s, bIsCandidateLambda, iJet));
        if(IsParticleInCone(v0, jet, fdDistanceV0JetMax)) // If good jet in event, find out whether V0 is in that jet
        {
          if(fDebug > 4) printf("%s %s::%s: %s\n", GetName(), ClassName(), __func__, Form("V0 %d %d found in jet cone %d", bIsCandidateK0s, bIsCandidateLambda, iJet));
          bIsInConeJet = kTRUE;
          dPtJetTrackLeading = jet->GetPtLeading();
          dPtJet = jet->Pt();
          break;
        }
      }
      // Selection of V0s in perp. cones
      if(fDebug > 4) printf("%s %s::%s: %s\n", GetName(), ClassName(), __func__, Form("Searching for V0 %d %d in %d perp. cones", bIsCandidateK0s, bIsCandidateLambda, iNJetPerp));
      for(Int_t iJet = 0; iJet < iNJetPerp; iJet++)
      {
        jetPerp = (AliAODJet*)arrayJetPerp->At(iJet); // load a jet in the list
        if(fDebug > 4) printf("%s %s::%s: %s\n", GetName(), ClassName(), __func__, Form("Checking if V0 %d %d in perp. cone %d", bIsCandidateK0s, bIsCandidateLambda, iJet));
        if(IsParticleInCone(v0, jetPerp, fdDistanceV0JetMax)) // V0 in perp. cone
        {
          if(fDebug > 4) printf("%s %s::%s: %s\n", GetName(), ClassName(), __func__, Form("V0 %d %d found in perp. cone %d", bIsCandidateK0s, bIsCandidateLambda, iJet));
          bIsInConePerp = kTRUE;
          break;
        }
      }
      // Selection of V0s in random cones
      if(jetRnd)
      {
        if(fDebug > 4) printf("%s %s::%s: %s\n", GetName(), ClassName(), __func__, Form("Searching for V0 %d %d in the rnd. cone", bIsCandidateK0s, bIsCandidateLambda));
        if(IsParticleInCone(v0, jetRnd, fdDistanceV0JetMax)) // V0 in rnd. cone?
        {
          if(fDebug > 4) printf("%s %s::%s: %s\n", GetName(), ClassName(), __func__, Form("V0 %d %d found in the rnd. cone", bIsCandidateK0s, bIsCandidateLambda));
          bIsInConeRnd = kTRUE;
        }
      }
      // Selection of V0s in median-cluster cones
      if(jetMed)
      {
        if(fDebug > 4) printf("%s %s::%s: %s\n", GetName(), ClassName(), __func__, Form("Searching for V0 %d %d in the med. cone", bIsCandidateK0s, bIsCandidateLambda));
        if(IsParticleInCone(v0, jetMed, fdDistanceV0JetMax)) // V0 in med. cone?
        {
          if(fDebug > 4) printf("%s %s::%s: %s\n", GetName(), ClassName(), __func__, Form("V0 %d %d found in the med. cone", bIsCandidateK0s, bIsCandidateLambda));
          bIsInConeMed = kTRUE;
        }
      }
      // Selection of V0s outside jet cones
      if(fDebug > 4) printf("%s %s::%s: %s\n", GetName(), ClassName(), __func__, Form("Searching for V0 %d %d outside jet cones", bIsCandidateK0s, bIsCandidateLambda));
      if(!OverlapWithJets(arrayJetSel, v0, dRadiusExcludeCone)) // V0 oustide jet cones
      {
        if(fDebug > 4) printf("%s %s::%s: %s\n", GetName(), ClassName(), __func__, Form("V0 %d %d found outside jet cones", bIsCandidateK0s, bIsCandidateLambda));
        bIsOutsideCones = kTRUE;
      }
    }

    // QA histograms after cuts
    FillQAHistogramV0(primVtx, v0, 1, bIsCandidateK0s, bIsCandidateLambda, bIsCandidateALambda, bIsInPeakK0s, bIsInPeakLambda, bIsInPeakALambda);
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
    }
    */

    //===== Start of filling V0 spectra =====

    Double_t dAngle = TMath::Pi(); // angle between V0 momentum and jet momentum
    if(bIsInConeJet)
    {
      dAngle = vecV0Momentum.Angle(vecJetMomentum);
    }

    // iCutIndex = 15
    if(bIsCandidateK0s)
    {
      // 15 K0s candidates after cuts
//          printf("K0S: i = %d, m = %g, pT = %g, eta = %g, phi = %g\n",iV0,dMassV0K0s,dPtV0,dEtaV0,dPhiV0);
      FillCandidates(dMassV0K0s, dMassV0Lambda, dMassV0ALambda, bIsCandidateK0s, kFALSE, kFALSE, iCutIndex, iCentIndex);
      Double_t valueKIncl[3] = {dMassV0K0s, dPtV0, dEtaV0};
      fhnV0InclusiveK0s[iCentIndex]->Fill(valueKIncl);
      fh1V0InvMassK0sCent[iCentIndex]->Fill(dMassV0K0s);

      fh1QACTau2D[1]->Fill(dMROverPtK0s / dCTauK0s);
      fh1QACTau3D[1]->Fill(dMLOverPK0s / dCTauK0s);
//      fh2Tau3DVs2D[1]->Fill(dPtV0, dLOverP / dROverPt);

      if(iNJetSel)
      {
        // 16 K0s in jet events
        FillCandidates(dMassV0K0s, dMassV0Lambda, dMassV0ALambda, bIsCandidateK0s, kFALSE, kFALSE, iCutIndex + 1, iCentIndex);
      }
      if(bIsInConeJet)
      {
        // 17 K0s in jets
        FillCandidates(dMassV0K0s, dMassV0Lambda, dMassV0ALambda, bIsCandidateK0s, kFALSE, kFALSE, iCutIndex + 2, iCentIndex);
        Double_t valueKInJC[4] = {dMassV0K0s, dPtV0, dEtaV0, jet->Pt()};
        fhnV0InJetK0s[iCentIndex]->Fill(valueKInJC);
        fh2V0PtJetAngleK0s[iCentIndex]->Fill(jet->Pt(), dAngle);
        Double_t valuesDaughter[5] = {dPtDaughterPos, dPtDaughterNeg, dPtV0, dPtJet, dPtJetTrackLeading};
        fhnPtDaughterK0s[iCentIndex]->Fill(valuesDaughter);
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
      // 15 Lambda candidates after cuts
//          printf("La: i = %d, m = %g, pT = %g, eta = %g, phi = %g\n",iV0,dMassV0Lambda,dPtV0,dEtaV0,dPhiV0);
      FillCandidates(dMassV0K0s, dMassV0Lambda, dMassV0ALambda, kFALSE, bIsCandidateLambda, kFALSE, iCutIndex, iCentIndex);
      Double_t valueLIncl[3] = {dMassV0Lambda, dPtV0, dEtaV0};
      fhnV0InclusiveLambda[iCentIndex]->Fill(valueLIncl);
      fh1V0InvMassLambdaCent[iCentIndex]->Fill(dMassV0Lambda);
      if(iNJetSel)
      {
        // 16 Lambda in jet events
        FillCandidates(dMassV0K0s, dMassV0Lambda, dMassV0ALambda, kFALSE, bIsCandidateLambda, kFALSE, iCutIndex + 1, iCentIndex);
      }
      if(bIsInConeJet)
      {
        // 17 Lambda in jets
        FillCandidates(dMassV0K0s, dMassV0Lambda, dMassV0ALambda, kFALSE, bIsCandidateLambda, kFALSE, iCutIndex + 2, iCentIndex);
        Double_t valueLInJC[4] = {dMassV0Lambda, dPtV0, dEtaV0, jet->Pt()};
        fhnV0InJetLambda[iCentIndex]->Fill(valueLInJC);
        fh2V0PtJetAngleLambda[iCentIndex]->Fill(jet->Pt(), dAngle);
        Double_t valuesDaughter[5] = {dPtDaughterPos, dPtDaughterNeg, dPtV0, dPtJet, dPtJetTrackLeading};
        fhnPtDaughterLambda[iCentIndex]->Fill(valuesDaughter);
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
      // 15 ALambda candidates after cuts
//          printf("AL: i = %d, m = %g, pT = %g, eta = %g, phi = %g\n",iV0,dMassV0ALambda,dPtV0,dEtaV0,dPhiV0);
      FillCandidates(dMassV0K0s, dMassV0Lambda, dMassV0ALambda, kFALSE, kFALSE, bIsCandidateALambda, iCutIndex, iCentIndex);
      Double_t valueALIncl[3] = {dMassV0ALambda, dPtV0, dEtaV0};
      fhnV0InclusiveALambda[iCentIndex]->Fill(valueALIncl);
      fh1V0InvMassALambdaCent[iCentIndex]->Fill(dMassV0ALambda);
      if(iNJetSel)
      {
        // 16 ALambda in jet events
        FillCandidates(dMassV0K0s, dMassV0Lambda, dMassV0ALambda, kFALSE, kFALSE, bIsCandidateALambda, iCutIndex + 1, iCentIndex);
      }
      if(bIsInConeJet)
      {
        // 17 ALambda in jets
        FillCandidates(dMassV0K0s, dMassV0Lambda, dMassV0ALambda, kFALSE, kFALSE, bIsCandidateALambda, iCutIndex + 2, iCentIndex);
        Double_t valueLInJC[4] = {dMassV0ALambda, dPtV0, dEtaV0, jet->Pt()};
        fhnV0InJetALambda[iCentIndex]->Fill(valueLInJC);
        fh2V0PtJetAngleALambda[iCentIndex]->Fill(jet->Pt(), dAngle);
        Double_t valuesDaughter[5] = {dPtDaughterPos, dPtDaughterNeg, dPtV0, dPtJet, dPtJetTrackLeading};
        fhnPtDaughterALambda[iCentIndex]->Fill(valuesDaughter);
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
    // V0-jet correlations
    if(fbCorrelations && iNJetSel)
    {
      Double_t dDPhi, dDEta;
      // Fill V0-jet correlations in same events
      for(Int_t iJet = 0; iJet < iNJetSel; iJet++)
      {
        AliAODJet* jetCorrel = (AliAODJet*)arrayJetSel->At(iJet); // load a jet in the list
        dDPhi = GetNormalPhi(dPhiV0 - jetCorrel->Phi());
        dDEta = dEtaV0 - jetCorrel->Eta();
        if(TMath::Abs(dDEta) > fdDeltaEtaMax)
          continue;
        if(bIsCandidateK0s)
        {
          Double_t valueKCorrel[6] = {dMassV0K0s, dPtV0, dEtaV0, jetCorrel->Pt(), dDPhi, dDEta};
          fhnV0CorrelSEK0s[iCentIndex]->Fill(valueKCorrel);
        }
        if(bIsCandidateLambda)
        {
          Double_t valueLCorrel[6] = {dMassV0Lambda, dPtV0, dEtaV0, jetCorrel->Pt(), dDPhi, dDEta};
          fhnV0CorrelSELambda[iCentIndex]->Fill(valueLCorrel);
        }
      }
      // Fill V0-jet correlations in mixed events
      if(bPoolReady)
      {
        for(Int_t iMix = 0; iMix < pool->GetCurrentNEvents(); iMix++)
        {
          TObjArray* arrayMixedEvent = pool->GetEvent(iMix);
          if(!arrayMixedEvent)
            continue;
          for(Int_t iJet = 0; iJet < arrayMixedEvent->GetEntriesFast(); iJet++)
          {
            TLorentzVector* jetMixed = (TLorentzVector*)arrayMixedEvent->At(iJet);
            if(!jetMixed)
              continue;
            dDPhi = GetNormalPhi(dPhiV0 - jetMixed->Phi());
            dDEta = dEtaV0 - jetMixed->Eta();
            if(TMath::Abs(dDEta) > fdDeltaEtaMax)
              continue;
            if(bIsCandidateK0s)
            {
              Double_t valueKCorrel[6] = {dMassV0K0s, dPtV0, dEtaV0, jetMixed->Pt(), dDPhi, dDEta};
              fhnV0CorrelMEK0s[iCentIndex]->Fill(valueKCorrel);
            }
            if(bIsCandidateLambda)
            {
              Double_t valueLCorrel[6] = {dMassV0Lambda, dPtV0, dEtaV0, jetMixed->Pt(), dDPhi, dDEta};
              fhnV0CorrelMELambda[iCentIndex]->Fill(valueLCorrel);
            }
          }
        }
      }
    }
    //===== End of filling V0 spectra =====


    //===== Association of reconstructed V0 candidates with MC particles =====
    if(fbMCAnalysis)
    {
      // Associate selected candidates only
//          if ( !(bIsCandidateK0s && bIsInPeakK0s) && !(bIsCandidateLambda && bIsInPeakLambda) ) // signal candidates
      if(!(bIsCandidateK0s) && !(bIsCandidateLambda)  && !(bIsCandidateALambda)) // chosen candidates with any mass
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
      Double_t dRapV0Gen = particleMCMother->Y();
      Double_t dEtaV0Gen = particleMCMother->Eta();
//      Double_t dPhiV0Gen = particleMCMother->Phi();

      // V0 pseudorapidity cut applied on generated particles
      if(fdCutEtaV0Max > 0.)
      {
        if(bPrintCuts) printf("Rec->Gen: Applying cut: V0 |eta|: < %g\n", fdCutEtaV0Max);
        if((TMath::Abs(dEtaV0Gen) > fdCutEtaV0Max))
          continue;
      }
      // V0 rapidity cut applied on generated particles
      if(fdCutRapV0Max > 0.)
      {
        if(bPrintCuts) printf("Rec->Gen: Applying cut: V0 |y|: < %g\n", fdCutRapV0Max);
        if((TMath::Abs(dRapV0Gen) > fdCutRapV0Max))
          continue;
      }

      // Select only particles from a specific generator
      if(!IsFromGoodGenerator(iIndexMotherPos))
        continue;

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

      // K0s
//          if (bIsCandidateK0s && bIsInPeakK0s) // selected candidates in peak
      if(bIsCandidateK0s) // selected candidates with any mass
      {
//              if (bV0MCIsK0s && bV0MCIsPrimary) // well reconstructed candidates
        if(bV0MCIsK0s && bV0MCIsPrimaryDist) // well reconstructed candidates
        {
          fh2V0K0sPtMassMCRec[iCentIndex]->Fill(dPtV0Gen, dMassV0K0s);
          Double_t valueEtaK[3] = {dMassV0K0s, dPtV0Gen, dEtaV0Gen};
          fh3V0K0sEtaPtMassMCRec[iCentIndex]->Fill(valueEtaK);

          Double_t valueEtaDKNeg[6] = {0, particleMCDaughterNeg->Eta(), particleMCDaughterNeg->Pt(), dEtaV0Gen, dPtV0Gen, 0};
          fhnV0K0sInclDaughterEtaPtPtMCRec[iCentIndex]->Fill(valueEtaDKNeg);
          Double_t valueEtaDKPos[6] = {1, particleMCDaughterPos->Eta(), particleMCDaughterPos->Pt(), dEtaV0Gen, dPtV0Gen, 0};
          fhnV0K0sInclDaughterEtaPtPtMCRec[iCentIndex]->Fill(valueEtaDKPos);

          fh2V0K0sMCResolMPt[iCentIndex]->Fill(dMassV0K0s - dMassPDGK0s, dPtV0);
          fh2V0K0sMCPtGenPtRec[iCentIndex]->Fill(dPtV0Gen, dPtV0);
          if(bIsInConeJet) // true V0 associated to a candidate in jet
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
        if(bV0MCIsK0s && !bV0MCIsPrimaryDist) // not primary K0s
        {
          fh1V0K0sPtMCRecFalse[iCentIndex]->Fill(dPtV0Gen);
        }
      }
      // Lambda
//          if (bIsCandidateLambda && bIsInPeakLambda) // selected candidates in peak
      if(bIsCandidateLambda) // selected candidates with any mass
      {
//              if (bV0MCIsLambda && bV0MCIsPrimaryLambda) // well reconstructed candidates
        if(bV0MCIsLambda && bV0MCIsPrimaryDist) // well reconstructed candidates
        {
          fh2V0LambdaPtMassMCRec[iCentIndex]->Fill(dPtV0Gen, dMassV0Lambda);
          Double_t valueEtaL[3] = {dMassV0Lambda, dPtV0Gen, dEtaV0Gen};
          fh3V0LambdaEtaPtMassMCRec[iCentIndex]->Fill(valueEtaL);

          Double_t valueEtaDLNeg[6] = {0, particleMCDaughterNeg->Eta(), particleMCDaughterNeg->Pt(), dEtaV0Gen, dPtV0Gen, 0};
          fhnV0LambdaInclDaughterEtaPtPtMCRec[iCentIndex]->Fill(valueEtaDLNeg);
          Double_t valueEtaDLPos[6] = {1, particleMCDaughterPos->Eta(), particleMCDaughterPos->Pt(), dEtaV0Gen, dPtV0Gen, 0};
          fhnV0LambdaInclDaughterEtaPtPtMCRec[iCentIndex]->Fill(valueEtaDLPos);

          fh2V0LambdaMCResolMPt[iCentIndex]->Fill(dMassV0Lambda - dMassPDGLambda, dPtV0);
          fh2V0LambdaMCPtGenPtRec[iCentIndex]->Fill(dPtV0Gen, dPtV0);
          if(bIsInConeJet) // true V0 associated to a reconstructed candidate in jet
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
        if(bV0MCIsLambda && !bV0MCIsPrimaryDist && !bV0MCComesFromXi) // not primary Lambda
        {
          fh1V0LambdaPtMCRecFalse[iCentIndex]->Fill(dPtV0Gen);
        }
      }
      // anti-Lambda
//          if (bIsCandidateALambda && bIsInPeakALambda) // selected candidates in peak
      if(bIsCandidateALambda) // selected candidates with any mass
      {
//              if (bV0MCIsALambda && bV0MCIsPrimaryALambda) // well reconstructed candidates
        if(bV0MCIsALambda && bV0MCIsPrimaryDist) // well reconstructed candidates
        {
          fh2V0ALambdaPtMassMCRec[iCentIndex]->Fill(dPtV0Gen, dMassV0ALambda);
          Double_t valueEtaAL[3] = {dMassV0ALambda, dPtV0Gen, dEtaV0Gen};
          fh3V0ALambdaEtaPtMassMCRec[iCentIndex]->Fill(valueEtaAL);

          Double_t valueEtaDALNeg[6] = {0, particleMCDaughterNeg->Eta(), particleMCDaughterNeg->Pt(), dEtaV0Gen, dPtV0Gen, 0};
          fhnV0ALambdaInclDaughterEtaPtPtMCRec[iCentIndex]->Fill(valueEtaDALNeg);
          Double_t valueEtaDALPos[6] = {1, particleMCDaughterPos->Eta(), particleMCDaughterPos->Pt(), dEtaV0Gen, dPtV0Gen, 0};
          fhnV0ALambdaInclDaughterEtaPtPtMCRec[iCentIndex]->Fill(valueEtaDALPos);

          fh2V0ALambdaMCResolMPt[iCentIndex]->Fill(dMassV0ALambda - dMassPDGLambda, dPtV0);
          fh2V0ALambdaMCPtGenPtRec[iCentIndex]->Fill(dPtV0Gen, dPtV0);
          if(bIsInConeJet) // true V0 associated to a reconstructed candidate in jet
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
        if(bV0MCIsALambda && !bV0MCIsPrimaryDist && !bV0MCComesFromAXi) // not primary anti-Lambda
        {
          fh1V0ALambdaPtMCRecFalse[iCentIndex]->Fill(dPtV0Gen);
        }
      }
    }
    //===== End Association of reconstructed V0 candidates with MC particles =====
  }
  //===== End of V0 loop =====
  if(fDebug > 0) printf("%s %s::%s: %s\n", GetName(), ClassName(), __func__, "End of V0 loop");

  if(fbCorrelations && iNJetSel)
  {
    if(bPoolReady)
      fh1VtxZME[iCentIndex]->Fill(dZVtxME); // fill z_vtx if event was used for event mixing
    pool->UpdatePool(arrayMixedEventAdd); // update the pool with jets from this event
  }

  fh1V0CandPerEvent->Fill(iNV0CandTot);
  fh1V0CandPerEventCentK0s[iCentIndex]->Fill(iNV0CandK0s);
  fh1V0CandPerEventCentLambda[iCentIndex]->Fill(iNV0CandLambda);
  fh1V0CandPerEventCentALambda[iCentIndex]->Fill(iNV0CandALambda);

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
      if((iPdgCodeParticleMC == 3312) && (TMath::Abs(particleMC->Y()) < 0.5) && IsFromGoodGenerator(iPartMC))
      {
        fh1V0XiPtMCGen[iCentIndex]->Fill(particleMC->Pt());
      }
      else if((iPdgCodeParticleMC == -3312) && (TMath::Abs(particleMC->Y()) < 0.5) && IsFromGoodGenerator(iPartMC))
      {
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

      // V0 pseudorapidity cut
      if(fdCutEtaV0Max > 0.)
      {
        if(bPrintCuts) printf("Gen: Applying cut: V0 |eta|: < %g\n", fdCutEtaV0Max);
        if((TMath::Abs(dEtaV0Gen) > fdCutEtaV0Max))
          continue;
      }
      // V0 rapidity cut
      if(fdCutRapV0Max > 0.)
      {
        if(bPrintCuts) printf("Gen: Applying cut: V0 |y|: < %g\n", fdCutRapV0Max);
        if((TMath::Abs(dRapV0Gen) > fdCutRapV0Max))
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

      // Get the distance between the production point of the MC V0 particle and the primary vertex
      Double_t dx = dPrimVtxMCX - particleMC->Xv();
      Double_t dy = dPrimVtxMCY - particleMC->Yv();
      Double_t dz = dPrimVtxMCZ - particleMC->Zv();
      Double_t dDistPrimary = TMath::Sqrt(dx * dx + dy * dy + dz * dz);
      Bool_t bV0MCIsPrimaryDist = (dDistPrimary < dDistPrimaryMax); // Is close enough to be considered primary-like?

      // Select only primary-like MC V0 particles
      if(!bV0MCIsPrimaryDist)
        continue;

      // Select only particles from a specific generator
      if(!IsFromGoodGenerator(iPartMC))
        continue;

      // Check whether the MC V0 particle is in a MC jet
      AliAODJet* jetMC = 0;
      Bool_t bIsMCV0InJet = kFALSE;
      if(iNJetSel)
      {
        if(fDebug > 4) printf("%s %s::%s: %s\n", GetName(), ClassName(), __func__, Form("Searching for gen V0 in %d MC jets", iNJetSel));
        for(Int_t iJet = 0; iJet < iNJetSel; iJet++)
        {
          jetMC = (AliAODJet*)arrayJetSel->At(iJet); // load a jet in the list
          if(fDebug > 4) printf("%s %s::%s: %s\n", GetName(), ClassName(), __func__, Form("Checking if gen V0 in MC jet %d", iJet));
          if(IsParticleInCone(particleMC, jetMC, fdDistanceV0JetMax)) // If good jet in event, find out whether V0 is in that jet
          {
            if(fDebug > 4) printf("%s %s::%s: %s\n", GetName(), ClassName(), __func__, Form("gen V0 found in MC jet %d", iJet));
            bIsMCV0InJet = kTRUE;
            break;
          }
        }
      }

      // K0s
//          if (bV0MCIsK0s && bV0MCIsPrimary) // well reconstructed candidates
      if(bV0MCIsK0s) // well reconstructed candidates
      {
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
      if(bV0MCIsLambda) // well reconstructed candidates
      {
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
      if(bV0MCIsALambda) // well reconstructed candidates
      {
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

  arrayJetSel->Delete();
  delete arrayJetSel;
  arrayJetPerp->Delete();
  delete arrayJetPerp;
  if(jetRnd)
    delete jetRnd;
  jetRnd = 0;

  PostData(1, fOutputListStd);
  PostData(2, fOutputListQA);
  PostData(3, fOutputListCuts);
  PostData(4, fOutputListMC);

  if(fDebug > 0) printf("%s %s::%s: %s\n", GetName(), ClassName(), __func__, "End");

  return kFALSE; // Must be false to avoid calling PostData from AliAnalysisTaskEmcal. Otherwise, slot 1 is not stored.
}

void AliAnalysisTaskV0sInJetsEmcal::FillQAHistogramV0(AliAODVertex* vtx, const AliAODv0* vZero, Int_t iIndexHisto, Bool_t IsCandK0s, Bool_t IsCandLambda, Bool_t IsCandALambda, Bool_t IsInPeakK0s, Bool_t IsInPeakLambda, Bool_t IsInPeakALambda)
{
  if(!IsCandK0s && !IsCandLambda && !IsCandALambda)
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
//    if((nCrossedRowsTPC > (160. / (250. - 85.) * (255.*TMath::Abs(tan(track->Theta())) - 85.)) + 20.) && (track->Eta() < 0) && (track->Pt() > 0.15))
//      if (IsCandK0s)
//    {
    fh2QAV0EtaRows[iIndexHisto]->Fill(track->Eta(), nCrossedRowsTPC);
    fh2QAV0PtRows[iIndexHisto]->Fill(track->Pt(), nCrossedRowsTPC);
    fh2QAV0PhiRows[iIndexHisto]->Fill(track->Phi(), nCrossedRowsTPC);
    fh2QAV0NClRows[iIndexHisto]->Fill(findable, nCrossedRowsTPC);
    fh2QAV0EtaNCl[iIndexHisto]->Fill(track->Eta(), findable);
    fh2QAV0PtNCls[iIndexHisto]->Fill(track->Pt(), track->GetTPCNcls());
    fh2QAV0PtChi[iIndexHisto]->Fill(track->Pt(), track->Chi2perNDF());
//    }

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
//      fh2QAV0EtaPtK0sPeak[iIndexHisto]->Fill(trackNeg->Eta(),vZero->Pt());
//      fh2QAV0EtaPtK0sPeak[iIndexHisto]->Fill(trackPos->Eta(),vZero->Pt());
      fh2QAV0EtaPtK0sPeak[iIndexHisto]->Fill(vZero->Eta(), vZero->Pt());
      fh2QAV0PtPtK0sPeak[iIndexHisto]->Fill(trackNeg->Pt(), trackPos->Pt());
      fh2ArmPodK0s[iIndexHisto]->Fill(dAlpha, dPtArm);
      fh2QAV0PhiPtK0sPeak[iIndexHisto]->Fill(vZero->Phi(), vZero->Pt());
    }
    fh2QAV0EtaEtaK0s[iIndexHisto]->Fill(trackNeg->Eta(), trackPos->Eta());
    fh2QAV0PhiPhiK0s[iIndexHisto]->Fill(trackNeg->Phi(), trackPos->Phi());
    fh1QAV0RapK0s[iIndexHisto]->Fill(vZero->RapK0Short());
  }

  if(IsCandLambda)
  {
    if(IsInPeakLambda)
    {
//      fh2QAV0EtaPtLambdaPeak[iIndexHisto]->Fill(trackNeg->Eta(),vZero->Pt());
//      fh2QAV0EtaPtLambdaPeak[iIndexHisto]->Fill(trackPos->Eta(),vZero->Pt());
      fh2QAV0EtaPtLambdaPeak[iIndexHisto]->Fill(vZero->Eta(), vZero->Pt());
      fh2QAV0PtPtLambdaPeak[iIndexHisto]->Fill(trackNeg->Pt(), trackPos->Pt());
      fh2ArmPodLambda[iIndexHisto]->Fill(dAlpha, dPtArm);
      fh2QAV0PhiPtLambdaPeak[iIndexHisto]->Fill(vZero->Phi(), vZero->Pt());
    }
    fh2QAV0EtaEtaLambda[iIndexHisto]->Fill(trackNeg->Eta(), trackPos->Eta());
    fh2QAV0PhiPhiLambda[iIndexHisto]->Fill(trackNeg->Phi(), trackPos->Phi());
    fh1QAV0RapLambda[iIndexHisto]->Fill(vZero->RapLambda());
  }

  if(IsCandALambda)
  {
    if(IsInPeakALambda)
    {
//      fh2QAV0EtaPtALambdaPeak[iIndexHisto]->Fill(trackNeg->Eta(),vZero->Pt());
//      fh2QAV0EtaPtALambdaPeak[iIndexHisto]->Fill(trackPos->Eta(),vZero->Pt());
      fh2QAV0EtaPtALambdaPeak[iIndexHisto]->Fill(vZero->Eta(), vZero->Pt());
      fh2QAV0PtPtALambdaPeak[iIndexHisto]->Fill(trackNeg->Pt(), trackPos->Pt());
      fh2ArmPodALambda[iIndexHisto]->Fill(dAlpha, dPtArm);
      fh2QAV0PhiPtALambdaPeak[iIndexHisto]->Fill(vZero->Phi(), vZero->Pt());
    }
    fh2QAV0EtaEtaALambda[iIndexHisto]->Fill(trackNeg->Eta(), trackPos->Eta());
    fh2QAV0PhiPhiALambda[iIndexHisto]->Fill(trackNeg->Phi(), trackPos->Phi());
    fh1QAV0RapALambda[iIndexHisto]->Fill(vZero->RapLambda());
  }

  fh2ArmPod[iIndexHisto]->Fill(dAlpha, dPtArm);
}

void AliAnalysisTaskV0sInJetsEmcal::FillCandidates(Double_t mK, Double_t mL, Double_t mAL, Bool_t isK, Bool_t isL, Bool_t isAL, Int_t iCut/*cut index*/, Int_t iCent/*cent index*/)
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

Bool_t AliAnalysisTaskV0sInJetsEmcal::IsParticleInCone(const AliVParticle* part1, const AliVParticle* part2, Double_t dRMax) const
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

Bool_t AliAnalysisTaskV0sInJetsEmcal::OverlapWithJets(const TClonesArray* array, const AliVParticle* part, Double_t dDistance) const
{
// decides whether a cone overlaps with other jets
  if(!part)
  {
    AliError("No particle!");
    return kFALSE;
  }
  if(!array)
  {
    AliError("No jet array!");
    return kFALSE;
  }
  Int_t iNJets = array->GetEntriesFast();
  if(iNJets <= 0)
  {
    if(fDebug > 0) printf("%s %s::%s: %s\n", GetName(), ClassName(), __func__, "Warning: No jets");
    return kFALSE;
  }
  AliVParticle* jet = 0;
  for(Int_t iJet = 0; iJet < iNJets; iJet++)
  {
    jet = (AliVParticle*)array->At(iJet);
    if(!jet)
    {
      AliError(Form("Failed to load jet %d/%d!", iJet, iNJets));
      continue;
    }
    if(IsParticleInCone(part, jet, dDistance))
      return kTRUE;
  }
  return kFALSE;
}

AliAODJet* AliAnalysisTaskV0sInJetsEmcal::GetRandomCone(const TClonesArray* array, Double_t dEtaConeMax, Double_t dDistance) const
{
// generate a random cone which does not overlap with selected jets
  if(fDebug > 3) printf("%s %s::%s: %s\n", GetName(), ClassName(), __func__, "Generating random cone...");
  TLorentzVector vecCone;
  AliAODJet* part = 0;
  Double_t dEta, dPhi;
  Int_t iNTrialsMax = 10;
  Bool_t bStatus = kFALSE;
  for(Int_t iTry = 0; iTry < iNTrialsMax; iTry++)
  {
    if(fDebug > 4) printf("%s %s::%s: %s\n", GetName(), ClassName(), __func__, Form("Try %d", iTry));
    dEta = dEtaConeMax * (2 * fRandom->Rndm() - 1.); // random eta in [-dEtaConeMax,+dEtaConeMax]
    dPhi = TMath::TwoPi() * fRandom->Rndm(); // random phi in [0,2*Pi]
    vecCone.SetPtEtaPhiM(1., dEta, dPhi, 0.);
    part = new AliAODJet(vecCone);
    if(!OverlapWithJets(array, part, dDistance))
    {
      bStatus = kTRUE;
      if(fDebug > 1) printf("%s %s::%s: %s\n", GetName(), ClassName(), __func__, "Random cone successfully generated");
      break;
    }
    else
      delete part;
  }
  if(!bStatus)
  {
    if(fDebug > 0) printf("%s %s::%s: %s\n", GetName(), ClassName(), __func__, "Failed to find a random cone");
    part = 0;
  }
  return part;
}

AliEmcalJet* AliAnalysisTaskV0sInJetsEmcal::GetMedianCluster(AliJetContainer* cont, Double_t dEtaConeMax) const
{
// sort kt clusters by pT/area and return the middle one, based on code in AliAnalysisTaskJetChem
  if(!cont)
  {
    AliError("No jet container!");
    return NULL;
  }
  Int_t iNClTot = cont->GetNJets(); // number of all clusters in the array
  Int_t iNCl = 0; // number of accepted clusters

  // get list of densities
  std::vector<std::vector<Double_t> > vecListClusters; // vector that contains pairs [ index, density ]
  if(fDebug > 3) printf("%s %s::%s: %s\n", GetName(), ClassName(), __func__, Form("Loop over %d clusters.", iNClTot));
  for(Int_t ij = 0; ij < iNClTot; ij++)
  {
    AliEmcalJet* clusterBg = (AliEmcalJet*)(cont->GetAcceptJet(ij));
    if(!clusterBg)
      continue;
    if(fDebug > 4) printf("%s %s::%s: %s\n", GetName(), ClassName(), __func__, Form("Cluster %d/%d used as accepted cluster %d.", ij, iNClTot, int(vecListClusters.size())));
    Double_t dPtBg = clusterBg->Pt();
    Double_t dAreaBg = clusterBg->Area();
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
    if(fDebug > 0) printf("%s %s::%s: %s\n", GetName(), ClassName(), __func__, "Warning: Too little clusters");
    return NULL;
  }

//  printf("%s %s::%s: %s\n", GetName(), ClassName(), __func__, "Original lists:");
//  for(Int_t i = 0; i < iNCl; i++)
//    printf("%g %g\n", (vecListClusters[i])[0], (vecListClusters[i])[1]);

  // sort list of indeces by density in descending order
  std::sort(vecListClusters.begin(), vecListClusters.end(), CompareClusters);

//  printf("%s %s::%s: %s\n", GetName(), ClassName(), __func__, "Sorted lists:");
//  for(Int_t i = 0; i < iNCl; i++)
//    printf("%g %g\n", (vecListClusters[i])[0], (vecListClusters[i])[1]);

  // get median cluster with median density
  AliEmcalJet* clusterMed = 0;
  Int_t iIndex = 0; // index of the median cluster in the sorted list
  Int_t iIndexMed = 0; // index of the median cluster in the original array
  if(TMath::Odd(iNCl)) // odd number of clusters
  {
    iIndex = (Int_t)(0.5 * (iNCl + 1)); // = (n - skip + 1)/2 + 1, skip = 2
//    printf("%s %s::%s: %s\n", GetName(), ClassName(), __func__, Form("Odd, median index = %d/%d", iIndex, iNCl));
  }
  else // even number: picking randomly one of the two closest to median
  {
    Int_t iIndex1 = (Int_t)(0.5 * iNCl); // = (n - skip)/2 + 1, skip = 2
    Int_t iIndex2 = (Int_t)(0.5 * iNCl + 1); // = (n - skip)/2 + 1 + 1, skip = 2
    iIndex = ((fRandom->Rndm() > 0.5) ? iIndex1 : iIndex2);
//    printf("%s %s::%s: %s\n", GetName(), ClassName(), __func__, Form("Even, median index = %d or %d -> %d/%d", iIndex1, iIndex2, iIndex, iNCl));
  }
  iIndexMed = Int_t((vecListClusters[iIndex])[0]);

  if(fDebug > 1) printf("%s %s::%s: %s\n", GetName(), ClassName(), __func__, Form("Getting median cluster %d/%d ok, rho = %g", iIndexMed, iNClTot, (vecListClusters[iIndex])[1]));
  clusterMed = (AliEmcalJet*)(cont->GetAcceptJet(iIndexMed));

  if(TMath::Abs(clusterMed->Eta()) > dEtaConeMax)
  {
    if(fDebug > 0) printf("%s %s::%s: %s\n", GetName(), ClassName(), __func__, Form("Chosen median cluster is out of window |eta| < %g, eta = %g", dEtaConeMax, clusterMed->Eta()));
    return NULL;
  }

  return clusterMed;
}

Double_t AliAnalysisTaskV0sInJetsEmcal::AreaCircSegment(Double_t dRadius, Double_t dDistance) const
{
// calculate area of a circular segment defined by the circle radius and the (oriented) distance between the secant line and the circle centre
  Double_t dEpsilon = 1e-2;
  Double_t dR = dRadius;
  Double_t dD = dDistance;
  if(TMath::Abs(dR) < dEpsilon)
  {
    AliError(Form("Too small radius: %g < %g!", dR, dEpsilon));
    return 0.;
  }
  if(dD >= dR)
    return 0.;
  if(dD <= -dR)
    return TMath::Pi() * dR * dR;
  return dR * dR * TMath::ACos(dD / dR) - dD * TMath::Sqrt(dR * dR - dD * dD);
}

Double_t AliAnalysisTaskV0sInJetsEmcal::GetD(const AliVParticle* part1, const AliVParticle* part2) const
{
// return D between two particles
  if(!part1 || !part2)
    return -1;

  TVector3 vecMom2(part2->Px(), part2->Py(), part2->Pz());
  TVector3 vecMom1(part1->Px(), part1->Py(), part1->Pz());
  Double_t dR = vecMom2.DeltaR(vecMom1); // = sqrt(dEta*dEta+dPhi*dPhi)
  return dR;
}

Bool_t AliAnalysisTaskV0sInJetsEmcal::IsSelectedForJets(AliAODEvent* fAOD, Double_t dVtxZCut, Double_t dVtxR2Cut, Double_t dCentCutLo, Double_t dCentCutUp, Double_t dDeltaZMax)
{
// event selection
  if(!fbIsPbPb)
  {
    if(fAOD->IsPileupFromSPD())
    {
      if(fDebug > 0) printf("%s %s::%s: %s\n", GetName(), ClassName(), __func__, "SPD pile-up");
      return kFALSE;
    }
    fh1EventCounterCut->Fill(2); // not pile-up from SPD
  }
  AliAODVertex* vertex = fAOD->GetPrimaryVertex();
  if(!vertex)
  {
    if(fDebug > 0) printf("%s %s::%s: %s\n", GetName(), ClassName(), __func__, "No vertex");
    return kFALSE;
  }
  Int_t iNContribMin = 3;
  if(vertex->GetNContributors() < iNContribMin)
  {
    if(fDebug > 0) printf("%s %s::%s: %s\n", GetName(), ClassName(), __func__, Form("Not enough contributors, %d", vertex->GetNContributors()));
    return kFALSE;
  }
  fh1EventCounterCut->Fill(3); // enough contributors
  TString vtxTitle(vertex->GetTitle());
  if(vtxTitle.Contains("TPCVertex"))
  {
    if(fDebug > 0) printf("%s %s::%s: %s\n", GetName(), ClassName(), __func__, "TPC vertex");
    return kFALSE;
  }
  fh1EventCounterCut->Fill(4); // not TPC vertex only
  Double_t zVertex = vertex->GetZ();
  if(TMath::Abs(zVertex) > dVtxZCut)
  {
    if(fDebug > 0) printf("%s %s::%s: %s\n", GetName(), ClassName(), __func__, Form("Cut on z, %g", zVertex));
    return kFALSE;
  }
  fh1EventCounterCut->Fill(5); // PV z coordinate within range
  if(dDeltaZMax > 0.) // cut on |delta z| between SPD vertex and nominal primary vertex
  {
    AliAODVertex* vertexSPD = fAOD->GetPrimaryVertexSPD();
    if(!vertexSPD)
    {
      if(fDebug > 0) printf("%s %s::%s: %s\n", GetName(), ClassName(), __func__, "No SPD vertex");
      return kFALSE;
    }
    Double_t zVertexSPD = vertexSPD->GetZ();
    if(TMath::Abs(zVertex - zVertexSPD) > dDeltaZMax)
    {
      if(fDebug > 0) printf("%s %s::%s: %s\n", GetName(), ClassName(), __func__, Form("Cut on Delta z = %g - %g = %g", zVertex, zVertexSPD, zVertex - zVertexSPD));
      return kFALSE;
    }
    fh1EventCounterCut->Fill(6); // delta z within range
  }
  Double_t xVertex = vertex->GetX();
  Double_t yVertex = vertex->GetY();
  Double_t radiusSq = yVertex * yVertex + xVertex * xVertex;
  if(radiusSq > dVtxR2Cut)
  {
    if(fDebug > 0) printf("%s %s::%s: %s\n", GetName(), ClassName(), __func__, Form("Cut on r, %g", radiusSq));
    return kFALSE;
  }
  fh1EventCounterCut->Fill(7); // radius within range
  if(fbIsPbPb)
  {
    fdCentrality = ((AliVAODHeader*)fAOD->GetHeader())->GetCentralityP()->GetCentralityPercentile("V0M");
    if(fdCentrality < 0)
    {
      if(fDebug > 0) printf("%s %s::%s: %s\n", GetName(), ClassName(), __func__, "Negative centrality");
      return kFALSE;
    }
    if((dCentCutUp < 0) || (dCentCutLo < 0) || (dCentCutUp > 100) || (dCentCutLo > 100) || (dCentCutLo > dCentCutUp))
    {
      if(fDebug > 0) printf("%s %s::%s: %s\n", GetName(), ClassName(), __func__, "Wrong centrality limits");
      return kFALSE;
    }
    if((fdCentrality < dCentCutLo) || (fdCentrality > dCentCutUp))
    {
      if(fDebug > 0) printf("%s %s::%s: %s\n", GetName(), ClassName(), __func__, Form("Centrality cut, %g", fdCentrality));
      return kFALSE;
    }
  }
  return kTRUE;
}

Int_t AliAnalysisTaskV0sInJetsEmcal::GetCentralityBinIndex(Double_t centrality)
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

Int_t AliAnalysisTaskV0sInJetsEmcal::GetCentralityBinEdge(Int_t index)
{
// returns the upper edge of the centrality bin corresponding to the provided value of index
  if(index < 0 || index >= fgkiNBinsCent)
    return -1;
  return fgkiCentBinRanges[index];
}

TString AliAnalysisTaskV0sInJetsEmcal::GetCentBinLabel(Int_t index)
{
// get string with centrality range for given bin
  TString lowerEdge = ((index == 0) ? "0" : Form("%d", GetCentralityBinEdge(index - 1)));
  TString upperEdge = Form("%d", GetCentralityBinEdge(index));
  return Form("%s-%s %%", lowerEdge.Data(), upperEdge.Data());
}

Double_t AliAnalysisTaskV0sInJetsEmcal::MassPeakSigmaOld(Double_t pt, Int_t particle)
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

bool AliAnalysisTaskV0sInJetsEmcal::CompareClusters(const std::vector<Double_t> cluster1, const std::vector<Double_t> cluster2)
{
  return (cluster1[1] > cluster2[1]);
}

Bool_t AliAnalysisTaskV0sInJetsEmcal::IsFromGoodGenerator(Int_t index)
{
  if(!fEventMC)
  {
    AliError("No MC event!");
    return kFALSE;
  }
  if(fsGeneratorName.Length())
  {
    TString sGenName = "";
    Bool_t bCocktailOK = fEventMC->GetCocktailGenerator(index, sGenName);
    if(!bCocktailOK || !sGenName.Contains(fsGeneratorName.Data()))
      return kFALSE;
  }
  return kTRUE;
}
