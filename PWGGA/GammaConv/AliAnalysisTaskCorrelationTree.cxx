/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *									  *
 * Authors: Svein Lindal, Daniel Lohner					  *
 * Version 1.0								  *
 *									  *
 * Permission to use, copy, modify and distribute this software and its	  *
 * documentation strictly for non-commercial purposes is hereby granted	  *
 * without fee, provided that the above copyright notice appears in all	  *
 * copies and that both the copyright notice and this permission notice	  *
 * appear in the supporting documentation. The authors make no claims	  *
 * about the suitability of this software for any purpose. It is	  *
 * provided "as is" without express or implied warranty.		  *
 **************************************************************************/

////////////////////////////////////////////////
//---------------------------------------------
// QA Task for V0 Reader V1
//---------------------------------------------
////////////////////////////////////////////////

#include "AliAnalysisTaskCorrelationTree.h"
#include "TChain.h"
#include "TRandom.h"
#include "AliAnalysisManager.h"
#include "TVectorF.h"
#include "AliPIDResponse.h"
#include "TFile.h"
#include "AliESDtrackCuts.h"
#include "AliAODMCParticle.h"
#include "AliAODMCHeader.h"
#include "AliAODEvent.h"

class iostream;

using namespace std;

ClassImp(AliAnalysisTaskCorrelationTree)

    //________________________________________________________________________
    AliAnalysisTaskCorrelationTree::AliAnalysisTaskCorrelationTree() : AliAnalysisTaskSE(),
                                                                       fV0Reader(NULL),
                                                                       fV0ReaderName("V0ReaderV1"),
                                                                       fConversionCuts(NULL),
                                                                       fEventCuts(NULL),
                                                                       fInputEvent(NULL),
                                                                       fMCEvent(NULL),
                                                                       fAODMCTrackArray(NULL),
                                                                       farrClustersProcess(NULL),
                                                                       fCorrTaskSetting(""),
                                                                       fClusterCuts(NULL),

                                                                       fAnalysisTree(NULL),
                                                                       fIsHeavyIon(false),
                                                                       fOutputList(NULL),
                                                                       fPidResponse(0),
                                                                       fHistoNEvents(NULL),

                                                                       fBuffer_NContributors(0),
                                                                       fBuffer_RunNumber(0),
                                                                       fBuffer_VertexZ(0),
                                                                       fBuffer_NEventTriggers(0),
                                                                       fBuffer_EventTrigger(0),
                                                                       fBuffer_NElectronCandidates(0),
                                                                       fBuffer_ElectronCandidate_E(0),
                                                                       fBuffer_ElectronCandidate_Px(0),
                                                                       fBuffer_ElectronCandidate_Py(0),
                                                                       fBuffer_ElectronCandidate_Pz(0),
                                                                       fBuffer_ElectronCandidate_PropEta(0),
                                                                       fBuffer_ElectronCandidate_PropPhi(0),
                                                                       fBuffer_ElectronCandidate_Charge(0),
                                                                       fBuffer_ElectronCandidate_NSigmaElecTPC(0),
                                                                       fBuffer_ElectronCandidate_NSigmaElecTOF(0),
                                                                       fBuffer_ElectronCandidate_NSigmaPionTPC(0),
                                                                       fBuffer_ElectronCandidate_NSigmaPionTOF(0),
                                                                       fBuffer_ElectronCandidate_NSigmaKaonTPC(0),
                                                                       fBuffer_ElectronCandidate_NSigmaProtonTPC(0),
                                                                       fBuffer_ElectronCandidate_MC_E(0),
                                                                       fBuffer_ElectronCandidate_MC_Px(0),
                                                                       fBuffer_ElectronCandidate_MC_Py(0),
                                                                       fBuffer_ElectronCandidate_MC_Pz(0),
                                                                       fBuffer_ElectronCandidate_MC_PDG(0),
                                                                       fBuffer_ElectronCandidate_MC_Mother_E(0),
                                                                       fBuffer_ElectronCandidate_MC_Mother_Px(0),
                                                                       fBuffer_ElectronCandidate_MC_Mother_Py(0),
                                                                       fBuffer_ElectronCandidate_MC_Mother_Pz(0),
                                                                       fBuffer_ElectronCandidate_MC_Mother_PDG(0),
                                                                       fBuffer_ElectronCandidate_MC_GrandMother_PDG(0),
                                                                       fBuffer_ElectronCandidate_MC_GrandMother_Pt(0),
                                                                       fBuffer_ElectronCandidate_MC_GGrandMother_PDG(0),
                                                                       fBuffer_ElectronCandidate_MC_GGrandMother_Pt(0),

                                                                       fBuffer_NMuonCandidates(0),
                                                                       fBuffer_MuonCandidate_E(0),
                                                                       fBuffer_MuonCandidate_Px(0),
                                                                       fBuffer_MuonCandidate_Py(0),
                                                                       fBuffer_MuonCandidate_Pz(0),
                                                                       fBuffer_MuonCandidate_RAbsEnd(0),
                                                                       fBuffer_MuonCandidate_Chi2NDF(0),
                                                                       fBuffer_MuonCandidate_DCA(0),
                                                                       fBuffer_MuonCandidate_Charge(0),
                                                                       fBuffer_MuonCandidate_MatchTrigger(0),
                                                                       fBuffer_MuonCandidate_MC_E(0),
                                                                       fBuffer_MuonCandidate_MC_Px(0),
                                                                       fBuffer_MuonCandidate_MC_Py(0),
                                                                       fBuffer_MuonCandidate_MC_Pz(0),
                                                                       fBuffer_MuonCandidate_MC_PDG(0),
                                                                       fBuffer_MuonCandidate_MC_Mother_E(0),
                                                                       fBuffer_MuonCandidate_MC_Mother_Px(0),
                                                                       fBuffer_MuonCandidate_MC_Mother_Py(0),
                                                                       fBuffer_MuonCandidate_MC_Mother_Pz(0),
                                                                       fBuffer_MuonCandidate_MC_Mother_PDG(0),
                                                                       fBuffer_MuonCandidate_MC_GrandMother_PDG(0),
                                                                       fBuffer_MuonCandidate_MC_GrandMother_Pt(0),
                                                                       fBuffer_MuonCandidate_MC_GGrandMother_PDG(0),
                                                                       fBuffer_MuonCandidate_MC_GGrandMother_Pt(0),

                                                                       fBuffer_NClusterCandidates(0),
                                                                       fBuffer_ClusterCandidate_E(0),
                                                                       fBuffer_ClusterCandidate_Px(0),
                                                                       fBuffer_ClusterCandidate_Py(0),
                                                                       fBuffer_ClusterCandidate_Pz(0),
                                                                       fBuffer_ClusterCandidate_Eta(0),
                                                                       fBuffer_ClusterCandidate_Phi(0),
                                                                       fBuffer_ClusterCandidate_M02(0),
                                                                       fBuffer_ClusterCandidate_MC_E(0),
                                                                       fBuffer_ClusterCandidate_MC_Px(0),
                                                                       fBuffer_ClusterCandidate_MC_Py(0),
                                                                       fBuffer_ClusterCandidate_MC_Pz(0),
                                                                       fBuffer_ClusterCandidate_MC_PDG(0),
                                                                       fBuffer_ClusterCandidate_MC_Mother_E(0),
                                                                       fBuffer_ClusterCandidate_MC_Mother_Px(0),
                                                                       fBuffer_ClusterCandidate_MC_Mother_Py(0),
                                                                       fBuffer_ClusterCandidate_MC_Mother_Pz(0),
                                                                       fBuffer_ClusterCandidate_MC_Mother_PDG(0),
                                                                       fBuffer_ClusterCandidate_MC_GrandMother_PDG(0),
                                                                       fIsMC(false)
{
  fBuffer_EventTrigger = new Int_t[kMaxTriggers];
  fBuffer_ElectronCandidate_E = new Float_t[kMaxTracks];
  fBuffer_ElectronCandidate_Px = new Float_t[kMaxTracks];
  fBuffer_ElectronCandidate_Py = new Float_t[kMaxTracks];
  fBuffer_ElectronCandidate_Pz = new Float_t[kMaxTracks];
  fBuffer_ElectronCandidate_PropEta = new Float_t[kMaxTracks];
  fBuffer_ElectronCandidate_PropPhi = new Float_t[kMaxTracks];
  fBuffer_ElectronCandidate_Charge = new Float_t[kMaxTracks];
  fBuffer_ElectronCandidate_NSigmaElecTPC = new Float_t[kMaxTracks];
  fBuffer_ElectronCandidate_NSigmaElecTOF = new Float_t[kMaxTracks];
  fBuffer_ElectronCandidate_NSigmaKaonTPC = new Float_t[kMaxTracks];
  fBuffer_ElectronCandidate_NSigmaPionTPC = new Float_t[kMaxTracks];
  fBuffer_ElectronCandidate_NSigmaPionTOF = new Float_t[kMaxTracks];
  fBuffer_ElectronCandidate_NSigmaProtonTPC = new Float_t[kMaxTracks];
  fBuffer_ElectronCandidate_MC_E = new Float_t[kMaxTracks];
  fBuffer_ElectronCandidate_MC_Px = new Float_t[kMaxTracks];
  fBuffer_ElectronCandidate_MC_Py = new Float_t[kMaxTracks];
  fBuffer_ElectronCandidate_MC_Pz = new Float_t[kMaxTracks];
  fBuffer_ElectronCandidate_MC_PDG = new Int_t[kMaxTracks];
  fBuffer_ElectronCandidate_MC_Mother_E = new Float_t[kMaxTracks];
  fBuffer_ElectronCandidate_MC_Mother_Px = new Float_t[kMaxTracks];
  fBuffer_ElectronCandidate_MC_Mother_Py = new Float_t[kMaxTracks];
  fBuffer_ElectronCandidate_MC_Mother_Pz = new Float_t[kMaxTracks];
  fBuffer_ElectronCandidate_MC_Mother_PDG = new Int_t[kMaxTracks];
  fBuffer_ElectronCandidate_MC_GrandMother_PDG = new Int_t[kMaxTracks];
  fBuffer_ElectronCandidate_MC_GrandMother_Pt = new Float_t[kMaxTracks];
  fBuffer_ElectronCandidate_MC_GGrandMother_PDG = new Int_t[kMaxTracks];
  fBuffer_ElectronCandidate_MC_GGrandMother_Pt = new Float_t[kMaxTracks];

  fBuffer_MuonCandidate_E = new Float_t[kMaxTracks];
  fBuffer_MuonCandidate_Px = new Float_t[kMaxTracks];
  fBuffer_MuonCandidate_Py = new Float_t[kMaxTracks];
  fBuffer_MuonCandidate_Pz = new Float_t[kMaxTracks];
  fBuffer_MuonCandidate_RAbsEnd = new Float_t[kMaxTracks];
  fBuffer_MuonCandidate_Chi2NDF = new Float_t[kMaxTracks];
  fBuffer_MuonCandidate_DCA = new Float_t[kMaxTracks];
  fBuffer_MuonCandidate_Charge = new Int_t[kMaxTracks];
  fBuffer_MuonCandidate_MatchTrigger = new Int_t[kMaxTracks];
  fBuffer_MuonCandidate_MC_E = new Float_t[kMaxTracks];
  fBuffer_MuonCandidate_MC_Px = new Float_t[kMaxTracks];
  fBuffer_MuonCandidate_MC_Py = new Float_t[kMaxTracks];
  fBuffer_MuonCandidate_MC_Pz = new Float_t[kMaxTracks];
  fBuffer_MuonCandidate_MC_PDG = new Int_t[kMaxTracks];
  fBuffer_MuonCandidate_MC_Mother_E = new Float_t[kMaxTracks];
  fBuffer_MuonCandidate_MC_Mother_Px = new Float_t[kMaxTracks];
  fBuffer_MuonCandidate_MC_Mother_Py = new Float_t[kMaxTracks];
  fBuffer_MuonCandidate_MC_Mother_Pz = new Float_t[kMaxTracks];
  fBuffer_MuonCandidate_MC_Mother_PDG = new Int_t[kMaxTracks];
  fBuffer_MuonCandidate_MC_GrandMother_PDG = new Int_t[kMaxTracks];
  fBuffer_MuonCandidate_MC_GrandMother_Pt = new Float_t[kMaxTracks];
  fBuffer_MuonCandidate_MC_GGrandMother_PDG = new Int_t[kMaxTracks];
  fBuffer_MuonCandidate_MC_GGrandMother_Pt = new Float_t[kMaxTracks];

  fBuffer_ClusterCandidate_E = new Float_t[kMaxTracks];
  fBuffer_ClusterCandidate_Px = new Float_t[kMaxTracks];
  fBuffer_ClusterCandidate_Py = new Float_t[kMaxTracks];
  fBuffer_ClusterCandidate_Pz = new Float_t[kMaxTracks];
  fBuffer_ClusterCandidate_Eta = new Float_t[kMaxTracks];
  fBuffer_ClusterCandidate_Phi = new Float_t[kMaxTracks];
  fBuffer_ClusterCandidate_M02 = new Float_t[kMaxTracks];
  fBuffer_ClusterCandidate_MC_E = new Float_t[kMaxTracks];
  fBuffer_ClusterCandidate_MC_Px = new Float_t[kMaxTracks];
  fBuffer_ClusterCandidate_MC_Py = new Float_t[kMaxTracks];
  fBuffer_ClusterCandidate_MC_Pz = new Float_t[kMaxTracks];
  fBuffer_ClusterCandidate_MC_PDG = new Int_t[kMaxTracks];
  fBuffer_ClusterCandidate_MC_Mother_E = new Float_t[kMaxTracks];
  fBuffer_ClusterCandidate_MC_Mother_Px = new Float_t[kMaxTracks];
  fBuffer_ClusterCandidate_MC_Mother_Py = new Float_t[kMaxTracks];
  fBuffer_ClusterCandidate_MC_Mother_Pz = new Float_t[kMaxTracks];
  fBuffer_ClusterCandidate_MC_Mother_PDG = new Int_t[kMaxTracks];
  fBuffer_ClusterCandidate_MC_GrandMother_PDG = new Int_t[kMaxTracks];
}

AliAnalysisTaskCorrelationTree::AliAnalysisTaskCorrelationTree(const char *name) : AliAnalysisTaskSE(name),
                                                                                   fV0Reader(NULL),
                                                                                   fV0ReaderName("V0ReaderV1"),
                                                                                   fConversionCuts(NULL),
                                                                                   fEventCuts(NULL),
                                                                                   fInputEvent(NULL),
                                                                                   fMCEvent(NULL),
                                                                                   fAODMCTrackArray(NULL),
                                                                                   farrClustersProcess(NULL),
                                                                                   fCorrTaskSetting(""),
                                                                                   fClusterCuts(NULL),

                                                                                   fAnalysisTree(NULL),
                                                                                   fIsHeavyIon(false),
                                                                                   fOutputList(NULL),
                                                                                   fPidResponse(0),
                                                                                   fHistoNEvents(NULL),

                                                                                   fBuffer_NContributors(0),
                                                                                   fBuffer_RunNumber(0),
                                                                                   fBuffer_VertexZ(0),
                                                                                   fBuffer_NEventTriggers(0),
                                                                                   fBuffer_EventTrigger(0),
                                                                                   fBuffer_NElectronCandidates(0),
                                                                                   fBuffer_ElectronCandidate_E(0),
                                                                                   fBuffer_ElectronCandidate_Px(0),
                                                                                   fBuffer_ElectronCandidate_Py(0),
                                                                                   fBuffer_ElectronCandidate_Pz(0),
                                                                                   fBuffer_ElectronCandidate_PropEta(0),
                                                                                   fBuffer_ElectronCandidate_PropPhi(0),
                                                                                   fBuffer_ElectronCandidate_Charge(0),
                                                                                   fBuffer_ElectronCandidate_NSigmaElecTPC(0),
                                                                                   fBuffer_ElectronCandidate_NSigmaElecTOF(0),
                                                                                   fBuffer_ElectronCandidate_NSigmaPionTPC(0),
                                                                                   fBuffer_ElectronCandidate_NSigmaPionTOF(0),
                                                                                   fBuffer_ElectronCandidate_NSigmaKaonTPC(0),
                                                                                   fBuffer_ElectronCandidate_NSigmaProtonTPC(0),
                                                                                   fBuffer_ElectronCandidate_MC_E(0),
                                                                                   fBuffer_ElectronCandidate_MC_Px(0),
                                                                                   fBuffer_ElectronCandidate_MC_Py(0),
                                                                                   fBuffer_ElectronCandidate_MC_Pz(0),
                                                                                   fBuffer_ElectronCandidate_MC_PDG(0),
                                                                                   fBuffer_ElectronCandidate_MC_Mother_E(0),
                                                                                   fBuffer_ElectronCandidate_MC_Mother_Px(0),
                                                                                   fBuffer_ElectronCandidate_MC_Mother_Py(0),
                                                                                   fBuffer_ElectronCandidate_MC_Mother_Pz(0),
                                                                                   fBuffer_ElectronCandidate_MC_Mother_PDG(0),
                                                                                   fBuffer_ElectronCandidate_MC_GrandMother_PDG(0),
                                                                                   fBuffer_ElectronCandidate_MC_GrandMother_Pt(0),
                                                                                   fBuffer_ElectronCandidate_MC_GGrandMother_PDG(0),
                                                                                   fBuffer_ElectronCandidate_MC_GGrandMother_Pt(0),

                                                                                   fBuffer_NMuonCandidates(0),
                                                                                   fBuffer_MuonCandidate_E(0),
                                                                                   fBuffer_MuonCandidate_Px(0),
                                                                                   fBuffer_MuonCandidate_Py(0),
                                                                                   fBuffer_MuonCandidate_Pz(0),
                                                                                   fBuffer_MuonCandidate_RAbsEnd(0),
                                                                                   fBuffer_MuonCandidate_Chi2NDF(0),
                                                                                   fBuffer_MuonCandidate_DCA(0),
                                                                                   fBuffer_MuonCandidate_Charge(0),
                                                                                   fBuffer_MuonCandidate_MatchTrigger(0),
                                                                                   fBuffer_MuonCandidate_MC_E(0),
                                                                                   fBuffer_MuonCandidate_MC_Px(0),
                                                                                   fBuffer_MuonCandidate_MC_Py(0),
                                                                                   fBuffer_MuonCandidate_MC_Pz(0),
                                                                                   fBuffer_MuonCandidate_MC_PDG(0),
                                                                                   fBuffer_MuonCandidate_MC_Mother_E(0),
                                                                                   fBuffer_MuonCandidate_MC_Mother_Px(0),
                                                                                   fBuffer_MuonCandidate_MC_Mother_Py(0),
                                                                                   fBuffer_MuonCandidate_MC_Mother_Pz(0),
                                                                                   fBuffer_MuonCandidate_MC_Mother_PDG(0),
                                                                                   fBuffer_MuonCandidate_MC_GrandMother_PDG(0),
                                                                                   fBuffer_MuonCandidate_MC_GrandMother_Pt(0),
                                                                                   fBuffer_MuonCandidate_MC_GGrandMother_PDG(0),
                                                                                   fBuffer_MuonCandidate_MC_GGrandMother_Pt(0),

                                                                                   fBuffer_NClusterCandidates(0),
                                                                                   fBuffer_ClusterCandidate_E(0),
                                                                                   fBuffer_ClusterCandidate_Px(0),
                                                                                   fBuffer_ClusterCandidate_Py(0),
                                                                                   fBuffer_ClusterCandidate_Pz(0),
                                                                                   fBuffer_ClusterCandidate_Eta(0),
                                                                                   fBuffer_ClusterCandidate_Phi(0),
                                                                                   fBuffer_ClusterCandidate_M02(0),
                                                                                   fBuffer_ClusterCandidate_MC_E(0),
                                                                                   fBuffer_ClusterCandidate_MC_Px(0),
                                                                                   fBuffer_ClusterCandidate_MC_Py(0),
                                                                                   fBuffer_ClusterCandidate_MC_Pz(0),
                                                                                   fBuffer_ClusterCandidate_MC_PDG(0),
                                                                                   fBuffer_ClusterCandidate_MC_Mother_E(0),
                                                                                   fBuffer_ClusterCandidate_MC_Mother_Px(0),
                                                                                   fBuffer_ClusterCandidate_MC_Mother_Py(0),
                                                                                   fBuffer_ClusterCandidate_MC_Mother_Pz(0),
                                                                                   fBuffer_ClusterCandidate_MC_Mother_PDG(0),
                                                                                   fBuffer_ClusterCandidate_MC_GrandMother_PDG(0),
                                                                                   fIsMC(false)
{
  // Default constructor
  fBuffer_EventTrigger = new Int_t[kMaxTriggers];
  fBuffer_ElectronCandidate_E = new Float_t[kMaxTracks];
  fBuffer_ElectronCandidate_Px = new Float_t[kMaxTracks];
  fBuffer_ElectronCandidate_Py = new Float_t[kMaxTracks];
  fBuffer_ElectronCandidate_Pz = new Float_t[kMaxTracks];
  fBuffer_ElectronCandidate_PropEta = new Float_t[kMaxTracks];
  fBuffer_ElectronCandidate_PropPhi = new Float_t[kMaxTracks];
  fBuffer_ElectronCandidate_Charge = new Float_t[kMaxTracks];
  fBuffer_ElectronCandidate_NSigmaElecTPC = new Float_t[kMaxTracks];
  fBuffer_ElectronCandidate_NSigmaElecTOF = new Float_t[kMaxTracks];
  fBuffer_ElectronCandidate_NSigmaKaonTPC = new Float_t[kMaxTracks];
  fBuffer_ElectronCandidate_NSigmaPionTPC = new Float_t[kMaxTracks];
  fBuffer_ElectronCandidate_NSigmaPionTOF = new Float_t[kMaxTracks];
  fBuffer_ElectronCandidate_NSigmaProtonTPC = new Float_t[kMaxTracks];

  fBuffer_ElectronCandidate_MC_E = new Float_t[kMaxTracks];
  fBuffer_ElectronCandidate_MC_Px = new Float_t[kMaxTracks];
  fBuffer_ElectronCandidate_MC_Py = new Float_t[kMaxTracks];
  fBuffer_ElectronCandidate_MC_Pz = new Float_t[kMaxTracks];
  fBuffer_ElectronCandidate_MC_PDG = new Int_t[kMaxTracks];
  fBuffer_ElectronCandidate_MC_Mother_E = new Float_t[kMaxTracks];
  fBuffer_ElectronCandidate_MC_Mother_Px = new Float_t[kMaxTracks];
  fBuffer_ElectronCandidate_MC_Mother_Py = new Float_t[kMaxTracks];
  fBuffer_ElectronCandidate_MC_Mother_Pz = new Float_t[kMaxTracks];
  fBuffer_ElectronCandidate_MC_Mother_PDG = new Int_t[kMaxTracks];
  fBuffer_ElectronCandidate_MC_GrandMother_PDG = new Int_t[kMaxTracks];
  fBuffer_ElectronCandidate_MC_GrandMother_Pt = new Float_t[kMaxTracks];
  fBuffer_ElectronCandidate_MC_GGrandMother_PDG = new Int_t[kMaxTracks];
  fBuffer_ElectronCandidate_MC_GGrandMother_Pt = new Float_t[kMaxTracks];

  fBuffer_MuonCandidate_E = new Float_t[kMaxTracks];
  fBuffer_MuonCandidate_Px = new Float_t[kMaxTracks];
  fBuffer_MuonCandidate_Py = new Float_t[kMaxTracks];
  fBuffer_MuonCandidate_Pz = new Float_t[kMaxTracks];
  fBuffer_MuonCandidate_RAbsEnd = new Float_t[kMaxTracks];
  fBuffer_MuonCandidate_Chi2NDF = new Float_t[kMaxTracks];
  fBuffer_MuonCandidate_DCA = new Float_t[kMaxTracks];
  fBuffer_MuonCandidate_Charge = new Int_t[kMaxTracks];
  fBuffer_MuonCandidate_MatchTrigger = new Int_t[kMaxTracks];
  fBuffer_MuonCandidate_MC_E = new Float_t[kMaxTracks];
  fBuffer_MuonCandidate_MC_Px = new Float_t[kMaxTracks];
  fBuffer_MuonCandidate_MC_Py = new Float_t[kMaxTracks];
  fBuffer_MuonCandidate_MC_Pz = new Float_t[kMaxTracks];
  fBuffer_MuonCandidate_MC_PDG = new Int_t[kMaxTracks];
  fBuffer_MuonCandidate_MC_Mother_E = new Float_t[kMaxTracks];
  fBuffer_MuonCandidate_MC_Mother_Px = new Float_t[kMaxTracks];
  fBuffer_MuonCandidate_MC_Mother_Py = new Float_t[kMaxTracks];
  fBuffer_MuonCandidate_MC_Mother_Pz = new Float_t[kMaxTracks];
  fBuffer_MuonCandidate_MC_Mother_PDG = new Int_t[kMaxTracks];
  fBuffer_MuonCandidate_MC_GrandMother_PDG = new Int_t[kMaxTracks];
  fBuffer_MuonCandidate_MC_GrandMother_Pt = new Float_t[kMaxTracks];
  fBuffer_MuonCandidate_MC_GGrandMother_PDG = new Int_t[kMaxTracks];
  fBuffer_MuonCandidate_MC_GGrandMother_Pt = new Float_t[kMaxTracks];

  fBuffer_ClusterCandidate_E = new Float_t[kMaxTracks];
  fBuffer_ClusterCandidate_Px = new Float_t[kMaxTracks];
  fBuffer_ClusterCandidate_Py = new Float_t[kMaxTracks];
  fBuffer_ClusterCandidate_Pz = new Float_t[kMaxTracks];
  fBuffer_ClusterCandidate_Eta = new Float_t[kMaxTracks];
  fBuffer_ClusterCandidate_Phi = new Float_t[kMaxTracks];
  fBuffer_ClusterCandidate_M02 = new Float_t[kMaxTracks];
  fBuffer_ClusterCandidate_MC_E = new Float_t[kMaxTracks];
  fBuffer_ClusterCandidate_MC_Px = new Float_t[kMaxTracks];
  fBuffer_ClusterCandidate_MC_Py = new Float_t[kMaxTracks];
  fBuffer_ClusterCandidate_MC_Pz = new Float_t[kMaxTracks];
  fBuffer_ClusterCandidate_MC_PDG = new Int_t[kMaxTracks];
  fBuffer_ClusterCandidate_MC_Mother_E = new Float_t[kMaxTracks];
  fBuffer_ClusterCandidate_MC_Mother_Px = new Float_t[kMaxTracks];
  fBuffer_ClusterCandidate_MC_Mother_Py = new Float_t[kMaxTracks];
  fBuffer_ClusterCandidate_MC_Mother_Pz = new Float_t[kMaxTracks];
  fBuffer_ClusterCandidate_MC_Mother_PDG = new Int_t[kMaxTracks];
  fBuffer_ClusterCandidate_MC_GrandMother_PDG = new Int_t[kMaxTracks];

  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
  DefineOutput(2, TTree::Class());
}

//________________________________________________________________________
AliAnalysisTaskCorrelationTree::~AliAnalysisTaskCorrelationTree()
{
  // default deconstructor
}
//________________________________________________________________________
void AliAnalysisTaskCorrelationTree::UserCreateOutputObjects()
{
  // Create User Output Objects

  if (fOutputList != NULL)
  {
    delete fOutputList;
    fOutputList = NULL;
  }
  if (fOutputList == NULL)
  {
    fOutputList = new TList();
    fOutputList->SetOwner(true);
    fHistoNEvents = new TH1F("NEvents", "NEvents", 14, -0.5, 13.5);
    fHistoNEvents->GetXaxis()->SetBinLabel(1, "Accepted");
    fHistoNEvents->GetXaxis()->SetBinLabel(2, "Centrality");
    fHistoNEvents->GetXaxis()->SetBinLabel(3, "Miss. MC or inc. ev.");
    fHistoNEvents->GetXaxis()->SetBinLabel(4, "Trigger");
    fHistoNEvents->GetXaxis()->SetBinLabel(5, "Vertex Z");
    fHistoNEvents->GetXaxis()->SetBinLabel(6, "Cont. Vertex");
    fHistoNEvents->GetXaxis()->SetBinLabel(7, "Pile-Up");
    fHistoNEvents->GetXaxis()->SetBinLabel(8, "no SDD");
    fHistoNEvents->GetXaxis()->SetBinLabel(9, "no V0AND");
    fHistoNEvents->GetXaxis()->SetBinLabel(10, "EMCAL/TPC problems");
    fHistoNEvents->GetXaxis()->SetBinLabel(11, "rejectedForJetJetMC");
    fHistoNEvents->GetXaxis()->SetBinLabel(12, "SPD hits vs tracklet");
    fHistoNEvents->GetXaxis()->SetBinLabel(13, "Out-of-Bunch pileup Past-Future");
    fHistoNEvents->GetXaxis()->SetBinLabel(14, "Pileup V0M-TPCout Tracks");
    fOutputList->Add(fHistoNEvents);
  }
  fAnalysisTree = new TTree(Form("CorrelationTree_%s_%s", (fEventCuts->GetCutNumber()).Data(), (fConversionCuts->GetCutNumber()).Data()), Form("CorrelationTree_%s_%s", (fEventCuts->GetCutNumber()).Data(), (fConversionCuts->GetCutNumber()).Data()));

  fAnalysisTree->Branch("NVertexContributors", &fBuffer_NContributors, "NVertexContributors/I"); // max 200 for now
  fAnalysisTree->Branch("RunNumber", &fBuffer_RunNumber, "RunNumber/I");                         // max 200 for now
  fAnalysisTree->Branch("Vertex_Z", &fBuffer_VertexZ, "Vertex_Z/F");                             // max 200 for now
  fAnalysisTree->Branch("NEventTriggers", &fBuffer_NEventTriggers, "NEventTriggers/I");          // max 200 for now
  fAnalysisTree->Branch("EventTrigger", fBuffer_EventTrigger, "EventTrigger[NEventTriggers]/I"); // max 200 for now

  fAnalysisTree->Branch("NElectronCandidates", &fBuffer_NElectronCandidates, "NElectronCandidates/I"); // max 200 for now
  fAnalysisTree->Branch("ElectronCandidate_E", fBuffer_ElectronCandidate_E, "ElectronCandidate_E[NElectronCandidates]/F");
  fAnalysisTree->Branch("ElectronCandidate_Px", fBuffer_ElectronCandidate_Px, "ElectronCandidate_Px[NElectronCandidates]/F");
  fAnalysisTree->Branch("ElectronCandidate_Py", fBuffer_ElectronCandidate_Py, "ElectronCandidate_Py[NElectronCandidates]/F");
  fAnalysisTree->Branch("ElectronCandidate_Pz", fBuffer_ElectronCandidate_Pz, "ElectronCandidate_Pz[NElectronCandidates]/F");
  fAnalysisTree->Branch("ElectronCandidate_PropEta", fBuffer_ElectronCandidate_PropEta, "ElectronCandidate_PropEta[NElectronCandidates]/F");
  fAnalysisTree->Branch("ElectronCandidate_PropPhi", fBuffer_ElectronCandidate_PropPhi, "ElectronCandidate_PropPhi[NElectronCandidates]/F");
  fAnalysisTree->Branch("ElectronCandidate_Charge", fBuffer_ElectronCandidate_Charge, "ElectronCandidate_Charge[NElectronCandidates]/F");
  fAnalysisTree->Branch("ElectronCandidate_NSigmaElecTPC", fBuffer_ElectronCandidate_NSigmaElecTPC, "ElectronCandidate_NSigmaElecTPC[NElectronCandidates]/F");
  fAnalysisTree->Branch("ElectronCandidate_NSigmaElecTOF", fBuffer_ElectronCandidate_NSigmaElecTOF, "ElectronCandidate_NSigmaElecTOF[NElectronCandidates]/F");
  fAnalysisTree->Branch("ElectronCandidate_NSigmaPionTPC", fBuffer_ElectronCandidate_NSigmaPionTPC, "ElectronCandidate_NSigmaPionTPC[NElectronCandidates]/F");
  fAnalysisTree->Branch("ElectronCandidate_NSigmaPionTOF", fBuffer_ElectronCandidate_NSigmaPionTOF, "ElectronCandidate_NSigmaPionTOF[NElectronCandidates]/F");
  fAnalysisTree->Branch("ElectronCandidate_NSigmaKaonTPC", fBuffer_ElectronCandidate_NSigmaKaonTPC, "ElectronCandidate_NSigmaKaonTPC[NElectronCandidates]/F");
  fAnalysisTree->Branch("ElectronCandidate_NSigmaProtonTPC", fBuffer_ElectronCandidate_NSigmaProtonTPC, "ElectronCandidate_NSigmaProtonTPC[NElectronCandidates]/F");

  if (fIsMC)
  {
    fAnalysisTree->Branch("ElectronCandidate_MC_E", fBuffer_ElectronCandidate_MC_E, "ElectronCandidate_MC_E[NElectronCandidates]/F");
    fAnalysisTree->Branch("ElectronCandidate_MC_Px", fBuffer_ElectronCandidate_MC_Px, "ElectronCandidate_MC_Px[NElectronCandidates]/F");
    fAnalysisTree->Branch("ElectronCandidate_MC_Py", fBuffer_ElectronCandidate_MC_Py, "ElectronCandidate_MC_Py[NElectronCandidates]/F");
    fAnalysisTree->Branch("ElectronCandidate_MC_Pz", fBuffer_ElectronCandidate_MC_Pz, "ElectronCandidate_MC_Pz[NElectronCandidates]/F");
    fAnalysisTree->Branch("ElectronCandidate_MC_PDG", fBuffer_ElectronCandidate_MC_PDG, "ElectronCandidate_MC_PDG[NElectronCandidates]/I");
    fAnalysisTree->Branch("ElectronCandidate_MC_Mother_E", fBuffer_ElectronCandidate_MC_Mother_E, "ElectronCandidate_MC_Mother_E[NElectronCandidates]/F");
    fAnalysisTree->Branch("ElectronCandidate_MC_Mother_Px", fBuffer_ElectronCandidate_MC_Mother_Px, "ElectronCandidate_MC_Mother_Px[NElectronCandidates]/F");
    fAnalysisTree->Branch("ElectronCandidate_MC_Mother_Py", fBuffer_ElectronCandidate_MC_Mother_Py, "ElectronCandidate_MC_Mother_Py[NElectronCandidates]/F");
    fAnalysisTree->Branch("ElectronCandidate_MC_Mother_Pz", fBuffer_ElectronCandidate_MC_Mother_Pz, "ElectronCandidate_MC_Mother_Pz[NElectronCandidates]/F");
    fAnalysisTree->Branch("ElectronCandidate_MC_Mother_PDG", fBuffer_ElectronCandidate_MC_Mother_PDG, "ElectronCandidate_MC_Mother_PDG[NElectronCandidates]/I");
    fAnalysisTree->Branch("ElectronCandidate_MC_GrandMother_PDG", fBuffer_ElectronCandidate_MC_GrandMother_PDG, "ElectronCandidate_MC_GrandMother_PDG[NElectronCandidates]/I");
    fAnalysisTree->Branch("ElectronCandidate_MC_GrandMother_Pt", fBuffer_ElectronCandidate_MC_GrandMother_Pt, "ElectronCandidate_MC_GrandMother_Pt[NElectronCandidates]/F");
    fAnalysisTree->Branch("ElectronCandidate_MC_GGrandMother_PDG", fBuffer_ElectronCandidate_MC_GGrandMother_PDG, "ElectronCandidate_MC_GGrandMother_PDG[NElectronCandidates]/I");
    fAnalysisTree->Branch("ElectronCandidate_MC_GGrandMother_Pt", fBuffer_ElectronCandidate_MC_GGrandMother_Pt, "ElectronCandidate_MC_GGrandMother_Pt[NElectronCandidates]/F");
  }
  fAnalysisTree->Branch("NMuonCandidates", &fBuffer_NMuonCandidates, "NMuonCandidates/I"); // max 200 for now
  fAnalysisTree->Branch("MuonCandidate_E", fBuffer_MuonCandidate_E, "MuonCandidate_E[NMuonCandidates]/F");
  fAnalysisTree->Branch("MuonCandidate_Px", fBuffer_MuonCandidate_Px, "MuonCandidate_Px[NMuonCandidates]/F");
  fAnalysisTree->Branch("MuonCandidate_Py", fBuffer_MuonCandidate_Py, "MuonCandidate_Py[NMuonCandidates]/F");
  fAnalysisTree->Branch("MuonCandidate_Pz", fBuffer_MuonCandidate_Pz, "MuonCandidate_Pz[NMuonCandidates]/F");
  fAnalysisTree->Branch("MuonCandidate_RAbsEnd", fBuffer_MuonCandidate_RAbsEnd, "MuonCandidate_RAbsEnd[NMuonCandidates]/F");
  fAnalysisTree->Branch("MuonCandidate_Chi2NDF", fBuffer_MuonCandidate_Chi2NDF, "MuonCandidate_Chi2NDF[NMuonCandidates]/F");
  fAnalysisTree->Branch("MuonCandidate_DCA", fBuffer_MuonCandidate_DCA, "MuonCandidate_DCA[NMuonCandidates]/F");
  fAnalysisTree->Branch("MuonCandidate_Charge", fBuffer_MuonCandidate_Charge, "MuonCandidate_Charge[NMuonCandidates]/I");
  fAnalysisTree->Branch("MuonCandidate_MatchTrigger", fBuffer_MuonCandidate_MatchTrigger, "MuonCandidate_MatchTrigger[NMuonCandidates]/I");

  if (fIsMC)
  {
    fAnalysisTree->Branch("MuonCandidate_MC_E", fBuffer_MuonCandidate_MC_E, "MuonCandidate_MC_E[NMuonCandidates]/F");
    fAnalysisTree->Branch("MuonCandidate_MC_Px", fBuffer_MuonCandidate_MC_Px, "MuonCandidate_MC_Px[NMuonCandidates]/F");
    fAnalysisTree->Branch("MuonCandidate_MC_Py", fBuffer_MuonCandidate_MC_Py, "MuonCandidate_MC_Py[NMuonCandidates]/F");
    fAnalysisTree->Branch("MuonCandidate_MC_Pz", fBuffer_MuonCandidate_MC_Pz, "MuonCandidate_MC_Pz[NMuonCandidates]/F");
    fAnalysisTree->Branch("MuonCandidate_MC_PDG", fBuffer_MuonCandidate_MC_PDG, "MuonCandidate_MC_PDG[NMuonCandidates]/I");
    fAnalysisTree->Branch("MuonCandidate_MC_Mother_E", fBuffer_MuonCandidate_MC_Mother_E, "MuonCandidate_MC_Mother_E[NMuonCandidates]/F");
    fAnalysisTree->Branch("MuonCandidate_MC_Mother_Px", fBuffer_MuonCandidate_MC_Mother_Px, "MuonCandidate_MC_Mother_Px[NMuonCandidates]/F");
    fAnalysisTree->Branch("MuonCandidate_MC_Mother_Py", fBuffer_MuonCandidate_MC_Mother_Py, "MuonCandidate_MC_Mother_Py[NMuonCandidates]/F");
    fAnalysisTree->Branch("MuonCandidate_MC_Mother_Pz", fBuffer_MuonCandidate_MC_Mother_Pz, "MuonCandidate_MC_Mother_Pz[NMuonCandidates]/F");
    fAnalysisTree->Branch("MuonCandidate_MC_Mother_PDG", fBuffer_MuonCandidate_MC_Mother_PDG, "MuonCandidate_MC_Mother_PDG[NMuonCandidates]/I");
    fAnalysisTree->Branch("MuonCandidate_MC_GrandMother_PDG", fBuffer_MuonCandidate_MC_GrandMother_PDG, "MuonCandidate_MC_GrandMother_PDG[NMuonCandidates]/I");
    fAnalysisTree->Branch("MuonCandidate_MC_GrandMother_Pt", fBuffer_MuonCandidate_MC_GrandMother_Pt, "MuonCandidate_MC_GrandMother_Pt[NMuonCandidates]/F");
    fAnalysisTree->Branch("MuonCandidate_MC_GGrandMother_PDG", fBuffer_MuonCandidate_MC_GGrandMother_PDG, "MuonCandidate_MC_GGrandMother_PDG[NMuonCandidates]/I");
    fAnalysisTree->Branch("MuonCandidate_MC_GGrandMother_Pt", fBuffer_MuonCandidate_MC_GGrandMother_Pt, "MuonCandidate_MC_GGrandMother_Pt[NMuonCandidates]/F");
  }
  fAnalysisTree->Branch("NClusterCandidates", &fBuffer_NClusterCandidates, "NClusterCandidates/I"); // max 200 for now
  fAnalysisTree->Branch("ClusterCandidate_E", fBuffer_ClusterCandidate_E, "ClusterCandidate_E[NClusterCandidates]/F");
  fAnalysisTree->Branch("ClusterCandidate_Px", fBuffer_ClusterCandidate_Px, "ClusterCandidate_Px[NClusterCandidates]/F");
  fAnalysisTree->Branch("ClusterCandidate_Py", fBuffer_ClusterCandidate_Py, "ClusterCandidate_Py[NClusterCandidates]/F");
  fAnalysisTree->Branch("ClusterCandidate_Pz", fBuffer_ClusterCandidate_Pz, "ClusterCandidate_Pz[NClusterCandidates]/F");
  fAnalysisTree->Branch("ClusterCandidate_Eta", fBuffer_ClusterCandidate_Eta, "ClusterCandidate_Eta[NClusterCandidates]/F");
  fAnalysisTree->Branch("ClusterCandidate_Phi", fBuffer_ClusterCandidate_Phi, "ClusterCandidate_Phi[NClusterCandidates]/F");
  fAnalysisTree->Branch("ClusterCandidate_M02", fBuffer_ClusterCandidate_M02, "ClusterCandidate_M02[NClusterCandidates]/F");

  if (fIsMC)
  {
    fAnalysisTree->Branch("ClusterCandidate_MC_E", fBuffer_ClusterCandidate_MC_E, "ClusterCandidate_MC_E[NClusterCandidates]/F");
    fAnalysisTree->Branch("ClusterCandidate_MC_Px", fBuffer_ClusterCandidate_MC_Px, "ClusterCandidate_MC_Px[NClusterCandidates]/F");
    fAnalysisTree->Branch("ClusterCandidate_MC_Py", fBuffer_ClusterCandidate_MC_Py, "ClusterCandidate_MC_Py[NClusterCandidates]/F");
    fAnalysisTree->Branch("ClusterCandidate_MC_Pz", fBuffer_ClusterCandidate_MC_Pz, "ClusterCandidate_MC_Pz[NClusterCandidates]/F");
    fAnalysisTree->Branch("ClusterCandidate_MC_PDG", fBuffer_ClusterCandidate_MC_PDG, "ClusterCandidate_MC_PDG[NClusterCandidates]/I");
    fAnalysisTree->Branch("ClusterCandidate_MC_Mother_E", fBuffer_ClusterCandidate_MC_Mother_E, "ClusterCandidate_MC_Mother_E[NClusterCandidates]/F");
    fAnalysisTree->Branch("ClusterCandidate_MC_Mother_Px", fBuffer_ClusterCandidate_MC_Mother_Px, "ClusterCandidate_MC_Mother_Px[NClusterCandidates]/F");
    fAnalysisTree->Branch("ClusterCandidate_MC_Mother_Py", fBuffer_ClusterCandidate_MC_Mother_Py, "ClusterCandidate_MC_Mother_Py[NClusterCandidates]/F");
    fAnalysisTree->Branch("ClusterCandidate_MC_Mother_Pz", fBuffer_ClusterCandidate_MC_Mother_Pz, "ClusterCandidate_MC_Mother_Pz[NClusterCandidates]/F");
    fAnalysisTree->Branch("ClusterCandidate_MC_Mother_PDG", fBuffer_ClusterCandidate_MC_Mother_PDG, "ClusterCandidate_MC_Mother_PDG[NClusterCandidates]/I");
    fAnalysisTree->Branch("ClusterCandidate_MC_GrandMother_PDG", fBuffer_ClusterCandidate_MC_GrandMother_PDG, "ClusterCandidate_MC_GrandMother_PDG[NClusterCandidates]/I");
  }
  fV0Reader = (AliV0ReaderV1 *)AliAnalysisManager::GetAnalysisManager()->GetTask(fV0ReaderName.Data());

  if (fV0Reader && fV0Reader->GetProduceV0FindingEfficiency())
    if (fV0Reader->GetV0FindingEfficiencyHistograms())
      fOutputList->Add(fV0Reader->GetV0FindingEfficiencyHistograms());

  PostData(1, fOutputList);
  OpenFile(2);
  PostData(2, fAnalysisTree);
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskCorrelationTree::Notify()
{
  if (fEventCuts->GetPeriodEnum() == AliConvEventCuts::kNoPeriod && ((AliConvEventCuts *)fV0Reader->GetEventCuts())->GetPeriodEnum() != AliConvEventCuts::kNoPeriod)
  {
    fEventCuts->SetPeriodEnumExplicit(((AliConvEventCuts *)fV0Reader->GetEventCuts())->GetPeriodEnum());
  }
  else if (fEventCuts->GetPeriodEnum() == AliConvEventCuts::kNoPeriod)
  {
    fEventCuts->SetPeriodEnum(fV0Reader->GetPeriodName());
  }

  if (!fEventCuts->GetDoEtaShift())
    return true; // No Eta Shift requested, continue

  if (fEventCuts->GetEtaShift() == 0.0)
  { // Eta Shift requested but not set, get shift automatically
    fEventCuts->GetCorrectEtaShiftFromPeriod();
    fEventCuts->DoEtaShift(false); // Eta Shift Set, make sure that it is called only once
    return true;
  }
  else
  {
    printf(" Gamma Conversion QA Task %s :: Eta Shift Manually Set to %f \n\n",
           (fEventCuts->GetCutNumber()).Data(), fEventCuts->GetEtaShift());
    fEventCuts->DoEtaShift(false); // Eta Shift Set, make sure that it is called only once
  }

  return true;
}
//________________________________________________________________________
void AliAnalysisTaskCorrelationTree::UserExec(Option_t *)
{
  ResetBuffer();
  fInputEvent = InputEvent();

  Double_t fMaxVertexZ = 10.0;

  Int_t eventQuality = ((AliConvEventCuts *)fV0Reader->GetEventCuts())->GetEventQuality();
  if (fInputEvent->IsIncompleteDAQ() == true)
    eventQuality = 2; // incomplete event


  if (eventQuality == 2)
  {
    // cout << "Event not accepted with event quality " << eventQuality << endl;
    return;
  }

  if (fIsMC)
  {
    fMCEvent = MCEvent();
    if (!fMCEvent)
    {
      AliError("ERROR: Could not retrieve MC event");
      return;
    }
    // if(!fAODMCTrackArray) fAODMCTrackArray = dynamic_cast<TClonesArray*>(fMCEvent->FindListObject(AliAODMCParticle::StdBranchName()));
    //   if (fAODMCTrackArray == NULL){
    //     AliError("No MC particle list available in AOD");
    //     return;
    //   }
  }

  Int_t eventNotAccepted =
      fEventCuts->IsEventAcceptedByCut(fV0Reader->GetEventCuts(), fInputEvent, fMCEvent, fIsHeavyIon, false);
  
  if(!eventNotAccepted && eventQuality!=3) fHistoNEvents->Fill(eventQuality);

  // if (eventNotAccepted)
  // {
  //   cout << "Event not accepted" << endl;
  //   return; // Check Centrality, PileUp, SDD and V0AND --> Not Accepted => eventQuality = 1
  // }
  // else
  // {
  //   cout << "Event accepted" << endl;
  // }
  // if (fInputEvent->GetPrimaryVertex()->GetNContributors() < 1)
  // {
  //   // cout << "vertex not reconstructed from tracks/tracklets" << endl;
  //   // vtx position is the mean vertex in x,y and z=0
  //   // return;
  // }
  fBuffer_NContributors = fInputEvent->GetPrimaryVertex()->GetNContributors();
  fBuffer_RunNumber = fInputEvent->GetRunNumber();
  if (TMath::Abs(fInputEvent->GetPrimaryVertex()->GetZ()) > fMaxVertexZ)
    return;
  fBuffer_VertexZ = fInputEvent->GetPrimaryVertex()->GetZ();
  fBuffer_NEventTriggers = 0;
  TString firedTrigClass = fInputEvent->GetFiredTriggerClasses();
  const int nTriggers = 13;
  TString triggerClasses[nTriggers] = {"CEMC7MUL-B-NOPF-ALLNOTRD", "CEMC7MSL-B-NOPF-ALLNOTRD", "CEMC7MSH-B-NOPF-ALLNOTRD", "CDMC7MUL-B-NOPF-ALLNOTRD", "CDMC7MSL-B-NOPF-ALLNOTRD", "CDMC7MSH-B-NOPF-ALLNOTRD", "CMUL7-B-NOPF-ALLNOTRD", "CMSL7-B-NOPF-ALLNOTRD", "CMSH7-B-NOPF-ALLNOTRD", "CEMC7-B-NOPF-ALLNOTRD", "CDMC7-B-NOPF-ALLNOTRD", "CPHI7MSL-B-NOPF-ALLNOTRD", "CPHI7-B-NOPF-ALLNOTRD"};
  if (fIsMC)
  {
    fBuffer_EventTrigger[fBuffer_NEventTriggers] = 0;
    fBuffer_NEventTriggers++;
  }
  else
  {
    for (Int_t i = 0; i < nTriggers; i++)
    {
      if (firedTrigClass.Contains(triggerClasses[i]))
      {
        fBuffer_EventTrigger[fBuffer_NEventTriggers] = i;
        fBuffer_NEventTriggers++;
        // cout << "\n\tTrigger fired: " << firedTrigClass << endl;
      }
    }
  }
  if (fBuffer_NEventTriggers == 0)
  {
    // cout << "No trigger fired" << endl;
    return;
  }

  fPidResponse = fInputHandler->GetPIDResponse();
  if (!fPidResponse)
  {
    AliError("No PID Response");
    return;
  }

  // reset tree buffer for each event
  ProcessClusters();
  ProcessElectrons();
  ProcessMuons();

  if (fBuffer_NMuonCandidates == 0)
    return;
  if (fBuffer_NElectronCandidates == 0 && fBuffer_NClusterCandidates == 0)
    return;

  // fill tree for each accepted event
  if (fAnalysisTree)
  {
    fAnalysisTree->Fill();
  }

  PostData(1, fOutputList);
}

//_________________________________________________
void AliAnalysisTaskCorrelationTree::ProcessElectrons()
{
  AliVEvent *event = (AliVEvent *)InputEvent();
  fBuffer_NElectronCandidates = 0;
  Int_t ntracks = event->GetNumberOfTracks();

  double absEtaCut = 0.9;
  // double fTPCnCrossedRows = -0.8;

  double fPtCutMainEle = 0.4;
  // cout << "ntracks: " << ntracks << endl;
  for (Int_t j = 0; j < ntracks; j++)
  {
    // AliVTrack *track = dynamic_cast<AliVTrack*>(Vtrack);
    AliAODTrack *track = dynamic_cast<AliAODTrack *>(event->GetTrack(j));
    if (!track)
      AliFatal("Not a standard AOD");

    double eta = track->Eta();
    if (abs(eta) > absEtaCut)
      continue;

    if (!track->TestFilterBit(16))
      continue; // standard loose DCA cuts

    if (track->GetTPCNcls() < 70)
      continue;
    if (track->GetTPCNCrossedRows() < 100)
      continue;
    Double_t nclsFindable = (Double_t)track->GetNcls(1) / (Double_t)track->GetTPCNclsF(); // Ncluster/Nfindablecluster
    if (nclsFindable < 0.8)
      continue;

    Int_t nClusterITS = track->GetITSNcls();
    if (nClusterITS < 2)
      continue;

    if (!track->IsHybridGlobalConstrainedGlobal())
      continue;

    // if (!((track->GetStatus() & AliVTrack::kITSrefit)) || (!(track->GetStatus()&AliESDtrack::kTPCrefit)))
    //   continue;

    // if (track->GetTPCCrossedRows() < fTPCnCrossedRows)
    //   continue;
    // if(fTPCandITSrefit){
    //     if((!(atrack->GetStatus()&AliESDtrack::kITSrefit))|| (!(atrack->GetStatus()&AliESDtrack::kTPCrefit))) continue;
    // }

    // kAny
    // if (fITSpixel == 1)
    // {
    // if (!(track->HasPointOnITSLayer(0) || track->HasPointOnITSLayer(1)))
    //   continue;
    // }
    // // kBoth
    // if (fITSpixel == 2)
    // {
    //   if (!(atrack->HasPointOnITSLayer(0) && atrack->HasPointOnITSLayer(1)))
    //     continue;
    // }
    // // kFirst
    // if (fITSpixel == 3)
    // {
    //   if (!(atrack->HasPointOnITSLayer(0)))
    //     continue;
    // }

    // // ITS Ncls
    // if ((atrack->GetITSNcls()) < fITSncls)
    //   continue;

    // // DCA cut
    // Double_t d0z0[2], cov[3];
    // AliAODVertex *pVtx = fAOD->GetPrimaryVertex();
    // if (atrack->PropagateToDCA(pVtx, fAOD->GetMagneticField(), 20., d0z0, cov))
    // {
    //   if (TMath::Abs(d0z0[0]) > fDCAxyCut || TMath::Abs(d0z0[1]) > fDCAzCut)
    //     continue;
    // }
    // // Pt cut
    Double_t fPt = track->Pt();
    if (fPt < fPtCutMainEle)
      continue;

    // // // reject kink mother
    // if (fRejectKinkMother)
    // {
    //   Bool_t kinkmotherpass = true;
    //   for (Int_t kinkmother = 0; kinkmother < fNumberOfMotherkink; kinkmother++)
    //   {
    //     if (track->GetID() == fListOfmotherkink[kinkmother])
    //     {
    //       kinkmotherpass = false;
    //       continue;
    //     }
    //   }
    //   if (!kinkmotherpass)
    //     continue;
    // }

    // // chi2 per cluster
    // if (((track->GetTPCchi2()) / (atrack->GetTPCNcls())) > fTPCchi2)
    // {
    //   continue;
    // }

    // // ITS Chi2
    // if (((atrack->GetITSchi2()) / (atrack->GetITSNcls())) > fITSchi2)
    // {
    //   continue;
    // }

    // if (fAODGlobalTracks)
    // {
    //   if (!atrack->TestFilterMask(AliAODTrack::kTrkGlobalNoDCA))
    //     continue; // mimimum cuts
    // }

    fBuffer_ElectronCandidate_E[fBuffer_NElectronCandidates] = track->E();
    fBuffer_ElectronCandidate_Px[fBuffer_NElectronCandidates] = track->Px();
    fBuffer_ElectronCandidate_Py[fBuffer_NElectronCandidates] = track->Py();
    fBuffer_ElectronCandidate_Pz[fBuffer_NElectronCandidates] = track->Pz();
    fBuffer_ElectronCandidate_Charge[fBuffer_NElectronCandidates] = track->Charge();
    fBuffer_ElectronCandidate_NSigmaElecTPC[fBuffer_NElectronCandidates] = fPidResponse->NumberOfSigmasTPC(track, AliPID::kElectron);
    fBuffer_ElectronCandidate_NSigmaElecTOF[fBuffer_NElectronCandidates] = fPidResponse->NumberOfSigmasTOF(track, AliPID::kElectron);
    fBuffer_ElectronCandidate_NSigmaKaonTPC[fBuffer_NElectronCandidates] = fPidResponse->NumberOfSigmasTPC(track, AliPID::kKaon);
    fBuffer_ElectronCandidate_NSigmaPionTPC[fBuffer_NElectronCandidates] = fPidResponse->NumberOfSigmasTPC(track, AliPID::kPion);
    fBuffer_ElectronCandidate_NSigmaPionTOF[fBuffer_NElectronCandidates] = fPidResponse->NumberOfSigmasTOF(track, AliPID::kPion);
    fBuffer_ElectronCandidate_NSigmaProtonTPC[fBuffer_NElectronCandidates] = fPidResponse->NumberOfSigmasTPC(track, AliPID::kProton);

    bool doPropagation = true;
    if (doPropagation)
    {
      AliExternalTrackParam *trackParam = 0;
      if (TMath::Abs(track->GetTrackEtaOnEMCal()) < 0.75)
      {
        if (!(track->GetTrackPhiOnEMCal() < 70 * TMath::DegToRad() || track->GetTrackPhiOnEMCal() > 190 * TMath::DegToRad()) && (track->GetTrackPhiOnEMCal() < 250 * TMath::DegToRad() || track->GetTrackPhiOnEMCal() > 340 * TMath::DegToRad()))
        {
          Double_t xyz[3] = {0}, pxpypz[3] = {0}, cv[21] = {0};
          track->GetPxPyPz(pxpypz);
          track->GetXYZ(xyz);
          track->GetCovarianceXYZPxPyPz(cv);
          trackParam = new AliExternalTrackParam(xyz, pxpypz, cv, track->Charge());
          Float_t eta, phi, pt;
          AliExternalTrackParam emcParam(*trackParam);

          if (!AliEMCALRecoUtils::ExtrapolateTrackToEMCalSurface(&emcParam, 440., 0.139, 20., eta, phi, pt))
          {
            delete trackParam;
          }
          else
          {
            if (TMath::Abs(eta) < 0.75 && ((phi > 70 * TMath::DegToRad() && phi < 190 * TMath::DegToRad()) || (phi > 250 * TMath::DegToRad() && phi < 340 * TMath::DegToRad())))
            {
              fBuffer_ElectronCandidate_PropEta[fBuffer_NElectronCandidates] = eta;
              fBuffer_ElectronCandidate_PropPhi[fBuffer_NElectronCandidates] = phi;
            }
            delete trackParam;
          }
        }
      }
    }
    // cout << "\ttrack " << j << "\tpx " << track->Px() << "\tpy " << track->Py() << "\tpz " << track->Pz() << "\tE " << track->E() << endl;

    if (fIsMC)
    {
      // cout << "\ttrack " << j << "\tpx " << track->Px() << "\tpy " << track->Py() << "\tpz " << track->Pz() << "\tE " << track->E() << endl;
      // AliAODMCParticle *mcTrack = (AliAODMCParticle *)fAODMCTrackArray->At(TMath::Abs(track->GetLabel()));
      AliAODMCParticle *mcTrack = dynamic_cast<AliAODMCParticle *>(fMCEvent->GetTrack(track->GetLabel()));

      if (!mcTrack)
      {
        AliError("Could not receive MC particle");
        continue;
      }

      fBuffer_ElectronCandidate_MC_E[fBuffer_NElectronCandidates] = mcTrack->E();
      fBuffer_ElectronCandidate_MC_Px[fBuffer_NElectronCandidates] = mcTrack->Px();
      fBuffer_ElectronCandidate_MC_Py[fBuffer_NElectronCandidates] = mcTrack->Py();
      fBuffer_ElectronCandidate_MC_Pz[fBuffer_NElectronCandidates] = mcTrack->Pz();
      fBuffer_ElectronCandidate_MC_PDG[fBuffer_NElectronCandidates] = mcTrack->GetPdgCode();
      AliAODMCParticle *mcMother = dynamic_cast<AliAODMCParticle *>(fMCEvent->GetTrack(mcTrack->GetMother()));
      if (mcMother)
      {
        fBuffer_ElectronCandidate_MC_Mother_E[fBuffer_NElectronCandidates] = mcMother->E();
        fBuffer_ElectronCandidate_MC_Mother_Px[fBuffer_NElectronCandidates] = mcMother->Px();
        fBuffer_ElectronCandidate_MC_Mother_Py[fBuffer_NElectronCandidates] = mcMother->Py();
        fBuffer_ElectronCandidate_MC_Mother_Pz[fBuffer_NElectronCandidates] = mcMother->Pz();
        fBuffer_ElectronCandidate_MC_Mother_PDG[fBuffer_NElectronCandidates] = mcMother->GetPdgCode();

        AliAODMCParticle *mcGrandMother = dynamic_cast<AliAODMCParticle *>(fMCEvent->GetTrack(mcMother->GetMother()));
        if (mcGrandMother)
        {
          fBuffer_ElectronCandidate_MC_GrandMother_PDG[fBuffer_NElectronCandidates] = mcGrandMother->GetPdgCode();
          fBuffer_ElectronCandidate_MC_GrandMother_Pt[fBuffer_NElectronCandidates] = mcGrandMother->Pt();

          AliAODMCParticle *mcGGrandMother = dynamic_cast<AliAODMCParticle *>(fMCEvent->GetTrack(mcGrandMother->GetMother()));
          if (mcGGrandMother)
          {
            fBuffer_ElectronCandidate_MC_GGrandMother_PDG[fBuffer_NElectronCandidates] = mcGGrandMother->GetPdgCode();
            fBuffer_ElectronCandidate_MC_GGrandMother_PDG[fBuffer_NElectronCandidates] = mcGGrandMother->Pt();
          }
        }
      }
      // fBuffer_ElectronCandidate_MCStatusCode[fBuffer_NElectronCandidates] = mcTrack->GetStatusCode();
      // fBuffer_ElectronCandidate_MCIsPrimary[fBuffer_NElectronCandidates] = mcTrack->IsPrimary();
      // fBuffer_ElectronCandidate_MCIsSecondaryFromWeakDecay[fBuffer_NElectronCandidates] = mcTrack->IsSecondaryFromWeakDecay();
      // fBuffer_ElectronCandidate_MCIsSecondaryFromMaterial[fBuffer_NElectronCandidates] = mcTrack->IsSecondaryFromMaterial();
      // fBuffer_ElectronCandidate_MCIsPhysicalPrimary[fBuffer_NElectronCandidates] = mcTrack->IsPhysicalPrimary();
      // fBuffer_ElectronCandidate_MCIsSecondaryFromWeakDecayOrMaterial[fBuffer_NElectronCandidates] = mcTrack->IsSecondaryFromWeakDecay() || mcTrack->IsSecondaryFromMaterial();
    }

    fBuffer_NElectronCandidates++;
    if (kMaxTracks <= fBuffer_NElectronCandidates)
    {
      AliFatal("Too many electrons");
    }
    // delete mutrack;
  }
}

//_________________________________________________
void AliAnalysisTaskCorrelationTree::ProcessMuons()
{
  AliVEvent *event = (AliVEvent *)InputEvent();
  fBuffer_NMuonCandidates = 0;
  Double_t fEtaCutMax = -2.5;
  Double_t fEtaCutMin = -4.0;
  Double_t fPtCutMainMuon = 0.1;

  Int_t ntracks = event->GetNumberOfTracks();

  enum
  {
    kMuEta = BIT(0),
    kMuThetaAbs = BIT(1),
    kMuPdca = BIT(2),
    kMuMatchApt = BIT(3),
    kMuMatchLpt = BIT(4),
    kMuMatchHpt = BIT(5),
    kMuTrackChiSquare = BIT(6)
  };

  for (Int_t j = 0; j < ntracks; j++)
  {
    AliAODTrack *mutrack = dynamic_cast<AliAODTrack *>(event->GetTrack(j));
    if (!mutrack)
      AliFatal("Not a standard AOD");
    if (!mutrack->IsMuonTrack())
      continue;

    double eta = mutrack->Eta();
    if (eta > fEtaCutMax || eta < fEtaCutMin)
      continue;

    // Int_t matchTrig = AliAnalysisMuonUtility::GetMatchTrigger(track);
    Int_t matchTrig = static_cast<const AliAODTrack *>(mutrack)->GetMatchTrigger();
    // Int_t cutLevel[3] = {kMuMatchApt, kMuMatchLpt, kMuMatchHpt};

    Double_t pt = mutrack->Pt();
    // for ( Int_t ilevel=0; ilevel<3; ilevel++ ) {
    //   if ( matchTrig < ilevel+1 ) break;
    //   if ( fSharpPtCut && pt < fOADBParam.GetSharpPtCut(ilevel) ) break;
    //   selectionMask |= cutLevel[ilevel];
    // }

    if (pt < fPtCutMainMuon)
      continue;

    fBuffer_MuonCandidate_E[fBuffer_NMuonCandidates] = mutrack->E();
    fBuffer_MuonCandidate_Px[fBuffer_NMuonCandidates] = mutrack->Px();
    fBuffer_MuonCandidate_Py[fBuffer_NMuonCandidates] = mutrack->Py();
    fBuffer_MuonCandidate_Pz[fBuffer_NMuonCandidates] = mutrack->Pz();
    fBuffer_MuonCandidate_RAbsEnd[fBuffer_NMuonCandidates] = mutrack->GetRAtAbsorberEnd();
    fBuffer_MuonCandidate_Charge[fBuffer_NMuonCandidates] = mutrack->Charge();
    fBuffer_MuonCandidate_Chi2NDF[fBuffer_NMuonCandidates] = mutrack->Chi2perNDF();
    fBuffer_MuonCandidate_DCA[fBuffer_NMuonCandidates] = mutrack->DCA();
    fBuffer_MuonCandidate_MatchTrigger[fBuffer_NMuonCandidates] = matchTrig;

    if (fIsMC)
    {
      // cout << "\ttrack " << j << "\tpx " << track->Px() << "\tpy " << track->Py() << "\tpz " << track->Pz() << "\tE " << track->E() << endl;
      // AliAODMCParticle *mcTrack = (AliAODMCParticle *)fAODMCTrackArray->At(TMath::Abs(track->GetLabel()));
      AliAODMCParticle *mcTrack = dynamic_cast<AliAODMCParticle *>(fMCEvent->GetTrack(mutrack->GetLabel()));

      if (!mcTrack)
      {
        AliError("Could not receive MC particle");
        continue;
      }

      fBuffer_MuonCandidate_MC_E[fBuffer_NMuonCandidates] = mcTrack->E();
      fBuffer_MuonCandidate_MC_Px[fBuffer_NMuonCandidates] = mcTrack->Px();
      fBuffer_MuonCandidate_MC_Py[fBuffer_NMuonCandidates] = mcTrack->Py();
      fBuffer_MuonCandidate_MC_Pz[fBuffer_NMuonCandidates] = mcTrack->Pz();
      fBuffer_MuonCandidate_MC_PDG[fBuffer_NMuonCandidates] = mcTrack->GetPdgCode();
      AliAODMCParticle *mcMother = dynamic_cast<AliAODMCParticle *>(fMCEvent->GetTrack(mcTrack->GetMother()));
      if (mcMother)
      {
        fBuffer_MuonCandidate_MC_Mother_E[fBuffer_NMuonCandidates] = mcMother->E();
        fBuffer_MuonCandidate_MC_Mother_Px[fBuffer_NMuonCandidates] = mcMother->Px();
        fBuffer_MuonCandidate_MC_Mother_Py[fBuffer_NMuonCandidates] = mcMother->Py();
        fBuffer_MuonCandidate_MC_Mother_Pz[fBuffer_NMuonCandidates] = mcMother->Pz();
        fBuffer_MuonCandidate_MC_Mother_PDG[fBuffer_NMuonCandidates] = mcMother->GetPdgCode();

        AliAODMCParticle *mcGrandMother = dynamic_cast<AliAODMCParticle *>(fMCEvent->GetTrack(mcMother->GetMother()));
        if (mcGrandMother)
        {
          fBuffer_MuonCandidate_MC_GrandMother_PDG[fBuffer_NMuonCandidates] = mcGrandMother->GetPdgCode();
          fBuffer_MuonCandidate_MC_GrandMother_Pt[fBuffer_NMuonCandidates] = mcGrandMother->Pt();

          AliAODMCParticle *mcGGrandMother = dynamic_cast<AliAODMCParticle *>(fMCEvent->GetTrack(mcGrandMother->GetMother()));
          if (mcGGrandMother)
          {
            fBuffer_MuonCandidate_MC_GGrandMother_PDG[fBuffer_NMuonCandidates] = mcGGrandMother->GetPdgCode();
            fBuffer_MuonCandidate_MC_GGrandMother_Pt[fBuffer_NMuonCandidates] = mcGGrandMother->Pt();
          }
        }
      }

      // fBuffer_MuonCandidate_MCStatusCode[fBuffer_NMuonCandidates] = mcTrack->GetStatusCode();
      // fBuffer_MuonCandidate_MCIsPrimary[fBuffer_NMuonCandidates] = mcTrack->IsPrimary();
      // fBuffer_MuonCandidate_MCIsSecondaryFromWeakDecay[fBuffer_NMuonCandidates] = mcTrack->IsSecondaryFromWeakDecay();
      // fBuffer_MuonCandidate_MCIsSecondaryFromMaterial[fBuffer_NMuonCandidates] = mcTrack->IsSecondaryFromMaterial();
      // fBuffer_MuonCandidate_MCIsPhysicalPrimary[fBuffer_NMuonCandidates] = mcTrack->IsPhysicalPrimary();
      // fBuffer_MuonCandidate_MCIsSecondaryFromWeakDecayOrMaterial[fBuffer_NMuonCandidates] = mcTrack->IsSecondaryFromWeakDecay() || mcTrack->IsSecondaryFromMaterial();
    }

    fBuffer_NMuonCandidates++;
    if (kMaxTracks <= fBuffer_NMuonCandidates)
    {
      AliFatal("Too many muons");
    }
    // delete mutrack;
  }
}
//_________________________________________________
void AliAnalysisTaskCorrelationTree::ProcessClusters()
{
  Int_t nclus = 0;
  fBuffer_NClusterCandidates = 0;
  float WeightJetJetMC = 1;
  Double_t tempClusterWeight = WeightJetJetMC;
  if (!fCorrTaskSetting.CompareTo(""))
  {
    nclus = fInputEvent->GetNumberOfCaloClusters();
  }
  else
  {
    if (!farrClustersProcess)
      farrClustersProcess = dynamic_cast<TClonesArray *>(fInputEvent->FindListObject(Form("%sClustersBranch", fCorrTaskSetting.Data())));
    if (!farrClustersProcess)
      AliFatal(Form("%sClustersBranch was not found in AliAnalysisTaskCorrelationTree! Check the correction framework settings!", fCorrTaskSetting.Data()));
    nclus = farrClustersProcess->GetEntries();
  }
  if (nclus == 0)
    return;
  // cout << "nclus: " << nclus << endl;

  fClusterCuts->MatchTracksToClusters(fInputEvent, WeightJetJetMC, true, fMCEvent);

  // vertex
  Double_t vertex[3] = {0};
  InputEvent()->GetPrimaryVertex()->GetXYZ(vertex);

  for (Long_t i = 0; i < nclus; i++)
  {

    AliVCluster *clus = NULL;
    if (farrClustersProcess)
      clus = new AliAODCaloCluster(*(AliAODCaloCluster *)farrClustersProcess->At(i));
    else
      clus = new AliAODCaloCluster(*(AliAODCaloCluster *)fInputEvent->GetCaloCluster(i));

    if (!clus)
      continue;

    // if open cluster cuts are not fullfilled I can abort
    if (!fClusterCuts->ClusterIsSelected(clus, fInputEvent, fMCEvent, fIsMC, tempClusterWeight, i))
    {
      delete clus;
      continue;
    }

    // // TLorentzvector with cluster
    TLorentzVector clusterVector;
    clus->GetMomentum(clusterVector, vertex);

    // TLorentzVector* tmpvec = new TLorentzVector();
    // tmpvec->SetPxPyPzE(clusterVector.Px(),clusterVector.Py(),clusterVector.Pz(),clusterVector.E());

    fBuffer_ClusterCandidate_E[fBuffer_NClusterCandidates] = clus->E();
    fBuffer_ClusterCandidate_Px[fBuffer_NClusterCandidates] = clusterVector.Px();
    fBuffer_ClusterCandidate_Py[fBuffer_NClusterCandidates] = clusterVector.Py();
    fBuffer_ClusterCandidate_Pz[fBuffer_NClusterCandidates] = clusterVector.Pz();
    fBuffer_ClusterCandidate_Eta[fBuffer_NClusterCandidates] = clusterVector.Eta();
    fBuffer_ClusterCandidate_Phi[fBuffer_NClusterCandidates] = clusterVector.Phi();
    fBuffer_ClusterCandidate_M02[fBuffer_NClusterCandidates] = clus->GetM02();

    if (fIsMC)
    {
      Int_t mcLabel = clus->GetLabel();
      if (mcLabel >= 0)
      {
        AliAODMCParticle *mcTrack = dynamic_cast<AliAODMCParticle *>(fMCEvent->GetTrack(mcLabel));
        if (mcTrack)
        {
          fBuffer_ClusterCandidate_MC_E[fBuffer_NClusterCandidates] = mcTrack->E();
          fBuffer_ClusterCandidate_MC_Px[fBuffer_NClusterCandidates] = mcTrack->Px();
          fBuffer_ClusterCandidate_MC_Py[fBuffer_NClusterCandidates] = mcTrack->Py();
          fBuffer_ClusterCandidate_MC_Pz[fBuffer_NClusterCandidates] = mcTrack->Pz();
          fBuffer_ClusterCandidate_MC_PDG[fBuffer_NClusterCandidates] = mcTrack->GetPdgCode();
          AliAODMCParticle *mcMother = dynamic_cast<AliAODMCParticle *>(fMCEvent->GetTrack(mcTrack->GetMother()));
          if (mcMother)
          {
            fBuffer_ClusterCandidate_MC_Mother_E[fBuffer_NClusterCandidates] = mcMother->E();
            fBuffer_ClusterCandidate_MC_Mother_Px[fBuffer_NClusterCandidates] = mcMother->Px();
            fBuffer_ClusterCandidate_MC_Mother_Py[fBuffer_NClusterCandidates] = mcMother->Py();
            fBuffer_ClusterCandidate_MC_Mother_Pz[fBuffer_NClusterCandidates] = mcMother->Pz();
            fBuffer_ClusterCandidate_MC_Mother_PDG[fBuffer_NClusterCandidates] = mcMother->GetPdgCode();

            AliAODMCParticle *mcGrandMother = dynamic_cast<AliAODMCParticle *>(fMCEvent->GetTrack(mcMother->GetMother()));
            if (mcGrandMother)
            {
              // fBuffer_ClusterCandidate_MC_GrandMother_E[fBuffer_NClusterCandidates] = mcGrandMother->E();
              // fBuffer_ClusterCandidate_MC_GrandMother_Px[fBuffer_NClusterCandidates] = mcGrandMother->Px();
              // fBuffer_ClusterCandidate_MC_GrandMother_Py[fBuffer_NClusterCandidates] = mcGrandMother->Py();
              // fBuffer_ClusterCandidate_MC_GrandMother_Pz[fBuffer_NClusterCandidates] = mcGrandMother->Pz();
              fBuffer_ClusterCandidate_MC_GrandMother_PDG[fBuffer_NClusterCandidates] = mcGrandMother->GetPdgCode();
            }
          }
        }
      }
    }

    fBuffer_NClusterCandidates++;
    if (kMaxTracks <= fBuffer_NClusterCandidates)
    {
      delete clus;
      AliFatal("Too many clusters");
    }
    // delete mutrack;
    delete clus;
  }
}

// ///________________________________________________________________________
// void AliAnalysisTaskCorrelationTree::ProcessQATree(AliAODConversionPhoton *gamma)
// {

//   // Fill Histograms for QA and MC
//   AliVEvent *event = (AliVEvent *)InputEvent();

//   fBuffer_ElectronCandidate_E[fBuffer_NElectronCandidates] = gamma->GetPhotonP();
//   fBuffer_ElectronCandidate_Px[fBuffer_NElectronCandidates] = gamma->GetPx();
//   fBuffer_ElectronCandidate_Py[fBuffer_NElectronCandidates] = gamma->GetPy();
//   fBuffer_ElectronCandidate_Pz[fBuffer_NElectronCandidates] = gamma->GetPz();
//   fBuffer_ElectronCandidate_NSigmaElecTPC[fBuffer_NElectronCandidates] = gamma->GetArmenterosQt();
//   fBuffer_ElectronCandidate_NSigmaElecTOF[fBuffer_NElectronCandidates] = gamma->GetArmenterosAlpha();
//   fBuffer_ElectronCandidate_PsiPair[fBuffer_NElectronCandidates] = gamma->GetPsiPair();
//   fBuffer_ElectronCandidate_Chi2[fBuffer_NElectronCandidates] = gamma->GetChi2perNDF();
//   fBuffer_ElectronCandidate_CosPA[fBuffer_NElectronCandidates] = fConversionCuts->GetCosineOfPointingAngle(gamma, event);
//   fBuffer_ElectronCandidate_Eta[fBuffer_NElectronCandidates] = gamma->GetPhotonEta();
//   fBuffer_ElectronCandidate_Phi[fBuffer_NElectronCandidates] = gamma->GetPhotonPhi();
//   fBuffer_ElectronCandidate_ConvPointX[fBuffer_NElectronCandidates] = gamma->GetConversionX();
//   fBuffer_ElectronCandidate_ConvPointY[fBuffer_NElectronCandidates] = gamma->GetConversionY();
//   fBuffer_ElectronCandidate_ConvPointZ[fBuffer_NElectronCandidates] = gamma->GetConversionZ();

//   fBuffer_ElectronCandidate_MC_E[fBuffer_NElectronCandidates] = 9;
//   fBuffer_ElectronCandidate_MC_Mother_Type[fBuffer_NElectronCandidates] = 9;
//   if (fMCEvent && fInputEvent->IsA() == AliESDEvent::Class())
//   {
//     fBuffer_ElectronCandidate_MC_E[fBuffer_NElectronCandidates] = IsTruePhotonESD(gamma);
//     fBuffer_ElectronCandidate_MC_Mother_Type[fBuffer_NElectronCandidates] = GetTrueMotherInfoESD(gamma);
//   }
//   else if (fMCEvent && fInputEvent->IsA() == AliAODEvent::Class())
//   {
//     fBuffer_ElectronCandidate_MC_E[fBuffer_NElectronCandidates] = IsTruePhotonAOD(gamma);
//     fBuffer_ElectronCandidate_MC_Mother_Type[fBuffer_NElectronCandidates] = GetTrueMotherInfoAOD(gamma);
//   }
//   fBuffer_NElectronCandidates++;
// }

// //________________________________________________________________________
// UInt_t AliAnalysisTaskCorrelationTree::IsTruePhotonAOD(AliAODConversionPhoton *TruePhotonCandidate)
// {

//   UInt_t kind = 9;
//   TClonesArray *AODMCTrackArray = dynamic_cast<TClonesArray *>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
//   if (AODMCTrackArray != NULL && TruePhotonCandidate != NULL)
//   {
//     AliAODMCParticle *posDaughter = (AliAODMCParticle *)AODMCTrackArray->At(TruePhotonCandidate->GetMCLabelPositive());
//     AliAODMCParticle *negDaughter = (AliAODMCParticle *)AODMCTrackArray->At(TruePhotonCandidate->GetMCLabelNegative());
//     Int_t pdgCodePos = 0;
//     Int_t pdgCodeNeg = 0;
//     Int_t pdgCode = 0;
//     if (posDaughter == NULL || negDaughter == NULL)
//     {
//       kind = 9;
//     }
//     else if (posDaughter->GetMother() != negDaughter->GetMother() || (posDaughter->GetMother() == negDaughter->GetMother() && posDaughter->GetMother() == -1))
//     {
//       kind = 1;
//       pdgCodePos = TMath::Abs(posDaughter->GetPdgCode());
//       pdgCodeNeg = TMath::Abs(negDaughter->GetPdgCode());
//       if (pdgCodePos == 11 && pdgCodeNeg == 11)
//         kind = 10; // Electron Combinatorial
//       if (pdgCodePos == 11 && pdgCodeNeg == 11 &&
//           (posDaughter->GetMother() == negDaughter->GetMother() && posDaughter->GetMother() == -1))
//         kind = 15; // direct Electron Combinatorial

//       if (pdgCodePos == 211 && pdgCodeNeg == 211)
//         kind = 11; // Pion Combinatorial
//       if ((pdgCodePos == 211 && pdgCodeNeg == 2212) || (pdgCodePos == 2212 && pdgCodeNeg == 211))
//         kind = 12; // Pion, Proton Combinatorics
//       if ((pdgCodePos == 11 && pdgCodeNeg == 2212) || (pdgCodePos == 2212 && pdgCodeNeg == 11))
//         kind = 16; // electron, Proton Combinatorics
//       if ((pdgCodePos == 11 && pdgCodeNeg == 321) || (pdgCodePos == 321 && pdgCodeNeg == 11))
//         kind = 17; // electron, kaon
//       if ((pdgCodePos == 211 && pdgCodeNeg == 321) || (pdgCodePos == 321 && pdgCodeNeg == 211))
//         kind = 18; // pion, kaon
//       if ((pdgCodePos == 211 && pdgCodeNeg == 11) || (pdgCodePos == 11 && pdgCodeNeg == 211))
//         kind = 13; // Pion, Electron Combinatorics
//       if (pdgCodePos == 321 && pdgCodeNeg == 321)
//         kind = 14; // Kaon,Kaon combinatorics
//     }
//     else
//     {
//       AliAODMCParticle *Photon = (AliAODMCParticle *)AODMCTrackArray->At(posDaughter->GetMother());
//       pdgCodePos = posDaughter->GetPdgCode();
//       pdgCodeNeg = negDaughter->GetPdgCode();

//       if (Photon->GetPdgCode())
//         pdgCode = Photon->GetPdgCode();
//       if (TMath::Abs(pdgCodePos) != 11 || TMath::Abs(pdgCodeNeg) != 11)
//         kind = 2; // true from hadronic decays
//       else if (!(pdgCodeNeg == pdgCodePos))
//       {
//         if (pdgCode == 111)
//           kind = 3; // pi0 Dalitz
//         else if (pdgCode == 221)
//           kind = 4; // eta Dalitz
//         else if (!(negDaughter->GetMCProcessCode() != 5 || posDaughter->GetMCProcessCode() != 5))
//         {
//           const AliVVertex *primVtxMC = fMCEvent->GetPrimaryVertex();
//           double mcProdVtxX = primVtxMC->GetX();
//           double mcProdVtxY = primVtxMC->GetY();
//           double mcProdVtxZ = primVtxMC->GetZ();
//           Bool_t isPrimary = fEventCuts->IsConversionPrimaryAOD(fInputEvent, Photon, mcProdVtxX, mcProdVtxY, mcProdVtxZ);

//           if (pdgCode == 22 && isPrimary)
//           {
//             kind = 0; // primary photons
//           }
//           else if (pdgCode == 22)
//           {
//             kind = 5; // secondary photons
//           }
//         }
//       }
//     }

//     return kind;
//   }
//   return kind;
// }

//________________________________________________________________________
void AliAnalysisTaskCorrelationTree::Terminate(Option_t *)
{
}

void AliAnalysisTaskCorrelationTree::ResetBuffer()
{

  fBuffer_NEventTriggers = 0;
  fBuffer_NContributors = 0;
  fBuffer_RunNumber = 0;
  fBuffer_VertexZ = 0;
  for (Int_t ntrig = 0; ntrig < kMaxTriggers; ntrig++)
  {
    fBuffer_EventTrigger[ntrig] = 0;
  }

  fBuffer_NElectronCandidates = 0;
  fBuffer_NMuonCandidates = 0;
  fBuffer_NClusterCandidates = 0;

  for (Int_t ccand = 0; ccand < kMaxTracks; ccand++)
  {
    fBuffer_ElectronCandidate_E[ccand] = 0;
    fBuffer_ElectronCandidate_Px[ccand] = 0;
    fBuffer_ElectronCandidate_Py[ccand] = 0;
    fBuffer_ElectronCandidate_Pz[ccand] = 0;
    fBuffer_ElectronCandidate_PropEta[ccand] = 0;
    fBuffer_ElectronCandidate_PropPhi[ccand] = 0;
    fBuffer_ElectronCandidate_Charge[ccand] = 0;
    fBuffer_ElectronCandidate_NSigmaElecTPC[ccand] = 0;
    fBuffer_ElectronCandidate_NSigmaElecTOF[ccand] = 0;
    fBuffer_ElectronCandidate_NSigmaPionTPC[ccand] = 0;
    fBuffer_ElectronCandidate_NSigmaPionTOF[ccand] = 0;
    fBuffer_ElectronCandidate_NSigmaKaonTPC[ccand] = 0;
    fBuffer_ElectronCandidate_NSigmaProtonTPC[ccand] = 0;

    if (fIsMC)
    {
      fBuffer_ElectronCandidate_MC_E[ccand] = 0;
      fBuffer_ElectronCandidate_MC_Px[ccand] = 0;
      fBuffer_ElectronCandidate_MC_Py[ccand] = 0;
      fBuffer_ElectronCandidate_MC_Pz[ccand] = 0;
      fBuffer_ElectronCandidate_MC_PDG[ccand] = 0;
      fBuffer_ElectronCandidate_MC_Mother_E[ccand] = 0;
      fBuffer_ElectronCandidate_MC_Mother_Px[ccand] = 0;
      fBuffer_ElectronCandidate_MC_Mother_Py[ccand] = 0;
      fBuffer_ElectronCandidate_MC_Mother_Pz[ccand] = 0;
      fBuffer_ElectronCandidate_MC_Mother_PDG[ccand] = 0;
      fBuffer_ElectronCandidate_MC_GrandMother_PDG[ccand] = 0;
      fBuffer_ElectronCandidate_MC_GrandMother_Pt[ccand] = 0;
      fBuffer_ElectronCandidate_MC_GGrandMother_PDG[ccand] = 0;
      fBuffer_ElectronCandidate_MC_GGrandMother_Pt[ccand] = 0;
    }
    fBuffer_MuonCandidate_E[ccand] = 0;
    fBuffer_MuonCandidate_Px[ccand] = 0;
    fBuffer_MuonCandidate_Py[ccand] = 0;
    fBuffer_MuonCandidate_Pz[ccand] = 0;
    fBuffer_MuonCandidate_RAbsEnd[ccand] = 0;
    fBuffer_MuonCandidate_Chi2NDF[ccand] = 0;
    fBuffer_MuonCandidate_DCA[ccand] = 0;
    fBuffer_MuonCandidate_Charge[ccand] = 0;
    fBuffer_MuonCandidate_MatchTrigger[ccand] = 0;

    if (fIsMC)
    {
      fBuffer_MuonCandidate_MC_E[ccand] = 0;
      fBuffer_MuonCandidate_MC_Px[ccand] = 0;
      fBuffer_MuonCandidate_MC_Py[ccand] = 0;
      fBuffer_MuonCandidate_MC_Pz[ccand] = 0;
      fBuffer_MuonCandidate_MC_PDG[ccand] = 0;
      fBuffer_MuonCandidate_MC_Mother_E[ccand] = 0;
      fBuffer_MuonCandidate_MC_Mother_Px[ccand] = 0;
      fBuffer_MuonCandidate_MC_Mother_Py[ccand] = 0;
      fBuffer_MuonCandidate_MC_Mother_Pz[ccand] = 0;
      fBuffer_MuonCandidate_MC_Mother_PDG[ccand] = 0;
      fBuffer_MuonCandidate_MC_GrandMother_PDG[ccand] = 0;
      fBuffer_MuonCandidate_MC_GrandMother_Pt[ccand] = 0;
      fBuffer_MuonCandidate_MC_GGrandMother_PDG[ccand] = 0;
      fBuffer_MuonCandidate_MC_GGrandMother_Pt[ccand] = 0;
    }
    fBuffer_ClusterCandidate_E[ccand] = 0;
    fBuffer_ClusterCandidate_Px[ccand] = 0;
    fBuffer_ClusterCandidate_Py[ccand] = 0;
    fBuffer_ClusterCandidate_Pz[ccand] = 0;
    fBuffer_ClusterCandidate_Eta[ccand] = 0;
    fBuffer_ClusterCandidate_Phi[ccand] = 0;
    fBuffer_ClusterCandidate_M02[ccand] = -10;

    if (fIsMC)
    {
      fBuffer_ClusterCandidate_MC_E[ccand] = 0;
      fBuffer_ClusterCandidate_MC_Px[ccand] = 0;
      fBuffer_ClusterCandidate_MC_Py[ccand] = 0;
      fBuffer_ClusterCandidate_MC_Pz[ccand] = 0;
      fBuffer_ClusterCandidate_MC_PDG[ccand] = 0;
      fBuffer_ClusterCandidate_MC_Mother_E[ccand] = 0;
      fBuffer_ClusterCandidate_MC_Mother_Px[ccand] = 0;
      fBuffer_ClusterCandidate_MC_Mother_Py[ccand] = 0;
      fBuffer_ClusterCandidate_MC_Mother_Pz[ccand] = 0;
      fBuffer_ClusterCandidate_MC_Mother_PDG[ccand] = 0;
      fBuffer_ClusterCandidate_MC_GrandMother_PDG[ccand] = 0;
    }
  }
}
