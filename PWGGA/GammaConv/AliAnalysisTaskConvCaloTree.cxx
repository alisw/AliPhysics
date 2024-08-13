/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Authors: Joshua Koenig                                              *
 * Version 1.0                                                            *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/////////////////////////////////////////////////////////////////
//---------------------------------------------------------------
// Class used to create a tree for meutral meson analysis studies
//---------------------------------------------------------------
/////////////////////////////////////////////////////////////////

#include "AliAnalysisTaskConvCaloTree.h"
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
#include "AliMultSelection.h"
#include "AliGenPythiaEventHeader.h"
#include "AliEmcalTriggerDecisionContainer.h"
#include <bitset>

class iostream;

using namespace std;

ClassImp(AliAnalysisTaskConvCaloTree)

  //________________________________________________________________________
  AliAnalysisTaskConvCaloTree::AliAnalysisTaskConvCaloTree() : AliAnalysisTaskSE(),
                                                               fV0Reader(NULL),
                                                               fV0ReaderName("V0ReaderV1"),
                                                               fAddNameConvJet("FullJet"),
                                                               fConvJetReader(NULL),
                                                               fCaloIsolation(NULL),
                                                               fCaloIsolationName("PhotonIsolation"),
                                                               fReaderGammas(NULL),
                                                               fPIDResponse(NULL),
                                                               fCorrTaskSetting(""),
                                                               fInputEvent(NULL),
                                                               fMCEvent(NULL),
                                                               fWeightJetJetMC(1),
                                                               fEventCuts(NULL),
                                                               fAliEventCuts(false),
                                                               fClusterCutsEMC(NULL),
                                                               fClusterCutsPHOS(NULL),
                                                               fClusterCuts(NULL),
                                                               fConversionCuts(NULL),
                                                               fMesonCuts(NULL),
                                                               fPhotonTree(NULL),
                                                               fIsHeavyIon(kFALSE),
                                                               fAODMCTrackArray(NULL),
                                                               farrClustersProcess(NULL),
                                                               fOutputList(NULL),
                                                               fIsMC(0),
                                                               fMCEventPos(),
                                                               fMCEventNeg(),
                                                               fESDArrayPos(),
                                                               fESDArrayNeg(),
                                                               fSaveMCInformation(0),
                                                               fSaveClusters(0),
                                                               fSaveConversions(0),
                                                               fSaveTracks(0),
                                                               fStoreTriggerCondition(0),
                                                               fUseClusterIsolation(0),
                                                               fSaveJets(0),
                                                               fDoDownscale(false),
                                                               fDownscaleFac(1.),
                                                               fMinTrackPt(0),
                                                               fInfoList(nullptr),
                                                               fHistoNEvents(nullptr),
                                                               fHistoMCGammaPt(nullptr),
                                                               fHistoMCElecPt(nullptr),
                                                               fHistoMCPi0Pt(nullptr),
                                                               fHistoMCEtaPt(nullptr),
                                                               fHistoMCOmegaPt(nullptr),
                                                               fHistoMCEtaPrimePt(nullptr),
                                                               fHistoMCK0sPt(nullptr),
                                                               fHistoMCSecPi0Pt(nullptr),
                                                               fHistoMCSecEtaPt(nullptr),
                                                               fBuffer_EventWeight(1),
                                                               fBuffer_Event_Vertex_Z(0),
                                                               fMCQ2(0.),
                                                               fV0Mult(0.),
                                                               fTriggerBits(0),
                                                               fVBuffer_Cluster_E(),
                                                               fVBuffer_Cluster_Eta(),
                                                               fVBuffer_Cluster_Phi(),
                                                               fVBuffer_Cluster_NCells(),
                                                               fVBuffer_Cluster_M02(),
                                                               fVBuffer_Cluster_Fcross(),
                                                               fVBuffer_Cluster_Time(),
                                                               fVBuffer_Cluster_Exotic(),
                                                               fVBuffer_Cluster_Isolated(),
                                                               fVBuffer_Cluster_CellTimes(),
                                                               fVBuffer_Cluster_CellEnergies(),
                                                               fVBuffer_Cluster_CellIDs(),
                                                               fVTrueClusterPi0DaughterIndex(),
                                                               fVTrueClusterEtaDaughterIndex(),
                                                               fVTrueClusterSecondaryIndex(),
                                                               fVTrueClusterMCId(),
                                                               fVTrueClusterMCId2(),
                                                               fVTrueClusterMCTrueEnergy(),
                                                               fVTrueClusterMCTrueEnergy2(),
                                                               fVTrueClusterMCIsMerged(),
                                                               fVTrueClusterConvRadius(),
                                                               fVBuffer_Conv_px(),
                                                               fVBuffer_Conv_py(),
                                                               fVBuffer_Conv_pz(),
                                                               fVBuffer_Elec1etaCalo(),
                                                               fVBuffer_Elec2etaCalo(),
                                                               fVBuffer_Elec1phiCalo(),
                                                               fVBuffer_Elec2phiCalo(),
                                                               fVBuffer_Conv_R(),
                                                               fVBuffer_Conv_PsiPair(),
                                                               fVBuffer_Conv_NTPCClusElec1(),
                                                               fVBuffer_Conv_NTPCClusElec2(),
                                                               fVBuffer_Conv_dEdxElec1(),
                                                               fVBuffer_Conv_dEdxElec2(),
                                                               fVBuffer_Conv_PElec1(),
                                                               fVBuffer_Conv_PElec2(),
                                                               fVBuffer_Conv_Quality(),
                                                               fVBuffer_Conv_CosPAngle(),
                                                               fVBuffer_Conv_Chi2(),
                                                               fVTrueConvPi0DaughterIndex(),
                                                               fVTrueConvEtaDaughterIndex(),
                                                               fVTrueConvSecDaughterIndex(),
                                                               fVTrueConvMCTruePx(),
                                                               fVTrueConvMCTruePy(),
                                                               fVTrueConvMCTruePz(),
                                                               fVTrueConvMCLabel(),
                                                               fVBuffer_Track_px(),
                                                               fVBuffer_Track_py(),
                                                               fVBuffer_Track_pz(),
                                                               fVBuffer_Track_P(),
                                                               fVBuffer_Track_dedx(),
                                                               fVBuffer_Track_TOFSignal(),
                                                               fVBuffer_Track_DCA(),
                                                               fVBuffer_Track_FracNClus(),
                                                               fVBuffer_Track_PDG(),
                                                               fVBuffer_Track_StackID(),
                                                               fVBuffer_Track_Calo_eta(),
                                                               fVBuffer_Track_Calo_phi(),
                                                               fVBuffer_Jet_Pt(),
                                                               fVBuffer_Jet_Eta(),
                                                               fVBuffer_Jet_Phi(),
                                                               fVBuffer_Jet_NNeutr(),
                                                               fVBuffer_Jet_NCh(),
                                                               fVBuffer_Jet_NEF(),
                                                               fVBuffer_TrueJet_Pt(),
                                                               fVBuffer_TrueJet_Eta(),
                                                               fVBuffer_TrueJet_Phi(),
                                                               fVBuffer_TrueJet_Parton_Pt(),
                                                               fVBuffer_MCGenID(),
                                                               fVBuffer_MCGenPDG(),
                                                               fVBuffer_MCGenPx(),
                                                               fVBuffer_MCGenPy(),
                                                               fVBuffer_MCGenPz(),
                                                               fVBuffer_MCGenMotherID()
{
}

AliAnalysisTaskConvCaloTree::AliAnalysisTaskConvCaloTree(const char* name) : AliAnalysisTaskSE(name),
                                                                             fV0Reader(NULL),
                                                                             fV0ReaderName("V0ReaderV1"),
                                                                             fAddNameConvJet("FullJet"),
                                                                             fConvJetReader(NULL),
                                                                             fCaloIsolation(NULL),
                                                                             fCaloIsolationName("PhotonIsolation"),
                                                                             fReaderGammas(NULL),
                                                                             fPIDResponse(NULL),
                                                                             fCorrTaskSetting(""),
                                                                             fInputEvent(NULL),
                                                                             fMCEvent(NULL),
                                                                             fWeightJetJetMC(1),
                                                                             fEventCuts(NULL),
                                                                             fAliEventCuts(false),
                                                                             fClusterCutsEMC(NULL),
                                                                             fClusterCutsPHOS(NULL),
                                                                             fClusterCuts(NULL),
                                                                             fConversionCuts(NULL),
                                                                             fMesonCuts(NULL),
                                                                             fPhotonTree(NULL),
                                                                             fIsHeavyIon(kFALSE),
                                                                             fAODMCTrackArray(NULL),
                                                                             farrClustersProcess(NULL),
                                                                             fOutputList(NULL),
                                                                             fIsMC(0),
                                                                             fMCEventPos(),
                                                                             fMCEventNeg(),
                                                                             fESDArrayPos(),
                                                                             fESDArrayNeg(),
                                                                             fSaveMCInformation(0),
                                                                             fSaveClusters(0),
                                                                             fSaveConversions(0),
                                                                             fSaveTracks(0),
                                                                             fStoreTriggerCondition(0),
                                                                             fUseClusterIsolation(0),
                                                                             fSaveJets(0),
                                                                             fDoDownscale(false),
                                                                             fDownscaleFac(1.),
                                                                             fMinTrackPt(0),
                                                                             fInfoList(nullptr),
                                                                             fHistoNEvents(nullptr),
                                                                             fHistoMCGammaPt(nullptr),
                                                                             fHistoMCElecPt(nullptr),
                                                                             fHistoMCPi0Pt(nullptr),
                                                                             fHistoMCEtaPt(nullptr),
                                                                             fHistoMCOmegaPt(nullptr),
                                                                             fHistoMCEtaPrimePt(nullptr),
                                                                             fHistoMCK0sPt(nullptr),
                                                                             fHistoMCSecPi0Pt(nullptr),
                                                                             fHistoMCSecEtaPt(nullptr),
                                                                             fBuffer_EventWeight(1),
                                                                             fBuffer_Event_Vertex_Z(0),
                                                                             fMCQ2(0.),
                                                                             fV0Mult(0.),
                                                                             fTriggerBits(0),
                                                                             fVBuffer_Cluster_E(),
                                                                             fVBuffer_Cluster_Eta(),
                                                                             fVBuffer_Cluster_Phi(),
                                                                             fVBuffer_Cluster_NCells(),
                                                                             fVBuffer_Cluster_M02(),
                                                                             fVBuffer_Cluster_Fcross(),
                                                                             fVBuffer_Cluster_Time(),
                                                                             fVBuffer_Cluster_Exotic(),
                                                                             fVBuffer_Cluster_Isolated(),
                                                                             fVBuffer_Cluster_CellTimes(),
                                                                             fVBuffer_Cluster_CellEnergies(),
                                                                             fVBuffer_Cluster_CellIDs(),
                                                                             fVTrueClusterPi0DaughterIndex(),
                                                                             fVTrueClusterEtaDaughterIndex(),
                                                                             fVTrueClusterSecondaryIndex(),
                                                                             fVTrueClusterMCId(),
                                                                             fVTrueClusterMCId2(),
                                                                             fVTrueClusterMCTrueEnergy(),
                                                                             fVTrueClusterMCTrueEnergy2(),
                                                                             fVTrueClusterMCIsMerged(),
                                                                             fVTrueClusterConvRadius(),
                                                                             fVBuffer_Conv_px(),
                                                                             fVBuffer_Conv_py(),
                                                                             fVBuffer_Conv_pz(),
                                                                             fVBuffer_Elec1etaCalo(),
                                                                             fVBuffer_Elec2etaCalo(),
                                                                             fVBuffer_Elec1phiCalo(),
                                                                             fVBuffer_Elec2phiCalo(),
                                                                             fVBuffer_Conv_R(),
                                                                             fVBuffer_Conv_PsiPair(),
                                                                             fVBuffer_Conv_NTPCClusElec1(),
                                                                             fVBuffer_Conv_NTPCClusElec2(),
                                                                             fVBuffer_Conv_dEdxElec1(),
                                                                             fVBuffer_Conv_dEdxElec2(),
                                                                             fVBuffer_Conv_PElec1(),
                                                                             fVBuffer_Conv_PElec2(),
                                                                             fVBuffer_Conv_Quality(),
                                                                             fVBuffer_Conv_CosPAngle(),
                                                                             fVBuffer_Conv_Chi2(),
                                                                             fVTrueConvPi0DaughterIndex(),
                                                                             fVTrueConvEtaDaughterIndex(),
                                                                             fVTrueConvSecDaughterIndex(),
                                                                             fVTrueConvMCTruePx(),
                                                                             fVTrueConvMCTruePy(),
                                                                             fVTrueConvMCTruePz(),
                                                                             fVTrueConvMCLabel(),
                                                                             fVBuffer_Track_px(),
                                                                             fVBuffer_Track_py(),
                                                                             fVBuffer_Track_pz(),
                                                                             fVBuffer_Track_P(),
                                                                             fVBuffer_Track_dedx(),
                                                                             fVBuffer_Track_TOFSignal(),
                                                                             fVBuffer_Track_DCA(),
                                                                             fVBuffer_Track_FracNClus(),
                                                                             fVBuffer_Track_PDG(),
                                                                             fVBuffer_Track_StackID(),
                                                                             fVBuffer_Track_Calo_eta(),
                                                                             fVBuffer_Track_Calo_phi(),
                                                                             fVBuffer_Jet_Pt(),
                                                                             fVBuffer_Jet_Eta(),
                                                                             fVBuffer_Jet_Phi(),
                                                                             fVBuffer_Jet_NNeutr(),
                                                                             fVBuffer_Jet_NCh(),
                                                                             fVBuffer_Jet_NEF(),
                                                                             fVBuffer_TrueJet_Pt(),
                                                                             fVBuffer_TrueJet_Eta(),
                                                                             fVBuffer_TrueJet_Phi(),
                                                                             fVBuffer_TrueJet_Parton_Pt(),
                                                                             fVBuffer_MCGenID(),
                                                                             fVBuffer_MCGenPDG(),
                                                                             fVBuffer_MCGenPx(),
                                                                             fVBuffer_MCGenPy(),
                                                                             fVBuffer_MCGenPz(),
                                                                             fVBuffer_MCGenMotherID()
{
  // Do not perform trigger selection in the AliEvent cuts but let the task do this before
  fAliEventCuts.OverrideAutomaticTriggerSelection(AliVEvent::kAny, true);

  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
  DefineOutput(2, TTree::Class());
}

//________________________________________________________________________
AliAnalysisTaskConvCaloTree::~AliAnalysisTaskConvCaloTree()
{
  // default deconstructor
}

//________________________________________________________________________
void AliAnalysisTaskConvCaloTree::UserCreateOutputObjects()
{
  fV0Reader = (AliV0ReaderV1*)AliAnalysisManager::GetAnalysisManager()->GetTask(fV0ReaderName.Data());
  if (!fV0Reader) {
    printf("Error: No V0 Reader");
    return;
  } // GetV0Reader

  if (fSaveJets) {
    fConvJetReader = (AliAnalysisTaskConvJet*)AliAnalysisManager::GetAnalysisManager()->GetTask(Form("AliAnalysisTaskConvJet%s", fAddNameConvJet.EqualTo("") == true ? "" : Form("_%s", fAddNameConvJet.Data())));
    if (!fConvJetReader) {
      AliFatal(Form("ERROR: No AliAnalysisTaskConvJet%s\n", fAddNameConvJet.EqualTo("") == true ? "" : Form("_%s", fAddNameConvJet.Data())));
    }
  }

  // retrieve pointer to CaloIsolation Instance
  if (fUseClusterIsolation)
    fCaloIsolation = (AliPhotonIsolation*)AliAnalysisManager::GetAnalysisManager()->GetTask(fCaloIsolationName.Data());
  if (!fCaloIsolation && fUseClusterIsolation) {
    AliFatal("AliPhotonIsolation instance could not be initialized!");
  }

  if (fOutputList != NULL) {
    delete fOutputList;
    fOutputList = NULL;
  }
  if (fOutputList == NULL) {
    fOutputList = new TList();
    fOutputList->SetOwner(kTRUE);
  }

  if (((AliConvEventCuts*)fEventCuts)->GetCutHistograms()) {
    fOutputList->Add(((AliConvEventCuts*)fEventCuts)->GetCutHistograms());
  }
  if (fClusterCutsEMC) {
    if (((AliCaloPhotonCuts*)fClusterCutsEMC)->GetCutHistograms()) {
      fOutputList->Add(((AliCaloPhotonCuts*)fClusterCutsEMC)->GetCutHistograms());
    }
  }
  if (fClusterCutsPHOS) {
    if (((AliCaloPhotonCuts*)fClusterCutsPHOS)->GetCutHistograms()) {
      fOutputList->Add(((AliCaloPhotonCuts*)fClusterCutsPHOS)->GetCutHistograms());
    }
  }

  fInfoList = new TList();
  fInfoList->SetName(Form("%s event histograms", (fEventCuts->GetCutNumber()).Data()));
  fInfoList->SetOwner(kTRUE);

  if (fInfoList) {
    fOutputList->Add(fInfoList);
  }
  PostData(1, fOutputList);

  fPhotonTree = new TTree(Form("ConvCaloPhotons_%s", (fEventCuts->GetCutNumber()).Data()), Form("ConvCaloPhotons_%s", (fEventCuts->GetCutNumber()).Data()));

  if (fIsMC > 1) {
    fPhotonTree->Branch("Event_Weight", &fBuffer_EventWeight, "Event_Weight/F");
  }
  if (fIsMC) {
    fPhotonTree->Branch("Event_Q2", &fMCQ2);
  }
  fPhotonTree->Branch("Event_VertexZ", &fBuffer_Event_Vertex_Z, "Event_VertexZ/F");
  fPhotonTree->Branch("Event_V0MultPerc", &fV0Mult);
  if (fStoreTriggerCondition) {
    fPhotonTree->Branch("Event_TriggerBits", &fTriggerBits);
  }

  if (fSaveClusters) {

    fPhotonTree->Branch("Cluster_E", &fVBuffer_Cluster_E);
    fPhotonTree->Branch("Cluster_Eta", &fVBuffer_Cluster_Eta);
    fPhotonTree->Branch("Cluster_Phi", &fVBuffer_Cluster_Phi);
    fPhotonTree->Branch("Cluster_NCells", &fVBuffer_Cluster_NCells);
    fPhotonTree->Branch("Cluster_M02", &fVBuffer_Cluster_M02);
    fPhotonTree->Branch("Cluster_Time", &fVBuffer_Cluster_Time);
    if (fClusterCutsEMC) {
      fPhotonTree->Branch("Cluster_Fcross", &fVBuffer_Cluster_Fcross);
      fPhotonTree->Branch("Cluster_isExotic", &fVBuffer_Cluster_Exotic);
      fPhotonTree->Branch("Cluster_CellEnergies", &fVBuffer_Cluster_CellEnergies);
      fPhotonTree->Branch("Cluster_CellIDs", &fVBuffer_Cluster_CellIDs);
      if (fIsMC == 0) {
        fPhotonTree->Branch("Cluster_CellTimes", &fVBuffer_Cluster_CellTimes);
      }
    }
    if (fUseClusterIsolation) {
      fPhotonTree->Branch("Cluster_Isolated", &fVBuffer_Cluster_Isolated);
    }
    if (fIsMC) {
      fPhotonTree->Branch("Cluster_MotherPi0", &fVTrueClusterPi0DaughterIndex);
      fPhotonTree->Branch("Cluster_MotherEta", &fVTrueClusterEtaDaughterIndex);
      fPhotonTree->Branch("Cluster_MotherSec", &fVTrueClusterSecondaryIndex);
      fPhotonTree->Branch("Cluster_MCId", &fVTrueClusterMCId);
      fPhotonTree->Branch("Cluster_MCId2", &fVTrueClusterMCId2);
      fPhotonTree->Branch("Cluster_MCTrueEnergy", &fVTrueClusterMCTrueEnergy);
      fPhotonTree->Branch("Cluster_MCTrueEnergy2", &fVTrueClusterMCTrueEnergy2);
      fPhotonTree->Branch("Cluster_MCTrueIsMerged", &fVTrueClusterMCIsMerged);
      fPhotonTree->Branch("Cluster_MCTrueConvR", &fVTrueClusterConvRadius);
    }
  }
  if (fSaveConversions) {
    fPhotonTree->Branch("Conv_px", &fVBuffer_Conv_px);
    fPhotonTree->Branch("Conv_py", &fVBuffer_Conv_py);
    fPhotonTree->Branch("Conv_pz", &fVBuffer_Conv_pz);
    fPhotonTree->Branch("Conv_Elec1etaCalo", &fVBuffer_Elec1etaCalo);
    fPhotonTree->Branch("Conv_Elec2etaCalo", &fVBuffer_Elec2etaCalo);
    fPhotonTree->Branch("Conv_Elec1phiCalo", &fVBuffer_Elec1phiCalo);
    fPhotonTree->Branch("Conv_Elec2phiCalo", &fVBuffer_Elec2phiCalo);
    if (fSaveConversions > 1) {
      fPhotonTree->Branch("Conv_Radius", &fVBuffer_Conv_R);
      fPhotonTree->Branch("Conv_PsiPair", &fVBuffer_Conv_PsiPair);
      fPhotonTree->Branch("Conv_NTPCClus_Elec1", &fVBuffer_Conv_NTPCClusElec1);
      fPhotonTree->Branch("Conv_NTPCClus_Elec2", &fVBuffer_Conv_NTPCClusElec2);
      fPhotonTree->Branch("Conv_dEdx_Elec1", &fVBuffer_Conv_dEdxElec1);
      fPhotonTree->Branch("Conv_dEdx_Elec2", &fVBuffer_Conv_dEdxElec2);
      fPhotonTree->Branch("Conv_P_Elec1", &fVBuffer_Conv_PElec1);
      fPhotonTree->Branch("Conv_P_Elec2", &fVBuffer_Conv_PElec2);
      fPhotonTree->Branch("Conv_Quality", &fVBuffer_Conv_Quality);
      fPhotonTree->Branch("Conv_CosPAngle", &fVBuffer_Conv_CosPAngle);
      fPhotonTree->Branch("Conv_Chi2NDF", &fVBuffer_Conv_Chi2);
    }
    if (fIsMC) {
      fPhotonTree->Branch("Conv_MotherPi0", &fVTrueConvPi0DaughterIndex);
      fPhotonTree->Branch("Conv_MotherEta", &fVTrueConvEtaDaughterIndex);
      fPhotonTree->Branch("Conv_MotherSec", &fVTrueConvSecDaughterIndex);
      fPhotonTree->Branch("Conv_MCTruePx", &fVTrueConvMCTruePx);
      fPhotonTree->Branch("Conv_MCTruePy", &fVTrueConvMCTruePy);
      fPhotonTree->Branch("Conv_MCTruePz", &fVTrueConvMCTruePz);
      fPhotonTree->Branch("Conv_MCLabel", &fVTrueConvMCLabel);
    }
  }
  if (fSaveTracks) {
    if (fSaveTracks == 1) {
      fPhotonTree->Branch("Track_P", &fVBuffer_Track_P);
    } else if (fSaveTracks == 2) {
      fPhotonTree->Branch("Track_px", &fVBuffer_Track_px);
      fPhotonTree->Branch("Track_py", &fVBuffer_Track_py);
      fPhotonTree->Branch("Track_pz", &fVBuffer_Track_pz);
      fPhotonTree->Branch("Track_dEdx", &fVBuffer_Track_dedx);
      fPhotonTree->Branch("Track_TOFSignal", &fVBuffer_Track_TOFSignal);
      fPhotonTree->Branch("Track_DCA", &fVBuffer_Track_DCA);
      fPhotonTree->Branch("Track_FracNClus_Chi2", &fVBuffer_Track_FracNClus);
      if (fIsMC) {
        fPhotonTree->Branch("Track_PDG", &fVBuffer_Track_PDG);
        fPhotonTree->Branch("Track_StackID", &fVBuffer_Track_StackID);
      }
    }
    fPhotonTree->Branch("Track_TracketaCalo", &fVBuffer_Track_Calo_eta);
    fPhotonTree->Branch("Track_TrackphiCalo", &fVBuffer_Track_Calo_phi);
  }
  if (fSaveJets) {
    fPhotonTree->Branch("Jet_Pt", &fVBuffer_Jet_Pt);
    fPhotonTree->Branch("Jet_Eta", &fVBuffer_Jet_Eta);
    fPhotonTree->Branch("Jet_Phi", &fVBuffer_Jet_Phi);
    fPhotonTree->Branch("Jet_NNeutral", &fVBuffer_Jet_NNeutr);
    fPhotonTree->Branch("Jet_NCharged", &fVBuffer_Jet_NCh);
    fPhotonTree->Branch("Jet_NEF", &fVBuffer_Jet_NEF);
    if (fIsMC) {
      fPhotonTree->Branch("Jet_Pt_True", &fVBuffer_TrueJet_Pt);
      fPhotonTree->Branch("Jet_Eta_True", &fVBuffer_TrueJet_Eta);
      fPhotonTree->Branch("Jet_Phi_True", &fVBuffer_TrueJet_Phi);
      fPhotonTree->Branch("Jet_Pt_Parton", &fVBuffer_TrueJet_Parton_Pt);
    }
  }

  if (fSaveMCInformation) {
    fPhotonTree->Branch("MC_gen_ID", &fVBuffer_MCGenID);
    fPhotonTree->Branch("MC_gen_PDG", &fVBuffer_MCGenPDG);
    fPhotonTree->Branch("MC_gen_Px", &fVBuffer_MCGenPx);
    fPhotonTree->Branch("MC_gen_Py", &fVBuffer_MCGenPy);
    fPhotonTree->Branch("MC_gen_Pz", &fVBuffer_MCGenPz);
    fPhotonTree->Branch("MC_gen_MotherID", &fVBuffer_MCGenMotherID);
  }

  fHistoNEvents = new TH1F("NEvents", "NEvents", 15, -0.5, 14.5);
  fHistoNEvents->GetXaxis()->SetBinLabel(1, "Accepted");
  fHistoNEvents->GetXaxis()->SetBinLabel(2, "Centrality");
  fHistoNEvents->GetXaxis()->SetBinLabel(3, "Miss. MC or inc. ev.");
  if (((AliConvEventCuts*)fEventCuts)->IsSpecialTrigger() > 1) {
    TString TriggerNames = "Not Trigger: ";
    TriggerNames = TriggerNames + ((AliConvEventCuts*)fEventCuts)->GetSpecialTriggerName();
    fHistoNEvents->GetXaxis()->SetBinLabel(4, TriggerNames.Data());
  } else {
    fHistoNEvents->GetXaxis()->SetBinLabel(4, "Trigger");
  }
  fHistoNEvents->GetXaxis()->SetBinLabel(5, "Vertex Z");
  fHistoNEvents->GetXaxis()->SetBinLabel(6, "Cont. Vertex");
  fHistoNEvents->GetXaxis()->SetBinLabel(7, "SPD Pile-Up");
  fHistoNEvents->GetXaxis()->SetBinLabel(8, "no SDD");
  fHistoNEvents->GetXaxis()->SetBinLabel(9, "no V0AND");
  fHistoNEvents->GetXaxis()->SetBinLabel(10, "EMCAL/TPC problems");
  fHistoNEvents->GetXaxis()->SetBinLabel(11, "rejectedForJetJetMC");
  fHistoNEvents->GetXaxis()->SetBinLabel(12, "SPD hits vs tracklet");
  fHistoNEvents->GetXaxis()->SetBinLabel(13, "Out-of-Bunch pileup Past-Future");
  fHistoNEvents->GetXaxis()->SetBinLabel(14, "Pileup V0M-TPCout Tracks");
  fHistoNEvents->GetXaxis()->SetBinLabel(15, "No Jet in event");
  fInfoList->Add(fHistoNEvents);

  std::vector<double> vecPt = {0.};
  double epsilon = 1e-12;
  double pt = 0;
  while (true) {
    if (pt < 20 - epsilon)
      pt += 0.1;
    else if (pt < 30 - epsilon)
      pt += 0.5;
    else if (pt < 50 - epsilon)
      pt += 1;
    else if (pt < 200 - epsilon)
      pt += 5;
    else if (pt < 400 - epsilon)
      pt += 10;
    else
      break;
    vecPt.push_back(pt);
  }

  if (fIsMC > 0) {
    fHistoMCGammaPt = new TH2F("MC_Gamma_Pt", "MC_Gamma_Pt", vecPt.size() - 1, vecPt.data(), 60, -1.5, 1.5);
    fInfoList->Add(fHistoMCGammaPt);

    fHistoMCElecPt = new TH2F("MC_Elec_Pt", "MC_Elec_Pt", vecPt.size() - 1, vecPt.data(), 60, -1.5, 1.5);
    fInfoList->Add(fHistoMCElecPt);

    fHistoMCPi0Pt = new TH2F("MC_Pi0_Pt", "MC_Pi0_Pt", vecPt.size() - 1, vecPt.data(), 60, -1.5, 1.5);
    fInfoList->Add(fHistoMCPi0Pt);

    fHistoMCEtaPt = new TH2F("MC_Eta_Pt", "MC_Eta_Pt", vecPt.size() - 1, vecPt.data(), 60, -1.5, 1.5);
    fInfoList->Add(fHistoMCEtaPt);

    fHistoMCOmegaPt = new TH2F("MC_Omega_Pt", "MC_Omega_Pt", vecPt.size() - 1, vecPt.data(), 60, -1.5, 1.5);
    fInfoList->Add(fHistoMCPi0Pt);

    fHistoMCEtaPrimePt = new TH2F("MC_EtaPrime_Pt", "MC_EtaPrime_Pt", vecPt.size() - 1, vecPt.data(), 60, -1.5, 1.5);
    fInfoList->Add(fHistoMCEtaPrimePt);

    fHistoMCK0sPt = new TH2F("MC_K0s_Pt", "MC_K0s_Pt", vecPt.size() - 1, vecPt.data(), 60, -1.5, 1.5);
    fInfoList->Add(fHistoMCK0sPt);

    fHistoMCSecPi0Pt = new TH2F("MC_SecPi0_Pt", "MC_SecPi0_Pt", vecPt.size() - 1, vecPt.data(), 60, -1.5, 1.5);
    fInfoList->Add(fHistoMCSecPi0Pt);

    fHistoMCSecEtaPt = new TH2F("MC_SecEta_Pt", "MC_SecEta_Pt", vecPt.size() - 1, vecPt.data(), 60, -1.5, 1.5);
    fInfoList->Add(fHistoMCSecEtaPt);
  }
  if (fIsMC > 1) {
    TIter iter(fInfoList->MakeIterator());
    while (TObject* obj = iter()) {
      TString className = obj->ClassName();
      if (className.Contains("TH1")) {
        static_cast<TH1*>(obj)->Sumw2();
      } else if (className.Contains("TH2")) {
        static_cast<TH2*>(obj)->Sumw2();
      }
    }
  }

  fAliEventCuts.AddQAplotsToList(fInfoList);

  fV0Reader = (AliV0ReaderV1*)AliAnalysisManager::GetAnalysisManager()->GetTask(fV0ReaderName.Data());

  OpenFile(2);
  PostData(2, fPhotonTree);
}

//_____________________________________________________________________________
Bool_t AliAnalysisTaskConvCaloTree::Notify()
{
  if (fEventCuts->GetPeriodEnum() == AliConvEventCuts::kNoPeriod && ((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetPeriodEnum() != AliConvEventCuts::kNoPeriod) {
    fEventCuts->SetPeriodEnumExplicit(((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetPeriodEnum());
  } else if (fEventCuts->GetPeriodEnum() == AliConvEventCuts::kNoPeriod) {
    fEventCuts->SetPeriodEnum(fV0Reader->GetPeriodName());
  }

  return kTRUE;
}

//________________________________________________________________________
void AliAnalysisTaskConvCaloTree::UserExec(Option_t*)
{

  if (fDoDownscale) {
    if (fDownscaleFac > gRandom->Uniform()) {
      return; // do not accept this event
    }
  }
  AliAnalysisManager* man = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = (AliInputEventHandler*)(man->GetInputEventHandler());
  fPIDResponse = (AliPIDResponse*)inputHandler->GetPIDResponse();

  fInputEvent = InputEvent();
  if (fIsMC > 0)
    fMCEvent = MCEvent();
  Int_t eventQuality = ((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetEventQuality();
  Int_t eventNotAccepted = fEventCuts->IsEventAcceptedByCut(fV0Reader->GetEventCuts(), fInputEvent, fMCEvent, fIsHeavyIon, kFALSE);

  if (eventNotAccepted) {
    fHistoNEvents->Fill(eventNotAccepted);
    return; // Check Centrality, PileUp, SDD and V0AND --> Not Accepted => eventQuality = 1
  }

  // trigger condition
  if (fStoreTriggerCondition) {
    std::bitset<32> triggerbits(0);
    triggerbits[0] = true; // always INT7 triggerd
    auto triggercont = static_cast<PWG::EMCAL::AliEmcalTriggerDecisionContainer*>(fInputEvent->FindListObject("EmcalTriggerDecision"));
    if (triggercont) {
      triggerbits[1] = triggercont->IsEventSelected("EMC7");
      triggerbits[2] = triggercont->IsEventSelected("DMC7");
      triggerbits[3] = triggercont->IsEventSelected("EG2");
      triggerbits[4] = triggercont->IsEventSelected("DG2");
      triggerbits[5] = triggercont->IsEventSelected("EG1");
      triggerbits[6] = triggercont->IsEventSelected("DG1");
      triggerbits[7] = triggercont->IsEventSelected("EJ2");
      triggerbits[8] = triggercont->IsEventSelected("DJ2");
      triggerbits[9] = triggercont->IsEventSelected("EJ1");
      triggerbits[10] = triggercont->IsEventSelected("DJ1");
    }
    // INEL>0
    triggerbits[11] = fEventCuts->IsEventINELgtZERO(fInputEvent);
    if (fIsMC > 0)
      triggerbits[11] = fEventCuts->IsEventTrueINELgtZERO(fInputEvent, fMCEvent);

    fTriggerBits = (unsigned int)(triggerbits.to_ulong());
    // cout << "triggerbits.to_string() " << triggerbits.to_string() << "  fTriggerBits " << fTriggerBits << endl;
  }

  const AliVVertex* primVtxMC = fInputEvent->GetPrimaryVertex();
  fBuffer_Event_Vertex_Z = primVtxMC->GetZ();

  fV0Mult = ((AliConvEventCuts*)fEventCuts)->GetCentrality(fInputEvent);

  if (fIsMC > 0 && fInputEvent->IsA() == AliAODEvent::Class() && !(fV0Reader->AreAODsRelabeled())) {
    RelabelAODPhotonCandidates(kTRUE); // In case of AODMC relabeling MC
    fV0Reader->RelabelAODs(kTRUE);
  }

  fWeightJetJetMC = 1;
  if (fIsMC > 1) {
    Float_t pthard = -1;
    Bool_t isMCJet = ((AliConvEventCuts*)fEventCuts)->IsJetJetMCEventAccepted(fMCEvent, fWeightJetJetMC, pthard, fInputEvent);
    if (!isMCJet) {
      return;
    }
    fBuffer_EventWeight = fWeightJetJetMC;
  }

  if (fIsMC) {
    // only works for pythia right now
    AliGenEventHeader* eventHeader = fMCEvent->GenEventHeader();
    fMCQ2 = dynamic_cast<AliGenPythiaEventHeader*>(eventHeader)->GetPtHard();
  }

  if (eventNotAccepted != 0) {
    fHistoNEvents->Fill(eventNotAccepted, fWeightJetJetMC); // Check Centrality, PileUp, SDD and V0AND --> Not Accepted => eventQuality = 1

    if (fIsMC > 0 && (eventQuality == 0 || eventQuality == 3 || eventQuality == 5)) {
      if (eventNotAccepted == 3) { // wrong trigger selected. However, we still want to count the MC particles fot these events! If MC particles should be rejected in addition, use IsMCTriggerSelected function
        if (fSaveMCInformation) {
          ProcessAODMCParticles();
        }
      }
    }
    fPhotonTree->Fill();
    return;
  }
  if (eventQuality != 0) { // Event Not Accepted
    // cout << "event rejected due to: " <<eventQuality << endl;
    fHistoNEvents->Fill(eventQuality, fWeightJetJetMC);

    if (fIsMC > 0) {
      if (eventQuality == 3 || eventQuality == 5) { // 3 = trigger, 5 = contr. to vertex
        ProcessAODMCParticles();
      }
    }
    fPhotonTree->Fill();
    return;
  }

  if (fSaveJets) {
    InitJets();
    if (fSaveJets == 2) {
      if (fVBuffer_Jet_Pt.size() == 0) {
        fHistoNEvents->Fill(14);
        return;
      }
    }
  }

  if (fIsMC > 0 && fSaveMCInformation) {
    ProcessAODMCParticles();
  }

  if (!fAliEventCuts.AcceptEvent(fInputEvent)) {
    return;
  }

  fHistoNEvents->Fill(0);

  if (fSaveClusters) {
    ProcessClustersAOD();
  }

  if (fSaveConversions) {
    ProcessConversionsAOD();
  }

  if (fSaveTracks > 0) {
    ProcessTracksAOD();
  }

  fPhotonTree->Fill();

  ResetBufferVectors();

  if (fIsMC > 0 && fInputEvent->IsA() == AliAODEvent::Class() && !(fV0Reader->AreAODsRelabeled())) {
    RelabelAODPhotonCandidates(kFALSE); // Back to ESDMC Label
    fV0Reader->RelabelAODs(kFALSE);
  }

  PostData(1, fOutputList);
}

///________________________________________________________________________
void AliAnalysisTaskConvCaloTree::ProcessClustersAOD()
{

  Int_t nclus = 0;
  if (!fCorrTaskSetting.CompareTo("")) {
    nclus = fInputEvent->GetNumberOfCaloClusters();
  } else {
    if (!farrClustersProcess)
      farrClustersProcess = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(Form("%sClustersBranch", fCorrTaskSetting.Data())));
    if (!farrClustersProcess)
      AliFatal(Form("%sClustersBranch was not found in AliAnalysisTaskGammaCalo! Check the correction framework settings!", fCorrTaskSetting.Data()));
    nclus = farrClustersProcess->GetEntries();
  }

  if (nclus == 0)
    return;

  AliVCluster* clus = NULL;
  for (Long_t i = 0; i < nclus; i++) {

    if (fInputEvent->IsA() == AliAODEvent::Class()) {
      if (farrClustersProcess)
        clus = new AliAODCaloCluster(*(AliAODCaloCluster*)farrClustersProcess->At(i));
      else
        clus = new AliAODCaloCluster(*(AliAODCaloCluster*)fInputEvent->GetCaloCluster(i));
    } else if (fInputEvent->IsA() == AliESDEvent::Class()) {
      AliFatal("ESD running is not supported");
    }
    if (!clus)
      continue;

    Bool_t clusIsAcc = kFALSE;
    if (fClusterCutsEMC && clus->IsEMCAL()) {
      if (((AliCaloPhotonCuts*)fClusterCutsEMC)->ClusterIsSelected(clus, fInputEvent, fMCEvent, fIsMC, fWeightJetJetMC, i)) {
        clusIsAcc = kTRUE;
        fClusterCuts = fClusterCutsEMC;
      }
    }
    if (fClusterCutsPHOS && clus->IsPHOS() && !clusIsAcc) {
      if (((AliCaloPhotonCuts*)fClusterCutsPHOS)->ClusterIsSelected(clus, fInputEvent, fMCEvent, fIsMC, fWeightJetJetMC, i)) {
        clusIsAcc = kTRUE;
        fClusterCuts = fClusterCutsPHOS;
      }
    }
    if (!clusIsAcc) {
      delete clus;
      continue;
    }

    Double_t vertex[3] = {0};
    InputEvent()->GetPrimaryVertex()->GetXYZ(vertex);
    TLorentzVector clusterVector;
    clus->GetMomentum(clusterVector, vertex);
    TLorentzVector tmpvec;
    tmpvec.SetPxPyPzE(clusterVector.Px(), clusterVector.Py(), clusterVector.Pz(), clusterVector.E());
    AliAODConversionPhoton* PhotonCandidate = new AliAODConversionPhoton(&tmpvec);

    Float_t etaCluster = (Float_t)PhotonCandidate->Eta();
    Float_t phiCluster = (Float_t)PhotonCandidate->Phi();
    if (phiCluster < 0)
      phiCluster += 2 * TMath::Pi();

    fVBuffer_Cluster_E.push_back(static_cast<float>(clus->E()));
    fVBuffer_Cluster_Eta.push_back(static_cast<float>(etaCluster));
    fVBuffer_Cluster_Phi.push_back(static_cast<float>(phiCluster));
    fVBuffer_Cluster_NCells.push_back(static_cast<short>((Short_t)clus->GetNCells()));
    fVBuffer_Cluster_M02.push_back((Short_t)(clus->GetM02() * 1000));
    fVBuffer_Cluster_Time.push_back((Short_t)(clus->GetTOF() * 1e9)); // time in ns
    // calculate Ecross
    if (fClusterCutsEMC) {
      AliVCaloCells* cells = fInputEvent->GetEMCALCells();
      Int_t largestCellID = ((AliCaloPhotonCuts*)fClusterCutsEMC)->FindLargestCellInCluster(clus, fInputEvent);
      Float_t ecell1 = cells->GetCellAmplitude(largestCellID);
      Float_t eCross = ((AliCaloPhotonCuts*)fClusterCutsEMC)->GetECross(largestCellID, cells);
      fVBuffer_Cluster_Fcross.push_back((Short_t)(1000 * (1 - eCross / ecell1)));
      fVBuffer_Cluster_Exotic.push_back(clus->GetIsExotic());

      Int_t nCellCluster = clus->GetNCells();
      for (Int_t iCell = 0; iCell < nCellCluster; iCell++) {
        Int_t cellID = clus->GetCellAbsId(iCell);
        Double_t cellAmp = cells->GetCellAmplitude(cellID);
        fVBuffer_Cluster_CellEnergies.push_back(cellAmp);
        fVBuffer_Cluster_CellIDs.push_back(cellID);
        if (fIsMC == 0) {
          Double_t cellTime = cells->GetCellTime(cellID);
          fVBuffer_Cluster_CellTimes.push_back(cellTime * 1e9);
        }
      }
    }
    if (fUseClusterIsolation && fCaloIsolation) {
      Float_t ClusterPt = PhotonCandidate->Pt();
      Float_t pTCone = 0.1 * ClusterPt; // momPercantage = 0.1
      fVBuffer_Cluster_Isolated.push_back(fCaloIsolation->GetIsolation(i, 0.4, pTCone));
    }

    if (fIsMC > 0) { // MC info

      const AliVVertex* primVtxMC = fMCEvent->GetPrimaryVertex();
      const double mcProdVtxX = primVtxMC->GetX();
      const double mcProdVtxY = primVtxMC->GetY();
      const double mcProdVtxZ = primVtxMC->GetZ();

      if (!PhotonCandidate) {
        fVTrueClusterPi0DaughterIndex.push_back(0);
        fVTrueClusterEtaDaughterIndex.push_back(0);
        fVTrueClusterSecondaryIndex.push_back(0);
        fVTrueClusterMCId.push_back(0);
        fVTrueClusterMCId2.push_back(0);
        fVTrueClusterMCTrueEnergy.push_back(0);
        fVTrueClusterMCTrueEnergy2.push_back(0);
        fVTrueClusterMCIsMerged.push_back(kFALSE);
        fVTrueClusterConvRadius.push_back(0);
        delete clus;
        continue;
      }
      PhotonCandidate->SetIsCaloPhoton(((AliCaloPhotonCuts*)fClusterCuts)->GetClusterType());
      PhotonCandidate->SetCaloClusterRef(i);
      PhotonCandidate->SetLeadingCellID(((AliCaloPhotonCuts*)fClusterCuts)->FindLargestCellInCluster(clus, fInputEvent));
      Int_t* mclabelsCluster = clus->GetLabels();
      PhotonCandidate->SetNCaloPhotonMCLabels(clus->GetNLabels());

      if (clus->GetNLabels() > 0) {
        for (Int_t k = 0; k < (Int_t)clus->GetNLabels(); k++) {
          PhotonCandidate->SetCaloPhotonMCLabel(k, mclabelsCluster[k]);
        }
      }

      if (!fAODMCTrackArray)
        fAODMCTrackArray = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
      if (!fAODMCTrackArray)
        return;

      PhotonCandidate->SetCaloPhotonMCFlagsAOD(fAODMCTrackArray, kTRUE);

      // get true MC id label
      const int NLabels = 2;
      std::vector<Short_t> arrMCID = {0, 0};
      std::vector<Float_t> arrMCEnergy = {0, 0};
      for (Int_t i = 0; i < NLabels; ++i) {
        Int_t gammaMCLabel = PhotonCandidate->GetCaloPhotonMCLabel(i); // get most probable MC label
        if (gammaMCLabel != -1) {
          AliAODMCParticle* gammaMC = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(gammaMCLabel));
          arrMCID[i] = (Short_t)gammaMC->GetPdgCode();
          arrMCEnergy[i] = (Float_t)gammaMC->E();
        }
      }
      fVTrueClusterMCId.push_back(static_cast<short>(arrMCID[0]));
      fVTrueClusterMCId2.push_back(static_cast<short>(arrMCID[1]));
      fVTrueClusterMCTrueEnergy.push_back(static_cast<short>(arrMCEnergy[0]));
      fVTrueClusterMCTrueEnergy2.push_back(static_cast<short>(arrMCEnergy[1]));
      Int_t gammaMCLabel = PhotonCandidate->GetCaloPhotonMCLabel(0); // get most probable MC label
      Int_t gammaMotherLabel = -1;
      Double_t convRadiusConvClus = 0;

      AliAODMCParticle* gammaMC = nullptr;
      if (gammaMCLabel != -1) {
        gammaMC = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(gammaMCLabel));
        if (PhotonCandidate->IsLargestComponentPhoton() || PhotonCandidate->IsLargestComponentElectron()) { // largest component is electron
          if (PhotonCandidate->IsLargestComponentPhoton()) {                                                // for photons its the direct mother
            gammaMotherLabel = gammaMC->GetMother();
          } else if (PhotonCandidate->IsLargestComponentElectron()) { // for electrons its either the direct mother or for conversions the grandmother
            if (PhotonCandidate->IsConversion()) {
              AliAODMCParticle* gammaGrandMotherMC = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(gammaMC->GetMother()));
              gammaMotherLabel = gammaGrandMotherMC->GetMother();
              // get the production vertex of the conversion electron
              convRadiusConvClus = sqrt(gammaMC->Xv() * gammaMC->Xv() + gammaMC->Yv() * gammaMC->Yv());
            } else
              gammaMotherLabel = gammaMC->GetMother();
          }
        }
        if (PhotonCandidate->IsMerged()) {
          fVTrueClusterMCIsMerged.push_back(true);
        } else {
          fVTrueClusterMCIsMerged.push_back(false);
        }

        // check that we really captured all pi0s and etas
        AliAODMCParticle* tmpMother = nullptr;
        int tmpMotherLabel = gammaMCLabel;
        for (int i = 0; i < 10; ++i) {
          tmpMother = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(tmpMotherLabel));
          int tmpPDGCode = tmpMother->GetPdgCode();
          if (tmpPDGCode == 111 || tmpPDGCode == 221) {
            if (tmpMotherLabel == gammaMotherLabel) {
              break; // in this case, all is fine :)
            } else { // in this case the chain was somehow longer... Setting the label to the new one
              gammaMotherLabel = tmpMotherLabel;
              break;
            }
          }
        }
      } else {
        fVTrueClusterMCIsMerged.push_back(false);
      }
      // conv radius of conv cluster
      fVTrueClusterConvRadius.push_back(static_cast<Short_t>(convRadiusConvClus));

      if (gammaMotherLabel > -1) {

        // get grandmother label to check for secondaries
        int grandmotherLabel = ((AliAODMCParticle*)fAODMCTrackArray->At(gammaMotherLabel))->GetMother();
        int labelSecondary = 0;
        if (grandmotherLabel > 0) {
          AliAODMCParticle* gammaGrandMotherMC = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(grandmotherLabel));
          if (gammaGrandMotherMC) {
            Bool_t isPrimary = ((AliConvEventCuts*)fEventCuts)->IsConversionPrimaryAOD(fInputEvent, gammaMC, mcProdVtxX, mcProdVtxY, mcProdVtxZ);
            if (!isPrimary) {
              labelSecondary = grandmotherLabel;
            }
          }
        }
        fVTrueClusterSecondaryIndex.push_back(static_cast<unsigned short>(labelSecondary));

        if (((AliAODMCParticle*)fAODMCTrackArray->At(gammaMotherLabel))->GetPdgCode() == 111) {
          fVTrueClusterPi0DaughterIndex.push_back(static_cast<unsigned short>(gammaMotherLabel));
          fVTrueClusterEtaDaughterIndex.push_back(0);
        } else if (((AliAODMCParticle*)fAODMCTrackArray->At(gammaMotherLabel))->GetPdgCode() == 221) {
          fVTrueClusterPi0DaughterIndex.push_back(0);
          fVTrueClusterEtaDaughterIndex.push_back(static_cast<unsigned short>(gammaMotherLabel));
        } else {
          fVTrueClusterPi0DaughterIndex.push_back(0);
          fVTrueClusterEtaDaughterIndex.push_back(0);
        }
      } else {
        fVTrueClusterPi0DaughterIndex.push_back(0);
        fVTrueClusterEtaDaughterIndex.push_back(0);
      }
    }
  }
}

///________________________________________________________________________
void AliAnalysisTaskConvCaloTree::ProcessConversionsAOD()
{

  fReaderGammas = fV0Reader->GetReconstructedGammas(); // Gammas from default Cut
  for (Int_t i = 0; i < fReaderGammas->GetEntriesFast(); i++) {
    AliAODConversionPhoton* PhotonCandidate = (AliAODConversionPhoton*)fReaderGammas->At(i);
    if (!PhotonCandidate)
      continue;

    // apply cuts
    if (!fConversionCuts)
      AliFatal("Conversion Cuts Not Initialized");
    if (!((AliConversionPhotonCuts*)fConversionCuts)->PhotonIsSelected(PhotonCandidate, fInputEvent))
      continue;

    fVBuffer_Conv_px.push_back(PhotonCandidate->GetPx());
    fVBuffer_Conv_py.push_back(PhotonCandidate->GetPy());
    fVBuffer_Conv_pz.push_back(PhotonCandidate->GetPz());

    if (fSaveConversions > 1) {
      fVBuffer_Conv_R.push_back(static_cast<UShort_t>(PhotonCandidate->GetConversionRadius() * 10));                                                                       // conv radius *10
      fVBuffer_Conv_PsiPair.push_back(static_cast<Short_t>(PhotonCandidate->GetPsiPair() * 1000));                                                                         // psi pair *1000
      fVBuffer_Conv_CosPAngle.push_back(static_cast<Short_t>(((AliConversionPhotonCuts*)fConversionCuts)->GetCosineOfPointingAngle(PhotonCandidate, fInputEvent) * 1000)); // pointing angle *1000
      fVBuffer_Conv_Chi2.push_back(static_cast<UShort_t>(PhotonCandidate->GetChi2perNDF() * 10));                                                                          // pointing angle *10
      fVBuffer_Conv_Quality.push_back(static_cast<Short_t>(PhotonCandidate->GetPhotonQuality()));
    }

    // get track eta and phi on Calo surface
    for (Int_t iElec = 0; iElec < 2; iElec++) {
      Int_t tracklabel = PhotonCandidate->GetLabel(iElec);
      AliVTrack* inTrack = 0x0;
      if (fInputEvent->IsA() == AliESDEvent::Class()) {
        AliFatal("ESD running is not supported");
      } else { /// AOD
        if (((AliV0ReaderV1*)AliAnalysisManager::GetAnalysisManager()->GetTask(fV0ReaderName.Data()))->AreAODsRelabeled()) {
          inTrack = dynamic_cast<AliVTrack*>(fInputEvent->GetTrack(tracklabel));
        } else {
          for (Int_t ii = 0; ii < fInputEvent->GetNumberOfTracks(); ii++) {
            inTrack = dynamic_cast<AliVTrack*>(fInputEvent->GetTrack(ii));
            if (inTrack) {
              if (inTrack->GetID() == tracklabel) {
                break;
              }
            }
          }
        }

        if (fSaveConversions > 1) {
          // Short_t Charge    = inTrack->Charge();
          Double_t electronNSigmaTPC = fPIDResponse->NumberOfSigmasTPC(inTrack, AliPID::kElectron);
          // Double_t P=inTrack->P();;
          // Double_t Eta=inTrack->Eta();
          // Double_t nSigdEdxCorr = ((AliConversionPhotonCuts*)fConversionCuts)->GetCorrectedElectronTPCResponse(Charge,electronNSigmaTPC,P,Eta,inTrack->GetTPCNcls(),PhotonCandidate->GetConversionRadius());
          if (iElec == 0) {
            fVBuffer_Conv_NTPCClusElec1.push_back(static_cast<UShort_t>(inTrack->GetTPCNcls())); // NTPC clus
            fVBuffer_Conv_dEdxElec1.push_back(static_cast<Short_t>(electronNSigmaTPC * 100));    // dedx
            fVBuffer_Conv_PElec1.push_back(static_cast<Short_t>(inTrack->P() * 100));            // momentum
          } else {
            fVBuffer_Conv_NTPCClusElec2.push_back(static_cast<UShort_t>(inTrack->GetTPCNcls())); // NTPC clus
            fVBuffer_Conv_dEdxElec2.push_back(static_cast<Short_t>(electronNSigmaTPC * 100));    // dedx
            fVBuffer_Conv_PElec2.push_back(static_cast<Short_t>(inTrack->P() * 100));            // momentum
          }
        }

        // track extrapolation
        Double_t xyz[3] = {0}, pxpypz[3] = {0}, cv[21] = {0};
        inTrack->GetPxPyPz(pxpypz);
        inTrack->GetXYZ(xyz);
        inTrack->GetCovarianceXYZPxPyPz(cv);

        AliExternalTrackParam emcParam(xyz, pxpypz, cv, inTrack->Charge());
        Float_t eta, phi, pt;

        // propagate tracks to emc surfaces
        if (!AliEMCALRecoUtils::ExtrapolateTrackToEMCalSurface(&emcParam, 440., 0.139, 20., eta, phi, pt)) {
          if (iElec == 0) {
            fVBuffer_Elec1etaCalo.push_back(eta);
            fVBuffer_Elec1phiCalo.push_back(phi);
          } else {
            fVBuffer_Elec2etaCalo.push_back(eta);
            fVBuffer_Elec2phiCalo.push_back(phi);
          }
          continue;
        }
        if (iElec == 0) {
          fVBuffer_Elec1etaCalo.push_back(eta);
          fVBuffer_Elec1phiCalo.push_back(phi);
        } else {
          fVBuffer_Elec2etaCalo.push_back(eta);
          fVBuffer_Elec2phiCalo.push_back(phi);
        }
      }
    }

    if (fIsMC > 0) {
      if (!fAODMCTrackArray)
        fAODMCTrackArray = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
      if (fAODMCTrackArray != NULL && PhotonCandidate != NULL) {

        AliAODMCParticle* posDaughter = (AliAODMCParticle*)fAODMCTrackArray->At(PhotonCandidate->GetMCLabelPositive());

        AliAODMCParticle* Photon = (AliAODMCParticle*)fAODMCTrackArray->At(posDaughter->GetMother());

        fVTrueConvMCLabel.push_back(static_cast<UShort_t>(Photon->GetLabel()));

        if (Photon->GetPdgCode() != 22) {
          fVTrueConvEtaDaughterIndex.push_back(0);
          fVTrueConvPi0DaughterIndex.push_back(0);
          fVTrueConvSecDaughterIndex.push_back(0);
          fVTrueConvMCTruePx.push_back(0.);
          fVTrueConvMCTruePy.push_back(0.);
          fVTrueConvMCTruePz.push_back(0.);
          return; // Mother is no Photon
        }
        // fill true px, py and pz
        fVTrueConvMCTruePx.push_back(Photon->Px());
        fVTrueConvMCTruePy.push_back(Photon->Py());
        fVTrueConvMCTruePz.push_back(Photon->Pz());

        Int_t gammaMotherLabel = Photon->GetMother(); // get the mother stack ID
        if (((AliAODMCParticle*)fAODMCTrackArray->At(Photon->GetMother()))->GetPdgCode() == 111) {
          fVTrueConvPi0DaughterIndex.push_back(static_cast<unsigned short>(gammaMotherLabel));
          fVTrueConvEtaDaughterIndex.push_back(0);
        } else if (((AliAODMCParticle*)fAODMCTrackArray->At(Photon->GetMother()))->GetPdgCode() == 221) {
          fVTrueConvEtaDaughterIndex.push_back(static_cast<unsigned short>(gammaMotherLabel));
          fVTrueConvPi0DaughterIndex.push_back(0);
        } else {
          fVTrueConvEtaDaughterIndex.push_back(0);
          fVTrueConvPi0DaughterIndex.push_back(0);
        }
      }
    }
  }
}

//________________________________________________________________________
void AliAnalysisTaskConvCaloTree::ProcessTracksAOD()
{

  AliAODEvent* aodev = dynamic_cast<AliAODEvent*>(fInputEvent);
  if (!aodev) {
    AliError("Task needs AODs , returning");
    return;
  }

  if (!fAODMCTrackArray)
        fAODMCTrackArray = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));

  for (Int_t itr = 0; itr < fInputEvent->GetNumberOfTracks(); itr++) {

    AliVTrack* inTrack = dynamic_cast<AliVTrack*>(aodev->GetTrack(itr));
    if (!inTrack)
      continue;
    if (inTrack->Pt() < fMinTrackPt)
      continue;
    AliAODTrack* aodt = dynamic_cast<AliAODTrack*>(inTrack);
    // check track cuts
    if (!aodt->IsHybridGlobalConstrainedGlobal())
      continue;

    if (fSaveTracks == 2) {
      fVBuffer_Track_px.push_back(static_cast<Short_t>(aodt->Px() * 100));
      fVBuffer_Track_py.push_back(static_cast<Short_t>(aodt->Py() * 100));
      fVBuffer_Track_pz.push_back(static_cast<Short_t>(aodt->Pz() * 100));
      fVBuffer_Track_dedx.push_back(static_cast<unsigned short>(aodt->GetTPCsignal() * 100));

      if ((aodt->GetStatus() & AliVTrack::kTOFout) && (aodt->GetStatus() & AliVTrack::kTIME)) {
        double t0 = fPIDResponse->GetTOFResponse().GetStartTime(aodt->P());
        double times[AliPID::kSPECIESC];
        aodt->GetIntegratedTimes(times, AliPID::kSPECIESC);
        double TOFsignal = aodt->GetTOFsignal();
        double dT = TOFsignal - t0 - times[0];

        fVBuffer_Track_TOFSignal.push_back(static_cast<unsigned short>(dT * 1000));
      } else {
        fVBuffer_Track_TOFSignal.push_back(0);
      }

      float b[2]; // dcaxy, dcaz
      float bCov[3];
      aodt->GetImpactParameters(b, bCov);

      fVBuffer_Track_DCA.push_back(static_cast<unsigned short>(b[0] * 100));

      // fraction of findable clusters
      const double maxNClus = aodt->GetTPCNclsF();
      double fracClsFound = maxNClus > 0 ? ((double)aodt->GetTPCncls())/maxNClus : 0;

      double chi2Track = aodt->GetTPCchi2perCluster();
      unsigned short chi2TrackShort = static_cast<unsigned short>(10 + (chi2Track * 10)); // times 10 to have it in granularity of 0.1
      if(chi2Track > 64) chi2TrackShort = 64000; // maximum size of unsigned short
      else chi2TrackShort*=100; // 

      fVBuffer_Track_FracNClus.push_back(static_cast<unsigned short>(fracClsFound * 100) + chi2TrackShort);
      
      if (fIsMC) {
        fVBuffer_Track_StackID.push_back(static_cast<unsigned short>(std::abs(aodt->GetLabel())));
        AliAODMCParticle* MCTrack = (AliAODMCParticle*)fAODMCTrackArray->At(std::abs(aodt->GetLabel()));
        fVBuffer_Track_PDG.push_back(static_cast<short>(MCTrack->GetPdgCode()));
      }
    }
    // This is not the proper extrapolation! Although for a rough matching of tracks and clusters it was found to save a lot of CPU time
    Float_t trackEta = aodt->GetTrackEtaOnEMCal();
    Float_t trackPhi = aodt->GetTrackPhiOnEMCal();
    if (trackEta == -999 || trackPhi == -999) {
      // maximum values for Short_t (--> not matched to Calo)
      // Store this information if all tracks are stored
      if (fSaveTracks == 2) {
        fVBuffer_Track_Calo_eta.push_back(static_cast<Short_t>(32767));
        fVBuffer_Track_Calo_phi.push_back(static_cast<UShort_t>(65535));
      }
      continue;
    }
    if (fSaveTracks == 1) {
      fVBuffer_Track_P.push_back(static_cast<Short_t>(aodt->P() * 100));
    }
    fVBuffer_Track_Calo_eta.push_back(static_cast<Short_t>(trackEta * 10000));
    fVBuffer_Track_Calo_phi.push_back(static_cast<UShort_t>(trackPhi * 10000));
  }
}

//________________________________________________________________________
bool AliAnalysisTaskConvCaloTree::IsMCParticleRelevant(AliAODMCParticle* particle)
{
  if (particle->Pt() < 0.1)
    return false;
  if (std::abs(particle->Y()) > 1.5)
    return false;
  const int pdgCode = particle->GetPdgCode();
  if (pdgCode == 111 ||          // pi0
      pdgCode == 221 ||          // eta
      pdgCode == 223 ||          // omega
      pdgCode == 331 ||          // eta prime
      pdgCode == 310 ||          // K0s
      pdgCode == 130 ||          // K0l
      pdgCode == 3122 ||         // Lambda
      std::abs(pdgCode) == 11 || // electron
      std::abs(pdgCode) == 22    // photon
  ) {
    return true;
  }
  if (particle->GetMother() <= 2) { // has to originate from protons or system
    return true;
  }
  if (particle->GetDaughterFirst() <= 0) { // has no daughter -> stable particle
    return true;
  }
  // check if particle decayed after a certain time and could still be a valid track
  auto labelDaughter = particle->GetDaughterFirst();
  AliAODMCParticle* daughter = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(labelDaughter));
  float rVtxDaughter = sqrt(daughter->Xv()*daughter->Xv() + daughter->Yv()*daughter->Yv());

  auto labelGrandDaughter = daughter->GetDaughterFirst();
  AliAODMCParticle* granddaughter = nullptr;
  float rVtxGrandDaughter = 0;
  if(labelGrandDaughter > 1){
    granddaughter = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(labelGrandDaughter));
    rVtxGrandDaughter = sqrt(granddaughter->Xv()*granddaughter->Xv() + granddaughter->Yv()*granddaughter->Yv());
  }

  if(rVtxDaughter > 0.05){ 
    return true;
  }
  float rVtxPart = sqrt(particle->Xv()*particle->Xv() + particle->Yv()*particle->Yv());
  if(rVtxPart > 0.03 && rVtxDaughter > 0.05){ // 3cm
    return true;
  }

  return false;
}

//________________________________________________________________________
void AliAnalysisTaskConvCaloTree::ProcessAODMCParticles()
{
  if (fIsMC == 0)
    return;
  const AliVVertex* primVtxMC = fMCEvent->GetPrimaryVertex();
  Double_t mcProdVtxX = primVtxMC->GetX();
  Double_t mcProdVtxY = primVtxMC->GetY();
  Double_t mcProdVtxZ = primVtxMC->GetZ();

  if (!fAODMCTrackArray)
    fAODMCTrackArray = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
  if (fAODMCTrackArray == NULL)
    return;
  // Loop over all primary MC particle
  for (Long_t i = 0; i < fAODMCTrackArray->GetEntriesFast(); i++) {
    Double_t tempParticleWeight = fWeightJetJetMC;
    AliAODMCParticle* particle = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(i));

    if (!particle)
      continue;

    if ((fSaveMCInformation == 2 && IsMCParticleRelevant(particle)) || fSaveMCInformation == 3) {
      // stack ID
      fVBuffer_MCGenID.push_back(static_cast<unsigned short>(i));
      // pdg code
      fVBuffer_MCGenPDG.push_back(static_cast<short>(particle->GetPdgCode()));
      // px
      fVBuffer_MCGenPx.push_back(static_cast<short>(particle->Px() * 100));
      // py
      fVBuffer_MCGenPy.push_back(static_cast<short>(particle->Py() * 100));
      // pz
      fVBuffer_MCGenPz.push_back(static_cast<short>(particle->Pz() * 100));
      // mother
      fVBuffer_MCGenMotherID.push_back(static_cast<unsigned short>(particle->GetMother()));
    }

    Bool_t isPrimary = fEventCuts->IsConversionPrimaryAOD(fInputEvent, particle, mcProdVtxX, mcProdVtxY, mcProdVtxZ);
    if (isPrimary) {
      // Int_t isMCFromMBHeader = -1;

      // if(fEventCuts->GetSignalRejection() != 0){
      //   isMCFromMBHeader = fEventCuts->IsParticleFromBGEvent(i, fMCEvent, fInputEvent);
      //   if(isMCFromMBHeader == 0 && fEventCuts->GetSignalRejection() != 3) continue;
      //   // Set the jetjet weight to 1 in case the particle orignated from the minimum bias header
      //   if(isMCFromMBHeader == 2 && fEventCuts->GetSignalRejection() == 4) tempParticleWeight = 1;
      // }

      if (fMesonCuts->MesonIsSelectedAODMC(particle, fAODMCTrackArray, fEventCuts->GetEtaShift())) {
        Float_t weighted = 1;

        // if(fEventCuts->IsParticleFromBGEvent(i, fMCEvent, fInputEvent)){
        //   if (particle->Pt()>0.005){
        //     weighted = fEventCuts->GetWeightForMeson(i, 0x0, fInputEvent);
        //   }
        // }

        if (particle->GetPdgCode() == 22) {
          fHistoMCGammaPt->Fill(particle->Pt(), particle->Y(), weighted * tempParticleWeight); // All MC Gamma
        } else if (std::abs(particle->GetPdgCode()) == 11) {
          fHistoMCElecPt->Fill(particle->Pt(), particle->Y(), weighted * tempParticleWeight); // All MC electrons
        } else if (particle->GetPdgCode() == 111) {
          fHistoMCPi0Pt->Fill(particle->Pt(), particle->Y(), weighted * tempParticleWeight); // All MC Pi0
        } else if (particle->GetPdgCode() == 221) {
          fHistoMCEtaPt->Fill(particle->Pt(), particle->Y(), weighted * tempParticleWeight); // All MC Eta
        } else if (particle->GetPdgCode() == 223) {
          fHistoMCOmegaPt->Fill(particle->Pt(), particle->Y(), weighted * tempParticleWeight); // All MC Omega
        } else if (particle->GetPdgCode() == 331) {
          fHistoMCEtaPrimePt->Fill(particle->Pt(), particle->Y(), weighted * tempParticleWeight); // All MC Eta Prime
        } else if (particle->GetPdgCode() == 310) {
          fHistoMCK0sPt->Fill(particle->Pt(), particle->Y(), weighted * tempParticleWeight); // All MC K0s
        }
      }
    } else {
      if (fMesonCuts->MesonIsSelectedAODMC(particle, fAODMCTrackArray, fEventCuts->GetEtaShift())) {
        Float_t weighted = 1;

        // if(fEventCuts->IsParticleFromBGEvent(i, fMCEvent, fInputEvent)){
        //   if (particle->Pt()>0.005){
        //     weighted = fEventCuts->GetWeightForMeson(i, 0x0, fInputEvent);
        //   }
        // }
        if (particle->GetPdgCode() == 111) {
          fHistoMCSecPi0Pt->Fill(particle->Pt(), particle->Y(), weighted * tempParticleWeight); // All MC sec Pi0
        } else if (particle->GetPdgCode() == 221) {
          fHistoMCSecEtaPt->Fill(particle->Pt(), particle->Y(), weighted * tempParticleWeight); // All MC sec. Eta
        }
      }
    }
  }
}

//________________________________________________________________________
void AliAnalysisTaskConvCaloTree::InitJets()
{

  if (fIsMC) {
    if (!fAODMCTrackArray)
      fAODMCTrackArray = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
    fConvJetReader->FindPartonsJet(fAODMCTrackArray);
  }
  for (unsigned int i = 0; i < (fConvJetReader->GetVectorJetPt()).size(); ++i) {
    fVBuffer_Jet_Pt.push_back((fConvJetReader->GetVectorJetPt())[i]);
    fVBuffer_Jet_Eta.push_back((static_cast<short>(fConvJetReader->GetVectorJetEta()[i] * 1000)));
    fVBuffer_Jet_Phi.push_back((static_cast<short>(fConvJetReader->GetVectorJetPhi()[i] * 1000)));
    fVBuffer_Jet_NNeutr.push_back(static_cast<unsigned short>(fConvJetReader->GetVectorJetNclus()[i]));
    fVBuffer_Jet_NCh.push_back(static_cast<unsigned short>(fConvJetReader->GetVectorJetNtracks()[i]));
    fVBuffer_Jet_NEF.push_back(static_cast<unsigned short>(fConvJetReader->GetVectorJetNEF()[i]));
  }
  if (fIsMC) {
    for (unsigned int i = 0; i < (fConvJetReader->GetTrueVectorJetPt()).size(); ++i) {
      fVBuffer_TrueJet_Pt.push_back((fConvJetReader->GetTrueVectorJetPt())[i]);
      fVBuffer_TrueJet_Eta.push_back((static_cast<short>(fConvJetReader->GetTrueVectorJetEta()[i] * 1000)));
      fVBuffer_TrueJet_Phi.push_back((static_cast<short>(fConvJetReader->GetTrueVectorJetPhi()[i] * 1000)));
      fVBuffer_TrueJet_Parton_Pt.push_back((fConvJetReader->GetTrueVectorJetPartonPt())[i]);
    }
  }
}

//________________________________________________________________________
void AliAnalysisTaskConvCaloTree::Terminate(Option_t*)
{
}

//________________________________________________________________________
void AliAnalysisTaskConvCaloTree::RelabelAODPhotonCandidates(Bool_t mode)
{

  // Relabeling For AOD Event
  // ESDiD -> AODiD
  // MCLabel -> AODMCLabel
  if (mode) {
    fMCEventPos.clear();
    fMCEventNeg.clear();
    fESDArrayPos.clear();
    fESDArrayNeg.clear();
    fMCEventPos.resize(fReaderGammas->GetEntries());
    fMCEventNeg.resize(fReaderGammas->GetEntries());
    fESDArrayPos.resize(fReaderGammas->GetEntries());
    fESDArrayNeg.resize(fReaderGammas->GetEntries());
  }

  for (Int_t iGamma = 0; iGamma < fReaderGammas->GetEntries(); iGamma++) {
    AliAODConversionPhoton* PhotonCandidate = (AliAODConversionPhoton*)fReaderGammas->At(iGamma);
    if (!PhotonCandidate)
      continue;
    if (!mode) { // Back to ESD Labels
      PhotonCandidate->SetMCLabelPositive(fMCEventPos[iGamma]);
      PhotonCandidate->SetMCLabelNegative(fMCEventNeg[iGamma]);
      PhotonCandidate->SetLabelPositive(fESDArrayPos[iGamma]);
      PhotonCandidate->SetLabelNegative(fESDArrayNeg[iGamma]);
      continue;
    }
    fMCEventPos[iGamma] = PhotonCandidate->GetMCLabelPositive();
    fMCEventNeg[iGamma] = PhotonCandidate->GetMCLabelNegative();
    fESDArrayPos[iGamma] = PhotonCandidate->GetTrackLabelPositive();
    fESDArrayNeg[iGamma] = PhotonCandidate->GetTrackLabelNegative();

    Bool_t AODLabelPos = kFALSE;
    Bool_t AODLabelNeg = kFALSE;

    for (Int_t i = 0; i < fInputEvent->GetNumberOfTracks(); i++) {
      AliAODTrack* tempDaughter = static_cast<AliAODTrack*>(fInputEvent->GetTrack(i));
      if (!AODLabelPos) {
        if (tempDaughter->GetID() == PhotonCandidate->GetTrackLabelPositive()) {
          PhotonCandidate->SetMCLabelPositive(TMath::Abs(tempDaughter->GetLabel()));
          PhotonCandidate->SetLabelPositive(i);
          AODLabelPos = kTRUE;
        }
      }
      if (!AODLabelNeg) {
        if (tempDaughter->GetID() == PhotonCandidate->GetTrackLabelNegative()) {
          PhotonCandidate->SetMCLabelNegative(TMath::Abs(tempDaughter->GetLabel()));
          PhotonCandidate->SetLabelNegative(i);
          AODLabelNeg = kTRUE;
        }
      }
      if (AODLabelNeg && AODLabelPos) {
        break;
      }
    }
    if (!AODLabelPos || !AODLabelNeg) {
      cout << "WARNING!!! AOD TRACKS NOT FOUND FOR" << endl;
      if (!AODLabelNeg) {
        PhotonCandidate->SetMCLabelNegative(-999999);
        PhotonCandidate->SetLabelNegative(-999999);
      }
      if (!AODLabelPos) {
        PhotonCandidate->SetMCLabelPositive(-999999);
        PhotonCandidate->SetLabelPositive(-999999);
      }
    }
  }
}

//________________________________________________________________________
void AliAnalysisTaskConvCaloTree::ResetBuffer()
{
  fBuffer_EventWeight = 1;
  fBuffer_Event_Vertex_Z = 0;
  fMCQ2 = 0.;
  fV0Mult = 0;
  fTriggerBits = 0;
}

//________________________________________________________________________
void AliAnalysisTaskConvCaloTree::ResetBufferVectors()
{

  fVBuffer_Cluster_E.clear();
  fVBuffer_Cluster_Eta.clear();
  fVBuffer_Cluster_Phi.clear();
  fVBuffer_Cluster_NCells.clear();
  fVBuffer_Cluster_M02.clear();
  fVBuffer_Cluster_Time.clear();
  fVBuffer_Cluster_Fcross.clear();
  fVBuffer_Cluster_Exotic.clear();
  fVBuffer_Cluster_CellTimes.clear();
  fVBuffer_Cluster_Isolated.clear();
  fVBuffer_Cluster_CellEnergies.clear();
  fVBuffer_Cluster_CellIDs.clear();
  fVTrueClusterPi0DaughterIndex.clear();
  fVTrueClusterEtaDaughterIndex.clear();
  fVTrueClusterSecondaryIndex.clear();
  fVTrueClusterMCId.clear();
  fVTrueClusterMCId2.clear();
  fVTrueClusterMCTrueEnergy.clear();
  fVTrueClusterMCTrueEnergy2.clear();
  fVTrueClusterMCIsMerged.clear();
  fVTrueClusterConvRadius.clear();

  fVBuffer_Conv_px.clear();
  fVBuffer_Conv_py.clear();
  fVBuffer_Conv_pz.clear();
  fVBuffer_Elec1etaCalo.clear();
  fVBuffer_Elec1phiCalo.clear();
  fVBuffer_Elec2etaCalo.clear();
  fVBuffer_Elec2phiCalo.clear();

  fVBuffer_Conv_R.clear();
  fVBuffer_Conv_PsiPair.clear();
  fVBuffer_Conv_NTPCClusElec1.clear();
  fVBuffer_Conv_NTPCClusElec2.clear();
  fVBuffer_Conv_dEdxElec1.clear();
  fVBuffer_Conv_dEdxElec2.clear();
  fVBuffer_Conv_PElec1.clear();
  fVBuffer_Conv_PElec2.clear();
  fVBuffer_Conv_Quality.clear();
  fVBuffer_Conv_CosPAngle.clear();
  fVBuffer_Conv_Chi2.clear();

  fVTrueConvPi0DaughterIndex.clear();
  fVTrueConvEtaDaughterIndex.clear();
  fVTrueConvSecDaughterIndex.clear();

  fVTrueConvMCTruePx.clear();
  fVTrueConvMCTruePy.clear();
  fVTrueConvMCTruePz.clear();
  fVTrueConvMCLabel.clear();

  fVBuffer_Track_px.clear();
  fVBuffer_Track_py.clear();
  fVBuffer_Track_pz.clear();
  fVBuffer_Track_P.clear();
  fVBuffer_Track_dedx.clear();
  fVBuffer_Track_TOFSignal.clear();
  fVBuffer_Track_DCA.clear();
  fVBuffer_Track_FracNClus.clear();
  fVBuffer_Track_PDG.clear();
  fVBuffer_Track_StackID.clear();
  fVBuffer_Track_Calo_eta.clear();
  fVBuffer_Track_Calo_phi.clear();

  fVBuffer_Jet_Pt.clear();
  fVBuffer_Jet_Eta.clear();
  fVBuffer_Jet_Phi.clear();
  fVBuffer_Jet_NNeutr.clear();
  fVBuffer_Jet_NCh.clear();
  fVBuffer_Jet_NEF.clear();

  fVBuffer_TrueJet_Pt.clear();
  fVBuffer_TrueJet_Eta.clear();
  fVBuffer_TrueJet_Phi.clear();
  fVBuffer_TrueJet_Parton_Pt.clear();

  fVBuffer_MCGenID.clear();
  fVBuffer_MCGenPDG.clear();
  fVBuffer_MCGenPx.clear();
  fVBuffer_MCGenPy.clear();
  fVBuffer_MCGenPz.clear();
  fVBuffer_MCGenMotherID.clear();
}
