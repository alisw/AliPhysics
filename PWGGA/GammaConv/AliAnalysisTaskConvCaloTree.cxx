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

class iostream;

using namespace std;

ClassImp(AliAnalysisTaskConvCaloTree)

//________________________________________________________________________
AliAnalysisTaskConvCaloTree::AliAnalysisTaskConvCaloTree() : AliAnalysisTaskSE(),
  fV0Reader(NULL),
  fV0ReaderName("V0ReaderV1"),
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
  fUseClusterIsolation(0),
  fSaveJets(0),
  fMinTrackPt(0),
  fInfoList(NULL),
  fHistoNEvents(NULL),
  fHistoMCPi0Pt(NULL),
  fHistoMCEtaPt(NULL),
  fBuffer_EventWeight(1),
  fBuffer_Event_Vertex_Z(0),
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
  fVBuffer_Conv_CosPAngle(),
  fVBuffer_Conv_Chi2(),
  fVTrueConvPi0DaughterIndex(),
  fVTrueConvEtaDaughterIndex(),
  fVTrueConvMCTruePx(),
  fVTrueConvMCTruePy(),
  fVTrueConvMCTruePz(),
  fVTrueConvMCLabel(),
  fVBuffer_Track_px(),
  fVBuffer_Track_py(),
  fVBuffer_Track_pz(),
  fVBuffer_Track_P(),
  fVBuffer_Track_Calo_eta(),
  fVBuffer_Track_Calo_phi(),
  fVBuffer_Jet_Pt(),
  fVBuffer_Jet_Eta(),
  fVBuffer_Jet_Phi()
{
}

AliAnalysisTaskConvCaloTree::AliAnalysisTaskConvCaloTree(const char *name) : AliAnalysisTaskSE(name),
  fV0Reader(NULL),
  fV0ReaderName("V0ReaderV1"),
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
  fUseClusterIsolation(0),
  fSaveJets(0),
  fMinTrackPt(0),
  fInfoList(0),
  fHistoNEvents(NULL),
  fHistoMCPi0Pt(NULL),
  fHistoMCEtaPt(NULL),
  fBuffer_EventWeight(1),
  fBuffer_Event_Vertex_Z(0),
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
  fVBuffer_Conv_CosPAngle(),
  fVBuffer_Conv_Chi2(),
  fVTrueConvPi0DaughterIndex(),
  fVTrueConvEtaDaughterIndex(),
  fVTrueConvMCTruePx(),
  fVTrueConvMCTruePy(),
  fVTrueConvMCTruePz(),
  fVTrueConvMCLabel(),
  fVBuffer_Track_px(),
  fVBuffer_Track_py(),
  fVBuffer_Track_pz(),
  fVBuffer_Track_P(),
  fVBuffer_Track_Calo_eta(),
  fVBuffer_Track_Calo_phi(),
  fVBuffer_Jet_Pt(),
  fVBuffer_Jet_Eta(),
  fVBuffer_Jet_Phi()
{
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
  if(!fV0Reader){printf("Error: No V0 Reader");return;}// GetV0Reader

  if(fSaveJets){
    fConvJetReader = (AliAnalysisTaskConvJet*)AliAnalysisManager::GetAnalysisManager()->GetTask("AliAnalysisTaskConvJet");
    if (!fConvJetReader) {
      printf("Error: No AliAnalysisTaskConvJet");
      return;
    }
  }

  //retrieve pointer to CaloIsolation Instance
  if(fUseClusterIsolation) fCaloIsolation = (AliPhotonIsolation*)AliAnalysisManager::GetAnalysisManager()->GetTask(fCaloIsolationName.Data());
  if(!fCaloIsolation && fUseClusterIsolation){ AliFatal("AliPhotonIsolation instance could not be initialized!");}

  if(fOutputList != NULL){
    delete fOutputList;
    fOutputList                         = NULL;
  }
  if(fOutputList == NULL){
    fOutputList                         = new TList();
    fOutputList->SetOwner(kTRUE);
  }

  if(((AliConvEventCuts*)fEventCuts)->GetCutHistograms()){
    fOutputList->Add(((AliConvEventCuts*)fEventCuts)->GetCutHistograms());
  }
  if(fClusterCutsEMC){
    if(((AliCaloPhotonCuts*)fClusterCutsEMC)->GetCutHistograms()){
      fOutputList->Add(((AliCaloPhotonCuts*)fClusterCutsEMC)->GetCutHistograms());
    }
  }
  if(fClusterCutsPHOS){
    if(((AliCaloPhotonCuts*)fClusterCutsPHOS)->GetCutHistograms()){
      fOutputList->Add(((AliCaloPhotonCuts*)fClusterCutsPHOS)->GetCutHistograms());
    }
  }

  fInfoList          = new TList();
  fInfoList->SetName(Form("%s event histograms", (fEventCuts->GetCutNumber()).Data()));
  fInfoList->SetOwner(kTRUE);


  if(fInfoList){
    fOutputList->Add(fInfoList);
  }
  PostData(1, fOutputList);

  fPhotonTree = new TTree(Form("ConvCaloPhotons_%s",(fEventCuts->GetCutNumber()).Data()),Form("ConvCaloPhotons_%s",(fEventCuts->GetCutNumber()).Data()));

  if(fSaveClusters)
  {
    if(fIsMC > 1){
      fPhotonTree->Branch("Event_Weight",                   &fBuffer_EventWeight,                     "Event_Weight/F");
    }
    fPhotonTree->Branch("Event_VertexZ",                    &fBuffer_Event_Vertex_Z,                  "Event_VertexZ/F");

    fPhotonTree->Branch("Cluster_E",                        &fVBuffer_Cluster_E);
    fPhotonTree->Branch("Cluster_Eta",                      &fVBuffer_Cluster_Eta);
    fPhotonTree->Branch("Cluster_Phi",                      &fVBuffer_Cluster_Phi);
    fPhotonTree->Branch("Cluster_NCells",                   &fVBuffer_Cluster_NCells);
    fPhotonTree->Branch("Cluster_M02",                      &fVBuffer_Cluster_M02);
    fPhotonTree->Branch("Cluster_Time",                     &fVBuffer_Cluster_Time);
    if(fClusterCutsEMC){
      fPhotonTree->Branch("Cluster_Fcross",                 &fVBuffer_Cluster_Fcross);
      fPhotonTree->Branch("Cluster_isExotic",               &fVBuffer_Cluster_Exotic);
      fPhotonTree->Branch("Cluster_CellEnergies",           &fVBuffer_Cluster_CellEnergies);
      fPhotonTree->Branch("Cluster_CellIDs",                &fVBuffer_Cluster_CellIDs);
      if(fIsMC == 0){
        fPhotonTree->Branch("Cluster_CellTimes",              &fVBuffer_Cluster_CellTimes);
      }
    }
    if(fUseClusterIsolation){
      fPhotonTree->Branch("Cluster_Isolated",               &fVBuffer_Cluster_Isolated);
    }
    if(fIsMC){
      fPhotonTree->Branch("Cluster_MotherPi0",      &fVTrueClusterPi0DaughterIndex);
      fPhotonTree->Branch("Cluster_MotherEta",      &fVTrueClusterEtaDaughterIndex);
      fPhotonTree->Branch("Cluster_MCId",           &fVTrueClusterMCId);
      fPhotonTree->Branch("Cluster_MCId2",          &fVTrueClusterMCId2);
      fPhotonTree->Branch("Cluster_MCTrueEnergy",   &fVTrueClusterMCTrueEnergy);
      fPhotonTree->Branch("Cluster_MCTrueEnergy2",  &fVTrueClusterMCTrueEnergy2);
      fPhotonTree->Branch("Cluster_MCTrueIsMerged", &fVTrueClusterMCIsMerged);
      fPhotonTree->Branch("Cluster_MCTrueConvR",    &fVTrueClusterConvRadius);
    }
  }
  if(fSaveConversions)
  {
    fPhotonTree->Branch("Conv_px",                          &fVBuffer_Conv_px);
    fPhotonTree->Branch("Conv_py",                          &fVBuffer_Conv_py);
    fPhotonTree->Branch("Conv_pz",                          &fVBuffer_Conv_pz);
    fPhotonTree->Branch("Conv_Elec1etaCalo",                &fVBuffer_Elec1etaCalo);
    fPhotonTree->Branch("Conv_Elec2etaCalo",                &fVBuffer_Elec2etaCalo);
    fPhotonTree->Branch("Conv_Elec1phiCalo",                &fVBuffer_Elec1phiCalo);
    fPhotonTree->Branch("Conv_Elec2phiCalo",                &fVBuffer_Elec2phiCalo);
    if(fSaveConversions > 1){
      fPhotonTree->Branch("Conv_Radius",                    &fVBuffer_Conv_R);
      fPhotonTree->Branch("Conv_PsiPair",                   &fVBuffer_Conv_PsiPair);
      fPhotonTree->Branch("Conv_NTPCClus_Elec1",            &fVBuffer_Conv_NTPCClusElec1);
      fPhotonTree->Branch("Conv_NTPCClus_Elec2",            &fVBuffer_Conv_NTPCClusElec2);
      fPhotonTree->Branch("Conv_dEdx_Elec1",                &fVBuffer_Conv_dEdxElec1);
      fPhotonTree->Branch("Conv_dEdx_Elec2",                &fVBuffer_Conv_dEdxElec2);
      fPhotonTree->Branch("Conv_P_Elec1",                   &fVBuffer_Conv_PElec1);
      fPhotonTree->Branch("Conv_P_Elec2",                   &fVBuffer_Conv_PElec2);
      fPhotonTree->Branch("Conv_CosPAngle",                 &fVBuffer_Conv_CosPAngle);
      fPhotonTree->Branch("Conv_Chi2NDF",                   &fVBuffer_Conv_Chi2);
    }
    if(fIsMC){
      fPhotonTree->Branch("Conv_MotherPi0",         &fVTrueConvPi0DaughterIndex);
      fPhotonTree->Branch("Conv_MotherEta",         &fVTrueConvEtaDaughterIndex);
      fPhotonTree->Branch("Conv_MCTruePx",          &fVTrueConvMCTruePx);
      fPhotonTree->Branch("Conv_MCTruePy",          &fVTrueConvMCTruePy);
      fPhotonTree->Branch("Conv_MCTruePz",          &fVTrueConvMCTruePz);
      fPhotonTree->Branch("Conv_MCLabel",           &fVTrueConvMCLabel);
    }
  }
  if(fSaveTracks)
  {
    if(fSaveTracks == 1){
      fPhotonTree->Branch("Track_P",                        &fVBuffer_Track_P);
    }
    else if(fSaveTracks == 2){
      fPhotonTree->Branch("Track_px",                       &fVBuffer_Track_px);
      fPhotonTree->Branch("Track_py",                       &fVBuffer_Track_py);
      fPhotonTree->Branch("Track_pz",                       &fVBuffer_Track_pz);
    }
    fPhotonTree->Branch("Track_TracketaCalo",               &fVBuffer_Track_Calo_eta);
    fPhotonTree->Branch("Track_TrackphiCalo",               &fVBuffer_Track_Calo_phi);
  }
  if(fSaveJets){
    fPhotonTree->Branch("Jet_Pt",                         &fVBuffer_Jet_Pt);
    fPhotonTree->Branch("Jet_Eta",                        &fVBuffer_Jet_Eta);
    fPhotonTree->Branch("Jet_Phi",                        &fVBuffer_Jet_Phi);
  }



  fHistoNEvents     = new TH1F("NEvents", "NEvents", 3, -0.5, 2.5);
  fHistoNEvents->GetXaxis()->SetBinLabel(1,"Accepted");
  fHistoNEvents->GetXaxis()->SetBinLabel(2,"Not Accepted");
  fHistoNEvents->GetXaxis()->SetBinLabel(3,"No jet in event");
  fInfoList->Add(fHistoNEvents);

  if(fIsMC > 0){
    fHistoMCPi0Pt = new TH1F("MC_Pi0_Pt", "MC_Pi0_Pt", 1000, 0, 100);
    fInfoList->Add(fHistoMCPi0Pt);

    fHistoMCEtaPt = new TH1F("MC_Eta_Pt", "MC_Eta_Pt", 1000, 0, 100);
    fInfoList->Add(fHistoMCEtaPt);
  }

  fV0Reader=(AliV0ReaderV1*)AliAnalysisManager::GetAnalysisManager()->GetTask(fV0ReaderName.Data());


  OpenFile(2);
  PostData(2, fPhotonTree);
}

//_____________________________________________________________________________
Bool_t AliAnalysisTaskConvCaloTree::Notify()
{
  if (fEventCuts->GetPeriodEnum() == AliConvEventCuts::kNoPeriod && ((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetPeriodEnum() != AliConvEventCuts::kNoPeriod){
    fEventCuts->SetPeriodEnumExplicit(((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetPeriodEnum());
  } else if (fEventCuts->GetPeriodEnum() == AliConvEventCuts::kNoPeriod ){
    fEventCuts->SetPeriodEnum(fV0Reader->GetPeriodName());
  }

  return kTRUE;
}

//________________________________________________________________________
void AliAnalysisTaskConvCaloTree::UserExec(Option_t *){

  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  fPIDResponse = (AliPIDResponse*)inputHandler->GetPIDResponse();

  fInputEvent                         = InputEvent();
  if(fIsMC>0) fMCEvent                = MCEvent();
  Int_t eventNotAccepted              = fEventCuts->IsEventAcceptedByCut(fV0Reader->GetEventCuts(),fInputEvent,fMCEvent,fIsHeavyIon,kFALSE);

  if(eventNotAccepted) {
    fHistoNEvents->Fill(1);
    return; // Check Centrality, PileUp, SDD and V0AND --> Not Accepted => eventQuality = 1
  }

  const AliVVertex* primVtxMC   = fInputEvent->GetPrimaryVertex();
  fBuffer_Event_Vertex_Z = primVtxMC->GetZ();

  if(fIsMC>0 && fInputEvent->IsA()==AliAODEvent::Class() && !(fV0Reader->AreAODsRelabeled())){
    RelabelAODPhotonCandidates(kTRUE);// In case of AODMC relabeling MC
    fV0Reader->RelabelAODs(kTRUE);
  }

  fWeightJetJetMC = 1;
  if(fIsMC > 1){
    Float_t pthard = -1;
    Bool_t isMCJet = ((AliConvEventCuts*)fEventCuts)->IsJetJetMCEventAccepted( fMCEvent, fWeightJetJetMC,pthard, fInputEvent );
    if (!isMCJet){
      return;
    }
    fBuffer_EventWeight = fWeightJetJetMC;
  }

  if(fSaveJets){
    InitJets();
    if(fSaveJets == 2){
      if(fVBuffer_Jet_Pt.size() == 0){
        fHistoNEvents->Fill(2);
        return;
      }
    }
  }

  fHistoNEvents->Fill(0);

  if(fSaveClusters)
  {
    ProcessClustersAOD();
  }
  if(fSaveConversions)
  {
    ProcessConversionsAOD();
  }
  if(fSaveTracks > 0)
  {
    ProcessTracksAOD();
  }
  if(fIsMC > 0 && fSaveMCInformation)
  {
    ProcessAODMCParticles();
  }
  fPhotonTree->Fill();

  ResetBufferVectors();

  if(fIsMC>0 && fInputEvent->IsA()==AliAODEvent::Class() && !(fV0Reader->AreAODsRelabeled())){
    RelabelAODPhotonCandidates(kFALSE); // Back to ESDMC Label
    fV0Reader->RelabelAODs(kFALSE);
  }

  PostData(1, fOutputList);
}

///________________________________________________________________________
void AliAnalysisTaskConvCaloTree::ProcessClustersAOD(){

  Int_t nclus                       = 0;
  if(!fCorrTaskSetting.CompareTo("")){
    nclus = fInputEvent->GetNumberOfCaloClusters();
  } else {
    if(!farrClustersProcess) farrClustersProcess = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(Form("%sClustersBranch",fCorrTaskSetting.Data())));
    if(!farrClustersProcess)
      AliFatal(Form("%sClustersBranch was not found in AliAnalysisTaskGammaCalo! Check the correction framework settings!",fCorrTaskSetting.Data()));
    nclus = farrClustersProcess->GetEntries();
  }

  if(nclus == 0)  return;

  AliVCluster* clus = NULL;
  for(Long_t i = 0; i < nclus; i++){

    if(fInputEvent->IsA()==AliAODEvent::Class()){
      if(farrClustersProcess)
        clus = new AliAODCaloCluster(*(AliAODCaloCluster*)farrClustersProcess->At(i));
      else
        clus = new AliAODCaloCluster(*(AliAODCaloCluster*)fInputEvent->GetCaloCluster(i));
    } else if(fInputEvent->IsA()==AliESDEvent::Class()){
      AliFatal("ESD running is not supported");
    }
    if(!clus) continue;

    Bool_t clusIsAcc = kFALSE;
    if(fClusterCutsEMC && clus->IsEMCAL()){
      if(((AliCaloPhotonCuts*)fClusterCutsEMC)->ClusterIsSelected(clus,fInputEvent,fMCEvent,fIsMC, fWeightJetJetMC,i)){
        clusIsAcc = kTRUE;
        fClusterCuts = fClusterCutsEMC;
      }
    }
    if(fClusterCutsPHOS && clus->IsPHOS() && !clusIsAcc){
      if(((AliCaloPhotonCuts*)fClusterCutsPHOS)->ClusterIsSelected(clus,fInputEvent,fMCEvent,fIsMC, fWeightJetJetMC,i)){
        clusIsAcc = kTRUE;
        fClusterCuts = fClusterCutsPHOS;
      }
    }
    if(!clusIsAcc) {
      delete clus;
      continue;
    }

    Double_t vertex[3] = {0};
    InputEvent()->GetPrimaryVertex()->GetXYZ(vertex);
    TLorentzVector clusterVector;
    clus->GetMomentum(clusterVector,vertex);
    TLorentzVector tmpvec;
    tmpvec.SetPxPyPzE(clusterVector.Px(),clusterVector.Py(),clusterVector.Pz(),clusterVector.E());
    AliAODConversionPhoton *PhotonCandidate=new AliAODConversionPhoton(&tmpvec);

    Float_t etaCluster                     = (Float_t)PhotonCandidate->Eta();
    Float_t phiCluster                     = (Float_t)PhotonCandidate->Phi();
    if (phiCluster < 0) phiCluster          += 2*TMath::Pi();

    fVBuffer_Cluster_E.push_back(clus->E());
    fVBuffer_Cluster_Eta.push_back(etaCluster);
    fVBuffer_Cluster_Phi.push_back(phiCluster);
    fVBuffer_Cluster_NCells.push_back((Short_t)clus->GetNCells());
    fVBuffer_Cluster_M02.push_back((Short_t)(clus->GetM02()*1000));
    fVBuffer_Cluster_Time.push_back((Short_t)(clus->GetTOF()*1e9)); // time in ns
    // calculate Ecross
    if(fClusterCutsEMC){
      AliVCaloCells* cells    = fInputEvent->GetEMCALCells();
      Int_t largestCellID     = ((AliCaloPhotonCuts*)fClusterCutsEMC)->FindLargestCellInCluster(clus,fInputEvent);
      Float_t ecell1          = cells->GetCellAmplitude(largestCellID);
      Float_t eCross          = ((AliCaloPhotonCuts*)fClusterCutsEMC)->GetECross(largestCellID,cells);
      fVBuffer_Cluster_Fcross.push_back((Short_t)(1000*(1-eCross/ecell1)));
      fVBuffer_Cluster_Exotic.push_back(clus->GetIsExotic());


      Int_t nCellCluster = clus->GetNCells();
      for(Int_t iCell=0;iCell<nCellCluster;iCell++){
        Int_t cellID = clus->GetCellAbsId(iCell);
        Double_t cellAmp = cells->GetCellAmplitude(cellID);
        fVBuffer_Cluster_CellEnergies.push_back(cellAmp);
        fVBuffer_Cluster_CellIDs.push_back(cellID);
        if(fIsMC == 0){
          Double_t cellTime = cells->GetCellTime(cellID);
          fVBuffer_Cluster_CellTimes.push_back(cellTime*1e9);
        }
      }

    }

    if(fUseClusterIsolation && fCaloIsolation){
      Float_t ClusterPt = PhotonCandidate->Pt();
      Float_t pTCone = 0.1*ClusterPt; // momPercantage = 0.1
      fVBuffer_Cluster_Isolated.push_back(fCaloIsolation->GetIsolation(i,0.4,pTCone));
    }



    if(fIsMC > 0){// MC info


      if(!PhotonCandidate){
        fVTrueClusterPi0DaughterIndex.push_back(0);
        fVTrueClusterEtaDaughterIndex.push_back(0);
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
      PhotonCandidate->SetLeadingCellID(((AliCaloPhotonCuts*)fClusterCuts)->FindLargestCellInCluster(clus,fInputEvent));
      Int_t* mclabelsCluster = clus->GetLabels();
      PhotonCandidate->SetNCaloPhotonMCLabels(clus->GetNLabels());

      if (clus->GetNLabels()>0){
        for (Int_t k =0; k< (Int_t)clus->GetNLabels(); k++){
          PhotonCandidate->SetCaloPhotonMCLabel(k,mclabelsCluster[k]);
        }
      }

      if(!fAODMCTrackArray) fAODMCTrackArray = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
      if(!fAODMCTrackArray) return;

      PhotonCandidate->SetCaloPhotonMCFlagsAOD(fAODMCTrackArray, kTRUE);

      // get true MC id labels

      const int NLabels = 2;
      std::vector<Short_t> arrMCID = {0, 0};
      std::vector<Float_t> arrMCEnergy = {0, 0};
      for(Int_t i = 0; i < NLabels; ++i){
        Int_t gammaMCLabel          = PhotonCandidate->GetCaloPhotonMCLabel(i);   // get most probable MC label
        if(gammaMCLabel != -1){
          AliAODMCParticle *gammaMC = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(gammaMCLabel));
          arrMCID[i] = (Short_t) gammaMC->GetPdgCode();
          arrMCEnergy[i] = (Float_t) gammaMC->E();
        }
      }
      fVTrueClusterMCId.push_back(arrMCID[0]);
      fVTrueClusterMCId2.push_back(arrMCID[1]);
      fVTrueClusterMCTrueEnergy.push_back(arrMCEnergy[0]);
      fVTrueClusterMCTrueEnergy2.push_back(arrMCEnergy[1]);

      Int_t gammaMCLabel          = PhotonCandidate->GetCaloPhotonMCLabel(0);   // get most probable MC label
      Int_t gammaMotherLabel      = -1;
      Double_t convRadiusConvClus = 0;

      AliAODMCParticle * gammaMC = 0x0;
      if(gammaMCLabel != -1){
        gammaMC = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(gammaMCLabel));
        if (PhotonCandidate->IsLargestComponentPhoton() || PhotonCandidate->IsLargestComponentElectron()){    // largest component is electron
          if (PhotonCandidate->IsLargestComponentPhoton()){                                        // for photons its the direct mother
              gammaMotherLabel=gammaMC->GetMother();
          } else if (PhotonCandidate->IsLargestComponentElectron()){                               // for electrons its either the direct mother or for conversions the grandmother
            if (PhotonCandidate->IsConversion()){
              AliAODMCParticle * gammaGrandMotherMC =  static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(gammaMC->GetMother()));
              gammaMotherLabel=gammaGrandMotherMC->GetMother();
              // get the production vertex of the conversion electron
              convRadiusConvClus = sqrt(gammaMC->Xv()*gammaMC->Xv() + gammaMC->Yv()*gammaMC->Yv());
            } else gammaMotherLabel=gammaMC->GetMother();
          }
        }
        if(PhotonCandidate->IsMerged()){
          fVTrueClusterMCIsMerged.push_back(kTRUE);
        } else {
          fVTrueClusterMCIsMerged.push_back(kFALSE);
        }
      } else {
        fVTrueClusterMCIsMerged.push_back(kFALSE);
      }
      // conv radius of conv cluster
      fVTrueClusterConvRadius.push_back(static_cast<Short_t>(convRadiusConvClus));

      if(gammaMotherLabel > -1){
        if(((AliAODMCParticle*)fAODMCTrackArray->At(gammaMotherLabel))->GetPdgCode() == 111){
          fVTrueClusterPi0DaughterIndex.push_back(gammaMotherLabel);
          fVTrueClusterEtaDaughterIndex.push_back(0);
        }
        else if(((AliAODMCParticle*)fAODMCTrackArray->At(gammaMotherLabel))->GetPdgCode() == 221){
          fVTrueClusterPi0DaughterIndex.push_back(0);
          fVTrueClusterEtaDaughterIndex.push_back(gammaMotherLabel);
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
void AliAnalysisTaskConvCaloTree::ProcessConversionsAOD(){

  fReaderGammas = fV0Reader->GetReconstructedGammas(); // Gammas from default Cut
  for(Int_t i = 0; i < fReaderGammas->GetEntriesFast(); i++){
    AliAODConversionPhoton* PhotonCandidate = (AliAODConversionPhoton*) fReaderGammas->At(i);
    if(!PhotonCandidate) continue;

    // apply cuts
    if(!fConversionCuts) AliFatal("Conversion Cuts Not Initialized");
    if(!((AliConversionPhotonCuts*)fConversionCuts)->PhotonIsSelected(PhotonCandidate,fInputEvent)) continue;

    fVBuffer_Conv_px.push_back(PhotonCandidate->GetPx());
    fVBuffer_Conv_py.push_back(PhotonCandidate->GetPy());
    fVBuffer_Conv_pz.push_back(PhotonCandidate->GetPz());

    if(fSaveConversions > 1){
      fVBuffer_Conv_R.push_back(static_cast<UShort_t>(PhotonCandidate->GetConversionRadius()*10));  // conv radius *10
      fVBuffer_Conv_PsiPair.push_back(static_cast<Short_t>(PhotonCandidate->GetPsiPair()*1000));    // psi pair *1000
      fVBuffer_Conv_CosPAngle.push_back(static_cast<Short_t>(((AliConversionPhotonCuts*)fConversionCuts)->GetCosineOfPointingAngle(PhotonCandidate, fInputEvent)*1000));  // pointing angle *1000
      fVBuffer_Conv_Chi2.push_back(static_cast<UShort_t>(PhotonCandidate->GetChi2perNDF()*10));  // pointing angle *1000
      
    }

    //get track eta and phi on Calo surface
    for (Int_t iElec = 0;iElec < 2;iElec++){
      Int_t tracklabel = PhotonCandidate->GetLabel(iElec);
      AliVTrack *inTrack = 0x0;
      if( fInputEvent->IsA()==AliESDEvent::Class() ) {
        AliFatal("ESD running is not supported");
      } else { ///AOD
        if( ((AliV0ReaderV1*)AliAnalysisManager::GetAnalysisManager()->GetTask(fV0ReaderName.Data()))->AreAODsRelabeled() ){
          inTrack = dynamic_cast<AliVTrack*>(fInputEvent->GetTrack(tracklabel));
        } else {
          for( Int_t ii=0;ii<fInputEvent->GetNumberOfTracks();ii++ ) {
            inTrack = dynamic_cast<AliVTrack*>(fInputEvent->GetTrack(ii));
            if(inTrack){
              if(inTrack->GetID() == tracklabel) {
                break;
              }
            }
          }
        }

        if(fSaveConversions > 1){
          // Short_t Charge    = inTrack->Charge();
          Double_t electronNSigmaTPC = fPIDResponse->NumberOfSigmasTPC(inTrack,AliPID::kElectron);
          // Double_t P=inTrack->P();;
          // Double_t Eta=inTrack->Eta();
          // Double_t nSigdEdxCorr = ((AliConversionPhotonCuts*)fConversionCuts)->GetCorrectedElectronTPCResponse(Charge,electronNSigmaTPC,P,Eta,inTrack->GetTPCNcls(),PhotonCandidate->GetConversionRadius());
          if(iElec == 0){
            fVBuffer_Conv_NTPCClusElec1.push_back(static_cast<UShort_t>(inTrack->GetTPCNcls())); // NTPC clus
            fVBuffer_Conv_dEdxElec1.push_back(static_cast<Short_t>(electronNSigmaTPC*100)); // dedx
            fVBuffer_Conv_PElec1.push_back(static_cast<Short_t>(inTrack->P()*100)); // momentum
          } else {
            fVBuffer_Conv_NTPCClusElec2.push_back(static_cast<UShort_t>(inTrack->GetTPCNcls())); // NTPC clus
            fVBuffer_Conv_dEdxElec2.push_back(static_cast<Short_t>(electronNSigmaTPC*100)); // dedx
            fVBuffer_Conv_PElec2.push_back(static_cast<Short_t>(inTrack->P()*100)); // momentum

          }
        }


        // track extrapolation
        Double_t xyz[3] = {0}, pxpypz[3] = {0}, cv[21] = {0};
        inTrack->GetPxPyPz(pxpypz);
        inTrack->GetXYZ(xyz);
        inTrack->GetCovarianceXYZPxPyPz(cv);

        AliExternalTrackParam emcParam(xyz,pxpypz,cv,inTrack->Charge());
        Float_t eta, phi, pt;

        //propagate tracks to emc surfaces
        if (!AliEMCALRecoUtils::ExtrapolateTrackToEMCalSurface(&emcParam, 440., 0.139, 20., eta, phi, pt)) {
          if(iElec == 0){
            fVBuffer_Elec1etaCalo.push_back(eta);
            fVBuffer_Elec1phiCalo.push_back(phi);
          } else {
            fVBuffer_Elec2etaCalo.push_back(eta);
            fVBuffer_Elec2phiCalo.push_back(phi);
          }
          continue;
        }
        if(iElec == 0){
          fVBuffer_Elec1etaCalo.push_back(eta);
          fVBuffer_Elec1phiCalo.push_back(phi);
        } else {
          fVBuffer_Elec2etaCalo.push_back(eta);
          fVBuffer_Elec2phiCalo.push_back(phi);
        }
      }
    }

    if(fIsMC > 0){
      if(!fAODMCTrackArray) fAODMCTrackArray = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
      if (fAODMCTrackArray != NULL && PhotonCandidate != NULL){

        AliAODMCParticle *posDaughter = (AliAODMCParticle*) fAODMCTrackArray->At(PhotonCandidate->GetMCLabelPositive());

        AliAODMCParticle *Photon = (AliAODMCParticle*) fAODMCTrackArray->At(posDaughter->GetMother());

        fVTrueConvMCLabel.push_back(Photon->GetLabel());

        if(Photon->GetPdgCode() != 22){
          fVTrueConvEtaDaughterIndex.push_back(-1);
          fVTrueConvPi0DaughterIndex.push_back(-1);
          fVTrueConvMCTruePx.push_back(0.);
          fVTrueConvMCTruePy.push_back(0.);
          fVTrueConvMCTruePz.push_back(0.);
          return; // Mother is no Photon
        }
        // fill true px, py and pz
        fVTrueConvMCTruePx.push_back(Photon->Px());
        fVTrueConvMCTruePy.push_back(Photon->Py());
        fVTrueConvMCTruePz.push_back(Photon->Pz());

        Int_t gammaMotherLabel          = Photon->GetMother();   // get the mother stack ID
        if ( ((AliAODMCParticle*)fAODMCTrackArray->At(Photon->GetMother()))->GetPdgCode() == 111) {
          fVTrueConvPi0DaughterIndex.push_back(gammaMotherLabel);
          fVTrueConvEtaDaughterIndex.push_back(-1);
        } else if ( ((AliAODMCParticle*)fAODMCTrackArray->At(Photon->GetMother()))->GetPdgCode() == 221) {
          fVTrueConvEtaDaughterIndex.push_back(gammaMotherLabel);
          fVTrueConvPi0DaughterIndex.push_back(-1);
        } else {
          fVTrueConvEtaDaughterIndex.push_back(-1);
          fVTrueConvPi0DaughterIndex.push_back(-1);
        }
      }
    }
  }
}

//________________________________________________________________________
void AliAnalysisTaskConvCaloTree::ProcessTracksAOD(){



  AliAODEvent *aodev = dynamic_cast<AliAODEvent*>(fInputEvent);
  if (!aodev) {
    AliError("Task needs AODs , returning");
    return;
  }

  for (Int_t itr=0;itr<fInputEvent->GetNumberOfTracks();itr++){


    AliVTrack* inTrack = dynamic_cast<AliVTrack*>(aodev->GetTrack(itr));
    if(!inTrack) continue;
    if(inTrack->Pt()<fMinTrackPt) continue;
    AliAODTrack *aodt = dynamic_cast<AliAODTrack*>(inTrack);
    // check track cuts
    if(!aodt->IsHybridGlobalConstrainedGlobal()) continue;

    if(fSaveTracks == 2){
      fVBuffer_Track_px.push_back(static_cast<Short_t>(aodt->Px()*100));
      fVBuffer_Track_py.push_back(static_cast<Short_t>(aodt->Py()*100));
      fVBuffer_Track_pz.push_back(static_cast<Short_t>(aodt->Pz()*100));
    }
    // This is not the proper extrapolation! Although for a rough matching of tracks and clusyers it was found to save a lot of CPU time
    Float_t trackEta = aodt->GetTrackEtaOnEMCal();
    Float_t trackPhi = aodt->GetTrackPhiOnEMCal();
    if(trackEta == -999 || trackPhi == -999){
      // maximum values for Short_t (--> not matched to Calo)
      // Store this information if all tracks are stored
      if(fSaveTracks == 2){
        fVBuffer_Track_Calo_eta.push_back(static_cast<Short_t>(32767));
        fVBuffer_Track_Calo_phi.push_back(static_cast<UShort_t>(65535));
      }
      continue;
    }
    if(fSaveTracks == 1){
      fVBuffer_Track_P.push_back(static_cast<Short_t>(aodt->P()*100));
    }
    fVBuffer_Track_Calo_eta.push_back(static_cast<Short_t>(trackEta*10000));
    fVBuffer_Track_Calo_phi.push_back(static_cast<UShort_t>(trackPhi*10000));

  }
}

//________________________________________________________________________
void AliAnalysisTaskConvCaloTree::ProcessAODMCParticles()
{
  if(fIsMC == 0) return;
  const AliVVertex* primVtxMC   = fMCEvent->GetPrimaryVertex();
  Double_t mcProdVtxX   = primVtxMC->GetX();
  Double_t mcProdVtxY   = primVtxMC->GetY();
  Double_t mcProdVtxZ   = primVtxMC->GetZ();

  if(!fAODMCTrackArray) fAODMCTrackArray = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
  if (fAODMCTrackArray == NULL) return;
  // Loop over all primary MC particle
  for(Long_t i = 0; i < fAODMCTrackArray->GetEntriesFast(); i++) {
    Double_t tempParticleWeight       = fWeightJetJetMC;
    AliAODMCParticle* particle = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(i));

    if (!particle) continue;

    Bool_t isPrimary = fEventCuts->IsConversionPrimaryAOD(fInputEvent, particle, mcProdVtxX, mcProdVtxY, mcProdVtxZ);
    if (isPrimary){
      Int_t isMCFromMBHeader = -1;

      if(fEventCuts->GetSignalRejection() != 0){
        isMCFromMBHeader = fEventCuts->IsParticleFromBGEvent(i, fMCEvent, fInputEvent);
        if(isMCFromMBHeader == 0 && fEventCuts->GetSignalRejection() != 3) continue;
        // Set the jetjet weight to 1 in case the particle orignated from the minimum bias header
        if(isMCFromMBHeader == 2 && fEventCuts->GetSignalRejection() == 4) tempParticleWeight = 1;
      }

      if(fMesonCuts->MesonIsSelectedAODMC(particle,fAODMCTrackArray,fEventCuts->GetEtaShift())){
        Float_t weighted= 1;

        if(fEventCuts->IsParticleFromBGEvent(i, fMCEvent, fInputEvent)){
          if (particle->Pt()>0.005){
            weighted = fEventCuts->GetWeightForMeson(i, 0x0, fInputEvent);
          }
        }

        if(particle->GetPdgCode() == 111){
          fHistoMCPi0Pt->Fill(particle->Pt(),weighted* tempParticleWeight); // All MC Pi0
        } else if(particle->GetPdgCode() == 221){
          fHistoMCEtaPt->Fill(particle->Pt(),weighted* tempParticleWeight); // All MC Eta
        }
      }
    }
  }
}

//________________________________________________________________________
void AliAnalysisTaskConvCaloTree::InitJets()
{
  for(unsigned int i = 0; i < (fConvJetReader->GetVectorJetPt()).size(); ++i){
    fVBuffer_Jet_Pt.push_back((fConvJetReader->GetVectorJetPt())[i]);
    fVBuffer_Jet_Eta.push_back((fConvJetReader->GetVectorJetEta())[i]);
    fVBuffer_Jet_Phi.push_back((fConvJetReader->GetVectorJetPhi())[i]);
  }
}

//________________________________________________________________________
void AliAnalysisTaskConvCaloTree::Terminate(Option_t *){

}

//________________________________________________________________________
void AliAnalysisTaskConvCaloTree::RelabelAODPhotonCandidates(Bool_t mode){

  // Relabeling For AOD Event
  // ESDiD -> AODiD
  // MCLabel -> AODMCLabel
  if(mode){
    fMCEventPos.clear();
    fMCEventNeg.clear();
    fESDArrayPos.clear();
    fESDArrayNeg.clear();
    fMCEventPos.resize(fReaderGammas->GetEntries());
    fMCEventNeg.resize(fReaderGammas->GetEntries());
    fESDArrayPos.resize(fReaderGammas->GetEntries());
    fESDArrayNeg.resize(fReaderGammas->GetEntries());
  }

  for(Int_t iGamma = 0;iGamma<fReaderGammas->GetEntries();iGamma++){
    AliAODConversionPhoton* PhotonCandidate = (AliAODConversionPhoton*) fReaderGammas->At(iGamma);
    if(!PhotonCandidate) continue;
    if(!mode){// Back to ESD Labels
    PhotonCandidate->SetMCLabelPositive(fMCEventPos[iGamma]);
    PhotonCandidate->SetMCLabelNegative(fMCEventNeg[iGamma]);
    PhotonCandidate->SetLabelPositive(fESDArrayPos[iGamma]);
    PhotonCandidate->SetLabelNegative(fESDArrayNeg[iGamma]);
    continue;
    }
    fMCEventPos[iGamma] =  PhotonCandidate->GetMCLabelPositive();
    fMCEventNeg[iGamma] =  PhotonCandidate->GetMCLabelNegative();
    fESDArrayPos[iGamma] = PhotonCandidate->GetTrackLabelPositive();
    fESDArrayNeg[iGamma] = PhotonCandidate->GetTrackLabelNegative();

    Bool_t AODLabelPos = kFALSE;
    Bool_t AODLabelNeg = kFALSE;

    for(Int_t i = 0; i<fInputEvent->GetNumberOfTracks();i++){
      AliAODTrack *tempDaughter = static_cast<AliAODTrack*>(fInputEvent->GetTrack(i));
      if(!AODLabelPos){
        if( tempDaughter->GetID() == PhotonCandidate->GetTrackLabelPositive() ){
        PhotonCandidate->SetMCLabelPositive(TMath::Abs(tempDaughter->GetLabel()));
        PhotonCandidate->SetLabelPositive(i);
        AODLabelPos = kTRUE;
        }
      }
      if(!AODLabelNeg){
        if( tempDaughter->GetID() == PhotonCandidate->GetTrackLabelNegative()){
        PhotonCandidate->SetMCLabelNegative(TMath::Abs(tempDaughter->GetLabel()));
        PhotonCandidate->SetLabelNegative(i);
        AODLabelNeg = kTRUE;
        }
      }
      if(AODLabelNeg && AODLabelPos){
        break;
      }
    }
    if(!AODLabelPos || !AODLabelNeg){
      cout<<"WARNING!!! AOD TRACKS NOT FOUND FOR"<<endl;
      if(!AODLabelNeg){
        PhotonCandidate->SetMCLabelNegative(-999999);
        PhotonCandidate->SetLabelNegative(-999999);
      }
      if(!AODLabelPos){
        PhotonCandidate->SetMCLabelPositive(-999999);
        PhotonCandidate->SetLabelPositive(-999999);
      }
    }
  }
}

//________________________________________________________________________
void AliAnalysisTaskConvCaloTree::ResetBuffer(){
  fBuffer_EventWeight                     = 1;
  fBuffer_Event_Vertex_Z                  = 0;
}

//________________________________________________________________________
void AliAnalysisTaskConvCaloTree::ResetBufferVectors(){

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
  fVBuffer_Conv_CosPAngle.clear();
  fVBuffer_Conv_Chi2.clear();

  fVTrueConvPi0DaughterIndex.clear();
  fVTrueConvEtaDaughterIndex.clear();

  fVTrueConvMCTruePx.clear();
  fVTrueConvMCTruePy.clear();
  fVTrueConvMCTruePz.clear();
  fVTrueConvMCLabel.clear();

  fVBuffer_Track_px.clear();
  fVBuffer_Track_py.clear();
  fVBuffer_Track_pz.clear();
  fVBuffer_Track_P.clear();
  fVBuffer_Track_Calo_eta.clear();
  fVBuffer_Track_Calo_phi.clear();

  fVBuffer_Jet_Pt.clear();
  fVBuffer_Jet_Eta.clear();
  fVBuffer_Jet_Phi.clear();

}
