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

/* AliAnalysisTaskAO2Dconverter
 *
 * Convert Run 2 ESDs to Run 3 prototype AODs (AliAO2D.root).
 */

#include <TChain.h>
#include <TTree.h>
#include <TMath.h>
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliEMCALGeometry.h"
#include "AliAnalysisTaskAO2Dconverter.h"
#include "AliVHeader.h"
#include "AliAnalysisManager.h"

#include "AliESDCaloCells.h"
#include "AliESDCaloTrigger.h"
#include "AliESDHeader.h"
#include "AliESDtrack.h"
#include "AliESDMuonTrack.h"
#include "AliESDMuonCluster.h"
#include "AliESDZDC.h"
#include "AliESDVZERO.h"
#include "AliESDv0.h"
#include "AliESDcascade.h"

#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
#include "AliPIDResponse.h"

#include "AliGenCocktailEventHeader.h"
#include "AliGenDPMjetEventHeader.h"
#include "AliGenEpos3EventHeader.h"
#include "AliGenEposEventHeader.h"
#include "AliGenEventHeader.h"
#include "AliGenEventHeaderTunedPbPb.h"
#include "AliGenGeVSimEventHeader.h"
#include "AliGenHepMCEventHeader.h"
#include "AliGenHerwigEventHeader.h"
#include "AliGenHijingEventHeader.h"
#include "AliGenPythiaEventHeader.h"
#include "AliGenToyEventHeader.h"

ClassImp(AliAnalysisTaskAO2Dconverter);

namespace
{

ULong64_t GetEventIdAsLong(AliVHeader *header)
{
  // return ((ULong64_t)header->GetBunchCrossNumber() +
  //         (ULong64_t)header->GetOrbitNumber() * 3564 +
  //         (ULong64_t)header->GetPeriodNumber() * 16777215 * 3564);

  ULong64_t lbc = (ULong64_t) header->GetBunchCrossNumber();      // The lowest 12 bits 
  ULong64_t lorb = ((ULong64_t) header->GetOrbitNumber()) << 12;  // The next 24 bits
  ULong64_t lper = ((ULong64_t) header->GetPeriodNumber()) << 36; // The last 28 bits 

  // In Run 3 we have only BC (12 bits) and orbit number (32 instead of 24 bits)
  ULong64_t id = lbc | lorb | lper;

  return id;
}

} // namespace

AliAnalysisTaskAO2Dconverter::AliAnalysisTaskAO2Dconverter(const char* name)
    : AliAnalysisTaskSE(name)
    , fTrackFilter(Form("AO2Dconverter%s", name), Form("fTrackFilter%s", name))
    , fEventCuts{}
{
  DefineInput(0, TChain::Class());
  for (Int_t i = 0; i < kTrees; i++) {
    fTreeStatus[i] = kTRUE;
    DefineOutput(1 + i, TTree::Class());
  }
}

AliAnalysisTaskAO2Dconverter::~AliAnalysisTaskAO2Dconverter()
{
  for (Int_t i = 0; i < kTrees; i++)
    if (fTree[i])
      delete fTree[i];
}

const TString AliAnalysisTaskAO2Dconverter::TreeName[kTrees] = { "O2collisions", "O2tracks", "O2calo",  "O2caloTrigger", "O2muon", "O2muoncls", "O2zdc", "O2vzero", "O2v0s", "O2cascades", "O2tof", "O2kine", "O2mcvtx", "O2range", "O2labels", "O2trigger" };

const TString AliAnalysisTaskAO2Dconverter::TreeTitle[kTrees] = { "Collision tree", "Barrel tracks", "Calorimeter cells", "Calorimeter triggers", "MUON tracks", "MUON clusters", "ZDC", "VZERO", "V0s", "Cascades", "TOF hits", "Kinematics", "MC vertex", "Range of MC labels", "MC labels", "Trigger info"};

const TClass* AliAnalysisTaskAO2Dconverter::Generator[kGenerators] = { AliGenEventHeader::Class(), AliGenCocktailEventHeader::Class(), AliGenDPMjetEventHeader::Class(), AliGenEpos3EventHeader::Class(), AliGenEposEventHeader::Class(), AliGenEventHeaderTunedPbPb::Class(), AliGenGeVSimEventHeader::Class(), AliGenHepMCEventHeader::Class(), AliGenHerwigEventHeader::Class(), AliGenHijingEventHeader::Class(), AliGenPythiaEventHeader::Class(), AliGenToyEventHeader::Class() };

TTree* AliAnalysisTaskAO2Dconverter::CreateTree(TreeIndex t)
{
  fTree[t] = new TTree(TreeName[t], TreeTitle[t]);
  // if (fTreeStatus[t])
  //   fTree[t]->Branch("fGlobalBC", &vtx.fGlobalBC, "fGlobalBC/l"); // Branch common to all trees
  return fTree[t];
}

void AliAnalysisTaskAO2Dconverter::PostTree(TreeIndex t)
{
  if (!fTreeStatus[t])
    return;
  PostData(t + 1, fTree[t]);
}

void AliAnalysisTaskAO2Dconverter::FillTree(TreeIndex t)
{
  if (!fTreeStatus[t])
    return;
  fTree[t]->Fill();
}

void AliAnalysisTaskAO2Dconverter::UserCreateOutputObjects()
{
  switch (fTaskMode) { // Setting active/inactive containers based on the TaskMode
  case kStandard:
    DisableTree(kKinematics);
    DisableTree(kMCvtx);
    break;
  default:
    break;
  }

  // Reset the offsets
  fOffsetMuTrackID = 0;
  fOffsetTrackID = 0;
  fOffsetV0ID = 0;
  fOffsetLabel = 0;

  // create output objects
  OpenFile(1); // Necessary for large outputs

  // Associate branches for fEventTree
  TTree* tEvents = CreateTree(kEvents);
  tEvents->SetAutoFlush(fNumberOfEventsPerCluster);
  if (fTreeStatus[kEvents]) {
    TString sstart = TString::Format("fStart[%d]/I", kTrees);
    TString sentries = TString::Format("fNentries[%d]/I", kTrees);
    tEvents->Branch("fRunNumber", &vtx.fRunNumber, "fRunNumber/I");
    tEvents->Branch("fStart", vtx.fStart, sstart.Data());
    tEvents->Branch("fNentries", vtx.fNentries, sentries.Data());
    tEvents->Branch("fGlobalBC", &vtx.fGlobalBC, "fGlobalBC/l");
    tEvents->Branch("fX", &vtx.fX, "fX/F");
    tEvents->Branch("fY", &vtx.fY, "fY/F");
    tEvents->Branch("fZ", &vtx.fZ, "fZ/F");
    tEvents->Branch("fCovXX", &vtx.fCovXX, "fCovXX/F");
    tEvents->Branch("fCovXY", &vtx.fCovXY, "fCovXY/F");
    tEvents->Branch("fCovXZ", &vtx.fCovXZ, "fCovXZ/F");
    tEvents->Branch("fCovYY", &vtx.fCovYY, "fCovYY/F");
    tEvents->Branch("fCovYZ", &vtx.fCovYZ, "fCovYZ/F");
    tEvents->Branch("fCovZZ", &vtx.fCovZZ, "fCovZZ/F");
    tEvents->Branch("fChi2", &vtx.fChi2, "fChi2/F");
    tEvents->Branch("fN", &vtx.fN, "fN/i");
    tEvents->Branch("fCollisionTime", &vtx.fCollisionTime, "fCollisionTime/F");
    tEvents->Branch("fCollisionTimeRes", &vtx.fCollisionTimeRes, "fCollisionTimeRes/F");
    tEvents->Branch("fCollisionTimeMask", &vtx.fCollisionTimeMask, "fCollisionTimeMask/b");
  }
  PostTree(kEvents);

  // Associate branches for fEventTree
  TTree* tTrigger = CreateTree(kTrigger);
  tTrigger->SetAutoFlush(fNumberOfEventsPerCluster);
  if (fTreeStatus[kTrigger]) {
    tTrigger->Branch("fGlobalBC", &trigger.fGlobalBC, "fGlobalBC/l");
    tTrigger->Branch("fTriggerMask", &trigger.fTriggerMask, "fTriggerMask/l");
  }
  PostTree(kTrigger);

  
  // Associate branches for fTrackTree
  TTree* tTracks = CreateTree(kTracks);
  tTracks->SetAutoFlush(fNumberOfEventsPerCluster);
  if (fTreeStatus[kTracks]) {
    tTracks->Branch("fCollisionsID", &tracks.fCollisionsID, "fCollisionsID/I");
//    tTracks->Branch("fTOFclsIndex", &tracks.fTOFclsIndex, "fTOFclsIndex/I");
//    tTracks->Branch("fNTOFcls", &tracks.fNTOFcls, "fNTOFcls/I");
    tTracks->Branch("fX", &tracks.fX, "fX/F");
    tTracks->Branch("fAlpha", &tracks.fAlpha, "fAlpha/F");
    tTracks->Branch("fY", &tracks.fY, "fY/F");
    tTracks->Branch("fZ", &tracks.fZ, "fZ/F");
    tTracks->Branch("fSnp", &tracks.fSnp, "fSnp/F");
    tTracks->Branch("fTgl", &tracks.fTgl, "fTgl/F");
    tTracks->Branch("fSigned1Pt", &tracks.fSigned1Pt, "fSigned1Pt/F");
    tTracks->Branch("fCYY", &tracks.fCYY, "fCYY/F");
    tTracks->Branch("fCZY", &tracks.fCZY, "fCZY/F");
    tTracks->Branch("fCZZ", &tracks.fCZZ, "fCZZ/F");
    tTracks->Branch("fCSnpY", &tracks.fCSnpY, "fCSnpY/F");
    tTracks->Branch("fCSnpZ", &tracks.fCSnpZ, "fCSnpZ/F");
    tTracks->Branch("fCSnpSnp", &tracks.fCSnpSnp, "fCSnpSnp/F");
    tTracks->Branch("fCTglY", &tracks.fCTglY, "fCTglY/F");
    tTracks->Branch("fCTglZ", &tracks.fCTglZ, "fCTglZ/F");
    tTracks->Branch("fCTglSnp", &tracks.fCTglSnp, "fCTglSnp/F");
    tTracks->Branch("fCTglTgl", &tracks.fCTglTgl, "fCTglTgl/F");
    tTracks->Branch("fC1PtY", &tracks.fC1PtY, "fC1PtY/F");
    tTracks->Branch("fC1PtZ", &tracks.fC1PtZ, "fC1PtZ/F");
    tTracks->Branch("fC1PtSnp", &tracks.fC1PtSnp, "fC1PtSnp/F");
    tTracks->Branch("fC1PtTgl", &tracks.fC1PtTgl, "fC1PtTgl/F");
    tTracks->Branch("fC1Pt21Pt2", &tracks.fC1Pt21Pt2, "fC1Pt21Pt2/F");
    tTracks->Branch("fTPCinnerP", &tracks.fTPCinnerP, "fTPCinnerP/F");
    tTracks->Branch("fFlags", &tracks.fFlags, "fFlags/l");
    tTracks->Branch("fITSClusterMap", &tracks.fITSClusterMap, "fITSClusterMap/b");
    tTracks->Branch("fTPCnclsFindable", &tracks.fTPCnclsFindable, "fTPCnclsFindable/b");
    tTracks->Branch("fTPCnclsFindableMinusFound",&tracks.fTPCnclsFindableMinusFound, "fTPCnclsFindableMinusFound/B");
    tTracks->Branch("fTPCnclsFindableMinusCrossedRows", &tracks.fTPCnclsFindableMinusCrossedRows, "fTPCnclsFindableMinusCrossedRows/B");
    tTracks->Branch("fTPCnclsShared", &tracks.fTPCnclsShared, "fTPCnclsShared/b");
    tTracks->Branch("fTRDntracklets", &tracks.fTRDntracklets, "fTRDntracklets/b");
    tTracks->Branch("fITSchi2Ncl", &tracks.fITSchi2Ncl, "fITSchi2Ncl/F");
    tTracks->Branch("fTPCchi2Ncl", &tracks.fTPCchi2Ncl, "fTPCchi2Ncl/F");
    tTracks->Branch("fTRDchi2", &tracks.fTRDchi2, "fTRDchi2/F");
    tTracks->Branch("fTOFchi2", &tracks.fTOFchi2, "fTOFchi2/F");
    tTracks->Branch("fTPCsignal", &tracks.fTPCsignal, "fTPCsignal/F");
    tTracks->Branch("fTRDsignal", &tracks.fTRDsignal, "fTRDsignal/F");
    tTracks->Branch("fTOFsignal", &tracks.fTOFsignal, "fTOFsignal/F");
    tTracks->Branch("fLength", &tracks.fLength, "fLength/F");
  }
  PostTree(kTracks);

  // Associate branches for Calo
  TTree* tCalo = CreateTree(kCalo);
  tCalo->SetAutoFlush(fNumberOfEventsPerCluster);
  if (fTreeStatus[kCalo]) {
    tCalo->Branch("fCollisionsID", &calo.fCollisionsID, "fCollisionsID/I");
    tCalo->Branch("fCellNumber", &calo.fCellNumber, "fCellNumber/S");
    tCalo->Branch("fAmplitude", &calo.fAmplitude, "fAmplitude/F");
    tCalo->Branch("fTime", &calo.fTime, "fTime/F");
    tCalo->Branch("fCellType", &calo.fCellType, "fCellType/C");
    tCalo->Branch("fType", &calo.fType, "fType/B");
  }
  PostTree(kCalo);

  TTree *tCaloTrigger = CreateTree(kCaloTrigger);
  tCaloTrigger->SetAutoFlush(fNumberOfEventsPerCluster);
  if (fTreeStatus[kCaloTrigger]) {
    tCaloTrigger->Branch("fCollisionsID", &calotrigger.fCollisionsID, "fCollisionsID/I");
    tCaloTrigger->Branch("fFastOrAbsID", &calotrigger.fFastorAbsID, "fFastorAbsID/S");
    tCaloTrigger->Branch("fL0Amplitude", &calotrigger.fL0Amplitude, "fL0Amplitude/F");
    tCaloTrigger->Branch("fL1TimeSum", &calotrigger.fL1TimeSum, "fL1TimeSum/F");
    tCaloTrigger->Branch("fNL0Times", &calotrigger.fNL0Times, "fNL0Times/C");
    tCaloTrigger->Branch("fTriggerBits", &calotrigger.fTriggerBits, "fTriggerBits/I");
    tCaloTrigger->Branch("fType", &calotrigger.fType, "fType/B");
  }
  PostTree(kCaloTrigger);

  // Associuate branches for MUON tracks
  TTree* tMuon = CreateTree(kMuon);
  tMuon->SetAutoFlush(fNumberOfEventsPerCluster);
  if (fTreeStatus[kMuon]) {
    tMuon->Branch("fCollisionsID", &muons.fCollisionsID, "fCollisionsID/I");
//    tMuon->Branch("fClusterIndex", &muons.fClusterIndex, "fClusterIndex/I");
//    tMuon->Branch("fNclusters", &muons.fNclusters, "fNclusters/I");
    tMuon->Branch("fInverseBendingMomentum", &muons.fInverseBendingMomentum, "fInverseBendingMomentum/F");
    tMuon->Branch("fThetaX", &muons.fThetaX, "fThetaX/F");
    tMuon->Branch("fThetaY", &muons.fThetaY, "fThetaY/F");
    tMuon->Branch("fZ", &muons.fZ, "fZ/F");
    tMuon->Branch("fBendingCoor", &muons.fBendingCoor, "fBendingCoor/F");
    tMuon->Branch("fNonBendingCoor", &muons.fNonBendingCoor, "fNonBendingCoor/F");
    tMuon->Branch("fCovariances", muons.fCovariances, "fCovariances[15]/F");
    tMuon->Branch("fChi2", &muons.fChi2, "fChi2/F");
    tMuon->Branch("fChi2MatchTrigger", &muons.fChi2MatchTrigger, "fChi2MatchTrigger/F");
  }
  PostTree(kMuon);

  // Associate branches for MUON tracks
  TTree* tMuonCls = CreateTree(kMuonCls);
  tMuonCls->SetAutoFlush(fNumberOfEventsPerCluster);
  if (fTreeStatus[kMuonCls]) {
    tMuonCls->Branch("fMuonsID",&mucls.fMuonsID,"fMuonsID/I");
    tMuonCls->Branch("fX",&mucls.fX,"fX/F");
    tMuonCls->Branch("fY",&mucls.fY,"fY/F");
    tMuonCls->Branch("fZ",&mucls.fZ,"fZ/F");
    tMuonCls->Branch("fErrX",&mucls.fErrX,"fErrX/F");
    tMuonCls->Branch("fErrY",&mucls.fErrY,"fErrY/F");
    tMuonCls->Branch("fCharge",&mucls.fCharge,"fCharge/F");
    tMuonCls->Branch("fChi2",&mucls.fChi2,"fChi2/F");
  }
  PostTree(kMuonCls);

  // Associuate branches for ZDC
  TTree* tZdc = CreateTree(kZdc);
  tZdc->SetAutoFlush(fNumberOfEventsPerCluster);
  if (fTreeStatus[kZdc]) {
    tZdc->Branch("fCollisionsID", &zdc.fCollisionsID, "fCollisionsID/I");
    tZdc->Branch("fZEM1Energy", &zdc.fZEM1Energy, "fZEM1Energy/F");
    tZdc->Branch("fZEM2Energy", &zdc.fZEM2Energy, "fZEM2Energy/F");
    tZdc->Branch("fZNCTowerEnergy", zdc.fZNCTowerEnergy, "fZNCTowerEnergy[5]/F");
    tZdc->Branch("fZNATowerEnergy", zdc.fZNATowerEnergy, "fZNATowerEnergy[5]/F");
    tZdc->Branch("fZPCTowerEnergy", zdc.fZPCTowerEnergy, "fZPCTowerEnergy[5]/F");
    tZdc->Branch("fZPATowerEnergy", zdc.fZPATowerEnergy, "fZPATowerEnergy[5]/F");
    tZdc->Branch("fZNCTowerEnergyLR", zdc.fZNCTowerEnergyLR, "fZNCTowerEnergyLR[5]/F");
    tZdc->Branch("fZNATowerEnergyLR", zdc.fZNATowerEnergyLR, "fZNATowerEnergyLR[5]/F");
    tZdc->Branch("fZPCTowerEnergyLR", zdc.fZPCTowerEnergyLR, "fZPCTowerEnergyLR[5]/F");
    tZdc->Branch("fZPATowerEnergyLR", zdc.fZPATowerEnergyLR, "fZPATowerEnergyLR[5]/F");
    tZdc->Branch("fZDCTDCCorrected", zdc.fZDCTDCCorrected, "fZDCTDCCorrected[32][4]/F");
    tZdc->Branch("fFired", &zdc.fFired, "fFired/b");
  }
  PostTree(kZdc);

  // Associuate branches for VZERO
  TTree* tVzero = CreateTree(kVzero);
  tVzero->SetAutoFlush(fNumberOfEventsPerCluster);
  if (fTreeStatus[kVzero]) {
    tVzero->Branch("fCollisionsID", &vzero.fCollisionsID, "fCollisionsID/I");
    tVzero->Branch("fAdc", vzero.fAdc, "fAdc[64]/F");
    tVzero->Branch("fTime", vzero.fTime, "fTime[64]/F");
    tVzero->Branch("fWidth", vzero.fWidth, "fWidth[64]/F");
    tVzero->Branch("fBBFlag", &vzero.fBBFlag, "fBBFlag/l");
    tVzero->Branch("fBGFlag", &vzero.fBGFlag, "fBGFlag/l");
  }
  PostTree(kVzero);

  // Associuate branches for V0s
  TTree* tV0s = CreateTree(kV0s);
  tV0s->SetAutoFlush(fNumberOfEventsPerCluster);
  if (fTreeStatus[kV0s]) {
    tV0s->Branch("fPosTrackID", &v0s.fPosTrackID, "fPosTrackID/I");
    tV0s->Branch("fNegTrackID", &v0s.fNegTrackID, "fNegTrackID/I");
  }
  PostTree(kV0s);

  // Associuate branches for cascades
  TTree* tCascades = CreateTree(kCascades);
  tCascades->SetAutoFlush(fNumberOfEventsPerCluster);
  if (fTreeStatus[kCascades]) {
    tCascades->Branch("fV0sID", &cascs.fV0sID, "fV0sID/I");
    tCascades->Branch("fTracksID", &cascs.fTracksID, "fTracksID/I");
  }
  PostTree(kCascades);

#ifdef USE_TOF_CLUST
  // Associate branches for TOF
  TTree* TOF = CreateTree(kTOF);
  TOF->SetAutoFlush(fNumberOfEventsPerCluster);
  if (fTreeStatus[kTOF]) {
    TOF->Branch("fTOFChannel", &tofClusters.fTOFChannel, "fTOFChannel/I");
    TOF->Branch("fTOFncls", &tofClusters.fTOFncls, "fTOFncls/S");
    TOF->Branch("fDx", &tofClusters.fDx, "fDx/F");
    TOF->Branch("fDz", &tofClusters.fDz, "fDz/F");
    TOF->Branch("fToT", &tofClusters.fToT, "fToT/F");
  }
  PostTree(kTOF);
#endif

  if (fTaskMode == kMC) {
    TTree * tMCvtx = CreateTree(kMCvtx);
    tMCvtx->SetAutoFlush(fNumberOfEventsPerCluster);
    if(fTreeStatus[kMCvtx]) {
      tMCvtx->Branch("fGeneratorsID", &mcvtx.fGeneratorsID, "fGeneratorsID/S");
      tMCvtx->Branch("fX", &mcvtx.fX, "fX/F");
      tMCvtx->Branch("fY", &mcvtx.fY, "fY/F");
      tMCvtx->Branch("fZ", &mcvtx.fZ, "fZ/F");
      tMCvtx->Branch("fT", &mcvtx.fT, "fT/F");
      tMCvtx->Branch("fWeight", &mcvtx.fWeight, "fWeight/F");
      tMCvtx->Branch("fNProduced", &mcvtx.fNProduced, "fNProduced/I");
    }
    PostTree(kMCvtx);

    // Associate branches for Kinematics
    TTree* Kinematics = CreateTree(kKinematics);
    Kinematics->SetAutoFlush(fNumberOfEventsPerCluster);
    if (fTreeStatus[kMC]) {
      Kinematics->Branch("fCollisionsID", &mcparticle.fCollisionsID, "fCollisionsID/I");

      Kinematics->Branch("fPdgCode", &mcparticle.fPdgCode, "fPdgCode/I");
      Kinematics->Branch("fStatusCode", &mcparticle.fStatusCode, "fStatusCode/I");
      Kinematics->Branch("fMother", &mcparticle.fMother, "fMother[2]/I");
      Kinematics->Branch("fDaughter", &mcparticle.fDaughter, "fDaughter[2]/I");
      Kinematics->Branch("fWeight", &mcparticle.fWeight, "fWeight/F");
      
      Kinematics->Branch("fPx", &mcparticle.fPx, "fPx/F");
      Kinematics->Branch("fPy", &mcparticle.fPy, "fPy/F");
      Kinematics->Branch("fPz", &mcparticle.fPz, "fPz/F");
      Kinematics->Branch("fE", &mcparticle.fE, "fE/F");
      
      Kinematics->Branch("fVx", &mcparticle.fVx, "fVx/F");
      Kinematics->Branch("fVy", &mcparticle.fVy, "fVy/F");
      Kinematics->Branch("fVz", &mcparticle.fVz, "fVz/F");
      Kinematics->Branch("fVt", &mcparticle.fVt, "fVt/F");
    }
    PostTree(kKinematics);

    // Range for the MC labels of each reconstructed track
    TTree* tRange = CreateTree(kRange);
    tRange->SetAutoFlush(fNumberOfEventsPerCluster);
    if (fTreeStatus[kRange]) {
      tRange->Branch("fRange", &range.fRange, "fRange/i");
      FillTree(kRange); // Put the begin of the first range to 0
    }
    PostTree(kRange);
    
    // MC labels of each reconstructed track
    TTree* tLabels = CreateTree(kLabels);
    tLabels->SetAutoFlush(fNumberOfEventsPerCluster);
    if (fTreeStatus[kLabels]) {
      tLabels->Branch("fLabel", &labels.fLabel, "fLabel/I");
    }
    PostTree(kLabels);
  }


  Prune(); //Removing all unwanted branches (if any)
}

void AliAnalysisTaskAO2Dconverter::Prune()
{
  if (fPruneList.IsNull() || fPruneList.IsWhitespace())
    return;
  TObjArray* arr = fPruneList.Tokenize(" ");
  for (Int_t i = 0; i < arr->GetEntries(); i++) {
    Bool_t found = kFALSE;
    for (Int_t j = 0; j < kTrees; j++) {
      TObjArray* branches = fTree[j]->GetListOfBranches();
      for (Int_t k = 0; k < branches->GetEntries(); k++) {
        TString bname = branches->At(k)->GetName();
        if (!bname.EqualTo(arr->At(i)->GetName()))
          continue;
        fTree[j]->SetBranchStatus(bname, 0);
        found = kTRUE;
      }
    }
    if (!found)
      AliFatal(Form("Did not find Branch %s", arr->At(i)->GetName()));
  }
  fPruneList = "";
}

void AliAnalysisTaskAO2Dconverter::UserExec(Option_t *)
{
  // Initialisation

  fESD = dynamic_cast<AliESDEvent *>(InputEvent());
  if (!fESD) {
    ::Fatal("AliAnalysisTaskAO2Dconverter::UserExec", "Something is wrong with the event handler");
  }

  // We can use event cuts to avoid cases where we have zero reconstructed tracks
  if (fUseEventCuts && !fEventCuts.AcceptEvent(fESD)) {
    return;
  }

  // Configuration of the PID response
  AliPIDResponse* PIDResponse = (AliPIDResponse*)((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->GetPIDResponse();
  PIDResponse->SetTOFResponse(fESD, AliPIDResponse::kBest_T0);
  AliTOFPIDResponse & TOFResponse = PIDResponse->GetTOFResponse();

  // Configuration of the MC event (if needed)
  AliMCEvent* MCEvt = nullptr;
  if (fTaskMode == kMC) {
    AliMCEventHandler* eventHandler = dynamic_cast<AliMCEventHandler*>(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()); //Get the MC handler

    if (!eventHandler) //Check on the MC handler
      AliFatal("Could not retrieve MC event handler");
    MCEvt = eventHandler->MCEvent(); //Get the MC Event

    if (!MCEvt) // Check on the MC Event
      AliFatal("Could not retrieve MC event");
    PIDResponse->SetCurrentMCEvent(MCEvt); //Set The PID response on the current MC event
  }

  // Selection of events with at least two contributors (GMI)
  // Can this be done using the physics selection? (PH)

  const AliESDVertex * pvtx = fESD->GetPrimaryVertex();
  if (!pvtx) {
    ::Fatal("AliAnalysisTaskAO2Dconverter::UserExec", "Vertex not defined");
  }
  TString title=pvtx->GetTitle();
  if(pvtx->IsFromVertexer3D() || pvtx->IsFromVertexerZ()) return;
  if(pvtx->GetNContributors()<2) return;

  // Get access to the current event number
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  Int_t eventID = fEventCount++;

  //---------------------------------------------------------------------------
  // Collision data

  // Adjust start indices for this event in all trees by adding the number of entries of the previous event
  for (auto i = 0; i < kTrees; ++i)
     vtx.fStart[i] += vtx.fNentries[i];

  vtx.fNentries[kEvents] = 1;  // one entry per vertex
  vtx.fRunNumber = fESD->GetRunNumber();
  ULong64_t evtid = GetEventIdAsLong(fESD->GetHeader());
  if(!evtid){
    evtid = (ULong64_t(fESD->GetTimeStamp())<<32) + ULong64_t((fESD->GetNumberOfTPCClusters()<<5)|(fESD->GetNumberOfTPCTracks()));
  }
  vtx.fGlobalBC = evtid;
  vtx.fX = pvtx->GetX();
  vtx.fY = pvtx->GetY();
  vtx.fZ = pvtx->GetZ();

  Double_t covmatrix[6];
  pvtx->GetCovMatrix(covmatrix);

  vtx.fCovXX = covmatrix[0];
  vtx.fCovXY = covmatrix[1];
  vtx.fCovXZ = covmatrix[2];
  vtx.fCovYY = covmatrix[3];
  vtx.fCovYZ = covmatrix[4];
  vtx.fCovZZ = covmatrix[5];

  vtx.fChi2 = pvtx->GetChi2();
  vtx.fN = (pvtx->GetNDF()+3)/2;

  Float_t eventTime[10];
  Float_t eventTimeRes[10];
  Double_t eventTimeWeight[10];
  
  for (Int_t i = 0; i < TOFResponse.GetNmomBins(); i++) {
    if (i >= 10)
      AliFatal("Index is too high!");
    Float_t mom = (TOFResponse.GetMinMom(i) + TOFResponse.GetMaxMom(i)) / 2.f;
    eventTime[i] = TOFResponse.GetStartTime(mom);
    eventTimeRes[i] = TOFResponse.GetStartTimeRes(mom);
    eventTimeWeight[i] = 1./(eventTimeRes[i]*eventTimeRes[i]);

    //PH The part below is just a place holder
    if (TOFResponse.GetStartTimeMask(mom) & 0x1)
      SETBIT(vtx.fCollisionTimeMask, 0);
    else
      CLRBIT(vtx.fCollisionTimeMask, 0);
    //
    if (TOFResponse.GetStartTimeMask(mom) & 0x2)
      SETBIT(vtx.fCollisionTimeMask, 1);
    else
      CLRBIT(vtx.fCollisionTimeMask, 1);
    //
    if (TOFResponse.GetStartTimeMask(mom) & 0x3)
      SETBIT(vtx.fCollisionTimeMask, 2);
    else
      CLRBIT(vtx.fCollisionTimeMask, 2);
  }

  // Recalculate unique event time and its resolution
  vtx.fCollisionTime = TMath::Mean(10,eventTime,eventTimeWeight); // Weighted mean of times per momentum interval
  vtx.fCollisionTimeRes = TMath::Sqrt(9./10.)*TMath::Mean(10,eventTimeRes); // PH bad approximation

  //---------------------------------------------------------------------------
  // Trigger data
  
  trigger.fGlobalBC = evtid;
  trigger.fTriggerMask = fESD->GetTriggerMask();
  FillTree(kTrigger);
  
  //---------------------------------------------------------------------------
  // Track data

  Int_t ntrk = fESD->GetNumberOfTracks();

  Int_t ntrk_filled = 0;     // total number of tracks filled per event
  Int_t ntofcls_filled = 0;  // total number of TOF clusters filled per event

  for (Int_t itrk = 0; itrk < ntrk; itrk++)
  {
    AliESDtrack *track = fESD->GetTrack(itrk);
//    if (!fTrackFilter.IsSelected(track))
//      continue;
    tracks.fCollisionsID = eventID;

    tracks.fX = track->GetX();
    tracks.fAlpha = track->GetAlpha();

    tracks.fY = track->GetY();
    tracks.fZ = track->GetZ();
    tracks.fSnp = track->GetSnp();
    tracks.fTgl = track->GetTgl();
    tracks.fSigned1Pt = track->GetSigned1Pt();

    tracks.fCYY = track->GetSigmaY2();
    tracks.fCZY = track->GetSigmaZY();
    tracks.fCZZ = track->GetSigmaZ2();
    tracks.fCSnpY = track->GetSigmaSnpY();
    tracks.fCSnpZ = track->GetSigmaSnpZ();
    tracks.fCSnpSnp = track->GetSigmaSnp2();
    tracks.fCTglY = track->GetSigmaTglY();
    tracks.fCTglZ = track->GetSigmaTglZ();
    tracks.fCTglSnp = track->GetSigmaTglSnp();
    tracks.fCTglTgl = track->GetSigmaTgl2();
    tracks.fC1PtY = track->GetSigma1PtY();
    tracks.fC1PtZ = track->GetSigma1PtZ();
    tracks.fC1PtSnp = track->GetSigma1PtSnp();
    tracks.fC1PtTgl = track->GetSigma1PtTgl();
    tracks.fC1Pt21Pt2 = track->GetSigma1Pt2();

    const AliExternalTrackParam *intp = track->GetTPCInnerParam();
    tracks.fTPCinnerP = (intp ? intp->GetP() : 0); // Set the momentum to 0 if the track did not reach TPC

    tracks.fFlags = track->GetStatus();

    tracks.fITSClusterMap = track->GetITSClusterMap();
    tracks.fTPCnclsFindable = track->GetTPCNclsF();
    tracks.fTPCnclsFindableMinusFound = tracks.fTPCnclsFindable - track->GetTPCNcls();
    tracks.fTPCnclsFindableMinusCrossedRows = tracks.fTPCnclsFindable - track->GetTPCCrossedRows();
    tracks.fTPCnclsShared = (track->GetTPCSharedMap()).CountBits();
    tracks.fTRDntracklets = track->GetTRDntracklets();

    tracks.fITSchi2Ncl = (track->GetITSNcls() ? track->GetITSchi2() / track->GetITSNcls() : 0);
    tracks.fTPCchi2Ncl = (track->GetTPCNcls() ? track->GetTPCchi2() / track->GetTPCNcls() : 0);
    tracks.fTRDchi2 = track->GetTRDchi2();
    tracks.fTOFchi2 = track->GetTOFchi2();

    tracks.fTPCsignal = track->GetTPCsignal();
    tracks.fTRDsignal = track->GetTRDsignal();
    tracks.fTOFsignal = track->GetTOFsignal();
    tracks.fLength = track->GetIntegratedLength();

    if (fTaskMode == kMC) {
      // Separate tables (trees) for the MC labels
      // Right now we have only one label, the data model is adapted to many labels
      // We expect a loop on the labels in Run3 MC
      // for (ilabel=0; ilabel<nlabels; ++ilabel) {
      Int_t alabel = track->GetLabel();
      labels.fLabel = TMath::Sign(TMath::Abs(alabel) + fOffsetLabel, alabel); // keep the sign of the label
      FillTree(kLabels);
      range.fRange++;
      // } // End of loop on labels
      FillTree(kRange);
  }

#ifdef USE_TOF_CLUST
    tofClusters.fTOFncls = track->GetNTOFclusters();

    if (fTreeStatus[kTOF] && tofClusters.fTOFncls > 0) {
      Int_t* TOFclsIndex = track->GetTOFclusterArray(); //Index of the matchable cluster (there are fNTOFClusters of them)
      for (Int_t icls = 0; icls < tofClusters.fTOFncls; icls++) {
        AliESDTOFCluster* TOFcls = (AliESDTOFCluster*)fESD->GetESDTOFClusters()->At(TOFclsIndex[icls]);
        tofClusters.fToT = TOFcls->GetTOFsignalToT(0);
        tofClusters.fTOFChannel = TOFcls->GetTOFchannel();
        for (Int_t mtchbl = 0; mtchbl < TOFcls->GetNMatchableTracks(); mtchbl++) {
          if (TOFcls->GetTrackIndex(mtchbl) != track->GetID())
            continue;
          tofClusters.fDx = TOFcls->GetDx(mtchbl);
          tofClusters.fDz = TOFcls->GetDz(mtchbl);
          tofClusters.fLengthRatio = tracks.fLength > 0 ? TOFcls->GetLength(mtchbl) / tracks.fLength : -1;
          break;
        }
        FillTree(kTOF);
	if (fTreeStatus[kTOF]) ntofcls_filled++;
      }
    }
#endif

    // In case we need connection to clusters, activate next lines
    // tracks.fTOFclsIndex += tracks.fNTOFcls;
    // tracks.fNTOFcls = ntofcls_filled;
    FillTree(kTracks);
    if (fTreeStatus[kTracks]) ntrk_filled++;
  } // end loop on tracks
  vtx.fNentries[kTOF]    = ntofcls_filled;
  vtx.fNentries[kTracks] = ntrk_filled;

  //---------------------------------------------------------------------------
  // Calorimeter data

  AliESDCaloCells *cells = fESD->GetEMCALCells();
  Short_t nCells = cells->GetNumberOfCells();
  Int_t ncalocells_filled = 0; // total number of calo cells filled per event
  for (Short_t ice = 0; ice < nCells; ++ice)
  {
    Short_t cellNumber;
    Double_t amplitude;
    Double_t time;
    Int_t mclabel;
    Double_t efrac;

    calo.fCollisionsID = eventID;
    
    cells->GetCell(ice, cellNumber, amplitude, time, mclabel, efrac);
    calo.fCellNumber = cellNumber;
    calo.fAmplitude = amplitude;
    calo.fTime = time;
    calo.fType = cells->GetType(); // common for all cells
    calo.fCellType = cells->GetHighGain(ice) ? 0. : 1.; 
    FillTree(kCalo);
    if (fTreeStatus[kCalo]) ncalocells_filled++;
  } // end loop on calo cells
  vtx.fNentries[kCalo] = ncalocells_filled;

  AliEMCALGeometry *geo = AliEMCALGeometry::GetInstanceFromRunNumber(fESD->GetRunNumber()); // Needed for EMCAL trigger mapping
  AliESDCaloTrigger *calotriggers = fESD->GetCaloTrigger("EMCAL");
  calotriggers->Reset();
  Int_t ncalotriggers_filled = 0; // total number of EMCAL triggers filled per event
  while(calotriggers->Next()){
    calotrigger.fCollisionsID = eventID;
    int col, row, fastorID;
    calotriggers->GetPosition(col, row);
    geo->GetTriggerMapping()->GetAbsFastORIndexFromPositionInEMCAL(col, row, fastorID);
    calotrigger.fFastorAbsID = fastorID;
    calotriggers->GetAmplitude(calotrigger.fL0Amplitude);
    calotriggers->GetTime(calotrigger.fL0Time);
    calotriggers->GetTriggerBits(calotrigger.fTriggerBits);
    Int_t nL0times;
    calotriggers->GetNL0Times(nL0times);
    calotrigger.fNL0Times = nL0times;
    calotriggers->GetL1TimeSum(calotrigger.fTriggerBits);
    calotrigger.fType = 1;
    FillTree(kCaloTrigger);
    if (fTreeStatus[kCaloTrigger]) ncalotriggers_filled++;
  }
  vtx.fNentries[kCaloTrigger] = ncalotriggers_filled;

  cells = fESD->GetPHOSCells();
  nCells = cells->GetNumberOfCells();
  Int_t nphoscells_filled = 0;
  for (Short_t icp = 0; icp < nCells; ++icp)
  {
    Short_t cellNumber;
    Double_t amplitude;
    Double_t time;
    Int_t mclabel;
    Double_t efrac;

    calo.fCollisionsID = eventID;
    
    cells->GetCell(icp, cellNumber, amplitude, time, mclabel, efrac);
    calo.fCellNumber = cellNumber;
    calo.fAmplitude = amplitude;
    calo.fTime = time;
    calo.fCellType = cells->GetHighGain(icp) ? 0. : 1.;     /// @TODO cell type value to be confirmed by PHOS experts
    calo.fType = cells->GetType(); // common for all cells

    FillTree(kCalo);
    if (fTreeStatus[kCalo]) nphoscells_filled++;
  } // end loop on PHOS cells
  vtx.fNentries[kCalo] = nphoscells_filled;

  //---------------------------------------------------------------------------
  // Muon tracks
  muons.fCollisionsID  = eventID;
  
  Int_t nmu = fESD->GetNumberOfMuonTracks();
  Int_t nmu_filled = 0;    // total number of muons filled per event
  Int_t nmucl_filled = 0;  // total number of clusters filled per event
  for (Int_t imu=0; imu<nmu; ++imu) {
    AliESDMuonTrack* mutrk = fESD->GetMuonTrack(imu);

    muons.fInverseBendingMomentum = mutrk->GetInverseBendingMomentum();
    muons.fThetaX = mutrk->GetThetaX();
    muons.fThetaY = mutrk->GetThetaY();
    muons.fZ = mutrk->GetZ();
    muons.fBendingCoor = mutrk->GetBendingCoor();
    muons.fNonBendingCoor = mutrk->GetNonBendingCoor();

    TMatrixD cov;
    mutrk->GetCovariances(cov);
    for (Int_t i = 0; i < 5; i++)
      for (Int_t j = 0; j <= i; j++)
	muons.fCovariances[i*(i+1)/2 + j] = cov(i,j);

    muons.fChi2 = mutrk->GetChi2();
    muons.fChi2MatchTrigger = mutrk->GetChi2MatchTrigger();

    // Now MUON clusters for the current track
    Int_t muTrackID = fOffsetMuTrackID + imu;
    Int_t nmucl = mutrk->GetNClusters();
    for (Int_t imucl=0; imucl<nmucl; ++imucl){
      AliESDMuonCluster *muCluster = fESD->FindMuonCluster(mutrk->GetClusterId(imucl));
      mucls.fMuonsID = muTrackID;
      mucls.fX = muCluster->GetX();
      mucls.fY = muCluster->GetY();
      mucls.fZ = muCluster->GetZ();
      mucls.fErrX = muCluster->GetErrX();
      mucls.fErrY = muCluster->GetErrY();
      mucls.fCharge = muCluster->GetCharge();
      mucls.fChi2   = muCluster->GetChi2();
      FillTree(kMuonCls);
      if (fTreeStatus[kMuonCls]) nmucl_filled++;
    } // End loop on muon clusters for the current muon track

    // In case we need connection to clusters, activate next lines
    // muons.fClusterIndex += muons.fNclusters;
    // muons.fNclusters = nmucl_filled;

    FillTree(kMuon);
    if (fTreeStatus[kMuon]) nmu_filled++;
  } // End loop on muon tracks
  vtx.fNentries[kMuon] = nmu_filled;
  vtx.fNentries[kMuonCls] = nmucl_filled;

  //---------------------------------------------------------------------------
  // ZDC
  AliESDZDC* esdzdc  =    fESD->GetESDZDC();
  zdc.fCollisionsID = eventID;
  // ZEM
  zdc.fZEM1Energy = esdzdc->GetZEM1Energy();
  zdc.fZEM2Energy = esdzdc->GetZEM2Energy();
  // ZDC (P,N) towers
  for (Int_t ich=0; ich<5; ++ich) {
    zdc.fZNCTowerEnergy[ich] = esdzdc->GetZNCTowerEnergy()[ich];
    zdc.fZNATowerEnergy[ich] = esdzdc->GetZNATowerEnergy()[ich];
    zdc.fZPCTowerEnergy[ich] = esdzdc->GetZPCTowerEnergy()[ich];
    zdc.fZPATowerEnergy[ich] = esdzdc->GetZPATowerEnergy()[ich];
    
    zdc.fZNCTowerEnergyLR[ich] = esdzdc->GetZNCTowerEnergyLR()[ich];
    zdc.fZNATowerEnergyLR[ich] = esdzdc->GetZNATowerEnergyLR()[ich];
    zdc.fZPCTowerEnergyLR[ich] = esdzdc->GetZPCTowerEnergyLR()[ich];
    zdc.fZPATowerEnergyLR[ich] = esdzdc->GetZPATowerEnergyLR()[ich];
  }
  // ZDC TDC
  for (Int_t ii=0; ii< 32; ++ii)
    for (Int_t jj=0; jj<4; ++jj)
      zdc.fZDCTDCCorrected[ii][jj] = esdzdc->GetZDCTDCCorrected(ii,jj);
  // ZDC flags
  zdc.fFired = 0x0;                  // Bits: 0 - ZNA, 1 - ZNC, 2 - ZPA, 3 - ZPC, 4 - ZEM1, 5 - ZEM2
  if (esdzdc->IsZNAhit()) zdc.fFired |= (0x1);
  if (esdzdc->IsZNChit()) zdc.fFired |= (0x1 << 1);
  if (esdzdc->IsZPAhit()) zdc.fFired |= (0x1 << 2);
  if (esdzdc->IsZPChit()) zdc.fFired |= (0x1 << 3);
  if (esdzdc->IsZEM1hit()) zdc.fFired |= (0x1 << 4);
  if (esdzdc->IsZEM2hit()) zdc.fFired |= (0x1 << 5);
  FillTree(kZdc);
  if (fTreeStatus[kZdc]) vtx.fNentries[kZdc] = 1;

  //---------------------------------------------------------------------------
  // VZERO
  AliESDVZERO * vz = fESD->GetVZEROData();
  vzero.fCollisionsID  = eventID;
  for (Int_t ich=0; ich<64; ++ich) {
    vzero.fAdc[ich] = vz->GetAdc(ich);
    vzero.fTime[ich] = vz->GetTime(ich);
    vzero.fWidth[ich] = vz->GetWidth(ich);
    vzero.fBBFlag = 0u;
    vzero.fBGFlag = 0u;
    ULong64_t mask = 1u;
    for (Int_t i=0; i<64; ++i) {
      if (vz->GetBBFlag(i))
	vzero.fBBFlag |= (mask << i);
      if (vz->GetBGFlag(i))
	vzero.fBGFlag |= (mask << i);
    }
  }
  FillTree(kVzero);
  if (fTreeStatus[kVzero]) vtx.fNentries[kVzero] = 1;

  //---------------------------------------------------------------------------
  // V0s (Lambda and KS)
  Int_t nv0 = fESD->GetNumberOfV0s();
  Int_t nv0_filled = 0; // total number of v0's filled per event
  for (Int_t iv0=0; iv0<nv0; ++iv0) {
    AliESDv0 * v0 = fESD->GetV0(iv0);
    // select only "offline" V0s, skip the "on-the-fly" ones
    if (v0 && !v0->GetOnFlyStatus()) {
      Int_t pidx = v0->GetPindex();
      Int_t nidx = v0->GetNindex();
      v0s.fPosTrackID = TMath::Sign(TMath::Abs(pidx) + fOffsetTrackID, pidx); // Positive track ID
      v0s.fNegTrackID = TMath::Sign(TMath::Abs(nidx) + fOffsetTrackID, nidx); // Negative track ID
      FillTree(kV0s);
      if (fTreeStatus[kV0s]) nv0_filled++;
    }
  } // End loop on V0s
  vtx.fNentries[kV0s] = nv0_filled;

  //---------------------------------------------------------------------------
  // Cascades
  // If we do not have V0s, we do not have cascades
  Int_t ncascades_filled = 0; // total number of cascades filled per event
  if (nv0>0) {
    // Combine the track indexes of V0 daughters in unique identifier
    ULong64_t * packedPosNeg = new ULong64_t[nv0];
    ULong64_t * sortedPosNeg = new ULong64_t[nv0];
    Int_t * sortIdx = new Int_t[nv0];

    for (Int_t iv0=0; iv0<nv0; ++iv0) {
      AliESDv0 * v0 = fESD->GetV0(iv0);
      packedPosNeg[iv0] = (v0->GetPindex() << 31) | (v0->GetNindex());
    }
    TMath::Sort(nv0,packedPosNeg,sortIdx,kFALSE);
    for (Int_t iv0=0; iv0<nv0; ++iv0) {
      sortedPosNeg[iv0] = packedPosNeg[sortIdx[iv0]];
    }
  
    Int_t ncas = fESD->GetNumberOfCascades();
    for (Int_t icas=0; icas<ncas; ++icas) {
      AliESDcascade *cas = fESD->GetCascade(icas);
      // Select only cascades containing "offline" V0s
      if (cas && !cas->GetOnFlyStatus()) {
	// Find the identifier of the V0 using the indexes of its daughters
	ULong64_t currV0 = (cas->GetPindex() << 31) | cas->GetNindex();
	// Use binary search in the sorted array
	Int_t v0idx = TMath::BinarySearch(nv0, sortedPosNeg, currV0);
	// Check if the match is exact
	if (sortedPosNeg[v0idx] == currV0) {
	  cascs.fV0sID = sortIdx[v0idx] + fOffsetV0ID;
	  cascs.fTracksID = cas->GetBindex() + fOffsetTrackID;
	  FillTree(kCascades);
	  if (fTreeStatus[kCascades]) ncascades_filled++;
	}
      }
    } // End loop on cascades

    delete [] packedPosNeg;
    delete [] sortedPosNeg;
    delete [] sortIdx;
  } // End if V0s
  vtx.fNentries[kCascades] = ncascades_filled;
  
  //---------------------------------------------------------------------------
  // MC data (to be modified)

  Int_t nkine_filled = 0; // Number of kine tracks filled
  if (MCEvt) {
    // Kinematics
    TParticle* particle = nullptr;
    Int_t nMCtracks = MCEvt->GetNumberOfTracks();
    for (Int_t i = 0; i < nMCtracks; ++i) { //loop on primary MC tracks Before Event Selection
      AliVParticle* vpt = MCEvt->GetTrack(i);
      particle = vpt->Particle();

      mcparticle.fCollisionsID = eventID;
      
      //Get the kinematic values of the particles
      mcparticle.fPdgCode = particle->GetPdgCode();
      mcparticle.fStatusCode = particle->GetStatusCode();
      mcparticle.fMother[0] = vpt->GetMother();
      if (mcparticle.fMother[0] > -1) mcparticle.fMother[0]+=fOffsetLabel;
      mcparticle.fMother[1] = vpt->GetMother();
      if (mcparticle.fMother[1] > -1) mcparticle.fMother[1]+=fOffsetLabel;
      mcparticle.fDaughter[0] = particle->GetFirstDaughter();
      if (mcparticle.fDaughter[0] > -1) mcparticle.fDaughter[0]+=fOffsetLabel;
      mcparticle.fDaughter[1] = particle->GetLastDaughter();
      if (mcparticle.fDaughter[1] > -1) mcparticle.fDaughter[1]+=fOffsetLabel;
      mcparticle.fWeight = particle->GetWeight();

      mcparticle.fPx = particle->Px();
      mcparticle.fPy = particle->Py();
      mcparticle.fPz = particle->Pz();
      mcparticle.fE  = particle->Energy();

      mcparticle.fVx = particle->Vx();
      mcparticle.fVy = particle->Vy();
      mcparticle.fVz = particle->Vz();
      mcparticle.fVt = particle->T();

      FillTree(kKinematics);
      if (fTreeStatus[kKinematics]) nkine_filled++;
    }
    fOffsetLabel += nMCtracks; // Offset for the labels of the next event
  }
  vtx.fNentries[kKinematics] = nkine_filled;

  if (MCEvt) {
    // MC vertex
    const AliVVertex* MCvtx = MCEvt->GetPrimaryVertex();
    if (!MCvtx) //Check on the MC vertex
      AliFatal("Could not retrieve MC vertex");
    mcvtx.fX = MCvtx->GetX();
    mcvtx.fY = MCvtx->GetY();
    mcvtx.fZ = MCvtx->GetZ();

    AliGenEventHeader* mcGenH = MCEvt->GenEventHeader();
    mcvtx.fT = mcGenH->InteractionTime();
    mcvtx.fWeight = mcGenH->EventWeight();
    mcvtx.fNProduced = mcGenH->NProduced();

    mcvtx.fGeneratorsID = 0;
    for (Int_t gen = 0; gen < kGenerators; gen++) {
      if (mcGenH->InheritsFrom(Generator[gen]))
        SETBIT(mcvtx.fGeneratorsID, gen);
      else
        CLRBIT(mcvtx.fGeneratorsID, gen);
    }
    if (mcGenH->InheritsFrom(Generator[kAliGenCocktailEventHeader])) {
      TList* headers = ((AliGenCocktailEventHeader*)mcGenH)->GetHeaders();
      for (Int_t cocktail = 0; cocktail < headers->GetEntries(); headers++) {
        for (Int_t gen = 0; gen < kGenerators; gen++) {
          if (mcGenH->InheritsFrom(Generator[gen]))
            SETBIT(mcvtx.fGeneratorsID, gen);
        }
      }
    }
    vtx.fNentries[kMCvtx] = 1;
  } else {
    vtx.fNentries[kMCvtx] = 0;
  }
  // Filling the tree of vertices has to be done last because it contains the
  // index data for the other trees
  FillTree(kMCvtx);


  // We can fill now the vertex + indexing data
  FillTree(kEvents);

  //---------------------------------------------------------------------------
  //Posting data
  for (Int_t i = 0; i < kTrees; i++)
    PostTree((TreeIndex)i);

  //---------------------------------------------------------------------------
  // Update the offsets at the end of each collision    
  fOffsetTrackID += ntrk;
  fOffsetMuTrackID += nmu;
  fOffsetV0ID += nv0;
}

void AliAnalysisTaskAO2Dconverter::Terminate(Option_t *)
{
  // terminate
  // called at the END of the analysis (when all events are processed)
}

AliAnalysisTaskAO2Dconverter *AliAnalysisTaskAO2Dconverter::AddTask(TString suffix)
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    return nullptr;
  }
  // get the input event handler, again via a static method.
  // this handler is part of the managing system and feeds events
  // to your task
  if (!mgr->GetInputEventHandler())
  {
    return nullptr;
  }
  // by default, a file is open for writing. here, we get the filename
  TString fileName = "AO2D.root";
  if (!suffix.IsNull())
    fileName += ":" + suffix; // create a subfolder in the file
  // now we create an instance of your task
  AliAnalysisTaskAO2Dconverter *task = new AliAnalysisTaskAO2Dconverter((TString("AO2D") + suffix).Data());
  if (!task)
    return nullptr;
  // add your task to the manager
  mgr->AddTask(task);
  // your task needs input: here we connect the manager to your task
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  // same for the output
  for (Int_t i = 0; i < kTrees; i++)
    mgr->ConnectOutput(task, 1 + i, mgr->CreateContainer(TreeName[i], TTree::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
  // in the end, this macro returns a pointer to your task. this will be convenient later on
  // when you will run your analysis in an analysis train on grid
  return task;
}
