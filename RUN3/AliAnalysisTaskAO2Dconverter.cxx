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

#include "AliMathBase.h"

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
    , collision()
    , eventextra()
    , bc()
    , tracks()
    , mccollision()
    , mctracklabel()
    , mccalolabel()
    , mccollisionlabel()
    , mcparticle()
#ifdef USE_TOF_CLUST
    , tofClusters()
#endif
    , calo()
    , calotrigger()
    , muons()
    , mucls()
    , zdc()
    , vzero()
    , fdd()
    , v0s()
    , cascs()
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

const TString AliAnalysisTaskAO2Dconverter::TreeName[kTrees] = { "O2collision", "DbgEventExtra", "O2track", "O2calo",  "O2calotrigger", "O2muon", "O2muoncluster", "O2zdc", "Run2v0", "O2fdd", "O2v0", "O2cascade", "O2tof", "O2mcparticle", "O2mccollision", "O2mctracklabel", "O2mccalolabel", "O2mccollisionlabel", "O2bc" };

const TString AliAnalysisTaskAO2Dconverter::TreeTitle[kTrees] = { "Collision tree", "Collision extra", "Barrel tracks", "Calorimeter cells", "Calorimeter triggers", "MUON tracks", "MUON clusters", "ZDC", "Run2 V0", "FDD", "V0s", "Cascades", "TOF hits", "Kinematics", "MC collisions", "MC track labels", "MC calo labels", "MC collision labels", "BC info" };

const TClass* AliAnalysisTaskAO2Dconverter::Generator[kGenerators] = { AliGenEventHeader::Class(), AliGenCocktailEventHeader::Class(), AliGenDPMjetEventHeader::Class(), AliGenEpos3EventHeader::Class(), AliGenEposEventHeader::Class(), AliGenEventHeaderTunedPbPb::Class(), AliGenGeVSimEventHeader::Class(), AliGenHepMCEventHeader::Class(), AliGenHerwigEventHeader::Class(), AliGenHijingEventHeader::Class(), AliGenPythiaEventHeader::Class(), AliGenToyEventHeader::Class() };

TTree* AliAnalysisTaskAO2Dconverter::CreateTree(TreeIndex t)
{
  fTree[t] = new TTree(TreeName[t], TreeTitle[t]);
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
    DisableTree(kMcParticle);
    DisableTree(kMcCollision);
    DisableTree(kMcTrackLabel);
    DisableTree(kMcCaloLabel);
    DisableTree(kMcCollisionLabel);
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
    tEvents->Branch("fBCsID", &collision.fBCsID, "fBCsID/I");
    tEvents->Branch("fPosX", &collision.fPosX, "fPosX/F");
    tEvents->Branch("fPosY", &collision.fPosY, "fPosY/F");
    tEvents->Branch("fPosZ", &collision.fPosZ, "fPosZ/F");
    tEvents->Branch("fCovXX", &collision.fCovXX, "fCovXX/F");
    tEvents->Branch("fCovXY", &collision.fCovXY, "fCovXY/F");
    tEvents->Branch("fCovXZ", &collision.fCovXZ, "fCovXZ/F");
    tEvents->Branch("fCovYY", &collision.fCovYY, "fCovYY/F");
    tEvents->Branch("fCovYZ", &collision.fCovYZ, "fCovYZ/F");
    tEvents->Branch("fCovZZ", &collision.fCovZZ, "fCovZZ/F");
    tEvents->Branch("fChi2", &collision.fChi2, "fChi2/F");
    tEvents->Branch("fNumContrib", &collision.fN, "fNumContrib/i");
    tEvents->Branch("fCollisionTime", &collision.fCollisionTime, "fCollisionTime/F");
    tEvents->Branch("fCollisionTimeRes", &collision.fCollisionTimeRes, "fCollisionTimeRes/F");
    tEvents->Branch("fCollisionTimeMask", &collision.fCollisionTimeMask, "fCollisionTimeMask/b");
  }
  PostTree(kEvents);
  
  // Extra information for debugging for event table
  TTree* tEventsExtra = CreateTree(kEventsExtra);
  tEventsExtra->SetAutoFlush(fNumberOfEventsPerCluster);
  if (fTreeStatus[kEventsExtra]) {
    TString sstart = TString::Format("fStart[%d]/I", kTrees);
    TString sentries = TString::Format("fNentries[%d]/I", kTrees);
    tEventsExtra->Branch("fStart", eventextra.fStart, sstart.Data());
    tEventsExtra->Branch("fNentries", eventextra.fNentries, sentries.Data());
  }
  PostTree(kEventsExtra);

  // Associate branches for fEventTree
  TTree* tBC = CreateTree(kBC);
  tBC->SetAutoFlush(fNumberOfEventsPerCluster);
  if (fTreeStatus[kBC]) {
    tBC->Branch("fRunNumber", &bc.fRunNumber, "fRunNumber/I");
    tBC->Branch("fGlobalBC", &bc.fGlobalBC, "fGlobalBC/l");
    tBC->Branch("fTriggerMask", &bc.fTriggerMask, "fTriggerMask/l");
  }
  PostTree(kBC);
  
  // Associate branches for fTrackTree
  TTree* tTracks = CreateTree(kTracks);
  tTracks->SetAutoFlush(fNumberOfEventsPerCluster);
  if (fTreeStatus[kTracks]) {
    tTracks->Branch("fCollisionsID", &tracks.fCollisionsID, "fCollisionsID/I");
    tTracks->Branch("fTrackType", &tracks.fTrackType, "fTrackType/b");
    //    tTracks->Branch("fTOFclsIndex", &tracks.fTOFclsIndex, "fTOFclsIndex/I");
    //    tTracks->Branch("fNTOFcls", &tracks.fNTOFcls, "fNTOFcls/I");
    tTracks->Branch("fX", &tracks.fX, "fX/F");
    tTracks->Branch("fAlpha", &tracks.fAlpha, "fAlpha/F");
    tTracks->Branch("fY", &tracks.fY, "fY/F");
    tTracks->Branch("fZ", &tracks.fZ, "fZ/F");
    tTracks->Branch("fSnp", &tracks.fSnp, "fSnp/F");
    tTracks->Branch("fTgl", &tracks.fTgl, "fTgl/F");
    tTracks->Branch("fSigned1Pt", &tracks.fSigned1Pt, "fSigned1Pt/F");
    // Modified covariance matrix
    tTracks->Branch("fSigmaY", &tracks.fSigmaY, "fSigmaY/F");
    tTracks->Branch("fSigmaZ", &tracks.fSigmaZ, "fSigmaZ/F");
    tTracks->Branch("fSigmaSnp", &tracks.fSigmaSnp, "fSigmaSnp/F");
    tTracks->Branch("fSigmaTgl", &tracks.fSigmaTgl, "fSigmaTgl/F");
    tTracks->Branch("fSigma1Pt", &tracks.fSigma1Pt, "fSigma1Pt/F");
    tTracks->Branch("fRhoZY", &tracks.fRhoZY, "fRhoZY/B");
    tTracks->Branch("fRhoSnpY", &tracks.fRhoSnpY, "fRhoSnpY/B");
    tTracks->Branch("fRhoSnpZ", &tracks.fRhoSnpZ, "fRhoSnpZ/B");
    tTracks->Branch("fRhoTglY", &tracks.fRhoTglY, "fRhoTglY/B");
    tTracks->Branch("fRhoTglZ", &tracks.fRhoTglZ, "fRhoTglZ/B");
    tTracks->Branch("fRhoTglSnp", &tracks.fRhoTglSnp, "fRhoTglSnp/B");
    tTracks->Branch("fRho1PtY", &tracks.fRho1PtY, "fRho1PtY/B");
    tTracks->Branch("fRho1PtZ", &tracks.fRho1PtZ, "fRho1PtZ/B");
    tTracks->Branch("fRho1PtSnp", &tracks.fRho1PtSnp, "fRho1PtSnp/B");
    tTracks->Branch("fRho1PtTgl", &tracks.fRho1PtTgl, "fRho1PtTgl/B");
    //
    tTracks->Branch("fTPCInnerParam", &tracks.fTPCinnerP, "fTPCInnerParam/F");
    tTracks->Branch("fFlags", &tracks.fFlags, "fFlags/l");
    tTracks->Branch("fITSClusterMap", &tracks.fITSClusterMap, "fITSClusterMap/b");
    tTracks->Branch("fTPCNClsFindable", &tracks.fTPCNClsFindable, "fTPCNClsFindable/b");
    tTracks->Branch("fTPCNClsFindableMinusFound",&tracks.fTPCNClsFindableMinusFound, "fTPCNClsFindableMinusFound/B");
    tTracks->Branch("fTPCNClsFindableMinusCrossedRows", &tracks.fTPCNClsFindableMinusCrossedRows, "fTPCNClsFindableMinusCrossedRows/B");
    tTracks->Branch("fTPCNClsShared", &tracks.fTPCNClsShared, "fTPCNClsShared/b");
    tTracks->Branch("fTRDTOFPattern", &tracks.fTRDTOFPattern, "fTRDTOFPattern/b");
    tTracks->Branch("fITSChi2NCl", &tracks.fITSChi2NCl, "fITSChi2NCl/F");
    tTracks->Branch("fTPCChi2NCl", &tracks.fTPCChi2NCl, "fTPCChi2NCl/F");
    tTracks->Branch("fTRDChi2", &tracks.fTRDChi2, "fTRDChi2/F");
    tTracks->Branch("fTOFChi2", &tracks.fTOFChi2, "fTOFChi2/F");
    tTracks->Branch("fTPCSignal", &tracks.fTPCSignal, "fTPCSignal/F");
    tTracks->Branch("fTRDSignal", &tracks.fTRDSignal, "fTRDSignal/F");
    tTracks->Branch("fTOFSignal", &tracks.fTOFSignal, "fTOFSignal/F");
    tTracks->Branch("fLength", &tracks.fLength, "fLength/F");
    tTracks->Branch("fTOFExpMom", &tracks.fTOFExpMom, "fTOFExpMom/F");
  }
  PostTree(kTracks);

  // Associate branches for Calo
  TTree* tCalo = CreateTree(kCalo);
  tCalo->SetAutoFlush(fNumberOfEventsPerCluster);
  if (fTreeStatus[kCalo]) {
    tCalo->Branch("fBCsID", &calo.fBCsID, "fBCsID/I");
    tCalo->Branch("fCellNumber", &calo.fCellNumber, "fCellNumber/S");
    tCalo->Branch("fAmplitude", &calo.fAmplitude, "fAmplitude/F");
    tCalo->Branch("fTime", &calo.fTime, "fTime/F");
    tCalo->Branch("fCellType", &calo.fCellType, "fCellType/C");
    tCalo->Branch("fCaloType", &calo.fCaloType, "fCaloType/B");
  }
  PostTree(kCalo);

  TTree *tCaloTrigger = CreateTree(kCaloTrigger);
  tCaloTrigger->SetAutoFlush(fNumberOfEventsPerCluster);
  if (fTreeStatus[kCaloTrigger]) {
    tCaloTrigger->Branch("fBCsID", &calotrigger.fBCsID, "fBCsID/I");
    tCaloTrigger->Branch("fFastOrAbsID", &calotrigger.fFastOrAbsID, "fFastOrAbsID/S");
    tCaloTrigger->Branch("fL0Amplitude", &calotrigger.fL0Amplitude, "fL0Amplitude/F");
    tCaloTrigger->Branch("fL1TimeSum", &calotrigger.fL1TimeSum, "fL1TimeSum/F");
    tCaloTrigger->Branch("fNL0Times", &calotrigger.fNL0Times, "fNL0Times/C");
    tCaloTrigger->Branch("fTriggerBits", &calotrigger.fTriggerBits, "fTriggerBits/I");
    tCaloTrigger->Branch("fCaloType", &calotrigger.fCaloType, "fCaloType/B");
  }
  PostTree(kCaloTrigger);

  // Associuate branches for MUON tracks
  TTree* tMuon = CreateTree(kMuon);
  tMuon->SetAutoFlush(fNumberOfEventsPerCluster);
  if (fTreeStatus[kMuon]) {
    tMuon->Branch("fBCsID", &muons.fBCsID, "fBCsID/I");
//    tMuon->Branch("fClusterIndex", &muons.fClusterIndex, "fClusterIndex/I");
//    tMuon->Branch("fNclusters", &muons.fNclusters, "fNclusters/I");
    tMuon->Branch("fInverseBendingMomentum", &muons.fInverseBendingMomentum, "fInverseBendingMomentum/F");
    tMuon->Branch("fThetaX", &muons.fThetaX, "fThetaX/F");
    tMuon->Branch("fThetaY", &muons.fThetaY, "fThetaY/F");
    tMuon->Branch("fZMu", &muons.fZMu, "fZMu/F");
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
    tZdc->Branch("fBCsID",           &zdc.fBCsID          , "fBCsID/I");
    tZdc->Branch("fEnergyZEM1",      &zdc.fEnergyZEM1     , "fEnergyZEM1/F");
    tZdc->Branch("fEnergyZEM2",      &zdc.fEnergyZEM2     , "fEnergyZEM2/F");
    tZdc->Branch("fEnergyCommonZNA", &zdc.fEnergyCommonZNA, "fEnergyCommonZNA/F");
    tZdc->Branch("fEnergyCommonZNC", &zdc.fEnergyCommonZNC, "fEnergyCommonZNC/F");
    tZdc->Branch("fEnergyCommonZPA", &zdc.fEnergyCommonZPA, "fEnergyCommonZPA/F");
    tZdc->Branch("fEnergyCommonZPC", &zdc.fEnergyCommonZPC, "fEnergyCommonZPC/F");
    tZdc->Branch("fEnergySectorZNA", &zdc.fEnergySectorZNA, "fEnergySectorZNA[4]/F");
    tZdc->Branch("fEnergySectorZNC", &zdc.fEnergySectorZNC, "fEnergySectorZNC[4]/F");
    tZdc->Branch("fEnergySectorZPA", &zdc.fEnergySectorZPA, "fEnergySectorZPA[4]/F");
    tZdc->Branch("fEnergySectorZPC", &zdc.fEnergySectorZPC, "fEnergySectorZPC[4]/F");
    tZdc->Branch("fTimeZEM1",        &zdc.fTimeZEM1       , "fTimeZEM1/F");
    tZdc->Branch("fTimeZEM2",        &zdc.fTimeZEM2       , "fTimeZEM2/F");
    tZdc->Branch("fTimeZNA",         &zdc.fTimeZNA        , "fTimeZNA/F");
    tZdc->Branch("fTimeZNC",         &zdc.fTimeZNC        , "fTimeZNC/F");
    tZdc->Branch("fTimeZPA",         &zdc.fTimeZPA        , "fTimeZPA/F");
    tZdc->Branch("fTimeZPC",         &zdc.fTimeZPC        , "fTimeZPC/F");
  }
  
  PostTree(kZdc);

  // Associuate branches for VZERO
  TTree* tVzero = CreateTree(kRun2V0);
  tVzero->SetAutoFlush(fNumberOfEventsPerCluster);
  if (fTreeStatus[kRun2V0]) {
    tVzero->Branch("fBCsID", &vzero.fBCsID, "fBCsID/I");
    tVzero->Branch("fAdc", vzero.fAdc, "fAdc[64]/F");
    tVzero->Branch("fTime", vzero.fTime, "fTime[64]/F");
    tVzero->Branch("fWidth", vzero.fWidth, "fWidth[64]/F");
    tVzero->Branch("fMultA", &vzero.fMultA, "fMultA/F");
    tVzero->Branch("fMultC", &vzero.fMultC, "fMultC/F");
    tVzero->Branch("fTimeA", &vzero.fTimeA, "fTimeA/F");
    tVzero->Branch("fTimeC", &vzero.fTimeC, "fTimeC/F");
    tVzero->Branch("fBBFlag", &vzero.fBBFlag, "fBBFlag/l");
    tVzero->Branch("fBGFlag", &vzero.fBGFlag, "fBGFlag/l");
  }
  PostTree(kRun2V0);

  // Associate branches for FDD (AD)
  TTree* tFDD = CreateTree(kFDD);
  tFDD->SetAutoFlush(fNumberOfEventsPerCluster);
  if (fTreeStatus[kFDD]) {
    tFDD->Branch("fBCsID", &fdd.fBCsID, "fBCsID/I");
    tFDD->Branch("fAmplitude", fdd.fAmplitude, "fAmplitude[8]/F");
    tFDD->Branch("fTimeA", &fdd.fTimeA, "fTimeA/F");
    tFDD->Branch("fTimeC", &fdd.fTimeC, "fTimeC/F");
    tFDD->Branch("fBCSignal", &fdd.fBCSignal, "fBCSignal/b");
  }
  PostTree(kFDD);

  
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
    TTree * tMCvtx = CreateTree(kMcCollision);
    tMCvtx->SetAutoFlush(fNumberOfEventsPerCluster);
    if(fTreeStatus[kMcCollision]) {
      tMCvtx->Branch("fBCsID", &mccollision.fBCsID, "fBCsID/I");
      tMCvtx->Branch("fGeneratorsID", &mccollision.fGeneratorsID, "fGeneratorsID/S");
      tMCvtx->Branch("fPosX", &mccollision.fPosX, "fPosX/F");
      tMCvtx->Branch("fPosY", &mccollision.fPosY, "fPosY/F");
      tMCvtx->Branch("fPosZ", &mccollision.fPosZ, "fPosZ/F");
      tMCvtx->Branch("fT", &mccollision.fT, "fT/F");
      tMCvtx->Branch("fWeight", &mccollision.fWeight, "fWeight/F");
    }
    PostTree(kMcCollision);

    // Associate branches for Kinematics
    TTree* Kinematics = CreateTree(kMcParticle);
    Kinematics->SetAutoFlush(fNumberOfEventsPerCluster);
    if (fTreeStatus[kMcParticle]) {
      Kinematics->Branch("fMcCollisionsID", &mcparticle.fMcCollisionsID, "fMcCollisionsID/I");

      Kinematics->Branch("fPdgCode", &mcparticle.fPdgCode, "fPdgCode/I");
      Kinematics->Branch("fStatusCode", &mcparticle.fStatusCode, "fStatusCode/I");
      Kinematics->Branch("fFlags", &mcparticle.fFlags, "fFlags/b");
      
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
    PostTree(kMcParticle);

    // MC labels of each reconstructed track
    TTree* tLabels = CreateTree(kMcTrackLabel);
    tLabels->SetAutoFlush(fNumberOfEventsPerCluster);
    if (fTreeStatus[kMcTrackLabel]) {
      tLabels->Branch("fLabel", &mctracklabel.fLabel, "fLabel/i");
      tLabels->Branch("fLabelMask", &mctracklabel.fLabelMask, "fLabelMask/s");
    }
    PostTree(kMcTrackLabel);

    // MC labels of each reconstructed calo cluster
    TTree* tCaloLabels = CreateTree(kMcCaloLabel);
    tCaloLabels->SetAutoFlush(fNumberOfEventsPerCluster);
    if (fTreeStatus[kMcCaloLabel]) {
      tCaloLabels->Branch("fLabel", &mccalolabel.fLabel, "fLabel/i");
      tCaloLabels->Branch("fLabelMask", &mccalolabel.fLabelMask, "fLabelMask/s");
    }
    PostTree(kMcCaloLabel);

    // MC labels of each reconstructed calo cluster
    TTree* tCollisionLabels = CreateTree(kMcCollisionLabel);
    tCollisionLabels->SetAutoFlush(fNumberOfEventsPerCluster);
    if (fTreeStatus[kMcCaloLabel]) {
      tCollisionLabels->Branch("fLabel", &mccollisionlabel.fLabel, "fLabel/i");
      tCollisionLabels->Branch("fLabelMask", &mccollisionlabel.fLabelMask, "fLabelMask/s");
    }
    PostTree(kMcCaloLabel);
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
  // Set up the precision masks used to truncate the corresponding float data members

  // Without truncation
  
  UInt_t mCollisionPosition = 0xFFFFFFFF;    // Position in x,y,z
  UInt_t mCollisionPositionCov = 0xFFFFFFFF; // Covariance matrix and chi2

  UInt_t mTrackX =  0xFFFFFFFF;
  UInt_t mTrackAlpha = 0xFFFFFFFF;
  UInt_t mtrackSnp = 0xFFFFFFFF;
  UInt_t mTrackTgl = 0xFFFFFFFF;
  UInt_t mTrack1Pt = 0xFFFFFFFF; // Including the momentun at the inner wall of TPC
  UInt_t mTrackCovDiag = 0xFFFFFFFF; // Including the chi2
  UInt_t mTrackCovOffDiag = 0xFFFFFFFF;
  UInt_t mTrackSignal = 0xFFFFFFFF; // PID signals and track length

  UInt_t mTracklets = 0xFFFFFFFF; // tracklet members

  UInt_t mMcParticleW   = 0xFFFFFFFF; // Precision for weight
  UInt_t mMcParticlePos = 0xFFFFFFFF; // Precision for (x,y,z,t)
  UInt_t mMcParticleMom = 0xFFFFFFFF; // Precision for (Px,Py,Pz,E)

  UInt_t mCaloAmp = 0xFFFFFFFF;
  UInt_t mCaloTime = 0xFFFFFFFF;

  UInt_t mMuonTr1P = 0xFFFFFFFF;
  UInt_t mMuonTrThetaX = 0xFFFFFFFF;
  UInt_t mMuonTrThetaY = 0xFFFFFFFF;
  UInt_t mMuonTrZmu = 0xFFFFFFFF;
  UInt_t mMuonTrBend = 0xFFFFFFFF;
  UInt_t mMuonTrNonBend = 0xFFFFFFFF;
  UInt_t mMuonTrCov = 0xFFFFFFFF; // Covariance matrix and chi2

  UInt_t mMuonCl = 0xFFFFFFFF; // Position and charge
  UInt_t mMuonClErr = 0xFFFFFFFF;

  UInt_t mADTime = 0xFFFFFFFF;
  
  // No compression for ZDC and Run2 VZERO for the moment

  if (fTruncate) {
    mCollisionPosition = 0xFFFFFFF0; // 19 bits mantissa
    mCollisionPositionCov = 0xFFFFE000; // 10 bits mantissa

    mTrackX =  0xFFFFFFF0; // 19 bits
    mTrackAlpha = 0xFFFFFFF0; // 19 bits
    mtrackSnp = 0xFFFFFF00; // 15 bits
    mTrackTgl = 0xFFFFFF00; // 15 bits
    mTrack1Pt = 0xFFFFFC00; // 13 bits
    mTrackCovDiag = 0xFFFFFF00; // 15 bits
    mTrackCovOffDiag = 0xFFFF0000; // 7 bits
    mTrackSignal = 0xFFFFFF00; // 15 bits

    mTracklets = 0xFFFFFF00; // 15 bits

    mMcParticleW   = 0xFFFFFFF0; // 19 bits
    mMcParticlePos = 0xFFFFFFF0; // 19 bits
    mMcParticleMom = 0xFFFFFFF0; // 19 bits

    mCaloAmp = 0xFFFFFF00; // 15 bits
    mCaloTime = 0xFFFFFF00; // 15 bits

    mMuonTr1P = 0xFFFFFC00; // 13 bits 
    mMuonTrThetaX = 0xFFFFFF00; // 15 bits
    mMuonTrThetaY = 0xFFFFFF00; // 15 bits
    mMuonTrZmu = 0xFFFFFFF0; // 19 bits
    mMuonTrBend = 0xFFFFFFF0; // 19 bits
    mMuonTrNonBend = 0xFFFFFFF0; // 19 bits
    mMuonTrCov = 0xFFFF0000; // 7 bits

    mMuonCl = 0xFFFFFF00; // 15 bits
    mMuonClErr = 0xFFFF0000; // 7 bits
    
    mADTime = 0xFFFFF000; // 11 bits
  }
  
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
  
  // bypass vertex selection for muon UPC triggers with at least one muon track
  if (!(fESD->GetFiredTriggerClasses().Contains("CMUP") && fESD->GetNumberOfMuonTracks()>0)) {
    if(pvtx->IsFromVertexer3D() || pvtx->IsFromVertexerZ()) return;
    if(pvtx->GetNContributors()<2) return;
  }
  Int_t eventID = fEventCount++;

  //---------------------------------------------------------------------------
  // Collision data

  // Adjust start indices for this event in all trees by adding the number of entries of the previous event
  for (auto i = 0; i < kTrees; ++i)
     eventextra.fStart[i] += eventextra.fNentries[i];

  eventextra.fNentries[kEvents] = 1;  // one entry per vertex
  collision.fBCsID = eventID;
  collision.fPosX = AliMathBase::TruncateFloatFraction(pvtx->GetX(), mCollisionPosition);
  collision.fPosY = AliMathBase::TruncateFloatFraction(pvtx->GetY(), mCollisionPosition);
  collision.fPosZ = AliMathBase::TruncateFloatFraction(pvtx->GetZ(), mCollisionPosition);

  Double_t covmatrix[6];
  pvtx->GetCovMatrix(covmatrix);

  collision.fCovXX = AliMathBase::TruncateFloatFraction(covmatrix[0], mCollisionPositionCov);
  collision.fCovXY = AliMathBase::TruncateFloatFraction(covmatrix[1], mCollisionPositionCov);
  collision.fCovXZ = AliMathBase::TruncateFloatFraction(covmatrix[2], mCollisionPositionCov);
  collision.fCovYY = AliMathBase::TruncateFloatFraction(covmatrix[3], mCollisionPositionCov);
  collision.fCovYZ = AliMathBase::TruncateFloatFraction(covmatrix[4], mCollisionPositionCov);
  collision.fCovZZ = AliMathBase::TruncateFloatFraction(covmatrix[5], mCollisionPositionCov);

  collision.fChi2 = AliMathBase::TruncateFloatFraction(pvtx->GetChi2(), mCollisionPositionCov);
  collision.fN = (pvtx->GetNDF()+3)/2;

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
      SETBIT(collision.fCollisionTimeMask, 0);
    else
      CLRBIT(collision.fCollisionTimeMask, 0);
    //
    if (TOFResponse.GetStartTimeMask(mom) & 0x2)
      SETBIT(collision.fCollisionTimeMask, 1);
    else
      CLRBIT(collision.fCollisionTimeMask, 1);
    //
    if (TOFResponse.GetStartTimeMask(mom) & 0x3)
      SETBIT(collision.fCollisionTimeMask, 2);
    else
      CLRBIT(collision.fCollisionTimeMask, 2);
  }

  // Recalculate unique event time and its resolution
  collision.fCollisionTime = AliMathBase::TruncateFloatFraction(TMath::Mean(10,eventTime,eventTimeWeight), mCollisionPosition); // Weighted mean of times per momentum interval
  collision.fCollisionTimeRes = AliMathBase::TruncateFloatFraction(TMath::Sqrt(9./10.)*TMath::Mean(10,eventTimeRes), mCollisionPositionCov); // PH bad approximation

  //---------------------------------------------------------------------------
  // BC data
  
  bc.fRunNumber = fESD->GetRunNumber();
  
  ULong64_t evtid = GetEventIdAsLong(fESD->GetHeader());
  if(!evtid){
    evtid = (ULong64_t(fESD->GetTimeStamp())<<32) + ULong64_t((fESD->GetNumberOfTPCClusters()<<5)|(fESD->GetNumberOfTPCTracks()));
  }
  bc.fGlobalBC = evtid;
  
  bc.fTriggerMask = fESD->GetTriggerMask();
  TString firedClasses = fESD->GetFiredTriggerClasses();
  if (firedClasses.Contains("CINT7-B-NOPF-CENTNOTRD")) bc.fTriggerMask |= 1ull << 50;
  if (firedClasses.Contains("CCUP8-B-NOPF-CENTNOTRD")) bc.fTriggerMask |= 1ull << 51;
  if (firedClasses.Contains("CCUP9-B-NOPF-CENTNOTRD")) bc.fTriggerMask |= 1ull << 52;
  if (firedClasses.Contains("CMUP10-B-NOPF-MUFAST"))   bc.fTriggerMask |= 1ull << 53;
  if (firedClasses.Contains("CMUP11-B-NOPF-MUFAST"))   bc.fTriggerMask |= 1ull << 54;
  if (firedClasses.Contains("CINT7-B-NOPF-MUFAST"))    bc.fTriggerMask |= 1ull << 55;
  if (firedClasses.Contains("CMSL7-B-NOPF-MUFAST"))    bc.fTriggerMask |= 1ull << 56;
  if (firedClasses.Contains("CMLL7-B-NOPF-MUFAST"))    bc.fTriggerMask |= 1ull << 57;
  if (firedClasses.Contains("CMUL7-B-NOPF-MUFAST"))    bc.fTriggerMask |= 1ull << 58;
  
  FillTree(kBC);
  
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
    tracks.fTrackType = TrackTypeEnum::GlobalTrack;

    tracks.fX = AliMathBase::TruncateFloatFraction(track->GetX(), mTrackX);
    tracks.fAlpha = AliMathBase::TruncateFloatFraction(track->GetAlpha(), mTrackAlpha);

    tracks.fY = track->GetY(); // no lossy compression
    tracks.fZ = track->GetZ();
    tracks.fSnp = AliMathBase::TruncateFloatFraction(track->GetSnp(), mtrackSnp);
    tracks.fTgl = AliMathBase::TruncateFloatFraction(track->GetTgl(), mTrackTgl);
    tracks.fSigned1Pt = AliMathBase::TruncateFloatFraction(track->GetSigned1Pt(), mTrack1Pt);

    // Modified covariance matrix
    // First sigmas on the diagonal
    tracks.fSigmaY = AliMathBase::TruncateFloatFraction(TMath::Sqrt(track->GetSigmaY2()), mTrackCovDiag);
    tracks.fSigmaZ = AliMathBase::TruncateFloatFraction(TMath::Sqrt(track->GetSigmaZ2()), mTrackCovDiag);
    tracks.fSigmaSnp = AliMathBase::TruncateFloatFraction(TMath::Sqrt(track->GetSigmaSnp2()), mTrackCovDiag);
    tracks.fSigmaTgl = AliMathBase::TruncateFloatFraction(TMath::Sqrt(track->GetSigmaTgl2()), mTrackCovDiag);
    tracks.fSigma1Pt = AliMathBase::TruncateFloatFraction(TMath::Sqrt(track->GetSigma1Pt2()), mTrackCovDiag);
    //
    tracks.fRhoZY = (Char_t)(128.*track->GetSigmaZY()/tracks.fSigmaZ/tracks.fSigmaY);
    tracks.fRhoSnpY = (Char_t)(128.*track->GetSigmaSnpY()/tracks.fSigmaSnp/tracks.fSigmaY);
    tracks.fRhoSnpZ = (Char_t)(128.*track->GetSigmaSnpZ()/tracks.fSigmaSnp/tracks.fSigmaZ);
    tracks.fRhoTglY = (Char_t)(128.*track->GetSigmaTglY()/tracks.fSigmaTgl/tracks.fSigmaY);
    tracks.fRhoTglZ = (Char_t)(128.*track->GetSigmaTglZ()/tracks.fSigmaTgl/tracks.fSigmaZ);
    tracks.fRhoTglSnp = (Char_t)(128.*track->GetSigmaTglSnp()/tracks.fSigmaTgl/tracks.fSigmaSnp);
    tracks.fRho1PtY = (Char_t)(128.*track->GetSigma1PtY()/tracks.fSigma1Pt/tracks.fSigmaY);
    tracks.fRho1PtZ = (Char_t)(128.*track->GetSigma1PtZ()/tracks.fSigma1Pt/tracks.fSigmaZ);
    tracks.fRho1PtSnp = (Char_t)(128.*track->GetSigma1PtSnp()/tracks.fSigma1Pt/tracks.fSigmaSnp);
    tracks.fRho1PtTgl = (Char_t)(128.*track->GetSigma1PtTgl()/tracks.fSigma1Pt/tracks.fSigmaTgl);

    const AliExternalTrackParam *intp = track->GetTPCInnerParam();
    tracks.fTPCinnerP = AliMathBase::TruncateFloatFraction((intp ? intp->GetP() : 0), mTrack1Pt); // Set the momentum to 0 if the track did not reach TPC

    tracks.fFlags = track->GetStatus();

    tracks.fITSClusterMap = track->GetITSClusterMap();
    tracks.fTPCNClsFindable = track->GetTPCNclsF();
    tracks.fTPCNClsFindableMinusFound = tracks.fTPCNClsFindable - track->GetTPCNcls();
    tracks.fTPCNClsFindableMinusCrossedRows = tracks.fTPCNClsFindable - track->GetTPCCrossedRows();
    tracks.fTPCNClsShared = (track->GetTPCSharedMap()).CountBits();
    tracks.fTRDTOFPattern = 0; // FIXME

    tracks.fITSChi2NCl = AliMathBase::TruncateFloatFraction((track->GetITSNcls() ? track->GetITSchi2() / track->GetITSNcls() : 0), mTrackCovOffDiag);
    tracks.fTPCChi2NCl = AliMathBase::TruncateFloatFraction((track->GetTPCNcls() ? track->GetTPCchi2() / track->GetTPCNcls() : 0), mTrackCovOffDiag);
    tracks.fTRDChi2 = AliMathBase::TruncateFloatFraction(track->GetTRDchi2(), mTrackCovOffDiag);
    tracks.fTOFChi2 = AliMathBase::TruncateFloatFraction(track->GetTOFchi2(), mTrackCovOffDiag);

    tracks.fTPCSignal = AliMathBase::TruncateFloatFraction(track->GetTPCsignal(), mTrackSignal);
    tracks.fTRDSignal = AliMathBase::TruncateFloatFraction(track->GetTRDsignal(), mTrackSignal);
    tracks.fTOFSignal = AliMathBase::TruncateFloatFraction(track->GetTOFsignal(), mTrackSignal);
    tracks.fLength = AliMathBase::TruncateFloatFraction(track->GetIntegratedLength(), mTrackSignal);

    // Speed of ligth in TOF units
    const Float_t cspeed = 0.029979246f;
    // PID hypothesis for the momentum extraction
    const AliPID::EParticleType tof_pid = AliPID::kPion;
    // Expected beta for such hypothesis
    const Float_t exp_beta =
        (track->GetIntegratedLength() /
         TOFResponse.GetExpectedSignal(track, tof_pid) / cspeed);

    tracks.fTOFExpMom = AliMathBase::TruncateFloatFraction(
        AliPID::ParticleMass(tof_pid) * exp_beta * cspeed /
            TMath::Sqrt(1. - (exp_beta * exp_beta)),
        mTrack1Pt);

    if (fTaskMode == kMC) {
      // Separate tables (trees) for the MC labels
      Int_t alabel = track->GetLabel();
      mctracklabel.fLabel = TMath::Abs(alabel) + fOffsetLabel;
      mctracklabel.fLabelMask = 0;
      // Use the ITS shared clusters to set the corresponding bits 0-6
      UChar_t itsMask = track->GetITSSharedMap() & 0x1F; // Normally only bits 0-5 are set in Run1/2
      mctracklabel.fLabelMask |= itsMask;
      // Use the number of TPC shared clusters as number of TPC mismatches
      // encode in bits 7-9 the values in the ranges 0, 1, 2-3, 4-7, 8-15, 16-31, 32-63, >64
      const TBits * tpcShared = track->GetTPCSharedMapPtr();
      UInt_t tpcCount = tpcShared->CountBits();
      UShort_t tpcMask = 0;
      while (tpcCount>0) {
	tpcCount = tpcCount >> 1;
	tpcMask++;
      }
      if (tpcMask>7) tpcMask = 7;
      mctracklabel.fLabelMask |= (tpcMask<<7);
      // TRD (bit 10)
      // We can also use labels per tracklet in the future
      Int_t trdLabel = track->GetTRDLabel();
      if (TMath::Abs(alabel)!=TMath::Abs(trdLabel)) mctracklabel.fLabelMask |= (0x1 << 10);
      // TOF (bit 11)
      Int_t tofLabel[3]={-1};
      track->GetTOFLabel(tofLabel);
      // Check if at least one of the TOF hits matches the track label
      if (!( TMath::Abs(alabel)==TMath::Abs(tofLabel[0])
	     || TMath::Abs(alabel)==TMath::Abs(tofLabel[1])
	     || TMath::Abs(alabel)==TMath::Abs(tofLabel[2])))
	mctracklabel.fLabelMask |= (0x1 << 11);

      if (alabel<0) mctracklabel.fLabelMask |= (0x1 << 15);

      FillTree(kMcTrackLabel);
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

  eventextra.fNentries[kTOF]    = ntofcls_filled;

  AliMultiplicity *mlt = fESD->GetMultiplicity();
  Int_t Ntracklets = mlt->GetNumberOfTracklets();

  Int_t ntracklet_filled = 0;
  Float_t theta, phi, dphi, dphiS, dist, x, tgl, alpha;

  for (Int_t itr = Ntracklets; itr--;) {
    dphi   = mlt->GetDeltaPhi(itr);
    dist   = mlt->CalcDist(itr);
    
    // on-the-fly filtering based on parameters tuned in Run2
    dphiS  = TMath::Abs(dphi) - 0.0045; 
    if (dphi<0) dphiS = -dphiS;
    if (dist<1. && dphiS<0.06) {
      theta  = mlt->GetTheta(itr);
      phi = mlt->GetPhi(itr);
      tracks.fTrackType = TrackTypeEnum::Run2Tracklet;
      
      // inversion formulas for snp and alpha
      tracks.fSnp = 0.;
      alpha = phi;
      tracks.fAlpha = AliMathBase::TruncateFloatFraction(alpha, mTracklets);

      // inversion formulas for tgl
      x = (TMath::Tan(theta/2.)-1.) / (TMath::Tan(theta/2.)+1.);
      if (TMath::Log(TMath::Tan(theta/2)) >= 0)
        tgl = TMath::Sqrt((TMath::Power((1.+TMath::Power(x,2))/(1.-TMath::Power(x,2)),2))-1.);
      else 
        tgl = - TMath::Sqrt((TMath::Power((1.+TMath::Power(x,2))/(1.-TMath::Power(x,2)),2))-1.);
      tracks.fTgl = AliMathBase::TruncateFloatFraction(tgl, mTracklets);
    
      // set global track parameters to NAN
      tracks.fX = NAN;
      tracks.fY = NAN;
      tracks.fZ = NAN; 
      tracks.fSigned1Pt = NAN;
      tracks.fSigmaY = NAN;
      tracks.fSigmaZ = NAN;
      tracks.fSigmaSnp = NAN;
      tracks.fSigmaTgl = NAN;
      tracks.fSigma1Pt = NAN;
      tracks.fRhoZY = NAN;
      tracks.fRhoSnpY = NAN;
      tracks.fRhoSnpZ = NAN;
      tracks.fRhoTglY = NAN;
      tracks.fRhoTglZ = NAN;
      tracks.fRhoTglSnp = NAN;
      tracks.fRho1PtY = NAN;
      tracks.fRho1PtZ = NAN;
      tracks.fRho1PtSnp = NAN;
      tracks.fRho1PtTgl = NAN;
      tracks.fTPCinnerP = NAN; 
      tracks.fFlags = 0;
      tracks.fITSClusterMap = 0;
      tracks.fTPCNClsFindable = 0;
      tracks.fTPCNClsFindableMinusFound = 0;
      tracks.fTPCNClsFindableMinusCrossedRows = 0;
      tracks.fTPCNClsShared = 0;
      tracks.fTRDTOFPattern = 0;
      tracks.fITSChi2NCl = NAN;
      tracks.fTPCChi2NCl = NAN;
      tracks.fTRDChi2 = NAN; 
      tracks.fTOFChi2 = NAN;
      tracks.fTPCSignal = NAN; 
      tracks.fTRDSignal = NAN;
      tracks.fTOFSignal = NAN;
      tracks.fLength = NAN;
      tracks.fTOFExpMom = NAN;

      if (fTaskMode == kMC) {
	// Separate tables (trees) for the MC labels: tracklets
	Int_t alabel = mlt->GetLabel(itr, 0); // Take the label of the first layer
	mctracklabel.fLabel = TMath::Abs(alabel) + fOffsetLabel;
	mctracklabel.fLabelMask = 0;
	// Mask fake tracklets
	if (alabel<0) mctracklabel.fLabelMask |= (0x1 << 15);
	if (mlt->GetLabel(itr, 0) != mlt->GetLabel(itr, 1)) mctracklabel.fLabelMask |= (0x1 << 15);

	FillTree(kMcTrackLabel);
      }

      FillTree(kTracks);
      if (fTreeStatus[kTracks]) ntracklet_filled++;
    }
  } // end loop on tracklets
  eventextra.fNentries[kTracks] = ntrk_filled + ntracklet_filled; 
  
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

    calo.fBCsID = eventID;
    
    cells->GetCell(ice, cellNumber, amplitude, time, mclabel, efrac);
    calo.fCellNumber = cellNumber;
    calo.fAmplitude = AliMathBase::TruncateFloatFraction(amplitude, mCaloAmp);
    calo.fTime = AliMathBase::TruncateFloatFraction(time, mCaloAmp);
    calo.fCaloType = cells->GetType(); // common for all cells
    calo.fCellType = cells->GetHighGain(ice) ? 0. : 1.; 
    FillTree(kCalo);
    if (fTreeStatus[kCalo]) ncalocells_filled++;
    if (fTaskMode == kMC) {
      mccalolabel.fLabel = TMath::Abs(mclabel) + fOffsetLabel;
      mccalolabel.fLabelMask = 0;
      if (mclabel<0) mccalolabel.fLabelMask |= (0x1 << 15);

      FillTree(kMcCaloLabel);
    }
  } // end loop on calo cells
  eventextra.fNentries[kCalo] = ncalocells_filled;

  AliEMCALGeometry *geo = AliEMCALGeometry::GetInstanceFromRunNumber(fESD->GetRunNumber()); // Needed for EMCAL trigger mapping
  AliESDCaloTrigger *calotriggers = fESD->GetCaloTrigger("EMCAL");
  calotriggers->Reset();
  Int_t ncalotriggers_filled = 0; // total number of EMCAL triggers filled per event
  while(calotriggers->Next()){
    calotrigger.fBCsID = eventID;
    int col, row, fastorID;
    calotriggers->GetPosition(col, row);
    geo->GetTriggerMapping()->GetAbsFastORIndexFromPositionInEMCAL(col, row, fastorID);
    calotrigger.fFastOrAbsID = fastorID;
    calotriggers->GetAmplitude(calotrigger.fL0Amplitude);
    calotrigger.fL0Amplitude = AliMathBase::TruncateFloatFraction(calotrigger.fL0Amplitude, mCaloAmp);
    calotriggers->GetTime(calotrigger.fL0Time);
    calotrigger.fL0Time = AliMathBase::TruncateFloatFraction(calotrigger.fL0Time, mCaloTime);
    calotriggers->GetTriggerBits(calotrigger.fTriggerBits);
    Int_t nL0times;
    calotriggers->GetNL0Times(nL0times);
    calotrigger.fNL0Times = nL0times;
    calotriggers->GetL1TimeSum(calotrigger.fTriggerBits);
    calotrigger.fCaloType = 1;
    FillTree(kCaloTrigger);
    if (fTreeStatus[kCaloTrigger]) ncalotriggers_filled++;
  }
  eventextra.fNentries[kCaloTrigger] = ncalotriggers_filled;

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

    calo.fBCsID = eventID;
    
    cells->GetCell(icp, cellNumber, amplitude, time, mclabel, efrac);
    calo.fCellNumber = cellNumber;
    calo.fAmplitude = AliMathBase::TruncateFloatFraction(amplitude, mCaloAmp);
    calo.fTime = AliMathBase::TruncateFloatFraction(time, mCaloTime);
    calo.fCellType = cells->GetHighGain(icp) ? 0. : 1.;     /// @TODO cell type value to be confirmed by PHOS experts
    calo.fCaloType = cells->GetType(); // common for all cells

    FillTree(kCalo);
    if (fTreeStatus[kCalo]) nphoscells_filled++;
    if (fTaskMode == kMC) {
      mccalolabel.fLabel = TMath::Abs(mclabel) + fOffsetLabel;
      mccalolabel.fLabelMask = 0;
      if (mclabel<0) mccalolabel.fLabelMask |= (0x1 << 15);

      FillTree(kMcCaloLabel);
    }
  } // end loop on PHOS cells
  eventextra.fNentries[kCalo] = nphoscells_filled;

  //---------------------------------------------------------------------------
  // Muon tracks
  muons.fBCsID  = eventID;
  
  Int_t nmu = fESD->GetNumberOfMuonTracks();
  Int_t nmu_filled = 0;    // total number of muons filled per event
  Int_t nmucl_filled = 0;  // total number of clusters filled per event
  for (Int_t imu=0; imu<nmu; ++imu) {
    AliESDMuonTrack* mutrk = fESD->GetMuonTrack(imu);

    muons.fInverseBendingMomentum = AliMathBase::TruncateFloatFraction(mutrk->GetInverseBendingMomentum(), mMuonTr1P);
    muons.fThetaX = AliMathBase::TruncateFloatFraction(mutrk->GetThetaX(), mMuonTrThetaX);
    muons.fThetaY = AliMathBase::TruncateFloatFraction(mutrk->GetThetaY(), mMuonTrThetaY);
    muons.fZMu = AliMathBase::TruncateFloatFraction(mutrk->GetZ(), mMuonTrZmu);
    muons.fBendingCoor = AliMathBase::TruncateFloatFraction(mutrk->GetBendingCoor(), mMuonTrBend);
    muons.fNonBendingCoor = AliMathBase::TruncateFloatFraction(mutrk->GetNonBendingCoor(), mMuonTrNonBend);

    TMatrixD cov;
    mutrk->GetCovariances(cov);
    for (Int_t i = 0; i < 5; i++)
      for (Int_t j = 0; j <= i; j++)
	muons.fCovariances[i*(i+1)/2 + j] = AliMathBase::TruncateFloatFraction(cov(i,j), mMuonTrCov);

    muons.fChi2 = AliMathBase::TruncateFloatFraction(mutrk->GetChi2(), mMuonTrCov);
    muons.fChi2MatchTrigger = AliMathBase::TruncateFloatFraction(mutrk->GetChi2MatchTrigger(), mMuonTrCov);

    // Now MUON clusters for the current track
    Int_t muTrackID = fOffsetMuTrackID + imu;
    Int_t nmucl = mutrk->GetNClusters();
    for (Int_t imucl=0; imucl<nmucl; ++imucl){
      AliESDMuonCluster *muCluster = fESD->FindMuonCluster(mutrk->GetClusterId(imucl));
      mucls.fMuonsID = muTrackID;
      mucls.fX = AliMathBase::TruncateFloatFraction(muCluster->GetX(), mMuonCl);
      mucls.fY = AliMathBase::TruncateFloatFraction(muCluster->GetY(), mMuonCl);
      mucls.fZ = AliMathBase::TruncateFloatFraction(muCluster->GetZ(), mMuonCl);
      mucls.fErrX = AliMathBase::TruncateFloatFraction(muCluster->GetErrX(), mMuonClErr);
      mucls.fErrY = AliMathBase::TruncateFloatFraction(muCluster->GetErrY(), mMuonClErr);
      mucls.fCharge = AliMathBase::TruncateFloatFraction(muCluster->GetCharge(), mMuonCl);
      mucls.fChi2   = AliMathBase::TruncateFloatFraction(muCluster->GetChi2(), mMuonClErr);
      FillTree(kMuonCls);
      if (fTreeStatus[kMuonCls]) nmucl_filled++;
    } // End loop on muon clusters for the current muon track

    // In case we need connection to clusters, activate next lines
    // muons.fClusterIndex += muons.fNclusters;
    // muons.fNclusters = nmucl_filled;

    FillTree(kMuon);
    if (fTreeStatus[kMuon]) nmu_filled++;
  } // End loop on muon tracks
  eventextra.fNentries[kMuon] = nmu_filled;
  eventextra.fNentries[kMuonCls] = nmucl_filled;

  //---------------------------------------------------------------------------
  // ZDC
  AliESDZDC* esdzdc  =    fESD->GetESDZDC();
  zdc.fBCsID = eventID;
  // ZEM
  zdc.fEnergyZEM1      = esdzdc->GetZEM1Energy();
  zdc.fEnergyZEM2      = esdzdc->GetZEM2Energy();
  zdc.fEnergyCommonZNA = esdzdc->GetZNATowerEnergy()[0];
  zdc.fEnergyCommonZNC = esdzdc->GetZNCTowerEnergy()[0];
  zdc.fEnergyCommonZPA = esdzdc->GetZPATowerEnergy()[0];
  zdc.fEnergyCommonZPC = esdzdc->GetZPCTowerEnergy()[0];
  
  // ZDC (P,N) sectors
  for (Int_t ich=0; ich<4; ++ich) {
    zdc.fEnergySectorZNA[ich] = esdzdc->GetZNATowerEnergy()[ich+1];
    zdc.fEnergySectorZNC[ich] = esdzdc->GetZNCTowerEnergy()[ich+1];
    zdc.fEnergySectorZPA[ich] = esdzdc->GetZPATowerEnergy()[ich+1];
    zdc.fEnergySectorZPC[ich] = esdzdc->GetZPCTowerEnergy()[ich+1];
  }
  // ZDC TDC
  Bool_t isHitFlagFilled = fESD->GetRunNumber()>=208502;
  Bool_t isZNAhit  = isHitFlagFilled ? esdzdc->IsZNAhit() : 1;
  Bool_t isZNChit  = isHitFlagFilled ? esdzdc->IsZNChit() : 1;
  Bool_t isZPAhit  = isHitFlagFilled ? esdzdc->IsZPAhit() : 1;
  Bool_t isZPChit  = isHitFlagFilled ? esdzdc->IsZPChit() : 1;
  Bool_t isZEM1hit = isHitFlagFilled ? esdzdc->IsZEM1hit() : 1;
  Bool_t isZEM2hit = isHitFlagFilled ? esdzdc->IsZEM2hit() : 1;
  
  zdc.fTimeZNA  = 999.f;
  zdc.fTimeZNC  = 999.f;
  zdc.fTimeZPA  = 999.f;
  zdc.fTimeZPC  = 999.f;
  zdc.fTimeZEM1 = 999.f;
  zdc.fTimeZEM2 = 999.f;

  // Storing first ZDC hit in +/-12.5 ns around 0
  for (Int_t i=0;i<4;i++) {
    Float_t tZNA  = isZNAhit  ? esdzdc->GetZDCTDCCorrected(esdzdc->GetZNATDCChannel(),i)  : 999.f;
    Float_t tZNC  = isZNChit  ? esdzdc->GetZDCTDCCorrected(esdzdc->GetZNCTDCChannel(),i)  : 999.f;
    Float_t tZPA  = isZPAhit  ? esdzdc->GetZDCTDCCorrected(esdzdc->GetZPATDCChannel(),i)  : 999.f;
    Float_t tZPC  = isZPChit  ? esdzdc->GetZDCTDCCorrected(esdzdc->GetZPCTDCChannel(),i)  : 999.f;
    Float_t tZEM1 = isZEM1hit ? esdzdc->GetZDCTDCCorrected(esdzdc->GetZEM1TDCChannel(),i) : 999.f;
    Float_t tZEM2 = isZEM2hit ? esdzdc->GetZDCTDCCorrected(esdzdc->GetZEM2TDCChannel(),i) : 999.f;
    if (tZNA >-12.5 && tZNA <12.5 && zdc.fTimeZNA >998) zdc.fTimeZNA  = tZNA;
    if (tZNC >-12.5 && tZNC <12.5 && zdc.fTimeZNC >998) zdc.fTimeZNC  = tZNC;
    if (tZPA >-12.5 && tZPA <12.5 && zdc.fTimeZPA >998) zdc.fTimeZPA  = tZPA;
    if (tZPC >-12.5 && tZPC <12.5 && zdc.fTimeZPC >998) zdc.fTimeZPC  = tZPC;
    if (tZEM1>-12.5 && tZEM1<12.5 && zdc.fTimeZEM1>998) zdc.fTimeZEM1 = tZEM1;
    if (tZEM2>-12.5 && tZEM2<12.5 && zdc.fTimeZEM2>998) zdc.fTimeZEM2 = tZEM2;
  }
  
  FillTree(kZdc);
  if (fTreeStatus[kZdc]) eventextra.fNentries[kZdc] = 1;

  //---------------------------------------------------------------------------
  // VZERO
  AliESDVZERO * vz = fESD->GetVZEROData();
  vzero.fBCsID  = eventID;
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
  vzero.fMultA = vz->GetMTotV0A();
  vzero.fMultC = vz->GetMTotV0C();
  vzero.fTimeA = vz->GetV0ATime();
  vzero.fTimeC = vz->GetV0CTime();
  FillTree(kRun2V0);
  if (fTreeStatus[kRun2V0]) eventextra.fNentries[kRun2V0] = 1;

  //---------------------------------------------------------------------------
  // AD (FDD)
  AliESDAD* esdad = fESD->GetADData();
  fdd.fBCsID = eventID;
  fdd.fTimeA = AliMathBase::TruncateFloatFraction(esdad->GetADATime(),mADTime);
  fdd.fTimeC = AliMathBase::TruncateFloatFraction(esdad->GetADCTime(),mADTime);
  FillTree(kFDD);
  if (fTreeStatus[kFDD]) eventextra.fNentries[kFDD] = 1;
  
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
  eventextra.fNentries[kV0s] = nv0_filled;

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
  eventextra.fNentries[kCascades] = ncascades_filled;
  
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

      mcparticle.fMcCollisionsID = eventID;
      
      //Get the kinematic values of the particles
      mcparticle.fPdgCode = particle->GetPdgCode();
      mcparticle.fStatusCode = particle->GetStatusCode();
      mcparticle.fFlags = 0;
      if (i >= MCEvt->Stack()->GetNprimary())
        mcparticle.fFlags |= MCParticleFlags::ProducedInTransport;
      mcparticle.fMother[0] = vpt->GetMother();
      if (mcparticle.fMother[0] > -1) mcparticle.fMother[0]+=fOffsetLabel;
      mcparticle.fMother[1] = -1;
      mcparticle.fDaughter[0] = particle->GetFirstDaughter();
      if (mcparticle.fDaughter[0] > -1) mcparticle.fDaughter[0]+=fOffsetLabel;
      mcparticle.fDaughter[1] = particle->GetLastDaughter();
      if (mcparticle.fDaughter[1] > -1) mcparticle.fDaughter[1]+=fOffsetLabel;
      mcparticle.fWeight = AliMathBase::TruncateFloatFraction(particle->GetWeight(), mMcParticleW);

      mcparticle.fPx = AliMathBase::TruncateFloatFraction(particle->Px(), mMcParticleMom);
      mcparticle.fPy = AliMathBase::TruncateFloatFraction(particle->Py(), mMcParticleMom);
      mcparticle.fPz = AliMathBase::TruncateFloatFraction(particle->Pz(), mMcParticleMom);
      mcparticle.fE  = AliMathBase::TruncateFloatFraction(particle->Energy(), mMcParticleMom);

      mcparticle.fVx = AliMathBase::TruncateFloatFraction(particle->Vx(), mMcParticlePos);
      mcparticle.fVy = AliMathBase::TruncateFloatFraction(particle->Vy(), mMcParticlePos);
      mcparticle.fVz = AliMathBase::TruncateFloatFraction(particle->Vz(), mMcParticlePos);
      mcparticle.fVt = AliMathBase::TruncateFloatFraction(particle->T(), mMcParticlePos);

      FillTree(kMcParticle);
      if (fTreeStatus[kMcParticle]) nkine_filled++;
    }
    fOffsetLabel += nMCtracks; // Offset for the labels of the next event
  }
  eventextra.fNentries[kMcParticle] = nkine_filled;

  if (MCEvt) {
    // MC vertex
    const AliVVertex* MCvtx = MCEvt->GetPrimaryVertex();
    if (!MCvtx) //Check on the MC vertex
      AliFatal("Could not retrieve MC vertex");

    mccollision.fBCsID = eventID;

    mccollision.fPosX = AliMathBase::TruncateFloatFraction(MCvtx->GetX(), mCollisionPosition);
    mccollision.fPosY = AliMathBase::TruncateFloatFraction(MCvtx->GetY(), mCollisionPosition);
    mccollision.fPosZ = AliMathBase::TruncateFloatFraction(MCvtx->GetZ(), mCollisionPosition);

    AliGenEventHeader* mcGenH = MCEvt->GenEventHeader();
    mccollision.fT = mcGenH->InteractionTime();
    mccollision.fWeight = mcGenH->EventWeight();

    mccollision.fGeneratorsID = 0;
    for (Int_t gen = 0; gen < kGenerators; gen++) {
      if (mcGenH->InheritsFrom(Generator[gen]))
        SETBIT(mccollision.fGeneratorsID, gen);
      else
        CLRBIT(mccollision.fGeneratorsID, gen);
    }
    if (mcGenH->InheritsFrom(Generator[kAliGenCocktailEventHeader])) {
      TList* headers = ((AliGenCocktailEventHeader*)mcGenH)->GetHeaders();
      for (Int_t cocktail = 0; cocktail < headers->GetEntries(); headers++) {
        for (Int_t gen = 0; gen < kGenerators; gen++) {
          if (mcGenH->InheritsFrom(Generator[gen]))
            SETBIT(mccollision.fGeneratorsID, gen);
        }
      }
    }
    eventextra.fNentries[kMcCollision] = 1;
  } else {
    eventextra.fNentries[kMcCollision] = 0;
  }
  // Filling the tree of vertices has to be done last because it contains the
  // index data for the other trees
  FillTree(kMcCollision);

  // MC collision label
  mccollisionlabel.fLabel = eventID;
  mccollisionlabel.fLabelMask = 0;
  FillTree(kMcCollisionLabel);

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
