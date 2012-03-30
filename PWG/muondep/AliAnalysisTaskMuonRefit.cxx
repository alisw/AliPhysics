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

/* $Id$ */

// ROOT includes
#include <TString.h>
#include <TList.h>
#include <TGeoManager.h>

// STEER includes
#include "AliESDEvent.h"
#include "AliESDMuonTrack.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliGeomManager.h"

// ANALYSIS includes
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliAnalysisTaskMuonRefit.h"

// MUON includes
#include "AliMUONCDB.h"
#include "AliMUONConstants.h"
#include "AliMUONRecoParam.h"
#include "AliMUONESDInterface.h"
#include "AliMUONRefitter.h"
#include "AliMUONTrack.h"
#include "AliMUONTrackParam.h"
#include "AliMUONVCluster.h"
#include "AliMUONVTrackStore.h"
#include "AliMUONLocalTrigger.h"
#include "AliMUONGeometryTransformer.h"

#ifndef SafeDelete
#define SafeDelete(x) if (x != NULL) { delete x; x = NULL; }
#endif

ClassImp(AliAnalysisTaskMuonRefit)

//________________________________________________________________________
AliAnalysisTaskMuonRefit::AliAnalysisTaskMuonRefit() :
AliAnalysisTaskSE(),
fDefaultStorage(""),
fImproveTracks(kFALSE),
fSigmaCut(-1.),
fSigmaCutForTrigger(-1.),
fReAlign(kFALSE),
fOldAlignStorage(""),
fNewAlignStorage(""),
fOldGeoTransformer(NULL),
fNewGeoTransformer(NULL),
fESDInterface(NULL),
fRefitter(NULL)
{
  /// Default constructor
  for (Int_t i = 0; i < 10; i++) ResetClusterResolution(i, -1., -1.);
}

//________________________________________________________________________
AliAnalysisTaskMuonRefit::AliAnalysisTaskMuonRefit(const char *name) :
AliAnalysisTaskSE(name),
fDefaultStorage("raw://"),
fImproveTracks(kFALSE),
fSigmaCut(-1.),
fSigmaCutForTrigger(-1.),
fReAlign(kFALSE),
fOldAlignStorage(""),
fNewAlignStorage(""),
fOldGeoTransformer(NULL),
fNewGeoTransformer(NULL),
fESDInterface(NULL),
fRefitter(NULL)
{
  /// Constructor
  for (Int_t i = 0; i < 10; i++) ResetClusterResolution(i, -1., -1.);
}

//________________________________________________________________________
AliAnalysisTaskMuonRefit::~AliAnalysisTaskMuonRefit()
{
  /// Destructor
  SafeDelete(fOldGeoTransformer);
  SafeDelete(fNewGeoTransformer);
  SafeDelete(fESDInterface);
  SafeDelete(fRefitter);
}

//___________________________________________________________________________
void AliAnalysisTaskMuonRefit::UserCreateOutputObjects()
{
}

//________________________________________________________________________
void AliAnalysisTaskMuonRefit::UserExec(Option_t *)
{
  /// Main event loop
  
  // check if refitter properly created
  if (!fRefitter) return;
  
  AliESDEvent* esd = dynamic_cast<AliESDEvent*>(InputEvent());
  if (!esd) return;
  
  Int_t nTracks = (Int_t)esd->GetNumberOfMuonTracks();
  if (nTracks < 1) return;
  
  TList newGhosts;
  newGhosts.SetOwner(kFALSE);
  UInt_t firstGhostId = 0xFFFFFFFF - 1;
  
  // load the current event
  fESDInterface->LoadEvent(*esd, kFALSE);
  
  // remove clusters from ESD (keep digits as they will not change, just eventually not used anymore)
  esd->FindListObject("MuonClusters")->Clear("C");
  
  // modify clusters
  AliMUONVCluster* cluster = 0x0;
  TIter nextCluster(fESDInterface->CreateClusterIterator());
  while ((cluster = static_cast<AliMUONVCluster*>(nextCluster()))) ModifyCluster(*cluster);
  
  // refit the tracks from clusters
  AliMUONVTrackStore* newTrackStore = fRefitter->ReconstructFromClusters();
  
  // loop over the list of ESD tracks
  TClonesArray *esdTracks = (TClonesArray*) esd->FindListObject("MuonTracks");
  for (Int_t iTrack = 0; iTrack < nTracks; iTrack++) {
    
    // get the ESD track
    AliESDMuonTrack* esdTrack = (AliESDMuonTrack*) esdTracks->UncheckedAt(iTrack);
    
    // skip ghost tracks (leave them unchanged)
    if (!esdTrack->ContainTrackerData()) {
      if (esdTrack->GetUniqueID() <= firstGhostId) firstGhostId = esdTrack->GetUniqueID()-1;
      continue;
    }
    
    // Find the corresponding re-fitted MUON track
    AliMUONTrack* newTrack = (AliMUONTrack*) newTrackStore->FindObject(esdTrack->GetUniqueID());
    
    // Find the corresponding locaTrigger if any
    AliMUONLocalTrigger *locTrg = (esdTrack->ContainTriggerData()) ? fESDInterface->FindLocalTrigger(esdTrack->LoCircuit()) : 0x0;
    if (locTrg && locTrg->IsNull()) locTrg = 0x0;
    
    // replace the content of the current ESD track or remove it
    if (newTrack && (!fImproveTracks || newTrack->IsImproved())) {
      
      // eventually remove the trigger part if matching chi2 do not pass the new cut
      if (locTrg && newTrack->GetChi2MatchTrigger() > fSigmaCutForTrigger*fSigmaCutForTrigger) {
	newTrack->SetMatchTrigger(0);
	newTrack->SetLocalTrigger(0,0,0,0,0,0,0);
	newTrack->SetChi2MatchTrigger(0.);
	newTrack->SetHitsPatternInTrigCh(0);
	newGhosts.AddLast(locTrg);
	locTrg = 0x0;
      }
      
      // fill the track info
      Double_t vertex[3] = {esdTrack->GetNonBendingCoor(), esdTrack->GetBendingCoor(), esdTrack->GetZ()};
      AliMUONESDInterface::MUONToESD(*newTrack, *esdTrack, vertex, locTrg);
      
      // add the clusters if not already there
      for (Int_t i = 0; i < newTrack->GetNClusters(); i++) {
	AliMUONVCluster *cl = static_cast<AliMUONTrackParam*>(newTrack->GetTrackParamAtCluster()->UncheckedAt(i))->GetClusterPtr();
	if (esd->FindMuonCluster(cl->GetUniqueID())) continue;
	AliESDMuonCluster *esdCl = esd->NewMuonCluster();
	AliMUONESDInterface::MUONToESD(*cl, *esdCl, kTRUE);
      }
      
    } else {
      
      // keep the trigger part if any
      if (locTrg) newGhosts.AddLast(locTrg);
      
      // remove the track
      esdTracks->Remove(esdTrack);
      
    }
    
  }
  
  // free memory
  delete newTrackStore;
  
  // compress the array of ESD tracks
  esdTracks->Compress();
  
  // add new ghosts if not already there
  TIter nextGhost(&newGhosts);
  AliMUONLocalTrigger *locTrg = 0x0;
  while ((locTrg = static_cast<AliMUONLocalTrigger*>(nextGhost()))) {
    Bool_t alreadyThere = kFALSE;
    for (Int_t iTrack = 0; iTrack < esdTracks->GetEntriesFast(); iTrack++) {
      AliESDMuonTrack* esdTrack = (AliESDMuonTrack*) esdTracks->UncheckedAt(iTrack);
      alreadyThere = (esdTrack->LoCircuit() == locTrg->LoCircuit());
      if (alreadyThere) break;
    }
    if (!alreadyThere) AliMUONESDInterface::MUONToESD(*locTrg, *esd, firstGhostId--);
  }
  
}

//________________________________________________________________________
void AliAnalysisTaskMuonRefit::NotifyRun()
{
  /// load necessary data from OCDB and create the refitter
  
  AliCDBManager* cdbm = AliCDBManager::Instance();
  cdbm->SetDefaultStorage(fDefaultStorage.Data());
  cdbm->SetRun(fCurrentRunNumber);
  
  if (!AliMUONCDB::LoadField()) return;
  
  if (!AliMUONCDB::LoadMapping(kTRUE)) return;
  
  AliMUONRecoParam* recoParam = AliMUONCDB::LoadRecoParam();
  if (!recoParam) return;
  
  if (fImproveTracks) {
    if (fSigmaCut > 0.) recoParam->ImproveTracks(kTRUE, fSigmaCut);
    else recoParam->ImproveTracks(kTRUE);
  } else recoParam->ImproveTracks(kFALSE);
  
  if (fSigmaCutForTrigger > 0.) recoParam->SetSigmaCutForTrigger(fSigmaCutForTrigger);
  else fSigmaCutForTrigger = recoParam->GetSigmaCutForTrigger();
  
  for (Int_t i = 0; i < AliMUONConstants::NTrackingCh(); i++) {
    if (fClusterResNB[i] < 0.) fClusterResNB[i] = recoParam->GetDefaultNonBendingReso(i);
    if (fClusterResB[i] < 0.) fClusterResB[i] = recoParam->GetDefaultBendingReso(i);
  }
  
  if (fReAlign) {
    
    // recover default storage full name (raw:// cannot be used to set specific storage)
    TString defaultStorage(cdbm->GetDefaultStorage()->GetType());
    if (defaultStorage == "alien") defaultStorage += Form("://folder=%s", cdbm->GetDefaultStorage()->GetBaseFolder().Data());
    else defaultStorage += Form("://%s", cdbm->GetDefaultStorage()->GetBaseFolder().Data());
    
    // reset existing geometry/alignment if any
    if (cdbm->GetEntryCache()->Contains("GRP/Geometry/Data")) cdbm->UnloadFromCache("GRP/Geometry/Data");
    if (cdbm->GetEntryCache()->Contains("MUON/Align/Data")) cdbm->UnloadFromCache("MUON/Align/Data");
    if (AliGeomManager::GetGeometry()) AliGeomManager::GetGeometry()->UnlockGeometry();
    
    // get original geometry transformer
    AliGeomManager::LoadGeometry();
    if (!AliGeomManager::GetGeometry()) return;
    if (fOldAlignStorage != "none") {
      if (!fOldAlignStorage.IsNull()) cdbm->SetSpecificStorage("MUON/Align/Data",fOldAlignStorage.Data());
      else cdbm->SetSpecificStorage("MUON/Align/Data",defaultStorage.Data());
      AliGeomManager::ApplyAlignObjsFromCDB("MUON");
    }
    fOldGeoTransformer = new AliMUONGeometryTransformer();
    fOldGeoTransformer->LoadGeometryData();
    
    // get new geometry transformer
    cdbm->UnloadFromCache("GRP/Geometry/Data");
    if (fOldAlignStorage != "none") cdbm->UnloadFromCache("MUON/Align/Data");
    AliGeomManager::GetGeometry()->UnlockGeometry();
    AliGeomManager::LoadGeometry();
    if (!AliGeomManager::GetGeometry()) return;
    if (!fNewAlignStorage.IsNull()) cdbm->SetSpecificStorage("MUON/Align/Data",fNewAlignStorage.Data());
    else cdbm->SetSpecificStorage("MUON/Align/Data",defaultStorage.Data());
    AliGeomManager::ApplyAlignObjsFromCDB("MUON");
    fNewGeoTransformer = new AliMUONGeometryTransformer();
    fNewGeoTransformer->LoadGeometryData();
    
  } else {
    
    // load geometry for track extrapolation to vertex
    if (cdbm->GetEntryCache()->Contains("GRP/Geometry/Data")) cdbm->UnloadFromCache("GRP/Geometry/Data");
    if (AliGeomManager::GetGeometry()) AliGeomManager::GetGeometry()->UnlockGeometry();
    AliGeomManager::LoadGeometry();
    if (!AliGeomManager::GetGeometry()) return;
    
  }
  
  fESDInterface = new AliMUONESDInterface();
  fRefitter = new AliMUONRefitter(recoParam);
  fRefitter->Connect(fESDInterface);
}

//________________________________________________________________________
void AliAnalysisTaskMuonRefit::Terminate(Option_t *)
{
}

//________________________________________________________________________
void AliAnalysisTaskMuonRefit::ModifyCluster(AliMUONVCluster& cl)
{
  /// Reset the cluster resolution to the one given to the task and change
  /// the cluster position according to the new alignment parameters if required
  
  Double_t gX,gY,gZ,lX,lY,lZ;
  
  // change their resolution
  cl.SetErrXY(fClusterResNB[cl.GetChamberId()], fClusterResB[cl.GetChamberId()]);
  
  // change their position
  if (fReAlign) {
    gX = cl.GetX();
    gY = cl.GetY();
    gZ = cl.GetZ();
    fOldGeoTransformer->Global2Local(cl.GetDetElemId(),gX,gY,gZ,lX,lY,lZ);
    fNewGeoTransformer->Local2Global(cl.GetDetElemId(),lX,lY,lZ,gX,gY,gZ);
    cl.SetXYZ(gX,gY,gZ);
  }
  
}

