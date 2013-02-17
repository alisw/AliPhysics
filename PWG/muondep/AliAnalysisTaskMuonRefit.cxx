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
#include <TGeoGlobalMagField.h>
#include <TPRegexp.h>

// STEER includes
#include "AliESDEvent.h"
#include "AliESDMuonTrack.h"
#include "AliCDBManager.h"
#include "AliGRPManager.h"
#include "AliGRPObject.h"
#include "AliCDBStorage.h"
#include "AliGeomManager.h"
#include "AliMagF.h"

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
#include "AliMUONTrackExtrap.h"
#include "AliMUONVCluster.h"
#include "AliMUONVDigit.h"
#include "AliMUONVDigitStore.h"
#include "AliMUONVTrackStore.h"
#include "AliMUONLocalTrigger.h"
#include "AliMUONTriggerTrack.h"
#include "AliMUONVTriggerStore.h"
#include "AliMUONVTriggerTrackStore.h"
#include "AliMUONGeometryTransformer.h"
#include "AliMUONTriggerCircuit.h"
#include "AliMUONTrackHitPattern.h"
#include "AliMUONVTrackReconstructor.h"

// MUON mapping includes
#include "AliMpSegmentation.h"
#include "AliMpVSegmentation.h"
#include "AliMpConstants.h"
#include "AliMpDDLStore.h"
#include "AliMpManuStore.h"
#include "AliMpPad.h"
#include "AliMpDetElement.h"
#include "AliMpCathodType.h"

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
fField(""),
fRemoveMonoCathCl(kFALSE),
fCheckAllPads(kFALSE),
fTagBadTracks(kFALSE),
fKeepOldParam(kFALSE),
fTriggerCircuit(NULL),
fESDInterface(NULL),
fRefitter(NULL),
fTrackHitPattern(NULL)
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
fField(""),
fRemoveMonoCathCl(kFALSE),
fCheckAllPads(kFALSE),
fTagBadTracks(kFALSE),
fKeepOldParam(kFALSE),
fTriggerCircuit(NULL),
fESDInterface(NULL),
fRefitter(NULL),
fTrackHitPattern(NULL)
{
  /// Constructor
  fBranchNames = "ESD:AliESDRun.,AliESDHeader.,MuonTracks,MuonClusters,MuonPads";
  for (Int_t i = 0; i < 10; i++) ResetClusterResolution(i, -1., -1.);
}

//________________________________________________________________________
AliAnalysisTaskMuonRefit::~AliAnalysisTaskMuonRefit()
{
  /// Destructor
  SafeDelete(fOldGeoTransformer);
  SafeDelete(fNewGeoTransformer);
  SafeDelete(fTriggerCircuit);
  SafeDelete(fESDInterface);
  SafeDelete(fRefitter);
  SafeDelete(fTrackHitPattern);
}

//___________________________________________________________________________
void AliAnalysisTaskMuonRefit::UserCreateOutputObjects()
{
}

//________________________________________________________________________
void AliAnalysisTaskMuonRefit::UserExec(Option_t *)
{
  /// Main event loop
  
  // check if refitter properly created (i.e. OCDB data properly loaded)
  if (!fRefitter) AliFatal("Problem occur while loading OCDB objects");
  
  AliESDEvent* esd = dynamic_cast<AliESDEvent*>(InputEvent());
  if (!esd) return;
  
  LoadBranches();
  
  Int_t nTracks = (Int_t)esd->GetNumberOfMuonTracks();
  if (nTracks < 1) return;
  
  TList oldGhosts;
  oldGhosts.SetOwner(kFALSE);
  UInt_t firstGhostId = 0xFFFFFFFF - 1;
  
  TList esdTracksToRemove;
  esdTracksToRemove.SetOwner(kFALSE);
  
  // load the current event
  fESDInterface->LoadEvent(*esd, kFALSE);
  
  // remove clusters from ESD (keep digits as they will not change, just eventually not used anymore)
  esd->FindListObject("MuonClusters")->Clear("C");
  
  // modify clusters
  Int_t firstClusterIndex = 0;
  AliMUONVCluster* cluster = 0x0;
  TIter nextCluster(fESDInterface->CreateClusterIterator());
  while ((cluster = static_cast<AliMUONVCluster*>(nextCluster()))) {
    Int_t clIndex = cluster->GetClusterIndex(cluster->GetUniqueID());
    if (clIndex >= firstClusterIndex) firstClusterIndex = clIndex+1;
    ModifyCluster(*cluster);
  }
  
  // to do not mix with old clusters in case we re-clusterize
  fRefitter->SetFirstClusterIndex(firstClusterIndex);
  
  // refit the tracks from clusters
  AliMUONVTrackStore* newTrackStore = fRefitter->ReconstructFromClusters();
  
  // reset the trigger part
  TIter nextNewTrack(newTrackStore->CreateIterator());
  AliMUONTrack *newTrack = 0x0;
  while ((newTrack = static_cast<AliMUONTrack*>(nextNewTrack()))) {
    newTrack->SetMatchTrigger(0);
    newTrack->SetLocalTrigger(0,0,0,0,0,0,0);
    newTrack->SetChi2MatchTrigger(0.);
    newTrack->SetHitsPatternInTrigCh(0);
    newTrack->SetHitsPatternInTrigChTrk(0);
  }
  
  // reconstruct trigger tracks
  AliMUONVTriggerStore *trigStore = fESDInterface->GetTriggers();
  AliMUONVTriggerTrackStore *trigTrackStore = AliMUONESDInterface::NewTriggerTrackStore();
  if (!trigTrackStore) return;
  AliMUONVTrackReconstructor* tracker = AliMUONESDInterface::GetTracker();
  tracker->EventReconstructTrigger(*fTriggerCircuit, *trigStore, *trigTrackStore);
  
  // recover the hit pattern for all trigger tracks
  TClonesArray *esdTracks = (TClonesArray*) esd->FindListObject("MuonTracks");
  for (Int_t iTrack = 0; iTrack < nTracks; iTrack++) {
    AliESDMuonTrack* esdTrack = (AliESDMuonTrack*) esdTracks->UncheckedAt(iTrack);
    if (!esdTrack->ContainTriggerData()) continue;
    // use the UniqueID of the local trigger in case several tracks match the same trigger
    AliMUONLocalTrigger *locTrg = fESDInterface->FindLocalTrigger(esdTrack->LoCircuit());
    AliMUONTriggerTrack *trigTrack = (AliMUONTriggerTrack*) trigTrackStore->FindObject(locTrg->GetUniqueID());
    trigTrack->SetHitsPatternInTrigCh(esdTrack->GetHitsPatternInTrigCh());
  }
  
  // match tracker/trigger tracks
  const Int_t kFirstTrigCh = AliMUONConstants::NTrackingCh();
  const Float_t kZFilterOut = AliMUONConstants::MuonFilterZEnd();
  const Float_t kFilterThickness = kZFilterOut-AliMUONConstants::MuonFilterZBeg();
  nextNewTrack.Reset();
  while ((newTrack = static_cast<AliMUONTrack*>(nextNewTrack()))) {
    AliMUONTrackParam trackParam(*((AliMUONTrackParam*) (newTrack->GetTrackParamAtCluster()->Last())));
    AliMUONTrackExtrap::ExtrapToZCov(&trackParam, kZFilterOut);
    AliMUONTrackExtrap::AddMCSEffect(&trackParam, kFilterThickness, AliMUONConstants::MuonFilterX0());
    AliMUONTrackExtrap::ExtrapToZCov(&trackParam, AliMUONConstants::DefaultChamberZ(kFirstTrigCh));
    AliMUONTriggerTrack *matchedTriggerTrack = fTrackHitPattern->MatchTriggerTrack(newTrack, trackParam, *trigTrackStore, *trigStore);
    if ( matchedTriggerTrack ) newTrack->SetHitsPatternInTrigCh(matchedTriggerTrack->GetHitsPatternInTrigCh());
  }
  
  // loop over the list of ESD tracks
  for (Int_t iTrack = 0; iTrack < nTracks; iTrack++) {
    
    // get the ESD track
    AliESDMuonTrack* esdTrack = (AliESDMuonTrack*) esdTracks->UncheckedAt(iTrack);
    
    // keep the memory of old ghosts to concerve their Id if their are still ghosts after refitting
    // and remember to remove them before adding the new ghosts
    if (!esdTrack->ContainTrackerData()) {
      AliMUONLocalTrigger *locTrg = (esdTrack->ContainTriggerData()) ? fESDInterface->FindLocalTrigger(esdTrack->LoCircuit()) : 0x0;
      if (locTrg && !locTrg->IsNull()) oldGhosts.AddLast(locTrg);
      if (esdTrack->GetUniqueID() <= firstGhostId) firstGhostId = esdTrack->GetUniqueID()-1;
      esdTracksToRemove.AddLast(esdTrack);
      continue;
    }
    
    // Find the corresponding re-fitted MUON track
    newTrack = (AliMUONTrack*) newTrackStore->FindObject(esdTrack->GetUniqueID());
    
    // replace the content of the current ESD track or remove it
    if (newTrack && (!fImproveTracks || newTrack->IsImproved())) {
      
      // eventually tag the tracks which do not match the trigger while they were before
      Bool_t noLongerMatched = (fTagBadTracks && esdTrack->ContainTriggerData() && newTrack->GetMatchTrigger() <= 0);
      
      // find the corresponding trigger info and eventually restore the old one
      AliMUONLocalTrigger *locTrg = 0x0;
      if (newTrack->GetMatchTrigger() > 0) locTrg = fESDInterface->FindLocalTrigger(newTrack->LoCircuit());
      else if (noLongerMatched && !fKeepOldParam) {
	locTrg = fESDInterface->FindLocalTrigger(esdTrack->LoCircuit());
	newTrack->SetMatchTrigger(esdTrack->GetMatchTrigger());
	newTrack->SetChi2MatchTrigger(esdTrack->GetChi2MatchTrigger());
	newTrack->SetHitsPatternInTrigCh(esdTrack->GetHitsPatternInTrigCh());
	newTrack->SetHitsPatternInTrigChTrk(esdTrack->GetHitsPatternInTrigChTrk());
	newTrack->SetLocalTrigger(esdTrack->LoCircuit(), esdTrack->LoStripX(), esdTrack->LoStripY(),
				  esdTrack->LoDev(), esdTrack->LoLpt(), esdTrack->LoHpt(),
				  esdTrack->GetTriggerWithoutChamber());
      }
      if (locTrg && locTrg->IsNull()) locTrg = 0x0;
      
      // fill the new track/cluster info
      if (!(noLongerMatched && fKeepOldParam)) {
	
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
	
	// restore the old clusters if not already done
	for (Int_t i = 0; i < esdTrack->GetNClusters(); i++) {
	  UInt_t clId = esdTrack->GetClusterId(i);
	  if (esd->FindMuonCluster(clId)) continue;
	  AliMUONVCluster *cl = fESDInterface->FindCluster(esdTrack->GetUniqueID(), clId);
	  if (!cl) continue;
	  AliESDMuonCluster *esdCl = esd->NewMuonCluster();
	  AliMUONESDInterface::MUONToESD(*cl, *esdCl, kTRUE);
	}
	
      }
      
      // eventually tag the trigger part as bad
      esdTrack->ResetBit(BIT(20));
      if (noLongerMatched) esdTrack->SetBit(BIT(21));
      else esdTrack->ResetBit(BIT(21));
      
    } else {
      
      // simply tag it or remember to remove that track (cannot remove it now as it will create a hole which
      // will eventually produce a crash in parallel track loop (AliESDEvent::MoveMuonObjects()))
      if (fTagBadTracks) {
	esdTrack->SetBit(BIT(20));
	if (esdTrack->ContainTriggerData()) esdTrack->SetBit(BIT(21));
	else esdTrack->ResetBit(BIT(21));
      }
      else esdTracksToRemove.AddLast(esdTrack);
      
    }
    
  }
  
  // free memory
  delete newTrackStore;
  
  // remove tracks to remove and compress the array of ESD tracks
  TIter nextTrackToRemove(&esdTracksToRemove);
  AliESDMuonTrack* esdTrack = 0x0;
  while ((esdTrack = static_cast<AliESDMuonTrack*>(nextTrackToRemove()))) esdTracks->Remove(esdTrack);
  esdTracks->Compress();
  
  // add new ghosts (ignoring tracks marked as bad or no longer matched)
  TIter nextGhost(trigStore->CreateIterator());
  AliMUONLocalTrigger *locTrg = 0x0;
  while ((locTrg = static_cast<AliMUONLocalTrigger*>(nextGhost()))) {
    Bool_t alreadyThere = kFALSE;
    for (Int_t iTrack = 0; iTrack < esdTracks->GetEntriesFast(); iTrack++) {
      esdTrack = (AliESDMuonTrack*) esdTracks->UncheckedAt(iTrack);
      alreadyThere = (esdTrack->LoCircuit() == locTrg->LoCircuit() && !esdTrack->TestBit(BIT(21)));
      if (alreadyThere) break;
    }
    if (!alreadyThere) {
      UInt_t ghostId = oldGhosts.Contains(locTrg) ? locTrg->GetUniqueID() : firstGhostId--;
      AliMUONTriggerTrack *trigTrack = (AliMUONTriggerTrack*) trigTrackStore->FindObject(locTrg->GetUniqueID());
      AliMUONESDInterface::MUONToESD(*locTrg, *esd, ghostId, trigTrack);
    }
  }
  
}

//________________________________________________________________________
void AliAnalysisTaskMuonRefit::NotifyRun()
{
  /// load necessary data from OCDB and create the refitter
  
  // do it only once
  if (fRefitter) return;
  
  // set OCDB location
  AliCDBManager* cdbm = AliCDBManager::Instance();
  if (cdbm->IsDefaultStorageSet()) printf("MCLabelAddition: CDB default storage already set!\n");
  else {
    cdbm->SetDefaultStorage(fDefaultStorage.Data());
    if (fOldAlignStorage != "none" && !fOldAlignStorage.IsNull())
      cdbm->SetSpecificStorage("MUON/Align/Data",fOldAlignStorage.Data());
  }
  if (cdbm->GetRun() > -1) printf("MCLabelAddition: run number already set!\n");
  else cdbm->SetRun(fCurrentRunNumber);
  
  // load magnetic field or create it for track extrapolation
  if (!TGeoGlobalMagField::Instance()->GetField()) {
    if (!fField.IsNull()) {
      if (!SetMagField()) return;
    } else if (!AliMUONCDB::LoadField()) return;
  }
  
  // load mapping
  if (!AliMpDDLStore::Instance(kFALSE) || !AliMpManuStore::Instance(kFALSE)) {
    if (!AliMUONCDB::LoadMapping()) return;
  }
  
  // load recoParam for refitting
  AliMUONRecoParam* recoParam = AliMUONCDB::LoadRecoParam();
  if (!recoParam) return;
  if (!AliMUONESDInterface::GetTracker()) AliMUONESDInterface::ResetTracker(recoParam);
  
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
  
  // load geometry for track extrapolation to vertex and mono-cathod cluster finding
  if (!AliGeomManager::GetGeometry()) {
    
    if (fReAlign) {
      
      // get original geometry transformer
      AliGeomManager::LoadGeometry();
      if (!AliGeomManager::GetGeometry()) return;
      if (fOldAlignStorage != "none" && !AliGeomManager::ApplyAlignObjsFromCDB("MUON")) return;
      fOldGeoTransformer = new AliMUONGeometryTransformer();
      fOldGeoTransformer->LoadGeometryData();
      
      // load the new geometry
      if (cdbm->GetEntryCache()->Contains("GRP/Geometry/Data")) cdbm->UnloadFromCache("GRP/Geometry/Data");
      AliGeomManager::GetGeometry()->UnlockGeometry();
      AliGeomManager::LoadGeometry();
      if (!AliGeomManager::GetGeometry()) return;
      if (cdbm->GetEntryCache()->Contains("MUON/Align/Data")) cdbm->UnloadFromCache("MUON/Align/Data");
      if (!fNewAlignStorage.IsNull()) cdbm->SetSpecificStorage("MUON/Align/Data",fNewAlignStorage.Data());
      else {
	// recover default storage full name (raw:// cannot be used to set specific storage)
	TString defaultStorage(cdbm->GetDefaultStorage()->GetType());
	if (defaultStorage == "alien") defaultStorage += Form("://folder=%s", cdbm->GetDefaultStorage()->GetBaseFolder().Data());
	else defaultStorage += Form("://%s", cdbm->GetDefaultStorage()->GetBaseFolder().Data());
	cdbm->SetSpecificStorage("MUON/Align/Data",defaultStorage.Data());
      }
      if (!AliGeomManager::ApplyAlignObjsFromCDB("MUON")) return;
      
    } else {
      
      AliGeomManager::LoadGeometry();
      if (!AliGeomManager::GetGeometry()) return;
      if (!AliGeomManager::ApplyAlignObjsFromCDB("MUON")) return;
      
    }
    
  } else fReAlign = kFALSE; // disable the realignment if the geometry was already loaded
  
  fNewGeoTransformer = new AliMUONGeometryTransformer();
  fNewGeoTransformer->LoadGeometryData();
  
  fTriggerCircuit = new AliMUONTriggerCircuit(fNewGeoTransformer);
  
  fESDInterface = new AliMUONESDInterface();
  
  fRefitter = new AliMUONRefitter(recoParam);
  fRefitter->Connect(fESDInterface);
  
  AliMUONVDigitStore *digitStore = AliMUONESDInterface::NewDigitStore();
  fTrackHitPattern = new AliMUONTrackHitPattern(recoParam, *fNewGeoTransformer, *digitStore, 0x0);
  delete digitStore; // assume the digitStore is never used when we use AliMUONTrackHitPattern
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
  
  // "remove" mono-cathod clusters on stations 3-4-5 if required
  // (to be done after moving clusters to the new position)
  if (fRemoveMonoCathCl && cl.GetChamberId() > 3) {
    Bool_t hasBending, hasNonBending;
    if (fCheckAllPads) CheckPads(&cl, hasBending, hasNonBending);
    else CheckPadsBelow(&cl, hasBending, hasNonBending);
    if (!hasNonBending) cl.SetErrXY(10., cl.GetErrY());
    if (!hasBending) cl.SetErrXY(cl.GetErrX(), 10.);
  }
  
}

//________________________________________________________________________
void AliAnalysisTaskMuonRefit::CheckPads(AliMUONVCluster *cl, Bool_t &hasBending, Bool_t &hasNonBending) const
{
  /// Check that this cluster contains pads on both cathods
  
  // reset
  hasBending = kFALSE;
  hasNonBending = kFALSE;
  
  // loop over digits contained in the cluster
  for (Int_t iDigit = 0; iDigit < cl->GetNDigits(); iDigit++) {
    
    Int_t manuId = AliMUONVDigit::ManuId(cl->GetDigitId(iDigit));
    
    // check the location of the manu the digit belongs to
    if (manuId > 0) {
      if (manuId & AliMpConstants::ManuMask(AliMp::kNonBendingPlane)) hasNonBending = kTRUE;
      else hasBending = kTRUE;
    }
    
  }
  
}

//________________________________________________________________________
void AliAnalysisTaskMuonRefit::CheckPadsBelow(AliMUONVCluster *cl, Bool_t &hasBending, Bool_t &hasNonBending) const
{
  /// Check that this cluster contains pads on both cathods just under its position
  
  // reset
  hasBending = kFALSE;
  hasNonBending = kFALSE;
  
  // get the cathod corresponding to the bending/non-bending plane
  Int_t deId = cl->GetDetElemId();
  AliMpDetElement* de = AliMpDDLStore::Instance()->GetDetElement(deId, kFALSE);
  if (!de) return;
  AliMp::CathodType cath1 = de->GetCathodType(AliMp::kBendingPlane); 
  AliMp::CathodType cath2 = de->GetCathodType(AliMp::kNonBendingPlane); 
  
  // get the corresponding segmentation
  const AliMpVSegmentation* seg1 = AliMpSegmentation::Instance()->GetMpSegmentation(deId, cath1);
  const AliMpVSegmentation* seg2 = AliMpSegmentation::Instance()->GetMpSegmentation(deId, cath2);
  if (!seg1 || !seg2) return;
  
  // get local coordinate of the cluster
  Double_t lX,lY,lZ;
  Double_t gX = cl->GetX();
  Double_t gY = cl->GetY();
  Double_t gZ = cl->GetZ();
  fNewGeoTransformer->Global2Local(deId,gX,gY,gZ,lX,lY,lZ);
  
  // find pads below the cluster
  AliMpPad pad1 = seg1->PadByPosition(lX, lY, kFALSE);
  AliMpPad pad2 = seg2->PadByPosition(lX, lY, kFALSE);
  
  // build their ID if pads are valid
  UInt_t padId1 = (pad1.IsValid()) ? AliMUONVDigit::BuildUniqueID(deId, pad1.GetManuId(), pad1.GetManuChannel(), cath1) : 0;
  UInt_t padId2 = (pad2.IsValid()) ? AliMUONVDigit::BuildUniqueID(deId, pad2.GetManuId(), pad2.GetManuChannel(), cath2) : 0;
  
  // check if the cluster contains these pads 
  for (Int_t iDigit = 0; iDigit < cl->GetNDigits(); iDigit++) {
    
    UInt_t digitId = cl->GetDigitId(iDigit);
    
    if (digitId == padId1) {
      
      hasBending = kTRUE;
      if (hasNonBending) break;
      
    } else if (digitId == padId2) {
      
      hasNonBending = kTRUE;
      if (hasBending) break;
      
    }
    
  }
  
}

//________________________________________________________________________
Bool_t AliAnalysisTaskMuonRefit::SetMagField() const
{
  // Dealing with the magnetic field map assuming nominal current
  // Construct the mag field map from the data in GRP and the given field map
  // Set the global mag field instance
  
  if (TGeoGlobalMagField::Instance()->IsLocked()) delete TGeoGlobalMagField::Instance();
  
  AliGRPManager grpMan;
  if (!grpMan.ReadGRPEntry()) {
    AliError("failed to load GRP Data from OCDB");
    return kFALSE;
  }
  
  const AliGRPObject *grpData = grpMan.GetGRPData();
  if (!grpData) {
    AliError("GRP Data is not loaded");
    return kFALSE;
  }
  
  Char_t l3Polarity = grpData->GetL3Polarity();
  if (l3Polarity == AliGRPObject::GetInvalidChar()) {
    AliError("GRP/GRP/Data entry:  missing value for the L3 polarity !");
    return kFALSE;
  }
  
  Char_t diPolarity = grpData->GetDipolePolarity();
  if (diPolarity == AliGRPObject::GetInvalidChar()) {
    AliError("GRP/GRP/Data entry:  missing value for the dipole polarity !");
    return kFALSE;
  }
  
  TString beamType = grpData->GetBeamType();
  if (beamType==AliGRPObject::GetInvalidString()) {
    AliError("GRP/GRP/Data entry:  missing value for the beam type !");
    return kFALSE;
  }
  
  Float_t beamEnergy = grpData->GetBeamEnergy();
  if (beamEnergy==AliGRPObject::GetInvalidFloat()) {
    AliError("GRP/GRP/Data entry:  missing value for the beam energy !");
    return kFALSE;
  }
  
  AliMagF::BeamType_t btype = AliMagF::kNoBeamField;
  TString btypestr = beamType;
  btypestr.ToLower();
  TPRegexp protonBeam("(proton|p)\\s*-?\\s*\\1");
  TPRegexp ionBeam("(lead|pb|ion|a|A)\\s*-?\\s*\\1");
  TPRegexp protonionBeam("(proton|p)\\s*-?\\s*(lead|pb|ion|a|A)");
  TPRegexp ionprotonBeam("(lead|pb|ion|a|A)\\s*-?\\s*(proton|p)");
  if (btypestr.Contains(ionBeam)) btype = AliMagF::kBeamTypeAA;
  else if (btypestr.Contains(protonBeam)) btype = AliMagF::kBeamTypepp;
  else if (btypestr.Contains(protonionBeam)) btype = AliMagF::kBeamTypepA;
  else if (btypestr.Contains(ionprotonBeam)) btype = AliMagF::kBeamTypeAp;
  else AliInfoGeneral("AliMagF",Form("Assume no LHC magnet field for the beam type %s !",beamType.Data()));
  
  TGeoGlobalMagField::Instance()->SetField(new AliMagF("fld", "fld", l3Polarity ? -1:1, diPolarity ? -1:1,
						       AliMagF::k5kG, btype, beamEnergy, 2, 15., fField.Data()));
  TGeoGlobalMagField::Instance()->Lock();
  AliInfo("Running with the B field constructed out of GRP and custom field map !");
  
  return kTRUE;
  
}

