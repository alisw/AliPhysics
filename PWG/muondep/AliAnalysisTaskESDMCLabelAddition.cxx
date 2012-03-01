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

// STEER includes
#include "AliESDEvent.h"
#include "AliESDMuonTrack.h"
#include "AliLog.h"
#include "AliMCEventHandler.h"
#include "AliGeomManager.h"
#include "AliCDBManager.h"

// ANALYSIS includes
#include "AliAnalysisManager.h"

// MUON includes
#include "AliMUONCDB.h"
#include "AliMUONRecoParam.h"
#include "AliMUONRecoCheck.h"
#include "AliMUONESDInterface.h"
#include "AliMUONTrack.h"
#include "AliMUONVTrackStore.h"
#include "AliMUONTriggerTrack.h"
#include "AliMUONVTriggerTrackStore.h"
#include "AliMUONLocalTrigger.h"

#include "AliAnalysisTaskESDMCLabelAddition.h"

ClassImp(AliAnalysisTaskESDMCLabelAddition)

//----------------------------------------------------------------------
AliAnalysisTaskESDMCLabelAddition::AliAnalysisTaskESDMCLabelAddition():
AliAnalysisTaskSE(),
fDefaultStorage(""),
fSigmaCut(-1.),
fSigmaCutTrig(-1.)
{
  /// Default constructor
}


//----------------------------------------------------------------------
AliAnalysisTaskESDMCLabelAddition::AliAnalysisTaskESDMCLabelAddition(const char* name):
AliAnalysisTaskSE(name),
fDefaultStorage("raw://"),
fSigmaCut(-1.),
fSigmaCutTrig(-1.)
{
  /// Constructor
}


//----------------------------------------------------------------------
void AliAnalysisTaskESDMCLabelAddition::UserCreateOutputObjects()
{
  /// Create output objects
}


//----------------------------------------------------------------------
void AliAnalysisTaskESDMCLabelAddition::NotifyRun()
{
  /// Load OCDB inputs
  
  // set OCDB location
  AliCDBManager* cdbm = AliCDBManager::Instance();
  cdbm->SetDefaultStorage(fDefaultStorage.Data());
  cdbm->SetRun(fCurrentRunNumber);
  
  // load mapping
  if (!AliMUONCDB::LoadMapping()) return;
  
  // load geometry
  if (!AliGeomManager::GetGeometry()) AliGeomManager::LoadGeometry();
  if (!AliGeomManager::GetGeometry()) return;
  
  // load recoParam
  AliMUONRecoParam* recoParam = AliMUONCDB::LoadRecoParam();
  if (!recoParam) {
    fSigmaCut = -1.;
    fSigmaCutTrig = -1.;
    return;
  }
  
  // get sigma cut from recoParam to associate clusters with TrackRefs in case the labels are not used
  fSigmaCut = (recoParam->ImproveTracks()) ? recoParam->GetSigmaCutForImprovement() : recoParam->GetSigmaCutForTracking();
  
  // get sigma cut from recoParam to associate trigger track to triggerable track
  fSigmaCutTrig = recoParam->GetSigmaCutForTrigger();
  
}


//----------------------------------------------------------------------
void AliAnalysisTaskESDMCLabelAddition::UserExec(Option_t */*option*/)
{
  /// Execute analysis for current event				
  
  AliDebug(1, Form("MCLabel Addition: Analysing event # %5d\n",(Int_t) Entry())); 
  
  // make sure necessary information from PCDB have been loaded
  if (fSigmaCut < 0) return;
  
  /// Load ESD event
  AliESDEvent* esd = dynamic_cast<AliESDEvent*>(InputEvent());
  if (!esd) {
    AliError("Cannot get input event");
    return;
  }      
  
  // Load MC event 
  AliMCEventHandler* mcH = static_cast<AliMCEventHandler*>(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
  if ( ! mcH ) {
    AliError ("MCH event handler not found. Nothing done!");
    return;
  }
  
  
  // Get reference tracks
  AliMUONRecoCheck rc(esd,mcH);
  AliMUONVTrackStore* trackRefStore = rc.TrackRefs(-1);
  AliMUONVTriggerTrackStore* triggerTrackRefStore = rc.TriggerableTracks(-1);
  
  // Loop over reconstructed tracks
  AliESDMuonTrack *esdTrack = 0x0;
  Int_t nMuTracks = esd->GetNumberOfMuonTracks();
  for (Int_t nMuTrack = 0; nMuTrack < nMuTracks; ++nMuTrack) {
    
    esdTrack = esd->GetMuonTrack(nMuTrack);
    
    // tracker tracks
    if (esdTrack->ContainTrackerData()) {
      
      // convert ESD track to MUON track (without recomputing track parameters at each clusters)
      AliMUONTrack muonTrack;
      AliMUONESDInterface::ESDToMUON(*esdTrack, muonTrack, kFALSE);
      
      // try to match the reconstructed track with a simulated one
      Int_t nMatchClusters = 0;
      AliMUONTrack* matchedTrackRef = rc.FindCompatibleTrack(muonTrack, *trackRefStore, nMatchClusters, kFALSE, fSigmaCut);
      
      // set the MC label
      if (matchedTrackRef) esdTrack->SetLabel(matchedTrackRef->GetUniqueID());
      else esdTrack->SetLabel(-1);
      
    } else { // ghosts
      
      // Convert ESD track to trigger track
      AliMUONLocalTrigger locTrg;
      AliMUONESDInterface::ESDToMUON(*esdTrack, locTrg);
      AliMUONTriggerTrack trigTrack;
      rc.TriggerToTrack(locTrg, trigTrack);
      
      // try to match the reconstructed track with a simulated one
      AliMUONTriggerTrack* matchedTrigTrackRef = rc.FindCompatibleTrack(trigTrack, *triggerTrackRefStore, fSigmaCutTrig);
      
      // set the MC label
      if (matchedTrigTrackRef) esdTrack->SetLabel(matchedTrigTrackRef->GetUniqueID());
      else esdTrack->SetLabel(-1);
      
    }
    
  }
  
}


//----------------------------------------------------------------------
void AliAnalysisTaskESDMCLabelAddition::Terminate(Option_t */*option*/)
{
  /// Terminate analysis
  AliDebug(2, "Terminate()");
}

