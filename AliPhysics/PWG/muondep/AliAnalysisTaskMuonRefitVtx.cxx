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
#include <TGeoGlobalMagField.h>

// STEER includes
#include "AliESDEvent.h"
#include "AliESDVertex.h"
#include "AliESDMuonTrack.h"
#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliGeomManager.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"

// ANALYSIS includes
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskMuonRefitVtx.h"

// MUON includes
#include "AliMUONCDB.h"
#include "AliMUONTrackParam.h"
#include "AliMUONTrackExtrap.h"
#include "AliMUONESDInterface.h"

/// \ingroup utils

ClassImp(AliAnalysisTaskMuonRefitVtx)

//________________________________________________________________________
AliAnalysisTaskMuonRefitVtx::AliAnalysisTaskMuonRefitVtx() :
AliAnalysisTaskSE(),
fDefaultStorage(""),
fOCDBLoaded(kFALSE),
fUseMeanVtxSPD(kFALSE),
fUseMCVtx(kFALSE),
fUseTrackVtx(kFALSE)
{
  /// Default constructor
  for (Int_t i = 0; i < 3; i++) {
    fVtxPos[i] = 0.;
    fVtxSig[i] = 0.;
    fVtxShift[i] = 0.;
  }
}

//________________________________________________________________________
AliAnalysisTaskMuonRefitVtx::AliAnalysisTaskMuonRefitVtx(const char *name) :
AliAnalysisTaskSE(name),
fDefaultStorage("raw://"),
fOCDBLoaded(kFALSE),
fUseMeanVtxSPD(kFALSE),
fUseMCVtx(kFALSE),
fUseTrackVtx(kFALSE)
{
  /// Constructor
  fBranchNames = "ESD:AliESDRun.,AliESDHeader.,MuonTracks,MuonClusters,MuonPads";
  for (Int_t i = 0; i < 3; i++) {
    fVtxPos[i] = 0.;
    fVtxSig[i] = 0.;
    fVtxShift[i] = 0.;
  }
}

//________________________________________________________________________
AliAnalysisTaskMuonRefitVtx::~AliAnalysisTaskMuonRefitVtx()
{
  /// Destructor
}

//___________________________________________________________________________
void AliAnalysisTaskMuonRefitVtx::UserCreateOutputObjects()
{
}

//________________________________________________________________________
void AliAnalysisTaskMuonRefitVtx::UserExec(Option_t *)
{
  /// Main event loop
  
  // check if OCDB info are properly loaded
  if (!fOCDBLoaded) AliFatal("Problem occur while loading OCDB objects");
  
  AliESDEvent* esd = dynamic_cast<AliESDEvent*>(InputEvent());
  if (!esd) return;
  
  LoadBranches();
  
  // load the MC vertex position
  if (fUseMCVtx) {
    AliMCEventHandler* mcH = static_cast<AliMCEventHandler*>(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
    if ( ! mcH ) {
      AliError("MCH event handler not found. Nothing done!");
      return;
    }
    mcH->MCEvent()->GetPrimaryVertex()->GetXYZ(fVtxPos);
    for (Int_t i = 0; i < 3; i++) fVtxPos[i] += fVtxShift[i];
    fVtxSig[0] = fVtxSig[1] = fVtxSig[2] = 0.;
  }
  
  // loop over the list of ESD tracks
  AliMUONTrackParam trackParam;
  Int_t nTracks = (Int_t)esd->GetNumberOfMuonTracks();
  for (Int_t iTrack = 0; iTrack < nTracks; iTrack++) {
    
    // get the ESD track
    AliESDMuonTrack* esdTrack = esd->GetMuonTrack(iTrack);
    if (!esdTrack->ContainTrackerData()) continue;
    
    // keep the same vertex (from track position)
    if (fUseTrackVtx) {
      fVtxPos[0] = esdTrack->GetNonBendingCoor() + fVtxShift[0];
      fVtxPos[1] = esdTrack->GetBendingCoor() + fVtxShift[1];
      fVtxPos[2] = esdTrack->GetZ() + fVtxShift[2];
      fVtxSig[0] = fVtxSig[1] = fVtxSig[2] = 0.;
    }
    
    // extrapolate to the new vertex position
    AliMUONESDInterface::GetParamAtFirstCluster(*esdTrack, trackParam);
    AliMUONTrackExtrap::ExtrapToVertex(&trackParam, fVtxPos[0], fVtxPos[1], fVtxPos[2], fVtxSig[0], fVtxSig[1]);
    AliMUONESDInterface::SetParamAtVertex(trackParam, *esdTrack);
    
  }
  
}

//________________________________________________________________________
void AliAnalysisTaskMuonRefitVtx::NotifyRun()
{
  /// load necessary data from OCDB and create the refitter
  
  if (fOCDBLoaded) return;
  
  // set OCDB location
  AliCDBManager* cdbm = AliCDBManager::Instance();
  if (cdbm->IsDefaultStorageSet()) printf("RefitVtxTask: CDB default storage already set!\n");
  else cdbm->SetDefaultStorage(fDefaultStorage.Data());
  if (cdbm->GetRun() > -1) printf("RefitVtxTask: run number already set!\n");
  else cdbm->SetRun(fCurrentRunNumber);
  
  // load magnetic field for track extrapolation
  if (!TGeoGlobalMagField::Instance()->GetField()) {
    if (!AliMUONCDB::LoadField()) return;
  }
  AliMUONTrackExtrap::SetField();
  
  // load geometry for track extrapolation to vertex
  if (!AliGeomManager::GetGeometry()) {
    AliGeomManager::LoadGeometry();
    if (!AliGeomManager::GetGeometry()) return;  
    if (!AliGeomManager::ApplyAlignObjsFromCDB("MUON")) return;
  }
  
  // load the mean SPD vertex position
  if (fUseMeanVtxSPD) {
    AliCDBEntry* entry = AliCDBManager::Instance()->Get("GRP/Calib/MeanVertexSPD");
    if (!entry) return;
    AliESDVertex* vertex = dynamic_cast<AliESDVertex*>(entry->GetObject());
    if (!vertex) return;
    if (vertex->GetXRes() < 2.8) {
      vertex->GetXYZ(fVtxPos);
      vertex->GetSigmaXYZ(fVtxSig);
      printf("SPD Vertex position from OCDB --> x = %g ± %g, y = %g ± %g, z = %g ± %g\n",
	     fVtxPos[0], fVtxSig[0], fVtxPos[1], fVtxSig[1], fVtxPos[2], fVtxSig[2]);
    } else {
      fVtxPos[0] = fVtxPos[1] = fVtxPos[2] = 0.;
      fVtxSig[0] = fVtxSig[1] = fVtxSig[2] = 0.;
      printf("SPD Vertex position not valid --> x = 0 ± 0, y = 0 ± 0, z = 0 ± 0\n");
    }
    for (Int_t i = 0; i < 3; i++) fVtxPos[i] += fVtxShift[i];
  }
  
  fOCDBLoaded = kTRUE;
  
}

//________________________________________________________________________
void AliAnalysisTaskMuonRefitVtx::Terminate(Option_t *)
{
}

