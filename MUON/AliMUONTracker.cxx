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

//-----------------------------------------------------------------------------
/// \class AliMUONTracker
///
/// Steering class for use in global tracking framework;
/// reconstruct tracks from recpoints
///
/// Actual tracking is performed by some AliMUONVTrackReconstructor children
///
/// \author Christian Finck and Laurent Aphecetche, SUBATECH Nantes
//-----------------------------------------------------------------------------

#include "AliMUONTracker.h"

#include "AliMUONReconstructor.h"
#include "AliMUONRecoParam.h"
#include "AliMUONTrack.h"
#include "AliMUONTrackExtrap.h"
#include "AliMUONTrackHitPattern.h"
#include "AliMUONTrackParam.h"
#include "AliMUONHitForRec.h"
#include "AliMUONTrackReconstructor.h"
#include "AliMUONTrackReconstructorK.h"
#include "AliMUONTrackStoreV1.h"
#include "AliMUONTriggerChamberEff.h"
#include "AliMUONTriggerTrackStoreV1.h"
#include "AliMUONVClusterStore.h"
#include "AliMUONVTriggerStore.h"

#include "AliESDEvent.h"
#include "AliESDMuonTrack.h"
#include "AliESDVertex.h"
#include "AliLog.h"
#include "AliCodeTimer.h"

#include <Riostream.h>
#include <TTree.h>

/// \cond CLASSIMP
ClassImp(AliMUONTracker)
/// \endcond


//_____________________________________________________________________________
AliMUONTracker::AliMUONTracker(const AliMUONDigitMaker* digitMaker,
                               const AliMUONGeometryTransformer* transformer,
                               const AliMUONTriggerCircuit* triggerCircuit,
                               AliMUONTriggerChamberEff* chamberEff)
: AliTracker(),
  fDigitMaker(digitMaker), // not owner
  fTransformer(transformer), // not owner
  fTriggerCircuit(triggerCircuit), // not owner
  fTrigChamberEff(chamberEff), // not owner
  fTrackHitPatternMaker(0x0),
  fTrackReco(0x0),
  fClusterStore(0x0),
  fTriggerStore(0x0)
{
  /// constructor
  if (fTransformer && fDigitMaker)
    fTrackHitPatternMaker = new AliMUONTrackHitPattern(*fTransformer,*fDigitMaker);
}

//_____________________________________________________________________________
AliMUONTracker::~AliMUONTracker()
{
  /// dtor
  delete fTrackReco;
  delete fTrackHitPatternMaker;
  delete fClusterStore;
  delete fTriggerStore;
}

//_____________________________________________________________________________
Int_t AliMUONTracker::LoadClusters(TTree* clustersTree)
{
  /// Load clusterStore and triggerStore from clustersTree
  delete fClusterStore;
  delete fTriggerStore;

  if ( ! clustersTree ) {
    AliFatal("No clustersTree");
    return 1;
  }

  fClusterStore = AliMUONVClusterStore::Create(*clustersTree);
  fTriggerStore = AliMUONVTriggerStore::Create(*clustersTree);
  
  if (!fClusterStore)
  {
    AliError("Could not get clusterStore");
    return 1;
  }
  if (!fTriggerStore)
  {
    AliError("Could not get triggerStore");
    return 2;
  }
  
  fClusterStore->Connect(*clustersTree,kFALSE);
  fTriggerStore->Connect(*clustersTree,kFALSE);
  
  clustersTree->GetEvent(0);

  return 0;
}

//_____________________________________________________________________________
Int_t AliMUONTracker::Clusters2Tracks(AliESDEvent* esd)
{
  /// Performs the tracking and store the resulting tracks in the ESD
  AliDebug(1,"");
  AliCodeTimerAuto("")
  
  if (!fTrackReco) CreateTrackReconstructor();
  
  // if the required tracking mode does not exist
  if  (!fTrackReco) return 1;
  
  if (!fClusterStore) {
    AliError("ClusterStore is NULL");
    return 2;
  }
  
  if (!fTriggerStore) {
    AliError("TriggerStore is NULL");
    return 3;
  }
  
  // Make tracker tracks
  AliMUONVTrackStore* trackStore = new AliMUONTrackStoreV1;
  fTrackReco->EventReconstruct(*fClusterStore,*trackStore);
  
  // Make trigger tracks
  AliMUONVTriggerTrackStore* triggerTrackStore(0x0);
  if ( fTriggerCircuit ) {
    triggerTrackStore = new AliMUONTriggerTrackStoreV1;
    fTrackReco->EventReconstructTrigger(*fTriggerCircuit,*fTriggerStore,*triggerTrackStore);
  }

  // Match tracker/trigger tracks
  if ( triggerTrackStore && fTrackHitPatternMaker ) {
    fTrackReco->ValidateTracksWithTrigger(*trackStore,*triggerTrackStore,*fTriggerStore,*fTrackHitPatternMaker);
  }
  
  // Compute trigger chamber efficiency
  if( triggerTrackStore && fTrigChamberEff){
      AliCodeTimerStart("EventChamberEff");
      fTrigChamberEff->EventChamberEff(*fTriggerStore,*triggerTrackStore,*trackStore);
      AliCodeTimerStop("EventChamberEff");
  }
  
  // Fill ESD
  FillESD(*trackStore,esd);
  
  // cleanup
  delete trackStore;
  delete triggerTrackStore;
  
  return 0;
}

//_____________________________________________________________________________
void AliMUONTracker::FillESD(AliMUONVTrackStore& trackStore, AliESDEvent* esd) const
{
  /// Fill the ESD from the trackStore
  AliDebug(1,"");
  AliCodeTimerAuto("")
  
  // Get vertex 
  Double_t vertex[3] = {0};
  const AliESDVertex* esdVert = esd->GetVertex(); 
  if (esdVert->GetNContributors()) 
  {
    esdVert->GetXYZ(vertex);
    AliDebug(1,Form("found vertex (%e,%e,%e)",vertex[0],vertex[1],vertex[2]));
  }
  
  // setting ESD MUON class
  AliESDMuonTrack esdTrack;
  
  AliMUONTrack* track;
  TIter next(trackStore.CreateIterator());
  
  while ( ( track = static_cast<AliMUONTrack*>(next()) ) )
  {
    AliMUONTrackParam* trackParam = static_cast<AliMUONTrackParam*>((track->GetTrackParamAtHit())->First());
    AliMUONTrackParam trackParamAtVtx(*trackParam);
    
    /// Extrapolate to vertex (which is set to (0,0,0) if not available, see above)
    AliMUONTrackExtrap::ExtrapToVertex(&trackParamAtVtx, vertex[0],vertex[1],vertex[2]);
    
    // setting data member of ESD MUON
    
    // at first station
    trackParam->SetParamForUncorrected(esdTrack);
    trackParam->SetCovFor(esdTrack);
    // at vertex
    trackParamAtVtx.SetParamFor(esdTrack);
    // global info
    esdTrack.SetChi2(track->GetFitFMin());
    esdTrack.SetNHit(track->GetNTrackHits());
    esdTrack.SetLocalTrigger(track->GetLocalTrigger());
    esdTrack.SetChi2MatchTrigger(track->GetChi2MatchTrigger());
    esdTrack.SetHitsPatternInTrigCh(track->GetHitsPatternInTrigCh());
    // muon cluster map
    AliMUONHitForRec* cluster = static_cast<AliMUONHitForRec*>((track->GetHitForRecAtHit())->First());
    while (cluster) {
      esdTrack.AddInMuonClusterMap(cluster->GetChamberNumber());
      cluster = static_cast<AliMUONHitForRec*>((track->GetHitForRecAtHit())->After(cluster));
    }
    
    // storing ESD MUON Track into ESD Event 
    esd->AddMuonTrack(&esdTrack);
  } // end of loop on tracks
}

//_____________________________________________________________________________
void AliMUONTracker::CreateTrackReconstructor()
{
  /// Create track reconstructor, depending on tracking mode set in RecoParam
  
  TString opt(AliMUONReconstructor::GetRecoParam()->GetTrackingMode());
  opt.ToUpper();
  
  if (strstr(opt,"ORIGINAL"))
  {
    fTrackReco = new AliMUONTrackReconstructor();
  }
  else if (strstr(opt,"KALMAN"))
  {
    fTrackReco = new AliMUONTrackReconstructorK();
  }
  else
  {
    AliError(Form("tracking mode \"%s\" does not exist",opt.Data()));
    return;
  }
  
  AliInfo(Form("Will use %s for tracking",fTrackReco->ClassName()));
}

//_____________________________________________________________________________
void AliMUONTracker::UnloadClusters()
{
  /// Delete internal clusterStore
  delete fClusterStore;
  fClusterStore = 0x0;
}

