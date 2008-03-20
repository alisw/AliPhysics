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
/// Tracking modes (ORIGINAL, KALMAN) and associated options and parameters
/// can be changed by using:
/// AliMUONRecoParam *muonRecoParam = AliMUONRecoParam::GetLow(High)FluxParam();
/// muonRecoParam->Set...(); // see methods in AliMUONRecoParam.h for details
/// AliMUONReconstructor::SetRecoParam(muonRecoParam);
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
#include "AliMUONVCluster.h"
#include "AliMUONVClusterServer.h"
#include "AliMUONVDigitStore.h"
#include "AliMUONTrackReconstructor.h"
#include "AliMUONTrackReconstructorK.h"
#include "AliMUONTrackStoreV1.h"
#include "AliMUONTriggerTrackStoreV1.h"
#include "AliMUONClusterStoreV2.h"
#include "AliMUONVTriggerStore.h"
#include "AliMUONESDInterface.h"

#include "AliESDEvent.h"
#include "AliESDVertex.h"
#include "AliESDMuonTrack.h"
#include "AliLog.h"
#include "AliCodeTimer.h"

#include <Riostream.h>
#include <TTree.h>
#include <TRandom.h>

/// \cond CLASSIMP
ClassImp(AliMUONTracker)
/// \endcond


//_____________________________________________________________________________
AliMUONTracker::AliMUONTracker(AliMUONVClusterServer& clusterServer,
			       const AliMUONVDigitStore& digitStore,
                               const AliMUONDigitMaker* digitMaker,
                               const AliMUONGeometryTransformer* transformer,
                               const AliMUONTriggerCircuit* triggerCircuit)
: AliTracker(),
  fDigitMaker(digitMaker), // not owner
  fTransformer(transformer), // not owner
  fTriggerCircuit(triggerCircuit), // not owner
  fTrackHitPatternMaker(0x0),
  fTrackReco(0x0),
  fClusterStore(0x0),
  fTriggerStore(0x0),
  fClusterServer(clusterServer), // not owner
  fDigitStore(digitStore) // not owner
{
  /// constructor
  if (fTransformer && fDigitMaker)
    fTrackHitPatternMaker = new AliMUONTrackHitPattern(*fTransformer,*fDigitMaker);
  
  TIter next(fDigitStore.CreateIterator());
  fClusterServer.UseDigits(next);
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
AliMUONVClusterStore*
AliMUONTracker::ClusterStore() const
{
  /// Return (and create if necessary) the cluster container
  if (!fClusterStore) 
  {
    fClusterStore = new AliMUONClusterStoreV2;
  }
  return fClusterStore;
}

//_____________________________________________________________________________
Int_t AliMUONTracker::LoadClusters(TTree* clustersTree)
{
  /// Load triggerStore from clustersTree

  ClusterStore()->Clear();
  
  delete fTriggerStore;

  if ( ! clustersTree ) {
    AliFatal("No clustersTree");
    return 1;
  }

  fTriggerStore = AliMUONVTriggerStore::Create(*clustersTree);
  
  if (!fTriggerStore)
  {
    AliError("Could not get triggerStore");
    return 2;
  }
  
  ClusterStore()->Connect(*clustersTree,kFALSE);
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
  
  if ( ! ClusterStore() ) 
  {
    AliError("ClusterStore is NULL");
    return 2;
  }
  
  if (!fTriggerStore) {
    AliError("TriggerStore is NULL");
    return 3;
  }

  // Make tracker tracks
  AliMUONVTrackStore* trackStore = new AliMUONTrackStoreV1;
  fTrackReco->EventReconstruct(*(ClusterStore()),*trackStore);
  
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
  
  AliMUONTrack* track;
  AliESDMuonTrack esdTrack;
  Double_t vertex[3] = {0., 0., 0.};
  TIter next(trackStore.CreateIterator());
  
  // get ITS vertex
  const AliESDVertex* esdVert = esd->GetVertex(); 
  if (esdVert->GetNContributors()) {
    esdVert->GetXYZ(vertex);
    AliDebug(1,Form("found vertex (%e,%e,%e)",vertex[0],vertex[1],vertex[2]));
  }
  
  // fill ESD event including all info in ESD cluster if required and only for the given fraction of events
  if (AliMUONReconstructor::GetRecoParam()->SaveFullClusterInESD() && 
      gRandom->Uniform(100.) <= AliMUONReconstructor::GetRecoParam()->GetPercentOfFullClusterInESD()) {
    
    while ( ( track = static_cast<AliMUONTrack*>(next()) ) ) {
      AliMUONESDInterface::MUONToESD(*track, esdTrack, vertex, &fDigitStore);
      esd->AddMuonTrack(&esdTrack);
    }
    
  } else {
    
    while ( ( track = static_cast<AliMUONTrack*>(next()) ) ) {
      AliMUONESDInterface::MUONToESD(*track, esdTrack, vertex);
      esd->AddMuonTrack(&esdTrack);
    }
    
  }
  
}

//_____________________________________________________________________________
void AliMUONTracker::CreateTrackReconstructor()
{
  /// Create track reconstructor, depending on tracking mode set in RecoParam
  
  TString opt(AliMUONReconstructor::GetRecoParam()->GetTrackingMode());
  opt.ToUpper();
  
  if (strstr(opt,"ORIGINAL"))
  {
    fTrackReco = new AliMUONTrackReconstructor(fClusterServer);
  }
  else if (strstr(opt,"KALMAN"))
  {
    fTrackReco = new AliMUONTrackReconstructorK(fClusterServer);
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
  /// Clear internal clusterStore
  
  ClusterStore()->Clear();
}
