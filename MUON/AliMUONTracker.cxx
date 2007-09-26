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
  {
    fTrackHitPatternMaker = new AliMUONTrackHitPattern(*fTransformer,*fDigitMaker);
  }
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
Int_t 
AliMUONTracker::LoadClusters(TTree* clustersTree)
{
  /// Load clusterStore and triggerStore from clustersTree
  delete fClusterStore;
  delete fTriggerStore;

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
Int_t
AliMUONTracker::Clusters2Tracks(AliESDEvent* esd)
{
  /// Performs the tracking and store the resulting tracks in both
  /// the TreeT and the ESD
  
  Int_t rv(0);
 
  TTree *tracksTree = new TTree;
 
  if (!fClusterStore)
  {
    AliError("ClusterStore is NULL");
    rv=2;
  }
  if (!fTriggerStore)
  {
    AliError("TriggerStore is NULL");
    rv=3;
  }
  if (!rv)
  {
    rv = Clusters2Tracks(*tracksTree,esd);
  }
  return rv;
}

//_____________________________________________________________________________
Int_t AliMUONTracker::Clusters2Tracks(TTree& tracksTree, AliESDEvent* esd)
{
  /// Performs the tracking
  
  AliDebug(1,"");
  
  AliMUONVTrackStore* trackStore(0x0);
  AliMUONVTriggerTrackStore* triggerTrackStore(0x0);
  
  // Make tracker tracks
  if ( fClusterStore ) 
  {
    trackStore = new AliMUONTrackStoreV1;
    Bool_t alone = ( ( fTriggerStore && fTriggerCircuit ) ? kFALSE : kTRUE );
    trackStore->Connect(tracksTree,alone);
    fTrackReco->EventReconstruct(*fClusterStore,*trackStore);
  }
  
  if ( fTriggerStore && fTriggerCircuit )
  {
    // Make trigger tracks
    triggerTrackStore = new AliMUONTriggerTrackStoreV1;
    Bool_t alone = ( fClusterStore ? kFALSE : kTRUE );
    triggerTrackStore->Connect(tracksTree,alone);
    fTrackReco->EventReconstructTrigger(*fTriggerCircuit,*fTriggerStore,*triggerTrackStore);
  }

  if ( trackStore && triggerTrackStore && fTriggerStore && fTrackHitPatternMaker )
  {
    fTrackReco->ValidateTracksWithTrigger(*trackStore,*triggerTrackStore,*fTriggerStore,*fTrackHitPatternMaker);
  }
  
  // Fills output TreeT 
  tracksTree.Fill();

  if( trackStore && triggerTrackStore && fTriggerStore && fTrigChamberEff){
      fTrigChamberEff->EventChamberEff(*fTriggerStore,*triggerTrackStore,*trackStore);
  }

  FillESD(*trackStore,esd);
  
  // cleanup
  delete trackStore;
  delete triggerTrackStore;
  
  return 0;
}

//_____________________________________________________________________________
void 
AliMUONTracker::FillESD(AliMUONVTrackStore& trackStore, AliESDEvent* esd) const
{
  /// Fill the ESD from the trackStore
  
  AliDebug(1,"");
  
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
void AliMUONTracker::SetOption(Option_t* option)
{
  /// set reconstructor class
  
  if (strstr(option,"Original")) 
  {
    fTrackReco = new AliMUONTrackReconstructor;
  }
  else 
  {
    fTrackReco = new AliMUONTrackReconstructorK();
  }
}

//_____________________________________________________________________________
void 
AliMUONTracker::UnloadClusters()
{
  /// Delete internal clusterStore
  delete fClusterStore;
  fClusterStore = 0x0;
}

