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
/// Tracking modes (ORIGINAL, KALMAN) and associated options and parameters can be changed
/// through the AliMUONRecoParam object set in the reconstruction macro or read from the CDB
/// (see methods in AliMUONRecoParam.h file for details)
///
/// \author Christian Finck and Laurent Aphecetche, SUBATECH Nantes
//-----------------------------------------------------------------------------

#include "AliMUONTracker.h"

#include "AliCodeTimer.h"
#include "AliESDEvent.h"
#include "AliESDMuonTrack.h"
#include "AliESDVertex.h"
#include "AliLog.h"
#include "AliMUONClusterStoreV2.h"
#include "AliMUONESDInterface.h"
#include "AliMUONLegacyClusterServer.h"
#include "AliMUONRecoParam.h"
#include "AliMUONReconstructor.h"
#include "AliMUONTrack.h"
#include "AliMUONTrackExtrap.h"
#include "AliMUONTrackHitPattern.h"
#include "AliMUONTrackParam.h"
#include "AliMUONTrackReconstructor.h"
#include "AliMUONTrackReconstructorK.h"
#include "AliMUONTrackStoreV1.h"
#include "AliMUONTriggerTrackStoreV1.h"
#include "AliMUONTriggerTrack.h"
#include "AliMUONLocalTrigger.h"
#include "AliMUONVClusterServer.h"
#include "AliMUONVDigitStore.h"
#include "AliMUONVTriggerStore.h"
#include "AliMUONTriggerUtilities.h"
#include <Riostream.h>
#include <TRandom.h>
#include <TTree.h>

/// \cond CLASSIMP
ClassImp(AliMUONTracker)
/// \endcond


//_____________________________________________________________________________
AliMUONTracker::AliMUONTracker(const AliMUONRecoParam* recoParam,
                               AliMUONVClusterServer* clusterServer,
                               AliMUONVDigitStore& digitStore,
                               const AliMUONGeometryTransformer* transformer,
                               const AliMUONTriggerCircuit* triggerCircuit,
                               const AliMUONTriggerUtilities* triggerUtilities)
: AliTracker(),
fkTransformer(transformer), // not owner
fkTriggerCircuit(triggerCircuit), // not owner
fTrackHitPatternMaker(0x0),
fTrackReco(0x0),
fClusterStore(0x0),
fTriggerStore(0x0),
fClusterServer(clusterServer), 
fIsOwnerOfClusterServer(kFALSE),
fkDigitStore(digitStore), // not owner
fInputClusterStore(0x0),
fTriggerTrackStore(0x0),
fkRecoParam(recoParam)
{
  /// constructor
  if (fkTransformer)
    fTrackHitPatternMaker = new AliMUONTrackHitPattern(recoParam,*fkTransformer,fkDigitStore,triggerUtilities);
  
  if (!fClusterServer)
  {
    AliDebug(1,"No cluster server given. Will use AliMUONLegacyClusterServer");
    fIsOwnerOfClusterServer = kTRUE;
  }
  else
  {
    TIter next(fkDigitStore.CreateIterator());
    fClusterServer->UseDigits(next,&digitStore);
    
    SetupClusterServer(*fClusterServer);
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
  if ( fIsOwnerOfClusterServer ) delete fClusterServer;
  delete fInputClusterStore;
  delete fTriggerTrackStore;
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
AliMUONVTriggerTrackStore*
AliMUONTracker::TriggerTrackStore() const
{
  /// Return (and create if necessary) the trigger track container
  if (!fTriggerTrackStore) 
  {
    fTriggerTrackStore = new AliMUONTriggerTrackStoreV1;
  }
  return fTriggerTrackStore;
}

//_____________________________________________________________________________
Int_t AliMUONTracker::LoadClusters(TTree* clustersTree)
{
  /// Load triggerStore from clustersTree

  delete fTriggerStore;
  delete fInputClusterStore;
  fInputClusterStore=0x0;

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
  
  if ( fIsOwnerOfClusterServer )
  {
    fInputClusterStore = AliMUONVClusterStore::Create(*clustersTree);
    if ( fInputClusterStore ) 
    {
      AliDebug(1,Form("Created %s from cluster tree",fInputClusterStore->ClassName()));
      fInputClusterStore->Clear();
      fInputClusterStore->Connect(*clustersTree,kFALSE);
    }
    delete fClusterServer;
    fClusterServer = new AliMUONLegacyClusterServer(*fkTransformer,fInputClusterStore,
																										GetRecoParam()->BypassSt4(),
																										GetRecoParam()->BypassSt5());
    SetupClusterServer(*fClusterServer);
  }
  
  fTriggerStore->Connect(*clustersTree,kFALSE);
  
  clustersTree->GetEvent(0);

  return 0;
}

//_____________________________________________________________________________
Int_t AliMUONTracker::Clusters2Tracks(AliESDEvent* esd)
{
  /// Performs the tracking and store the resulting tracks in the ESD
  AliDebug(1,"");
  AliCodeTimerAuto("",0)
  
  if (!fTrackReco) 
  {
    fTrackReco = CreateTrackReconstructor(GetRecoParam(),fClusterServer);
  }
  
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

  // Make trigger tracks
  if ( fkTriggerCircuit ) 
  {
    TriggerTrackStore()->Clear();
    fTrackReco->EventReconstructTrigger(*fkTriggerCircuit,*fTriggerStore,*(TriggerTrackStore()));
  }
  
  if ( TriggerTrackStore()->GetSize() > GetRecoParam()->GetMaxTriggerTracks() ) 
  {
    // cut to reject shower events
    
    AliCodeTimerAuto("MUON Shower events",1);

    AliWarning(Form("Probably got a shower event (%d trigger tracks). Will not reconstruct tracks.",
                    TriggerTrackStore()->GetSize()));
    
    return 0;
  }
       
  // Make tracker tracks
  AliMUONVTrackStore* trackStore = new AliMUONTrackStoreV1;
  fTrackReco->EventReconstruct(*(ClusterStore()),*trackStore);
  
  // Match tracker/trigger tracks
  if ( fTrackHitPatternMaker ) 
  {
    fTrackReco->ValidateTracksWithTrigger(*trackStore,*(TriggerTrackStore()),*fTriggerStore,*fTrackHitPatternMaker);
  }
  
  // Fill ESD
  FillESD(*trackStore,esd);
  
  // cleanup
  delete trackStore;
  
  return 0;
}

//_____________________________________________________________________________
void AliMUONTracker::FillESD(const AliMUONVTrackStore& trackStore, AliESDEvent* esd) const
{
  /// Fill the ESD from the trackStore
  AliDebug(1,"");
  AliCodeTimerAuto("",0)
  
  // get ITS vertex
  Double_t vertex[3] = {0., 0., 0.};
  const AliESDVertex* esdVert = esd->GetVertex(); 
  if (esdVert->GetNContributors() > 0 || !strcmp(esdVert->GetTitle(),"vertexer: smearMC")) {
    esdVert->GetXYZ(vertex);
    AliDebug(1,Form("found vertex (%e,%e,%e)",vertex[0],vertex[1],vertex[2]));
  }
  
  // fill ESD event including all info in ESD cluster if required and only for the given fraction of events
  AliMUONTrack* track;
  AliMUONLocalTrigger* locTrg;
  AliESDMuonTrack esdTrack;
  TIter next(trackStore.CreateIterator());
  if (GetRecoParam()->SaveFullClusterInESD() && 
      gRandom->Uniform(100.) <= GetRecoParam()->GetPercentOfFullClusterInESD()) {
    
    while ( ( track = static_cast<AliMUONTrack*>(next()) ) ) {
      
      if (track->GetMatchTrigger() > 0) {
	locTrg = static_cast<AliMUONLocalTrigger*>(fTriggerStore->FindLocal(track->LoCircuit()));
	AliMUONESDInterface::MUONToESD(*track, esdTrack, vertex, &fkDigitStore, locTrg);
      } else AliMUONESDInterface::MUONToESD(*track, esdTrack, vertex, &fkDigitStore);
      
      esd->AddMuonTrack(&esdTrack);
    }
    
  } else {
    
    while ( ( track = static_cast<AliMUONTrack*>(next()) ) ) {
      
      if (track->GetMatchTrigger() > 0) {
	locTrg = static_cast<AliMUONLocalTrigger*>(fTriggerStore->FindLocal(track->LoCircuit()));
	AliMUONESDInterface::MUONToESD(*track, esdTrack, vertex, 0x0, locTrg);
      } else AliMUONESDInterface::MUONToESD(*track, esdTrack, vertex);
      
      esd->AddMuonTrack(&esdTrack);
    }
    
  }
  
  // fill the local trigger decisions not matched with tracks (associate them to "ghost" tracks)
  UInt_t ghostId = 0xFFFFFFFF - 1;
  Bool_t matched = kFALSE;
  AliMUONTriggerTrack *triggerTrack;
  TIter itTriggerTrack(fTriggerTrackStore->CreateIterator());
  while ( ( triggerTrack = static_cast<AliMUONTriggerTrack*>(itTriggerTrack()) ) ) {
    
    locTrg = static_cast<AliMUONLocalTrigger*>(fTriggerStore->FindLocal(triggerTrack->GetLoTrgNum()));
    
    // check if this local trigger has already been matched
    TIter itTrack(trackStore.CreateIterator());
    while ( ( track = static_cast<AliMUONTrack*>(itTrack()) ) ) {
      matched = (track->LoCircuit() == locTrg->LoCircuit());
      if (matched) break;
    }
    if (matched) continue;

    AliMUONESDInterface::MUONToESD(*locTrg, esdTrack, ghostId, triggerTrack);
    
    esd->AddMuonTrack(&esdTrack);
    ghostId -= 1;
  }
  
}

//_____________________________________________________________________________
AliMUONVTrackReconstructor* AliMUONTracker::CreateTrackReconstructor(const AliMUONRecoParam* recoParam, AliMUONVClusterServer* clusterServer)
{
  /// Create track reconstructor, depending on tracking mode set in RecoParam
  
  AliMUONVTrackReconstructor* trackReco(0x0);
  
  TString opt(recoParam->GetTrackingMode());
  opt.ToUpper();
  
  if (strstr(opt,"ORIGINAL"))
  {
    trackReco = new AliMUONTrackReconstructor(recoParam,clusterServer);
  }
  else if (strstr(opt,"KALMAN"))
  {
    trackReco = new AliMUONTrackReconstructorK(recoParam,clusterServer);
  }
  else
  {
    AliErrorClass(Form("tracking mode \"%s\" does not exist",opt.Data()));
    return 0x0;
  }
  
  AliDebugClass(1,Form("Will use %s for tracking",trackReco->ClassName()));
  
  return trackReco;
}

//_____________________________________________________________________________
void AliMUONTracker::UnloadClusters()
{
  /// Clear internal clusterStore
  
  delete fInputClusterStore;
  fInputClusterStore = 0x0;
}


//_____________________________________________________________________________
void
AliMUONTracker::SetupClusterServer(AliMUONVClusterServer& clusterServer)
{
  /// Setup the cluster server
  
  if ( GetRecoParam()->BypassSt4() ||
			 GetRecoParam()->BypassSt5() )
  {
    Bool_t ok = clusterServer.UseTriggerTrackStore(TriggerTrackStore());
  
		TString msg1;
		TString msg2;
		
		if ( GetRecoParam()->BypassSt45() )
		{
			msg1 = "STATIONS 4 AND 5";
			msg2 = "THOSE TWO STATIONS";
		}
		else if ( GetRecoParam()->BypassSt4() )
		{
			msg1 = "STATION 4";
			msg2 = "THAT STATION";
		}
		else if ( GetRecoParam()->BypassSt5() )
		{
			msg1 = "STATION 5";
			msg2 = "THAT STATION";
		}
		
    if ( ok ) 
    {
      AliWarning(Form("WILL USE TRIGGER TRACKS TO GENERATE CLUSTERS IN %s, "
											"THUS BYPASSING REAL CLUSTERS IN %s!!!",msg1.Data(),msg2.Data()));    
    }
    else
    {
      AliWarning("BYPASSING OF ST4 AND/OR 5 REQUESTED, BUT CLUSTERSERVER DOES NOT SEEM TO SUPPORT IT !!!");    
    }
  }
}


