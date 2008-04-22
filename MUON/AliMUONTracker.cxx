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
#include "AliMUONVCluster.h"
#include "AliMUONVClusterServer.h"
#include "AliMUONVDigitStore.h"
#include "AliMUONVTriggerStore.h"
#include <Riostream.h>
#include <TRandom.h>
#include <TTree.h>

/// \cond CLASSIMP
ClassImp(AliMUONTracker)
/// \endcond


//_____________________________________________________________________________
AliMUONTracker::AliMUONTracker(AliMUONVClusterServer* clusterServer,
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
  fClusterServer(clusterServer), 
  fIsOwnerOfClusterServer(kFALSE),
  fDigitStore(digitStore), // not owner
  fInputClusterStore(0x0),
  fTriggerTrackStore(0x0)
{
  /// constructor
  if (fTransformer && fDigitMaker)
    fTrackHitPatternMaker = new AliMUONTrackHitPattern(*fTransformer,*fDigitMaker);
  
  if (!fClusterServer)
  {
    AliInfo("No cluster server given. Will use AliMUONLegacyClusterServer");
    fIsOwnerOfClusterServer = kTRUE;
  }
  else
  {
    TIter next(fDigitStore.CreateIterator());
    fClusterServer->UseDigits(next);
    
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
      AliInfo(Form("Created %s from cluster tree",fInputClusterStore->ClassName()));
      fInputClusterStore->Clear();
      fInputClusterStore->Connect(*clustersTree,kFALSE);
    }
    delete fClusterServer;
    fClusterServer = new AliMUONLegacyClusterServer(*fTransformer,fInputClusterStore);
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
  AliCodeTimerAuto("")
  
  if (!fTrackReco) 
  {
    fTrackReco = CreateTrackReconstructor(AliMUONReconstructor::GetRecoParam()->GetTrackingMode(),fClusterServer);
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
  if ( fTriggerCircuit ) 
  {
    TriggerTrackStore()->Clear();
    fTrackReco->EventReconstructTrigger(*fTriggerCircuit,*fTriggerStore,*(TriggerTrackStore()));
  }
  
  if ( AliMUONReconstructor::GetRecoParam()->BypassSt45() && TriggerTrackStore()->GetSize() > 5 ) 
  {
    // Hard cut to reject shower events
    
    AliCodeTimerAuto("MUON Shower events");

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
      // set the trigger x/y strips pattern
      if (esdTrack.GetMatchTrigger()) {
	AliMUONLocalTrigger* locTrg = static_cast<AliMUONLocalTrigger*>(fTriggerStore->FindLocal(esdTrack.LoCircuit()));
	esdTrack.SetTriggerX1Pattern(locTrg->GetX1Pattern());
	esdTrack.SetTriggerY1Pattern(locTrg->GetY1Pattern());
	esdTrack.SetTriggerX2Pattern(locTrg->GetX2Pattern());
	esdTrack.SetTriggerY2Pattern(locTrg->GetY2Pattern());
	esdTrack.SetTriggerX3Pattern(locTrg->GetX3Pattern());
	esdTrack.SetTriggerY3Pattern(locTrg->GetY3Pattern());
	esdTrack.SetTriggerX4Pattern(locTrg->GetX4Pattern());
	esdTrack.SetTriggerY4Pattern(locTrg->GetY4Pattern());
      }
      esd->AddMuonTrack(&esdTrack);
    }
    
  } else {
    
    while ( ( track = static_cast<AliMUONTrack*>(next()) ) ) {
      AliMUONESDInterface::MUONToESD(*track, esdTrack, vertex);
      // set the trigger x/y strips pattern
      if (esdTrack.GetMatchTrigger()) {
	AliMUONLocalTrigger* locTrg = static_cast<AliMUONLocalTrigger*>(fTriggerStore->FindLocal(esdTrack.LoCircuit()));
	esdTrack.SetTriggerX1Pattern(locTrg->GetX1Pattern());
	esdTrack.SetTriggerY1Pattern(locTrg->GetY1Pattern());
	esdTrack.SetTriggerX2Pattern(locTrg->GetX2Pattern());
	esdTrack.SetTriggerY2Pattern(locTrg->GetY2Pattern());
	esdTrack.SetTriggerX3Pattern(locTrg->GetX3Pattern());
	esdTrack.SetTriggerY3Pattern(locTrg->GetY3Pattern());
	esdTrack.SetTriggerX4Pattern(locTrg->GetX4Pattern());
	esdTrack.SetTriggerY4Pattern(locTrg->GetY4Pattern());
      }
      esd->AddMuonTrack(&esdTrack);
    }
    
  }

  // fill the local trigger decisions not matched with tracks 
  // associate them to "ghost" tracks

  Int_t loTrgNum(-1);
  Bool_t matched = kFALSE;

  AliMUONTriggerTrack *triggerTrack;
  AliMUONTrack muonTrack;
  AliESDMuonTrack esdGhostTrack;
  TIter itTriggerTrack(fTriggerTrackStore->CreateIterator());
  while ( ( triggerTrack = static_cast<AliMUONTriggerTrack*>(itTriggerTrack() )
) )
    {
      loTrgNum = triggerTrack->GetLoTrgNum();
      AliMUONLocalTrigger* locTrg = static_cast<AliMUONLocalTrigger*>(fTriggerStore->FindLocal(loTrgNum));

      /* verify if this local trigger has been already matched */
      TIter itTrack(trackStore.CreateIterator());
      while ( ( track = static_cast<AliMUONTrack*>(itTrack()) ) )
        {
          if (matched = (track->LoCircuit() == locTrg->LoCircuit())) break;
        }
      if (matched) continue;

      muonTrack.SetLocalTrigger(locTrg->LoCircuit(),
                                locTrg->LoStripX(),
                                locTrg->LoStripY(),
                                locTrg->LoDev(),
                                locTrg->LoLpt(),
                                locTrg->LoHpt());

      /* make the AliESDMuonTrack from the "track" object */

      esdGhostTrack.SetLocalTrigger(muonTrack.GetLocalTrigger());
      // set the transverse momentum of this track to "zero"
      esdGhostTrack.SetInverseBendingMomentum(1.E+10);
      esdGhostTrack.SetInverseBendingMomentumAtDCA(1.E+10);
      esdGhostTrack.SetInverseBendingMomentumUncorrected(1.E+10);
      // set the trigger x/y strips pattern
      esdGhostTrack.SetTriggerX1Pattern(locTrg->GetX1Pattern());
      esdGhostTrack.SetTriggerY1Pattern(locTrg->GetY1Pattern());
      esdGhostTrack.SetTriggerX2Pattern(locTrg->GetX2Pattern());
      esdGhostTrack.SetTriggerY2Pattern(locTrg->GetY2Pattern());
      esdGhostTrack.SetTriggerX3Pattern(locTrg->GetX3Pattern());
      esdGhostTrack.SetTriggerY3Pattern(locTrg->GetY3Pattern());
      esdGhostTrack.SetTriggerX4Pattern(locTrg->GetX4Pattern());
      esdGhostTrack.SetTriggerY4Pattern(locTrg->GetY4Pattern());
      esd->AddMuonTrack(&esdGhostTrack);

    }
  
}

//_____________________________________________________________________________
AliMUONVTrackReconstructor* AliMUONTracker::CreateTrackReconstructor(const char* trackingMode, AliMUONVClusterServer* clusterServer)
{
  /// Create track reconstructor, depending on tracking mode set in RecoParam
  
  AliMUONVTrackReconstructor* trackReco(0x0);
  
  TString opt(trackingMode);
  opt.ToUpper();
  
  if (strstr(opt,"ORIGINAL"))
  {
    trackReco = new AliMUONTrackReconstructor(*clusterServer);
  }
  else if (strstr(opt,"KALMAN"))
  {
    trackReco = new AliMUONTrackReconstructorK(*clusterServer);
  }
  else
  {
    AliErrorClass(Form("tracking mode \"%s\" does not exist",opt.Data()));
    return 0x0;
  }
  
  AliInfoClass(Form("Will use %s for tracking",trackReco->ClassName()));
  
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
  
  if ( AliMUONReconstructor::GetRecoParam()->BypassSt45() )
  {
    Bool_t ok = clusterServer.UseTriggerTrackStore(TriggerTrackStore());
  
    if ( ok ) 
    
    {
      AliWarning("WILL USE TRIGGER TRACKS TO GENERATE CLUSTERS IN STATIONS 4 AND 5, THUS BYPASSING REAL CLUSTERS IN THOSE TWO STATIONS !!!");    
    }
    else
    {
      AliWarning("BYPASSING OF ST45 REQUESTED, BUT CLUSTERSERVER DOES NOT SEEM TO SUPPORT IT !!!");    
    }
  }
}


