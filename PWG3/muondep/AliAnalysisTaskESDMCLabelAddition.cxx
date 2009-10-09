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

#include <TChain.h>
#include <TFile.h>
#include <TParticle.h>

#include "AliAnalysisTaskESDMCLabelAddition.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliESDInputHandler.h"
#include "AliAODHandler.h"
#include "AliAnalysisFilter.h"
#include "AliESDtrack.h"
#include "AliESDMuonTrack.h"
#include "AliESDMuonCluster.h"
#include "AliESDVertex.h"
#include "AliMultiplicity.h"
#include "AliLog.h"
#include "AliStack.h"
#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
#include "AliAODMCParticle.h"

#include "AliMUONRecoCheck.h"
#include "AliMUONESDInterface.h"
#include "AliMUONTrack.h"
#include "AliMUONTrackParam.h"
#include "AliMUONVTrackStore.h"
#include "AliMUONVCluster.h"
#include "AliMUONVClusterStore.h"

ClassImp(AliAnalysisTaskESDMCLabelAddition)

// sigma cut applied to match a reconstructed cluster with a trackref
const Double_t AliAnalysisTaskESDMCLabelAddition::fgkSigmaCut = 10.;

//----------------------------------------------------------------------
AliAnalysisTaskESDMCLabelAddition::AliAnalysisTaskESDMCLabelAddition():
  AliAnalysisTaskSE()
{
  // Default constructor
}


//----------------------------------------------------------------------
AliAnalysisTaskESDMCLabelAddition::AliAnalysisTaskESDMCLabelAddition(const char* name):
  AliAnalysisTaskSE(name)
{
  // Constructor
}


//----------------------------------------------------------------------
void AliAnalysisTaskESDMCLabelAddition::UserCreateOutputObjects()
{
}


//----------------------------------------------------------------------
void AliAnalysisTaskESDMCLabelAddition::Init()
{
  AliDebug(2, "Init()");
}


//----------------------------------------------------------------------
void AliAnalysisTaskESDMCLabelAddition::UserExec(Option_t */*option*/)
{
  // Execute analysis for current event					    
  Long64_t ientry = Entry();
  AliDebug(1, Form("MCLabel Addition: Analysing event # %5d\n",(Int_t) ientry)); 
  AddMCLabel();
}


//----------------------------------------------------------------------
void AliAnalysisTaskESDMCLabelAddition::AddMCLabel() 
{
  // Load ESD event
  AliESDEvent* esd = dynamic_cast<AliESDEvent*>(InputEvent());
  
  // Load MC event 
  AliMCEventHandler *mcH = 0;
  if(MCEvent()) mcH = (AliMCEventHandler*) ((AliAnalysisManager::GetAnalysisManager())->GetMCtruthEventHandler()); 
  
  // Get reference tracks
  AliMUONRecoCheck rc(esd,mcH);
  AliMUONVTrackStore* trackRefStore = rc.TrackRefs(-1);
  
  // Loop over reconstructed tracks
  AliESDMuonTrack *esdTrack = 0x0;
  Int_t nMuTracks = esd->GetNumberOfMuonTracks();
  for (Int_t nMuTrack = 0; nMuTrack < nMuTracks; ++nMuTrack) {
    
    esdTrack = esd->GetMuonTrack(nMuTrack);
    esdTrack->SetLabel(-1);
    
    // skip ghosts
    if (!esdTrack->ContainTrackerData()) continue;
    
    // try to match the reconstructed track with a simulated one
    AliMUONTrack* matchedTrackRef = MatchWithTrackRef(*esdTrack, *trackRefStore);
    
    if (matchedTrackRef) {
      
      // set the MC label
      esdTrack->SetLabel(matchedTrackRef->GetUniqueID());
       
      // remove already matched trackRefs
      trackRefStore->Remove(*matchedTrackRef);
      
    }
    
  }
  
}


//----------------------------------------------------------------------
void AliAnalysisTaskESDMCLabelAddition::Terminate(Option_t */*option*/)
{
  // Terminate analysis
  //
  AliDebug(2, "Terminate()");
}


//----------------------------------------------------------------------
AliMUONTrack* AliAnalysisTaskESDMCLabelAddition::ESDToMUON(AliESDMuonTrack &esdTrack)
{
  /// Convert an ESD track into a MUON track with dummy parameters (just to fill the clusters).
  /// It is the responsability of the user to delete the track afterward
  
  AliMUONTrack *track = new AliMUONTrack();
  
  // ckeck whether the ESD track contains clusters
  if(!esdTrack.ClustersStored()) return track;
  
  // track parameters at first cluster
  AliMUONTrackParam param;
  AliMUONESDInterface::GetParamAtFirstCluster(esdTrack, param);
  AliMUONESDInterface::GetParamCov(esdTrack, param);
  
  // create empty cluster
  AliMUONVClusterStore* cStore = AliMUONESDInterface::NewClusterStore();
  AliMUONVCluster* cluster = cStore->CreateCluster(0,0,0);
  
  // loop over ESD clusters
  AliESDMuonCluster *esdCluster = (AliESDMuonCluster*) esdTrack.GetClusters().First();
  while (esdCluster) {
    
    // copy cluster information
    AliMUONESDInterface::ESDToMUON(*esdCluster, *cluster);
    
    // only set the Z parameter to avoid error in the AddTrackParamAtCluster(...) method
    param.SetZ(cluster->GetZ());
    
    // add common track parameters at current cluster
    track->AddTrackParamAtCluster(param, *cluster, kTRUE);
    
    esdCluster = (AliESDMuonCluster*) esdTrack.GetClusters().After(esdCluster);
  }
  
  // clean memory
  delete cluster;
  delete cStore;
  
  return track;
  
}


//----------------------------------------------------------------------
AliMUONTrack* AliAnalysisTaskESDMCLabelAddition::MatchWithTrackRef(AliESDMuonTrack &esdTrack,
								   AliMUONVTrackStore &trackRefStore)
{
  /// Return the trackRef matched with the reconstructed track and the fraction of matched clusters
  
  AliMUONTrack *matchedTrackRef = 0x0;
  
  // convert ESD track to MUON track
  AliMUONTrack *track = ESDToMUON(esdTrack);
  
  // look for the corresponding simulated track if any
  TIter next(trackRefStore.CreateIterator());
  AliMUONTrack* trackRef;
  while ( ( trackRef = static_cast<AliMUONTrack*>(next()) ) ) {
    
    // check compatibility
    if (TrackMatched(*track, *trackRef)) {
      matchedTrackRef = trackRef;
      break;
    }
    
  }
  
  // clean memory
  delete track;
  
  return matchedTrackRef;
  
}


//----------------------------------------------------------------------
Bool_t AliAnalysisTaskESDMCLabelAddition::TrackMatched(AliMUONTrack &track, AliMUONTrack &trackRef)
{
  /// Try to match 2 tracks
  
  Bool_t compTrack[10];
  Int_t nMatchClusters = track.CompatibleTrack(&trackRef, fgkSigmaCut, compTrack);
  Double_t fractionOfMatchCluster = ((Double_t)nMatchClusters) / ((Double_t)track.GetNClusters());
  
  if ((compTrack[0] || compTrack[1] || compTrack[2] || compTrack[3]) && // at least 1 cluster matched in st 1 & 2
      (compTrack[6] || compTrack[7] || compTrack[8] || compTrack[9]) && // at least 1 cluster matched in st 4 & 5
      fractionOfMatchCluster > 0.5) return kTRUE;                       // more than 50% of clusters matched
  else return kFALSE;
  
}

