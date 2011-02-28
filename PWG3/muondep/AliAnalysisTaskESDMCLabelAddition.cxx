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
  if (!esd) {
    AliError("Cannot get input event");
    return;
  }      
  
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
    
    // skip ghosts
    if (!esdTrack->ContainTrackerData()) continue;
    
    // convert ESD track to MUON track (without recomputing track parameters at each clusters)
    AliMUONTrack muonTrack;
    AliMUONESDInterface::ESDToMUON(*esdTrack, muonTrack, kFALSE);
    
    // try to match the reconstructed track with a simulated one
    Int_t nMatchClusters = 0;
    AliMUONTrack* matchedTrackRef = rc.FindCompatibleTrack(muonTrack, *trackRefStore, nMatchClusters, kFALSE, fgkSigmaCut);
    
    // set the MC label
    if (matchedTrackRef) esdTrack->SetLabel(matchedTrackRef->GetUniqueID());
    else esdTrack->SetLabel(-1);
    
  }
  
}


//----------------------------------------------------------------------
void AliAnalysisTaskESDMCLabelAddition::Terminate(Option_t */*option*/)
{
  // Terminate analysis
  //
  AliDebug(2, "Terminate()");
}

