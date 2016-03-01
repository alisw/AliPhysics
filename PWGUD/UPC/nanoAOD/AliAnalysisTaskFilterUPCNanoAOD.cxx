/**************************************************************************
 * Copyright(c) 1998-2013, ALICE Experiment at CERN, All rights reserved. *
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

#include "AliAnalysisTaskFilterUPCNanoAOD.h"

#include "AliAnalysisManager.h"
#include "AliAnalysisTaskESDMuonFilter.h"
#include "AliAODEvent.h"
#include "AliAODHandler.h"
#include "AliAODUPCReplicator.h"
#include "AliLog.h"

///
/// \brief AliAnalysisTaskFilterUPCNanoAOD : a class to convert full AODs to muon ones
///
/// The actual work is done, like for the ESD->muon AOD case, by the
/// AliAODUPCReplicator class. This very class is just steering it.
///

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskFilterUPCNanoAOD)
/// \endcond

//_____________________________________________________________________________
AliAnalysisTaskFilterUPCNanoAOD::AliAnalysisTaskFilterUPCNanoAOD(Bool_t withSPDTracklets)
: AliAnalysisTaskSE("AliAnalysisTaskFilterUPCNanoAOD"),
fBranchReplicator(new AliAODUPCReplicator("UPCReplicator",
                                           "remove non muon tracks and non primary or pileup vertices",
                                           kTRUE,
                                           withSPDTracklets))
{
  /// ctor. For the parameters \see AliAODUPCReplicator::AliAODUPCReplicator
}

//_____________________________________________________________________________
AliAnalysisTaskFilterUPCNanoAOD::~AliAnalysisTaskFilterUPCNanoAOD()
{
  /// dtor. Delete our internal worker class
  delete fBranchReplicator;
}

//_____________________________________________________________________________
void AliAnalysisTaskFilterUPCNanoAOD::UserCreateOutputObjects()
{
  /// Create the output AOD branches, from the list of object
  /// we have in the replicator
  
  TList* list = fBranchReplicator->GetList();
  
  TIter next(list);
  TObject* o;
    
  while ( (o = next()) )
  {
    AddAODBranch(o->ClassName(),list->GetObjectRef(o));
  }  
}

//_____________________________________________________________________________
void AliAnalysisTaskFilterUPCNanoAOD::UserExec(Option_t*)
{
  /// Main method doing the actual filtering (delegating it to
  /// the muon replicator)
  
  AliAODEvent* aod = dynamic_cast<AliAODEvent*>(InputEvent());
  if (!aod)
  {
    AliError("Event is not a the AOD type");
    return;
  }
  
 //Trigger
  Bool_t isTriggered = kFALSE;
  TString trigger = aod->GetFiredTriggerClasses();

  if(trigger.Contains("CCUP4-B")) isTriggered = kTRUE; // *0VBA *0VBC *0UBA *0UBC 0SM2 0OMU
  if(trigger.Contains("CCUP2-B")) isTriggered = kTRUE; // *0VBA *0VBC *0UBA *0UBC 0STP 0OM2
  if(trigger.Contains("CCUP7-B")) isTriggered = kTRUE; // *0VBA *0VBC *0UBA *0UBC 0STP 0OMU
  if(trigger.Contains("CINT1-B")) isTriggered = kTRUE; // 0VBA || 0VBC || 0SMB
  if(trigger.Contains("CTEST58-B")) isTriggered = kTRUE; // *0VBA *0VBC *0UBA *0UBC 0SH1
  if(trigger.Contains("CTEST59-B")) isTriggered = kTRUE; // *0VBA *0VBC *0UBA *0UBC 0STP
  if(trigger.Contains("CTEST60-B")) isTriggered = kTRUE; // *0VBA *0VBC *0UBA *0UBC 0OM2
  if(trigger.Contains("CTEST61-B")) isTriggered = kTRUE; // *0VBA *0VBC *0UBA *0UBC 0OMU
  if(trigger.Contains("CCUP8-B")) isTriggered = kTRUE; //*0VBA *0VBC *0UBA *0UBC 0STP 0OMU
  if(trigger.Contains("CCUP9-B")) isTriggered = kTRUE; //*0VBA *0VBC *0UBA *0UBC 0STP
  if(trigger.Contains("CCUP10-B")) isTriggered = kTRUE; //*0VBA *0VBC *0UBA *0UBC 0SH1
  if(trigger.Contains("CCUP11-B")) isTriggered = kTRUE; //*0UBA *0UBC 0STP 0OMU
  if(trigger.Contains("CCUP12-B")) isTriggered = kTRUE; //*0UBA *0UBC 0STP
  
  //Vertex
  Bool_t hasGoodVertex = kFALSE;
  AliAODVertex *aodVertex = aod->GetPrimaryVertex();
  if(aodVertex->GetNContributors() > 1 && TMath::Abs(aodVertex->GetZ()) < 15) hasGoodVertex = kTRUE;
  
  UInt_t nGoodTracks = 0;
  for(Int_t itr=0; itr<aod ->GetNumberOfTracks(); itr++) {
  	AliAODTrack *trk = dynamic_cast<AliAODTrack*>(aod->GetTrack(itr));
  	if( !trk ) continue;
  	if( trk->TestFilterBit(1<<0) || trk->TestFilterBit(1<<1)) nGoodTracks++;
	}
  
  if(!isTriggered || !hasGoodVertex || nGoodTracks == 0) return;
  //if(!isTriggered) return;
  //AliInfo("Good UPC event");
  
  
  // not sure why but the fillAOD flag seems to be lost somewhere, so must
  // set it per event.
  AliAODHandler* aodH = static_cast<AliAODHandler*>(AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler());
  
  aodH->SetFillAOD(kTRUE);
  
  fBranchReplicator->ReplicateAndFilter(*aod);
}
