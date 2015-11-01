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

#include "AliAnalysisTaskAOD2MuonAOD.h"

#include "AliAnalysisManager.h"
#include "AliAnalysisNonMuonTrackCuts.h"
#include "AliAnalysisNonPrimaryVertices.h"
#include "AliAnalysisTaskESDMuonFilter.h"
#include "AliAODEvent.h"
#include "AliAODHandler.h"
#include "AliAODMuonReplicator.h"
#include "AliLog.h"
#include "AliMuonEventCuts.h"

///
/// \brief AliAnalysisTaskAOD2MuonAOD : a class to convert full AODs to muon ones
///
/// The actual work is done, like for the ESD->muon AOD case, by the
/// AliAODMuonReplicator class. This very class is just steering it.
///

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskAOD2MuonAOD)
/// \endcond

//_____________________________________________________________________________
AliAnalysisTaskAOD2MuonAOD::AliAnalysisTaskAOD2MuonAOD(Int_t mcMode, Bool_t withSPDTracklets, AliMuonEventCuts* muonEventCuts)
: AliAnalysisTaskSE("AliAnalysisTaskAOD2MuonAOD"),
fBranchReplicator(new AliAODMuonReplicator("MuonReplicator",
                                           "remove non muon tracks and non primary or pileup vertices",
                                           new AliAnalysisNonMuonTrackCuts,
                                           new AliAnalysisNonPrimaryVertices,
                                           mcMode,
                                           kTRUE,
                                           withSPDTracklets)),
fMuonEventCuts(muonEventCuts)
{
  /// ctor. For the parameters \see AliAODMuonReplicator::AliAODMuonReplicator
}

//_____________________________________________________________________________
AliAnalysisTaskAOD2MuonAOD::~AliAnalysisTaskAOD2MuonAOD()
{
  /// dtor. Delete our internal worker class
  delete fBranchReplicator;
}

//_____________________________________________________________________________
void AliAnalysisTaskAOD2MuonAOD::UserCreateOutputObjects()
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
Bool_t AliAnalysisTaskAOD2MuonAOD::HasMuonInformation(const AliAODEvent& event) const
{
  /// Determine whether or not this event has any information related to muon
  /// spectrometer. This means either some muon tracks or just some L0 inputs (or both)
  
  AliAODHeader* header = static_cast<AliAODHeader*>(event.GetHeader());
  if ( header->GetNumberOfMuons()>0 )
  {
    return kTRUE;
  }
  
  const TObjArray* selectTrigClasses = fMuonEventCuts->GetSelectedTrigClassesInEvent(InputEvent());
  
  if ( selectTrigClasses->GetEntries() > 0 )
  {
    return kTRUE;
  }
  return kFALSE;
}

//_____________________________________________________________________________
void AliAnalysisTaskAOD2MuonAOD::UserExec(Option_t*)
{
  /// Main method doing the actual filtering (delegating it to
  /// the muon replicator)
  
  AliAODEvent* event = dynamic_cast<AliAODEvent*>(InputEvent());
  if (!event)
  {
    AliError("Event is not a the AOD type");
    return;
  }
  
  if ( fMuonEventCuts && !HasMuonInformation(*event) )
  {
    return;
  }

  // not sure why but the fillAOD flag seems to be lost somewhere, so must
  // set it per event.
  AliAODHandler* aodH = static_cast<AliAODHandler*>(AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler());
  
  aodH->SetFillAOD(kTRUE);
  
  fBranchReplicator->ReplicateAndFilter(*event);
}
