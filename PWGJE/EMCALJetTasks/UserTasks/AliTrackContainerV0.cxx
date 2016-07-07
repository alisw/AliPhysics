/*************************************************************************
* Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
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

#include <TClonesArray.h>
#include <AliLog.h>
#include <AliAODTrack.h>
#include <AliAODv0.h>

#include "AliTrackContainerV0.h"

/// \cond CLASSIMP
ClassImp(AliTrackContainerV0)
/// \endcond

/// This is the default constructor, used for ROOT I/O purposes.
AliTrackContainerV0::AliTrackContainerV0() :
  AliTrackContainer(),
  fFilterDaughterTracks(0),
  fV0s(0),
  fDaughterList(10)
{
  // Constructor.

  fBaseClassName = "AliAODTrack";
  SetClassName("AliAODTrack");
  fDaughterList.SetOwner(kFALSE);
}

/// This is the standard named constructor.
///
/// \param name Name of the particle collection
AliTrackContainerV0::AliTrackContainerV0(const char *name) :
  AliTrackContainer(name),
  fFilterDaughterTracks(0),
  fV0s(0),
  fDaughterList(10)
{
  // Constructor.

  fBaseClassName = "AliAODTrack";
  SetClassName("AliAODTrack");
  fDaughterList.SetOwner(kFALSE);
}

/// Get list of V0 candidates from AOD event.
/// (pointer itself does not change, only the array it points to)
///
/// \param event Pointer to the event.
void AliTrackContainerV0::SetArray(const AliVEvent *event)
{
  AliTrackContainer::SetArray(event);

  if(!fFilterDaughterTracks)
    return;
  
  //fEvent = dynamic_cast<AliAODEvent*>(const_cast<AliVEvent*>(event));  // alternative?
  fEvent = dynamic_cast<AliAODEvent*>((AliVEvent*)event->Clone()); 
}

/// Preparation for next event. 
/// Run in each event (via AliAnalysisTaskEmcal::RetrieveEventObjects)
void AliTrackContainerV0::NextEvent()
{
  AliTrackContainer::NextEvent();
  
  if(fFilterDaughterTracks)
  {
    // V0s daughter tracks will be removed from track sample

    if(!fEvent)
    {
      ::Warning("AliTrackContainerV0::NextEvent","fEvent is not valid!"); 
      return;
    }

    fV0s = dynamic_cast<TClonesArray*>(fEvent->GetV0s());
    if(!fV0s)
    {
      ::Warning("AliTrackContainerV0::NextEvent","fV0s is not valid!"); 
      return;
    }

    fDaughterList.Clear();
    Int_t iNumV0s = fV0s->GetEntriesFast();
    
    for(Int_t iV0 = 0; iV0 < iNumV0s; iV0++)
      ExtractDaughters(dynamic_cast<AliAODv0*>(fV0s->At(iV0)));
  }
}

/// Extract the V0 candidate daughter tracks and add them to daughter track list. 
/// 
/// \param cand Pointer to V0 candidate to be set.
void AliTrackContainerV0::ExtractDaughters(AliAODv0* cand)
{
  for(Int_t i = 0; i < 2; i++)
  {   
    AliAODTrack* track = dynamic_cast<AliAODTrack*>(cand->GetDaughter(i));
    
    if(!track)
    {
      ::Warning("AliTrackContainerV0::ExtractDaughters","Track is not valid! Skipping candidate.");
      return;
    }

    fDaughterList.AddLast(track);
  }
}

/// Check whether the particle is a daughter of
/// V0 candidate (in which case the particle is rejected);
/// then calls the base class AcceptTrack(AliVTrack*) method.
///
/// \param vp Pointer to track to be checked.
/// \param rejectionReason Rejection reason.
///
/// \return kTRUE if the particle is accepted, kFALSE otherwise.
Bool_t AliTrackContainerV0::ApplyTrackCuts(const AliVTrack* vp, UInt_t &rejectionReason) const
{
  const AliAODTrack* track = dynamic_cast<const AliAODTrack*>(vp);

  if(IsV0Daughter(track))
  {
    // track is one of the V0 daughter - will be rejected
    return kFALSE;
  } 
  else
  {
    // track is NOT a V0 daughter - regular cuts will be applied
    return AliTrackContainer::ApplyTrackCuts(vp, rejectionReason);
  }
}

/// Check whether the particle is among the daughter tracks of the V0 candidate 
/// 
/// \param track Pointer to input track to be checked.
///
/// \return kTRUE if the particle is a daughter of V0 candidate, kFALSE otherwise
Bool_t AliTrackContainerV0::IsV0Daughter(const AliAODTrack* track) const
{
	if(fDaughterList.FindObject(track)) 
		return kTRUE;

	return kFALSE;
}

