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
#include <AliAODMCParticle.h>
#include <AliLog.h>
#include <AliAODRecoDecay.h>
#include <AliAODTrack.h>

#include "AliHFTrackContainer.h"

/// \cond CLASSIMP
ClassImp(AliHFAODMCParticleContainer)
/// \endcond

/// This is the default constructor, used for ROOT I/O purposes.
AliHFTrackContainer::AliHFTrackContainer() :
  AliTrackContainer(),
  fDMesonCandidate(0),
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
AliHFTrackContainer::AliHFTrackContainer(const char *name) :
  AliTrackContainer(name),
  fDMesonCandidate(0),
  fDaughterList(10)
{
  // Constructor.

  fBaseClassName = "AliAODTrack";
  SetClassName("AliAODTrack");
  fDaughterList.SetOwner(kFALSE);
}

/// First check whether the particle is a daughter of the
/// D meson candidate (in which case the particle is rejected);
/// then calls the base class AcceptTrack(AliVTrack*) method.
///
/// \param i Index of the particle to be checked.
///
/// \return kTRUE if the particle is accepted, kFALSE otherwise.
Bool_t AliHFTrackContainer::ApplyTrackCuts(const AliVTrack* vp, UInt_t &rejectionReason) const
{
  const AliAODTrack* part = static_cast<const AliAODTrack*>(vp);

  if (IsDMesonDaughter(part)) {
    rejectionReason = kHFCut;
    return kFALSE;  // daughter the D meson candidate
  }

  // Not a daughter of the D meson. Apply regular cuts.
  return AliTrackContainer::ApplyTrackCuts(vp, rejectionReason);
}

/// Check if particle it's a daughter of the D meson candidate
/// \param track Pointer to a valid AliAODTrack object
///
/// \result kTRUE if it is a daughter of the D meson candidate, kFALSE otherwise
Bool_t AliHFTrackContainer::IsDMesonDaughter(const AliAODTrack* track) const
{
  if (!fDMesonCandidate) return kFALSE;

  if (fDaughterList.FindObject(track)) return kTRUE;

  return kFALSE;
}

/// Set the D meson candidate pointer and generate the list of daughters
///
/// \param cand A pointer to an AliAODRecoDecay object
void AliHFTrackContainer::SetDMesonCandidate(AliAODRecoDecay* c)
{
  fDMesonCandidate = c;
  GenerateDaughterList();
}


/// Generate the list of the daughters of the D meson candidate
void AliHFTrackContainer::GenerateDaughterList()
{
  fDaughterList.Clear();
  AddDaughters(fDMesonCandidate);
}

// Add all the daughters of a D meson candidate in a TObjArray. Follows all the decay cascades.
///
/// \param cand A pointer to an AliAODRecoDecay object
void AliHFTrackContainer::AddDaughters(const AliAODRecoDecay* cand)
{
  if (!cand) return;

   Int_t n = cand->GetNDaughters();

   for (Int_t i = 0; i < n; i++) {
     AliVTrack* track = dynamic_cast<AliVTrack*>(cand->GetDaughter(i));
     if (!track) continue;

     AliAODRecoDecay* cand2 = dynamic_cast<AliAODRecoDecay*>(track);

     if (cand2) {
       AddDaughters(cand2);
     }
     else {
       if (!track->InheritsFrom("AliAODTrack")) {
         ::Warning("AliHFTrackContainer::AddDaughters", "One of the daughters is not of type 'AliAODTrack' nor 'AliAODRecoDecay'.");
         continue;
       }
       fDaughterList.AddLast(track);
     }
   }
}
