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
#include <TClonesArray.h>

#include "AliVEvent.h"
#include "AliLog.h"
#include "AliVCuts.h"
#include "AliESDtrack.h"

#include "AliTLorentzVector.h"
#include "AliEmcalTrackSelectionAOD.h"
#include "AliEmcalTrackSelectionESD.h"
#include "AliTrackContainer.h"

/// \cond CLASSIMP
ClassImp(AliTrackContainer);
/// \endcond

TString AliTrackContainer::fgDefTrackCutsPeriod = "";

/**
 * Default constructor.
 */
AliTrackContainer::AliTrackContainer():
  AliParticleContainer(),
  fTrackFilterType(AliEmcalTrackSelection::kHybridTracks),
  fListOfCuts(0),
  fSelectionModeAny(kFALSE),
  fAODFilterBits(0),
  fTrackCutsPeriod(),
  fEmcalTrackSelection(0),
  fFilteredTracks(0),
  fTrackTypes(5000)
{
  fBaseClassName = "AliVTrack";
  SetClassName("AliVTrack");
  fMassHypothesis = 0.139;
}

/**
 * Standard constructor.
 * @param[in] name Name of the container (= name of the array operated on)
 * @param[in] period Name of the period, needed for track cut selection
 */
AliTrackContainer::AliTrackContainer(const char *name, const char *period):
  AliParticleContainer(name),
  fTrackFilterType(AliEmcalTrackSelection::kHybridTracks),
  fListOfCuts(0),
  fSelectionModeAny(kFALSE),
  fAODFilterBits(0),
  fTrackCutsPeriod(period),
  fEmcalTrackSelection(0),
  fFilteredTracks(0),
  fTrackTypes(5000)
{
  fBaseClassName = "AliVTrack";
  SetClassName("AliVTrack");

  if (fTrackCutsPeriod.IsNull() && !AliTrackContainer::fgDefTrackCutsPeriod.IsNull()) {
    AliInfo(Form("Default track cuts period is %s", AliTrackContainer::fgDefTrackCutsPeriod.Data()));
    fTrackCutsPeriod = AliTrackContainer::fgDefTrackCutsPeriod;
  }
  fMassHypothesis = 0.139;
}

/**
 * Get array from event. Also
 * creating the virtual track selection
 * for the period provided in the constructor.
 * @param[in] event Event from which the data is read
 */
void AliTrackContainer::SetArray(const AliVEvent *event)
{
  AliParticleContainer::SetArray(event);

  if (fTrackFilterType == AliEmcalTrackSelection::kNoTrackFilter) {
    if (fEmcalTrackSelection) delete fEmcalTrackSelection;
    fEmcalTrackSelection = 0;
  }
  else {
    if (fTrackFilterType == AliEmcalTrackSelection::kCustomTrackFilter) {

      AliInfo("Using custom track cuts");

      if (fLoadedClass) {
        if (fLoadedClass->InheritsFrom("AliAODTrack")) {
          AliInfo(Form("Objects are of type %s: AOD track selection will be done.", fLoadedClass->GetName()));
          fEmcalTrackSelection = new AliEmcalTrackSelectionAOD(0, fAODFilterBits);
        }
        else if (fLoadedClass->InheritsFrom("AliESDtrack")) {
          AliInfo(Form("Objects are of type %s: ESD track selection will be done.", fLoadedClass->GetName()));
          fEmcalTrackSelection = new AliEmcalTrackSelectionESD(0);
        }
        else {
          AliWarning(Form("Objects are of type %s: no track filtering will be done!!", fLoadedClass->GetName()));
        }
      }

      if (fEmcalTrackSelection) {
        if (fSelectionModeAny) {
          fEmcalTrackSelection->SetSelectionModeAny();
        }
        else {
          fEmcalTrackSelection->SetSelectionModeAll();
        }

        fEmcalTrackSelection->AddTrackCuts(fListOfCuts);
      }
    }
    else {
      if (!fTrackCutsPeriod.IsNull()) {
        AliInfo(Form("Using track cuts %d for period %s", fTrackFilterType, fTrackCutsPeriod.Data()));
      }
      else {
        AliInfo(Form("Using track cuts %d (no data period was provided!)", fTrackFilterType));
      }

      if (fLoadedClass->InheritsFrom("AliAODTrack")) {
        AliInfo(Form("Objects are of type %s: AOD track selection will be done.", fLoadedClass->GetName()));
        fEmcalTrackSelection = new AliEmcalTrackSelectionAOD(fTrackFilterType, fTrackCutsPeriod);
      }
      else if (fLoadedClass->InheritsFrom("AliESDtrack")) {
        AliInfo(Form("Objects are of type %s: ESD track selection will be done.", fLoadedClass->GetName()));
        fEmcalTrackSelection = new AliEmcalTrackSelectionESD(fTrackFilterType, fTrackCutsPeriod);
      }
      else {
        AliWarning(Form("Objects are of type %s: no track filtering will be done!!", fLoadedClass->GetName()));
      }
    }
  }
}

/**
 * Preparation for the next event: Run the track
 * selection of all bit and store the pointers to
 * selected tracks in a separate array.
 */
void AliTrackContainer::NextEvent()
{
  fTrackTypes.Reset(kUndefined);
  if (fEmcalTrackSelection) {
    fFilteredTracks = fEmcalTrackSelection->GetAcceptedTracks(fClArray);

    const TClonesArray* trackBitmaps = fEmcalTrackSelection->GetAcceptedTrackBitmaps();
    TIter nextBitmap(trackBitmaps);
    TBits* bits = 0;
    Int_t i = 0;
    while ((bits = static_cast<TBits*>(nextBitmap()))) {
      if (i >= fTrackTypes.GetSize()) fTrackTypes.Set((i+1)*2);
      AliVTrack* vTrack = static_cast<AliVTrack*>(fFilteredTracks->At(i));
      if (!vTrack) {
        fTrackTypes[i] = kRejected;
      }
      else if (fTrackFilterType == AliEmcalTrackSelection::kHybridTracks) {
        if (bits->FirstSetBit() == 0) {
          fTrackTypes[i] = kHybridGlobal;
        }
        else if (bits->FirstSetBit() == 1) {
          if ((vTrack->GetStatus()&AliVTrack::kITSrefit) != 0) {
            fTrackTypes[i] = kHybridConstrained;
          }
          else {
            fTrackTypes[i] = kHybridConstrainedNoITSrefit;
          }
        }
      }
     i++;
    }
  }
  else {
    fFilteredTracks = fClArray;
  }
}

/**
 * Get track at index in the container
 * @param[in] i Index of the particle in the container
 * @return pointer to particle if particle is accepted, NULL otherwise
 */
AliVTrack* AliTrackContainer::GetTrack(Int_t i) const
{
  //Get i^th jet in array

  if (i < 0 || i >= fFilteredTracks->GetEntriesFast()) return 0;
  AliVTrack *vp = static_cast<AliVTrack*>(fFilteredTracks->At(i));
  return vp;
}

/**
 * Get track at index in the container if accepted by the track selection provided
 * @param[in] i Index of the particle in the container
 * @return pointer to particle if particle is accepted, NULL otherwise
 */
AliVTrack* AliTrackContainer::GetAcceptTrack(Int_t i) const
{
  UInt_t rejectionReason;
  if (i == -1) i = fCurrentID;
  if (AcceptTrack(i, rejectionReason)) {
      return GetTrack(i);
  }
  else {
    AliDebug(2,"Track not accepted.");
    return 0;
  }
}

/**
 * Get next accepted particle in the container selected using the track cuts provided.
 * @deprecated Old style iterator - for compatibility reasons, use AliParticleContainer::accept_iterator instead
 * @return Next accepted particle (NULL if the end of the array is reached)
 */
AliVTrack* AliTrackContainer::GetNextAcceptTrack()
{
  const Int_t n = GetNEntries();
  AliVTrack *p = 0;
  do {
    fCurrentID++;
    if (fCurrentID >= n) break;
    p = GetAcceptTrack(fCurrentID);
  } while (!p);

  return p;
}

/**
 * Get next particle in the container
 * @deprecated Old style iterator - for compatibility reasons, use AliParticleContainer::all_iterator instead
 * @return Next track in the container (NULL if end of the container is reached)
 */
AliVTrack* AliTrackContainer::GetNextTrack()
{
  const Int_t n = GetNEntries();
  AliVTrack *p = 0;
  do {
    fCurrentID++;
    if (fCurrentID >= n) break;
    p = GetTrack(fCurrentID);
  } while (!p);

  return p;
}

/**
 * Retrieve the track type using the IndexOf method of the TClonesArray
 * to retrieve the index of the provided track.
 * @param track Track for which the type is requested.
 * @return The track type from fTrackTypes
 */
Char_t AliTrackContainer::GetTrackType(const AliVTrack* track) const
{
  Int_t id = fFilteredTracks->IndexOf(track);
  if (id >= 0) {
    return fTrackTypes[id];
  }
  else {
    return kUndefined;
  }
}

/**
 * Retrieve momentum information of a track and fill a TLorentzVector
 * with it. In case the optional parameter mass is provided, it is used as mass
 * hypothesis, otherwise the mass hypothesis from the particle itself is used.
 * @param[out] mom Momentum vector to be filled
 * @param[in] track Track from which the momentum information is obtained.
 * @param[in] mass (Optional) Mass hypothesis
 * @return
 */
Bool_t AliTrackContainer::GetMomentumFromTrack(TLorentzVector &mom, const AliVTrack* track, Double_t mass) const
{
  if (track) {
    if (mass < 0) mass = track->M();

    Bool_t useConstrainedParams = kFALSE;
    if (fLoadedClass->InheritsFrom("AliESDtrack") && fTrackFilterType == AliEmcalTrackSelection::kHybridTracks) {
      Char_t trackType = GetTrackType(track);
      if (trackType == kHybridConstrained || trackType == kHybridConstrainedNoITSrefit) {
        useConstrainedParams = kTRUE;
      }
    }

    if (useConstrainedParams) {
      const AliESDtrack *esdtrack = static_cast<const AliESDtrack*>(track);
      mom.SetPtEtaPhiM(esdtrack->GetConstrainedParam()->Pt(), esdtrack->GetConstrainedParam()->Eta(), esdtrack->GetConstrainedParam()->Phi(), mass);
    }
    else {
      mom.SetPtEtaPhiM(track->Pt(), track->Eta(), track->Phi(), mass);
    }
    return kTRUE;
  }
  else {
    mom.SetPtEtaPhiM(0, 0, 0, 0);
    return kFALSE;
  }
}

/**
 * Fills a TLorentzVector with the momentum information of the track provided
 * under a global mass hypothesis.
 * @param[out] mom Momentum vector of the particle provided
 * @param[in] track Track from which to obtain the momentum information
 * @return Always true
 */
Bool_t AliTrackContainer::GetMomentumFromTrack(TLorentzVector &mom, const AliVTrack* part) const
{
  return GetMomentumFromTrack(mom,part,fMassHypothesis);
}

/**
 * Fills a TLorentzVector with the momentum information of the
 * \f$ i^{th} \f$ particle in the container, using a global
 * mass hypothesis. In case the provided index is out of
 * range, false is returned as return value.
 * @param[out] mom Momentum vector of the \f$ i^{th} \f$ particle in the array
 * @param[in] i Index of the particle to check
 * @return True if the request was successful, false otherwise
 */
Bool_t AliTrackContainer::GetMomentum(TLorentzVector &mom, Int_t i) const
{
  Double_t mass = fMassHypothesis;

  if (i == -1) i = fCurrentID;
  AliVTrack *vp = GetTrack(i);
  if (vp) {
    if (mass < 0) mass = vp->M();

    if (fLoadedClass->InheritsFrom("AliESDtrack") && fTrackFilterType == AliEmcalTrackSelection::kHybridTracks &&
        (fTrackTypes[i] == kHybridConstrained || fTrackTypes[i] == kHybridConstrainedNoITSrefit)) {
      AliESDtrack *track = static_cast<AliESDtrack*>(vp);
      mom.SetPtEtaPhiM(track->GetConstrainedParam()->Pt(), track->GetConstrainedParam()->Eta(), track->GetConstrainedParam()->Phi(), mass);
    }
    else {
      mom.SetPtEtaPhiM(vp->Pt(), vp->Eta(), vp->Phi(), mass);
    }
    return kTRUE;
  }
  else {
    mom.SetPtEtaPhiM(0, 0, 0, 0);
    return kFALSE;
  }
}

/**
 * Fills a TLorentzVector with the momentum information of the
 * next particle in the container, using a global mass hypothesis.
 * In case the iterator reached the end of the array, false
 * is returned as return value.
 * @deprecated Old style iterator - use all_iterator instead
 * @param[out] mom Momentum vector of the next particle
 * @return True if the request was successful, false otherwise
 */
Bool_t AliTrackContainer::GetNextMomentum(TLorentzVector &mom)
{
  Double_t mass = fMassHypothesis;

  AliVTrack *vp = GetNextTrack();
  if (vp) {
    if (mass < 0) mass = vp->M();

    if (fLoadedClass->InheritsFrom("AliESDtrack") && fTrackFilterType == AliEmcalTrackSelection::kHybridTracks &&
        (fTrackTypes[fCurrentID] == kHybridConstrained || fTrackTypes[fCurrentID] == kHybridConstrainedNoITSrefit)) {
      AliESDtrack *track = static_cast<AliESDtrack*>(vp);
      mom.SetPtEtaPhiM(track->GetConstrainedParam()->Pt(), track->GetConstrainedParam()->Eta(), track->GetConstrainedParam()->Phi(), mass);
    }
    else {
      mom.SetPtEtaPhiM(vp->Pt(), vp->Eta(), vp->Phi(), mass);
    }
    return kTRUE;
  }
  else {
    mom.SetPtEtaPhiM(0, 0, 0, 0);
    return kFALSE;
  }
}

/**
 * Fills a TLorentzVector with the monentum infomation of the
 * \f$ i^{th} \f$ accepted particle in the container, using a
 * global mass hypothesis. In case the provided index is out of
 * range, or the particle under the index is not accepted, false
 * is returned as return value.
 * @param[out] mom Momentum vector of the accepted particle
 * @param[in] i Index to check
 * @return True if the request was successfull, false otherwise
 */
Bool_t AliTrackContainer::GetAcceptMomentum(TLorentzVector &mom, Int_t i) const
{

  Double_t mass = fMassHypothesis;

  if (i == -1) i = fCurrentID;
  AliVTrack *vp = GetAcceptTrack(i);
  if (vp) {
    if (mass < 0) mass = vp->M();

    if (fLoadedClass->InheritsFrom("AliESDtrack") && fTrackFilterType == AliEmcalTrackSelection::kHybridTracks &&
        (fTrackTypes[i] == kHybridConstrained || fTrackTypes[i] == kHybridConstrainedNoITSrefit)) {
      AliESDtrack *track = static_cast<AliESDtrack*>(vp);
      mom.SetPtEtaPhiM(track->GetConstrainedParam()->Pt(), track->GetConstrainedParam()->Eta(), track->GetConstrainedParam()->Phi(), mass);
    }
    else {
      mom.SetPtEtaPhiM(vp->Pt(), vp->Eta(), vp->Phi(), mass);
    }

    return kTRUE;
  }
  else {
    mom.SetPtEtaPhiM(0, 0, 0, 0);
    return kFALSE;
  }
}

/**
 * Fills a TLorentzVector with the monentum infomation of the
 * next accepted particle in the container, using a global
 * mass hypothesis. In case the iteration reached the end of
 * the array, false is returned as return value.
 * @deprecated Old style iterator - use accept_iterator instead
 * @param[out] mom Momentum vector of the next particle in the array
 * @return True if the request was successfull, false (no more entries) otherwise
 */
Bool_t AliTrackContainer::GetNextAcceptMomentum(TLorentzVector &mom)
{
  Double_t mass = fMassHypothesis;

  AliVTrack *vp = GetNextAcceptTrack();
  if (vp) {
    if (mass < 0) mass = vp->M();

    if (fLoadedClass->InheritsFrom("AliESDtrack") && fTrackFilterType == AliEmcalTrackSelection::kHybridTracks &&
        (fTrackTypes[fCurrentID] == kHybridConstrained || fTrackTypes[fCurrentID] == kHybridConstrainedNoITSrefit)) {
      AliESDtrack *track = static_cast<AliESDtrack*>(vp);
      mom.SetPtEtaPhiM(track->GetConstrainedParam()->Pt(), track->GetConstrainedParam()->Eta(), track->GetConstrainedParam()->Phi(), mass);
    }
    else {
      mom.SetPtEtaPhiM(vp->Pt(), vp->Eta(), vp->Phi(), mass);
    }

    return kTRUE;
  }
  else {
    mom.SetPtEtaPhiM(0, 0, 0, 0);
    return kFALSE;
  }
}

/**
 * Perform full track selection for the particle vp, consisting
 * of kinematical track selection and track quality
 * cut provided by the cuts assigned to this container
 * @param[in] vp Track to be checked
 * @param[in] rejectionReason Bitmap encoding the reason why the
 * track was rejected. Note: The variable is not set to NULL
 * inside this function before changing its value.
 * @return True if the track is accepted, false otherwise
 */
Bool_t AliTrackContainer::AcceptTrack(const AliVTrack *vp, UInt_t &rejectionReason) const
{
  Bool_t r = ApplyTrackCuts(vp, rejectionReason);
  if (!r) return kFALSE;

  AliTLorentzVector mom;
  GetMomentumFromTrack(mom, vp);

  return ApplyKinematicCuts(mom, rejectionReason);
}

/**
 * Perform full track selection for the particle in
 * the container stored at position i, consisting
 * of kinematical track selection and track quality
 * cut provided by the cuts assigned to this container
 * @param[in] i Index of the track to check
 * @param[in] rejectionReason Bitmap encoding the reason why the
 * track was rejected. Note: The variable is not set to NULL
 * inside this function before changing its value.
 * @return True if the track is accepted, false otherwise
 */
Bool_t AliTrackContainer::AcceptTrack(Int_t i, UInt_t &rejectionReason) const
{
  Bool_t r = ApplyTrackCuts(GetTrack(i), rejectionReason);
  if (!r) return kFALSE;

  AliTLorentzVector mom;
  GetMomentum(mom, i);

  return ApplyKinematicCuts(mom, rejectionReason);
}

/**
 * Perform track quality selection of the track vp using
 * the track cuts assigned to this container.
 * @param[in] vp Track to be checked
 * @param[out] rejectionReason Bitmap encoding the reason why the
 * track was rejected. Note: The variable is not set to NULL
 * inside this function before changing its value.
 * @return True if the particle was accepted, false otherwise
 */
Bool_t AliTrackContainer::ApplyTrackCuts(const AliVTrack* vp, UInt_t &rejectionReason) const
{
  return ApplyParticleCuts(vp, rejectionReason);
}

/**
 * Add new track cuts to the container.
 * @param[in] cuts Cuts to be  added
 */
void AliTrackContainer::AddTrackCuts(AliVCuts *cuts)
{
  if (!fListOfCuts) {
    fListOfCuts = new TObjArray;
    fListOfCuts->SetOwner(true);
  }
  fListOfCuts->Add(cuts);
}

/**
 * Get number of track cut objects assigned to this container.
 * @return Number of track cut objects
 */
Int_t AliTrackContainer::GetNumberOfCutObjects() const
{
  if (!fListOfCuts) return 0;
  return fListOfCuts->GetEntries();
}

/**
 * Get the cut object at index (icut) assigned to this container.
 * @param[in] icut Index of the cut in the container
 * @return Cut object at the index if existing, NULL otherwise
 */
AliVCuts* AliTrackContainer::GetTrackCuts(Int_t icut)
{
  if (!fListOfCuts) return NULL;
  if (icut < fListOfCuts->GetEntries()) {
    return static_cast<AliVCuts *>(fListOfCuts->At(icut));
  }
  return NULL;
}

/**
 * Create an iterable container interface over all objects in the
 * EMCAL container.
 * @return iterable container over all objects in the EMCAL container
 */
const AliTrackIterableContainer AliTrackContainer::all() const {
  return AliTrackIterableContainer(this, false);
}

/**
 * Create an iterable container interface over accepted objects in the
 * EMCAL container.
 * @return iterable container over accepted objects in the EMCAL container
 */
const AliTrackIterableContainer AliTrackContainer::accepted() const {
  return AliTrackIterableContainer(this, true);
}

/**
 * Create an iterable container interface over all objects in the
 * EMCAL container.
 * @return iterable container over all objects in the EMCAL container
 */
const AliTrackIterableMomentumContainer AliTrackContainer::all_momentum() const {
  return AliTrackIterableMomentumContainer(this, false);
}

/**
 * Create an iterable container interface over accepted objects in the
 * EMCAL container.
 * @return iterable container over accepted objects in the EMCAL container
 */
const AliTrackIterableMomentumContainer AliTrackContainer::accepted_momentum() const {
  return AliTrackIterableMomentumContainer(this, true);
}

/**
 * Build title of the container consisting of the container name
 * and a string encoding the minimum \f$ p_{t} \f$ cut applied
 * in the kinematic track selection.
 * @return Title of the container
 */
const char* AliTrackContainer::GetTitle() const
{
  static TString trackString;

  if (GetMinPt() == 0) {
    trackString = TString::Format("%s_pT0000", GetArrayName().Data());
  }
  else if (GetMinPt() < 1.0) {
    trackString = TString::Format("%s_pT0%3.0f", GetArrayName().Data(), GetMinPt()*1000.0);
  }
  else {
    trackString = TString::Format("%s_pT%4.0f", GetArrayName().Data(), GetMinPt()*1000.0);
  }

  return trackString.Data();
}
