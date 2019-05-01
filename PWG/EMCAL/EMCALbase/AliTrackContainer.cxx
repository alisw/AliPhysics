/************************************************************************************
 * Copyright (C) 2017, Copyright Holders of the ALICE Collaboration                 *
 * All rights reserved.                                                             *
 *                                                                                  *
 * Redistribution and use in source and binary forms, with or without               *
 * modification, are permitted provided that the following conditions are met:      *
 *     * Redistributions of source code must retain the above copyright             *
 *       notice, this list of conditions and the following disclaimer.              *
 *     * Redistributions in binary form must reproduce the above copyright          *
 *       notice, this list of conditions and the following disclaimer in the        *
 *       documentation and/or other materials provided with the distribution.       *
 *     * Neither the name of the <organization> nor the                             *
 *       names of its contributors may be used to endorse or promote products       *
 *       derived from this software without specific prior written permission.      *
 *                                                                                  *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND  *
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED    *
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE           *
 * DISCLAIMED. IN NO EVENT SHALL ALICE COLLABORATION BE LIABLE FOR ANY              *
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES       *
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;     *
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND      *
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT       *
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS    *
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                     *
 ************************************************************************************/
#include <bitset>
#include <iostream>
#include <TClonesArray.h>

#include "AliAODEvent.h"
#include "AliESDEvent.h"
#include "AliVEvent.h"
#include "AliLog.h"
#include "AliVCuts.h"
#include "AliESDtrack.h"

#include "AliTLorentzVector.h"
#include "AliEmcalTrackSelectionAOD.h"
#include "AliEmcalTrackSelectionESD.h"
#include "AliEmcalTrackSelResultPtr.h"
#include "AliEmcalTrackSelResultCombined.h"
#include "AliEmcalTrackSelResultHybrid.h"
#include "AliTrackContainer.h"

/// \cond CLASSIMP
ClassImp(AliTrackContainer);
/// \endcond

TString AliTrackContainer::fgDefTrackCutsPeriod = "";

// string to enum map for use with the %YAML config
const std::map <std::string, AliEmcalTrackSelection::ETrackFilterType_t> AliTrackContainer::fgkTrackFilterTypeMap = {
  {"kNoTrackFilter", AliEmcalTrackSelection::kNoTrackFilter },
  {"kCustomTrackFilter", AliEmcalTrackSelection::kCustomTrackFilter },
  {"kHybridTracks",  AliEmcalTrackSelection::kHybridTracks },
  {"kTPCOnlyTracks", AliEmcalTrackSelection::kTPCOnlyTracks },
  {"kITSPureTracks", AliEmcalTrackSelection::kITSPureTracks },
  {"kHybridTracks2010wNoRefit", AliEmcalTrackSelection::kHybridTracks2010wNoRefit },
  {"kHybridTracks2010woNoRefit", AliEmcalTrackSelection::kHybridTracks2010woNoRefit },
  {"kHybridTracks2011wNoRefit", AliEmcalTrackSelection::kHybridTracks2011wNoRefit },
  {"kHybridTracks2011woNoRefit", AliEmcalTrackSelection::kHybridTracks2011woNoRefit }
};

/**
 * Default constructor.
 */
AliTrackContainer::AliTrackContainer():
  AliParticleContainer(),
  fTrackFilterType(AliEmcalTrackSelection::kHybridTracks),
  fListOfCuts(0),
  fSelectionModeAny(kFALSE),
  fITSHybridTrackDistinction(kFALSE),
  fAODFilterBits(0),
  fTrackCutsPeriod(),
  fEmcalTrackSelection(0),
  fFilteredTracks(),
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
  fITSHybridTrackDistinction(kFALSE),
  fAODFilterBits(0),
  fTrackCutsPeriod(period),
  fEmcalTrackSelection(0),
  fFilteredTracks(),
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

        if(fListOfCuts) fEmcalTrackSelection->AddTrackCuts(fListOfCuts);
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
void AliTrackContainer::NextEvent(const AliVEvent * event)
{
  AliParticleContainer::NextEvent(event);

  fTrackTypes.Reset(kUndefined);
  if (fEmcalTrackSelection) {
    auto acceptedTracks = fEmcalTrackSelection->GetAcceptedTracks(fClArray);

    TObjArray *trackarray(fFilteredTracks.GetData());
    if(!trackarray){
      trackarray = new TObjArray;
      trackarray->SetOwner(false);
      fFilteredTracks.SetObject(trackarray);
      fFilteredTracks.SetOwner(true);
    } else {
      trackarray->Clear();
    }

    int naccepted(0), nrejected(0), nhybridTracks1(0), nhybridTracks2a(0), nhybridTracks2b(0), nhybridTracks3(0);
    Int_t i = 0;
    for(auto accresult : *acceptedTracks) {
      if (i >= fTrackTypes.GetSize()) fTrackTypes.Set((i+1)*2);
      PWG::EMCAL::AliEmcalTrackSelResultPtr *selectionResult = static_cast<PWG::EMCAL::AliEmcalTrackSelResultPtr *>(accresult);
      AliVTrack *vTrack = selectionResult->GetTrack();
      trackarray->AddLast(vTrack);
      if (!(*selectionResult) || !vTrack) {
        nrejected++;
        fTrackTypes[i] = kRejected;
      }
      else{ 
        // track is accepted;
        naccepted++;
        if (IsHybridTrackSelection()) {
          switch(GetHybridDefinition(*selectionResult)) {
            case PWG::EMCAL::AliEmcalTrackSelResultHybrid::kHybridGlobal:
              fTrackTypes[i] = kHybridGlobal;
              nhybridTracks1++;
              break;
            case PWG::EMCAL::AliEmcalTrackSelResultHybrid::kHybridConstrainedTrue:
              fTrackTypes[i] = kHybridConstrainedTrue;
              nhybridTracks2a++;
              break;
            case PWG::EMCAL::AliEmcalTrackSelResultHybrid::kHybridConstrainedFake:
              fTrackTypes[i] = kHybridConstrainedFake;
              nhybridTracks2b++;
              break;
            case PWG::EMCAL::AliEmcalTrackSelResultHybrid::kHybridConstrainedNoITSrefit:
              fTrackTypes[i] = kHybridConstrainedNoITSrefit;
              nhybridTracks3++;
              break;
            case PWG::EMCAL::AliEmcalTrackSelResultHybrid::kUndefined:
              fTrackTypes[i] = kRejected; // should in principle never happen
              break;
          };
        }
      }
     i++;
    }
    AliDebugStream(1) << "Accepted: " << naccepted << ", Rejected: " << nrejected << ", hybrid: (" << nhybridTracks1 << " | [" << nhybridTracks2a << " | " << nhybridTracks2b  << "] | " << nhybridTracks3 << ")" << std::endl;
  }
  else {
    fFilteredTracks.SetOwner(false);
    fFilteredTracks.SetObject(fClArray);
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

  if (i < 0 || i >= fFilteredTracks.GetData()->GetEntriesFast()) return 0;
  AliVTrack *vp = static_cast<AliVTrack*>(fFilteredTracks.GetData()->At(i));
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
  Int_t id = fFilteredTracks.GetData()->IndexOf(track);
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
    if (fLoadedClass->InheritsFrom("AliESDtrack") && IsHybridTrackSelection()) {
      Char_t trackType = GetTrackType(track);
      if (trackType == kHybridConstrainedTrue || trackType == kHybridConstrainedNoITSrefit) {
        AliDebugStream(2) << "Found a constrained track" << std::endl;
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

    if (fLoadedClass->InheritsFrom("AliESDtrack") && IsHybridTrackSelection() &&
        (fTrackTypes[i] == kHybridConstrainedTrue || fTrackTypes[i] == kHybridConstrainedNoITSrefit)) {
      AliDebugStream(2) << "Found a constrained track" << std::endl;
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

    if (fLoadedClass->InheritsFrom("AliESDtrack") && IsHybridTrackSelection() &&
        (fTrackTypes[fCurrentID] == kHybridConstrainedTrue || fTrackTypes[fCurrentID] == kHybridConstrainedNoITSrefit)) {
      AliDebugStream(2) << "Found a constrained track" << std::endl; 
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

    if (fLoadedClass->InheritsFrom("AliESDtrack") && IsHybridTrackSelection() &&
        (fTrackTypes[i] == kHybridConstrainedTrue || fTrackTypes[i] == kHybridConstrainedNoITSrefit)) {
      AliDebugStream(2) << "Found a constrained track" << std::endl;
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

    if (fLoadedClass->InheritsFrom("AliESDtrack") && IsHybridTrackSelection() &&
        (fTrackTypes[fCurrentID] == kHybridConstrainedTrue || fTrackTypes[fCurrentID] == kHybridConstrainedNoITSrefit)) {
      AliDebugStream(2) << "Found a constrained track" << std::endl;
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
  if(!GetMomentumFromTrack(mom, vp)) return false;

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
  if(fTrackTypes[i] == kRejected) return false; // track was rejected by the track selection
  Bool_t r = ApplyTrackCuts(GetTrack(i), rejectionReason);
  if (!r) return kFALSE;

  AliTLorentzVector mom;
  if(!GetMomentum(mom, i)) return false;

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

bool AliTrackContainer::IsHybridTrackSelection() const {
  return (fTrackFilterType == AliEmcalTrackSelection::kHybridTracks) ||
         (fTrackFilterType == AliEmcalTrackSelection::kHybridTracks2010wNoRefit) ||
         (fTrackFilterType == AliEmcalTrackSelection::kHybridTracks2010woNoRefit) ||
         (fTrackFilterType == AliEmcalTrackSelection::kHybridTracks2011wNoRefit) ||
         (fTrackFilterType == AliEmcalTrackSelection::kHybridTracks2011woNoRefit) ||
         (fTrackFilterType == AliEmcalTrackSelection::kHybridTracks2018TRD);
}

PWG::EMCAL::AliEmcalTrackSelResultHybrid::HybridType_t AliTrackContainer::GetHybridDefinition(const PWG::EMCAL::AliEmcalTrackSelResultPtr &selectionResult) const {
  PWG::EMCAL::AliEmcalTrackSelResultHybrid::HybridType_t hybridDefinition = PWG::EMCAL::AliEmcalTrackSelResultHybrid::kUndefined;
  if(auto hybriddata = dynamic_cast<const PWG::EMCAL::AliEmcalTrackSelResultHybrid *>(selectionResult.GetUserInfo())) {
    hybridDefinition = hybriddata->GetHybridTrackType();
  } else {
    if(auto combineddata = dynamic_cast<const PWG::EMCAL::AliEmcalTrackSelResultCombined *>(selectionResult.GetUserInfo())){
      for(int icut = 0; icut < combineddata->GetNumberOfSelectionResults(); icut++){
        try{
          auto cutresult = GetHybridDefinition((*combineddata)[icut]);
          if(cutresult != PWG::EMCAL::AliEmcalTrackSelResultHybrid::kUndefined) hybridDefinition = cutresult;
        } catch(PWG::EMCAL::AliEmcalTrackSelResultCombined::IndexException &e) {
          AliErrorStream() << "Index error: " << e.what() << std::endl;
        }
      }
    }
  }
  return hybridDefinition;
}

Bool_t AliTrackContainer::CheckArrayConsistency() const {
  bool teststatus = true;
  auto selected = fFilteredTracks.GetData();
  if(selected->GetEntries() != fClArray->GetEntries()) {
    std::cout << "Mismatch array size: selected " << selected->GetEntries() << ", input " << fClArray->GetEntries() << std::endl; 
    teststatus = false;
  } else {
    std::cout << "Array size consistent: " << selected->GetEntries() << std::endl;
  }

  // Now check whehter tracks are stored in the same indices:
  if(teststatus){
    int nmatch(0), nfail(0);
    for(int i = 0; i < fClArray->GetEntries(); i++) {
      AliVTrack *trackAll = dynamic_cast<AliVTrack *>(fClArray->At(i));
      AliVTrack *trackSel = dynamic_cast<AliVTrack *>(selected->At(i));
      if(trackSel == trackAll) {
        nmatch++;    
      } else {
        std::cout << "Mismatch in array position: " << i << std::endl;
        nfail++;
      }
    }
    if(nfail) {
      std::cout << "found " << nfail << "mismatches" << std::endl;
      teststatus = false;
    }
  }
  return teststatus;
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
  trackString = TString::Format("%s_pT%04d", GetArrayName().Data(), static_cast<int>(GetMinPt()*1000.0));
  return trackString.Data();
}

TString AliTrackContainer::GetDefaultArrayName(const AliVEvent *const ev) const {
  if(ev->IsA() == AliAODEvent::Class()) return "tracks";
  else if(ev->IsA() == AliESDEvent::Class()) return "Tracks";
  else return "";
}

AliTrackContainer::TrackOwnerHandler::TrackOwnerHandler():
  TObject(),
  fManagedObject(nullptr),
  fOwnership(false)
{
  
}

AliTrackContainer::TrackOwnerHandler::TrackOwnerHandler(TObjArray *managedobject, Bool_t ownership):
  TObject(),
  fManagedObject(managedobject),
  fOwnership(ownership)
{

}

AliTrackContainer::TrackOwnerHandler::TrackOwnerHandler(const AliTrackContainer::TrackOwnerHandler &other) :
  TObject(other),
  fManagedObject(other.fManagedObject),
  fOwnership(false)
{

}

AliTrackContainer::TrackOwnerHandler &AliTrackContainer::TrackOwnerHandler::operator=(const AliTrackContainer::TrackOwnerHandler &other) {
  TObject::operator=(other);
  if(this != &other) {
    if(fOwnership && fManagedObject) delete fManagedObject;
    fManagedObject = other.fManagedObject;
    fOwnership = false;
  }
  return *this;
}

AliTrackContainer::TrackOwnerHandler::~TrackOwnerHandler() {
  if(fOwnership && fManagedObject) delete fManagedObject;
}

void AliTrackContainer::TrackOwnerHandler::SetObject(TObjArray *obj) {
  if(fOwnership && fManagedObject) delete fManagedObject;
  fManagedObject = obj;
}

void AliTrackContainer::TrackOwnerHandler::SetOwner(Bool_t owner) {
  fOwnership = owner;  
}

void AliTrackContainer::TrackOwnerHandler::TransferOwnershipTo(AliTrackContainer::TrackOwnerHandler &target) {
  if(fOwnership) target.SetOwner(fOwnership);
  fOwnership = false;
}

void AliTrackContainer::TrackOwnerHandler::ReceiveOwnershipFrom(AliTrackContainer::TrackOwnerHandler &source) {
  if(source.IsOwner()) {
    fOwnership = true;
    source.SetOwner(false);
  }
}
