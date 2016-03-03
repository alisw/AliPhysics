//
// Container with name, TClonesArray and cuts for particles
//
// Author: M. Verweij, S. Aiola

#include <TClonesArray.h>

#include "AliVEvent.h"
#include "AliLog.h"
#include "AliVCuts.h"
#include "AliESDtrack.h"

#include "AliTLorentzVector.h"
#include "AliEmcalTrackSelectionAOD.h"
#include "AliEmcalTrackSelectionESD.h"
#include "AliTrackContainer.h"

ClassImp(AliTrackContainer)

TString AliTrackContainer::fgDefTrackCutsPeriod = "";

//________________________________________________________________________
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
  // Default constructor.

  fClassName = "AliVTrack";
  fMassHypothesis = 0.139;
}

//________________________________________________________________________
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
  // Standard constructor.

  fClassName = "AliVTrack";

  if (fTrackCutsPeriod.IsNull() && !AliTrackContainer::fgDefTrackCutsPeriod.IsNull()) {
    AliInfo(Form("Default track cuts period is %s", AliTrackContainer::fgDefTrackCutsPeriod.Data()));
    fTrackCutsPeriod = AliTrackContainer::fgDefTrackCutsPeriod;
  }
  fMassHypothesis = 0.139;
}

//________________________________________________________________________
void AliTrackContainer::SetArray(AliVEvent *event)
{
  // Get array from event.

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

//________________________________________________________________________
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

//________________________________________________________________________
AliVTrack* AliTrackContainer::GetTrack(Int_t i) const
{
  //Get i^th jet in array

  if (i == -1) i = fCurrentID;
  if (i < 0 || i >= fFilteredTracks->GetEntriesFast()) return 0;
  AliVTrack *vp = static_cast<AliVTrack*>(fFilteredTracks->At(i));
  return vp;
}

//________________________________________________________________________
AliVTrack* AliTrackContainer::GetAcceptTrack(Int_t i)
{
  //return pointer to particle if particle is accepted

  if (i == -1) i = fCurrentID;
  if (AcceptTrack(i)) {
      return GetTrack(i);
  }
  else {
    AliDebug(2,"Track not accepted.");
    return 0;
  }
}

//________________________________________________________________________
AliVTrack* AliTrackContainer::GetNextAcceptTrack()
{
  //Get next accepted particle

  const Int_t n = GetNEntries();
  AliVTrack *p = 0;
  do {
    fCurrentID++;
    if (fCurrentID >= n) break;
    p = GetAcceptTrack(fCurrentID);
  } while (!p);

  return p;
}

//________________________________________________________________________
AliVTrack* AliTrackContainer::GetNextTrack()
{
  //Get next particle

  const Int_t n = GetNEntries();
  AliVTrack *p = 0;
  do {
    fCurrentID++;
    if (fCurrentID >= n) break;
    p = GetTrack(fCurrentID);
  } while (!p);

  return p;
}

//________________________________________________________________________
Bool_t AliTrackContainer::GetMomentum(TLorentzVector &mom, const AliVTrack* part, Double_t mass)
{
  if (part) {
    if (mass < 0) mass = part->M();
    mom.SetPtEtaPhiM(part->Pt(), part->Eta(), part->Phi(), mass);
    return kTRUE;
  }
  else {
    mom.SetPtEtaPhiM(0, 0, 0, 0);
    return kFALSE;
  }
}

//________________________________________________________________________
Bool_t AliTrackContainer::GetMomentum(TLorentzVector &mom, const AliVTrack* part)
{
  return GetMomentum(mom,part,fMassHypothesis);
}

//________________________________________________________________________
Bool_t AliTrackContainer::GetMomentum(TLorentzVector &mom, Int_t i)
{
  //Get momentum of the i^th particle in array

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

//________________________________________________________________________
Bool_t AliTrackContainer::GetNextMomentum(TLorentzVector &mom)
{
  //Get momentum of the next particle in array

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

//________________________________________________________________________
Bool_t AliTrackContainer::GetAcceptMomentum(TLorentzVector &mom, Int_t i)
{
  //Get momentum of the i^th particle in array

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

//________________________________________________________________________
Bool_t AliTrackContainer::GetNextAcceptMomentum(TLorentzVector &mom)
{
  //Get momentum of the next accepted particle in array

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

//________________________________________________________________________
Bool_t AliTrackContainer::AcceptTrack(const AliVTrack *vp)
{
  // Return true if vp is accepted.
  Bool_t r = ApplyTrackCuts(vp);
  if (!r) return kFALSE;

  AliTLorentzVector mom;

  Int_t id = fFilteredTracks->IndexOf(vp);
  if (id >= 0) {
    GetMomentum(mom, id);
  }
  else {
    GetMomentum(mom, vp);
  }

  return ApplyKinematicCuts(mom);
}

//________________________________________________________________________
Bool_t AliTrackContainer::AcceptTrack(Int_t i)
{
  // Return true if vp is accepted.
  Bool_t r = ApplyTrackCuts(GetTrack(i));
  if (!r) return kFALSE;

  AliTLorentzVector mom;
  GetMomentum(mom, i);

  return ApplyKinematicCuts(mom);
}

//________________________________________________________________________
Bool_t AliTrackContainer::ApplyTrackCuts(const AliVTrack* vp)
{
  // Return true if i^th particle is accepted.

  return ApplyParticleCuts(vp);
}

//________________________________________________________________________
void AliTrackContainer::SetClassName(const char *clname)
{
  // Set the class name

  TClass cls(clname);
  if (cls.InheritsFrom("AliVTrack")) fClassName = clname;
  else AliError(Form("Unable to set class name %s for a AliTrackContainer, it must inherits from AliVTrack!",clname));
}

//________________________________________________________________________
void AliTrackContainer::AddTrackCuts(AliVCuts *cuts)
{
  if (!fListOfCuts) {
    fListOfCuts = new TObjArray;
    fListOfCuts->SetOwner(true);
  }
  fListOfCuts->Add(cuts);
}

//________________________________________________________________________
Int_t AliTrackContainer::GetNumberOfCutObjects() const
{
  if (!fListOfCuts) return 0;
  return fListOfCuts->GetEntries();
}

//________________________________________________________________________
AliVCuts* AliTrackContainer::GetTrackCuts(Int_t icut)
{
  if (!fListOfCuts) return NULL;
  if (icut < fListOfCuts->GetEntries()) {
    return static_cast<AliVCuts *>(fListOfCuts->At(icut));
  }
  return NULL;
}

//________________________________________________________________________
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
