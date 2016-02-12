//
// Container with name, TClonesArray and cuts for particles
//
// Author: M. Verweij, S. Aiola

#include <TClonesArray.h>

#include "AliVEvent.h"
#include "AliLog.h"
#include "AliVCuts.h"
#include "AliVTrack.h"
#include "AliESDtrack.h"

#include "AliTLorentzVector.h"
#include "AliEmcalTrackSelectionAOD.h"
#include "AliEmcalTrackSelectionESD.h"
#include "AliParticleContainer.h"

ClassImp(AliParticleContainer)


TString AliParticleContainer::fgDefTrackCutsPeriod = "";

//________________________________________________________________________
AliParticleContainer::AliParticleContainer():
  AliEmcalContainer(),
  fParticlePtCut(0.15),
  fParticleMinEta(-0.9),
  fParticleMaxEta(0.9),
  fParticleMinPhi(-10),
  fParticleMaxPhi(10),
  fPhiOffset(0.),
  fMinDistanceTPCSectorEdge(-1),
  fMinMCLabel(-1),
  fMaxMCLabel(-1),
  fMCFlag(0),
  fGeneratorIndex(-1),
  fCharge(-1),
  fTrackFilterType(AliEmcalTrackSelection::kNoTrackFilter),
  fListOfCuts(0),
  fSelectionModeAny(kFALSE),
  fAODFilterBits(0),
  fTrackCutsPeriod(),
  fEmcalTrackSelection(0),
  fFilteredTracks(0),
  fTrackTypes(5000)
{
  // Default constructor.

  fClassName = "AliVParticle";
}

//________________________________________________________________________
AliParticleContainer::AliParticleContainer(const char *name, const char *period):
  AliEmcalContainer(name),
  fParticlePtCut(0.15),
  fParticleMinEta(-0.9),
  fParticleMaxEta(0.9),
  fParticleMinPhi(-10),
  fParticleMaxPhi(10),
  fPhiOffset(0.),
  fMinDistanceTPCSectorEdge(-1),
  fMinMCLabel(-1),
  fMaxMCLabel(-1),
  fMCFlag(0),
  fGeneratorIndex(-1),
  fCharge(-1),
  fTrackFilterType(AliEmcalTrackSelection::kNoTrackFilter),
  fListOfCuts(0),
  fSelectionModeAny(kFALSE),
  fAODFilterBits(0),
  fTrackCutsPeriod(period),
  fEmcalTrackSelection(0),
  fFilteredTracks(0),
  fTrackTypes(5000)
{
  // Standard constructor.

  fClassName = "AliVParticle";

  if (fTrackCutsPeriod.IsNull() && !AliParticleContainer::fgDefTrackCutsPeriod.IsNull()) {
    AliInfo(Form("Default track cuts period is %s", AliParticleContainer::fgDefTrackCutsPeriod.Data()));
    fTrackCutsPeriod = AliParticleContainer::fgDefTrackCutsPeriod;
  }
}

//________________________________________________________________________
void AliParticleContainer::SetArray(AliVEvent *event)
{
  // Get array from event.

  AliEmcalContainer::SetArray(event);

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
void AliParticleContainer::NextEvent()
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
AliVParticle* AliParticleContainer::GetLeadingParticle(const char* opt) 
{
  // Get the leading particle; use p if "p" is contained in opt

  TString option(opt);
  option.ToLower();

  Int_t tempID = fCurrentID;

  AliVParticle *partMax = GetNextAcceptParticle(0);
  AliVParticle *part = 0;

  if (option.Contains("p")) {
    while ((part = GetNextAcceptParticle())) {
      if (part->P() > partMax->P()) partMax = part;
    }
  }
  else {
    while ((part = GetNextAcceptParticle())) {
      if (part->Pt() > partMax->Pt()) partMax = part;
    }
  }

  fCurrentID = tempID;

  return partMax;
}

//________________________________________________________________________
AliVParticle* AliParticleContainer::GetParticle(Int_t i) const 
{
  //Get i^th jet in array

  if (i == -1) i = fCurrentID;
  if (i < 0 || i >= fFilteredTracks->GetEntriesFast()) return 0;
  AliVParticle *vp = static_cast<AliVParticle*>(fFilteredTracks->At(i));
  return vp;
}

//________________________________________________________________________
AliVParticle* AliParticleContainer::GetAcceptParticle(Int_t i)
{
  //return pointer to particle if particle is accepted

  if (i == -1) i = fCurrentID;
  if (AcceptParticle(i)) {
      return GetParticle(i);
  }
  else {
    AliDebug(2,"Particle not accepted.");
    return 0;
  }
}

//________________________________________________________________________
AliVParticle* AliParticleContainer::GetParticleWithLabel(Int_t lab) const 
{
  //Get particle with label lab in array
  
  Int_t i = GetIndexFromLabel(lab);
  return GetParticle(i);
}

//________________________________________________________________________
AliVParticle* AliParticleContainer::GetAcceptParticleWithLabel(Int_t lab)  
{
  //Get particle with label lab in array
  
  Int_t i = GetIndexFromLabel(lab);
  return GetAcceptParticle(i);
}


//________________________________________________________________________
AliVParticle* AliParticleContainer::GetNextAcceptParticle(Int_t i) 
{
  //Get next accepted particle; if i >= 0 (re)start counter from i; return 0 if no accepted particle could be found

  if (i >= 0) fCurrentID = i;

  const Int_t n = GetNEntries();

  AliVParticle *p = 0;
  do {
    fCurrentID++;
    if (fCurrentID >= n) break;
    p = GetAcceptParticle(fCurrentID);
  } while (!p);

  return p;
}

//________________________________________________________________________
AliVParticle* AliParticleContainer::GetNextParticle(Int_t i) 
{
  //Get next particle; if i >= 0 (re)start counter from i; return 0 if no particle could be found

  if (i >= 0) fCurrentID = i;

  const Int_t n = GetNEntries();
  AliVParticle *p = 0;
  do {
    fCurrentID++;
    if (fCurrentID >= n) break;
    p = GetParticle(fCurrentID);
  } while (!p);

  return p;
}

//________________________________________________________________________
Bool_t AliParticleContainer::GetMomentum(TLorentzVector &mom, Int_t i)
{
  //Get momentum of the i^th particle in array

  if (i == -1) i = fCurrentID;
  AliVParticle *vp = GetParticle(i);
  if (vp) {
    if (fLoadedClass->InheritsFrom("AliESDtrack") && fTrackFilterType == AliEmcalTrackSelection::kHybridTracks &&
        (fTrackTypes[i] == kHybridConstrained || fTrackTypes[i] == kHybridConstrainedNoITSrefit)) {
      AliESDtrack *track = static_cast<AliESDtrack*>(vp);
      mom.SetPtEtaPhiM(track->GetConstrainedParam()->Pt(), track->GetConstrainedParam()->Eta(), track->GetConstrainedParam()->Phi(), 0.139);
    }
    else {
      mom.SetPtEtaPhiM(vp->Pt(), vp->Eta(), vp->Phi(), 0.139);
    }
    return kTRUE;
  }
  else {
    mom.SetPtEtaPhiM(0, 0, 0, 0.139);
    return kFALSE;
  }
}

//________________________________________________________________________
Bool_t AliParticleContainer::GetNextMomentum(TLorentzVector &mom, Int_t i)
{
  //Get momentum of the i^th particle in array

  AliVParticle *vp = GetNextParticle(i);
  if (vp) {
    if (fLoadedClass->InheritsFrom("AliESDtrack") && fTrackFilterType == AliEmcalTrackSelection::kHybridTracks &&
        (fTrackTypes[fCurrentID] == kHybridConstrained || fTrackTypes[fCurrentID] == kHybridConstrainedNoITSrefit)) {
      AliESDtrack *track = static_cast<AliESDtrack*>(vp);
      mom.SetPtEtaPhiM(track->GetConstrainedParam()->Pt(), track->GetConstrainedParam()->Eta(), track->GetConstrainedParam()->Phi(), 0.139);
    }
    else {
      mom.SetPtEtaPhiM(vp->Pt(), vp->Eta(), vp->Phi(), 0.139);
    }
    return kTRUE;
  }
  else {
    mom.SetPtEtaPhiM(0, 0, 0, 0.139);
    return kFALSE;
  }
}

//________________________________________________________________________
Bool_t AliParticleContainer::GetAcceptMomentum(TLorentzVector &mom, Int_t i)
{
  //Get momentum of the i^th particle in array

  if (i == -1) i = fCurrentID;
  AliVParticle *vp = GetAcceptParticle(i);
  if (vp) {
    if (fLoadedClass->InheritsFrom("AliESDtrack") && fTrackFilterType == AliEmcalTrackSelection::kHybridTracks &&
        (fTrackTypes[i] == kHybridConstrained || fTrackTypes[i] == kHybridConstrainedNoITSrefit)) {
      AliESDtrack *track = static_cast<AliESDtrack*>(vp);
      mom.SetPtEtaPhiM(track->GetConstrainedParam()->Pt(), track->GetConstrainedParam()->Eta(), track->GetConstrainedParam()->Phi(), 0.139);
    }
    else {
      mom.SetPtEtaPhiM(vp->Pt(), vp->Eta(), vp->Phi(), 0.139);
    }

    return kTRUE;
  }
  else {
    mom.SetPtEtaPhiM(0, 0, 0, 0.139);
    return kFALSE;
  }
}

//________________________________________________________________________
Bool_t AliParticleContainer::GetNextAcceptMomentum(TLorentzVector &mom, Int_t i)
{
  //Get momentum of the i^th particle in array

  AliVParticle *vp = GetNextAcceptParticle(i);
  if (vp) {
    if (fLoadedClass->InheritsFrom("AliESDtrack") && fTrackFilterType == AliEmcalTrackSelection::kHybridTracks &&
        (fTrackTypes[fCurrentID] == kHybridConstrained || fTrackTypes[fCurrentID] == kHybridConstrainedNoITSrefit)) {
      AliESDtrack *track = static_cast<AliESDtrack*>(vp);
      mom.SetPtEtaPhiM(track->GetConstrainedParam()->Pt(), track->GetConstrainedParam()->Eta(), track->GetConstrainedParam()->Phi(), 0.139);
    }
    else {
      mom.SetPtEtaPhiM(vp->Pt(), vp->Eta(), vp->Phi(), 0.139);
    }

    return kTRUE;
  }
  else {
    mom.SetPtEtaPhiM(0, 0, 0, 0.139);
    return kFALSE;
  }
}

//________________________________________________________________________
Bool_t AliParticleContainer::AcceptParticle(AliVParticle *vp)
{
  // Return true if vp is accepted.
  Int_t id = fFilteredTracks->IndexOf(vp);
  if (id >= 0) {
    return AcceptParticle(id);
  }
  else {
    return kFALSE;
  }
}

//________________________________________________________________________
Bool_t AliParticleContainer::AcceptParticle(Int_t i)
{
  // Return true if i^th particle is accepted.

  fRejectionReason = 0;

  // Cuts on the particle properties
  AliVParticle* vp = GetParticle(i);

  if (!vp) {
    fRejectionReason |= kNullObject;
    return kFALSE;
  }

  if (vp->TestBits(fBitMap) != (Int_t)fBitMap) {
    fRejectionReason |= kBitMapCut;
    return kFALSE;
  }

  if (fMinMCLabel >= 0 && TMath::Abs(vp->GetLabel()) > fMinMCLabel) {
    fRejectionReason |= kMCLabelCut;
    return kFALSE;
  }

  if (fMaxMCLabel >= 0 && TMath::Abs(vp->GetLabel()) < fMaxMCLabel) {
    fRejectionReason |= kMCLabelCut;
    return kFALSE;
  }

  if ((vp->GetFlag() & fMCFlag) != fMCFlag) {
    fRejectionReason |= kMCFlag;
    return kFALSE;
  }

  if (fGeneratorIndex >= 0 && fGeneratorIndex != vp->GetGeneratorIndex()) {
    fRejectionReason |= kMCGeneratorCut;
    return kFALSE;
  }

  if (fCharge>=0 && fCharge != vp->Charge()) {
    fRejectionReason |= kChargeCut;
    return kFALSE;
  }

  // Cuts on the 4-momentum
  AliTLorentzVector mom;
  GetMomentum(mom, i);

  if (mom.Pt() < fParticlePtCut) {
    fRejectionReason |= kPtCut;
    return kFALSE;
  }

  Double_t phi = mom.Phi_0_2pi() + fPhiOffset;

  if (mom.Eta() < fParticleMinEta || mom.Eta() > fParticleMaxEta ||
      phi < fParticleMinPhi       || phi > fParticleMaxPhi) {
    fRejectionReason |= kAcceptanceCut;
    return kFALSE;
  }

  if(fMinDistanceTPCSectorEdge>0.) {
    const Double_t pi = TMath::Pi();
    const Double_t kSector = pi/9;
    Double_t phiDist = TMath::Abs(mom.Phi() - TMath::FloorNint(mom.Phi()/kSector)*kSector);
    if(phiDist<fMinDistanceTPCSectorEdge) {
      fRejectionReason |= kMinDistanceTPCSectorEdgeCut;
      return kFALSE;
    }
  }
  
  return kTRUE;
}

//________________________________________________________________________
Int_t AliParticleContainer::GetNAcceptedParticles()
{
  // Get number of accepted particles

  Int_t nPart = 0;

  AliVParticle *vp = GetNextAcceptParticle(0);
  if(vp) nPart = 1;
  while (GetNextAcceptParticle())
    nPart++;

  return nPart;
}

//________________________________________________________________________
void AliParticleContainer::SetClassName(const char *clname)
{
  // Set the class name

  TClass cls(clname);
  if (cls.InheritsFrom("AliVParticle")) fClassName = clname;
  else AliError(Form("Unable to set class name %s for a AliParticleContainer, it must inherits from AliVParticle!",clname));
}

//________________________________________________________________________
void AliParticleContainer::AddTrackCuts(AliVCuts *cuts)
{
  if (!fListOfCuts) {
    fListOfCuts = new TObjArray;
    fListOfCuts->SetOwner(true);
  }
  fListOfCuts->Add(cuts);
}

//________________________________________________________________________
Int_t AliParticleContainer::GetNumberOfCutObjects() const
{
  if (!fListOfCuts) return 0;
  return fListOfCuts->GetEntries();
}

//________________________________________________________________________
AliVCuts* AliParticleContainer::GetTrackCuts(Int_t icut)
{
  if (!fListOfCuts) return NULL;
  if (icut < fListOfCuts->GetEntries()) {
    return static_cast<AliVCuts *>(fListOfCuts->At(icut));
  }
  return NULL;
}

//________________________________________________________________________
const char* AliParticleContainer::GetTitle() const
{
  static TString trackString;

  if (GetParticlePtCut() == 0) {
    trackString = TString::Format("%s_pT0000", GetArrayName().Data());
  }
  else if (GetParticlePtCut() < 1.0) {
    trackString = TString::Format("%s_pT0%3.0f", GetArrayName().Data(), GetParticlePtCut()*1000.0);
  }
  else {
    trackString = TString::Format("%s_pT%4.0f", GetArrayName().Data(), GetParticlePtCut()*1000.0);
  }

  return trackString.Data();
}
