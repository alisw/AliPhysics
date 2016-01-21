//
// Container with name, TClonesArray and cuts for particles
//
// Author: M. Verweij, S. Aiola

#include <TClonesArray.h>

#include "AliVEvent.h"
#include "AliLog.h"
#include "AliVCuts.h"
#include "AliVTrack.h"

#include "AliEmcalTrackSelectionAOD.h"
#include "AliEmcalTrackSelectionESD.h"
#include "AliParticleContainer.h"

ClassImp(AliParticleContainer)


TString AliParticleContainer::fgDefTrackCutsPeriod = "";

//________________________________________________________________________
AliParticleContainer::AliParticleContainer():
  AliEmcalContainer("AliParticleContainer"),
  fParticlePtCut(0.15),
  fParticleMinEta(-0.9),
  fParticleMaxEta(0.9),
  fParticleMinPhi(-10),
  fParticleMaxPhi(10),
  fPhiOffset(0.),
  fMinDistanceTPCSectorEdge(-1),
  fTrackBitMap(0),
  fMCTrackBitMap(0),
  fMinMCLabel(0),
  fMinMCLabelAccept(-1),
  fMCFlag(0),
  fGeneratorIndex(-1),
  fCharge(-1),
  fTrackFilterType(AliEmcalTrackSelection::kNoTrackFilter),
  fListOfCuts(0),
  fSelectionModeAny(kFALSE),
  fAODFilterBits(0),
  fTrackCutsPeriod(),
  fEmcalTrackSelection(0),
  fTrackType(kUndefined)
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
  fTrackBitMap(0),
  fMCTrackBitMap(0),
  fMinMCLabel(0),
  fMinMCLabelAccept(-1),
  fMCFlag(0),
  fGeneratorIndex(-1),
  fCharge(-1),
  fTrackFilterType(AliEmcalTrackSelection::kNoTrackFilter),
  fListOfCuts(0),
  fSelectionModeAny(kFALSE),
  fAODFilterBits(0),
  fTrackCutsPeriod(period),
  fEmcalTrackSelection(0),
  fTrackType(kUndefined)
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

  TClass* particleClass = fClArray->GetClass();

  if (fTrackFilterType == AliEmcalTrackSelection::kNoTrackFilter) {
    if (fEmcalTrackSelection) delete fEmcalTrackSelection;
    fEmcalTrackSelection = 0;
  }
  else {
    if (fTrackFilterType == AliEmcalTrackSelection::kCustomTrackFilter) {

      AliInfo("Using custom track cuts");

      if (particleClass->InheritsFrom("AliAODTrack")) {
        AliInfo(Form("Objects are of type %s: AOD track selection will be done.", particleClass->GetName()));
        fEmcalTrackSelection = new AliEmcalTrackSelectionAOD(0, fAODFilterBits);
      }
      else if (particleClass->InheritsFrom("AliESDtrack")) {
        AliInfo(Form("Objects are of type %s: ESD track selection will be done.", particleClass->GetName()));
        fEmcalTrackSelection = new AliEmcalTrackSelectionESD(0);
      }
      else {
        AliWarning(Form("Objects are of type %s: no track filtering will be done!!", particleClass->GetName()));
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

      if (particleClass->InheritsFrom("AliAODTrack")) {
        AliInfo(Form("Objects are of type %s: AOD track selection will be done.", particleClass->GetName()));
        fEmcalTrackSelection = new AliEmcalTrackSelectionAOD(fTrackFilterType, fTrackCutsPeriod);
      }
      else if (particleClass->InheritsFrom("AliESDtrack")) {
        AliInfo(Form("Objects are of type %s: ESD track selection will be done.", particleClass->GetName()));
        fEmcalTrackSelection = new AliEmcalTrackSelectionESD(fTrackFilterType, fTrackCutsPeriod);
      }
      else {
        AliWarning(Form("Objects are of type %s: no track filtering will be done!!", particleClass->GetName()));
      }
    }
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

  if(i<0 || i>=fClArray->GetEntriesFast()) return 0;
  AliVParticle *vp = static_cast<AliVParticle*>(fClArray->At(i));
  return vp;

}

//________________________________________________________________________
AliVParticle* AliParticleContainer::GetAcceptParticle(Int_t i) {
  //return pointer to particle if particle is accepted

  AliVParticle *vp = GetParticle(i);
  if(!vp) return 0;

  if(AcceptParticle(vp))
      return vp;
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

  if (i>=0) fCurrentID = i;

  const Int_t n = GetNEntries();
  AliVParticle *p = 0;
  while (fCurrentID < n && !p) { 
    p = GetAcceptParticle(fCurrentID);
    fCurrentID++;
  }

  return p;
}

//________________________________________________________________________
AliVParticle* AliParticleContainer::GetNextParticle(Int_t i) 
{
  //Get next particle; if i >= 0 (re)start counter from i; return 0 if no particle could be found

  if (i>=0) fCurrentID = i;

  const Int_t n = GetNEntries();
  AliVParticle *p = 0;
  while (fCurrentID < n && !p) { 
    p = GetParticle(fCurrentID);
    fCurrentID++;
  }

  return p;
}

//________________________________________________________________________
void AliParticleContainer::GetMomentum(TLorentzVector &mom, Int_t i) const
{
  //Get momentum of the i^th particle in array

  AliVParticle *vp = GetParticle(i);
  if (vp) mom.SetPtEtaPhiM(vp->Pt(),vp->Eta(),vp->Phi(),0.139);
}

//________________________________________________________________________
Bool_t AliParticleContainer::AcceptParticle(AliVParticle *vp)
{
  // Return true if vp is accepted.

  fRejectionReason = 0;
  fTrackType = kUndefined;

  if (!vp) {
    fRejectionReason |= kNullObject;
    return kFALSE;
  }

  if (fEmcalTrackSelection) {
    AliVTrack* vTrack = static_cast<AliVTrack*>(vp);
    if (!fEmcalTrackSelection->IsTrackAccepted(vTrack)) {
      fRejectionReason |= kNotHybridTrack;
      return kFALSE;
    }
    else {
      if (fTrackFilterType == AliEmcalTrackSelection::kHybridTracks) {
        if (fEmcalTrackSelection->GetTrackBitmap().FirstSetBit() == 0) {
          fTrackType = kHybridGlobal;
        }
        else if (fEmcalTrackSelection->GetTrackBitmap().FirstSetBit() == 1) {
          if (vTrack->GetStatus()&AliVTrack::kITSrefit != 0) {
            fTrackType = kHybridConstrained;
          }
          else {
            fTrackType = kHybridConstrainedNoITSrefit;
          }
        }
      }
    }
  }

  if (vp->Pt() < fParticlePtCut) {
    fRejectionReason |= kPtCut;
    return kFALSE;
  }

  Double_t phi = vp->Phi() + fPhiOffset;
  Double_t tpi = TMath::TwoPi();
  if(phi<0.)  phi+=tpi;
  if(phi>tpi) phi-=tpi;

  if (vp->Eta() < fParticleMinEta || vp->Eta() > fParticleMaxEta || 
      phi < fParticleMinPhi       || phi > fParticleMaxPhi) {
    fRejectionReason |= kAcceptanceCut;
    return kFALSE;
  }

  if(fMinDistanceTPCSectorEdge>0.) {
    const Double_t pi = TMath::Pi();
    const Double_t kSector = pi/9;
    Double_t phiDist = TMath::Abs(vp->Phi() - TMath::FloorNint(vp->Phi()/kSector)*kSector);
    if(phiDist<fMinDistanceTPCSectorEdge) {
      fRejectionReason |= kMinDistanceTPCSectorEdgeCut;
      return kFALSE;
    }
  }

  if (TMath::Abs(vp->GetLabel()) < fMinMCLabelAccept) {
    AliDebug(2,"Particle not accepted because label too small.");
    fRejectionReason |= kMinMCLabelAccept;
    return kFALSE;
  }

  if (TMath::Abs(vp->GetLabel()) > fMinMCLabel) {
    if(vp->TestBits(fMCTrackBitMap) != (Int_t)fMCTrackBitMap) {
      AliDebug(2,"MC particle not accepted because of MC bit map.");
      fRejectionReason |= kBitMapCut;
      return kFALSE;
    }
  }
  else {
    if(vp->TestBits(fTrackBitMap) != (Int_t)fTrackBitMap) {
      AliDebug(2,"Track not accepted because of bit map.");
      fRejectionReason |= kBitMapCut;
      return kFALSE;
    }
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
