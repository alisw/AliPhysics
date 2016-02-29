//
// Container with name, TClonesArray and cuts for particles
//
// Author: M. Verweij, S. Aiola

#include <TClonesArray.h>

#include "AliVEvent.h"
#include "AliLog.h"

#include "AliTLorentzVector.h"
#include "AliParticleContainer.h"

ClassImp(AliParticleContainer)

//________________________________________________________________________
AliParticleContainer::AliParticleContainer():
  AliEmcalContainer(),
  fMinDistanceTPCSectorEdge(-1),
  fCharge(-1),
  fGeneratorIndex(-1)
{
  // Default constructor.

  fClassName = "AliVParticle";
}

//________________________________________________________________________
AliParticleContainer::AliParticleContainer(const char *name) :
  AliEmcalContainer(name),
  fMinDistanceTPCSectorEdge(-1),
  fCharge(-1),
  fGeneratorIndex(-1)
{
  // Standard constructor.

  fClassName = "AliVParticle";
}

//________________________________________________________________________
AliVParticle* AliParticleContainer::GetLeadingParticle(const char* opt) 
{
  // Get the leading particle; use p if "p" is contained in opt

  TString option(opt);
  option.ToLower();

  Int_t tempID = fCurrentID;
  ResetCurrentID();

  AliVParticle *partMax = GetNextAcceptParticle();
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
  if (i < 0 || i >= this->fClArray->GetEntriesFast()) return 0;
  AliVParticle *vp = static_cast<AliVParticle*>(fClArray->At(i));
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
AliVParticle* AliParticleContainer::GetNextAcceptParticle()
{
  //Get next accepted particle

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
AliVParticle* AliParticleContainer::GetNextParticle()
{
  //Get next particle

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
Bool_t AliParticleContainer::GetMomentum(TLorentzVector &mom, const AliVParticle* part, Double_t mass)
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
Bool_t AliParticleContainer::GetMomentum(TLorentzVector &mom, const AliVParticle* part)
{
  return GetMomentum(mom,part,fMassHypothesis);
}

//________________________________________________________________________
Bool_t AliParticleContainer::GetMomentum(TLorentzVector &mom, Int_t i)
{
  //Get momentum of the i^th particle in array

  Double_t mass = fMassHypothesis;

  if (i == -1) i = fCurrentID;
  AliVParticle *vp = GetParticle(i);
  if (vp) {
    if (mass < 0) mass = vp->M();
    mom.SetPtEtaPhiM(vp->Pt(), vp->Eta(), vp->Phi(), mass);
    return kTRUE;
  }
  else {
    mom.SetPtEtaPhiM(0, 0, 0, 0);
    return kFALSE;
  }
}

//________________________________________________________________________
Bool_t AliParticleContainer::GetNextMomentum(TLorentzVector &mom)
{
  //Get momentum of the next particle in array

  Double_t mass = fMassHypothesis;

  AliVParticle *vp = GetNextParticle();
  if (vp) {
    if (mass < 0) mass = vp->M();
    mom.SetPtEtaPhiM(vp->Pt(), vp->Eta(), vp->Phi(), mass);
    return kTRUE;
  }
  else {
    mom.SetPtEtaPhiM(0, 0, 0, 0);
    return kFALSE;
  }
}

//________________________________________________________________________
Bool_t AliParticleContainer::GetAcceptMomentum(TLorentzVector &mom, Int_t i)
{
  //Get momentum of the i^th particle in array

  Double_t mass = fMassHypothesis;

  if (i == -1) i = fCurrentID;
  AliVParticle *vp = GetAcceptParticle(i);
  if (vp) {
    if (mass < 0) mass = vp->M();
    mom.SetPtEtaPhiM(vp->Pt(), vp->Eta(), vp->Phi(), mass);
    return kTRUE;
  }
  else {
    mom.SetPtEtaPhiM(0, 0, 0, 0);
    return kFALSE;
  }
}

//________________________________________________________________________
Bool_t AliParticleContainer::GetNextAcceptMomentum(TLorentzVector &mom)
{
  //Get momentum of the next accepted particle in array

  Double_t mass = fMassHypothesis;

  AliVParticle *vp = GetNextAcceptParticle();
  if (vp) {
    if (mass < 0) mass = vp->M();
    mom.SetPtEtaPhiM(vp->Pt(), vp->Eta(), vp->Phi(), mass);
    return kTRUE;
  }
  else {
    mom.SetPtEtaPhiM(0, 0, 0, 0);
    return kFALSE;
  }
}

//________________________________________________________________________
Bool_t AliParticleContainer::AcceptParticle(const AliVParticle *vp)
{
  // Return true if vp is accepted.
  Bool_t r = ApplyParticleCuts(vp);
  if (!r) return kFALSE;

  AliTLorentzVector mom;

  Int_t id = fClArray->IndexOf(vp);
  if (id >= 0) {
    GetMomentum(mom, id);
  }
  else {
    GetMomentum(mom, vp);
  }

  return ApplyKinematicCuts(mom);
}

//________________________________________________________________________
Bool_t AliParticleContainer::AcceptParticle(Int_t i)
{
  // Return true if vp is accepted.
  Bool_t r = ApplyParticleCuts(GetParticle(i));
  if (!r) return kFALSE;

  AliTLorentzVector mom;
  GetMomentum(mom, i);

  return ApplyKinematicCuts(mom);
}

//________________________________________________________________________
Bool_t AliParticleContainer::ApplyParticleCuts(const AliVParticle* vp)
{
  // Return true if i^th particle is accepted.

  fRejectionReason = 0;

  // Cuts on the particle properties

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

  if (fCharge>=0 && fCharge != vp->Charge()) {
    fRejectionReason |= kChargeCut;
    return kFALSE;
  }

  if (fGeneratorIndex >= 0 && fGeneratorIndex != vp->GetGeneratorIndex()) {
    fRejectionReason |= kMCGeneratorCut;
    return kFALSE;
  }

  return kTRUE;
}

//________________________________________________________________________
Bool_t AliParticleContainer::ApplyKinematicCuts(const AliTLorentzVector& mom)
{
  if(fMinDistanceTPCSectorEdge>0.) {
    const Double_t pi = TMath::Pi();
    const Double_t kSector = pi/9;
    Double_t phiDist = TMath::Abs(mom.Phi() - TMath::FloorNint(mom.Phi()/kSector)*kSector);
    if(phiDist<fMinDistanceTPCSectorEdge) {
      fRejectionReason |= kMinDistanceTPCSectorEdgeCut;
      return kFALSE;
    }
  }

  return AliEmcalContainer::ApplyKinematicCuts(mom);
}

//________________________________________________________________________
Int_t AliParticleContainer::GetNAcceptedParticles()
{
  // Get number of accepted particles

  Int_t nPart = 0;
  Int_t tempID = fCurrentID;
  ResetCurrentID();

  AliVParticle *vp = GetNextAcceptParticle();
  if(vp) nPart = 1;
  while (GetNextAcceptParticle())
    nPart++;

  fCurrentID = tempID;

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
const char* AliParticleContainer::GetTitle() const
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
