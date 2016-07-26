//
// Container with name, TClonesArray and cuts for particles
//
// Author: M. Verweij, S. Aiola

#include <TClonesArray.h>

#include "AliVEvent.h"
#include "AliLog.h"

#include "AliTLorentzVector.h"
#include "AliMCParticleContainer.h"

ClassImp(AliParticleContainer)

//________________________________________________________________________
AliMCParticleContainer::AliMCParticleContainer():
  AliParticleContainer(),
  fMCFlag(AliAODMCParticle::kPhysicalPrim)
{
  // Default constructor.

  fClassName = "AliAODMCParticle";
}

//________________________________________________________________________
AliMCParticleContainer::AliMCParticleContainer(const char *name):
  AliParticleContainer(name),
  fMCFlag(AliAODMCParticle::kPhysicalPrim)
{
  // Standard constructor.

  fClassName = "AliAODMCParticle";
}

//________________________________________________________________________
AliAODMCParticle* AliMCParticleContainer::GetMCParticleWithLabel(Int_t lab) const
{
  //Get particle with label lab in array

  Int_t i = GetIndexFromLabel(lab);
  return GetMCParticle(i);
}

//________________________________________________________________________
AliAODMCParticle* AliMCParticleContainer::GetAcceptMCParticleWithLabel(Int_t lab)
{
  //Get particle with label lab in array

  Int_t i = GetIndexFromLabel(lab);
  return GetAcceptMCParticle(i);
}

//________________________________________________________________________
AliAODMCParticle* AliMCParticleContainer::GetMCParticle(Int_t i) const
{
  //Get i^th jet in array

  if (i == -1) i = fCurrentID;
  if (i < 0 || i >= fClArray->GetEntriesFast()) return 0;
  AliAODMCParticle *vp = static_cast<AliAODMCParticle*>(fClArray->At(i));
  return vp;
}

//________________________________________________________________________
AliAODMCParticle* AliMCParticleContainer::GetAcceptMCParticle(Int_t i)
{
  //return pointer to particle if particle is accepted

  if (i == -1) i = fCurrentID;
  if (AcceptMCParticle(i)) {
      return GetMCParticle(i);
  }
  else {
    AliDebug(2,"Particle not accepted.");
    return 0;
  }
}

//________________________________________________________________________
AliAODMCParticle* AliMCParticleContainer::GetNextAcceptMCParticle()
{
  //Get next accepted particle

  const Int_t n = GetNEntries();
  AliAODMCParticle *p = 0;
  do {
    fCurrentID++;
    if (fCurrentID >= n) break;
    p = GetAcceptMCParticle(fCurrentID);
  } while (!p);

  return p;
}

//________________________________________________________________________
AliAODMCParticle* AliMCParticleContainer::GetNextMCParticle()
{
  //Get next particle

  const Int_t n = GetNEntries();
  AliAODMCParticle *p = 0;
  do {
    fCurrentID++;
    if (fCurrentID >= n) break;
    p = GetMCParticle(fCurrentID);
  } while (!p);

  return p;
}

//________________________________________________________________________
Bool_t AliMCParticleContainer::GetMomentum(TLorentzVector &mom, const AliAODMCParticle* part, Double_t mass)
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
Bool_t AliMCParticleContainer::GetMomentum(TLorentzVector &mom, const AliAODMCParticle* part)
{
  return GetMomentum(mom,part,fMassHypothesis);
}

//________________________________________________________________________
Bool_t AliMCParticleContainer::GetMomentum(TLorentzVector &mom, Int_t i)
{
  //Get momentum of the i^th particle in array

  Double_t mass = fMassHypothesis;

  if (i == -1) i = fCurrentID;
  AliAODMCParticle *vp = GetMCParticle(i);
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
Bool_t AliMCParticleContainer::GetNextMomentum(TLorentzVector &mom)
{
  //Get momentum of the next particle in array

  Double_t mass = fMassHypothesis;

  AliAODMCParticle *vp = GetNextMCParticle();
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
Bool_t AliMCParticleContainer::GetAcceptMomentum(TLorentzVector &mom, Int_t i)
{
  //Get momentum of the i^th particle in array

  Double_t mass = fMassHypothesis;

  if (i == -1) i = fCurrentID;
  AliAODMCParticle *vp = GetAcceptMCParticle(i);
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
Bool_t AliMCParticleContainer::GetNextAcceptMomentum(TLorentzVector &mom)
{
  //Get momentum of the next accepted particle in array

  Double_t mass = fMassHypothesis;

  AliAODMCParticle *vp = GetNextAcceptMCParticle();
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
Bool_t AliMCParticleContainer::AcceptMCParticle(const AliAODMCParticle *vp)
{
  // Return true if vp is accepted.
  Bool_t r = ApplyMCParticleCuts(vp);
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
Bool_t AliMCParticleContainer::AcceptMCParticle(Int_t i)
{
  // Return true if vp is accepted.
  Bool_t r = ApplyMCParticleCuts(GetMCParticle(i));
  if (!r) return kFALSE;

  AliTLorentzVector mom;
  GetMomentum(mom, i);

  return ApplyKinematicCuts(mom);
}

//________________________________________________________________________
Bool_t AliMCParticleContainer::ApplyMCParticleCuts(const AliAODMCParticle* vp)
{
  // Return true if i^th particle is accepted.

  fRejectionReason = 0;

  // Cuts on the particle properties

  if ((vp->GetFlag() & fMCFlag) != fMCFlag) {
    fRejectionReason |= kMCFlag;
    return kFALSE;
  }

  return ApplyParticleCuts(vp);
}

//________________________________________________________________________
void AliMCParticleContainer::SetClassName(const char *clname)
{
  // Set the class name

  TClass cls(clname);
  if (cls.InheritsFrom("AliAODMCParticle")) fClassName = clname;
  else AliError(Form("Unable to set class name %s for a AliMCParticleContainer, it must inherits from AliAODMCParticle!",clname));
}

//________________________________________________________________________
const char* AliMCParticleContainer::GetTitle() const
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
