// $Id$
//
// Container with name, TClonesArray and cuts for particles
//
// Author: M. Verweij

#include <TClonesArray.h>

#include "AliVEvent.h"
#include "AliLog.h"

#include "AliParticleContainer.h"

ClassImp(AliParticleContainer)

//________________________________________________________________________
AliParticleContainer::AliParticleContainer():
  AliEmcalContainer("AliParticleContainer"),
  fParticlePtCut(0.15),
  fParticleMinEta(-0.9),
  fParticleMaxEta(0.9),
  fParticleMinPhi(-10),
  fParticleMaxPhi(10),
  fTrackBitMap(0),
  fMCTrackBitMap(0),
  fMinMCLabel(0)
{
  // Default constructor.

  fClassName = "AliVParticle";
}

//________________________________________________________________________
AliParticleContainer::AliParticleContainer(const char *name):
  AliEmcalContainer(name),
  fParticlePtCut(0.15),
  fParticleMinEta(-0.9),
  fParticleMaxEta(0.9),
  fParticleMinPhi(-10),
  fParticleMaxPhi(10),
  fTrackBitMap(0),
  fMCTrackBitMap(0),
  fMinMCLabel(0)
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

  if(i<0 || i>fClArray->GetEntriesFast()) return 0;
  AliVParticle *vp = static_cast<AliVParticle*>(fClArray->At(i));
  return vp;

}

//________________________________________________________________________
AliVParticle* AliParticleContainer::GetAcceptParticle(Int_t i) const {
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
AliVParticle* AliParticleContainer::GetAcceptParticleWithLabel(Int_t lab) const 
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
  if(vp) mom.SetPtEtaPhiM(vp->Pt(),vp->Eta(),vp->Phi(),0.139);
}

//________________________________________________________________________
Bool_t AliParticleContainer::AcceptParticle(AliVParticle *vp) const
{
  // Return true if vp is accepted.

  if (!vp)
    return kFALSE;

  if (TMath::Abs(vp->GetLabel()) > fMinMCLabel) {
    if(vp->TestBits(fMCTrackBitMap) != (Int_t)fMCTrackBitMap) {
      AliDebug(2,"MC particle not accepted because of MC bit map.");
      return kFALSE;
    }
  }
  else {
    if(vp->TestBits(fTrackBitMap) != (Int_t)fTrackBitMap) {
      AliDebug(2,"Track not accepted because of bit map.");
      return kFALSE;
    }
  }

  if (vp->Pt() < fParticlePtCut)
    return kFALSE;

  if (vp->Eta() < fParticleMinEta || vp->Eta() > fParticleMaxEta || 
      vp->Phi() < fParticleMinPhi || vp->Phi() > fParticleMaxPhi)
    return kFALSE;
  
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
