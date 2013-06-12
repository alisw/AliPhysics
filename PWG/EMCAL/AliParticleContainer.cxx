//
// Particle Container
//
// Author: M. Verweij

#include <TROOT.h>
#include <TSystem.h>
#include <TInterpreter.h>

#include <TChain.h>
#include <TClonesArray.h>
#include <TObject.h>
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

}

//________________________________________________________________________
AliParticleContainer::~AliParticleContainer()
{
  // Destructor.
}

//________________________________________________________________________
void AliParticleContainer::SetParticleArray(AliVEvent *event) 
{
  // Set jet array

  SetArray(event, "AliVParticle");

}

//________________________________________________________________________
AliVParticle* AliParticleContainer::GetParticle(Int_t i) const {

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
