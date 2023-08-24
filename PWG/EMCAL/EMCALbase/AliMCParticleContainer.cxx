
/************************************************************************************
 * Copyright (C) 2013, Copyright Holders of the ALICE Collaboration                 *
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
#include <iostream>

#include <TClonesArray.h>

#include "AliVEvent.h"
#include "AliLog.h"

#include "AliTLorentzVector.h"
#include "AliMCParticleContainer.h"

ClassImp(AliMCParticleContainer);

AliMCParticleContainer::AliMCParticleContainer():
  AliParticleContainer(),
  fMCFlag(AliAODMCParticle::kPhysicalPrim),
  fMinPtCharged(0.),
  fMinPtNeutral(0.),
  fMinECharged(0.),
  fMinENeutral(0.)
{
  fBaseClassName = "AliAODMCParticle";
  SetClassName("AliAODMCParticle");
}

AliMCParticleContainer::AliMCParticleContainer(const char *name):
  AliParticleContainer(name),
  fMCFlag(AliAODMCParticle::kPhysicalPrim),
  fMinPtCharged(0.),
  fMinPtNeutral(0.),
  fMinECharged(0.),
  fMinENeutral(0.)
{
  fBaseClassName = "AliAODMCParticle";
  SetClassName("AliAODMCParticle");
}

AliAODMCParticle* AliMCParticleContainer::GetMCParticleWithLabel(Int_t lab) const
{
  Int_t i = GetIndexFromLabel(lab);
  if (i >= 0) {
    return GetMCParticle(i);
  }
  else {
    return nullptr;
  }
}

AliAODMCParticle* AliMCParticleContainer::GetAcceptMCParticleWithLabel(Int_t lab)
{
  Int_t i = GetIndexFromLabel(lab);
  if (i >= 0) {
    return GetAcceptMCParticle(i);
  }
  else {
    return nullptr;
  }
}

AliAODMCParticle* AliMCParticleContainer::GetMCParticle(Int_t i) const
{
  if (i == -1) i = fCurrentID;
  if (i < 0 || i >= fClArray->GetEntriesFast()) return 0;
  AliAODMCParticle *vp = static_cast<AliAODMCParticle*>(fClArray->At(i));
  return vp;
}

AliAODMCParticle* AliMCParticleContainer::GetAcceptMCParticle(Int_t i) const
{
  //return pointer to particle if particle is accepted

  UInt_t rejectionReason = 0;
  if (i == -1) i = fCurrentID;
  if (AcceptMCParticle(i, rejectionReason)) {
      return GetMCParticle(i);
  }
  else {
    AliDebugStream(2) << "Particle " << i << " not accepted." << std::endl;
    return nullptr;
  }
}

AliAODMCParticle* AliMCParticleContainer::GetNextAcceptMCParticle()
{
  const Int_t n = GetNEntries();
  AliAODMCParticle *p = 0;
  do {
    fCurrentID++;
    if (fCurrentID >= n) break;
    p = GetAcceptMCParticle(fCurrentID);
  } while (!p);

  return p;
}

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

Bool_t AliMCParticleContainer::AcceptMCParticle(const AliAODMCParticle *vp, UInt_t &rejectionReason) const
{
  // Return true if vp is accepted.
  Bool_t r = ApplyMCParticleCuts(vp, rejectionReason);
  if (!r) return kFALSE;

  AliTLorentzVector mom;
  if(!GetMomentumFromParticle(mom, vp)) return false;

  return ApplyKinematicCuts(mom, rejectionReason);
}

Bool_t AliMCParticleContainer::AcceptMCParticle(Int_t i, UInt_t &rejectionReason) const
{
  // Return true if vp is accepted.

  Bool_t r = ApplyMCParticleCuts(GetMCParticle(i), rejectionReason);
  if (!r) return kFALSE;

  AliTLorentzVector mom;
  if (!GetMomentum(mom, i)) return kFALSE;

  return ApplyKinematicCuts(mom, rejectionReason);
}

Bool_t AliMCParticleContainer::ApplyMCParticleCuts(const AliAODMCParticle* vp, UInt_t &rejectionReason) const
{
  // Return true if i^th particle is accepted.
  // Cuts on the particle properties

  if ((vp->GetFlag() & fMCFlag) != fMCFlag) {
    rejectionReason |= kMCFlag;
    return kFALSE;
  }
  
  if(fMinPtCharged > 0.) {
    AliDebugStream(1) << "Applying pt-cut on charged particles: " << fMinPtCharged << std::endl;
    if((vp->Charge() != 0) && (TMath::Abs(vp->Pt()) < fMinPtCharged)) {
      rejectionReason |= kPtCut;
      return kFALSE;
    }
  }
  if(fMinPtNeutral > 0.) {
    AliDebugStream(1) << "Applying pt-cut on neutral particles: " << fMinPtNeutral << std::endl;
    if((vp->Charge() == 0) && (TMath::Abs(vp->Pt()) < fMinPtNeutral)) {
      rejectionReason |= kPtCut;
      return kFALSE;      
    }
  }
  if(fMinECharged > 0.) {
    AliDebugStream(1) << "Applying E-cut on charged particles: " << fMinECharged << std::endl;
    if((vp->Charge() != 0) && (vp->E() < fMinECharged)) {
      rejectionReason |= kPtCut;
      return kFALSE;
    }
  }
  if(fMinENeutral > 0.) {
    AliDebugStream(1) << "Applying E-cut on neutral particles: " << fMinENeutral << std::endl;
    if((vp->Charge() == 0) && (vp->E() < fMinENeutral)) {
      rejectionReason |= kPtCut;
      return kFALSE;      
    }
  }

  return ApplyParticleCuts(vp, rejectionReason);
}

const AliMCParticleIterableContainer AliMCParticleContainer::all() const {
  return AliMCParticleIterableContainer(this, false);
}

const AliMCParticleIterableContainer AliMCParticleContainer::accepted() const {
  return AliMCParticleIterableContainer(this, true);
}

const AliMCParticleIterableMomentumContainer AliMCParticleContainer::all_momentum() const {
  return AliMCParticleIterableMomentumContainer(this, false);
}

const AliMCParticleIterableMomentumContainer AliMCParticleContainer::accepted_momentum() const {
  return AliMCParticleIterableMomentumContainer(this, true);
}

const char* AliMCParticleContainer::GetTitle() const
{
  static TString trackString;
  trackString = TString::Format("%s_pT%04d", GetArrayName().Data(), static_cast<int>(GetMinPt()*1000.0));
  if(fMinPtCharged > 0.) {
    trackString += TString::Format("_pTCh%04d", static_cast<int>(fMinPtCharged*1000.0));
  }
  if(fMinPtNeutral > 0.) {
    trackString += TString::Format("_pTNe%04d", static_cast<int>(fMinPtNeutral*1000.0));
  }
  if(fMinECharged > 0.) {
    trackString += TString::Format("_ECh%04d", static_cast<int>(fMinECharged*1000.0));
  }
  if(fMinENeutral > 0.) {
    trackString += TString::Format("_ENe%04d", static_cast<int>(fMinENeutral*1000.0));
  }
  return trackString.Data();
}
