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
#include <algorithm>
#include <iostream>
#include <vector>
#include <TClonesArray.h>

#include "AliVEvent.h"
#include "AliLog.h"

#include "AliTLorentzVector.h"
#include "AliParticleContainer.h"

ClassImp(AliParticleContainer);

// Properly instantiate the object
AliEmcalContainerIndexMap <TClonesArray, AliVParticle> AliParticleContainer::fgEmcalContainerIndexMap;

AliParticleContainer::AliParticleContainer():
  AliEmcalContainer(),
  fMinDistanceTPCSectorEdge(-1),
  fChargeCut(kNoChargeCut),
  fGeneratorIndex(-1)
{
  fBaseClassName = "AliVParticle";
  SetClassName("AliVParticle");
}

AliParticleContainer::AliParticleContainer(const char *name) :
  AliEmcalContainer(name),
  fMinDistanceTPCSectorEdge(-1),
  fChargeCut(kNoChargeCut),
  fGeneratorIndex(-1)
{
  fBaseClassName = "AliVParticle";
  SetClassName("AliVParticle");
}

AliVParticle* AliParticleContainer::GetLeadingParticle(const char* opt) 
{
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

AliVParticle* AliParticleContainer::GetParticle(Int_t i) const 
{
  //
  if (i == -1) i = fCurrentID;
  if (i < 0 || i >= this->fClArray->GetEntriesFast()) return 0;
  AliVParticle *vp = static_cast<AliVParticle*>(fClArray->At(i));
  return vp;
}

AliVParticle* AliParticleContainer::GetAcceptParticle(Int_t i) const
{
  UInt_t rejectionReason = 0;
  if (i == -1) i = fCurrentID;
  if (AcceptParticle(i, rejectionReason)) {
      return GetParticle(i);
  }
  else {
    AliDebug(2,"Particle not accepted.");
    return 0;
  }
}

AliVParticle* AliParticleContainer::GetNextAcceptParticle()
{
  const Int_t n = GetNEntries();
  AliVParticle *p = 0;
  do {
    fCurrentID++;
    if (fCurrentID >= n) break;
    p = GetAcceptParticle(fCurrentID);
  } while (!p);

  return p;
}

AliVParticle* AliParticleContainer::GetNextParticle()
{
  const Int_t n = GetNEntries();
  AliVParticle *p = 0;
  do {
    fCurrentID++;
    if (fCurrentID >= n) break;
    p = GetParticle(fCurrentID);
  } while (!p);

  return p;
}

Bool_t AliParticleContainer::GetMomentumFromParticle(TLorentzVector &mom, const AliVParticle* part, Double_t mass) const
{
  if (part && part->Eta() < 1e6 && part->Eta() > -1e6) { // protection against FPE in sinh(eta)
    if (mass < 0) mass = part->M();
    mom.SetPtEtaPhiM(part->Pt(), part->Eta(), part->Phi(), mass);
    return kTRUE;
  }
  else {
    mom.SetPtEtaPhiM(0, 0, 0, 0);
    return kFALSE;
  }
}

Bool_t AliParticleContainer::GetMomentumFromParticle(TLorentzVector &mom, const AliVParticle* part) const
{
  return GetMomentumFromParticle(mom,part,fMassHypothesis);
}

Bool_t AliParticleContainer::GetMomentum(TLorentzVector &mom, Int_t i) const
{
  if (i == -1) i = fCurrentID;
  AliVParticle *vp = GetParticle(i);
  return GetMomentumFromParticle(mom, vp);
}

Bool_t AliParticleContainer::GetNextMomentum(TLorentzVector &mom)
{
  AliVParticle *vp = GetNextParticle();
  return GetMomentumFromParticle(mom, vp);
}

Bool_t AliParticleContainer::GetAcceptMomentum(TLorentzVector &mom, Int_t i) const
{
  if (i == -1) i = fCurrentID;
  AliVParticle *vp = GetAcceptParticle(i);
  return GetMomentumFromParticle(mom, vp);
}

Bool_t AliParticleContainer::GetNextAcceptMomentum(TLorentzVector &mom)
{
  AliVParticle *vp = GetNextAcceptParticle();
  return GetMomentumFromParticle(mom, vp);
}

Bool_t AliParticleContainer::AcceptParticle(const AliVParticle *vp, UInt_t &rejectionReason) const
{
  Bool_t r = ApplyParticleCuts(vp, rejectionReason);
  if (!r) return kFALSE;

  AliTLorentzVector mom;

  Int_t id = fClArray->IndexOf(vp);
  bool status(true);
  if (id >= 0) {
    status = GetMomentum(mom, id);
  }
  else {
    status = GetMomentumFromParticle(mom, vp);
  }
  if(!status) return false;

  return ApplyKinematicCuts(mom, rejectionReason);
}

Bool_t AliParticleContainer::AcceptParticle(Int_t i, UInt_t &rejectionReason) const
{
  Bool_t r = ApplyParticleCuts(GetParticle(i), rejectionReason);
  if (!r) return kFALSE;

  AliTLorentzVector mom;
  if(!GetMomentum(mom, i)) return false;

  return ApplyKinematicCuts(mom, rejectionReason);
}

Bool_t AliParticleContainer::ApplyParticleCuts(const AliVParticle* vp, UInt_t &rejectionReason) const
{
  if (!vp) {
    rejectionReason |= kNullObject;
    return kFALSE;
  }

  if (vp->TestBits(fBitMap) != (Int_t)fBitMap) {
    rejectionReason |= kBitMapCut;
    return kFALSE;
  }

  if (fMinMCLabel >= 0 && TMath::Abs(vp->GetLabel()) < fMinMCLabel) {
    rejectionReason |= kMCLabelCut;
    return kFALSE;
  }

  if (fMaxMCLabel >= 0 && TMath::Abs(vp->GetLabel()) > fMaxMCLabel) {
    rejectionReason |= kMCLabelCut;
    return kFALSE;
  }

  switch (fChargeCut) {
  case kCharged:
    if (vp->Charge() == 0) {
      rejectionReason |= kChargeCut;
      return kFALSE;
    }
    break;

  case kNeutral:
    if (vp->Charge() != 0) {
      rejectionReason |= kChargeCut;
      return kFALSE;
    }
    break;

  case kPositiveCharge:
    if (vp->Charge() <= 0) {
      rejectionReason |= kChargeCut;
      return kFALSE;
    }
    break;

  case kNegativeCharge:
    if (vp->Charge() >= 0) {
      rejectionReason |= kChargeCut;
      return kFALSE;
    }
    break;

  default:
    break;
  }

  if (fGeneratorIndex >= 0 && fGeneratorIndex != vp->GetGeneratorIndex()) {
    rejectionReason |= kMCGeneratorCut;
    return kFALSE;
  }

  return kTRUE;
}

Bool_t AliParticleContainer::ApplyKinematicCuts(const AliTLorentzVector& mom, UInt_t &rejectionReason) const
{
  if(fMinDistanceTPCSectorEdge>0.) {
    const Double_t pi = TMath::Pi();
    const Double_t kSector = pi/9;
    Double_t phiDist = TMath::Abs(mom.Phi() - TMath::FloorNint(mom.Phi()/kSector)*kSector);
    if(phiDist<fMinDistanceTPCSectorEdge) {
      rejectionReason |= kMinDistanceTPCSectorEdgeCut;
      return kFALSE;
    }
  }

  return AliEmcalContainer::ApplyKinematicCuts(mom, rejectionReason);
}

Int_t AliParticleContainer::GetNAcceptedParticles() const
{
  Int_t nPart = 0;
  for(int ipart = 0; ipart < this->GetNParticles(); ipart++){
    UInt_t rejectionReason = 0;
    if(this->AcceptParticle(ipart, rejectionReason)) nPart++;
  }
  return nPart;
}

const char* AliParticleContainer::GetTitle() const
{
  static TString trackString;
  trackString = TString::Format("%s_pT%04d", GetArrayName().Data(), static_cast<int>(GetMinPt()*1000.0));
  return trackString.Data();
}

void AliParticleContainer::SetArray(const AliVEvent * event)
{
  AliEmcalContainer::SetArray(event);

  // Register TClonesArray in index map
  fgEmcalContainerIndexMap.RegisterArray(GetArray());
}

const AliParticleIterableContainer AliParticleContainer::all() const {
  return AliParticleIterableContainer(this, false);
}

const AliParticleIterableContainer AliParticleContainer::accepted() const {
  return AliParticleIterableContainer(this, true);
}

const AliParticleIterableMomentumContainer AliParticleContainer::all_momentum() const {
  return AliParticleIterableMomentumContainer(this, false);
}

const AliParticleIterableMomentumContainer AliParticleContainer::accepted_momentum() const {
  return AliParticleIterableMomentumContainer(this, true);
}

/******************************************
 * Unit tests                             *
 ******************************************/

int TestParticleContainerIterator(const AliParticleContainer *const cont, int iteratorType, bool verbose){
  std::vector<AliVParticle *> reference, variation;
  AliVParticle *test = NULL;
  for(int iclust = 0; iclust < cont->GetNParticles(); iclust++){
    test = cont->GetParticle(iclust);
    if(!iteratorType){
      UInt_t rejectionReason = 0;
      if(!cont->AcceptParticle(test, rejectionReason)) continue;
    }
    reference.push_back(test);
  }

  if(!iteratorType){
    // test accept iterator
    for(auto part : cont->accepted()){
      variation.push_back(part);
    }
  } else {
    // test all iterator
    for(auto part : cont->all()){
      variation.push_back(part);
    }
  }

  if(verbose) {
    // Printing cluster addresses:
    if(reference.size() < 30){
      std::cout << "Paritcles in reference container: " << std::endl;
      std::cout << "===========================================" << std::endl;
      for(std::vector<AliVParticle *>::iterator refit = reference.begin(); refit != reference.end(); ++refit){
        std::cout << "Address: " << *refit << std::endl;
      }
    }
    if(variation.size() < 30){
      std::cout << "Paritcles in test container: " << std::endl;
      std::cout << "===========================================" << std::endl;
      for(std::vector<AliVParticle *>::iterator varit = variation.begin(); varit != variation.end(); ++varit){
        std::cout << "Address: " << *varit << std::endl;
      }
    }
  }

  int testresult = 0;
  // compare distributions: all particles in one vector needs to be found in the other vector and vice versa
  bool failure = false;
  // first test: cleck if all particles are found by the iterator
  for(std::vector<AliVParticle *>::iterator clit = reference.begin(); clit != reference.end(); ++clit){
    if(std::find(variation.begin(), variation.end(), *clit) == variation.end()) {
      if(verbose)
        std::cout << "Could not find particle with address " << *clit << " in test container" << std::endl;
      failure = true;
      break;
    }
  }
  if(!failure){
    // second test: check if there are no extra particles found
    for(std::vector<AliVParticle *>::iterator clit = variation.begin(); clit != variation.end(); ++clit){
      if(std::find(reference.begin(), reference.end(), *clit) == reference.end()) {
       if(verbose)
          std::cout << "Could not find particle with address " << *clit << " in reference container" << std::endl;
        failure = true;
        break;
      }
    }
    if(failure) testresult = 2;
  } else testresult = 1;
  if(verbose) {
    std::cout << "Unit test particle container, iterator type " << iteratorType << std::endl;
    std::cout << "Number of expected particles: " << reference.size() << ", number of found particles: " << variation.size() << std::endl;
    std::cout << "Test result: " << testresult << std::endl;
    std::cout << "-----------------------------------------------------------------------" << std::endl;
  }
  return testresult;
}
