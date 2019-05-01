/************************************************************************************
 * Copyright (C) 2016, Copyright Holders of the ALICE Collaboration                 *
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
#include <TClonesArray.h>
#include "AliVEvent.h"
#include "AliLog.h"
#include "AliNamedArrayI.h"
#include "AliVParticle.h"
#include "AliTLorentzVector.h"

#include "AliEmcalContainerUtils.h"

#include "AliEmcalContainer.h"

ClassImp(AliEmcalContainer);

AliEmcalContainer::AliEmcalContainer():
  TObject(),
  fName(),
  fClArrayName(),
  fBaseClassName(),
  fIsParticleLevel(kFALSE),
  fBitMap(0),
  fMinPt(0.15),
  fMaxPt(1000.),
  fMaxE(1000.),
  fMinE(0.),
  fMinEta(-0.9),
  fMaxEta(0.9),
  fMinPhi(-10),
  fMaxPhi(10),
  fMinMCLabel(-1),
  fMaxMCLabel(-1),
  fMassHypothesis(-1),
  fIsEmbedding(kFALSE),
  fClArray(0),
  fCurrentID(0),
  fLabelMap(0),
  fLoadedClass(0),
  fClassName()
{
  fVertex[0] = 0;
  fVertex[1] = 0;
  fVertex[2] = 0;
}

AliEmcalContainer::AliEmcalContainer(const char *name):
  TObject(),
  fName(name),
  fClArrayName(name),
  fBaseClassName(),
  fIsParticleLevel(kFALSE),
  fBitMap(0),
  fMinPt(0.15),
  fMaxPt(1000.),
  fMaxE(1000.),
  fMinE(0.),
  fMinEta(-0.9),
  fMaxEta(0.9),
  fMinPhi(-10),
  fMaxPhi(10),
  fMinMCLabel(-1),
  fMaxMCLabel(-1),
  fMassHypothesis(-1),
  fIsEmbedding(kFALSE),
  fClArray(0),
  fCurrentID(0),
  fLabelMap(0),
  fLoadedClass(0),
  fClassName()
{
  fVertex[0] = 0;
  fVertex[1] = 0;
  fVertex[2] = 0;
}

TObject *AliEmcalContainer::operator[](int index) const {
  if(index >= 0 && index < GetNEntries()) return fClArray->At(index);
  return NULL;
}

void AliEmcalContainer::SetClassName(const char *clname)
{
  TClass cls(clname);
  if (cls.InheritsFrom(fBaseClassName)) {
    fClassName = clname;
  }
  else {
    AliError(Form("Unable to set class name %s for this container, it must inherits from %s!",clname,fBaseClassName.Data()));
  }
}

void AliEmcalContainer::GetVertexFromEvent(const AliVEvent * event)
{
  const AliVVertex *vertex = event->GetPrimaryVertex();
  if (vertex) vertex->GetXYZ(fVertex);
}

void AliEmcalContainer::SetArray(const AliVEvent *event)
{
  // Handling of default containers
  if(fClArrayName == "usedefault"){
    fClArrayName = GetDefaultArrayName(event);
  }

  // Get the right event (either the current event of the embedded event)
  event = AliEmcalContainerUtils::GetEvent(event, fIsEmbedding);

  if (!event) return;

  GetVertexFromEvent(event);

  if (!fClArrayName.IsNull() && !fClArray) {
    fClArray = dynamic_cast<TClonesArray*>(event->FindListObject(fClArrayName));
    if (!fClArray) {
      AliError(Form("%s: Could not retrieve array with name %s!", GetName(), fClArrayName.Data())); 
      return;
    }
  } else {
    return;
  }

  fLoadedClass = fClArray->GetClass();

  if (!fClassName.IsNull()) {
    if (!fLoadedClass->InheritsFrom(fClassName)) {
      AliError(Form("%s: Objects of type %s in %s are not inherited from %s!", 
		    GetName(), fLoadedClass->GetName(), fClArrayName.Data(), fClassName.Data()));
      fClArray = 0;
      fLoadedClass = 0;
    }
  }

  fLabelMap = dynamic_cast<AliNamedArrayI*>(event->FindListObject(fClArrayName + "_Map"));
}

void AliEmcalContainer::NextEvent(const AliVEvent * event)
{
  // Get the right event (either the current event of the embedded event)
  event = AliEmcalContainerUtils::GetEvent(event, fIsEmbedding);

  if (!event) return;

  GetVertexFromEvent(event);
}

Int_t AliEmcalContainer::GetNAcceptEntries() const{
  Int_t result = 0;
  for(int index = 0; index < GetNEntries(); index++){
    UInt_t rejectionReason = 0;
    if(AcceptObject(index, rejectionReason)) result++;
  }
  return result;
}

Int_t AliEmcalContainer::GetIndexFromLabel(Int_t lab) const
{ 
  if (fLabelMap) {
    if (lab < fLabelMap->GetSize()) {
      return fLabelMap->At(lab); 
    }
    else {
      AliDebug(3,Form("%s_AliEmcalContainer::GetIndexFromLabel - Label not found in the map, returning -1...",fClArrayName.Data()));
      return -1;
    }
  }
  else {
    AliDebug(3,Form("%s_AliEmcalContainer::GetIndexFromLabel - No index-label map found, returning label...",fClArrayName.Data()));
    return lab; 
  }
}

UShort_t AliEmcalContainer::GetRejectionReasonBitPosition(UInt_t rejectionReason)
{ 
  UInt_t rs = rejectionReason;
  UShort_t p = 0;
  while (rs >>= 1) { p++; }
  return p;
}

Bool_t AliEmcalContainer::SamePart(const AliVParticle* part1, const AliVParticle* part2, Double_t dist)
{
  if(!part1) return kFALSE;
  if(!part2) return kFALSE;
  Double_t dPhi = TMath::Abs(part1->Phi() - part2->Phi());
  Double_t dEta = TMath::Abs(part1->Eta() - part2->Eta());
  Double_t dpT  = TMath::Abs(part1->Pt() - part2->Pt());
  dPhi = TVector2::Phi_mpi_pi(dPhi);
  if (dPhi > dist) return kFALSE;
  if (dEta > dist) return kFALSE;
  if (dpT  > dist) return kFALSE;
  return kTRUE;
}

Bool_t AliEmcalContainer::ApplyKinematicCuts(const AliTLorentzVector& mom, UInt_t &rejectionReason) const
{
  if (mom.Pt() < fMinPt || mom.Pt() > fMaxPt) {
    rejectionReason |= kPtCut;
    return kFALSE;
  }

  if (mom.E() < fMinE || mom.E() > fMaxE) {
    rejectionReason |= kPtCut;
    return kFALSE;
  }

  Double_t eta = mom.Eta();
  Double_t phi = mom.Phi_0_2pi();

  if (fMinEta < fMaxEta && (eta < fMinEta || eta > fMaxEta)) {
    rejectionReason |= kAcceptanceCut;
    return kFALSE;
  }

  if (fMinPhi < fMaxPhi && (phi < fMinPhi || phi > fMaxPhi)) {
    rejectionReason |= kAcceptanceCut;
    return kFALSE;
  }

  return kTRUE;
}

const AliEmcalIterableContainer AliEmcalContainer::all() const {
  return AliEmcalIterableContainer(this, false);
}

const AliEmcalIterableContainer AliEmcalContainer::accepted() const {
  return AliEmcalIterableContainer(this, true);
}

const AliEmcalIterableMomentumContainer AliEmcalContainer::all_momentum() const {
  return AliEmcalIterableMomentumContainer(this, false);
}

const AliEmcalIterableMomentumContainer AliEmcalContainer::accepted_momentum() const {
  return AliEmcalIterableMomentumContainer(this, true);
}

Double_t AliEmcalContainer::RelativePhi(Double_t mphi, Double_t vphi)
{
  vphi = TVector2::Phi_0_2pi(vphi);
  mphi = TVector2::Phi_0_2pi(mphi);

  Double_t dphi = TVector2::Phi_mpi_pi(mphi - vphi);
  return dphi;
}
