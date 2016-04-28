/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
#include <TClonesArray.h>
#include "AliVEvent.h"
#include "AliLog.h"
#include "AliNamedArrayI.h"
#include "AliVParticle.h"
#include "AliTLorentzVector.h"

#include "AliAnalysisTaskEmcalEmbeddingHelper.h"

#include "AliEmcalContainer.h"

/// \cond CLASSIMP
ClassImp(AliEmcalContainer);
/// \endcond

/**
 * Default constructor. This constructor is only for ROOT I/O and
 * not to be used by users. The container will not connect to an
 * array in the input event.
 */
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

/**
 * Standard (named) constructor. The name provided must match the name of the
 * array inside the list objects in the input event the EMCAL container connects
 * to. The EMCAL container can get a different name, to be specified in the function
 * SetEvent.
 * @param name Name of the container in the input event.
 */
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

/**
 * Index operator, accessing object in the container at a given index.
 * Operates on all objects inside the container.
 * @param index Index of the object to access
 * @return Object at the given index (NULL if out of range)
 */
TObject *AliEmcalContainer::operator[](int index) const {
  if(index >= 0 && index < GetNEntries()) return fClArray->At(index);
  return NULL;
}

/**
 * Set the name of the class of the objets inside the underlying array.
 * @param[in] clname Name of the class of the object inside the underlying array.
 */
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

/**
 * Connect the container to the array with content stored inside the virtual event.
 * The object name in the event must match the name given in the constructor
 * @param event Input event containing the array with content.
 */
void AliEmcalContainer::SetArray(const AliVEvent *event)
{
  if (fIsEmbedding) {
    // this is an embedding container
    // will ignore the provided event and use the event
    // from the embedding helper class

    const AliAnalysisTaskEmcalEmbeddingHelper* embedding = AliAnalysisTaskEmcalEmbeddingHelper::GetInstance();
    if (!embedding) return;

    event = embedding->GetExternalEvent();
  }

  if (!event) return;

  const AliVVertex *vertex = event->GetPrimaryVertex();
  if (vertex) vertex->GetXYZ(fVertex);

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

/**
 * Count accepted entries in the container
 * @return Number of accepted events in the container
 */
Int_t AliEmcalContainer::GetNAcceptEntries() const{
  Int_t result = 0;
  for(int index = 0; index < GetNEntries(); index++){
    UInt_t rejectionReason = 0;
    if(AcceptObject(index, rejectionReason)) result++;
  }
  return result;
}

/**
 * Get the index in the container from a given label
 * @param lab Label to check
 * @return Index (-1 if not found)
 */
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

/**
 * Returns the highest bit in the rejection map as reason why the object
 * was rejected.
 * @return
 */
UShort_t AliEmcalContainer::GetRejectionReasonBitPosition(UInt_t rejectionReason)
{ 
  UInt_t rs = rejectionReason;
  UShort_t p = 0;
  while (rs >>= 1) { p++; }
  return p;
}

/**
 * Helper function to calculate the distance between two jets or a jet and a particle
 * @param part1 First particle in the check
 * @param part2 Second particle to compare to
 * @param dist Maximum distance under which partices are considered as "same" in \f$ p_{t} \f$, \f$ \eta \f$ and \f$ \phi \f$
 * @return True if the particles are considered as the same, false otherwise
 */
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

/**
 * Apply kinematical selection to the momentum vector provided. Selection is done in
 * - \f$ p_{t} \f$ (E)
 * - \f$ \eta \f$
 * - \f$ \phi \f$
 * @param[in] mom Momentum vector to select
 * @param[out] rejectionReason Bitmap for reason why object is rejected
 * @return True if the momentum vector is selected, false otherwise
 */
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

/**
 * Create an iterable container interface over all objects in the
 * EMCAL container.
 * @return iterable container over all objects in the EMCAL container
 */
const AliEmcalIterableContainer AliEmcalContainer::all() const {
  return AliEmcalIterableContainer(this, false);
}

/**
 * Create an iterable container interface over accepted objects in the
 * EMCAL container.
 * @return iterable container over accepted objects in the EMCAL container
 */
const AliEmcalIterableContainer AliEmcalContainer::accepted() const {
  return AliEmcalIterableContainer(this, true);
}

/**
 * Create an iterable container interface over all objects in the
 * EMCAL container.
 * @return iterable container over all objects in the EMCAL container
 */
const AliEmcalIterableMomentumContainer AliEmcalContainer::all_momentum() const {
  return AliEmcalIterableMomentumContainer(this, false);
}

/**
 * Create an iterable container interface over accepted objects in the
 * EMCAL container.
 * @return iterable container over accepted objects in the EMCAL container
 */
const AliEmcalIterableMomentumContainer AliEmcalContainer::accepted_momentum() const {
  return AliEmcalIterableMomentumContainer(this, true);
}

/**
 * Calculates the relative phi between two angle values and returns it in [-Pi, +Pi] range.
 * @param mphi First angle value
 * @param vphi Second angle value
 * @return Difference between mphi and vphi
 */
Double_t AliEmcalContainer::RelativePhi(Double_t mphi, Double_t vphi)
{
  vphi = TVector2::Phi_0_2pi(vphi);
  mphi = TVector2::Phi_0_2pi(mphi);

  Double_t dphi = TVector2::Phi_mpi_pi(mphi - vphi);
  return dphi;
}
