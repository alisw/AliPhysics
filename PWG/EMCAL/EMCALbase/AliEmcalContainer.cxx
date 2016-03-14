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

#include "AliEmcalContainer.h"

/// \cond CLASSIMP
ClassImp(AliEmcalContainer)
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
  fClassName(),
  fIsParticleLevel(kFALSE),
  fBitMap(0),
  fMinPt(0.15),
  fMaxPt(1000.),
  fMinE(0.),
  fMaxE(1000.),
  fMinEta(-0.9),
  fMaxEta(0.9),
  fMinPhi(-10),
  fMaxPhi(10),
  fMinMCLabel(-1),
  fMaxMCLabel(-1),
  fMassHypothesis(-1),
  fClArray(0),
  fCurrentID(0),
  fLabelMap(0),
  fRejectionReason(0),
  fLoadedClass(0)
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
  fClassName(),
  fIsParticleLevel(kFALSE),
  fBitMap(0),
  fMinPt(0.15),
  fMaxPt(1000.),
  fMinE(0.),
  fMaxE(1000.),
  fMinEta(-0.9),
  fMaxEta(0.9),
  fMinPhi(-10),
  fMaxPhi(10),
  fMinMCLabel(-1),
  fMaxMCLabel(-1),
  fMassHypothesis(-1),
  fClArray(0),
  fCurrentID(0),
  fLabelMap(0),
  fRejectionReason(0),
  fLoadedClass(0)
{
  fVertex[0] = 0;
  fVertex[1] = 0;
  fVertex[2] = 0;
}

/**
 * Connect the container to the array with content stored inside the virtual event.
 * The object name in the event must match the name given in the constructor
 * @param event Input event containing the array with content.
 */
void AliEmcalContainer::SetArray(AliVEvent *event) 
{

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
UShort_t AliEmcalContainer::GetRejectionReasonBitPosition() const
{ 
  UInt_t rs = fRejectionReason;
  UShort_t p = 0;
  while (rs >>= 1) { p++; }
  return p;
}

/**
 * Helper function to calculate the distance between two jets or a jet and a particle
 * @param part1 First particle in the check
 * @param part2 Second particle to compare to
 * @param dist Maximum distance under which partices are considered as "same" in \f$ p_{t} \f$, \f$ \eta \f$ and \$ \phi \f$
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
 * @param mom Momentum vector to select
 * @return True if the momentum vector is selected, false otherwise
 */
Bool_t AliEmcalContainer::ApplyKinematicCuts(const AliTLorentzVector& mom)
{
  if (mom.Pt() < fMinPt || mom.Pt() > fMaxPt) {
    fRejectionReason |= kPtCut;
    return kFALSE;
  }

  if (mom.E() < fMinE || mom.E() > fMaxE) {
    fRejectionReason |= kPtCut;
    return kFALSE;
  }

  Double_t eta = mom.Eta();
  Double_t phi = mom.Phi_0_2pi();

  if (fMinEta < fMaxEta && (eta < fMinEta || eta > fMaxEta)) {
    fRejectionReason |= kAcceptanceCut;
    return kFALSE;
  }

  if (fMinPhi < fMaxPhi && (phi < fMinPhi || phi > fMaxPhi)) {
    fRejectionReason |= kAcceptanceCut;
    return kFALSE;
  }

  return kTRUE;
}
