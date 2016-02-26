//
// Emcal Container Base class
//
// Author: M. Verweij


#include <TClonesArray.h>
#include "AliVEvent.h"
#include "AliLog.h"
#include "AliNamedArrayI.h"
#include "AliVParticle.h"
#include "AliTLorentzVector.h"

#include "AliEmcalContainer.h"

ClassImp(AliEmcalContainer)

//________________________________________________________________________
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
  // Default constructor.

  fVertex[0] = 0;
  fVertex[1] = 0;
  fVertex[2] = 0;
}

//________________________________________________________________________
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
  // Standard constructor.

  fVertex[0] = 0;
  fVertex[1] = 0;
  fVertex[2] = 0;
}

//________________________________________________________________________
void AliEmcalContainer::SetArray(AliVEvent *event) 
{
  // Get array from event.

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

//________________________________________________________________________
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

//________________________________________________________________________
UShort_t AliEmcalContainer::GetRejectionReasonBitPosition() const
{ 
  // Returns the highest bit in the rejection map.

  UInt_t rs = fRejectionReason;
  UShort_t p = 0;
  while (rs >>= 1) { p++; }
  return p;
}


//__________________________________________________________________________________________________
Bool_t AliEmcalContainer::SamePart(const AliVParticle* part1, const AliVParticle* part2, Double_t dist)
{
  // Helper function to calculate the distance between two jets or a jet and a particle
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

//________________________________________________________________________
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
