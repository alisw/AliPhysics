#include "AliSigma0EventContainer.h"

#include <iostream>
#include "AliAODEvent.h"
#include "AliAODTracklets.h"
#include "TString.h"

ClassImp(AliSigma0EventContainer)

    //____________________________________________________________________________________________________
    AliSigma0EventContainer::AliSigma0EventContainer()
    : TObject(),
      fHistograms(nullptr),
      fHistogramsMC(nullptr),
      fHistogramsExtended(nullptr),
      fIsMC(false),
      fIsExtendedQA(false),
      fUseEmptyEvents(false),
      fInputEvent(nullptr),
      fMCEvent(nullptr),
      fMixingDepthPhoton(10),
      fMixingDepthProton(10),
      fMixingDepthLambda(10),
      fMixingDepthSigma(25),
      fZVertex(-1),
      fMultiplicity(-1),
      fMassSigma(1.192642),
      fSigmaMassCut(0.f),
      fSigmaSidebandLow(0.f),
      fSigmaSidebandUp(999.f),
      fProton(),
      fAntiProton(),
      fLambda(),
      fAntiLambda(),
      fSigma(),
      fAntiSigma(),
      fPhoton(),
      fProtonMixed(),
      fAntiProtonMixed(),
      fLambdaMixed(),
      fAntiLambdaMixed(),
      fSigmaMixed(),
      fAntiSigmaMixed(),
      fPhotonMixed(),
      fHistCuts(nullptr),
      fHistKProtonProton(nullptr),
      fHistKProtonProtonMixed(nullptr),
      fHistKProtonProtonMulti(),
      fHistKProtonProtonMixedMulti(),
      fHistKAntiProtonAntiProton(nullptr),
      fHistKAntiProtonAntiProtonMixed(nullptr),
      fHistKAntiProtonAntiProtonMulti(),
      fHistKAntiProtonAntiProtonMixedMulti(),
      fHistKProtonLambda(nullptr),
      fHistKProtonLambdaMixed(nullptr),
      fHistKProtonLambdaMulti(),
      fHistKProtonLambdaMixedMulti(),
      fHistKAntiProtonAntiLambda(nullptr),
      fHistKAntiProtonAntiLambdaMixed(nullptr),
      fHistKAntiProtonAntiLambdaMulti(),
      fHistKAntiProtonAntiLambdaMixedMulti(),
      fHistKLambdaLambda(nullptr),
      fHistKLambdaLambdaMixed(nullptr),
      fHistKLambdaLambdaMulti(),
      fHistKLambdaLambdaMixedMulti(),
      fHistKAntiLambdaAntiLambda(nullptr),
      fHistKAntiLambdaAntiLambdaMixed(nullptr),
      fHistKAntiLambdaAntiLambdaMulti(),
      fHistKAntiLambdaAntiLambdaMixedMulti(),
      fHistKProtonSigma(nullptr),
      fHistKProtonSigmaMixed(nullptr),
      fHistKProtonSigmaMulti(),
      fHistKProtonSigmaMixedMulti(),
      fHistKAntiProtonAntiSigma(nullptr),
      fHistKAntiProtonAntiSigmaMixed(nullptr),
      fHistKAntiProtonAntiSigmaMulti(),
      fHistKAntiProtonAntiSigmaMixedMulti(),
      fHistKProtonSigmaSidebandLower(nullptr),
      fHistKProtonSigmaMixedSidebandLower(nullptr),
      fHistKAntiProtonAntiSigmaSidebandLower(nullptr),
      fHistKAntiProtonAntiSigmaMixedSidebandLower(nullptr),
      fHistKProtonSigmaSidebandUpper(nullptr),
      fHistKProtonSigmaMixedSidebandUpper(nullptr),
      fHistKAntiProtonAntiSigmaSidebandUpper(nullptr),
      fHistKAntiProtonAntiSigmaMixedSidebandUpper(nullptr),
      fHistKProtonSigmaVsLambda(nullptr),
      fHistKAntiProtonAntiSigmaVsLambda(nullptr),
      fHistKProtonSigmaLambda(nullptr),
      fHistKProtonSigmaMixedLambda(nullptr),
      fHistKAntiProtonAntiSigmaLambda(nullptr),
      fHistKAntiProtonAntiSigmaMixedLambda(nullptr),
      fHistMCProtonProtonMomentum(nullptr),
      fHistMCAntiProtonAntiProtonMomentum(nullptr),
      fHistMCProtonLambdaMomentum(nullptr),
      fHistMCAntiProtonAntiLambdaMomentum(nullptr),
      fHistMCLambdaLambdaMomentum(nullptr),
      fHistMCAntiLambdaAntiLambdaMomentum(nullptr),
      fHistMCProtonSigmaMomentum(nullptr),
      fHistMCAntiProtonAntiSigmaMomentum(nullptr),
      fHistProtonProtonPhiStar(),
      fHistAntiProtonAntiProtonPhiStar(),
      fHistProtonLambdaProtonPhiStar(),
      fHistAntiProtonAntiLambdaProtonPhiStar(),
      fHistProtonLambdaPionPhiStar(),
      fHistAntiProtonAntiLambdaPionPhiStar(),
      fHistProtonProtonPhiStarMixed(),
      fHistAntiProtonAntiProtonPhiStarMixed(),
      fHistProtonLambdaProtonPhiStarMixed(),
      fHistAntiProtonAntiLambdaProtonPhiStarMixed(),
      fHistProtonLambdaPionPhiStarMixed(),
      fHistAntiProtonAntiLambdaPionPhiStarMixed(),
      fHistProtonLambdaNShared(nullptr),
      fHistProtonAntiLambdaNShared(nullptr),
      fHistAntiProtonLambdaNShared(nullptr),
      fHistAntiProtonAntiLambdaNShared(nullptr),
      fHistLambdaNShared(nullptr),
      fHistAntiLambdaNShared(nullptr),
      fHistSigmaMixedPt(nullptr),
      fHistSigmaMixedInvMass(nullptr),
      fHistSigmaMixedInvMassPt(nullptr),
      fHistSigmaMixedInvMassEta(nullptr),
      fHistAntiSigmaMixedPt(nullptr),
      fHistAntiSigmaMixedInvMass(nullptr),
      fHistAntiSigmaMixedInvMassPt(nullptr),
      fHistAntiSigmaMixedInvMassEta(nullptr),
      fHistProtonPt(nullptr),
      fHistAntiProtonPt(nullptr),
      fHistLambdaPt(nullptr),
      fHistLambdaInvMass(nullptr),
      fHistLambdaInvMassRec(nullptr),
      fHistAntiLambdaPt(nullptr),
      fHistAntiLambdaInvMass(nullptr),
      fHistAntiLambdaInvMassRec(nullptr),
      fHistSigmaPt(nullptr),
      fHistSigmaInvMass(nullptr),
      fHistSigmaInvMassRec(nullptr),
      fHistAntiSigmaPt(nullptr),
      fHistAntiSigmaInvMass(nullptr),
      fHistAntiSigmaInvMassRec(nullptr) {}

//____________________________________________________________________________________________________
AliSigma0EventContainer::AliSigma0EventContainer(
    const AliSigma0EventContainer &ref)
    : TObject(),
      fHistograms(nullptr),
      fHistogramsMC(nullptr),
      fHistogramsExtended(nullptr),
      fIsMC(false),
      fIsExtendedQA(false),
      fUseEmptyEvents(false),
      fInputEvent(nullptr),
      fMCEvent(nullptr),
      fMixingDepthPhoton(10),
      fMixingDepthProton(10),
      fMixingDepthLambda(10),
      fMixingDepthSigma(25),
      fZVertex(-1),
      fMultiplicity(-1),
      fMassSigma(1.192642),
      fSigmaMassCut(0.f),
      fSigmaSidebandLow(0.f),
      fSigmaSidebandUp(999.f),
      fProton(),
      fAntiProton(),
      fLambda(),
      fAntiLambda(),
      fSigma(),
      fAntiSigma(),
      fPhoton(),
      fProtonMixed(),
      fAntiProtonMixed(),
      fLambdaMixed(),
      fAntiLambdaMixed(),
      fSigmaMixed(),
      fAntiSigmaMixed(),
      fPhotonMixed(),
      fHistCuts(nullptr),
      fHistKProtonProton(nullptr),
      fHistKProtonProtonMixed(nullptr),
      fHistKProtonProtonMulti(),
      fHistKProtonProtonMixedMulti(),
      fHistKAntiProtonAntiProton(nullptr),
      fHistKAntiProtonAntiProtonMixed(nullptr),
      fHistKAntiProtonAntiProtonMulti(),
      fHistKAntiProtonAntiProtonMixedMulti(),
      fHistKProtonLambda(nullptr),
      fHistKProtonLambdaMixed(nullptr),
      fHistKProtonLambdaMulti(),
      fHistKProtonLambdaMixedMulti(),
      fHistKAntiProtonAntiLambda(nullptr),
      fHistKAntiProtonAntiLambdaMixed(nullptr),
      fHistKAntiProtonAntiLambdaMulti(),
      fHistKAntiProtonAntiLambdaMixedMulti(),
      fHistKLambdaLambda(nullptr),
      fHistKLambdaLambdaMixed(nullptr),
      fHistKLambdaLambdaMulti(),
      fHistKLambdaLambdaMixedMulti(),
      fHistKAntiLambdaAntiLambda(nullptr),
      fHistKAntiLambdaAntiLambdaMixed(nullptr),
      fHistKAntiLambdaAntiLambdaMulti(),
      fHistKAntiLambdaAntiLambdaMixedMulti(),
      fHistKProtonSigma(nullptr),
      fHistKProtonSigmaMixed(nullptr),
      fHistKProtonSigmaMulti(),
      fHistKProtonSigmaMixedMulti(),
      fHistKAntiProtonAntiSigma(nullptr),
      fHistKAntiProtonAntiSigmaMixed(nullptr),
      fHistKAntiProtonAntiSigmaMulti(),
      fHistKAntiProtonAntiSigmaMixedMulti(),
      fHistKProtonSigmaSidebandLower(nullptr),
      fHistKProtonSigmaMixedSidebandLower(nullptr),
      fHistKAntiProtonAntiSigmaSidebandLower(nullptr),
      fHistKAntiProtonAntiSigmaMixedSidebandLower(nullptr),
      fHistKProtonSigmaSidebandUpper(nullptr),
      fHistKProtonSigmaMixedSidebandUpper(nullptr),
      fHistKAntiProtonAntiSigmaSidebandUpper(nullptr),
      fHistKAntiProtonAntiSigmaMixedSidebandUpper(nullptr),
      fHistKProtonSigmaVsLambda(nullptr),
      fHistKAntiProtonAntiSigmaVsLambda(nullptr),
      fHistKProtonSigmaLambda(nullptr),
      fHistKProtonSigmaMixedLambda(nullptr),
      fHistKAntiProtonAntiSigmaLambda(nullptr),
      fHistKAntiProtonAntiSigmaMixedLambda(nullptr),
      fHistMCProtonProtonMomentum(nullptr),
      fHistMCAntiProtonAntiProtonMomentum(nullptr),
      fHistMCProtonLambdaMomentum(nullptr),
      fHistMCAntiProtonAntiLambdaMomentum(nullptr),
      fHistMCLambdaLambdaMomentum(nullptr),
      fHistMCAntiLambdaAntiLambdaMomentum(nullptr),
      fHistMCProtonSigmaMomentum(nullptr),
      fHistMCAntiProtonAntiSigmaMomentum(nullptr),
      fHistProtonProtonPhiStar(),
      fHistAntiProtonAntiProtonPhiStar(),
      fHistProtonLambdaProtonPhiStar(),
      fHistAntiProtonAntiLambdaProtonPhiStar(),
      fHistProtonLambdaPionPhiStar(),
      fHistAntiProtonAntiLambdaPionPhiStar(),
      fHistProtonProtonPhiStarMixed(),
      fHistAntiProtonAntiProtonPhiStarMixed(),
      fHistProtonLambdaProtonPhiStarMixed(),
      fHistAntiProtonAntiLambdaProtonPhiStarMixed(),
      fHistProtonLambdaPionPhiStarMixed(),
      fHistAntiProtonAntiLambdaPionPhiStarMixed(),
      fHistProtonLambdaNShared(nullptr),
      fHistProtonAntiLambdaNShared(nullptr),
      fHistAntiProtonLambdaNShared(nullptr),
      fHistAntiProtonAntiLambdaNShared(nullptr),
      fHistLambdaNShared(nullptr),
      fHistAntiLambdaNShared(nullptr),
      fHistSigmaMixedPt(nullptr),
      fHistSigmaMixedInvMass(nullptr),
      fHistSigmaMixedInvMassPt(nullptr),
      fHistSigmaMixedInvMassEta(nullptr),
      fHistAntiSigmaMixedPt(nullptr),
      fHistAntiSigmaMixedInvMass(nullptr),
      fHistAntiSigmaMixedInvMassPt(nullptr),
      fHistAntiSigmaMixedInvMassEta(nullptr),
      fHistProtonPt(nullptr),
      fHistAntiProtonPt(nullptr),
      fHistLambdaPt(nullptr),
      fHistLambdaInvMass(nullptr),
      fHistLambdaInvMassRec(nullptr),
      fHistAntiLambdaPt(nullptr),
      fHistAntiLambdaInvMass(nullptr),
      fHistAntiLambdaInvMassRec(nullptr),
      fHistSigmaPt(nullptr),
      fHistSigmaInvMass(nullptr),
      fHistSigmaInvMassRec(nullptr),
      fHistAntiSigmaPt(nullptr),
      fHistAntiSigmaInvMass(nullptr),
      fHistAntiSigmaInvMassRec(nullptr) {}

//____________________________________________________________________________________________________
AliSigma0EventContainer &AliSigma0EventContainer::operator=(
    const AliSigma0EventContainer &ref) {
  // Assignment operator
  if (this == &ref) return *this;
  return (*this);
}

//____________________________________________________________________________________________________
AliSigma0EventContainer::~AliSigma0EventContainer() {}

//____________________________________________________________________________________________________
void AliSigma0EventContainer::ProcessEvent(
    AliVEvent *inputEvent, AliMCEvent *mcEvent,
    vector<AliSigma0ParticleBase> &ParticleContainer,
    vector<AliSigma0ParticleBase> &AntiParticleContainer,
    vector<AliSigma0ParticleV0> &V0Container,
    vector<AliSigma0ParticleV0> &AntiV0Container,
    vector<AliSigma0ParticlePhotonMother> &SigmaContainer,
    vector<AliSigma0ParticlePhotonMother> &AntiSigmaContainer,
    vector<AliAODConversionPhoton> &PhotonContainer) {
  fInputEvent = inputEvent;
  fMCEvent = mcEvent;

  fZVertex = GetZvertexBins(inputEvent->GetPrimaryVertex()->GetZ());
  if (fZVertex < 0) return;

  fMultiplicity = GetMultiplicityBins(
      fInputEvent->GetMultiplicity()->GetNumberOfTracklets());
  if (fMultiplicity < 0) return;

  fProton = ParticleContainer;
  fAntiProton = AntiParticleContainer;
  fLambda = V0Container;
  fAntiLambda = AntiV0Container;
  fSigma = SigmaContainer;
  fAntiSigma = AntiSigmaContainer;
  fPhoton = PhotonContainer;

  CorrelationFunction();
  PhiStar();
  EventMixing();
  FillEventBuffer();
}

//____________________________________________________________________________________________________
void AliSigma0EventContainer::CorrelationFunction() {
  // PROTON - PROTON
  // Same Event
  for (auto ProtonCandidate1 = fProton.begin();
       ProtonCandidate1 != fProton.end(); ++ProtonCandidate1) {
    if (!ProtonCandidate1->GetIsUse()) continue;
    fHistProtonPt->Fill(ProtonCandidate1->GetPt());
    for (auto ProtonCandidate2 = ProtonCandidate1 + 1;
         ProtonCandidate2 != fProton.end(); ++ProtonCandidate2) {
      if (!ProtonCandidate2->GetIsUse()) continue;
      const float krel = ProtonCandidate1->ComputeRelK(*ProtonCandidate2);
      fHistKProtonProton->Fill(krel);
      if (fIsExtendedQA) fHistKProtonProtonMulti[fMultiplicity]->Fill(krel);
    }
  }

  // Mixed event
  for (const auto &ProtonContainer : fProtonMixed[fZVertex][fMultiplicity]) {
    for (const auto &ProtonLastEvent : ProtonContainer) {
      if (!ProtonLastEvent.GetIsUse()) continue;
      for (const auto &ProtonCandidate : fProton) {
        if (!ProtonCandidate.GetIsUse()) continue;
        const float krel = ProtonCandidate.ComputeRelK(ProtonLastEvent);
        fHistKProtonProtonMixed->Fill(krel);
        if (fIsExtendedQA)
          fHistKProtonProtonMixedMulti[fMultiplicity]->Fill(krel);
        if (fIsMC) {
          if (std::abs(ProtonCandidate.GetPDGcode()) == 2212 &&
              std::abs(ProtonLastEvent.GetPDGcode()) == 2212) {
            const float krelMC = ProtonCandidate.ComputeRelKMC(ProtonLastEvent);
            fHistMCProtonProtonMomentum->Fill(krelMC, krel);
          }
        }
      }
    }
  }

  // ANTI-PROTON - ANTI-PROTON
  // Same Event
  for (auto ProtonCandidate1 = fAntiProton.begin();
       ProtonCandidate1 != fAntiProton.end(); ++ProtonCandidate1) {
    if (!ProtonCandidate1->GetIsUse()) continue;
    fHistAntiProtonPt->Fill(ProtonCandidate1->GetPt());
    for (auto ProtonCandidate2 = ProtonCandidate1 + 1;
         ProtonCandidate2 != fAntiProton.end(); ++ProtonCandidate2) {
      if (!ProtonCandidate2->GetIsUse()) continue;
      const float krel = ProtonCandidate1->ComputeRelK(*ProtonCandidate2);
      fHistKAntiProtonAntiProton->Fill(krel);
      if (fIsExtendedQA)
        fHistKAntiProtonAntiProtonMulti[fMultiplicity]->Fill(krel);
    }
  }

  // Mixed event
  for (const auto &AntiProtonContainer :
       fAntiProtonMixed[fZVertex][fMultiplicity]) {
    for (const auto &ProtonLastEvent : AntiProtonContainer) {
      if (!ProtonLastEvent.GetIsUse()) continue;
      for (const auto &ProtonCandidate : fAntiProton) {
        if (!ProtonCandidate.GetIsUse()) continue;
        const float krel = ProtonCandidate.ComputeRelK(ProtonLastEvent);
        fHistKAntiProtonAntiProtonMixed->Fill(krel);
        if (fIsExtendedQA)
          fHistKAntiProtonAntiProtonMixedMulti[fMultiplicity]->Fill(krel);
        if (fIsMC) {
          if (std::abs(ProtonCandidate.GetPDGcode()) == 2212 &&
              std::abs(ProtonLastEvent.GetPDGcode()) == 2212) {
            const float krelMC = ProtonCandidate.ComputeRelKMC(ProtonLastEvent);
            fHistMCAntiProtonAntiProtonMomentum->Fill(krelMC, krel);
          }
        }
      }
    }
  }

  // PROTON - LAMBDA
  // Same Event
  for (const auto &ProtonCandidate : fProton) {
    if (!ProtonCandidate.GetIsUse()) continue;
    for (const auto &LambdaCandidate : fLambda) {
      if (!LambdaCandidate.GetIsUse()) continue;
      const float krel = ProtonCandidate.ComputeRelK(LambdaCandidate);
      fHistKProtonLambda->Fill(krel);
      if (fIsExtendedQA) fHistKProtonLambdaMulti[fMultiplicity]->Fill(krel);
    }
  }

  // Mixed event
  for (const auto &LambdaCandidate : fLambda) {
    if (!LambdaCandidate.GetIsUse()) continue;
    for (const auto &ProtonContainer : fProtonMixed[fZVertex][fMultiplicity]) {
      for (const auto &ProtonLastEvent : ProtonContainer) {
        if (!ProtonLastEvent.GetIsUse()) continue;
        const float krel = ProtonLastEvent.ComputeRelK(LambdaCandidate);
        fHistKProtonLambdaMixed->Fill(krel);
        if (fIsExtendedQA)
          fHistKProtonLambdaMixedMulti[fMultiplicity]->Fill(krel);
        if (fIsMC) {
          if (std::abs(LambdaCandidate.GetPDGcode()) == 3122 &&
              std::abs(ProtonLastEvent.GetPDGcode()) == 2212) {
            const float krelMC = ProtonLastEvent.ComputeRelKMC(LambdaCandidate);
            fHistMCProtonLambdaMomentum->Fill(krelMC, krel);
          }
        }
      }
    }
  }

  // ANTI-PROTON - ANTI-LAMBDA
  // Same Event
  for (const auto &ProtonCandidate : fAntiProton) {
    if (!ProtonCandidate.GetIsUse()) continue;
    for (const auto &LambdaCandidate : fAntiLambda) {
      if (!LambdaCandidate.GetIsUse()) continue;
      const float krel = ProtonCandidate.ComputeRelK(LambdaCandidate);
      fHistKAntiProtonAntiLambda->Fill(krel);
      if (fIsExtendedQA)
        fHistKAntiProtonAntiLambdaMulti[fMultiplicity]->Fill(krel);
    }
  }

  // Mixed event
  for (const auto &LambdaCandidate : fAntiLambda) {
    if (!LambdaCandidate.GetIsUse()) continue;
    for (const auto &ProtonContainer :
         fAntiProtonMixed[fZVertex][fMultiplicity]) {
      for (const auto &ProtonLastEvent : ProtonContainer) {
        if (!ProtonLastEvent.GetIsUse()) continue;
        const float krel = ProtonLastEvent.ComputeRelK(LambdaCandidate);
        fHistKAntiProtonAntiLambdaMixed->Fill(krel);
        if (fIsExtendedQA)
          fHistKAntiProtonAntiLambdaMixedMulti[fMultiplicity]->Fill(krel);
        if (fIsMC) {
          if (std::abs(LambdaCandidate.GetPDGcode()) == 3122 &&
              std::abs(ProtonLastEvent.GetPDGcode()) == 2212) {
            const float krelMC = ProtonLastEvent.ComputeRelKMC(LambdaCandidate);
            fHistMCAntiProtonAntiLambdaMomentum->Fill(krelMC, krel);
          }
        }
      }
    }
  }

  // LAMBDA - LAMBDA
  // same event
  for (auto LambdaCandidate1 = fLambda.begin();
       LambdaCandidate1 != fLambda.end(); ++LambdaCandidate1) {
    if (!LambdaCandidate1->GetIsUse()) continue;
    fHistLambdaPt->Fill(LambdaCandidate1->GetPt());
    fHistLambdaInvMass->Fill(LambdaCandidate1->GetMass());
    fHistLambdaInvMassRec->Fill(LambdaCandidate1->GetRecMass());
    for (auto LambdaCandidate2 = LambdaCandidate1 + 1;
         LambdaCandidate2 != fLambda.end(); ++LambdaCandidate2) {
      if (!LambdaCandidate2->GetIsUse()) continue;
      const float krel = LambdaCandidate1->ComputeRelK(*LambdaCandidate2);
      fHistKLambdaLambda->Fill(krel);
      if (fIsExtendedQA) fHistKLambdaLambdaMulti[fMultiplicity]->Fill(krel);
    }
  }

  // Mixed event
  for (const auto &LambdaContainer : fLambdaMixed[fZVertex][fMultiplicity]) {
    for (const auto &LambdaCandidate1 : LambdaContainer) {
      if (!LambdaCandidate1.GetIsUse()) continue;
      for (const auto &LambdaCandidate2 : fLambda) {
        if (!LambdaCandidate1.GetIsUse() || !LambdaCandidate2.GetIsUse())
          continue;
        const float krel = LambdaCandidate1.ComputeRelK(LambdaCandidate2);
        fHistKLambdaLambdaMixed->Fill(krel);
        if (fIsExtendedQA)
          fHistKLambdaLambdaMixedMulti[fMultiplicity]->Fill(krel);
        if (fIsMC) {
          if (std::abs(LambdaCandidate1.GetPDGcode()) == 3122 &&
              std::abs(LambdaCandidate2.GetPDGcode()) == 3122) {
            const float krelMC =
                LambdaCandidate1.ComputeRelKMC(LambdaCandidate2);
            fHistMCLambdaLambdaMomentum->Fill(krelMC, krel);
          }
        }
      }
    }
  }

  // ANTI-LAMBDA - ANTI-LAMBDA
  // same event
  for (auto LambdaCandidate1 = fAntiLambda.begin();
       LambdaCandidate1 != fAntiLambda.end(); ++LambdaCandidate1) {
    if (!LambdaCandidate1->GetIsUse()) continue;
    fHistAntiLambdaPt->Fill(LambdaCandidate1->GetPt());
    fHistAntiLambdaInvMass->Fill(LambdaCandidate1->GetMass());
    fHistAntiLambdaInvMassRec->Fill(LambdaCandidate1->GetRecMass());
    for (auto LambdaCandidate2 = LambdaCandidate1 + 1;
         LambdaCandidate2 != fAntiLambda.end(); ++LambdaCandidate2) {
      if (!LambdaCandidate2->GetIsUse()) continue;
      const float krel = LambdaCandidate1->ComputeRelK(*LambdaCandidate2);
      fHistKAntiLambdaAntiLambda->Fill(krel);
      if (fIsExtendedQA)
        fHistKAntiLambdaAntiLambdaMulti[fMultiplicity]->Fill(krel);
    }
  }

  // Mixed event
  for (const auto &AntiLambdaContainer :
       fAntiLambdaMixed[fZVertex][fMultiplicity]) {
    for (const auto &LambdaCandidate1 : AntiLambdaContainer) {
      if (!LambdaCandidate1.GetIsUse()) continue;
      for (const auto &LambdaCandidate2 : fAntiLambda) {
        if (!LambdaCandidate2.GetIsUse()) continue;
        const float krel = LambdaCandidate1.ComputeRelK(LambdaCandidate2);
        fHistKAntiLambdaAntiLambdaMixed->Fill(krel);
        if (fIsExtendedQA)
          fHistKAntiLambdaAntiLambdaMixedMulti[fMultiplicity]->Fill(krel);
        if (fIsMC) {
          if (std::abs(LambdaCandidate1.GetPDGcode()) == 3122 &&
              std::abs(LambdaCandidate2.GetPDGcode()) == 3122) {
            const float krelMC =
                LambdaCandidate1.ComputeRelKMC(LambdaCandidate2);
            fHistMCAntiLambdaAntiLambdaMomentum->Fill(krelMC, krel);
          }
        }
      }
    }
  }

  // PROTON - SIGMA
  // Same Event
  for (const auto &ProtonCandidate : fProton) {
    if (!ProtonCandidate.GetIsUse()) continue;
    for (const auto &SigmaCandidate : fSigma) {
      if (!SigmaCandidate.GetIsUse()) continue;
      const float invMass = SigmaCandidate.GetMass();
      if (invMass < fMassSigma + fSigmaMassCut &&
          invMass > fMassSigma - fSigmaMassCut) {
        fHistSigmaPt->Fill(SigmaCandidate.GetPt());
        fHistSigmaInvMass->Fill(invMass);
      }
      const float krel = ProtonCandidate.ComputeRelK(SigmaCandidate);
      if (invMass < fMassSigma + fSigmaMassCut &&
          invMass > fMassSigma - fSigmaMassCut) {
        fHistKProtonSigma->Fill(krel);
        if (fIsExtendedQA) fHistKProtonSigmaMulti[fMultiplicity]->Fill(krel);
        auto lambda = SigmaCandidate.GetV0();
        const float krelLambda = ProtonCandidate.ComputeRelK(lambda);
        fHistKProtonSigmaLambda->Fill(krelLambda);
        fHistKProtonSigmaVsLambda->Fill(krel, krelLambda);
      } else if (invMass > fMassSigma + fSigmaSidebandLow &&
                 invMass < fMassSigma + fSigmaSidebandUp) {
        fHistKProtonSigmaSidebandUpper->Fill(krel);
      } else if (invMass < fMassSigma - fSigmaSidebandLow &&
                 invMass > fMassSigma - fSigmaSidebandUp) {
        fHistKProtonSigmaSidebandLower->Fill(krel);
      }
    }
  }

  // Plotting only
  for (const auto &SigmaCandidate : fSigma) {
    if (!SigmaCandidate.GetIsUse()) continue;
    const float invMass = SigmaCandidate.GetMass();
    if (invMass < fMassSigma + fSigmaMassCut &&
        invMass > fMassSigma - fSigmaMassCut) {
      fHistSigmaPt->Fill(SigmaCandidate.GetPt());
      fHistSigmaInvMass->Fill(invMass);
      fHistSigmaInvMassRec->Fill(SigmaCandidate.GetRecMass());
    }
  }

  // Mixed event
  // Mixed sigmas with current protons
  for (const auto &SigmaContainer : fSigmaMixed[fZVertex]) {
    for (const auto &SigmaCandidate : SigmaContainer) {
      if (!SigmaCandidate.GetIsUse()) continue;
      const float invMass = SigmaCandidate.GetMass();
      auto lambda = SigmaCandidate.GetV0();
      for (const auto &ProtonCandidate : fProton) {
        if (!ProtonCandidate.GetIsUse()) continue;
        const float krel = ProtonCandidate.ComputeRelK(SigmaCandidate);
        const float krelLambda = ProtonCandidate.ComputeRelK(lambda);
        if (invMass < fMassSigma + fSigmaMassCut &&
            invMass > fMassSigma - fSigmaMassCut) {
          fHistKProtonSigmaMixed->Fill(krel);
          if (fIsExtendedQA)
            fHistKProtonSigmaMixedMulti[fMultiplicity]->Fill(krel);
          fHistKProtonSigmaMixedLambda->Fill(krelLambda);
          if (fIsMC) {
            const float krelMC = ProtonCandidate.ComputeRelKMC(SigmaCandidate);
            if (SigmaCandidate.IsTrueSigma(fMCEvent))
              fHistMCProtonSigmaMomentum->Fill(krelMC, krel);
          }
        } else if (invMass > fMassSigma + fSigmaSidebandLow &&
                   invMass < fMassSigma + fSigmaSidebandUp) {
          fHistKProtonSigmaMixedSidebandUpper->Fill(krel);
        } else if (invMass < fMassSigma - fSigmaSidebandLow &&
                   invMass > fMassSigma - fSigmaSidebandUp) {
          fHistKProtonSigmaMixedSidebandLower->Fill(krel);
        }
      }
    }
  }

  // Mixed protons with current sigmas
  for (const auto &ProtonContainer : fProtonMixed[fZVertex][fMultiplicity]) {
    for (const auto &ProtonCandidate : ProtonContainer) {
      if (!ProtonCandidate.GetIsUse()) continue;
      for (const auto &SigmaCandidate : fSigma) {
        if (!SigmaCandidate.GetIsUse()) continue;
        const float invMass = SigmaCandidate.GetMass();
        auto lambda = SigmaCandidate.GetV0();
        const float krel = ProtonCandidate.ComputeRelK(SigmaCandidate);
        const float krelLambda = ProtonCandidate.ComputeRelK(lambda);
        if (invMass < fMassSigma + fSigmaMassCut &&
            invMass > fMassSigma - fSigmaMassCut) {
          fHistKProtonSigmaMixed->Fill(krel);
          if (fIsExtendedQA)
            fHistKProtonSigmaMixedMulti[fMultiplicity]->Fill(krel);
          fHistKProtonSigmaMixedLambda->Fill(krelLambda);
          if (fIsMC) {
            const float krelMC = ProtonCandidate.ComputeRelKMC(SigmaCandidate);
            if (SigmaCandidate.IsTrueSigma(fMCEvent))
              fHistMCProtonSigmaMomentum->Fill(krelMC, krel);
          }
        } else if (invMass > fMassSigma + fSigmaSidebandLow &&
                   invMass < fMassSigma + fSigmaSidebandUp) {
          fHistKProtonSigmaMixedSidebandUpper->Fill(krel);
        } else if (invMass < fMassSigma - fSigmaSidebandLow &&
                   invMass > fMassSigma - fSigmaSidebandUp) {
          fHistKProtonSigmaMixedSidebandLower->Fill(krel);
        }
      }
    }
  }

  // ANTI-PROTON - ANTI-SIGMA
  // Same Event
  for (const auto &ProtonCandidate : fAntiProton) {
    if (!ProtonCandidate.GetIsUse()) continue;
    for (const auto &SigmaCandidate : fAntiSigma) {
      if (!SigmaCandidate.GetIsUse()) continue;
      const float invMass = SigmaCandidate.GetMass();
      const float krel = ProtonCandidate.ComputeRelK(SigmaCandidate);
      if (invMass < fMassSigma + fSigmaMassCut &&
          invMass > fMassSigma - fSigmaMassCut) {
        fHistKAntiProtonAntiSigma->Fill(krel);
        if (fIsExtendedQA)
          fHistKAntiProtonAntiSigmaMulti[fMultiplicity]->Fill(krel);
        auto lambda = SigmaCandidate.GetV0();
        const float krelLambda = ProtonCandidate.ComputeRelK(lambda);
        fHistKAntiProtonAntiSigmaLambda->Fill(krelLambda);
        fHistKAntiProtonAntiSigmaVsLambda->Fill(krel, krelLambda);
      } else if (invMass > fMassSigma + fSigmaSidebandLow &&
                 invMass < fMassSigma + fSigmaSidebandUp) {
        fHistKAntiProtonAntiSigmaSidebandUpper->Fill(krel);
      } else if (invMass < fMassSigma - fSigmaSidebandLow &&
                 invMass > fMassSigma - fSigmaSidebandUp) {
        fHistKAntiProtonAntiSigmaSidebandLower->Fill(krel);
      }
    }
  }

  // Plotting only
  for (const auto &SigmaCandidate : fAntiSigma) {
    if (!SigmaCandidate.GetIsUse()) continue;
    const float invMass = SigmaCandidate.GetMass();
    if (invMass < fMassSigma + fSigmaMassCut &&
        invMass > fMassSigma - fSigmaMassCut) {
      fHistAntiSigmaPt->Fill(SigmaCandidate.GetPt());
      fHistAntiSigmaInvMass->Fill(invMass);
      fHistAntiSigmaInvMassRec->Fill(SigmaCandidate.GetRecMass());
    }
  }

  // Mixed event
  // Mixed sigmas with current protons
  for (const auto &SigmaContainer : fAntiSigmaMixed[fZVertex]) {
    for (const auto &SigmaCandidate : SigmaContainer) {
      if (!SigmaCandidate.GetIsUse()) continue;
      auto lambda = SigmaCandidate.GetV0();
      const float invMass = SigmaCandidate.GetMass();
      for (const auto &ProtonCandidate : fAntiProton) {
        if (!ProtonCandidate.GetIsUse()) continue;
        const float krel = ProtonCandidate.ComputeRelK(SigmaCandidate);
        const float krelLambda = ProtonCandidate.ComputeRelK(lambda);
        if (invMass < fMassSigma + fSigmaMassCut &&
            invMass > fMassSigma - fSigmaMassCut) {
          fHistKAntiProtonAntiSigmaMixed->Fill(krel);
          if (fIsExtendedQA)
            fHistKAntiProtonAntiSigmaMixedMulti[fMultiplicity]->Fill(krel);
          fHistKAntiProtonAntiSigmaMixedLambda->Fill(krelLambda);
          if (fIsMC) {
            const float krelMC = ProtonCandidate.ComputeRelKMC(SigmaCandidate);
            if (SigmaCandidate.IsTrueSigma(fMCEvent))
              fHistMCAntiProtonAntiSigmaMomentum->Fill(krelMC, krel);
          }
        } else if (invMass > fMassSigma + fSigmaSidebandLow &&
                   invMass < fMassSigma + fSigmaSidebandUp) {
          fHistKAntiProtonAntiSigmaMixedSidebandUpper->Fill(krel);
        } else if (invMass < fMassSigma - fSigmaSidebandLow &&
                   invMass > fMassSigma - fSigmaSidebandUp) {
          fHistKAntiProtonAntiSigmaMixedSidebandLower->Fill(krel);
        }
      }
    }
  }

  // Mixed protons with current sigmas
  for (const auto &ProtonContainer :
       fAntiProtonMixed[fZVertex][fMultiplicity]) {
    for (const auto &ProtonCandidate : ProtonContainer) {
      if (!ProtonCandidate.GetIsUse()) continue;
      for (const auto &SigmaCandidate : fAntiSigma) {
        if (!SigmaCandidate.GetIsUse()) continue;
        const float invMass = SigmaCandidate.GetMass();
        auto lambda = SigmaCandidate.GetV0();
        const float krel = ProtonCandidate.ComputeRelK(SigmaCandidate);
        const float krelLambda = ProtonCandidate.ComputeRelK(lambda);
        if (invMass < fMassSigma + fSigmaMassCut &&
            invMass > fMassSigma - fSigmaMassCut) {
          fHistKAntiProtonAntiSigmaMixed->Fill(krel);
          if (fIsExtendedQA)
            fHistKAntiProtonAntiSigmaMixedMulti[fMultiplicity]->Fill(krel);
          fHistKAntiProtonAntiSigmaMixedLambda->Fill(krelLambda);
          if (fIsMC) {
            const float krelMC = ProtonCandidate.ComputeRelKMC(SigmaCandidate);
            if (SigmaCandidate.IsTrueSigma(fMCEvent))
              fHistMCAntiProtonAntiSigmaMomentum->Fill(krelMC, krel);
          }
        } else if (invMass > fMassSigma + fSigmaSidebandLow &&
                   invMass < fMassSigma + fSigmaSidebandUp) {
          fHistKAntiProtonAntiSigmaMixedSidebandUpper->Fill(krel);
        } else if (invMass < fMassSigma - fSigmaSidebandLow &&
                   invMass > fMassSigma - fSigmaSidebandUp) {
          fHistKAntiProtonAntiSigmaMixedSidebandLower->Fill(krel);
        }
      }
    }
  }
}

//____________________________________________________________________________________________________
void AliSigma0EventContainer::EventMixing() {
  // MIXED EVENT INV. MASS LAMBDA PHOTON
  // lambda
  // photons from this event with mixed v0s
  for (const auto &LambdaContainer : fLambdaMixed[fZVertex][fMultiplicity]) {
    for (const auto LambdaCandidate : LambdaContainer) {
      if (!LambdaCandidate.GetIsUse()) continue;
      for (auto PhotonCandidate : fPhoton) {
        AliSigma0ParticlePhotonMother *sigma =
            new AliSigma0ParticlePhotonMother(LambdaCandidate, PhotonCandidate,
                                              fInputEvent);
        fHistSigmaMixedPt->Fill(sigma->GetPt());
        fHistSigmaMixedInvMass->Fill(sigma->GetMass());
        fHistSigmaMixedInvMassPt->Fill(sigma->GetPt(), sigma->GetMass());
        fHistSigmaMixedInvMassEta->Fill(sigma->GetEta(), sigma->GetMass());
        delete sigma;
      }
    }
  }

  // lambdas from this event with mixed photons
  for (const auto PhotonContainer : fPhotonMixed[fZVertex][fMultiplicity]) {
    for (const auto PhotonCandidate : PhotonContainer) {
      for (const auto LambdaCandidate : fLambda) {
        if (!LambdaCandidate.GetIsUse()) continue;
        AliSigma0ParticlePhotonMother *sigma =
            new AliSigma0ParticlePhotonMother(LambdaCandidate, PhotonCandidate,
                                              fInputEvent);
        fHistSigmaMixedPt->Fill(sigma->GetPt());
        fHistSigmaMixedInvMass->Fill(sigma->GetMass());
        fHistSigmaMixedInvMassPt->Fill(sigma->GetPt(), sigma->GetMass());
        fHistSigmaMixedInvMassEta->Fill(sigma->GetEta(), sigma->GetMass());
        delete sigma;
      }
    }
  }

  // anti-lambda
  // photons from this event with mixed v0s
  for (const auto &AntiLambdaContainer :
       fAntiLambdaMixed[fZVertex][fMultiplicity]) {
    for (const auto LambdaCandidate : AntiLambdaContainer) {
      if (!LambdaCandidate.GetIsUse()) continue;
      for (const auto PhotonCandidate : fPhoton) {
        AliSigma0ParticlePhotonMother *antiSigma =
            new AliSigma0ParticlePhotonMother(LambdaCandidate, PhotonCandidate,
                                              fInputEvent);
        fHistAntiSigmaMixedPt->Fill(antiSigma->GetPt());
        fHistAntiSigmaMixedInvMass->Fill(antiSigma->GetMass());
        fHistAntiSigmaMixedInvMassPt->Fill(antiSigma->GetPt(),
                                           antiSigma->GetMass());
        fHistAntiSigmaMixedInvMassEta->Fill(antiSigma->GetEta(),
                                            antiSigma->GetMass());
        delete antiSigma;
      }
    }
  }

  // lambdas from this event with mixed photons
  for (const auto PhotonContainer : fPhotonMixed[fZVertex][fMultiplicity]) {
    for (const auto PhotonCandidate : PhotonContainer) {
      for (const auto LambdaCandidate : fAntiLambda) {
        if (!LambdaCandidate.GetIsUse()) continue;
        AliSigma0ParticlePhotonMother *antiSigma =
            new AliSigma0ParticlePhotonMother(LambdaCandidate, PhotonCandidate,
                                              fInputEvent);
        fHistAntiSigmaMixedPt->Fill(antiSigma->GetPt());
        fHistAntiSigmaMixedInvMass->Fill(antiSigma->GetMass());
        fHistAntiSigmaMixedInvMassPt->Fill(antiSigma->GetPt(),
                                           antiSigma->GetMass());
        fHistAntiSigmaMixedInvMassEta->Fill(antiSigma->GetEta(),
                                            antiSigma->GetMass());
        delete antiSigma;
      }
    }
  }
}

//____________________________________________________________________________________________________
void AliSigma0EventContainer::PhiStar() {
  // PROTON - PROTON
  // Same Event
  for (auto ProtonCandidate1 = fProton.begin();
       ProtonCandidate1 != fProton.end(); ++ProtonCandidate1) {
    if (!ProtonCandidate1->GetIsUse()) continue;
    for (auto ProtonCandidate2 = ProtonCandidate1 + 1;
         ProtonCandidate2 != fProton.end(); ++ProtonCandidate2) {
      if (!ProtonCandidate2->GetIsUse()) continue;
      for (int i = 0; i < 9; ++i) {
        const float deltaEta =
            ProtonCandidate1->GetEta() - ProtonCandidate2->GetEta();
        const float deltaPhiStar =
            ProtonCandidate1->GetPhiStar(i) - ProtonCandidate2->GetPhiStar(i);
        fHistProtonProtonPhiStar[i]->Fill(std::fabs(deltaEta),
                                          std::fabs(deltaPhiStar));
      }
    }
  }

  // Mixed event
  for (const auto &ProtonContainer : fProtonMixed[fZVertex][fMultiplicity]) {
    for (const auto &ProtonLastEvent : ProtonContainer) {
      if (!ProtonLastEvent.GetIsUse()) continue;
      for (const auto &ProtonCandidate : fProton) {
        if (!ProtonCandidate.GetIsUse()) continue;
        for (int i = 0; i < 9; ++i) {
          const float deltaEta =
              ProtonCandidate.GetEta() - ProtonLastEvent.GetEta();
          const float deltaPhiStar =
              ProtonCandidate.GetPhiStar(i) - ProtonLastEvent.GetPhiStar(i);
          fHistProtonProtonPhiStarMixed[i]->Fill(std::fabs(deltaEta),
                                                 std::fabs(deltaPhiStar));
        }
      }
    }
  }

  // ANTI-PROTON - ANTI-PROTON
  // Same Event
  for (auto ProtonCandidate1 = fAntiProton.begin();
       ProtonCandidate1 != fAntiProton.end(); ++ProtonCandidate1) {
    if (!ProtonCandidate1->GetIsUse()) continue;
    for (auto ProtonCandidate2 = ProtonCandidate1 + 1;
         ProtonCandidate2 != fAntiProton.end(); ++ProtonCandidate2) {
      if (!ProtonCandidate2->GetIsUse()) continue;
      for (int i = 0; i < 9; ++i) {
        const float deltaEta =
            ProtonCandidate1->GetEta() - ProtonCandidate2->GetEta();
        const float deltaPhiStar =
            ProtonCandidate1->GetPhiStar(i) - ProtonCandidate2->GetPhiStar(i);
        fHistAntiProtonAntiProtonPhiStar[i]->Fill(std::fabs(deltaEta),
                                                  std::fabs(deltaPhiStar));
      }
    }
  }

  // Mixed event
  for (const auto &AntiProtonContainer :
       fAntiProtonMixed[fZVertex][fMultiplicity]) {
    for (const auto &ProtonLastEvent : AntiProtonContainer) {
      if (!ProtonLastEvent.GetIsUse()) continue;
      for (const auto &ProtonCandidate : fAntiProton) {
        if (!ProtonCandidate.GetIsUse()) continue;
        for (int i = 0; i < 9; ++i) {
          const float deltaEta =
              ProtonCandidate.GetEta() - ProtonLastEvent.GetEta();
          const float deltaPhiStar =
              ProtonCandidate.GetPhiStar(i) - ProtonLastEvent.GetPhiStar(i);
          fHistAntiProtonAntiProtonPhiStarMixed[i]->Fill(
              std::fabs(deltaEta), std::fabs(deltaPhiStar));
        }
      }
    }
  }

  // PROTON - LAMBDA
  // Same Event
  for (const auto &ProtonCandidate : fProton) {
    if (!ProtonCandidate.GetIsUse()) continue;
    for (const auto &LambdaCandidate : fLambda) {
      if (!LambdaCandidate.GetIsUse()) continue;
      auto proton = LambdaCandidate.GetPosDaughter();
      auto pion = LambdaCandidate.GetNegDaughter();
      for (int i = 0; i < 9; ++i) {
        const float deltaEtaPion = ProtonCandidate.GetEta() - pion.GetEta();
        const float deltaPhiStarPion =
            ProtonCandidate.GetPhiStar(i) - pion.GetPhiStar(i);
        const float deltaEtaProton = ProtonCandidate.GetEta() - proton.GetEta();
        const float deltaPhiStarProton =
            ProtonCandidate.GetPhiStar(i) - proton.GetPhiStar(i);
        fHistProtonLambdaPionPhiStar[i]->Fill(std::fabs(deltaEtaPion),
                                              std::fabs(deltaPhiStarPion));
        fHistProtonLambdaProtonPhiStar[i]->Fill(std::fabs(deltaEtaProton),
                                                std::fabs(deltaPhiStarProton));
      }
    }
  }

  // Mixed event
  for (const auto &LambdaCandidate : fLambda) {
    if (!LambdaCandidate.GetIsUse()) continue;
    for (const auto &ProtonContainer : fProtonMixed[fZVertex][fMultiplicity]) {
      for (const auto &ProtonLastEvent : ProtonContainer) {
        if (!ProtonLastEvent.GetIsUse()) continue;
        auto proton = LambdaCandidate.GetPosDaughter();
        auto pion = LambdaCandidate.GetNegDaughter();
        for (int i = 0; i < 9; ++i) {
          const float deltaEtaPion = ProtonLastEvent.GetEta() - pion.GetEta();
          const float deltaPhiStarPion =
              ProtonLastEvent.GetPhiStar(i) - pion.GetPhiStar(i);
          const float deltaEtaProton =
              ProtonLastEvent.GetEta() - proton.GetEta();
          const float deltaPhiStarProton =
              ProtonLastEvent.GetPhiStar(i) - proton.GetPhiStar(i);
          fHistProtonLambdaPionPhiStarMixed[i]->Fill(
              std::fabs(deltaEtaPion), std::fabs(deltaPhiStarPion));
          fHistProtonLambdaProtonPhiStarMixed[i]->Fill(
              std::fabs(deltaEtaProton), std::fabs(deltaPhiStarProton));
        }
      }
    }
  }

  // ANTI-PROTON - ANTI-LAMBDA
  // Same Event
  for (const auto &ProtonCandidate : fAntiProton) {
    if (!ProtonCandidate.GetIsUse()) continue;
    for (const auto &LambdaCandidate : fAntiLambda) {
      if (!LambdaCandidate.GetIsUse()) continue;
      auto proton = LambdaCandidate.GetNegDaughter();
      auto pion = LambdaCandidate.GetPosDaughter();
      for (int i = 0; i < 9; ++i) {
        const float deltaEtaPion = ProtonCandidate.GetEta() - pion.GetEta();
        const float deltaPhiStarPion =
            ProtonCandidate.GetPhiStar(i) - pion.GetPhiStar(i);
        const float deltaEtaProton = ProtonCandidate.GetEta() - proton.GetEta();
        const float deltaPhiStarProton =
            ProtonCandidate.GetPhiStar(i) - proton.GetPhiStar(i);
        fHistAntiProtonAntiLambdaPionPhiStar[i]->Fill(
            std::fabs(deltaEtaPion), std::fabs(deltaPhiStarPion));
        fHistAntiProtonAntiLambdaProtonPhiStar[i]->Fill(
            std::fabs(deltaEtaProton), std::fabs(deltaPhiStarProton));
      }
    }
  }

  // Mixed event
  for (const auto &LambdaCandidate : fAntiLambda) {
    if (!LambdaCandidate.GetIsUse()) continue;
    for (const auto &ProtonContainer :
         fAntiProtonMixed[fZVertex][fMultiplicity]) {
      for (const auto &ProtonLastEvent : ProtonContainer) {
        if (!ProtonLastEvent.GetIsUse()) continue;
        auto proton = LambdaCandidate.GetNegDaughter();
        auto pion = LambdaCandidate.GetPosDaughter();
        for (int i = 0; i < 9; ++i) {
          const float deltaEtaPion = ProtonLastEvent.GetEta() - pion.GetEta();
          const float deltaPhiStarPion =
              ProtonLastEvent.GetPhiStar(i) - pion.GetPhiStar(i);
          const float deltaEtaProton =
              ProtonLastEvent.GetEta() - proton.GetEta();
          const float deltaPhiStarProton =
              ProtonLastEvent.GetPhiStar(i) - proton.GetPhiStar(i);
          fHistAntiProtonAntiLambdaPionPhiStarMixed[i]->Fill(
              std::fabs(deltaEtaPion), std::fabs(deltaPhiStarPion));
          fHistAntiProtonAntiLambdaProtonPhiStarMixed[i]->Fill(
              std::fabs(deltaEtaProton), std::fabs(deltaPhiStarProton));
        }
      }
    }
  }
}

//____________________________________________________________________________________________________
void AliSigma0EventContainer::FillEventBuffer() {
  // ++++++++++++++
  // Protons
  if (static_cast<int>(fProton.size()) > 0) {
    if (static_cast<int>(fProtonMixed[fZVertex][fMultiplicity].size()) <
        fMixingDepthProton) {
      fProtonMixed[fZVertex][fMultiplicity].push_back(fProton);
    } else {
      fProtonMixed[fZVertex][fMultiplicity].pop_front();
      fProtonMixed[fZVertex][fMultiplicity].push_back(fProton);
    }
  }

  if (static_cast<int>(fAntiProton.size()) > 0) {
    if (static_cast<int>(fAntiProtonMixed[fZVertex][fMultiplicity].size()) <
        fMixingDepthProton) {
      fAntiProtonMixed[fZVertex][fMultiplicity].push_back(fAntiProton);
    } else {
      fAntiProtonMixed[fZVertex][fMultiplicity].pop_front();
      fAntiProtonMixed[fZVertex][fMultiplicity].push_back(fAntiProton);
    }
  }

  // ++++++++++++++
  // Lambdas
  if (static_cast<int>(fLambda.size()) > 0) {
    if (static_cast<int>(fLambdaMixed[fZVertex][fMultiplicity].size()) <
        fMixingDepthLambda) {
      fLambdaMixed[fZVertex][fMultiplicity].push_back(fLambda);
    } else {
      fLambdaMixed[fZVertex][fMultiplicity].pop_front();
      fLambdaMixed[fZVertex][fMultiplicity].push_back(fLambda);
    }
  }

  if (static_cast<int>(fAntiLambda.size()) > 0) {
    if (static_cast<int>(fAntiLambdaMixed[fZVertex][fMultiplicity].size()) <
        fMixingDepthLambda) {
      fAntiLambdaMixed[fZVertex][fMultiplicity].push_back(fAntiLambda);
    } else {
      fAntiLambdaMixed[fZVertex][fMultiplicity].pop_front();
      fAntiLambdaMixed[fZVertex][fMultiplicity].push_back(fAntiLambda);
    }
  }

  // ++++++++++++++
  // Photons
  if (static_cast<int>(fPhoton.size()) > 0) {
    if (static_cast<int>(fPhotonMixed[fZVertex][fMultiplicity].size()) <
        fMixingDepthPhoton) {
      fPhotonMixed[fZVertex][fMultiplicity].push_back(fPhoton);
    } else {
      fPhotonMixed[fZVertex][fMultiplicity].pop_front();
      fPhotonMixed[fZVertex][fMultiplicity].push_back(fPhoton);
    }
  }

  // ++++++++++++++
  // Sigmas
  if (static_cast<int>(fSigma.size()) > 0) {
    if (static_cast<int>(fSigmaMixed[fZVertex].size()) < fMixingDepthSigma) {
      fSigmaMixed[fZVertex].push_back(fSigma);
    } else {
      fSigmaMixed[fZVertex].pop_front();
      fSigmaMixed[fZVertex].push_back(fSigma);
    }
  }

  if (static_cast<int>(fAntiSigma.size()) > 0) {
    if (static_cast<int>(fAntiSigmaMixed[fZVertex].size()) <
        fMixingDepthSigma) {
      fAntiSigmaMixed[fZVertex].push_back(fAntiSigma);
    } else {
      fAntiSigmaMixed[fZVertex].pop_front();
      fAntiSigmaMixed[fZVertex].push_back(fAntiSigma);
    }
  }
}

//____________________________________________________________________________________________________
void AliSigma0EventContainer::TrackCleaner(
    vector<AliSigma0ParticleBase> &ParticleContainer,
    vector<AliSigma0ParticleBase> &AntiParticleContainer,
    vector<AliSigma0ParticleV0> &V0Container,
    vector<AliSigma0ParticleV0> &AntiV0Container) {
  // +++++++++++++++++++++++++++++++++++++++++++
  // clean proton - lambda
  int nShared = 0;
  for (auto &LambdaCandidate : V0Container) {
    if (!LambdaCandidate.GetIsUse()) continue;
    for (const auto ProtonCandidate : ParticleContainer) {
      if ((LambdaCandidate.GetTrackLabelNeg() ==
              ProtonCandidate.GetTrackLabel()) ||
          (LambdaCandidate.GetTrackLabelPos() ==
              ProtonCandidate.GetTrackLabel())) {
        LambdaCandidate.SetUse(false);
        ++nShared;
      }
    }
  }
  fHistProtonLambdaNShared->Fill(nShared);

  // +++++++++++++++++++++++++++++++++++++++++++
  // clean proton - anti-lambda
  nShared = 0;
  for (auto &LambdaCandidate : AntiV0Container) {
    if (!LambdaCandidate.GetIsUse()) continue;
    for (const auto ProtonCandidate : ParticleContainer) {
      if ((LambdaCandidate.GetTrackLabelNeg() ==
              ProtonCandidate.GetTrackLabel()) ||
          (LambdaCandidate.GetTrackLabelPos() ==
              ProtonCandidate.GetTrackLabel())) {
        LambdaCandidate.SetUse(false);
        ++nShared;
      }
    }
  }
  fHistProtonAntiLambdaNShared->Fill(nShared);

  // +++++++++++++++++++++++++++++++++++++++++++
  // clean anti-proton - lambda
  nShared = 0;
  for (auto &LambdaCandidate : V0Container) {
    if (!LambdaCandidate.GetIsUse()) continue;
    for (const auto ProtonCandidate : AntiParticleContainer) {
      if ((LambdaCandidate.GetTrackLabelNeg() ==
              ProtonCandidate.GetTrackLabel()) ||
          (LambdaCandidate.GetTrackLabelPos() ==
              ProtonCandidate.GetTrackLabel())) {
        LambdaCandidate.SetUse(false);
        ++nShared;
      }
    }
  }
  fHistAntiProtonLambdaNShared->Fill(nShared);

  // +++++++++++++++++++++++++++++++++++++++++++
  // clean anti-proton - anti-lambda
  nShared = 0;
  for (auto &LambdaCandidate : AntiV0Container) {
    if (!LambdaCandidate.GetIsUse()) continue;
    for (const auto ProtonCandidate : AntiParticleContainer) {
      if ((LambdaCandidate.GetTrackLabelNeg() ==
              ProtonCandidate.GetTrackLabel()) ||
          (LambdaCandidate.GetTrackLabelPos() ==
              ProtonCandidate.GetTrackLabel())) {
        LambdaCandidate.SetUse(false);
        ++nShared;
      }
    }
  }
  fHistAntiProtonAntiLambdaNShared->Fill(nShared);

  // +++++++++++++++++++++++++++++++++++++++++++
  // clean lambda - lambda
  nShared = 0;
  for (auto iLambdaCandidate = V0Container.begin();
       iLambdaCandidate != V0Container.end(); ++iLambdaCandidate) {
    if (!iLambdaCandidate->GetIsUse()) continue;
    for (auto jLambdaCandidate = iLambdaCandidate + 1;
         jLambdaCandidate != V0Container.end(); ++jLambdaCandidate) {
      if (!jLambdaCandidate->GetIsUse()) continue;
      if ((iLambdaCandidate->GetTrackLabelNeg() ==
              jLambdaCandidate->GetTrackLabelNeg()) ||
          (iLambdaCandidate->GetTrackLabelPos() ==
              jLambdaCandidate->GetTrackLabelPos()) ||
          (iLambdaCandidate->GetTrackLabelPos() ==
              jLambdaCandidate->GetTrackLabelNeg()) ||
          (iLambdaCandidate->GetTrackLabelNeg() ==
              jLambdaCandidate->GetTrackLabelPos())) {
        ++nShared;
        if (iLambdaCandidate->GetCosineAlpha() >=
            jLambdaCandidate->GetCosineAlpha()) {
          // keep i, remove j
          jLambdaCandidate->SetUse(false);
        } else {
          // remove i, keep j
          iLambdaCandidate->SetUse(false);
          break;
        }
      }
    }
  }
  fHistLambdaNShared->Fill(nShared);

  //   +++++++++++++++++++++++++++++++++++++++++++
  //   clean anti-lambda - anti-lambda
  nShared = 0;
  for (auto iLambdaCandidate = AntiV0Container.begin();
       iLambdaCandidate != AntiV0Container.end(); ++iLambdaCandidate) {
    if (!iLambdaCandidate->GetIsUse()) continue;
    for (auto jLambdaCandidate = iLambdaCandidate + 1;
         jLambdaCandidate != AntiV0Container.end(); ++jLambdaCandidate) {
      if (!jLambdaCandidate->GetIsUse()) continue;
      if ((iLambdaCandidate->GetTrackLabelNeg() ==
              jLambdaCandidate->GetTrackLabelNeg()) ||
          (iLambdaCandidate->GetTrackLabelPos() ==
              jLambdaCandidate->GetTrackLabelPos()) ||
          (iLambdaCandidate->GetTrackLabelPos() ==
              jLambdaCandidate->GetTrackLabelNeg()) ||
          (iLambdaCandidate->GetTrackLabelNeg() ==
              jLambdaCandidate->GetTrackLabelPos())) {
        ++nShared;
        if (iLambdaCandidate->GetCosineAlpha() >=
            jLambdaCandidate->GetCosineAlpha()) {
          // keep i, remove j
          jLambdaCandidate->SetUse(false);
        } else {
          // remove i, keep j
          iLambdaCandidate->SetUse(false);
          break;
        }
      }
    }
  }
  fHistAntiLambdaNShared->Fill(nShared);
}

//____________________________________________________________________________________________________
void AliSigma0EventContainer::SetZvertexBins(float max, float min, float step) {
  //  float nSteps = max-min;
  //  fNZVertexBins = int(nSteps/step);
  //  if (std::fabs(int(nSteps/step) - (nSteps/step)) > 0.01) {
  //    AliFatal("Non integer number of bins");
  //  }

  //  float currentMax = max;
  //  for (int i=0; i<nSteps; ++i) {
  //    if(currentMax < min) break;
  //    fZVertexBins.push_back(currentMax);
  //    currentMax -= step;
  //  }

  //  fEvents.resize(fNZVertexBins);
  //  fSigmaMixed.resize(fNZVertexBins);
  //  fAntiSigmaMixed.resize(fNZVertexBins);
}

//____________________________________________________________________________________________________
int AliSigma0EventContainer::GetZvertexBins(float zVertex) const {
  //  // check in the first bin - necessary to consistently use >= and <
  //  if(fZVertexBins[0] > zVertex && fZVertexBins[1] < zVertex) return 0;

  //  // check the rest of the vector
  //  for( int i=1; i<static_cast<int>(fZVertexBins.size()-1); ++i ) {
  //    if(fZVertexBins[i] >= zVertex && fZVertexBins[i+1] < zVertex) return i;
  //  }

  //  // if not within the specified range return -1
  //  return -1;

  if (zVertex >= 8 && zVertex < 10)
    return 0;
  else if (zVertex >= 6 && zVertex < 8)
    return 1;
  else if (zVertex >= 4 && zVertex < 6)
    return 2;
  else if (zVertex >= 2 && zVertex < 4)
    return 3;
  else if (zVertex >= 0 && zVertex < 2)
    return 4;
  else if (zVertex >= -2 && zVertex < 0)
    return 5;
  else if (zVertex >= -4 && zVertex < -2)
    return 6;
  else if (zVertex >= -6 && zVertex < -4)
    return 7;
  else if (zVertex >= -8 && zVertex < -6)
    return 8;
  else if (zVertex > -10 && zVertex < -8)
    return 9;
  else
    return -1;  // should not happen
}

//____________________________________________________________________________________________________
int AliSigma0EventContainer::GetMultiplicityBins(int multiplicity) const {
  if (multiplicity > 0 && multiplicity < 5)
    return 0;
  else if (multiplicity >= 5 && multiplicity < 9)
    return 1;
  else if (multiplicity >= 9 && multiplicity < 13)
    return 2;
  else if (multiplicity >= 13 && multiplicity < 17)
    return 3;
  else if (multiplicity >= 17 && multiplicity < 21)
    return 4;
  else if (multiplicity >= 21 && multiplicity < 25)
    return 5;
  else if (multiplicity >= 25 && multiplicity < 29)
    return 6;
  else if (multiplicity >= 29 && multiplicity < 33)
    return 7;
  else if (multiplicity >= 33 && multiplicity < 37)
    return 8;
  else if (multiplicity >= 37 && multiplicity < 41)
    return 9;
  else if (multiplicity >= 41 && multiplicity < 61)
    return 10;
  else if (multiplicity >= 61 && multiplicity < 80)
    return 11;
  else if (multiplicity > 80)
    return 12;
  else
    return -1;  // should not happen
}

//____________________________________________________________________________________________________
void AliSigma0EventContainer::InitCutHistograms() {
  std::cout << "============================\n"
            << " EVENT CONTAINER CONFIGURATION \n"
            << " Sigma0 mass     " << fMassSigma << "\n"
            << " Sigma0 select   " << fSigmaMassCut << "\n"
            << " Sigma0 sideb    " << fSigmaSidebandLow << " "
            << fSigmaSidebandUp << "\n"
            << " Mixing depth p  " << fMixingDepthProton << "\n"
            << " Mixing depth L  " << fMixingDepthLambda << "\n"
            << " Mixing depth Ph " << fMixingDepthProton << "\n"
            << " Mixing depth Si " << fMixingDepthSigma << "\n"
            << "============================\n";

  TH1::AddDirectory(kFALSE);

  if (fHistograms != nullptr) {
    delete fHistograms;
    fHistograms = nullptr;
  }
  if (fHistograms == nullptr) {
    fHistograms = new TList();
    fHistograms->SetOwner(kTRUE);
    fHistograms->SetName("EventCollection");
  }

  fHistCuts = new TH1F("fHistCuts", ";;Cut values", 10, 0, 10);
  fHistCuts->GetXaxis()->SetBinLabel(1, "Mixing depth p");
  fHistCuts->GetXaxis()->SetBinLabel(2, "Mixing depth #Lambda");
  fHistCuts->GetXaxis()->SetBinLabel(3, "Mixing depth #gamma");
  fHistCuts->GetXaxis()->SetBinLabel(4, "Mixing depth #Sigma^{0}");
  fHistCuts->GetXaxis()->SetBinLabel(5, "#Sigma^{0} selection");
  fHistCuts->GetXaxis()->SetBinLabel(6, "#Sigma^{0} sideband up");
  fHistCuts->GetXaxis()->SetBinLabel(7, "#Sigma^{0} sideband low");
  fHistograms->Add(fHistCuts);

  fHistCuts->Fill(0.f, fMixingDepthProton);
  fHistCuts->Fill(1.f, fMixingDepthLambda);
  fHistCuts->Fill(2.f, fMixingDepthPhoton);
  fHistCuts->Fill(3.f, fMixingDepthSigma);
  fHistCuts->Fill(4.f, fSigmaMassCut);
  fHistCuts->Fill(5.f, fSigmaSidebandUp);
  fHistCuts->Fill(6.f, fSigmaSidebandLow);

  fHistKProtonProton =
      new TH1F("fHistKProtonProton", "; k* [GeV/#it{c}]; Entries", 3000, 0, 3);
  fHistKProtonProtonMixed = new TH1F("fHistKProtonProtonMixed",
                                     "; k* [GeV/#it{c}]; Entries", 3000, 0, 3);
  fHistKAntiProtonAntiProton = new TH1F(
      "fHistKAntiProtonAntiProton", "; k* [GeV/#it{c}]; Entries", 3000, 0, 3);
  fHistKAntiProtonAntiProtonMixed =
      new TH1F("fHistKAntiProtonAntiProtonMixed", "; k* [GeV/#it{c}]; Entries",
               3000, 0, 3);
  fHistograms->Add(fHistKProtonProton);
  fHistograms->Add(fHistKProtonProtonMixed);
  fHistograms->Add(fHistKAntiProtonAntiProton);
  fHistograms->Add(fHistKAntiProtonAntiProtonMixed);
  fHistKProtonLambda =
      new TH1F("fHistKProtonLambda", "; k* [GeV/#it{c}]; Entries", 3000, 0, 3);
  fHistKProtonLambdaMixed = new TH1F("fHistKProtonLambdaMixed",
                                     "; k* [GeV/#it{c}]; Entries", 3000, 0, 3);
  fHistKAntiProtonAntiLambda = new TH1F(
      "fHistKAntiProtonAntiLambda", "; k* [GeV/#it{c}]; Entries", 3000, 0, 3);
  fHistKAntiProtonAntiLambdaMixed =
      new TH1F("fHistKAntiProtonAntiLambdaMixed", "; k* [GeV/#it{c}]; Entries",
               3000, 0, 3);
  fHistograms->Add(fHistKProtonLambda);
  fHistograms->Add(fHistKProtonLambdaMixed);
  fHistograms->Add(fHistKAntiProtonAntiLambda);
  fHistograms->Add(fHistKAntiProtonAntiLambdaMixed);
  fHistKLambdaLambda =
      new TH1F("fHistKLambdaLambda", "; k* [GeV/#it{c}]; Entries", 3000, 0, 3);
  fHistKLambdaLambdaMixed = new TH1F("fHistKLambdaLambdaMixed",
                                     "; k* [GeV/#it{c}]; Entries", 3000, 0, 3);
  fHistKAntiLambdaAntiLambda = new TH1F(
      "fHistKAntiLambdaAntiLambda", "; k* [GeV/#it{c}]; Entries", 3000, 0, 3);
  fHistKAntiLambdaAntiLambdaMixed =
      new TH1F("fHistKAntiLambdaAntiLambdaMixed", "; k* [GeV/#it{c}]; Entries",
               3000, 0, 3);
  fHistograms->Add(fHistKLambdaLambda);
  fHistograms->Add(fHistKLambdaLambdaMixed);
  fHistograms->Add(fHistKAntiLambdaAntiLambda);
  fHistograms->Add(fHistKAntiLambdaAntiLambdaMixed);
  fHistKProtonSigma =
      new TH1F("fHistKProtonSigma", "; k* [GeV/#it{c}]; Entries", 3000, 0, 3);
  fHistKProtonSigmaMixed = new TH1F("fHistKProtonSigmaMixed",
                                    "; k* [GeV/#it{c}]; Entries", 3000, 0, 3);
  fHistKAntiProtonAntiSigma = new TH1F(
      "fHistKAntiProtonAntiSigma", "; k* [GeV/#it{c}]; Entries", 3000, 0, 3);
  fHistKAntiProtonAntiSigmaMixed =
      new TH1F("fHistKAntiProtonAntiSigmaMixed", "; k* [GeV/#it{c}]; Entries",
               3000, 0, 3);
  fHistograms->Add(fHistKProtonSigma);
  fHistograms->Add(fHistKProtonSigmaMixed);
  fHistograms->Add(fHistKAntiProtonAntiSigma);
  fHistograms->Add(fHistKAntiProtonAntiSigmaMixed);

  fHistKProtonSigmaSidebandLower =
      new TH1F("fHistKProtonSigmaSidebandLower", "; k* [GeV/#it{c}]; Entries",
               3000, 0, 3);
  fHistKProtonSigmaMixedSidebandLower =
      new TH1F("fHistKProtonSigmaMixedSidebandLower",
               "; k* [GeV/#it{c}]; Entries", 3000, 0, 3);
  fHistKAntiProtonAntiSigmaSidebandLower =
      new TH1F("fHistKAntiProtonAntiSigmaSidebandLower",
               "; k* [GeV/#it{c}]; Entries", 3000, 0, 3);
  fHistKAntiProtonAntiSigmaMixedSidebandLower =
      new TH1F("fHistKAntiProtonAntiSigmaMixedSidebandLower",
               "; k* [GeV/#it{c}]; Entries", 3000, 0, 3);
  fHistograms->Add(fHistKProtonSigmaSidebandLower);
  fHistograms->Add(fHistKProtonSigmaMixedSidebandLower);
  fHistograms->Add(fHistKAntiProtonAntiSigmaSidebandLower);
  fHistograms->Add(fHistKAntiProtonAntiSigmaMixedSidebandLower);

  fHistKProtonSigmaSidebandUpper =
      new TH1F("fHistKProtonSigmaSidebandUpper", "; k* [GeV/#it{c}]; Entries",
               3000, 0, 3);
  fHistKProtonSigmaMixedSidebandUpper =
      new TH1F("fHistKProtonSigmaMixedSidebandUpper",
               "; k* [GeV/#it{c}]; Entries", 3000, 0, 3);
  fHistKAntiProtonAntiSigmaSidebandUpper =
      new TH1F("fHistKAntiProtonAntiSigmaSidebandUpper",
               "; k* [GeV/#it{c}]; Entries", 3000, 0, 3);
  fHistKAntiProtonAntiSigmaMixedSidebandUpper =
      new TH1F("fHistKAntiProtonAntiSigmaMixedSidebandUpper",
               "; k* [GeV/#it{c}]; Entries", 3000, 0, 3);
  fHistograms->Add(fHistKProtonSigmaSidebandUpper);
  fHistograms->Add(fHistKProtonSigmaMixedSidebandUpper);
  fHistograms->Add(fHistKAntiProtonAntiSigmaSidebandUpper);
  fHistograms->Add(fHistKAntiProtonAntiSigmaMixedSidebandUpper);

  fHistKProtonSigmaVsLambda =
      new TH2F("fHistKProtonSigmaVsLambda",
               "; k*_{p#Sigma} [GeV/#it{c}]; k*_{p#Lambda} [GeV/#it{c}]", 3000,
               0, 3, 3000, 0, 3);
  fHistKAntiProtonAntiSigmaVsLambda = new TH2F(
      "fHistKAntiProtonAntiSigmaVsLambda",
      "; k*_{p#bar{#Sigma}} [GeV/#it{c}]; k*_{p#bar{#Lambda}} [GeV/#it{c}]",
      3000, 0, 3, 3000, 0, 3);
  fHistograms->Add(fHistKProtonSigmaVsLambda);
  fHistograms->Add(fHistKAntiProtonAntiSigmaVsLambda);

  fHistKProtonSigmaLambda = new TH1F("fHistKProtonSigmaLambda",
                                     "; k* [GeV/#it{c}]; Entries", 3000, 0, 3);
  fHistKProtonSigmaMixedLambda = new TH1F(
      "fHistKProtonSigmaMixedLambda", "; k* [GeV/#it{c}]; Entries", 3000, 0, 3);
  fHistKAntiProtonAntiSigmaLambda =
      new TH1F("fHistKAntiProtonAntiSigmaLambda", "; k* [GeV/#it{c}]; Entries",
               3000, 0, 3);
  fHistKAntiProtonAntiSigmaMixedLambda =
      new TH1F("fHistKAntiProtonAntiSigmaMixedLambda",
               "; k* [GeV/#it{c}]; Entries", 3000, 0, 3);
  fHistograms->Add(fHistKProtonSigmaLambda);
  fHistograms->Add(fHistKProtonSigmaMixedLambda);
  fHistograms->Add(fHistKAntiProtonAntiSigmaLambda);
  fHistograms->Add(fHistKAntiProtonAntiSigmaMixedLambda);

  fHistProtonLambdaNShared =
      new TH1F("fHistProtonLambdaNShared",
               "; # shared p-V0 pairs per event; Entries", 50, 0, 50);
  fHistProtonAntiLambdaNShared =
      new TH1F("fHistProtonAntiLambdaNShared",
               "; # shared p-#bar{V0} pairs per event; Entries", 50, 0, 50);
  fHistAntiProtonLambdaNShared =
      new TH1F("fHistAntiProtonLambdaNShared",
               "; # shared #bar{p}-V0 pairs per event; Entries", 50, 0, 50);
  fHistAntiProtonAntiLambdaNShared = new TH1F(
      "fHistProtonAntiAntiLambdaNShared",
      "; # shared #bar{p}-#bar{V0} pairs per event; Entries", 50, 0, 50);
  fHistLambdaNShared =
      new TH1F("fHistLambdaNShared", "; # shared V0 pairs per event; Entries",
               50, 0, 50);
  fHistAntiLambdaNShared =
      new TH1F("fHistAntiLambdaNShared",
               "; # shared #bar{V0} pairs per event; Entries", 50, 0, 50);
  fHistograms->Add(fHistProtonLambdaNShared);
  fHistograms->Add(fHistProtonAntiLambdaNShared);
  fHistograms->Add(fHistAntiProtonLambdaNShared);
  fHistograms->Add(fHistAntiProtonAntiLambdaNShared);
  fHistograms->Add(fHistLambdaNShared);
  fHistograms->Add(fHistAntiLambdaNShared);

  fHistSigmaMixedPt =
      new TH1F("fHistSigmaMixedPt",
               "; #it{p}_{T} #Lambda#gamma [GeV/c]; Entries", 1000, 0, 20);
  fHistSigmaMixedInvMass = new TH1F(
      "fHistSigmaMixedInvMass",
      "; Invariant mass #Lambda#gamma [GeV/#it{c}^{2}]; Entries", 2000, 1., 2.);
  fHistSigmaMixedInvMassPt =
      new TH2F("fHistSigmaMixedInvMassPt",
               "; #it{p}_{T} #Lambda#gamma [GeV/c]; Invariant mass "
               "#Lambda#gamma [GeV/#it{c}^{2}]",
               1000, 0, 20, 2000, 1., 2.);
  fHistSigmaMixedInvMassEta =
      new TH2F("fHistSigmaMixedInvMassEta",
               "; #eta; Invariant mass #Lambda#gamma [GeV/#it{c}^{2}]", 1000, 0,
               1, 2000, 1., 2.);
  fHistograms->Add(fHistSigmaMixedPt);
  fHistograms->Add(fHistSigmaMixedInvMass);
  fHistograms->Add(fHistSigmaMixedInvMassPt);
  fHistograms->Add(fHistSigmaMixedInvMassEta);

  fHistAntiSigmaMixedPt = new TH1F(
      "fHistAntiSigmaMixedPt",
      "; #it{p}_{T} #bar{#Lambda}#gamma [GeV/c]; Entries", 1000, 0, 20);
  fHistAntiSigmaMixedInvMass = new TH1F(
      "fHistAntiSigmaMixedInvMass",
      "; Invariant mass  #bar{#Lambda}#gamma [GeV/#it{c}^{2}]; Entries", 2000,
      1., 2.);
  fHistAntiSigmaMixedInvMassPt =
      new TH2F("fHistAntiSigmaMixedInvMassPt",
               "; #it{p}_{T}  #bar{#Lambda}#gamma [GeV/c]; Invariant mass "
               "#Lambda#gamma [GeV/#it{c}^{2}]",
               1000, 0, 20, 2000, 1., 2.);
  fHistAntiSigmaMixedInvMassEta =
      new TH2F("fHistAntiSigmaMixedInvMassEta",
               "; #eta; Invariant mass #bar{#Lambda}#gamma [GeV/#it{c}^{2}]",
               1000, 0, 1, 2000, 1., 2.);
  fHistograms->Add(fHistAntiSigmaMixedPt);
  fHistograms->Add(fHistAntiSigmaMixedInvMass);
  fHistograms->Add(fHistAntiSigmaMixedInvMassPt);
  fHistograms->Add(fHistAntiSigmaMixedInvMassEta);

  fHistProtonPt = new TH1F("fHistProtonPt", ";#it{p}_{T} [GeV/#it{c}]; Entries",
                           500, 0, 10);
  fHistAntiProtonPt = new TH1F("fHistAntiProtonPt",
                               ";#it{p}_{T} [GeV/#it{c}]; Entries", 500, 0, 10);
  fHistLambdaPt = new TH1F("fHistLambdaPt",
                           "; #it{p}_{T} [GeV/#it{c}]; Entries", 500, 0, 10);
  fHistLambdaInvMass = new TH1F(
      "fHistLambdaInvMass", "; Invariant mass p#pi [GeV/#it{c}^{2}]; Entries",
      600, 1.1, 1.13);
  fHistLambdaInvMassRec = new TH1F(
      "fHistLambdaInvMassRec",
      "; Invariant mass rec p#pi [GeV/#it{c}^{2}]; Entries", 600, 1.1, 1.13);
  fHistAntiLambdaPt = new TH1F(
      "fHistAntiLambdaPt", "; #it{p}_{T} [GeV/#it{c}]; Entries", 500, 0, 10);
  fHistAntiLambdaInvMass = new TH1F(
      "fHistAntiLambdaInvMass",
      "; Invariant mass p#pi [GeV/#it{c}^{2}]; Entries", 600, 1.1, 1.13);
  fHistAntiLambdaInvMassRec = new TH1F(
      "fHistAntiLambdaInvMassRec",
      "; Invariant mass rec p#pi [GeV/#it{c}^{2}]; Entries", 600, 1.1, 1.13);
  fHistSigmaPt =
      new TH1F("fHistSigmaPt", ";#it{p}_{T} [GeV/#it{c}]; Entries", 500, 0, 10);
  fHistSigmaInvMass =
      new TH1F("fHistSigmaInvMass",
               ";Invariant mass #Lambda#gamma [GeV/#it{c}^{2}]; Entries", 600,
               1.18, 1.21);
  fHistSigmaInvMassRec =
      new TH1F("fHistSigmaInvMassRec",
               ";Invariant mass rec #Lambda#gamma [GeV/#it{c}^{2}]; Entries",
               600, 1.18, 1.21);
  fHistAntiSigmaPt = new TH1F("fHistAntiSigmaPt",
                              ";#it{p}_{T} [GeV/#it{c}]; Entries", 500, 0, 10);
  fHistAntiSigmaInvMass =
      new TH1F("fHistAntiSigmaInvMass",
               ";Invariant mass #Lambda#gamma [GeV/#it{c}^{2}]; Entries", 600,
               1.18, 1.21);
  fHistAntiSigmaInvMassRec =
      new TH1F("fHistAntiSigmaInvMassRec",
               ";Invariant mass rec #Lambda#gamma [GeV/#it{c}^{2}]; Entries",
               600, 1.18, 1.21);
  fHistograms->Add(fHistProtonPt);
  fHistograms->Add(fHistAntiProtonPt);
  fHistograms->Add(fHistLambdaPt);
  fHistograms->Add(fHistLambdaInvMass);
  fHistograms->Add(fHistLambdaInvMassRec);
  fHistograms->Add(fHistAntiLambdaPt);
  fHistograms->Add(fHistAntiLambdaInvMass);
  fHistograms->Add(fHistAntiLambdaInvMassRec);
  fHistograms->Add(fHistSigmaPt);
  fHistograms->Add(fHistSigmaInvMass);
  fHistograms->Add(fHistSigmaInvMassRec);
  fHistograms->Add(fHistAntiSigmaPt);
  fHistograms->Add(fHistAntiSigmaInvMass);
  fHistograms->Add(fHistAntiSigmaInvMassRec);

  if (fIsExtendedQA) {
    if (fHistogramsExtended != nullptr) {
      delete fHistogramsExtended;
      fHistogramsExtended = nullptr;
    }
    if (fHistogramsExtended == nullptr) {
      fHistogramsExtended = new TList();
      fHistogramsExtended->SetOwner(kTRUE);
      fHistogramsExtended->SetName("EventCollection_Extended");
    }
    std::vector<int> fMultiplicityBins = {
        {1, 5, 9, 13, 17, 21, 25, 29, 33, 37, 41, 61, 80, 999}};

    for (int j = 0; j < static_cast<int>(fMultiplicityBins.size() - 1); ++j) {
      fHistKProtonProtonMulti[j] =
          new TH1F(Form("fHistKProtonProton_Multi_%i_%i", fMultiplicityBins[j],
                        fMultiplicityBins[j + 1] - 1),
                   Form("Multiplicity [%i, %i]; k* [GeV/#it{c}]; Entries",
                        fMultiplicityBins[j], fMultiplicityBins[j + 1] - 1),
                   3000, 0, 3);
      fHistKProtonProtonMixedMulti[j] =
          new TH1F(Form("fHistKProtonProtonMixed_Multi_%i_%i",
                        fMultiplicityBins[j], fMultiplicityBins[j + 1] - 1),
                   Form("Multiplicity [%i, %i]; k* [GeV/#it{c}]; Entries",
                        fMultiplicityBins[j], fMultiplicityBins[j + 1] - 1),
                   3000, 0, 3);
      fHistKAntiProtonAntiProtonMulti[j] =
          new TH1F(Form("fHistKAntiProtonAntiProton_Multi_%i_%i",
                        fMultiplicityBins[j], fMultiplicityBins[j + 1] - 1),
                   Form("Multiplicity [%i, %i]; k* [GeV/#it{c}]; Entries",
                        fMultiplicityBins[j], fMultiplicityBins[j + 1] - 1),
                   3000, 0, 3);
      fHistKAntiProtonAntiProtonMixedMulti[j] =
          new TH1F(Form("fHistKAntiProtonAntiProtonMixed_Multi_%i_%i",
                        fMultiplicityBins[j], fMultiplicityBins[j + 1] - 1),
                   Form("Multiplicity [%i, %i]; k* [GeV/#it{c}]; Entries",
                        fMultiplicityBins[j], fMultiplicityBins[j + 1] - 1),
                   3000, 0, 3);
      fHistogramsExtended->Add(fHistKProtonProtonMulti[j]);
      fHistogramsExtended->Add(fHistKProtonProtonMixedMulti[j]);
      fHistogramsExtended->Add(fHistKAntiProtonAntiProtonMulti[j]);
      fHistogramsExtended->Add(fHistKAntiProtonAntiProtonMixedMulti[j]);
      fHistKProtonLambdaMulti[j] =
          new TH1F(Form("fHistKProtonLambda_Multi_%i_%i", fMultiplicityBins[j],
                        fMultiplicityBins[j + 1] - 1),
                   Form("Multiplicity [%i, %i]; k* [GeV/#it{c}]; Entries",
                        fMultiplicityBins[j], fMultiplicityBins[j + 1] - 1),
                   3000, 0, 3);
      fHistKProtonLambdaMixedMulti[j] =
          new TH1F(Form("fHistKProtonLambdaMixed_Multi_%i_%i",
                        fMultiplicityBins[j], fMultiplicityBins[j + 1] - 1),
                   Form("Multiplicity [%i, %i]; k* [GeV/#it{c}]; Entries",
                        fMultiplicityBins[j], fMultiplicityBins[j + 1] - 1),
                   3000, 0, 3);
      fHistKAntiProtonAntiLambdaMulti[j] =
          new TH1F(Form("fHistKAntiProtonAntiLambda_Multi_%i_%i",
                        fMultiplicityBins[j], fMultiplicityBins[j + 1] - 1),
                   Form("Multiplicity [%i, %i]; k* [GeV/#it{c}]; Entries",
                        fMultiplicityBins[j], fMultiplicityBins[j + 1] - 1),
                   3000, 0, 3);
      fHistKAntiProtonAntiLambdaMixedMulti[j] =
          new TH1F(Form("fHistKAntiProtonAntiLambdaMixed_Multi_%i_%i",
                        fMultiplicityBins[j], fMultiplicityBins[j + 1] - 1),
                   Form("Multiplicity [%i, %i]; k* [GeV/#it{c}]; Entries",
                        fMultiplicityBins[j], fMultiplicityBins[j + 1] - 1),
                   3000, 0, 3);
      fHistogramsExtended->Add(fHistKProtonLambdaMulti[j]);
      fHistogramsExtended->Add(fHistKProtonLambdaMixedMulti[j]);
      fHistogramsExtended->Add(fHistKAntiProtonAntiLambdaMulti[j]);
      fHistogramsExtended->Add(fHistKAntiProtonAntiLambdaMixedMulti[j]);
      fHistKLambdaLambdaMulti[j] =
          new TH1F(Form("fHistKLambdaLambda_Multi_%i_%i", fMultiplicityBins[j],
                        fMultiplicityBins[j + 1] - 1),
                   Form("Multiplicity [%i, %i]; k* [GeV/#it{c}]; Entries",
                        fMultiplicityBins[j], fMultiplicityBins[j + 1] - 1),
                   3000, 0, 3);
      fHistKLambdaLambdaMixedMulti[j] =
          new TH1F(Form("fHistKLambdaLambdaMixed_Multi_%i_%i",
                        fMultiplicityBins[j], fMultiplicityBins[j + 1] - 1),
                   Form("Multiplicity [%i, %i]; k* [GeV/#it{c}]; Entries",
                        fMultiplicityBins[j], fMultiplicityBins[j + 1] - 1),
                   3000, 0, 3);
      fHistKAntiLambdaAntiLambdaMulti[j] =
          new TH1F(Form("fHistKAntiLambdaAntiLambda_Multi_%i_%i",
                        fMultiplicityBins[j], fMultiplicityBins[j + 1] - 1),
                   Form("Multiplicity [%i, %i]; k* [GeV/#it{c}]; Entries",
                        fMultiplicityBins[j], fMultiplicityBins[j + 1] - 1),
                   3000, 0, 3);
      fHistKAntiLambdaAntiLambdaMixedMulti[j] =
          new TH1F(Form("fHistKAntiLambdaAntiLambdaMixed_Multi_%i_%i",
                        fMultiplicityBins[j], fMultiplicityBins[j + 1] - 1),
                   Form("Multiplicity [%i, %i]; k* [GeV/#it{c}]; Entries",
                        fMultiplicityBins[j], fMultiplicityBins[j + 1] - 1),
                   3000, 0, 3);
      fHistogramsExtended->Add(fHistKLambdaLambdaMulti[j]);
      fHistogramsExtended->Add(fHistKLambdaLambdaMixedMulti[j]);
      fHistogramsExtended->Add(fHistKAntiLambdaAntiLambdaMulti[j]);
      fHistogramsExtended->Add(fHistKAntiLambdaAntiLambdaMixedMulti[j]);
      fHistKProtonSigmaMulti[j] =
          new TH1F(Form("fHistKProtonSigma_Multi_%i_%i", fMultiplicityBins[j],
                        fMultiplicityBins[j + 1] - 1),
                   Form("Multiplicity [%i, %i]; k* [GeV/#it{c}]; Entries",
                        fMultiplicityBins[j], fMultiplicityBins[j + 1] - 1),
                   3000, 0, 3);
      fHistKProtonSigmaMixedMulti[j] =
          new TH1F(Form("fHistKProtonSigmaMixed_Multi_%i_%i",
                        fMultiplicityBins[j], fMultiplicityBins[j + 1] - 1),
                   Form("Multiplicity [%i, %i]; k* [GeV/#it{c}]; Entries",
                        fMultiplicityBins[j], fMultiplicityBins[j + 1] - 1),
                   3000, 0, 3);
      fHistKAntiProtonAntiSigmaMulti[j] =
          new TH1F(Form("fHistKAntiProtonAntiSigma_Multi_%i_%i",
                        fMultiplicityBins[j], fMultiplicityBins[j + 1] - 1),
                   Form("Multiplicity [%i, %i]; k* [GeV/#it{c}]; Entries",
                        fMultiplicityBins[j], fMultiplicityBins[j + 1] - 1),
                   3000, 0, 3);
      fHistKAntiProtonAntiSigmaMixedMulti[j] =
          new TH1F(Form("fHistKAntiProtonAntiSigmaMixed_Multi_%i_%i",
                        fMultiplicityBins[j], fMultiplicityBins[j + 1] - 1),
                   Form("Multiplicity [%i, %i]; k* [GeV/#it{c}]; Entries",
                        fMultiplicityBins[j], fMultiplicityBins[j + 1] - 1),
                   3000, 0, 3);
      fHistogramsExtended->Add(fHistKProtonSigmaMulti[j]);
      fHistogramsExtended->Add(fHistKProtonSigmaMixedMulti[j]);
      fHistogramsExtended->Add(fHistKAntiProtonAntiSigmaMulti[j]);
      fHistogramsExtended->Add(fHistKAntiProtonAntiSigmaMixedMulti[j]);
    }

    for (int i = 0; i < 9; ++i) {
      fHistProtonProtonPhiStar[i] = new TH2F(
          Form("fHistProtonProtonPhiStar_%i", 85 + i * 20),
          Form("r_{TPC} = %i cm; #Delta #eta; #Delta #phi*", 85 + i * 20), 200,
          0, 0.5, 200., 0, 0.5);
      fHistAntiProtonAntiProtonPhiStar[i] = new TH2F(
          Form("fHistAntiProtonAntiProtonPhiStar_%i", 85 + i * 20),
          Form("r_{TPC} = %i cm; #Delta #eta; #Delta #phi*", 85 + i * 20), 200,
          0, 0.5, 200., 0, 0.5);
      fHistProtonLambdaProtonPhiStar[i] = new TH2F(
          Form("fHistProtonLambdaProtonPhiStar_%i", 85 + i * 20),
          Form("r_{TPC} = %i cm; #Delta #eta; #Delta #phi*", 85 + i * 20), 200,
          0, 0.5, 200., 0, 0.5);
      fHistAntiProtonAntiLambdaProtonPhiStar[i] = new TH2F(
          Form("fHistAntiProtonAntiLambdaProtonPhiStar_%i", 85 + i * 20),
          Form("r_{TPC} = %i cm; #Delta #eta; #Delta #phi*", 85 + i * 20), 200,
          0, 0.5, 200., 0, 0.5);
      fHistProtonLambdaPionPhiStar[i] = new TH2F(
          Form("fHistProtonLambdaPionPhiStar_%i", 85 + i * 20),
          Form("r_{TPC} = %i cm; #Delta #eta; #Delta #phi*", 85 + i * 20), 200,
          0, 0.5, 200., 0, 0.5);
      fHistAntiProtonAntiLambdaPionPhiStar[i] = new TH2F(
          Form("fHistAntiProtonAntiLambdaPionPhiStar_%i", 85 + i * 20),
          Form("r_{TPC} = %i cm; #Delta #eta; #Delta #phi*", 85 + i * 20), 200,
          0, 0.5, 200., 0, 0.5);
      fHistProtonProtonPhiStarMixed[i] = new TH2F(
          Form("fHistProtonProtonPhiStarMixed_%i", 85 + i * 20),
          Form("r_{TPC} = %i cm; #Delta #eta; #Delta #phi*", 85 + i * 20), 200,
          0, 0.5, 200., 0, 0.5);
      fHistAntiProtonAntiProtonPhiStarMixed[i] = new TH2F(
          Form("fHistAntiProtonAntiProtonPhiStarMixed_%i", 85 + i * 20),
          Form("r_{TPC} = %i cm; #Delta #eta; #Delta #phi*", 85 + i * 20), 200,
          0, 0.5, 200., 0, 0.5);
      fHistProtonLambdaProtonPhiStarMixed[i] = new TH2F(
          Form("fHistProtonLambdaProtonPhiStarMixed_%i", 85 + i * 20),
          Form("r_{TPC} = %i cm; #Delta #eta; #Delta #phi*", 85 + i * 20), 200,
          0, 0.5, 200., 0, 0.5);
      fHistAntiProtonAntiLambdaProtonPhiStarMixed[i] = new TH2F(
          Form("fHistAntiProtonAntiLambdaProtonPhiStarMixed_%i", 85 + i * 20),
          Form("r_{TPC} = %i cm; #Delta #eta; #Delta #phi*", 85 + i * 20), 200,
          0, 0.5, 200., 0, 0.5);
      fHistProtonLambdaPionPhiStarMixed[i] = new TH2F(
          Form("fHistProtonLambdaPionPhiStarMixed_%i", 85 + i * 20),
          Form("r_{TPC} = %i cm; #Delta #eta; #Delta #phi*", 85 + i * 20), 200,
          0, 0.5, 200., 0, 0.5);
      fHistAntiProtonAntiLambdaPionPhiStarMixed[i] = new TH2F(
          Form("fHistAntiProtonAntiLambdaPionPhiStarMixed_%i", 85 + i * 20),
          Form("r_{TPC} = %i cm; #Delta #eta; #Delta #phi*", 85 + i * 20), 200,
          0, 0.5, 200., 0, 0.5);

      fHistogramsExtended->Add(fHistProtonProtonPhiStar[i]);
      fHistogramsExtended->Add(fHistAntiProtonAntiProtonPhiStar[i]);
      fHistogramsExtended->Add(fHistProtonLambdaProtonPhiStar[i]);
      fHistogramsExtended->Add(fHistAntiProtonAntiLambdaProtonPhiStar[i]);
      fHistogramsExtended->Add(fHistProtonLambdaPionPhiStar[i]);
      fHistogramsExtended->Add(fHistAntiProtonAntiLambdaPionPhiStar[i]);
      fHistogramsExtended->Add(fHistProtonProtonPhiStarMixed[i]);
      fHistogramsExtended->Add(fHistAntiProtonAntiProtonPhiStarMixed[i]);
      fHistogramsExtended->Add(fHistProtonLambdaProtonPhiStarMixed[i]);
      fHistogramsExtended->Add(fHistAntiProtonAntiLambdaProtonPhiStarMixed[i]);
      fHistogramsExtended->Add(fHistProtonLambdaPionPhiStarMixed[i]);
      fHistogramsExtended->Add(fHistAntiProtonAntiLambdaPionPhiStarMixed[i]);
    }
    fHistograms->Add(fHistogramsExtended);
  }

  if (fIsMC) {
    if (fHistogramsMC != nullptr) {
      delete fHistogramsMC;
      fHistogramsMC = nullptr;
    }
    if (fHistogramsMC == nullptr) {
      fHistogramsMC = new TList();
      fHistogramsMC->SetOwner(kTRUE);
      fHistogramsMC->SetName("EventCollection_MC");
    }

    fHistMCProtonProtonMomentum =
        new TH2F("fHistMCProtonProtonMomentum",
                 "; k*_{truth} [GeV/#it{c}]; k*_{rec} [GeV/#it{c}]", 750, 0.,
                 3., 750, 0., 3.);
    fHistMCAntiProtonAntiProtonMomentum =
        new TH2F("fHistMCAntiProtonAntiProtonMomentum",
                 "; k*_{truth} [GeV/#it{c}]; k*_{rec} [GeV/#it{c}]", 750, 0.,
                 3., 750, 0., 3.);
    fHistMCProtonLambdaMomentum =
        new TH2F("fHistMCProtonLambdaMomentum",
                 "; k*_{truth} [GeV/#it{c}]; k*_{rec} [GeV/#it{c}]", 600, 0.,
                 3., 600, 0., 3.);
    fHistMCAntiProtonAntiLambdaMomentum =
        new TH2F("fHistMCAntiProtonAntiLambdaMomentum",
                 "; k*_{truth} [GeV/#it{c}]; k*_{rec} [GeV/#it{c}]", 600, 0.,
                 3., 600, 0., 3.);
    fHistMCLambdaLambdaMomentum =
        new TH2F("fHistMCLambdaLambdaMomentum",
                 "; k*_{truth} [GeV/#it{c}]; k*_{rec} [GeV/#it{c}]", 600, 0.,
                 3., 600, 0., 3.);
    fHistMCAntiLambdaAntiLambdaMomentum =
        new TH2F("fHistMCAntiLambdaAntiLambdaMomentum",
                 "; k*_{truth} [GeV/#it{c}]; k*_{rec} [GeV/#it{c}]", 600, 0.,
                 3., 600, 0., 3.);
    fHistMCProtonSigmaMomentum =
        new TH2F("fHistMCProtonSigmaMomentum",
                 "; k*_{truth} [GeV/#it{c}]; k*_{rec} [GeV/#it{c}]", 600, 0.,
                 3., 600, 0., 3.);
    fHistMCAntiProtonAntiSigmaMomentum =
        new TH2F("fHistMCAntiProtonAntiSigmaMomentum",
                 "; k*_{truth} [GeV/#it{c}]; k*_{rec} [GeV/#it{c}]", 600, 0.,
                 3., 600, 0., 3.);
    fHistogramsMC->Add(fHistMCProtonProtonMomentum);
    fHistogramsMC->Add(fHistMCAntiProtonAntiProtonMomentum);
    fHistogramsMC->Add(fHistMCProtonLambdaMomentum);
    fHistogramsMC->Add(fHistMCAntiProtonAntiLambdaMomentum);
    fHistogramsMC->Add(fHistMCLambdaLambdaMomentum);
    fHistogramsMC->Add(fHistMCAntiLambdaAntiLambdaMomentum);
    fHistogramsMC->Add(fHistMCProtonSigmaMomentum);
    fHistogramsMC->Add(fHistMCAntiProtonAntiSigmaMomentum);

    fHistograms->Add(fHistogramsMC);
  }
}
