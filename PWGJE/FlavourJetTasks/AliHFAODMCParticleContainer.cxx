/*************************************************************************
* Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
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
#include <AliAODMCParticle.h>
#include <AliLog.h>

#include "AliHFAODMCParticleContainer.h"

/// \cond CLASSIMP
ClassImp(AliHFAODMCParticleContainer)
/// \endcond

/// This is the default constructor, used for ROOT I/O purposes.
AliHFAODMCParticleContainer::AliHFAODMCParticleContainer() :
  AliMCParticleContainer(),
  fSpecialPDG(0),
  fRejectedOrigin(0),
  fAcceptedDecay(0),
  fRejectISR(kFALSE),
  fHistOrigin(0)
{
  // Constructor.
}

/// This is the standard named constructor.
///
/// \param name Name of the particle collection
AliHFAODMCParticleContainer::AliHFAODMCParticleContainer(const char *name) :
  AliMCParticleContainer(name),
  fSpecialPDG(0),
  fRejectedOrigin(AliAnalysisTaskDmesonJets::kUnknownQuark | AliAnalysisTaskDmesonJets::kFromBottom),
  fAcceptedDecay(AliAnalysisTaskDmesonJets::kAnyDecay),
  fRejectISR(kFALSE),
  fHistOrigin(0)
{
  // Constructor.
}


/// Automatically sets parameters to select only the decay chain c->D0->Kpi
void AliHFAODMCParticleContainer::SelectCharmtoD0toKpi()
{
  SetSpecialPDG(421);
  SetKeepOnlyD0toKpi();
  SetRejectDfromB(kTRUE);
  SetRejectQuarkNotFound(kTRUE);
  SetKeepOnlyDfromB(kFALSE);
}

/// Automatically sets parameters to select only the decay chain c->D*->Kpipi
void AliHFAODMCParticleContainer::SelectCharmtoDStartoKpipi()
{
  SetSpecialPDG(413);
  SetKeepOnlyDStartoKpipi();
  SetRejectDfromB(kTRUE);
  SetRejectQuarkNotFound(kTRUE);
  SetKeepOnlyDfromB(kFALSE);
}

/// First check whether the particle is a daughter of a "special" PDG particle (in which case the particle is rejected);
/// if the particle itself is "special" PDG particle, skips the regular MC particle cuts and apply only the general particle and
/// kinematic cuts.
///
/// \param i Index of the particle to be checked.
///
/// \return kTRUE if the particle is accepted, kFALSE otherwise.
Bool_t AliHFAODMCParticleContainer::AcceptMCParticle(const AliAODMCParticle *vp, UInt_t &rejectionReason) const
{
  // Return true if vp is accepted.

  if (IsSpecialPDGDaughter(vp)) {
    rejectionReason |= kHFCut;
    return kFALSE;  // daughter of a special PDG particle, reject it without any other check.
  }

  if (IsSpecialPDG(vp)) {
    // Special PDG particle, skip regular MC particle cuts and apply kinematic cuts.
    // For the future, may want to implement special kinematic cuts for D mesons
    AliTLorentzVector mom;
    GetMomentumFromParticle(mom, vp);
    return ApplyKinematicCuts(mom, rejectionReason);
  }
  else {
    return AliMCParticleContainer::AcceptMCParticle(vp, rejectionReason);
  }
}

/// First check whether the particle is a daughter of a "special" PDG particle (in which case the particle is rejected);
/// if the particle itself is "special" PDG particle, skips the regular MC particle cuts and apply only the general particle and
/// kinematic cuts.
///
/// \param i Index of the particle to be checked.
///
/// \return kTRUE if the particle is accepted, kFALSE otherwise.
Bool_t AliHFAODMCParticleContainer::AcceptMCParticle(Int_t i, UInt_t &rejectionReason) const
{
  // Return true if vp is accepted.

  AliAODMCParticle* vp = GetMCParticle(i);

  if (IsSpecialPDGDaughter(vp)) {
    rejectionReason = kHFCut;
    return kFALSE;  // daughter of a special PDG particle, reject it without any other check.
  }

  if (IsSpecialPDG(vp, fHistOrigin)) {
    // Special PDG particle, skip regular MC particle cuts and apply particle cuts.
    // For the future, may want to implement special kinematic cuts for D mesons
    AliTLorentzVector mom;
    GetMomentum(mom, i);
    return ApplyKinematicCuts(mom, rejectionReason);
  }
  else {
    return AliMCParticleContainer::AcceptMCParticle(i, rejectionReason);
  }
}

/// Check if particle it's a "special" PDG particle: AOD mode
/// \param part Pointer to a valid AliAODMCParticle object
///
/// \result kTRUE if it is a "special" PDG particle, kFALSE otherwise
Bool_t AliHFAODMCParticleContainer::IsSpecialPDGDaughter(const AliAODMCParticle* part) const
{
  if (fSpecialPDG == 0) return kFALSE;
  
  const AliAODMCParticle* pm = part;
  Int_t imo = -1;
  while (pm != 0) {
    imo = pm->GetMother();
    if (imo < 0) break;
    pm = static_cast<const AliAODMCParticle*>(fClArray->At(imo));
    if (IsSpecialPDG(pm)) {
      AliDebug(2, Form("Rejecting particle (PDG = %d, pT = %.3f, eta = %.3f, phi = %.3f) daughter of %d (PDG = %d, pT = %.3f, eta = %.3f, phi = %.3f)",
                       part->PdgCode(), part->Pt(), part->Eta(), part->Phi(), imo, pm->PdgCode(), pm->Pt(), pm->Eta(), pm->Phi()));
      return kTRUE;
    }
  }
  return kFALSE;
}

/// Check if particle it's a daughter of a "special" PDG particle: AOD mode
/// \param part Pointer to a valid AliAODMCParticle object
///
/// \result kTRUE if it is a daughter of the "special" PDG particle, kFALSE otherwise
Bool_t AliHFAODMCParticleContainer::IsSpecialPDG(const AliAODMCParticle* part, TH1* histOrigin) const
{
  if (fSpecialPDG == 0) return kFALSE;

  Int_t partPdgCode = TMath::Abs(part->PdgCode());

  if (partPdgCode != fSpecialPDG) return kFALSE;

  if (!part->IsPrimary()) return kFALSE;

  if (fRejectISR) {
    // proton has PDG code 2212
    std::set<UInt_t> pdgSet = {2212};
    auto origin = AliAnalysisTaskDmesonJets::AnalysisEngine::FindParticleOrigin(part, fClArray, AliAnalysisTaskDmesonJets::AnalysisEngine::kFindFirst, pdgSet);
    if (origin) return kFALSE;
  }

  auto origin = AliAnalysisTaskDmesonJets::AnalysisEngine::IsPromptCharm(part, fClArray);

  if (histOrigin) {
    UInt_t rs = origin.first;
    UShort_t p = 0;
    while (rs >>= 1) { p++; }
    histOrigin->Fill(p);
  }

  if ((origin.first & fRejectedOrigin) != 0) return kFALSE;

  AliAnalysisTaskDmesonJets::EMesonDecayChannel_t decayChannel = AliAnalysisTaskDmesonJets::AnalysisEngine::CheckDecayChannel(part, fClArray);

  if (fAcceptedDecay && (decayChannel & fAcceptedDecay) == 0) return kFALSE;

  AliDebug(2, Form("Including particle %d (PDG = %d, pT = %.3f, eta = %.3f, phi = %.3f)",
      part->Label(), partPdgCode, part->Pt(), part->Eta(), part->Phi()));

  // Special PDG particle
  return kTRUE;
}

/// Checks whether the special PDG particle is found in the container
///
/// \result kTRUE if the special PDG particle is found
Bool_t AliHFAODMCParticleContainer::IsSpecialPDGFound() const
{
  for (auto part : accepted()) {
    if (IsSpecialPDG(part)) return kTRUE;
  }
  return kFALSE;
}
