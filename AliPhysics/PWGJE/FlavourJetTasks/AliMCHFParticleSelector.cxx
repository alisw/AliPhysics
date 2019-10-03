//
// Class to select particles in MC events.
//
// Author: S. Aiola

#include "AliMCHFParticleSelector.h"

#include <TClonesArray.h>

#include <AliVEvent.h>
#include <AliMCEvent.h>
#include <AliAODMCParticle.h>
#include <AliMCParticle.h>
#include <AliStack.h>
#include <AliLog.h>

#include "AliAnalysisTaskSEDmesonsFilterCJ.h"

ClassImp(AliMCHFParticleSelector)

//________________________________________________________________________
AliMCHFParticleSelector::AliMCHFParticleSelector() : 
  AliEmcalMCTrackSelector("AliMCHFParticleSelector"),
  fSpecialPDG(0),
  fRejectQuarkNotFound(kTRUE),
  fRejectDfromB(kTRUE),
  fKeepOnlyDfromB(kFALSE),
  fKeepOnlyD0toKpi(kFALSE),
  fKeepOnlyDStartoKpipi(kFALSE)
{
  // Constructor.
}

//________________________________________________________________________
AliMCHFParticleSelector::AliMCHFParticleSelector(const char *name) : 
  AliEmcalMCTrackSelector(name),
  fSpecialPDG(0),
  fRejectQuarkNotFound(kTRUE),
  fRejectDfromB(kTRUE),
  fKeepOnlyDfromB(kFALSE),
  fKeepOnlyD0toKpi(kFALSE),
  fKeepOnlyDStartoKpipi(kFALSE)
{
  // Constructor.
}

//________________________________________________________________________
AliMCHFParticleSelector::~AliMCHFParticleSelector()
{
  // Destructor.
}

//________________________________________________________________________
void AliMCHFParticleSelector::SelectCharmtoD0toKpi()
{
  SetSpecialPDG(421);
  SetKeepOnlyD0toKpi(kTRUE);
  SetKeepOnlyDStartoKpipi(kFALSE);
  SetRejectDfromB(kTRUE);
  SetRejectQuarkNotFound(kTRUE);
  SetKeepOnlyDfromB(kFALSE);
}

//________________________________________________________________________
void AliMCHFParticleSelector::SelectCharmtoDStartoKpipi()
{
  SetSpecialPDG(413);
  SetKeepOnlyD0toKpi(kFALSE);
  SetKeepOnlyDStartoKpipi(kTRUE);
  SetRejectDfromB(kTRUE);
  SetRejectQuarkNotFound(kTRUE);
  SetKeepOnlyDfromB(kFALSE);
}

//________________________________________________________________________
Bool_t AliMCHFParticleSelector::AcceptParticle(AliAODMCParticle* part) const
{
  // Determine whether the MC particle is accepted.

  if (!part) return kFALSE;

  Int_t partPdgCode = TMath::Abs(part->PdgCode());

  if (fOnlyHIJING && (part->GetGeneratorIndex() != 0)) return kFALSE;

  if (fEtaMax > 0. && TMath::Abs(part->Eta()) > fEtaMax) return kFALSE;
    
  if (fRejectNK && (partPdgCode == 130 || partPdgCode == 2112)) return kFALSE;

  Bool_t isSpecialPdg = (fSpecialPDG != 0 && partPdgCode == fSpecialPDG && part->IsPrimary());

  if (isSpecialPdg) {
    Int_t origin = -1;

    if (fIsESD) {
      origin = AliAnalysisTaskSEDmesonsFilterCJ::CheckOrigin(part->Label(), fMC->Stack());
    }
    else {
      origin = AliAnalysisTaskSEDmesonsFilterCJ::CheckOrigin(part, fParticlesIn);
    }

    if (origin < 0) isSpecialPdg = kFALSE;
    
    if (fRejectQuarkNotFound && origin == AliAnalysisTaskSEDmesonsFilterCJ::kQuarkNotFound) {
      isSpecialPdg = kFALSE;
    }
    else if (fRejectDfromB && origin == AliAnalysisTaskSEDmesonsFilterCJ::kFromBottom) {
      isSpecialPdg = kFALSE;
    }
    else if (fKeepOnlyDfromB && origin != AliAnalysisTaskSEDmesonsFilterCJ::kFromBottom) {
      isSpecialPdg = kFALSE;
    }

    Int_t decayChannel = -1;

    if (fIsESD) {
      decayChannel = AliAnalysisTaskSEDmesonsFilterCJ::CheckDecayChannel(part->Label(), fMC->Stack());
    }
    else {
      decayChannel = AliAnalysisTaskSEDmesonsFilterCJ::CheckDecayChannel(part, fParticlesIn);
    }

    if (fKeepOnlyD0toKpi && decayChannel != AliAnalysisTaskSEDmesonsFilterCJ::kDecayD0toKpi) {
      isSpecialPdg = kFALSE;
    }
    else if (fKeepOnlyDStartoKpipi && decayChannel != AliAnalysisTaskSEDmesonsFilterCJ::kDecayDStartoKpipi) {
      isSpecialPdg = kFALSE;
    }
  }

  if (isSpecialPdg) {
    AliDebug(2, Form("Including particle %d (PDG = %d, pT = %.3f, eta = %.3f, phi = %.3f)",
                     part->Label(), partPdgCode, part->Pt(), part->Eta(), part->Phi()));
  }

  if (!isSpecialPdg) {
    if (fChargedMC && part->Charge() == 0) return kFALSE;
    if (fOnlyPhysPrim && !part->IsPhysicalPrimary()) return kFALSE;
  }

  if (fIsESD) {
    if (IsSpecialPDGDaughter(part->Label())) return kFALSE;
  }
  else {
    if (IsSpecialPDGDaughter(part)) return kFALSE;
  }

  return kTRUE;
}

//________________________________________________________________________
Bool_t AliMCHFParticleSelector::IsSpecialPDGDaughter(AliAODMCParticle* part) const
{
  // Check if particle it's a daughter of a "special" PDG particle: AOD mode

  if (fSpecialPDG == 0) return kTRUE;
  
  AliAODMCParticle* pm = part;
  Int_t imo = -1;
  while (pm != 0) {
    imo = pm->GetMother();
    if (imo < 0) break;
    pm = static_cast<AliAODMCParticle*>(fParticlesIn->At(imo));
    if (TMath::Abs(pm->GetPdgCode()) == fSpecialPDG && pm->IsPrimary()) {
      AliDebug(2, Form("Rejecting particle (PDG = %d, pT = %.3f, eta = %.3f, phi = %.3f) daughter of %d (PDG = %d, pT = %.3f, eta = %.3f, phi = %.3f)",
                       part->PdgCode(), part->Pt(), part->Eta(), part->Phi(), imo, pm->PdgCode(), pm->Pt(), pm->Eta(), pm->Phi()));
      return kTRUE;
    }
  }
  return kFALSE;
}

//________________________________________________________________________
Bool_t AliMCHFParticleSelector::IsSpecialPDGDaughter(Int_t iPart) const
{
  // Check if particle it's a daughter of a "special" PDG particle: ESD mode
  
  if (fSpecialPDG == 0) return kTRUE;
  
  AliStack* stack = fMC->Stack();
  TParticle* part = stack->Particle(iPart);
  TParticle* pm = part;
  Int_t imo = -1;
  while (pm != 0) {
    imo = pm->GetFirstMother();
    if (imo < 0) break;
    pm = stack->Particle(imo);
    if (TMath::Abs(pm->GetPdgCode()) == fSpecialPDG && imo < stack->GetNprimary()) {
      AliDebug(2, Form("Rejecting particle (PDG = %d, pT = %.3f, eta = %.3f, phi = %.3f) daughter of %d (PDG = %d, pT = %.3f, eta = %.3f, phi = %.3f)",
                       part->GetPdgCode(), part->Pt(), part->Eta(), part->Phi(), imo, pm->GetPdgCode(), pm->Pt(), pm->Eta(), pm->Phi()));
      return kTRUE;
    } 
  }
  return kFALSE;
}
