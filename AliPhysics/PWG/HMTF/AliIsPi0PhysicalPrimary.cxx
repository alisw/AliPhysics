#include "TParticle.h"

#include "AliIsPi0PhysicalPrimary.h"

Bool_t AliIsPi0PhysicalPrimary(Int_t index, AliStack *stack)
{
  // pi0 are considered unstable and thus not primary in AliStack.
  // If the given track index is a pi0 and if that pi0 has no ancestors, it is a primary
  // Return kFALSE if the index points to anything but a pi0

  TParticle* p = stack->Particle(index);
  Int_t ist = p->GetStatusCode();

  // Initial state particle
  //if (ist > 1) return kFALSE; // pi0 are not initial state (?)
  // pi0 are unstable and thus this returned false in the original implementation.
  //if (!stack->IsStable(pdg)) return kFALSE;

  Int_t pdg = TMath::Abs(p->GetPdgCode());

  // The function is only for pi0's so I'm out of here if its not a pi'0 we are looking at
  if (pdg != 111) return kFALSE;

  if (index < stack->GetNprimary()) {
    // Particle produced by generator
    return kTRUE;
  }
  else {
    // Particle produced during transport
    Int_t imo =  p->GetFirstMother();
    TParticle* pm  = stack->Particle(imo);
    Int_t mpdg = TMath::Abs(pm->GetPdgCode());
    // Check for Sigma0
    if ((mpdg == 3212) &&  (imo <  stack->GetNprimary())) return kTRUE;

    // Check if it comes from a pi0 decay
    if ((mpdg == 111) && (imo < stack->GetNprimary()))   return kTRUE;

    // Check if this is a heavy flavor decay product
    Int_t mfl  = Int_t (mpdg / TMath::Power(10, Int_t(TMath::Log10(mpdg))));

    // Light hadron
    if (mfl < 4) return kFALSE;

    // Heavy flavor hadron produced by generator
    if (imo <  stack->GetNprimary()) {
      return kTRUE;
    }

    // To be sure that heavy flavor has not been produced in a secondary interaction
    // Loop back to the generated mother
    while (imo >=  stack->GetNprimary()) {
      imo = pm->GetFirstMother();
      pm  =  stack->Particle(imo);
    }
    mpdg = TMath::Abs(pm->GetPdgCode());
    mfl  = Int_t (mpdg / TMath::Power(10, Int_t(TMath::Log10(mpdg))));

    if (mfl < 4) {
      return kFALSE;
    } else {
      return kTRUE;
    }
    } // produced by generator ?
}
