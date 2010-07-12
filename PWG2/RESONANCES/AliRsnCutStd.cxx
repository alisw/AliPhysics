//
// Class AliRsnCutStd
//
// General implementation of a single cut strategy, which can be:
// - a value contained in a given interval  [--> IsBetween()   ]
// - a value equal to a given reference     [--> MatchesValue()]
//
// In all cases, the reference value(s) is (are) given as data members
// and each kind of cut requires a given value type (Int, UInt, Double),
// but the cut check procedure is then automatized and chosen thanks to
// an enumeration of the implemented cut types.
// At the end, the user (or any other point which uses this object) has
// to use the method IsSelected() to check if this cut has been passed.
//
// authors: Martin Vala (martin.vala@cern.ch)
//          Alberto Pulvirenti (alberto.pulvirenti@ct.infn.it)
//

#include "TParticle.h"
#include "TMath.h"

#include "AliRsnDaughter.h"
#include "AliRsnEvent.h"
#include "AliRsnPairParticle.h"

#include "AliRsnCutStd.h"

ClassImp(AliRsnCutStd)

//_________________________________________________________________________________________________
AliRsnCutStd::AliRsnCutStd() :
    AliRsnCut(),
    fType(kLastType),
    fUseMC(kFALSE),
    fMass(0.0)
{
//
// Default constructor.
//
}

//_________________________________________________________________________________________________
AliRsnCutStd::AliRsnCutStd
(const char *name, EType type, Int_t val1, Int_t val2, Bool_t useMC) :
    AliRsnCut(name, val1, val2),
    fType(type),
    fUseMC(useMC),
    fMass(0.0)
{
//
// Main constructor.
// Checks also that cut values are given in the correct type,
// in order to avoid that one passes, for example, a value which should be double
// but is interpreted as integer due to the overloading of constructors
//

  switch (fType) {
    // int
  case kMult:
  case kMultDiff:
  case kKink:
  case kKinkMother:
  case kAssignedPID:
  case kTruePID:
  case kRequiredPID:
  case kCharge:
    break;
    // ulong
  case kStatus:
    if (fVarType != kULong) {
      AliWarning(Form("[INT CONSTRUCTOR] Cut '%s' is based on ULONG. Casting values to ULONG", GetName()));
      SetRange((ULong_t)val1, (ULong_t)val2);
      AliWarning(Form("[INT CONSTRUCTOR] Cut '%s' ULONG range = %lu, %lu", GetName(), fMinU, fMaxU));
    }
    break;
    // double
  case kP:
  case kPt:
  case kEta:
  case kY:
  case kThetaDeg:
  case kDr:
  case kDz:
  case kTPCsignal:
  case kMultDiffRel:
  case kVzDiff:
    if (fVarType != kDouble) {
      AliWarning(Form("[INT CONSTRUCTOR] Cut '%s' is based on DOUBLE. Casting values to DOUBLE", GetName()));
      SetRange((Double_t)val1, (Double_t)val2);
      AliWarning(Form("[INT CONSTRUCTOR] Cut '%s' DOUBLE range = %f, %f", GetName(), fMinD, fMaxD));
    }
    break;
    // other cuts are not based on a value, so no problem
  default:
    break;
  }
}

//_________________________________________________________________________________________________
AliRsnCutStd::AliRsnCutStd
(const char *name, EType type, ULong_t val1, ULong_t val2, Bool_t useMC) :
    AliRsnCut(name, val1, val2),
    fType(type),
    fUseMC(useMC),
    fMass(0.0)
{
//
// Main constructor.
// Checks also that cut values are given in the correct type,
// in order to avoid that one passes, for example, a value which should be double
// but is interpreted as integer due to the overloading of constructors
//

  switch (fType) {
    // int
  case kMult:
  case kMultDiff:
  case kKink:
  case kKinkMother:
  case kAssignedPID:
  case kTruePID:
  case kRequiredPID:
  case kCharge:
    if (fVarType != kInt) {
      AliWarning(Form("[ULONG CONSTRUCTOR] Cut '%s' is based on INT. Casting values to INT", GetName()));
      SetRange((Int_t)val1, (Int_t)val2);
      AliWarning(Form("[ULONG CONSTRUCTOR] Cut '%s' INT range = %d, %d", GetName(), fMinI, fMaxI));
    }
    break;
    // ulong
  case kStatus:
    break;
    // double
  case kP:
  case kPt:
  case kEta:
  case kY:
  case kThetaDeg:
  case kDr:
  case kDz:
  case kTPCsignal:
  case kMultDiffRel:
  case kVzDiff:
    if (fVarType != kDouble) {
      AliWarning(Form("[ULONG CONSTRUCTOR] Cut '%s' is based on DOUBLE. Casting values to DOUBLE", GetName()));
      SetRange((Double_t)val1, (Double_t)val2);
      AliWarning(Form("[ULONG CONSTRUCTOR] Cut '%s' DOUBLE range = %f, %f", GetName(), fMinD, fMaxD));
    }
    break;
    // other cuts are not based on a value, so no problem
  default:
    break;
  }
}

//_________________________________________________________________________________________________
AliRsnCutStd::AliRsnCutStd
(const char *name, EType type, Double_t val1, Double_t val2, Bool_t useMC) :
    AliRsnCut(name, val1, val2),
    fType(type),
    fUseMC(useMC),
    fMass(0.0)
{
//
// Main constructor.
// Checks also that cut values are given in the correct type,
// in order to avoid that one passes, for example, a value which should be double
// but is interpreted as integer due to the overloading of constructors
//

  switch (fType) {
    // int
  case kMult:
  case kMultDiff:
  case kKink:
  case kKinkMother:
  case kAssignedPID:
  case kTruePID:
  case kRequiredPID:
  case kCharge:
    if (fVarType != kInt) {
      AliWarning(Form("[DOUBLE CONSTRUCTOR] Cut '%s' is based on INT. Casting values to INT", GetName()));
      SetRange((Int_t)val1, (Int_t)val2);
      AliWarning(Form("[DOUBLE CONSTRUCTOR] Cut '%s' INT range = %d, %d", GetName(), fMinI, fMaxI));
    }
    break;
    // ulong
  case kStatus:
    if (fVarType != kULong) {
      AliWarning(Form("[DOUBLE CONSTRUCTOR] Cut '%s' is based on ULONG. Casting values to ULONG", GetName()));
      SetRange((ULong_t)val1, (ULong_t)val2);
      AliWarning(Form("[DOUBLE CONSTRUCTOR] Cut '%s' ULONG range = %lu, %lu", GetName(), fMinU, fMaxU));
    }
    break;
    // double
  case kP:
  case kPt:
  case kEta:
  case kY:
  case kThetaDeg:
  case kDr:
  case kDz:
  case kTPCsignal:
  case kMultDiffRel:
  case kVzDiff:
    break;
    // other cuts are not based on a value, so no problem
  default:
    break;
  }
}

//_________________________________________________________________________________________________
Bool_t AliRsnCutStd::IsSelected(ETarget tgt, AliRsnDaughter * const track)
{
//
// Cut checker.
//

  // coherence check
  if (tgt != AliRsnCut::kParticle) {
    AliError(Form("[%s] Wrong target. Skipping cut", GetName()));
    return kTRUE;
  }


  // if MC info is required, try to retrieve the TParticle:
  // if it is not present, the cut is skipped
  // this avoids to raise a seg-fault when accessing the NULL TParticle object
  TParticle *part = 0x0;
  if (fUseMC) {
    part = track->GetParticle();
    if (!part) {
      AliError("Required MC info but it is not present. Cut Skipped.");
      return kTRUE;
    }
  }

  // loop on allowed possibilities
  switch (fType) {
  case kP:
    fCutValueD = fUseMC ? part->P() : track->P();
    return OkRange();
  case kPt:
    fCutValueD = fUseMC ? part->Pt() : track->Pt();
    return OkRange();
  case kThetaDeg:
    fCutValueD = track->ThetaDeg();
    return OkRange();
  case kEta:
    fCutValueD = fUseMC ? part->Eta() : track->Eta();
    return OkRange();
  case kDr:
    fCutValueD = track->Dr();
    return OkRange();
  case kDz:
    fCutValueD = track->Dz();
    return OkRange();
  case kStatus:
    fCutValueU = track->GetStatus();
    return OkValue();
  case kKink:
    fCutValueI = track->IsKink();
    return OkValue();
  case kKinkMother:
    fCutValueI = track->IsKinkMother();
    return OkValue();
  case kCharge:
    fCutValueI = (Int_t)track->Charge();
    return OkValue();
  case kTruePID:
    fCutValueI = (Int_t)track->PerfectPID();
    return OkValue();
  case kAssignedPID:
    fCutValueI = (Int_t)track->AssignedPID();
    return OkValue();
  case kRequiredPID:
    fCutValueI = (Int_t)track->RequiredPID();
    return OkValue();
  case kRealisticPID:
    fCutValueI = (Int_t)track->RealisticPID();
    return OkValue();
  case kPairIndex:
    fCutValueI = track->PairIndex();
    return OkValue();
  case kTruePIDMatch:
    return (track->PerfectPID() == track->RequiredPID());
  case kRealisticPIDMatch:
    return (track->RealisticPID() == track->RequiredPID());
  default:
    AliWarning(Form("Value %d is not included in available cuts for DAUGHTER. Cut skipped.", fType));
    return kTRUE;
  }
}

//_________________________________________________________________________________________________
Bool_t AliRsnCutStd::IsSelected(ETarget tgt, AliRsnPairParticle * const pair)
{
//
// Cut checker
//

  // coherence check
  if (tgt != AliRsnCut::kPair) {
    AliError(Form("[%s] Wrong target. Skipping cut", GetName()));
    return kTRUE;
  }

  // loop on allowed possibilities
  switch (fType) {
  case kP:
    fCutValueD = (fUseMC ? pair->GetPMC() : pair->GetP());
    return OkRange();
  case kPt:
    fCutValueD = (fUseMC ? pair->GetPtMC() : pair->GetPt());
    return OkRange();
  case kEta:
    fCutValueD = (fUseMC ? pair->GetEtaMC() : pair->GetEta());
    return OkRange();
  case kY:
    fCutValueD = (fUseMC ? pair->GetYMC(fMass) : pair->GetY(fMass));
    return OkRange();
  case kSameLabel:
    return pair->IsLabelEqual();
  case kTruePair:
    fCutValueI = TMath::Abs(pair->CommonMother());
    return OkValue();
  default:
    AliWarning(Form("Value %d is not included in available cuts for PAIR. Cut skipped.", fType));
    return kTRUE;
  }
}

//_________________________________________________________________________________________________
Bool_t AliRsnCutStd::IsSelected(ETarget tgt, AliRsnEvent * const event)
{
//
// Cut checker
//

  // coherence check
  if (tgt != AliRsnCut::kEvent) {
    AliError(Form("[%s] Wrong target. Skipping cut", GetName()));
    return kTRUE;
  }

  // loop on allowed possibilities
  switch (fType) {
  case kMult:
    fCutValueI = event->GetMultiplicity();
    return OkRange();
  default:
    AliWarning(Form("Value %d is not included in available cuts for EVENT. Cut skipped.", fType));
    return kTRUE;
  }
}

//_________________________________________________________________________________________________
Bool_t AliRsnCutStd::IsSelected(ETarget tgt, AliRsnEvent * const ev1, AliRsnEvent * const ev2)
{
//
// Cut checker
//

  // coherence check
  if (tgt != AliRsnCut::kMixEvent) {
    AliError(Form("[%s] Wrong target. Skipping cut", GetName()));
    return kTRUE;
  }

  // loop on allowed possibilities
  Double_t mult1, mult2;
  switch (fType) {
  case kMultDiff:
    fCutValueI = TMath::Abs(ev1->GetMultiplicity() - ev2->GetMultiplicity());
    return OkRange();
  case kMultDiffRel:
    mult1 = (Double_t)ev1->GetMultiplicity();
    mult2 = (Double_t)ev2->GetMultiplicity();
    if (mult1 == 0.0  && mult2 == 0.0) return kTRUE;
    fCutValueD = 100.0 * TMath::Abs(mult1 - mult2) / TMath::Max(mult1, mult2); // in %
    return OkRange();
  case kVzDiff:
    fCutValueD = TMath::Abs(ev1->GetVz() - ev2->GetVz());
    return OkRange();
  default:
    AliWarning(Form("Value %d is not included in available cuts for MIXEVENT. Cut skipped.", fType));
    return kTRUE;
  }
}
