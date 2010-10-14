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

#include <TMath.h>
#include <TLorentzVector.h>

#include "AliStack.h"
#include "AliMCEvent.h"
#include "AliRsnDaughter.h"
#include "AliRsnMother.h"
#include "AliRsnEvent.h"

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
(const char *name, ETarget target, EType type, Int_t val1, Int_t val2, Bool_t useMC) :
  AliRsnCut(name, target, val1, val2),
  fType(type),
  fUseMC(useMC),
  fMass(0.0)
{
//
// Main constructor.
// Checks also that cut values are given in the correct type,
// in order to avoid that one passes, for example, a value which should be double
// but is interpreted as integer due to the overloading of constructors.
//

  switch (fType) 
  {
    // int
    case kMult:
    case kCharge:
      break;
    // double
    case kP:
    case kPt:
    case kPtLeading:
    case kEta:
    case kY:
    case kDipAngle:
    case kThetaDeg:
      if (fVarType != kDouble) 
      {
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
(const char *name, ETarget target, EType type, Double_t val1, Double_t val2, Bool_t useMC) :
  AliRsnCut(name, target, val1, val2),
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

  switch (fType) 
  {
    // int
    case kMult:
    case kCharge:
      if (fVarType != kInt) 
      {
        AliWarning(Form("[DOUBLE CONSTRUCTOR] Cut '%s' is based on INT. Casting values to INT", GetName()));
        SetRange((Int_t)val1, (Int_t)val2);
        AliWarning(Form("[DOUBLE CONSTRUCTOR] Cut '%s' INT range = %d, %d", GetName(), fMinI, fMaxI));
      }
      break;
    // double
    case kP:
    case kPt:
    case kPtLeading:
    case kEta:
    case kY:
    case kDipAngle:
    case kThetaDeg:
      break;
    // other cuts are not based on a value, so no problem
    default:
      break;
  }
}

//_________________________________________________________________________________________________
AliRsnCut::EVarType AliRsnCutStd::CheckType()
{
//
// Returns the variable type expected for the selected cut type
//

  switch (fType) 
  {
    // integer cuts
    case kMult:
    case kCharge:
      return kInt;
    // double couts
    case kP:
    case kPt:
    case kPtLeading:
    case kEta:
    case kY:
    case kDipAngle:
    case kThetaDeg:
      return kDouble;
    // other cuts are not based on a value, so no problem
    default:
      return kNoVar;
  }
}

//_________________________________________________________________________________________________
Bool_t AliRsnCutStd::IsSelected(TObject *obj1, TObject *obj2)
{
  // coherence check using the mother class method
  if (!TargetOK(obj1, obj2)) 
  {
    AliError("Wrong target. Skipping cut");
    return kTRUE;
  }
  
  // if coherence check is OK, only one of the following
  // dynamic casts will be successful, and it will trigger
  // the correct internal method to check the cut
  AliRsnDaughter *objD = dynamic_cast<AliRsnDaughter*>(obj1);
  AliRsnMother   *objM = dynamic_cast<AliRsnMother*>(obj1);
  AliRsnEvent    *objE = dynamic_cast<AliRsnEvent*>(obj1);  
  if (objD) return IsDaughterSelected(objD);
  else if (objM) return IsMotherSelected(objM);
  else if (objE) return IsEventSelected(objE);
  else return kTRUE;
}

//_________________________________________________________________________________________________
Bool_t AliRsnCutStd::IsDaughterSelected(AliRsnDaughter *daughter)
{
//
// Cut checker.
//

  // get the correct reference object for kinematic cuts
  // and check that this reference is not NULL, to avoid segfaults
  AliVParticle *ref = fUseMC ? daughter->GetRefMC() : daughter->GetRef();

  // loop on allowed possibilities
  switch (fType) 
  {
    case kP:
      fCutValueD = ref->P();
      return OkRange();
    case kPt:
      fCutValueD = ref->Pt();
      return OkRange();
    case kThetaDeg:
      fCutValueD = ref->Theta() * TMath::RadToDeg();
      return OkRange();
    case kEta:
      fCutValueD = ref->Eta();
      return OkRange();
    case kCharge:
      fCutValueI = (Int_t)ref->Charge();
      return OkValue();
    case kPhysPrimary:
      if (!fEvent->GetRefMC()) return kFALSE;
      else
      {
        return fEvent->GetRefMC()->Stack()->IsPhysicalPrimary(TMath::Abs(((AliVTrack*)ref)->GetLabel()));
      }
    default:
      AliWarning(Form("Value %d is not included in available cuts for DAUGHTER. Cut skipped.", fType));
      return kTRUE;
  }
}

//_________________________________________________________________________________________________
Bool_t AliRsnCutStd::IsMotherSelected(AliRsnMother * const mother)
{
//
// Cut checker
//
  
  // use MC flag if required
  TLorentzVector &sum = fUseMC ? mother->SumMC() : mother->Sum();
  TLorentzVector &ref = fUseMC ? mother->RefMC() : mother->Ref();
  
  // loop on allowed possibilities
  switch (fType) 
  {
    case kP:
      fCutValueD = sum.P();
      return OkRange();
    case kPt:
      fCutValueD = sum.Perp();
      return OkRange();
    case kEta:
      fCutValueD = sum.Eta();
      return OkRange();
    case kY:
      fCutValueD = ref.Rapidity();
      return OkRange();
    case kDipAngle:
      fCutValueD = mother->GetDaughter(0)->P().Angle(mother->GetDaughter(1)->P().Vect());
      fCutValueD = TMath::Abs(TMath::ACos(fCutValueD));
      return OkRangeD();
    case kSameLabel:
      return mother->IsLabelEqual();
    default:
      AliWarning(Form("Value %d is not included in available cuts for PAIR. Cut skipped.", fType));
      return kTRUE;
  }
}

//_________________________________________________________________________________________________
Bool_t AliRsnCutStd::IsEventSelected(AliRsnEvent * const event)
{
//
// Cut checker
//

  // loop on allowed possibilities
  switch (fType) 
  {
    case kMult:
      fCutValueI = event->GetMultiplicity();
      return OkRange();
    case kPtLeading:
    {
      int leadingID = event->SelectLeadingParticle(0);
      if(leadingID >= 0) {
    	  AliRsnDaughter leadingPart = event->GetDaughter(leadingID);
    	  AliVParticle *ref = fUseMC ? leadingPart.GetRefMC() : leadingPart.GetRef();
    	  fCutValueD = ref->Pt();
      }
      else fCutValueD = 0;
      return OkRange();
    }
    default:
      AliWarning(Form("Value %d is not included in available cuts for EVENT. Cut skipped.", fType));
      return kTRUE;
  }
}
