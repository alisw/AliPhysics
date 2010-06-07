//
// Class AliRsnValue
//
// Definition of a single value which can be computed
// from any of the defined input objects implemented
// in the resonance package.
//

#include "AliRsnEvent.h"
#include "AliRsnDaughter.h"
#include "AliRsnPairParticle.h"
#include "AliRsnPairDef.h"

#include "AliRsnValue.h"

ClassImp(AliRsnValue)

//_____________________________________________________________________________
AliRsnValue::AliRsnValue(EAxisType type) : fType(type)
{
//
// Constructor
//
}

//_____________________________________________________________________________
const char* AliRsnValue::GetName() const
{
//
// Return the name of this object defined by the type
//

  switch (fType)
  {
    case kTrackPt:        return "PT";
    case kTrackEta:       return "ETA";
    case kTrack1P:        return "P1";
    case kTrack2P:        return "P2";
    case kTrack1Pt:       return "PT1";
    case kTrack2Pt:       return "PT2";
    case kPairInvMass:    return "IM";
    case kPairInvMassMC:  return "IMMC";
    case kPairInvMassRes: return "IMRES";
    case kPairPt:         return "PT";
    case kPairEta:        return "ETA";
    case kPairMt:         return "MT";
    case kPairY:          return "Y";
    case kEventMult:      return "MULT";
    default:              return "UNDEF";
  }
}

//_____________________________________________________________________________
AliRsnValue::EAxisObject AliRsnValue::GetAxisObject() const
{
//
// Tells what kind of object must be evaluated for this axis
//

  switch (fType)
  {
    // values coming from single track
    case kTrackPt:
    case kTrackEta:
      return kParticle;
    // values coming from pair
    case kTrack1P:
    case kTrack2P:
    case kTrack1Pt:
    case kTrack2Pt:
    case kPairInvMass:
    case kPairInvMassMC:
    case kPairInvMassRes:
    case kPairPt:
    case kPairEta:
    case kPairMt:
    case kPairY:
      return kPair;
    // values coming from event
    case kEventMult:
      return kEvent;
    default:
      return kNone;
  }
}

//_____________________________________________________________________________
Double_t AliRsnValue::Eval(TObject * const obj, const AliRsnPairDef *pairDef) const
{
//
// Evaluation of the required value.
// Checks that the passed object is of the right type
// and if this check is successful, returns the required value.
// The output of the function tells if it was successful,
// and the values must be taken with GetValue().
//

  // dynamic casting
  AliRsnDaughter     *track = dynamic_cast<AliRsnDaughter*>(obj);
  AliRsnPairParticle *pair  = dynamic_cast<AliRsnPairParticle*>(obj);
  AliRsnEvent        *event = dynamic_cast<AliRsnEvent*>(obj);

  // check that the object type and required input match
  EAxisObject requiredInput = GetAxisObject();
  Double_t    mass          = pairDef->GetMotherMass();
  if (requiredInput == kParticle && !track)
    return 0.0;
  else if (requiredInput == kPair)
  {
    if (!pair)
      return 0.0;
    else
    {
      // if we are filling a pair, we need the pair def
      if (pairDef)
        mass = pairDef->GetMotherMass();
      else
      {
        AliError("Cannot compute a 'pair' value if I don't have a vaild PairDef");
        return 0.0;
      }
    }
  }
  else if (requiredInput == kEvent && !event)
  {
    return 0.0;
  }
  else
  {
    AliError(Form("Failed computation: expected type %d, passed a '%s'", (Int_t)requiredInput, obj->ClassName()));
    return 0.0;
  }

  switch (fType)
  {
    case kTrackPt:
      return track->Pt();
    case kTrackEta:
      return track->Eta();
    case kTrack1P:
      return pair->GetDaughter(0)->P();
    case kTrack2P:
      return pair->GetDaughter(1)->P();
    case kTrack1Pt:
      return pair->GetDaughter(0)->Pt();
    case kTrack2Pt:
      return pair->GetDaughter(1)->Pt();
    case kPairInvMass:
      return pair->GetInvMass(pairDef->GetMass(0), pairDef->GetMass(1));
    case kPairInvMassMC:
      return pair->GetInvMassMC(pairDef->GetMass(0), pairDef->GetMass(1));
    case kPairInvMassRes:
      {
        Double_t value;
        value  = pair->GetInvMass  (pairDef->GetMass(0), pairDef->GetMass(1));
        value -= pair->GetInvMassMC(pairDef->GetMass(0), pairDef->GetMass(1));
        value /= pair->GetInvMassMC(pairDef->GetMass(0), pairDef->GetMass(1));
        return value;
      }
    case kPairPt:
      return pair->GetPt();
    case kPairEta:
      return pair->GetEta();
    case kPairMt:
      if (TMath::Abs(mass) < 1E-5) AliWarning(Form("Suspicious mass value specified: %f", mass));
      return (TMath::Sqrt(pair->GetPt()*pair->GetPt() + mass*mass) - mass);
    case kPairY:
      if (TMath::Abs(mass) < 1E-5) AliWarning(Form("Suspicious mass value specified: %f", mass));
      return pair->GetY(mass);
    case kEventMult:
      return (Double_t)event->GetMultiplicity();
    default:
      AliWarning("Invalid value type");
      return 0.0;
  }
}
