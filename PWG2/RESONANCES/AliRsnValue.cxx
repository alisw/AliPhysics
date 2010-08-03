
//
// Class AliRsnValue
//
// Definition of a single value which can be computed
// from any of the defined input objects implemented
// in the resonance package.
//

#include "AliRsnEvent.h"
#include "AliRsnDaughter.h"
#include "AliRsnMother.h"
#include "AliRsnPairDef.h"

#include "AliRsnValue.h"

ClassImp(AliRsnValue)

//_____________________________________________________________________________
AliRsnValue::AliRsnValue() :
  TNamed(),
  fType(kValueTypes),
  fNBins(0),
  fMin(0.0),
  fMax(0.0),
  fValue(0.0)
{
//
// Main constructor (version 1)
// This can also be created without any argument.
//
}

//_____________________________________________________________________________
AliRsnValue::AliRsnValue
(const char *name, EValueType type, Int_t nbins, Double_t min, Double_t max) :
  TNamed(name, ""),
  fType(type),
  fNBins(0),
  fMin(0.0),
  fMax(0.0),
  fValue(0.0)
{
//
// Main constructor (version 1)
// This can also be created without any argument.
//

  SetBins(nbins, min, max);
}

//_____________________________________________________________________________
AliRsnValue::AliRsnValue
(const char *name, EValueType type, Double_t min, Double_t max, Double_t step) :
  TNamed(name, ""),
  fType(type),
  fNBins(0),
  fMin(0.0),
  fMax(0.0),
  fValue(0.0)
{
//
// Main constructor (version 2)
//

  SetBins(min, max, step);
}

/*
//_____________________________________________________________________________
const char* AliRsnValue::GetName() const
{
//
// Return the name of this object defined by the type
//

  switch (fType)
  {
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
*/

//_____________________________________________________________________________
TArrayD AliRsnValue::GetArray() const
{
//
// Creates an array with all bin edges
//

  TArrayD out(fNBins + 1);

  Int_t i;
  Double_t step = (fMax - fMin) / (Double_t)fNBins;

  for (i = 0; i <= fNBins; i++) out[i] = fMin + step * (Double_t)i;

  return out;
}

//_____________________________________________________________________________
void AliRsnValue::SetBins(Int_t n, Double_t min, Double_t max)
{
//
// Set binning for histogram.
//

  fNBins = n;

  if (min < max) 
  {
    fMin = min;
    fMax = max;
  } 
  else 
  {
    fMin = max;
    fMax = min;
  }
}

//_____________________________________________________________________________
void AliRsnValue::SetBins(Double_t min, Double_t max, Double_t step)
{
//
// Binning for histogram.
//

  if (min < max) 
  {
    fMin = min;
    fMax = max;
  } 
  else 
  {
    fMin = max;
    fMax = min;
  }

  fNBins = (Int_t)((fMax - fMin) / (step)) + 1;
}

//_____________________________________________________________________________
Bool_t AliRsnValue::Eval(AliRsnMother *mother, const AliRsnPairDef *pairDef, AliRsnEvent *event)
{
//
// Evaluation of the required value.
// Checks that the passed object is of the right type
// and if this check is successful, returns the required value.
// The output of the function tells if it was successful,
// and the values must be taken with GetValue().
//

  // avoid segfaults
  if (!mother) return kFALSE;
  if (!pairDef) return kFALSE;

  Double_t mass = pairDef->GetMotherMass();

  switch (fType)
  {
    case kTrack1P:
      fValue = mother->GetDaughter(0)->P().Mag();
      break;
    case kTrack2P:
      fValue = mother->GetDaughter(1)->P().Mag();
      break;
    case kTrack1Pt:
      fValue = mother->GetDaughter(0)->P().Perp();
      break;
    case kTrack2Pt:
      fValue = mother->GetDaughter(1)->P().Perp();
      break;
    case kPairInvMass:
      fValue = mother->Sum().M();
      break;
    case kPairInvMassMC:
      fValue = mother->SumMC().M();
      break;
    case kPairInvMassRes:
      fValue = (mother->SumMC().M() - mother->Sum().M()) / mother->SumMC().M();
      break;
    case kPairPt:
      fValue = mother->Sum().Perp();
      break;
    case kPairEta:
      fValue = mother->Sum().Eta();
      break;
    case kPairMt:
      if (TMath::Abs(mass) < 1E-5) AliWarning(Form("Suspicious mass value specified: %f", mass));
      fValue = (TMath::Sqrt(mother->Sum().Perp2() + mass*mass) - mass);
      break;
    case kPairY:
      if (TMath::Abs(mass) < 1E-5) AliWarning(Form("Suspicious mass value specified: %f", mass));
      mother->SetDefaultMass(mass);
      fValue = mother->Ref().Rapidity();
      break;
    case kPairCosThetaStar1:
      fValue = TMath::Cos(mother->ThetaStar(kTRUE, kFALSE));
      break;
    case kPairCosThetaStar2:
      fValue = TMath::Cos(mother->ThetaStar(kFALSE, kFALSE));
      break;
    case kPairCosThetaStarMC1:
      fValue = TMath::Cos(mother->ThetaStar(kTRUE, kTRUE));
      break;
    case kPairCosThetaStarMC2:
      fValue = TMath::Cos(mother->ThetaStar(kFALSE, kTRUE));
      break;
    case kEventMult:
      if (!event) fValue = 0.0;
      fValue = (Double_t)event->GetMultiplicity();
      break;
    default:
      AliWarning("Invalid value type");
      return kFALSE;
  }
  
  return kTRUE;
}
