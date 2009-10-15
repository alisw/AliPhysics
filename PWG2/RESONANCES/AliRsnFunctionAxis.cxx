//
// Class AliRsnFunctionAxis
//
// Definition for a histogram type.
// Since one could do an analysis which is not an invariant mass
// the histogram definition should be more flexible, and it is stored
// separately in a new class.
// This class considers the possibility of a 1D or 2D histograms
// with its related binning, and can create a new histo from his definitions
//

#include <TArrayD.h>

#include "AliRsnEvent.h"
#include "AliRsnPairParticle.h"
#include "AliRsnPairDef.h"
#include "AliRsnFunctionAxis.h"

ClassImp(AliRsnFunctionAxis)

//_____________________________________________________________________________
AliRsnFunctionAxis::AliRsnFunctionAxis() :
    fType(kAxisTypes),
    fNBins(0),
    fMin(0.0),
    fMax(0.0),
    fMass(0.0)
{
//
// Default constructor
//
}

//_____________________________________________________________________________
AliRsnFunctionAxis::AliRsnFunctionAxis
(EAxisType type, Int_t nbins, Double_t min, Double_t max) :
    fType(type),
    fNBins(0),
    fMin(0.0),
    fMax(0.0),
    fMass(0.0)
{
//
// Main constructor (version 1)
//

  SetBins(nbins, min, max);
}

//_____________________________________________________________________________
AliRsnFunctionAxis::AliRsnFunctionAxis
(EAxisType type, Double_t min, Double_t max, Double_t step) :
    fType(type),
    fNBins(0),
    fMin(0.0),
    fMax(0.0),
    fMass(0.0)
{
//
// Main constructor (version 2)
//

  SetBins(min, max, step);
}

//_____________________________________________________________________________
const char* AliRsnFunctionAxis::GetName() const
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
TArrayD AliRsnFunctionAxis::GetArray() const
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
AliRsnFunctionAxis::EAxisObject AliRsnFunctionAxis::GetAxisObject() const
{
//
// Tells what kind of object must be evaluated for this axis
//

  switch (fType)
  {
    case kTrackPt:
    case kTrackEta:
      return kParticle;
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
    case kEventMult:
      return kEvent;
    default:
      return kNone;
  }
}

//_____________________________________________________________________________
void AliRsnFunctionAxis::SetBins(Int_t n, Double_t min, Double_t max)
{
//
// Set binning for histogram.
//

  fNBins = n;

  if (min < max) {
    fMin = min;
    fMax = max;
  } else {
    fMin = max;
    fMax = min;
  }
}

//_____________________________________________________________________________
void AliRsnFunctionAxis::SetBins(Double_t min, Double_t max, Double_t step)
{
//
// Binning for histogram.
//

  if (min < max) {
    fMin = min;
    fMax = max;
  } else {
    fMin = max;
    fMax = min;
  }

  fNBins = (Int_t)((fMax - fMin) / (step)) + 1;
}

//_____________________________________________________________________________
Double_t AliRsnFunctionAxis::Eval(AliRsnDaughter* daughter) const
{
//
// EValuation method for single tracks
// (currently disabled)
//

  switch (fType)
  {
    case kTrackPt:
      return daughter->Pt();
    case kTrackEta:
      return daughter->Eta();
    default:
      AliWarning("This axis type cannot be applied to particles");
      return -999.0;
  }

  return 0.0;
}

//_____________________________________________________________________________
Double_t AliRsnFunctionAxis::Eval(AliRsnPairParticle * const pair, AliRsnPairDef *const pairDef) const
{
//
// EValuation method for pairs.
// Requires also the pair definitions, in order to retrieve mass.
//

  switch (fType)
  {
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
    case kPairInvMassRes: {
        Double_t value;
        value  = pair->GetInvMass(pairDef->GetMass(0), pairDef->GetMass(1));
        value -= pair->GetInvMassMC(pairDef->GetMass(0), pairDef->GetMass(1));
        value /= pair->GetInvMassMC(pairDef->GetMass(0), pairDef->GetMass(1));
        return value;
      }
    case kPairPt:
      return pair->GetPt();
    case kPairEta:
      return pair->GetEta();
    case kPairMt:
      if (TMath::Abs(fMass) < 1E-5) AliWarning(Form("Suspicious mass value specified: %f", fMass));
      return TMath::Sqrt(pair->GetPt()*pair->GetPt() + fMass*fMass);
    case kPairY:
      if (TMath::Abs(fMass) < 1E-5) AliWarning(Form("Suspicious mass value specified: %f", fMass));
      return pair->GetY(fMass);
    default:
      AliWarning("This axis type cannot be applied to pairs");
      return -999.0;
  }
}

//_____________________________________________________________________________
Double_t AliRsnFunctionAxis::Eval(AliRsnEvent *const event) const
{
//
// EValuation method for events
//

  switch (fType)
  {
    case kEventMult:
      return (Double_t)event->GetMultiplicity();
    default:
      AliWarning("This axis type cannot be applied to events");
      return 0.0;
  }
}
