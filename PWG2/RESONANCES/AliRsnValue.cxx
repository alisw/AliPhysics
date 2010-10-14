
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
  fValue(0.0),
  fType(kValueTypes),
  fArray(0)
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
  fValue(0.0),
  fType(type),
  fArray(0)
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
  fValue(0.0),
  fType(type),
  fArray(0)
{
//
// Main constructor (version 2)
//

  SetBins(min, max, step);
}

//_____________________________________________________________________________
AliRsnValue::AliRsnValue
(const char *name, EValueType type, Int_t nbins, Double_t *array) :
  TNamed(name, ""),
  fValue(0.0),
  fType(type),
  fArray(0)
{
//
// Main constructor (version 2)
//

  SetBins(nbins, array);
}

//_____________________________________________________________________________
void AliRsnValue::SetBins(Int_t nbins, Double_t min, Double_t max)
{
//
// Set binning for the axis in equally spaced bins
// where the number of bins, minimum and maximum are given.
//

  fArray.Set(nbins + 1);
  
  Double_t mymax = TMath::Max(min, max);
  Double_t mymin = TMath::Min(min, max);
  
  Int_t    k = 0;
  Double_t binSize = (mymax - mymin) / ((Double_t)nbins);
  
  fArray[0] = mymin;
  for (k = 1; k <= nbins; k++) fArray[k] = fArray[k-1] + binSize;
  for (k = 0; k < fArray.GetSize() - 1; k++) AliDebug(AliLog::kDebug + 3, Form("Bin #%d: %f - %f", k, fArray[k], fArray[k+1]));
}

//_____________________________________________________________________________
void AliRsnValue::SetBins(Double_t min, Double_t max, Double_t step)
{
//
// Set binning for the axis in equally spaced bins
// where the bin size, minimum and maximum are given.
//

  Double_t dblNbins = TMath::Abs(max - min) / step;
  Int_t    intNbins = ((Int_t)dblNbins) + 1;
  
  SetBins(intNbins, min, max);
}

//_____________________________________________________________________________
void AliRsnValue::SetBins(Int_t nbins, Double_t *array)
{
//
// Set binning for the axis in unequally spaced bins
// using the same way it is done in TAxis
//

  fArray.Adopt(nbins, array);
  for (Int_t k = 0; k < fArray.GetSize() - 1; k++) AliDebug(AliLog::kDebug + 3, Form("Bin #%d: %f - %f", k, fArray[k], fArray[k+1]));
}

//_____________________________________________________________________________
Bool_t AliRsnValue::Eval(AliRsnMother * const mother, AliRsnPairDef * const pairDef, AliRsnEvent * const event)
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
    case kTrack1Px:
      fValue = mother->GetDaughter(0)->P().X();
      break;
    case kTrack1Py:
      fValue = mother->GetDaughter(0)->P().Y();
      break;
    case kTrack1Pz:
      fValue = mother->GetDaughter(0)->P().Z();
      break;
    case kTrack2Px:
      fValue = mother->GetDaughter(1)->P().X();
      break;
    case kTrack2Py:
      fValue = mother->GetDaughter(1)->P().Y();
      break;
    case kTrack2Pz:
      fValue = mother->GetDaughter(1)->P().Z();
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
    case kPairPhi:
      fValue = mother->Sum().Phi();
      break;
    case kPairPhiMC:
      fValue = mother->SumMC().Phi();
      break;
    case kPairPtRatio:
      fValue  = TMath::Abs(mother->GetDaughter(0)->P().Perp() - mother->GetDaughter(1)->P().Perp());
      fValue /= TMath::Abs(mother->GetDaughter(0)->P().Perp() + mother->GetDaughter(1)->P().Perp());
      break;
    case kPairDipAngle:
      fValue = mother->GetDaughter(0)->P().Angle(mother->GetDaughter(1)->P().Vect());
      fValue = TMath::Abs(TMath::ACos(fValue));
      break;
    case kPairCosThetaStar:
      fValue = mother->CosThetaStar();
      break;
    case kPairCosThetaStar1:
      //fValue = TMath::Cos(mother->ThetaStar(kTRUE, kFALSE));
      break;
    case kPairCosThetaStar2:
      //fValue = TMath::Cos(mother->ThetaStar(kFALSE, kFALSE));
      break;
    case kPairCosThetaStarMC1:
      //fValue = TMath::Cos(mother->ThetaStar(kTRUE, kTRUE));
      break;
    case kPairCosThetaStarMC2:
      //fValue = TMath::Cos(mother->ThetaStar(kFALSE, kTRUE));
      break;
    case kAngleToLeading:
      {
    	  int ID1 = (mother->GetDaughter(0))->GetID();
    	  int ID2 = (mother->GetDaughter(1))->GetID();
    	  int leadingID = event->SelectLeadingParticle(0);
    	  if(leadingID == ID1 || leadingID == ID2) return kFALSE;
    	  AliRsnDaughter  leadingPart = event->GetDaughter(leadingID);
    	  AliVParticle *ref = leadingPart.GetRef();

    	  fValue = ref->Phi() - mother->Sum().Phi();
    	  //return angle w.r.t. leading particle in the range -pi/2, 3/2pi
    	  while(fValue >= TMath::Pi()) fValue -= 2*TMath::Pi();
    	  while(fValue < -0.5*TMath::Pi()) fValue += 2*TMath::Pi();
    	  //Printf("%g", fValue);
      }
      break;
    case kEventMult:
      if (!event) 
      {
        fValue = 0.0;
        return kFALSE;
      }
      else fValue = (Double_t)event->GetMultiplicity();
      break;
    case kLeadingPt:
      if (!event) 
      {
        fValue = 0.0;
        return kFALSE;
      }
      else
      {
    	  int leadingID = event->SelectLeadingParticle(0);
    	  if(leadingID >= 0) {
    		  AliRsnDaughter leadingPart = event->GetDaughter(leadingID);
    		  AliVParticle *ref = leadingPart.GetRef();
    		  fValue = ref->Pt();
    	  }
    	  else fValue = 0;
      }
      break;
    case kQInv:
      {
        TLorentzVector diff = mother->GetDaughter(0)->P() - mother->GetDaughter(1)->P();
        fValue = diff.M();
      }
      break;
    default:
      AliWarning("Invalid value type");
      return kFALSE;
  }
  
  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliRsnValue::Eval(AliRsnDaughter * const daughter, AliRsnEvent * const event)
{
//
// Evaluation of the required value.
// Checks that the passed object is of the right type
// and if this check is successful, returns the required value.
// The output of the function tells if it was successful,
// and the values must be taken with GetValue().
//

  // avoid segfaults
  if (!daughter) return kFALSE;

  switch (fType)
  {
    case kEventMult:
      if (!event) 
      {
        fValue = 0.0;
        return kFALSE;
      }
      else fValue = (Double_t)event->GetMultiplicity();
      break;
    case kLeadingPt:
      if (!event) 
      {
        fValue = 0.0;
        return kFALSE;
      }
      else
      {
    	  int leadingID = event->SelectLeadingParticle(0);
    	  if(leadingID >= 0) {
    		  AliRsnDaughter leadingPart = event->GetDaughter(leadingID);
    		  AliVParticle *ref = leadingPart.GetRef();
    		  fValue = ref->Pt();
    	  }
    	  else fValue = 0;
      }
      break;
    default:
      AliWarning("Invalid value type");
      return kFALSE;
  }
  
  return kTRUE;
}

//_____________________________________________________________________________
void AliRsnValue::Print(Option_t *) const
{
//
// Print all bins
//

  Int_t   i, n = fArray.GetSize();
  TString msg("Array values: ");
  
  for (i = 0; i < n; i++) msg += Form("%f, ", fArray[i]);
  
  AliInfo(Form("Axis name: %s", GetName()));
  AliInfo(msg.Data());
}
