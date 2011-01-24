
//
// Class AliRsnValue
//
// Definition of a single value which can be computed
// from any of the defined input objects implemented
// in the resonance package.
//

#include <Riostream.h>
#include "AliESDtrackCuts.h"
#include "AliRsnEvent.h"
#include "AliRsnDaughter.h"
#include "AliRsnMother.h"
#include "AliRsnPairDef.h"

#include "AliRsnValue.h"

ClassImp(AliRsnValue)

//_____________________________________________________________________________
AliRsnValue::AliRsnValue() :
  AliRsnTarget(),
  fComputedValue(0),
  fValueType(kValueTypes),
  fBinArray(0),
  fSupportObject(0x0)
{
//
// Default constructor without arguments.
// Initialize data members to meaningless values.
// This method is provided for ROOT streaming, 
// but should never be used directly by a user.
//

  AssignTarget();
}

//_____________________________________________________________________________
AliRsnValue::AliRsnValue
(const char *name, EValueType type, Int_t nbins, Double_t min, Double_t max) :
  AliRsnTarget(name, AliRsnTarget::kTargetTypes),
  fComputedValue(0.0),
  fValueType(type),
  fBinArray(0),
  fSupportObject(0x0)
{
//
// Main constructor (version 1).
// This constructor defines in meaningful way all data members,
// and defined a fixed binnings, subdividing the specified interval
// into that many bins as specified in the integer argument.
// ---
// This method is also the entry point for all instances
// of this class which don't need to do binning (e.g.: TNtuple inputs),
// since arguments 3 to 5 have default values which don't create any
// binning array, in order not to allocate memory when this is useless.
//

  AssignTarget();
  SetBins(nbins, min, max);
}

//_____________________________________________________________________________
AliRsnValue::AliRsnValue
(const char *name, EValueType type, Double_t min, Double_t max, Double_t step) :
  AliRsnTarget(name, AliRsnTarget::kTargetTypes),
  fComputedValue(0.0),
  fValueType(type),
  fBinArray(0),
  fSupportObject(0x0)
{
//
// Main constructor (version 2).
// This constructor defines in meaningful way all data members
// and creates enough equal bins of the specified size to cover
// the required interval.
//

  AssignTarget();
  SetBins(min, max, step);
}

//_____________________________________________________________________________
AliRsnValue::AliRsnValue
(const char *name, EValueType type, Int_t nbins, Double_t *array) :
  AliRsnTarget(name, AliRsnTarget::kTargetTypes),
  fComputedValue(0.0),
  fValueType(type),
  fBinArray(0),
  fSupportObject(0x0)
{
//
// Main constructor (version 3).
// This constructor defines in meaningful way all data members
// and creates a set of variable bins delimited by the passed array.
//

  AssignTarget();
  SetBins(nbins, array);
}

//_____________________________________________________________________________
AliRsnValue::AliRsnValue(const AliRsnValue& copy) : 
  AliRsnTarget(copy),
  fComputedValue(copy.fComputedValue),
  fValueType(copy.fValueType),
  fBinArray(copy.fBinArray),
  fSupportObject(copy.fSupportObject)
{
//
// Copy constructor.
// Duplicates the binning array and copies all settings.
// Calls also the function that assigns properly the 
// expected target, depending on the computation type.
//

  AssignTarget();
}

//_____________________________________________________________________________
AliRsnValue& AliRsnValue::operator=(const AliRsnValue& copy)
{
//
// Assignment operator.
// Works like copy constructor.
//

  AliRsnTarget::operator=(copy);
  
  fComputedValue = copy.fComputedValue;
  fBinArray = copy.fBinArray;
  fSupportObject = copy.fSupportObject;
  
  AssignTarget();
  
  return (*this);
}

//_____________________________________________________________________________
void AliRsnValue::SetBins(Int_t nbins, Double_t min, Double_t max)
{
//
// Set binning for the axis in equally spaced bins
// where the number of bins, minimum and maximum are given.
//

  if (!nbins)
  {
    fBinArray.Set(0);
    return;
  }

  fBinArray.Set(nbins + 1);
  
  Double_t mymax = TMath::Max(min, max);
  Double_t mymin = TMath::Min(min, max);
  
  Int_t    k = 0;
  Double_t binSize = (mymax - mymin) / ((Double_t)nbins);
  
  fBinArray[0] = mymin;
  for (k = 1; k <= nbins; k++) fBinArray[k] = fBinArray[k-1] + binSize;
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

  if (!nbins)
  {
    fBinArray.Set(0);
    return;
  }
  
  fBinArray.Adopt(nbins, array);
}

//_____________________________________________________________________________
const char* AliRsnValue::GetValueTypeName() const
{
//
// This method returns a string to give a name to each possible
// computation value.
//
  
  switch (fValueType)
  {
    case kTrackP:             return "SingleTrackPtot";
    case kTrackPt:            return "SingleTrackPt";
    case kTrackEta:           return "SingleTrackEta";
    case kPairP1:             return "PairPtotDaughter1";
    case kPairP2:             return "PairPtotDaughter2";
    case kPairP1t:            return "PairPtDaughter1";
    case kPairP2t:            return "PairPtDaughter2";
    case kPairP1z:            return "PairPzDaughter1";
    case kPairP2z:            return "PairPzDaughter2";
    case kPairInvMass:        return "PairInvMass";
    case kPairInvMassMC:      return "PairInvMassMC";
    case kPairInvMassRes:     return "PairInvMassResolution";
    case kPairPt:             return "PairPt";
    case kPairPz:             return "PairPz";
    case kPairEta:            return "PairEta";
    case kPairMt:             return "PairMt";
    case kPairY:              return "PairY";
    case kPairPhi:            return "PairPhi";
    case kPairPhiMC:          return "PairPhiMC";
    case kPairPtRatio:        return "PairPtRatio";
    case kPairDipAngle:       return "PairDipAngle";
    case kPairCosThetaStar:   return "PairCosThetaStar";
    case kPairQInv:           return "PairQInv";
    case kPairAngleToLeading: return "PairAngleToLeading";
    case kEventLeadingPt:     return "EventLeadingPt";
    case kEventMult:          return "EventMult";
    case kEventMultESDCuts:   return "EventMultESDCuts";
    case kEventVz:            return "EventVz";
    default:                  return "Undefined";
  }
}

//_____________________________________________________________________________
void AliRsnValue::AssignTarget()
{
//
// This method assigns the target to be expected by this object
// in the computation, depending on its type chosen in the enum.
//
  
  switch (fValueType)
  {
    // track related values
    case kTrackP:
    case kTrackPt:
    case kTrackEta:
      SetTargetType(AliRsnTarget::kDaughter); // end of track-related values
      break;
    // pair related values
    case kPairP1:
    case kPairP2:
    case kPairP1t:
    case kPairP2t:
    case kPairP1z:
    case kPairP2z:
    case kPairInvMass:
    case kPairInvMassMC:
    case kPairInvMassRes:
    case kPairPt:
    case kPairPz:
    case kPairEta:
    case kPairMt:
    case kPairY:
    case kPairPhi:
    case kPairPhiMC:
    case kPairPtRatio:
    case kPairDipAngle:
    case kPairCosThetaStar:
    case kPairQInv:
    case kPairAngleToLeading:
      SetTargetType(AliRsnTarget::kMother); // end of pair-related values
      break;
    // event related values
    case kEventLeadingPt:
    case kEventMult:
    case kEventMultESDCuts:
    case kEventVz:
      SetTargetType(AliRsnTarget::kEvent); // end of event-related values
      break;
    // undefined value
    default:
      SetTargetType(AliRsnTarget::kTargetTypes); // undefined targets
  }
}

//_____________________________________________________________________________
Bool_t AliRsnValue::Eval(TObject *object, Bool_t useMC)
{
//
// Evaluation of the required value.
// Checks that the passed object is of the right type
// and if this check is successful, computes the required value.
// The output of the function tells if computing was successful,
// and the values must be taken with GetValue().
//

  // cast the input to the allowed types
  AliRsnDaughter *daughter = dynamic_cast<AliRsnDaughter*>(object);
  AliRsnMother   *mother   = dynamic_cast<AliRsnMother*>(object);
  
  // common variables
  TLorentzVector pRec;   // 4-momentum for single track or pair sum (reco)
  TLorentzVector pSim;   // 4-momentum for single track or pair sum (MC)
  TLorentzVector pRec0;  // 4-momentum of first daughter (reco)
  TLorentzVector pSim0;  // 4-momentum of first daughter (MC)
  TLorentzVector pRec1;  // 4-momentum of second daughter (reco)
  TLorentzVector pSim1;  // 4-momentum of second daughter (MC)
  
  // check that the input object is the correct class type
  switch (fTargetType)
  {
    case AliRsnTarget::kDaughter:
      if (daughter)
      {
        pRec = daughter->Psim();
        pSim = daughter->Prec();
      }
      else
      {
        AliError(Form("[%s] expected: AliRsnDaughter, passed: [%s]", GetName(), object->ClassName()));
        return kFALSE;
      }
      break;
    case AliRsnTarget::kMother:
      if (mother)
      {
        pRec  = mother->Sum();
        pSim  = mother->SumMC();
        pRec0 = mother->GetDaughter(0)->Prec();
        pRec1 = mother->GetDaughter(1)->Prec();
        pSim0 = mother->GetDaughter(0)->Psim();
        pSim1 = mother->GetDaughter(1)->Psim();
      }
      else
      {
        AliError(Form("[%s] expected: AliRsnMother, passed: [%s]", GetName(), object->ClassName()));
        return kFALSE;
      }
      break;
    case AliRsnTarget::kEvent:
      if (!AliRsnTarget::GetCurrentEvent())
      {
        AliError(Form("[%s] expected: AliRsnEvent, passed: [%s]", GetName(), object->ClassName()));
        return kFALSE;
      }
      break;
    default:
      AliError(Form("[%s] Wrong type", GetName()));
      return kFALSE;
  }

  // cast the support object to the types which could be needed
  AliESDtrackCuts *esdCuts = dynamic_cast<AliESDtrackCuts*>(fSupportObject);
  AliRsnPairDef   *pairDef = dynamic_cast<AliRsnPairDef*>(fSupportObject);
    
  // compute value depending on type
  switch (fValueType)
  {
    case kTrackP:
      fComputedValue = useMC ? pSim.Mag() : pRec.Mag();
      break;
    case kTrackPt:
      fComputedValue = useMC ? pSim.Perp() : pRec.Perp();
      break;
    case kTrackEta:
      fComputedValue = useMC ? pSim.Eta() : pRec.Eta();
      break;
    case kPairP1:
      fComputedValue = useMC ? pSim0.Mag() : pRec0.Mag();
      break;
    case kPairP2:
      fComputedValue = useMC ? pSim1.Mag() : pRec1.Mag();
      break;
    case kPairP1t:
      fComputedValue = useMC ? pSim0.Perp() : pRec0.Perp();
      break;
    case kPairP2t:
      fComputedValue = useMC ? pSim1.Perp() : pRec1.Perp();
      break;
    case kPairP1z:
      fComputedValue = useMC ? pSim0.Z() : pRec0.Z();
      break;
    case kPairP2z:
      fComputedValue = useMC ? pSim1.Z() : pRec1.Z();
      break;
    case kPairInvMass:
      fComputedValue = useMC ? pSim.M() : pRec.M();
      break;
    case kPairInvMassRes:
      fComputedValue = (pSim.M() - pRec.M()) / pSim.M();
      break;
    case kPairPt:
      fComputedValue = useMC ? pSim.Perp() : pRec.Perp();
      break;
    case kPairEta:
      fComputedValue = useMC ? pSim.Eta() : pRec.Eta();
      break;
    case kPairMt:
      // for this computation, replace the computed mass with the default mass
      // for doing this, an initialized pairDef is required to get the mass
      if (!pairDef)
      {
        AliError(Form("[%s] Required a correctly initialized PairDef to compute this value", GetName()));
        fComputedValue = 1E+10;
        return kFALSE;
      }
      else
      {
        pRec.SetXYZM(pRec.X(), pRec.Y(), pRec.Z(), pairDef->GetMotherMass());
        pSim.SetXYZM(pSim.X(), pSim.Y(), pSim.Z(), pairDef->GetMotherMass());
        fComputedValue = useMC ? pSim.Mt() : pRec.Mt();
      }
      break;
    case kPairY:
      // for this computation, replace the computed mass with the default mass
      // for doing this, an initialized pairDef is required to get the mass
      if (!pairDef)
      {
        AliError(Form("[%s] Required a correctly initialized PairDef to compute this value", GetName()));
        fComputedValue = 1E+10;
        return kFALSE;
      }
      else
      {
        pRec.SetXYZM(pRec.X(), pRec.Y(), pRec.Z(), pairDef->GetMotherMass());
        pSim.SetXYZM(pSim.X(), pSim.Y(), pSim.Z(), pairDef->GetMotherMass());
        fComputedValue = useMC ? pSim.Rapidity() : pRec.Rapidity();
      }
      break;
    case kPairPhi:
      fComputedValue = useMC ? pSim.Phi() : pRec.Phi();
      break;
    case kPairPtRatio:
      if (useMC)
      {
        fComputedValue  = TMath::Abs(pSim0.Perp() - pSim1.Perp());
        fComputedValue /= TMath::Abs(pSim0.Perp() + pSim1.Perp());
      }
      else
      {
        fComputedValue  = TMath::Abs(pRec0.Perp() - pRec1.Perp());
        fComputedValue /= TMath::Abs(pRec0.Perp() + pRec1.Perp());
      }
      break;
    case kPairDipAngle:
      fComputedValue = useMC ? pSim0.Angle(pSim1.Vect()) : pRec0.Angle(pRec1.Vect());
      fComputedValue = TMath::Abs(TMath::Cos(fComputedValue));
      break;
    case kPairCosThetaStar:
      fComputedValue = mother->CosThetaStar(useMC);
      break;
    case kPairQInv:
      pSim0 -= pSim1;
      pRec0 -= pRec1;
      fComputedValue = useMC ? pSim0.M() : pRec0.M();
      break;
    case kPairAngleToLeading:
      {
        AliRsnEvent *event = AliRsnTarget::GetCurrentEvent();
    	  int ID1 = (mother->GetDaughter(0))->GetID();
    	  int ID2 = (mother->GetDaughter(1))->GetID();
    	  //int leadingID = event->SelectLeadingParticle(0);
        Int_t leadingID = event->GetLeadingParticleID();
    	  if (leadingID == ID1 || leadingID == ID2) return kFALSE;
    	  AliRsnDaughter leadingPart = event->GetDaughter(leadingID);
    	  AliVParticle  *ref = leadingPart.GetRef();
    	  fComputedValue = ref->Phi() - mother->Sum().Phi();
    	  //return angle w.r.t. leading particle in the range -pi/2, 3/2pi
    	  while(fComputedValue >= TMath::Pi()) fComputedValue -= 2*TMath::Pi();
    	  while(fComputedValue < -0.5*TMath::Pi()) fComputedValue += 2*TMath::Pi();
    	  //Printf("%g", fComputedValue);
      }
      break;
    case kEventMult:
      fComputedValue = (Double_t)AliRsnTarget::GetCurrentEvent()->GetMultiplicity(0x0);
      break;
    case kEventMultESDCuts:
      // this value requires an initialized ESDtrackCuts
      if (!esdCuts)
      {
        AliError(Form("[%s] Required a correctly initialized ESDtrackCuts to compute this value", GetName()));
        fComputedValue = 1E+10;
        return kFALSE;
      }
      fComputedValue = (Double_t)AliRsnTarget::GetCurrentEvent()->GetMultiplicity(esdCuts);
      break;
    case kEventLeadingPt:
      {
    	  int leadingID = AliRsnTarget::GetCurrentEvent()->GetLeadingParticleID(); //fEvent->SelectLeadingParticle(0);
    	  if(leadingID >= 0) 
        {
    		  AliRsnDaughter leadingPart = AliRsnTarget::GetCurrentEvent()->GetDaughter(leadingID);
    		  AliVParticle *ref = leadingPart.GetRef();
    		  fComputedValue = ref->Pt();
    	  }
    	  else fComputedValue = 0;
      }
      break;
    case kEventVz:
      fComputedValue = AliRsnTarget::GetCurrentEvent()->GetRef()->GetPrimaryVertex()->GetZ();
      break;
    default:
      AliError(Form("[%s] Invalid value type for this computation", GetName()));
      return kFALSE;
  }
  
  return kTRUE;
}

//_____________________________________________________________________________
void AliRsnValue::Print(Option_t * /*option */) const
{
//
// Print informations about this object
//

  AliInfo("=== VALUE INFO =================================================");
  AliInfo(Form(" Name                  : %s", GetName()));
  AliInfo(Form(" Type                  : %s", GetValueTypeName()));
  AliInfo(Form(" Current computed value: %f", fComputedValue));
  Int_t i;
  for (i = 0; i < fBinArray.GetSize(); i++)
  {
    AliInfo(Form(" Bin limit #%d         = %f", i, fBinArray[i]));
  }
  AliInfo(Form(" Support object        : %s", (fSupportObject ? fSupportObject->ClassName() : " NO SUPPORT")));
  AliInfo("=== END VALUE INFO =============================================");
}
