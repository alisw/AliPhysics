//
// Class AliRsnCut
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

#include "AliLog.h"

#include "AliRsnDaughter.h"
#include "AliRsnMCInfo.h"
#include "AliRsnPairParticle.h"
#include "AliRsnPairDef.h"
#include "AliRsnEvent.h"
#include "AliRsnCut.h"

const Double_t AliRsnCut::fgkDSmallNumber = 1e-100;
const Double_t AliRsnCut::fgkDBigNumber = 1e10;
const Int_t    AliRsnCut::fgkIBigNumber = 32767;

ClassImp(AliRsnCut)

//________________________________________________________________________________________________________________
AliRsnCut::AliRsnCut() :
    TNamed(),
    fDMin(-fgkDBigNumber),
    fDMax(fgkDBigNumber),
    fIMin(-fgkIBigNumber),
    fIMax(fgkIBigNumber),
    fUIMin(0),
    fUIMax(2 * (UInt_t) fgkIBigNumber),
    fULMin(0),
    fULMax(2 * (ULong_t) fgkIBigNumber),
    fType(kLastCutType),
    fVarType(kDouble_t)
{
//
// Constructor
//
}

//________________________________________________________________________________________________________________
AliRsnCut::AliRsnCut(const char *name, const char *title, EType type) :
    TNamed(name,title),
    fDMin(-fgkDBigNumber),
    fDMax(fgkDBigNumber),
    fIMin(-fgkIBigNumber),
    fIMax(fgkIBigNumber),
    fUIMin(0),
    fUIMax(2 * (UInt_t) fgkIBigNumber),
    fULMin(0),
    fULMax(2 * (ULong_t) fgkIBigNumber),
    fType(type),
    fVarType(kDouble_t)
{
//
// Constructor with arguments but not limits
//
}

//________________________________________________________________________________________________________________
AliRsnCut::AliRsnCut(const char *name, const char *title, EType type, Double_t min, Double_t max) :
    TNamed(name,title),
    fDMin(min),
    fDMax(max),
    fIMin(-fgkIBigNumber),
    fIMax(fgkIBigNumber),
    fUIMin(0),
    fUIMax(2 * (UInt_t) fgkIBigNumber),
    fULMin(min),
    fULMax(max),
    fType(type),
    fVarType(kDouble_t)
{
//
// Constructor with arguments and limits
//
}

//________________________________________________________________________________________________________________
AliRsnCut::AliRsnCut(const char * name, const char * title, EType type, Int_t min, Int_t max) :
    TNamed(name,title),
    fDMin(-fgkDBigNumber),
    fDMax(fgkDBigNumber),
    fIMin(min),
    fIMax(max),
    fUIMin(0),
    fUIMax(2 * (UInt_t) fgkIBigNumber),
    fULMin(min),
    fULMax(max),
    fType(type),
    fVarType(kInt_t)
{
//
// Constructor with arguments and limits
//
}

//________________________________________________________________________________________________________________
AliRsnCut::AliRsnCut(const char * name, const char * title, EType type, UInt_t min, UInt_t max) :
    TNamed(name,title),
    fDMin(-fgkDBigNumber),
    fDMax(fgkDBigNumber),
    fIMin(-fgkIBigNumber),
    fIMax(fgkIBigNumber),
    fUIMin(min),
    fUIMax(max),
    fULMin(min),
    fULMax(max),
    fType(type),
    fVarType(kUInt_t)
{
//
// Constructor with arguments and limits
//
}

//________________________________________________________________________________________________________________
AliRsnCut::AliRsnCut(const char * name, const char * title, EType type, ULong_t min, ULong_t max) :
    TNamed(name,title),
    fDMin(-fgkDBigNumber),
    fDMax(fgkDBigNumber),
    fIMin(-fgkIBigNumber),
    fIMax(fgkIBigNumber),
    fUIMin(min),
    fUIMax(max),
    fULMin(min),
    fULMax(max),
    fType(type),
    fVarType(kUInt_t)
{
//
// Constructor with arguments and limits
//
}

//________________________________________________________________________________________________________________
AliRsnCut::~AliRsnCut()
{
//
// Destructor.
// Does absolutely nothing.
//
}

//________________________________________________________________________________________________________________
Bool_t AliRsnCut::IsBetween(const Double_t &theValue)
{
//
// Interval check.
// Question: "Is the argument included between fDMin and fDMax?"
// (not implemented for integer values because usually it is not used with them)
//
  return ((theValue >= fDMin) && (theValue <= fDMax));
}

//________________________________________________________________________________________________________________
Bool_t AliRsnCut::IsBetween(const Int_t &theValue)
{
//
// Interval check.
// Question: "Is the argument included between fDMin and fDMax?"
// (not implemented for integer values because usually it is not used with them)
//
  return ((theValue >= fIMin) && (theValue <= fIMax));
}

//________________________________________________________________________________________________________________
Bool_t AliRsnCut::MatchesValue(const Int_t &theValue)
{
//
// Reference check.
// Question: "Is the argument equal to fIMin?" (fIMax is assumed never used)
//
  return (theValue == fIMin);
}

//________________________________________________________________________________________________________________
Bool_t AliRsnCut::MatchesValue(const UInt_t &theValue)
{
//
// Reference check.
// Question: "Is the argument equal to fUIMin?" (fUIMax is assumed never used)
//
  return (theValue == fUIMin);
}

//________________________________________________________________________________________________________________
Bool_t AliRsnCut::MatchesValue(const ULong_t &theValue)
{
  //
// Reference check.
// Question: "Is the argument equal to fUIMin?" (fUIMax is assumed never used)
  //
  return (theValue == fULMin);
}

//________________________________________________________________________________________________________________
Bool_t AliRsnCut::MatchesValue(const Double_t &theValue)
{
//
// Reference check.
// Question: "Is the argument reasonably close to fDMin?" (fDMax is assumed never used)
// Here, "reasonably close" means that the difference is smaller than the
// 'fgkSmallNumber' global static data member of this class
//
  return (TMath::Abs(theValue - fDMin) < fgkDSmallNumber);
}

//________________________________________________________________________________________________________________
void AliRsnCut::SetCutValues(EType type, const Double_t &theValue, const Double_t &theValue2)
{
//
// (Re)assignment of cut values
//
  fType = type;
  fDMin = theValue;
  fDMax = theValue2;
}

//________________________________________________________________________________________________________________
void AliRsnCut::SetCutValues(EType type, const Int_t &theValue, const Int_t &theValue2)
{
//
// (Re)assignment of cut values
//
  fType = type;
  fIMin = theValue;
  fIMax = theValue2;
}

//________________________________________________________________________________________________________________
void AliRsnCut::SetCutValues(EType type, const UInt_t &theValue, const UInt_t &theValue2)
{
//
// (Re)assignment of cut values
//
  fType = type;
  fUIMin = theValue;
  fUIMax = theValue2;
}

//________________________________________________________________________________________________________________
void AliRsnCut::SetCutValues(EType type, const ULong_t &theValue, const ULong_t &theValue2)
{
  //
// (Re)assignment of cut values
  //
  fType = type;
  fULMin = theValue;
  fULMax = theValue2;
}

//________________________________________________________________________________________________________________
Bool_t AliRsnCut::IsSelected(ETarget type, AliRsnDaughter *daughter)
{
//
// Core of the whole class.
// According to the kind of cut selected in the enumeration,
// checks the cut taking the right values from the argument.
// Depending on the second argument type, only some cuts are checked
// (the ones for that type of object), otherwise kTRUE is returned in order
// not to act as a cleaning factor for an AND with other cuts.
//
  AliDebug(AliLog::kDebug, "<-");
  AliRsnMCInfo *mcinfo = daughter->GetMCInfo();

  // check type
  if (type != kParticle)
  {
    AliWarning(Form("Mismatch: type = %d (expected %d), class type = %s (expected AliRsnDaughter)", type, kParticle, daughter->ClassName()));
    return kTRUE;
  }

  // utility variables
  AliRsnPID::EType pidType;
  Double_t prob;
  Int_t pdg;
  Bool_t cut;

  switch (fType)
  {
    case kMomentum:
      return IsBetween(daughter->P());
    case kTransMomentum:
      return IsBetween(daughter->Pt());
    case kEta:
      return IsBetween(daughter->Eta());
    case kRadialImpactParam:
      return IsBetween(daughter->Dr());
    case kMomentumMC:
      if (mcinfo) return IsBetween(mcinfo->P());
      else return kTRUE;
    case kTransMomentumMC:
      if (mcinfo) return IsBetween(mcinfo->P());
      else return kTRUE;
    case kStatus:
      return daughter->CheckFlag(fUIMin);
    case kChargePos:
      return (daughter->Charge() > 0);
    case kChargeNeg:
      return (daughter->Charge() < 0);
    case kPIDType:
    case kPIDProb:
      pidType = daughter->PIDType(prob);
      if (fType == kPIDType) return MatchesValue((Int_t) pidType);
      if (fType == kPIDProb) return IsBetween(prob);
    case kTruePID:
      pdg = TMath::Abs(mcinfo->PDG());
      cut = MatchesValue(pdg);
      //AliError(Form("PDG = %d -- CUT = %s", pdg, (cut ? "passed" : "not passed")));
      if (mcinfo) return MatchesValue((Int_t) TMath::Abs(mcinfo->PDG()));
      else return kTRUE;
    case kEtaMC:
      if (mcinfo) return IsBetween(mcinfo->Eta());
      else return kTRUE;
    case kIsPrimary:
      if (mcinfo) return (mcinfo->Mother() < 0);
      else return kTRUE;
    case kNSigma:
      return IsBetween(daughter->NSigmaToVertex());
      /*
      case kEsdNSigmaCalculate:
      return IsBetween (daughter->GetESDInfo()->GetNSigmaCalculate());
      */
    default:
      AliWarning("Requested a cut which cannot be applied to a single track");
      return kTRUE;
  }

  return kTRUE;
}

//________________________________________________________________________________________________________________
Bool_t AliRsnCut::IsSelected(ETarget type, AliRsnPairParticle * pair)
{
  AliDebug(AliLog::kDebug, "<-");

  // check type
  if (type != kPair)
  {
    AliWarning(Form("Mismatch: type = %d (expected %d), class type = %s (expected AliRsnPairParticle)", type, kPair, pair->ClassName()));
    return kTRUE;
  }

  switch (fType)
  {
    case kMomentum:
      return IsBetween(pair->GetP());
    case kTransMomentum:
      return IsBetween(pair->GetPt());
    case kEta:
      return IsBetween(pair->GetEta());
    case kEtaMC:
      return IsBetween(pair->GetEtaMC());
    case kMomentumMC:
      return IsBetween(pair->GetPMC());
    case kTransMomentumMC:
      return IsBetween(pair->GetPtMC());
    case kIsLabelEqual:
      return pair->IsLabelEqual();
    case kIsTruePair:
      return pair->IsTruePair(fIMin);
    default:
      AliWarning("Requested a cut which cannot be applied to a pair");
      return kTRUE;
  }

  return kTRUE;
}

//________________________________________________________________________________________________________________
Bool_t AliRsnCut::IsSelected(ETarget type, AliRsnEvent * event)
{
  AliDebug(AliLog::kDebug, "<-");

  // check type
  if (type != kEvent)
  {
    AliWarning(Form("Mismatch: type = %d (expected %d), class type = %s (expected AliRsnEvent)", type, kEvent, event->ClassName()));
    return kTRUE;
  }

  switch (fType)
  {
    case kMultiplicity:
      return IsBetween((Int_t) event->GetMultiplicity());
    default:
      AliWarning("Requested a cut which cannot be applied to an event");
      return kTRUE;
  }

  return kTRUE;
}

//________________________________________________________________________________________________________________
Bool_t AliRsnCut::IsSelected(ETarget type, AliRsnEvent * ev1, AliRsnEvent * ev2)
{
  AliDebug(AliLog::kDebug, "<-");

  // check type
  if (type != kMixEvent)
  {
    AliWarning(Form("Mismatch: type = %d (expected %d)", type, kMixEvent));
    return kTRUE;
  }

  Double_t valueD, mult1, mult2;
  Int_t    valueI;

  switch (fType)
  {
    case kMultiplicityDifference:
      valueI = TMath::Abs(ev1->GetMultiplicity() - ev2->GetMultiplicity());
      return IsBetween((Int_t)valueI);
    case kMultiplicityRatio:
      mult1 = (Double_t)ev1->GetMultiplicity();
      mult2 = (Double_t)ev2->GetMultiplicity();
      if (mult1 == 0.0  && mult2 == 0.0) return kTRUE;
      valueD = 100.0 * TMath::Abs(mult1 - mult2) / TMath::Max(mult1, mult2); // in %
      return IsBetween((Double_t)valueD);
    case kVzDifference:
      valueD = TMath::Abs(ev1->GetVz() - ev2->GetVz());
      return IsBetween((Double_t)valueD);
    case kPhiMeanDifference:
      valueD = TMath::Abs(ev1->GetPhiMean() - ev2->GetPhiMean());
      if (valueD > 180.0) valueD = 360.0 - valueD;
      return IsBetween((Double_t)valueD);
    default:
      AliWarning("Requested a cut which cannot be applied to an event");
      return kTRUE;
  }

  return kTRUE;
}

//________________________________________________________________________________________________________________
void AliRsnCut::PrintAllValues()
{
  AliInfo(Form("fType=%d fVarType=%d",fType,fVarType));
  AliInfo(Form("fDMin=%.2e fDMax=%.2e",fDMin,fDMax));
  AliInfo(Form("fIMin=%d fIMax=%d",fIMin,fIMax));
  AliInfo(Form("fUIMin=%d fUIMax=%d",fUIMin,fUIMax));
  AliInfo(Form("fULMin=%d fULMax=%d",fULMin,fULMax));
}
