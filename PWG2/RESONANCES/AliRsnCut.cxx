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
#include "AliRsnMother.h"
#include "AliRsnEvent.h"

#include "AliRsnCut.h"

ClassImp(AliRsnCut)

//_________________________________________________________________________________________________
AliRsnCut::AliRsnCut(ETarget target) :
  TNamed(),
  fVarType(kInt),
  fTarget(target),
  fMinI(0),
  fMaxI(0),
  fMinD(0.0),
  fMaxD(0.0),
  fCutValueI(0),
  fCutValueD(0.0),
  fCutResult(kTRUE),
  fEvent(0x0)
{
//
// Default constructor.
//
}

//_________________________________________________________________________________________________
AliRsnCut::AliRsnCut(const AliRsnCut& copy) :
  TNamed(copy),
  fVarType(copy.fVarType),
  fTarget(copy.fTarget),
  fMinI(copy.fMinI),
  fMaxI(copy.fMaxI),
  fMinD(copy.fMinD),
  fMaxD(copy.fMaxD),
  fCutValueI(copy.fCutValueI),
  fCutValueD(copy.fCutValueD),
  fCutResult(copy.fCutResult),
  fEvent(copy.fEvent)
{
//
// Copy constructor.
//
}

//_________________________________________________________________________________________________
AliRsnCut::AliRsnCut
(const char *name, ETarget target, Int_t min, Int_t max) :
  TNamed(name, ""),
  fVarType(kInt),
  fTarget(target),
  fMinI(min),
  fMaxI(max),
  fMinD(0.0),
  fMaxD(0.0),
  fCutValueI(0),
  fCutValueD(0.0),
  fCutResult(kTRUE),
  fEvent(0x0)
{
//
// Constructor with integer values.
// If the cut must check values inside a range,
// both 'value' arguments must be used, and they are, in the order,
// the minimum and maximum of the allowed range.
// If the cut must check a value, the second 'value' argument will never be used.
//
}

//_________________________________________________________________________________________________
AliRsnCut::AliRsnCut
(const char *name, ETarget target, Double_t min, Double_t max) :
  TNamed(name, ""),
  fVarType(kDouble),
  fTarget(target),
  fMinI(0),
  fMaxI(0),
  fMinD(min),
  fMaxD(max),
  fCutValueI(0),
  fCutValueD(0.0),
  fCutResult(kTRUE),
  fEvent(0x0)
{
//
// Constructor with double values.
// If the cut must check values inside a range,
// both 'value' arguments must be used, and they are, in the order,
// the minimum and maximum of the allowed range.
// If the cut must check a value, the second 'value' argument will never be used.
//
}

//_________________________________________________________________________________________________
AliRsnCut& AliRsnCut::operator=(const AliRsnCut& copy)
{
//
// Assignment operator
// don't duplicate memory occupancy for pointer
//

  fVarType   = copy.fVarType;
  fTarget    = copy.fTarget;
  fMinI      = copy.fMinI;
  fMaxI      = copy.fMaxI;
  fMinD      = copy.fMinD;
  fMaxD      = copy.fMaxD;
  fCutValueI = copy.fCutValueI;
  fCutValueD = copy.fCutValueD;
  fCutResult = copy.fCutResult;
  fEvent     = copy.fEvent;

  return (*this);
}

//_________________________________________________________________________________________________
Bool_t AliRsnCut::TargetOK(TObject *obj1, TObject *obj2)
{
//
// This method checks if the expected target and the passed object match.
//

  if (!obj1)
  {
    AliError("Cannot cut on a NULL object!");
    return kFALSE;
  }

  switch (fTarget)
  {
    case kDaughter:
      if (dynamic_cast<AliRsnDaughter*>(obj1) == 0x0)
      {
        AliError(Form("[%s] Target mismatch (obj #1): expected  'AliRsnDaughter', passed '%s'", GetName(), obj1->ClassName()));
        Print();
        return kFALSE;
      }
      break;
    case kMother:
      if (dynamic_cast<AliRsnMother*>(obj1) == 0x0)
      {
        AliError(Form("[%s] Target mismatch (obj #1): expected  'AliRsnMother', passed '%s'", GetName(), obj1->ClassName()));
        Print();
        return kFALSE;
      }
      break;
    case kEvent:
      if (dynamic_cast<AliRsnEvent*>(obj1) == 0x0)
      {
        AliError(Form("[%s] Target mismatch (obj #1): expected  'AliRsnEvent', passed '%s'", GetName(), obj1->ClassName()));
        Print();
        return kFALSE;
      }
      break;
    case kMixEvent:
      if (dynamic_cast<AliRsnEvent*>(obj1) == 0x0)
      {
        AliError(Form("[%s] Target mismatch (obj #1): expected  'AliRsnEvent', passed '%s' an", GetName(), obj1->ClassName()));
        Print();
        return kFALSE;
      }
      if (obj2)
      {
        if (dynamic_cast<AliRsnEvent*>(obj2) == 0x0)
        {
          AliError(Form("[%s] Target mismatch (obj #2): expected  'AliRsnEvent', passed '%s' an", GetName(), obj2->ClassName()));
          Print();
          return kFALSE;
        }
      }
      else
      {
        AliError("Mix-event cuts require 2 not NULL objects");
        Print();
        return kFALSE;
      }
      break;
    default:
      return kTRUE;
  }
  
  return kTRUE;
}

//_________________________________________________________________________________________________
Bool_t AliRsnCut::IsSelected(TObject* /*obj1*/, TObject* /*obj2*/)
{
//
// Virtual cut-checking method for event mixing.
// This method checks only that the target is the oner for mixing.
//

  AliWarning("Single-object cuts are not implemented here.");
  return kTRUE;
}

//_________________________________________________________________________________________________
Bool_t AliRsnCut::OkValue()
{
//
// This method is used when the cut consists in comparing the cut value
// with a reference value to which it must be equal (in case of doubles, 'almost' equal).
// Then, the cut result is kTRUE if the cut value is equal to this reference value.
//

  switch (fVarType) 
  {
    case kInt:
      return OkValueI();
    case kDouble:
      return OkValueD();
    default:
      AliError(Form("fVarType = %d --> not allowed", fVarType));
      return kFALSE;
  }
}

//_________________________________________________________________________________________________
Bool_t AliRsnCut::OkRange()
{
//
// This method is used when the cut consists in an allowed range
// where the cut value must be included to pass the cut.
// Then, the cut result is kTRUE if the cut value is inside this range.
//

  switch (fVarType) 
  {
    case kInt:
      return OkRangeI();
    case kDouble:
      return OkRangeD();
    default:
      AliError(Form("fVarType = %d --> not allowed", fVarType));
      return kFALSE;
  }
}

//_________________________________________________________________________________________________
Bool_t AliRsnCut::OkValueI()
{
//
// This method is used when the cut consists in comparing the cut value
// with a reference integer value to which it must be equal.
//

  // eval result
  fCutResult = (fCutValueI == fMinI);
  
  // print debug message
  AliDebug(AliLog::kDebug + 2, "=== CUT DEBUG ====================================");
  AliDebug(AliLog::kDebug + 2, Form("Cut name     : %s", GetName()));
  AliDebug(AliLog::kDebug + 2, Form("Checked value: %d", fCutValueI));
  AliDebug(AliLog::kDebug + 2, Form("Cut value    : %d", fMinI));
  AliDebug(AliLog::kDebug + 2, Form("Cut result   : %s", (fCutResult ? "PASSED" : "NOT PASSED")));
  AliDebug(AliLog::kDebug + 2, "=== END CUT DEBUG ================================");
  
  return fCutResult;
}

//_________________________________________________________________________________________________
Bool_t AliRsnCut::OkValueD()
{
//
// This method is used when the cut consists in comparing the cut value
// with a reference double value to which it must be equal (or at least, almost).
//

  // eval result
  fCutResult = (TMath::Abs(fCutValueD - fMinD) < 1E-6);
  
  // print debug message
  AliDebug(AliLog::kDebug + 2, "=== CUT DEBUG ====================================");
  AliDebug(AliLog::kDebug + 2, Form("Cut name     : %s", GetName()));
  AliDebug(AliLog::kDebug + 2, Form("Checked value: %f", fCutValueD));
  AliDebug(AliLog::kDebug + 2, Form("Cut value    : %f", fMinD));
  AliDebug(AliLog::kDebug + 2, Form("Cut result   : %s", (fCutResult ? "PASSED" : "NOT PASSED")));
  AliDebug(AliLog::kDebug + 2, "=== END CUT DEBUG ================================");
  
  return fCutResult;
}

//_________________________________________________________________________________________________
Bool_t AliRsnCut::OkRangeI()
{
//
// This method is used when the cut consists in an allowed range
// where the cut value must be included to pass the cut.
// Then, the cut result is kTRUE if the cut value is inside this range.
//

  // eval result
  fCutResult = ((fCutValueI >= fMinI) && (fCutValueI <= fMaxI));
  
  // print debug message
  AliDebug(AliLog::kDebug + 2, "=== CUT DEBUG ====================================");
  AliDebug(AliLog::kDebug + 2, Form("Cut name     : %s", GetName()));
  AliDebug(AliLog::kDebug + 2, Form("Checked value: %d", fCutValueI));
  AliDebug(AliLog::kDebug + 2, Form("Cut range    : %d , %d", fMinI, fMaxI));
  AliDebug(AliLog::kDebug + 2, Form("Cut result   : %s", (fCutResult ? "PASSED" : "NOT PASSED")));
  AliDebug(AliLog::kDebug + 2, "=== END CUT DEBUG ================================");
  
  return fCutResult;
}

//_________________________________________________________________________________________________
Bool_t AliRsnCut::OkRangeD()
{
//
// This method is used when the cut consists in an allowed range
// where the cut value must be included to pass the cut.
// Then, the cut result is kTRUE if the cut value is inside this range.
//

  // eval result
  fCutResult = ((fCutValueD >= fMinD) && (fCutValueD <= fMaxD));
   
  // print debug message
  AliDebug(AliLog::kDebug + 2, "=== CUT DEBUG ====================================");
  AliDebug(AliLog::kDebug + 2, Form("Cut name     : %s", GetName()));
  AliDebug(AliLog::kDebug + 2, Form("Checked value: %f", fCutValueD));
  AliDebug(AliLog::kDebug + 2, Form("Cut range    : %f , %f", fMinD, fMaxD));
  AliDebug(AliLog::kDebug + 2, Form("Cut result   : %s", (fCutResult ? "PASSED" : "NOT PASSED")));
  AliDebug(AliLog::kDebug + 2, "=== END CUT DEBUG ================================");

  return fCutResult;
}

//_________________________________________________________________________________________________
void AliRsnCut::Print(Option_t*) const
{
//
// Override TObject::Print() method
//

  Char_t target[100];
  switch (fTarget)
  {
    case kDaughter: snprintf(target, strlen("DAUGHTER") , "DAUGHTER") ; break;
    case kMother  : snprintf(target, strlen("MOTHER")   , "MOTHER")   ; break;
    case kEvent   : snprintf(target, strlen("EVENT")    , "EVENT")    ; break;
    case kMixEvent: snprintf(target, strlen("MIX EVENT"), "MIX EVENT"); break;
    default       : snprintf(target, strlen("UNDEFINED"), "UNDEFINED"); break;
  }

  AliInfo("=== CUT DETAILS ====================================");
  AliInfo(Form("Cut name     : [%s]", GetName()));
  AliInfo(Form("Cut target   : [%s]", target));
  AliInfo(Form("Cut edges [D]: [%f - %f]", fMinD, fMaxD));
  AliInfo(Form("Cut edges [I]: [%d - %d]", fMinI, fMaxI));
  AliInfo("====================================================");
}

//_________________________________________________________________________________________________
void AliRsnCut::SetEvent(AliRsnEvent *event)
{
//
// Sets the reference event
//

  fEvent = event;
}
