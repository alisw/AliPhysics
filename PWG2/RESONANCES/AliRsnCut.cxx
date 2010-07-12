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

#include "AliRsnCut.h"

ClassImp(AliRsnCut)

//_________________________________________________________________________________________________
AliRsnCut::AliRsnCut() :
    TNamed(),
    fVarType(kInt),
    fMinI(0),
    fMaxI(0),
    fMinU(0),
    fMaxU(0),
    fMinD(0.0),
    fMaxD(0.0),
    fCutValueI(0),
    fCutValueU(0),
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
    fMinI(copy.fMinI),
    fMaxI(copy.fMaxI),
    fMinU(copy.fMinU),
    fMaxU(copy.fMaxU),
    fMinD(copy.fMinD),
    fMaxD(copy.fMaxD),
    fCutValueI(copy.fCutValueI),
    fCutValueU(copy.fCutValueU),
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
(const char *name, Int_t min, Int_t max) :
    TNamed(name, ""),
    fVarType(kInt),
    fMinI(min),
    fMaxI(max),
    fMinU(0),
    fMaxU(0),
    fMinD(0.0),
    fMaxD(0.0),
    fCutValueI(0),
    fCutValueU(0),
    fCutValueD(0.0),
    fCutResult(kTRUE),
    fEvent(0x0)
{
//
// Constructor.
// If the cut must check values inside a range,
// both 'value' arguments must be used, and they are, in the order,
// the minimum and maximum of the allowed range.
// If the cut must check a value, the second 'value' argument will never be used.
//
}

//_________________________________________________________________________________________________
AliRsnCut::AliRsnCut
(const char *name, ULong_t min, ULong_t max) :
    TNamed(name, ""),
    fVarType(kULong),
    fMinI(0),
    fMaxI(0),
    fMinU(min),
    fMaxU(max),
    fMinD(0.0),
    fMaxD(0.0),
    fCutValueI(0),
    fCutValueU(0),
    fCutValueD(0.0),
    fCutResult(kTRUE),
    fEvent(0x0)
{
//
// Constructor.
// If the cut must check values inside a range,
// both 'value' arguments must be used, and they are, in the order,
// the minimum and maximum of the allowed range.
// If the cut must check a value, the second 'value' argument will never be used.
//
}

//_________________________________________________________________________________________________
AliRsnCut::AliRsnCut
(const char *name, Double_t min, Double_t max) :
    TNamed(name, ""),
    fVarType(kDouble),
    fMinI(0),
    fMaxI(0),
    fMinU(0),
    fMaxU(0),
    fMinD(min),
    fMaxD(max),
    fCutValueI(0),
    fCutValueU(0),
    fCutValueD(0.0),
    fCutResult(kTRUE),
    fEvent(0x0)
{
//
// Constructor.
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
  fMinI      = copy.fMinI;
  fMaxI      = copy.fMaxI;
  fMinD      = copy.fMinD;
  fMaxD      = copy.fMaxD;
  fMinU      = copy.fMinU;
  fMaxU      = copy.fMaxU;
  fCutValueI = copy.fCutValueI;
  fCutValueD = copy.fCutValueD;
  fCutValueU = copy.fCutValueU;
  fCutResult = copy.fCutResult;
  fEvent     = copy.fEvent;

  return (*this);
}

//_________________________________________________________________________________________________
Bool_t AliRsnCut::IsSelected(ETarget /*tgt*/, AliRsnDaughter* /*track*/)
{
//
// Virtual cut-checking method.
// In base class, these methods compare the argument type
// with the defined target, in order to detect a mismatch
//

  AliWarning("This cut does not provide checks on AliRsnDaughter. This function will return kTRUE");
  return kTRUE;
}

//_________________________________________________________________________________________________
Bool_t AliRsnCut::IsSelected(ETarget /*tgt*/, AliRsnPairParticle* /*pair*/)
{
//
// Virtual cut-checking method.
// In base class, these methods compare the argument type
// with the defined target, in order to detect a mismatch
//

  AliWarning("This cut does not provide checks on AliRsnPairParticle. This function will return kTRUE");
  return kTRUE;
}

//_________________________________________________________________________________________________
Bool_t AliRsnCut::IsSelected(ETarget /*tgt*/, AliRsnEvent* /*event*/)
{
//
// Virtual cut-checking method.
// In base class, these methods compare the argument type
// with the defined target, in order to detect a mismatch
//

  AliWarning("This cut does not provide checks on AliRsnEvent. This function will return kTRUE");
  return kTRUE;
}

//_________________________________________________________________________________________________
Bool_t AliRsnCut::IsSelected(ETarget /*tgt*/, AliRsnEvent* /*ev1*/, AliRsnEvent* /*ev2*/)
{
//
// Virtual cut-checking method.
// In base class, these methods compare the argument type
// with the defined target, in order to detect a mismatch
//

  AliWarning("This cut does not provide checks on two AliRsnEvent's. This function will return kTRUE");
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

  switch (fVarType) {
  case kInt:
    // eval result
    fCutResult = (fCutValueI == fMinI);
    // print debug message
    AliDebug(AliLog::kDebug + 3, "=== CUT DEBUG ====================================");
    AliDebug(AliLog::kDebug + 3, Form("Cut name     : %s", GetName()));
    AliDebug(AliLog::kDebug + 3, Form("Checked value: %d", fCutValueI));
    AliDebug(AliLog::kDebug + 3, Form("Cut value    : %d", fMinI));
    AliDebug(AliLog::kDebug + 3, Form("Cut result   : %s", (fCutResult ? "PASSED" : "NOT PASSED")));
    AliDebug(AliLog::kDebug + 3, "=== END CUT DEBUG ================================");
    break;
  case kULong:
    // eval result
    fCutResult = (fCutValueU == fMinU);
    // print debug message
    AliDebug(AliLog::kDebug + 3, "=== CUT DEBUG ====================================");
    AliDebug(AliLog::kDebug + 3, Form("Cut name     : %s" , GetName()));
    AliDebug(AliLog::kDebug + 3, Form("Checked value: %lu", fCutValueU));
    AliDebug(AliLog::kDebug + 3, Form("Cut value    : %lu", fMinU));
    AliDebug(AliLog::kDebug + 3, Form("Cut result   : %s" , (fCutResult ? "PASSED" : "NOT PASSED")));
    AliDebug(AliLog::kDebug + 3, "=== END CUT DEBUG ================================");
    break;
  case kDouble:
    // eval result
    fCutResult = (TMath::Abs(fCutValueD - fMinD) < 1E-6);
    // print debug message
    AliDebug(AliLog::kDebug + 3, "=== CUT DEBUG ====================================");
    AliDebug(AliLog::kDebug + 3, Form("Cut name     : %s", GetName()));
    AliDebug(AliLog::kDebug + 3, Form("Checked value: %f", fCutValueD));
    AliDebug(AliLog::kDebug + 3, Form("Cut value    : %f", fMinD));
    AliDebug(AliLog::kDebug + 3, Form("Cut result   : %s", (fCutResult ? "PASSED" : "NOT PASSED")));
    AliDebug(AliLog::kDebug + 3, "=== END CUT DEBUG ================================");
    break;
  default:
    AliError(Form("fVarType = %d --> not allowed", fVarType));
    return kFALSE;
  }

  return fCutResult;
}

//_________________________________________________________________________________________________
Bool_t AliRsnCut::OkRange()
{
//
// This method is used when the cut consists in an allowed range
// where the cut value must be included to pass the cut.
// Then, the cut result is kTRUE if the cut value is inside this range.
//

  switch (fVarType) {
  case kInt:
    // eval result
    fCutResult = ((fCutValueI >= fMinI) && (fCutValueI <= fMaxI));
    // print debug message
    AliDebug(AliLog::kDebug + 3, "=== CUT DEBUG ====================================");
    AliDebug(AliLog::kDebug + 2, Form("Cut name     : %s", GetName()));
    AliDebug(AliLog::kDebug + 2, Form("Checked value: %d", fCutValueI));
    AliDebug(AliLog::kDebug + 2, Form("Cut range    : %d , %d", fMinI, fMaxI));
    AliDebug(AliLog::kDebug + 2, Form("Cut result   : %s", (fCutResult ? "PASSED" : "NOT PASSED")));
    AliDebug(AliLog::kDebug + 3, "=== END CUT DEBUG ================================");
    break;
  case kULong:
    // eval result
    fCutResult = ((fCutValueU >= fMinU) && (fCutValueU <= fMaxU));
    // print debug message
    AliDebug(AliLog::kDebug + 3, "=== CUT DEBUG ====================================");
    AliDebug(AliLog::kDebug + 2, Form("Cut name     : %s"       , GetName()));
    AliDebug(AliLog::kDebug + 2, Form("Checked value: %lu"      , fCutValueU));
    AliDebug(AliLog::kDebug + 2, Form("Cut range    : %lu , %lu", fMinU, fMaxU));
    AliDebug(AliLog::kDebug + 2, Form("Cut result   : %s"       , (fCutResult ? "PASSED" : "NOT PASSED")));
    AliDebug(AliLog::kDebug + 3, "=== END CUT DEBUG ================================");
    break;
  case kDouble:
    // eval result
    fCutResult = ((fCutValueD >= fMinD) && (fCutValueD <= fMaxD));
    // print debug message
    AliDebug(AliLog::kDebug + 3, "=== CUT DEBUG ====================================");
    AliDebug(AliLog::kDebug + 2, Form("Cut name     : %s", GetName()));
    AliDebug(AliLog::kDebug + 2, Form("Checked value: %f", fCutValueD));
    AliDebug(AliLog::kDebug + 2, Form("Cut range    : %f , %f", fMinD, fMaxD));
    AliDebug(AliLog::kDebug + 2, Form("Cut result   : %s", (fCutResult ? "PASSED" : "NOT PASSED")));
    AliDebug(AliLog::kDebug + 3, "=== END CUT DEBUG ================================");
    break;
  default:
    AliError(Form("fVarType = %d --> not allowed", fVarType));
    return kFALSE;
  }

  return fCutResult;
}

