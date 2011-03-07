//
// *** Class AliRsnCut ***
//
// Cut base class: all other cuts inherit from it.
// The 'core' of the class is the method "IsSelected()" which
// must be overloaded by any specific cut implementation.
//
// This class provides some default instruments to check values
// agains a reference or an allowed range, in order to permit
// a unique way to execute such kind of checks.
//
// authors: Alberto Pulvirenti (alberto.pulvirenti@ct.infn.it)
//          Martin Vala (martin.vala@cern.ch)
//

#include "AliRsnCut.h"

ClassImp(AliRsnCut)

//______________________________________________________________________________
AliRsnCut::AliRsnCut(const char *name, RSNTARGET target) :
   AliRsnTarget(name, target),
   fMinI(0),
   fMaxI(0),
   fMinD(0.),
   fMaxD(0.),
   fCutValueI(0),
   fCutValueD(0.0),
   fCutResult(kTRUE)
{
//
// Default constructor.
//
}

//______________________________________________________________________________
AliRsnCut::AliRsnCut
(const char *name, RSNTARGET target, Int_t imin, Int_t imax, Double_t dmin, Double_t dmax) :
   AliRsnTarget(name, target),
   fMinI(imin),
   fMaxI(imax),
   fMinD(dmin),
   fMaxD(dmax),
   fCutValueI(0),
   fCutValueD(0.0),
   fCutResult(kTRUE)
{
//
// Constructor with arguments.
// This is provided to allow a quick setting of all data members.
//
}

//______________________________________________________________________________
AliRsnCut::AliRsnCut
(const char *name, RSNTARGET target, Double_t dmin, Double_t dmax, Int_t imin, Int_t imax) :
   AliRsnTarget(name, target),
   fMinI(imin),
   fMaxI(imax),
   fMinD(dmin),
   fMaxD(dmax),
   fCutValueI(0),
   fCutValueD(0.0),
   fCutResult(kTRUE)
{
//
// Constructor with arguments.
// This is provided to allow a quick setting of all data members.
//
}

//______________________________________________________________________________
AliRsnCut::AliRsnCut(const AliRsnCut& copy) :
   AliRsnTarget(copy),
   fMinI(copy.fMinI),
   fMaxI(copy.fMaxI),
   fMinD(copy.fMinD),
   fMaxD(copy.fMaxD),
   fCutValueI(copy.fCutValueI),
   fCutValueD(copy.fCutValueD),
   fCutResult(copy.fCutResult)
{
//
// Copy constructor.
// Don't duplicate memory occupancy for pointer
//
}

//______________________________________________________________________________
AliRsnCut& AliRsnCut::operator=(const AliRsnCut& copy)
{
//
// Assignment operator.
// Don't duplicate memory occupancy for pointer
//

   AliRsnTarget::operator=(copy);

   fMinI      = copy.fMinI;
   fMaxI      = copy.fMaxI;
   fMinD      = copy.fMinD;
   fMaxD      = copy.fMaxD;
   fCutValueI = copy.fCutValueI;
   fCutValueD = copy.fCutValueD;
   fCutResult = copy.fCutResult;

   return (*this);
}

//______________________________________________________________________________
Bool_t AliRsnCut::IsSelected(TObject* /*object*/)
{
//
// Virtual cut-checking method.
// In this implementation, it does nothing, and all classes
// inheriting from this, should provide a proper implementation
// which must return kTRUE if the cut is passed, and kFALSE otherwise.
//

   AliWarning("This virtual function must be implemented properly");
   return kTRUE;
}

//______________________________________________________________________________
Bool_t AliRsnCut::OkValueI()
{
//
// This method is used to compare a value with a reference.
// In the case of integers, the equality must be exact.
//

   // eval result
   fCutResult = (fCutValueI == fMinI);

   // print debug message
   AliDebug(AliLog::kDebug + 2, "=== CUT DEBUG ========================================================");
   AliDebug(AliLog::kDebug + 2, Form("Cut name     : %s", GetName()));
   AliDebug(AliLog::kDebug + 2, Form("Checked value: %d", fCutValueI));
   AliDebug(AliLog::kDebug + 2, Form("Cut value    : %d", fMinI));
   AliDebug(AliLog::kDebug + 2, Form("Cut result   : %s", (fCutResult ? "PASSED" : "NOT PASSED")));
   AliDebug(AliLog::kDebug + 2, "=== END CUT DEBUG ====================================================");

   return fCutResult;
}

//______________________________________________________________________________
Bool_t AliRsnCut::OkValueD()
{
//
// This method is used to compare a value with a reference.
// In the case of doubles, the equality consists in being very close.
//

   // eval result
   fCutResult = (TMath::Abs(fCutValueD - fMinD) < 1E-6);

   // print debug message
   AliDebug(AliLog::kDebug + 2, "=== CUT DEBUG =======================================================");
   AliDebug(AliLog::kDebug + 2, Form("Cut name     : %s", GetName()));
   AliDebug(AliLog::kDebug + 2, Form("Checked value: %f", fCutValueD));
   AliDebug(AliLog::kDebug + 2, Form("Cut value    : %f", fMinD));
   AliDebug(AliLog::kDebug + 2, Form("Cut result   : %s", (fCutResult ? "PASSED" : "NOT PASSED")));
   AliDebug(AliLog::kDebug + 2, "=== END CUT DEBUG ===================================================");

   return fCutResult;
}

//______________________________________________________________________________
Bool_t AliRsnCut::OkRangeI()
{
//
// This method is used to compare a value with an integer range.
//

   // eval result
   fCutResult = ((fCutValueI >= fMinI) && (fCutValueI <= fMaxI));

   // print debug message
   AliDebug(AliLog::kDebug + 2, "=== CUT DEBUG ========================================================");
   AliDebug(AliLog::kDebug + 2, Form("Cut name     : %s", GetName()));
   AliDebug(AliLog::kDebug + 2, Form("Checked value: %d", fCutValueI));
   AliDebug(AliLog::kDebug + 2, Form("Cut range    : %d , %d", fMinI, fMaxI));
   AliDebug(AliLog::kDebug + 2, Form("Cut result   : %s", (fCutResult ? "PASSED" : "NOT PASSED")));
   AliDebug(AliLog::kDebug + 2, "=== END CUT DEBUG ====================================================");

   return fCutResult;
}

//______________________________________________________________________________
Bool_t AliRsnCut::OkRangeD()
{
//
// This method is used to compare a value with a double-float range.
//

   // eval result
   fCutResult = ((fCutValueD >= fMinD) && (fCutValueD <= fMaxD));

   // print debug message
   AliDebug(AliLog::kDebug + 2, "=== CUT DEBUG ========================================================");
   AliDebug(AliLog::kDebug + 2, Form("Cut name     : %s", GetName()));
   AliDebug(AliLog::kDebug + 2, Form("Checked value: %f", fCutValueD));
   AliDebug(AliLog::kDebug + 2, Form("Cut range    : %f , %f", fMinD, fMaxD));
   AliDebug(AliLog::kDebug + 2, Form("Cut result   : %s", (fCutResult ? "PASSED" : "NOT PASSED")));
   AliDebug(AliLog::kDebug + 2, "=== END CUT DEBUG ====================================================");

   return fCutResult;
}

//______________________________________________________________________________
void AliRsnCut::Print(Option_t*) const
{
//
// Override TObject::Print() method,
// and print some useful info about the cut general parameters.
//

   AliInfo("=== CUT DETAILS ====================================");
   AliInfo(Form("Cut name     : [%s]", GetName()));
   AliInfo(Form("Cut target   : [%s]", GetTargetTypeName()));
   AliInfo(Form("Cut edges [D]: [%f - %f]", fMinD, fMaxD));
   AliInfo(Form("Cut edges [I]: [%d - %d]", fMinI, fMaxI));
   AliInfo("====================================================");
}
