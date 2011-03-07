//
// *** Class AliRsnPairNtuple ***
//
// "Core" method for defining the work on a pari of particles.
// For one analysis, one must setup one of this for each pair he wants to analyze,
// adding to it all analysis which he desires to do.
// Here he defines the cuts, and the particle types and charges, and can add
// functions which do different operations on the same pair, and some binning
// with respect to some kinematic variables (eta, momentum)
//
// authors: A. Pulvirenti (email: alberto.pulvirenti@ct.infn.it)
//          M. Vala (email: martin.vala@cern.ch)
//

#include <TList.h>
#include <TNtuple.h>

#include "AliLog.h"

#include "AliRsnMother.h"
#include "AliRsnEvent.h"
#include "AliRsnFunction.h"
#include "AliRsnCutSet.h"
#include "AliRsnValue.h"

#include "AliRsnPairNtuple.h"

ClassImp(AliRsnPairNtuple)

//_____________________________________________________________________________
AliRsnPairNtuple::AliRsnPairNtuple(const char *name, AliRsnPairDef *def) :
   AliRsnPair(name, def),
   fValues("AliRsnValue", 0),
   fNtuple(0x0)
{
//
// Default constructor
//
}

//_____________________________________________________________________________
AliRsnPairNtuple::AliRsnPairNtuple(const AliRsnPairNtuple& copy) :
   AliRsnPair(copy),
   fValues(copy.fValues),
   fNtuple(copy.fNtuple)
{
//
// Default constructor
//
}

//_____________________________________________________________________________
AliRsnPairNtuple& AliRsnPairNtuple::operator=(const AliRsnPairNtuple& copy)
{
   AliRsnPair::operator=(copy);

   Int_t i, n = copy.fValues.GetEntries();
   for (i = 0; i < n; i++) {
      AliRsnValue *fcn = (AliRsnValue*)copy.fValues[i];
      if (fcn) AddValue(fcn);
   }

   fNtuple = copy.fNtuple;

   return (*this);
}

//_____________________________________________________________________________
AliRsnPairNtuple::~AliRsnPairNtuple()
{
//
// Destructor
//
}

//_____________________________________________________________________________
void AliRsnPairNtuple::Compute()
{
//
// Makes computations using the two passed daughter objects.
// Checks all cuts and then computes the corresponding pair object
// and then fill the list of required values using it.
//

   AliDebug(AliLog::kDebug + 2, "<-");

   // compute all values
   Int_t        i, n = fValues.GetEntries();
   TArrayF      values(n);
   AliRsnValue *value = 0x0;
   Bool_t       computeOK = kFALSE;
   for (i = 0; i < n; i++) {
      values[i] = -1E10;
      value = (AliRsnValue*)fValues[i];
      switch (value->GetTargetType()) {
         case AliRsnTarget::kMother:
            value->SetSupportObject(fPairDef);
            computeOK = value->Eval(&fMother);
            break;
         case AliRsnTarget::kEvent:
            computeOK = value->Eval(fMother.GetRefEvent());
            break;
         default:
            computeOK = kFALSE;
      }
      if (computeOK) values[i] = ((Float_t)value->GetComputedValue());
   }

   fNtuple->Fill(values.GetArray());

   AliDebug(AliLog::kDebug + 2, "->");
}

//_____________________________________________________________________________
void AliRsnPairNtuple::Init(const char *prefix, TList *list)
{
//
// Generates needed histograms, giving them a name based on
// the flags defined here, on the pair definition, and attaches
// a prefix to it, according to the argument.
//
// All generated histograms are stored into the output TList.
//

   TString nameList("");

   Int_t        i, n = fValues.GetEntries();
   AliRsnValue *val = 0;
   for (i = 0; i < n; i++) {
      val = (AliRsnValue*)fValues.At(i);
      nameList += val->GetName();
      if (i < n - 1) nameList += ':';
   }

   if (fNtuple) delete fNtuple;
   fNtuple = new TNtuple(Form("%sntp%s", prefix, GetName()), "", nameList.Data());
   if (list) list->Add(fNtuple);
}

//_____________________________________________________________________________
Bool_t AliRsnPairNtuple::AddValue(AliRsnValue* const val)
{
//
// Adds a new computing function.
//

   RSNTARGET target = val->GetTargetType();
   if (target != AliRsnTarget::kMother && target != AliRsnTarget::kEvent) {
      AliError(Form("Allowed targets are mothers and events; cannot use axis '%s' which has target '%s'", val->GetName(), val->GetTargetTypeName()));
      return kFALSE;
   }

   Int_t size = fValues.GetEntries();
   new (fValues[size]) AliRsnValue(*val);

   return kTRUE;
}
