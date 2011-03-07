//
// *** Class AliRsnMonitorNtuple ***
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

#include "AliRsnMonitorNtuple.h"

ClassImp(AliRsnMonitorNtuple)

//_____________________________________________________________________________
AliRsnMonitorNtuple::AliRsnMonitorNtuple(const char *name, AliRsnDaughterDef *def) :
   AliRsnMonitor(name, def),
   fValues("AliRsnValue", 0),
   fNtuple(0x0)
{
//
// Default constructor
//

   AliDebug(AliLog::kDebug + 2, "<-");
   AliDebug(AliLog::kDebug + 2, "->");
}

//_____________________________________________________________________________
AliRsnMonitorNtuple::AliRsnMonitorNtuple(const AliRsnMonitorNtuple& copy) :
   AliRsnMonitor(copy),
   fValues(copy.fValues),
   fNtuple(copy.fNtuple)
{
//
// Default constructor
//

   AliDebug(AliLog::kDebug + 2, "<-");
   AliDebug(AliLog::kDebug + 2, "->");
}

//_____________________________________________________________________________
AliRsnMonitorNtuple& AliRsnMonitorNtuple::operator=(const AliRsnMonitorNtuple& copy)
{
   AliRsnMonitor::operator=(copy);

   Int_t i, n = copy.fValues.GetEntries();
   for (i = 0; i < n; i++) {
      AliRsnValue *fcn = (AliRsnValue*)copy.fValues[i];
      if (fcn) AddValue(fcn);
   }

   fNtuple = copy.fNtuple;

   return (*this);
}

//_____________________________________________________________________________
AliRsnMonitorNtuple::~AliRsnMonitorNtuple()
{
//
// Destructor
//

   AliDebug(AliLog::kDebug + 2, "<-");
   AliDebug(AliLog::kDebug + 2, "->");
}

//_____________________________________________________________________________
void AliRsnMonitorNtuple::Compute()
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
   Bool_t       computeOK = kFALSE, globalOK = kTRUE;
   for (i = 0; i < n; i++) {
      values[i] = -1E10;
      value = (AliRsnValue*)fValues[i];
      switch (value->GetTargetType()) {
         case AliRsnTarget::kDaughter:
            computeOK = value->Eval(fDaughter);
            break;
         case AliRsnTarget::kEvent:
            computeOK = value->Eval(fDaughter->GetOwnerEvent());
            break;
         default:
            AliError(Form("Allowed targets are mothers and events; cannot use axis '%s' which has target '%s'", value->GetName(), value->GetTargetTypeName()));
            computeOK = kFALSE;
      }
      if (computeOK) 
         values[i] = ((Float_t)value->GetComputedValue());
      else 
         globalOK = kFALSE;
   }

   if (globalOK) fNtuple->Fill(values.GetArray());

   AliDebug(AliLog::kDebug + 2, "->");
}

//_____________________________________________________________________________
void AliRsnMonitorNtuple::Init(const char *prefix, TList *list)
{
//
// Generates needed histograms, giving them a name based on
// the flags defined here, on the pair definition, and attaches
// a prefix to it, according to the argument.
//
// All generated histograms are stored into the output TList.
//

   AliDebug(AliLog::kDebug + 2, "<-");

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

   AliDebug(AliLog::kDebug + 2, "->");
}

//_____________________________________________________________________________
Bool_t AliRsnMonitorNtuple::AddValue(AliRsnValue *const val)
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
   new(fValues[size]) AliRsnValue(*val);

   return kTRUE;
}
