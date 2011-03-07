//
// *** Class AliRsnMonitorFunctions ***
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

#include "AliLog.h"

#include "AliRsnMother.h"
#include "AliRsnEvent.h"
#include "AliRsnFunction.h"
#include "AliRsnCutSet.h"
#include "AliRsnValue.h"

#include "AliRsnMonitorFunctions.h"

ClassImp(AliRsnMonitorFunctions)

//_____________________________________________________________________________
AliRsnMonitorFunctions::AliRsnMonitorFunctions(const char *name, AliRsnDaughterDef *def) :
   AliRsnMonitor(name, def),
   fFunctions("AliRsnFunction", 0)
{
//
// Default constructor
//

   AliDebug(AliLog::kDebug + 2, "<-");
   AliDebug(AliLog::kDebug + 2, "->");
}

//_____________________________________________________________________________
AliRsnMonitorFunctions::AliRsnMonitorFunctions(const AliRsnMonitorFunctions& copy) :
   AliRsnMonitor(copy),
   fFunctions(copy.fFunctions)
{
//
// Default constructor
//

   AliDebug(AliLog::kDebug + 2, "<-");
   AliDebug(AliLog::kDebug + 2, "->");
}

//_____________________________________________________________________________
AliRsnMonitorFunctions& AliRsnMonitorFunctions::operator=(const AliRsnMonitorFunctions& copy)
{
   AliRsnMonitor::operator=(copy);

   Int_t i, n = copy.fFunctions.GetEntries();
   for (i = 0; i < n; i++) {
      AliRsnFunction *fcn = (AliRsnFunction*)copy.fFunctions[i];
      if (fcn) AddFunction(fcn);
   }

   return (*this);
}

//_____________________________________________________________________________
AliRsnMonitorFunctions::~AliRsnMonitorFunctions()
{
//
// Destructor
//

   AliDebug(AliLog::kDebug + 2, "<-");
   AliDebug(AliLog::kDebug + 2, "->");
}

//_____________________________________________________________________________
void AliRsnMonitorFunctions::Compute()
{
//
// Makes computations using the two passed daughter objects.
// Checks all cuts and then computes the corresponding pair object
// and then fill the list of required values using it.
//

   AliDebug(AliLog::kDebug + 2, "<-");

   TObjArrayIter   nextFcn(&fFunctions);
   AliRsnFunction *fcn = 0x0;

   while ((fcn = (AliRsnFunction*)nextFcn())) {
      fcn->Fill(fDaughter);
   }

   AliDebug(AliLog::kDebug + 2, "->");
}

//_____________________________________________________________________________
void AliRsnMonitorFunctions::Init(const char *prefix, TList *list)
{
//
// Generates needed histograms, giving them a name based on
// the flags defined here, on the pair definition, and attaches
// a prefix to it, according to the argument.
//
// All generated histograms are stored into the output TList.
//

   AliDebug(AliLog::kDebug + 2, "<-");

   Int_t  i;
   TString hName("");
   AliRsnFunction *fcn = 0;
   for (i = 0; i < fFunctions.GetEntries(); i++) {
      fcn = (AliRsnFunction*)fFunctions.At(i);
      hName  = prefix;
      hName += '_';
      hName += GetName();
      hName += '_';
      hName += fcn->GetName();
      if (fcn->IsUsingTH1()) list->Add(fcn->CreateHistogram(hName.Data(), ""));
      else list->Add(fcn->CreateHistogramSparse(hName.Data(), ""));
   }

   AliDebug(AliLog::kDebug + 2, "->");
}

//_____________________________________________________________________________
void AliRsnMonitorFunctions::AddFunction(AliRsnFunction *const fcn)
{
//
// Adds a new computing function
//

   AliDebug(AliLog::kDebug + 2, "<-");

   Int_t size = fFunctions.GetEntries();
   new(fFunctions[size]) AliRsnFunction(*fcn);
   fFunctions.Print();

   AliDebug(AliLog::kDebug + 2, "->");
}
