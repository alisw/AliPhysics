//
// *** Class AliRsnPairFunctions ***
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

#include "AliRsnPairFunctions.h"

ClassImp(AliRsnPairFunctions)

//_____________________________________________________________________________
AliRsnPairFunctions::AliRsnPairFunctions(const char *name, AliRsnPairDef *def) :
   AliRsnPair(name, def),
   fFunctions("AliRsnFunction", 0)
{
//
// Default constructor
//
}

//_____________________________________________________________________________
AliRsnPairFunctions::AliRsnPairFunctions(const AliRsnPairFunctions& copy) :
   AliRsnPair(copy),
   fFunctions(copy.fFunctions)
{
//
// Default constructor
//
}

//_____________________________________________________________________________
AliRsnPairFunctions& AliRsnPairFunctions::operator=(const AliRsnPairFunctions& copy)
{
   AliRsnPair::operator=(copy);

   Int_t i, n = copy.fFunctions.GetEntries();
   for (i = 0; i < n; i++) {
      AliRsnFunction *fcn = (AliRsnFunction*)copy.fFunctions[i];
      if (fcn) AddFunction(fcn);
   }

   return (*this);
}

//_____________________________________________________________________________
AliRsnPairFunctions::~AliRsnPairFunctions()
{
//
// Destructor
//
}

//_____________________________________________________________________________
void AliRsnPairFunctions::Compute()
{
//
// Makes computations using the two passed daughter objects.
// Checks all cuts and then computes the corresponding pair object
// and then fill the list of required values using it.
//

   TObjArrayIter   nextFcn(&fFunctions);
   AliRsnFunction *fcn = 0x0;

   while ((fcn = (AliRsnFunction*)nextFcn())) {
      fcn->Fill(&fMother);
   }
}

//_____________________________________________________________________________
void AliRsnPairFunctions::Init(const char *prefix, TList *list)
{
//
// Generates needed histograms, giving them a name based on
// the flags defined here, on the pair definition, and attaches
// a prefix to it, according to the argument.
//
// All generated histograms are stored into the output TList.
//

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
}

//_____________________________________________________________________________
void AliRsnPairFunctions::AddFunction(AliRsnFunction *const fcn)
{
//
// Adds a new computing function
//

   Int_t size = fFunctions.GetEntries();
   new (fFunctions[size]) AliRsnFunction(*fcn);
   fFunctions.Print();
}
