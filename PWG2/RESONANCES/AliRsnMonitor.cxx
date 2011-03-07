//
// *** Class AliRsnMonitor ***
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

#include "AliLog.h"

#include "AliRsnTarget.h"
#include "AliRsnEvent.h"
#include "AliRsnValue.h"

#include "AliRsnMonitor.h"

ClassImp(AliRsnMonitor)

//_____________________________________________________________________________
AliRsnMonitor::AliRsnMonitor(const char *name, AliRsnDaughterDef *def) :
   TNamed(name, ""),
   fOnlyTrue(kFALSE),
   fCount(0),
   fDaughterDef(def),
   fCuts(Form("cuts_%s", name), AliRsnTarget::kDaughter),
   fDaughter(0)
{
//
// Default constructor
//
}

//_____________________________________________________________________________
AliRsnMonitor::AliRsnMonitor(const AliRsnMonitor& copy) :
   TNamed(copy),
   fOnlyTrue(copy.fOnlyTrue),
   fCount(copy.fCount),
   fDaughterDef(copy.fDaughterDef),
   fCuts(copy.fCuts),
   fDaughter(copy.fDaughter)
{
//
// Default constructor
//
}

//_____________________________________________________________________________
AliRsnMonitor& AliRsnMonitor::operator=(const AliRsnMonitor& copy)
{
   fOnlyTrue = copy.fOnlyTrue;
   fCount = copy.fCount;
   fCuts = copy.fCuts;
   fDaughter = copy.fDaughter;
   fDaughterDef = copy.fDaughterDef;

   return (*this);
}

//_____________________________________________________________________________
AliRsnMonitor::~AliRsnMonitor()
{
//
// Destructor
//
}

//_____________________________________________________________________________
void AliRsnMonitor::Print(Option_t* /*option*/) const
{
//
// Prints info about pair
//
}

//_____________________________________________________________________________
Bool_t AliRsnMonitor::Fill(AliRsnDaughter *daughter)
{
//
// Sets the two passed daughters to the AliRsnMother data member of this object
// which is used to perform all computations to fill the value list.
// This operation is done successfully only when the first passed object matches
// the required object type (track/V0) and the required charge for first element in pair def,
// and the second passed object does the same w.r. to the second element in pair def.
// Moreover, all cuts are checked and the operation fails if a cut check is unsuccessful.
// Finally, if a true pair is required, this is checked at the end.
//

   AliDebug(AliLog::kDebug + 2, "<-");
   
   // check match with prototype
   // include check on true PID if required
   if (!fDaughterDef->MatchesDaughter(daughter, fOnlyTrue)) return kFALSE;
   
   // if matching is successful
   // update track data member and assigh default mass
   fDaughter = daughter;
   fDaughter->SetMass(fDaughterDef->GetMass());

   // check the cuts
   if (!fCuts.IsSelected(daughter)) return kFALSE;
   
   // if track is accepted increment counter   
   ++fCount;

   return kTRUE;
}

//_____________________________________________________________________________
void AliRsnMonitor::Compute()
{
//
// Virtual method to compute pair quantities of interest
//

   AliWarning("Implement this method in derived classes");
}

//_____________________________________________________________________________
void AliRsnMonitor::Init(const char* /*prefix*/, TList* /*list*/)
{
//
// Virtual method to compute pair quantities of interest
//

   AliWarning("Implement this method in derived classes");
}
