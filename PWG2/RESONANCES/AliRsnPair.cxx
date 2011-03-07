//
// *** Class AliRsnPair ***
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
#include "AliRsnCutSet.h"

#include "AliRsnPair.h"

ClassImp(AliRsnPair)

//_____________________________________________________________________________
AliRsnPair::AliRsnPair(const char *name, AliRsnPairDef *def) :
   TNamed(name, ""),
   fOnlyTrue(kFALSE),
   fCheckDecay(kFALSE),
   fIsMixed(kFALSE),
   fCount(0),
   fPairDef(def),
   fCutManager(Form("cutMgr_%s", name)),
   fMother()
{
//
// Default constructor
//
}

//_____________________________________________________________________________
AliRsnPair::AliRsnPair(const AliRsnPair& copy) :
   TNamed(copy),
   fOnlyTrue(copy.fOnlyTrue),
   fCheckDecay(copy.fCheckDecay),
   fIsMixed(copy.fIsMixed),
   fCount(copy.fCount),
   fPairDef(copy.fPairDef),
   fCutManager(copy.fCutManager),
   fMother(copy.fMother)
{
//
// Default constructor
//
}

//_____________________________________________________________________________
AliRsnPair& AliRsnPair::operator=(const AliRsnPair& copy)
{
   fOnlyTrue = copy.fOnlyTrue;
   fCheckDecay = copy.fCheckDecay;
   fIsMixed = copy.fIsMixed;
   fCount = copy.fCount;
   fPairDef = copy.fPairDef;
   fMother = copy.fMother;
   fCutManager = copy.fCutManager;

   return (*this);
}

//_____________________________________________________________________________
AliRsnPair::~AliRsnPair()
{
//
// Destructor
//
}

//_____________________________________________________________________________
void AliRsnPair::Print(Option_t* /*option*/) const
{
//
// Prints info about pair
//
}

//_____________________________________________________________________________
Bool_t AliRsnPair::Fill
(AliRsnDaughter *daughter1, AliRsnDaughter *daughter2, Bool_t refFirst)
{
//
// Checks that first argument matches definitions for first daughter
// and the same for second argument, where the order is defined by
// the AliRsnPairDef data member.
// If the matching is successful, the AliRsnMother data member is 
// initialized using the mass hypotheses defined here and the momenta
// in the passed daughters.
// The third argument is necessary to choose which one of the possible two
// events owning the two daughter will be used as reference.
//
   
   // check matching and exit if one of them fails
   // if true pair is required, this is taken into account:
   // if both true pairs and correct decay tree is required,
   // then we must be sure that also the true PID of daughters matches,
   // instead if correct decay tree is not required this additional check is not done
   if (!fPairDef->GetDef1()->MatchesDaughter(daughter1, fOnlyTrue && fCheckDecay)) return kFALSE;
   if (!fPairDef->GetDef2()->MatchesDaughter(daughter2, fOnlyTrue && fCheckDecay)) return kFALSE;
   
   // if matching is successful
   // compute 4-momenta of daughters and mother
   fMother.SetDaughter(0, daughter1);
   fMother.SetDaughter(1, daughter2);
   fMother.ComputeSum(fPairDef->GetMass1(), fPairDef->GetMass2());
   
   // if required a true pair, check this here and eventually return a fail message
   // this is done using the method AliRsnMother::CommonMother with 2 arguments
   // passed by reference, where the real GEANT label of the particle is stored
   // and one can check if these tracks are both really secondaries (ID >= 0)
   if (fOnlyTrue) {
      Int_t m0, m1, common;
      common = fMother.CommonMother(m0, m1);
      if (m0 < 0 || m1 < 0) return kFALSE;
      if (common != fPairDef->GetMotherPDG()) return kFALSE;
   }
   
   // point to first event as reference
   // and checks the pair cuts,
   // (done first because it is more likely 
   // that it is not passed and execution is faster)
   if (!fCutManager.PassMotherCuts(&fMother)) return kFALSE;

   // cuts on track #1 & common
   if (!fCutManager.PassDaughter1Cuts(daughter1)) {
      AliDebug(AliLog::kDebug + 2, "Specific cuts for track #1 not passed");
      return kFALSE;
   }
   if (!fCutManager.PassCommonDaughterCuts(daughter1)) {
      AliDebug(AliLog::kDebug + 2, "Common cuts for track #1 not passed");
      return kFALSE;
   }

   // cuts on track #2 & common
   if (!fCutManager.PassDaughter2Cuts(daughter2)) {
      AliDebug(AliLog::kDebug + 2, "Specific cuts for track #2 not passed");
      return kFALSE;
   }
   if (!fCutManager.PassCommonDaughterCuts(daughter2)) {
      AliDebug(AliLog::kDebug + 2, "Common cuts for track #2 not passed");
      return kFALSE;
   }

   // if pair is accepted, increment counter
   ++fCount;
   
   // assign reference event
   if (refFirst) fMother.SetRefEvent(daughter1->GetOwnerEvent()); else fMother.SetRefEvent(daughter2->GetOwnerEvent());

   return kTRUE;
}

//_____________________________________________________________________________
void AliRsnPair::Compute()
{
//
// Virtual method to compute pair quantities of interest
//

   AliWarning("Implement this method in derived classes");
}

//_____________________________________________________________________________
void AliRsnPair::Init(const char* /*prefix*/, TList* /*list*/)
{
//
// Virtual method to compute pair quantities of interest
//

   AliWarning("Implement this method in derived classes");
}
