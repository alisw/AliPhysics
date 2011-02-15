//
// Class AliRsnCutMomentumComparison
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
#include "TMath.h"

#include "AliRsnDaughter.h"
#include "AliRsnMother.h"
#include "AliRsnCutMomentumComparison.h"

ClassImp(AliRsnCutMomentumComparison)

//_________________________________________________________________________________________________
AliRsnCutMomentumComparison::AliRsnCutMomentumComparison(const char *name, EMode mode) :
   AliRsnCut(name, AliRsnCut::kMother),
   fMode(mode)
{
//
// Default constructor.
//
}

//_________________________________________________________________________________________________
AliRsnCutMomentumComparison::AliRsnCutMomentumComparison(const AliRsnCutMomentumComparison& copy) :
   AliRsnCut(copy),
   fMode(copy.fMode)
{
//
// Copy constructor
//
}

//_________________________________________________________________________________________________
AliRsnCutMomentumComparison& AliRsnCutMomentumComparison::operator=(const AliRsnCutMomentumComparison& copy)
{
//
// Assignment operator
//

   AliRsnCut::operator=(copy);
   fMode = copy.fMode;
}

//_________________________________________________________________________________________________
Bool_t AliRsnCutMomentumComparison::IsSelected(TObject *object)
{
//
// Cut checker.
//

   // convert the object into the unique correct type
   if (!TargetOK(object)) {
      AliError(Form("[%s]: this cut works only with AliRsnMother objects", GetName()));
      return kTRUE;
   }

   // compare momenta
   AliRsnMother *mother = dynamic_cast<AliRsnMother*>(object);
   Double_t p1  = mother->GetDaughter(0)->GetRef()->P();
   Double_t p2  = mother->GetDaughter(1)->GetRef()->P();
   Double_t pt1 = mother->GetDaughter(0)->GetRef()->Pt();
   Double_t pt2 = mother->GetDaughter(1)->GetRef()->Pt();
   
   switch (fMode)
   {
      case kFirstLargerP  : return (p1 > p2);
      case kFirstSmallerP : return (p1 < p2);
      case kFirstLargerPt : return (pt1 > pt2);
      case kFirstSmallerPt: return (pt1 < pt2);
      default:
         AliError("Invalid mode selected. Cut is skipped");
         return kTRUE;
   }
}
