//
// This cut definition works on mini-pairs for the 'mini' subpackage.
// Since cuts on mini-pairs can usually be just kinematic selections
// or kinematic comparisons between 4-momenta of daughters, they are all
// implemented in this class, by the use of an enumeration which allows
// the user to choose between all possibilities.
//

#include <Riostream.h>

#include "AliRsnMiniPair.h"
#include "AliRsnCutMiniPair.h"

ClassImp(AliRsnCutMiniPair)

//__________________________________________________________________________________________________
AliRsnCutMiniPair::AliRsnCutMiniPair(const char *name, EType type) :
   AliRsnCut(name, AliRsnTarget::kMother),
   fType(type)
{
//
// Constructor
//
}

//__________________________________________________________________________________________________
Bool_t AliRsnCutMiniPair::IsSelected(TObject *obj)
{
//
// Global check
//

   AliRsnMiniPair *pair = dynamic_cast<AliRsnMiniPair*>(obj);
   if (!pair) {
      AliError("This cut applies only to mini-pairs");
      return kFALSE;
   }
   
   switch (fType) {
      case kRapidityRange:
         fCutValueD = pair->Y(0);
         return OkRangeD();
      case kRapidityRangeMC:
         fCutValueD = pair->Y(1);
         return OkRangeD();
      case kMomentumComparison:
         AliWarning("TODO: implement this");
         return kTRUE;
      default:
         AliWarning("Undefined enum value");
         return kTRUE;
   }
}
