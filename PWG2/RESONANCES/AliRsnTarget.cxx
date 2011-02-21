//
// *** Class AliRsnTarget ***
//
// Base class used wherever it is needed to check the class type of
// an object (daughter, mother, event) which could be used for
// cut checking or value computing.
// Since most of these operation are implemented into classes that
// operate on any of such objects, then this class helps in making sure
// that the object being processed corresponds to what is expected.
//
// authors: Alberto Pulvirenti (alberto.pulvirenti@ct.infn.it)
//          Martin Vala (martin.vala@cern.ch)
//

#include "AliLog.h"

#include "AliRsnDaughter.h"
#include "AliRsnMother.h"

#include "AliRsnTarget.h"

ClassImp(AliRsnTarget)

AliRsnEvent*   AliRsnTarget::fgCurrentEvent = 0x0;
const Double_t AliRsnTarget::fgkVeryBig     = 1E+10;
const Double_t AliRsnTarget::fgkVerySmall   = 1E-10;

//_____________________________________________________________________________
Bool_t AliRsnTarget::TargetOK(TObject *object)
{
//
// This method doew two things:
// 1) check if the object class matches the required target type
// 2) if (1) is successful, set the built-in pointer data member
//    in order to point to it, after being casted accordingly
//

   // fails by default if a NULL pointer is passed
   if (!object) return kFALSE;
   
   // reset local pointers and then initialize
   // only the right one by static cast, if found
   fDaughter = 0x0;
   fMother   = 0x0;
   fEvent    = 0x0;
   if (object->IsA() == AliRsnDaughter::Class() && fTargetType == kDaughter) {
      fDaughter = static_cast<AliRsnDaughter*>(object);
      return kTRUE;
   }
   else if (object->IsA() == AliRsnMother::Class() && fTargetType == kMother) {
      fMother = static_cast<AliRsnMother*>(object);
      return kTRUE;
   }
   else if (object->IsA() == AliRsnEvent::Class() && fTargetType == kEvent) {
      fEvent = static_cast<AliRsnEvent*>(object);
      return kTRUE;
   }
   else {
      AliError(Form("[%s] Target mismatch: expected '%s', passed '%s'", GetName(), GetTargetTypeName(), object->ClassName()));
      return kFALSE;
   }
}

//______________________________________________________________________________
Char_t AliRsnTarget::GetTargetTypeChar() const
{
//
// Returns a single character identifying the cut target type.
//

   switch (fTargetType) {
      case kDaughter: return 'D';
      case kMother: return 'M';
      case kEvent: return 'E';
      default: return 'X';
   }
}

//______________________________________________________________________________
const char* AliRsnTarget::GetTargetTypeName() const
{
//
// Returns a string with the name of the cut target type-
//

   switch (fTargetType) {
      case kDaughter: return "Daughter";
      case kMother: return "Mother";
      case kEvent: return "Event";
      default: return "Undefined";
   }
}
