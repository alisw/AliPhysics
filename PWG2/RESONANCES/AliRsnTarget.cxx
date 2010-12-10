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

#include "AliLog.h"

#include "AliRsnEvent.h"
#include "AliRsnDaughter.h"
#include "AliRsnMother.h"

#include "AliRsnTarget.h"

ClassImp(AliRsnTarget)

AliRsnEvent* AliRsnTarget::fgCurrentEvent = 0x0;

//_____________________________________________________________________________
Bool_t AliRsnTarget::TargetOK(TObject *object)
{
//
// This method compares the target type stored as data member
// with the type of the object passed as argument, and returns
// kTRUE or kFALSE depending if they match or not.
//

  // fails by default if a NULL pointer is passed
  if (!object)
  {
    AliError("Object is NULL");
    return kFALSE;
  }

  // checks if the object is correct by dynamic casting
  switch (fTargetType)
  {
    case kDaughter:
      if (dynamic_cast<AliRsnDaughter*>(object) == 0x0)
      {
        AliError(Form("[%s] Target mismatch: expected 'AliRsnDaughter', passed '%s'", GetName(), object->ClassName()));
        Print();
        return kFALSE;
      }
      break;
    case kMother:
      if (dynamic_cast<AliRsnMother*>(object) == 0x0)
      {
        AliError(Form("[%s] Target mismatch: expected 'AliRsnMother', passed '%s'", GetName(), object->ClassName()));
        Print();
        return kFALSE;
      }
      break;
    case kEvent:
      if (dynamic_cast<AliRsnEvent*>(object) == 0x0)
      {
        AliError(Form("[%s] Target mismatch: expected 'AliRsnEvent', passed '%s'", GetName(), object->ClassName()));
        Print();
        return kFALSE;
      }
      break;
    default:
      return kTRUE;
  }
  
  return kTRUE;
}

//______________________________________________________________________________
Char_t AliRsnTarget::GetTargetTypeChar() const
{
//
// Returns a single character identifying the cut target type.
//

  switch (fTargetType) 
  {
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

  switch (fTargetType) 
  {
    case kDaughter: return "Daughter";
    case kMother: return "Mother";
    case kEvent: return "Event";
    default: return "Undefined";
  }
}
