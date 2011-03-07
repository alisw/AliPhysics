/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

////////////////////////////////////////////////////////////////////////////////
//
//  Base class used wherever it is needed to check the class type of
//  an object w.r. to the RSN framework analysis (daughter, mother, event) 
//  which could be used for cut checking or value computing.
//  Since most of these operation are implemented into classes that
//  operate on any of such objects, then this class helps in making sure
//  that the object being processed corresponds to what is expected.
//  It also contains three pointers to which any passed object is casted
//  in order to have a quick reference to any allowed object type from
//  an appropriate pointer, which is propagated to all inheriting classes.
//
//  authors: A. Pulvirenti (alberto.pulvirenti@ct.infn.it)
//           M. Vala (martin.vala@cern.ch)
//
////////////////////////////////////////////////////////////////////////////////

#include "AliLog.h"

#include "AliRsnDaughter.h"
#include "AliRsnMother.h"

#include "AliRsnTarget.h"

ClassImp(AliRsnTarget)

//_____________________________________________________________________________
Bool_t AliRsnTarget::TargetOK(TObject *object)
{
//
// This method doew the following things:
// 1) check if the object class matches the required target type
// 2) if (1) is successful, set the built-in pointer data member
//    in order to point to it, after being casted accordingly
// 3) if the target is a daughter or a mother, adjust the pointer 
//    to reference event accordingly
//

   // reset local pointers and then initialize
   // only the right one by static cast, if found
   fDaughter = 0x0;
   fMother   = 0x0;
   fEvent    = 0x0;

   // fails by default if a NULL pointer is passed
   if (!object) {
      AliError("Passed a NULL object");
      return kFALSE;
   }
   
   // checks matching between argument type and expected object type
   // the passed object is also casted into the appropriate pointer
   // when a daughter or mother is passed and the target type expects
   // and event, the check is yes successful and the owner event is used
   if (object->IsA() == AliRsnDaughter::Class() && (fTargetType == kDaughter || fTargetType == kEvent)) {
      fDaughter = static_cast<AliRsnDaughter*>(object);
      fEvent = fDaughter->GetOwnerEvent();
      return kTRUE;
   }
   else if (object->IsA() == AliRsnMother::Class() && (fTargetType == kMother || fTargetType == kEvent)) {
      fMother = static_cast<AliRsnMother*>(object);
      fEvent = fMother->GetRefEvent();
      return kTRUE;
   }
   else if (object->IsA() == AliRsnEvent::Class() && fTargetType == kEvent) {
      fEvent = static_cast<AliRsnEvent*>(object);
      return kTRUE;
   }
   else {
      if (object) {
         AliError(Form("[%s] Target mismatch: expected '%s', passed '%s'", GetName(), GetTargetTypeName(), object->ClassName()));
      } else {
         AliError(Form("[%s] Target mismatch: expected '%s', passed 'NULL'", GetName(), GetTargetTypeName()));
      }
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
