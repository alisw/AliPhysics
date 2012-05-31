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

/* $Id$ */
// Author: Andrei Gheata, 28/07/2009

#include "AliTrigEvent.h"

#include <TClass.h>
#include <TBits.h>
#include <TROOT.h>

ClassImp(AliTrigEvent)

//==============================================================================
//   AliTrigEvent - Base class for generic information exchanged by a trigger
//                  device. Trigger inputs and outputs are represented and
//                  handled via AliTrigEvent objects. Trigger events are typically
//                  wrappers for the information exchanged on a single I/O slot
//                  or a group of correlated inputs.
//==============================================================================

//______________________________________________________________________________
AliTrigEvent &AliTrigEvent::operator=(const AliTrigEvent &other)
{
// Assignment operator.
   if (&other == this) return *this;
   TNamed::operator=(other);
   return *this;
}
//______________________________________________________________________________
void AliTrigEvent::Activate(Bool_t flag)
{
// Activate/deactivate this signal.
  TObject::SetBit(kActive, flag);
}                        

ClassImp(AliTrigEventWithMask)

//______________________________________________________________________________
Bool_t AliTrigEventWithMask::ImportData(AliTrigEvent *source)
{
// Import data from likewise signal.
  AliTrigEventWithMask *src = (AliTrigEventWithMask *)source;
  SetValue(src->GetValue());
  return kTRUE;
}

//______________________________________________________________________________
void AliTrigEventWithMask::SetValue(TBits *value)
{
// Set the mask value.
  *fValue = *value;
}

//______________________________________________________________________________
AliTrigEventWithMask::AliTrigEventWithMask(const AliTrigEventWithMask &other)
                       :AliTrigEvent(other),
                        fValue(NULL)
{
// Copy constructor.   
  *fValue = *other.fValue;
}

//______________________________________________________________________________
AliTrigEventWithMask &AliTrigEventWithMask::operator=(const AliTrigEventWithMask &other)
{
// Assignment operator.
   if (&other == this) return *this;
   AliTrigEvent::operator=(other);
   *fValue = *other.fValue;
   return *this;
}

ClassImp(AliTrigEventWithObject)

//______________________________________________________________________________
AliTrigEventWithObject::AliTrigEventWithObject(const char *name,const char *classname)
                       :AliTrigEvent(name),
                        fValue(0),
                        fType("")
{
// Normal constructor where a class name is provided for the embedded object.
// If the event is created in this way one will only be able to connect to 
// events embedding the same object type (via connectors). If empty string the type
// will be set upon the first call of SetValue.
  SetType(classname);
}   

//______________________________________________________________________________
AliTrigEventWithObject::AliTrigEventWithObject(const AliTrigEventWithObject &other)
                       :AliTrigEvent(other),
                        fValue(NULL),
                        fType(other.fType)
{
// Copy constructor.   
  TClass* pClass=TClass::GetClass(fType);
  if (!pClass) return;
  fValue = (TObject*)pClass->New();
  fValue->Copy(*other.fValue);
}

//______________________________________________________________________________
AliTrigEventWithObject &AliTrigEventWithObject::operator=(const AliTrigEventWithObject &other)
{
// Assignment operator.
  if (&other == this) return *this;
  AliTrigEvent::operator=(other);
  fType = other.fType;
  TClass* pClass=TClass::GetClass(fType);
  if (!pClass) return *this;
  fValue = (TObject*)pClass->New();
  fValue->Copy(*other.fValue);
  return *this;
}

//______________________________________________________________________________
Bool_t AliTrigEventWithObject::ImportData(AliTrigEvent *source)
{
// Import data from likewise signal.
  AliTrigEventWithObject *src = (AliTrigEventWithObject *)source;
  Bool_t done = SetValue(src->GetValue());
  if (!done) Error("ImportData", "Cannot import object <%s> of class <%s> since event type was set to: <%s>",
                   src->GetValue()->GetName(), src->GetValue()->ClassName(), fType.Data());
  return done;
}

//______________________________________________________________________________
Bool_t AliTrigEventWithObject::SetType(const char *classname)
{
// Set the type of this event. Can be done only once.
  if (!strlen(classname)) return kFALSE;
  if (!fType.IsNull()) {
    Error("SetType", "Type for this trigger event already set to: %s", fType.Data());
    return kFALSE;
  }  
  TClass *type = gROOT->GetClass(classname);
  if (!type) {
    Error("SetType", "Unknown class %s", classname);
    return kFALSE;
  }
  fType = classname;
  return kTRUE;
}

//______________________________________________________________________________
Bool_t AliTrigEventWithObject::SetValue(TObject *value)
{
// Set the current event content. Checks consistency with event type.
  if (!value) {
    // Reset current value.
    fValue = NULL;
    return kTRUE;
  }
  // Set the type if used for the first time.
  if (!fType) fType = value->ClassName();
  fValue = value;
  return kTRUE;
}  
