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

//==============================================================================
//   AliTrigEvent - Base class for generic information exchanged by a trigger
//                  device. Trigger inputs and outputs are represented and
//                  handled via AliTrigEvent objects. Trigger events are typically
//                  wrappers for the information exchanged on a single I/O slot
//                  or a group of correlated inputs.
//==============================================================================

#include "AliTrigEvent.h"

#include <TClass.h>
#include <TBits.h>
#include <TROOT.h>


ClassImp(AliTrigEvent)

//______________________________________________________________________________
AliTrigEvent::Activate()
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
  fValue = src->GetValue();
  return kTRUE;
}

ClassImp(AliTrigEventWithObject)

//______________________________________________________________________________
AliTrigEventWithObject::AliTrigEventWithObject(const char *classname)
                       :AliTrigEvent(),
                        fValue(0),
                        fType(0)
{
// Normal constructor where a class name is provided for the embedded object.
// If the event is created in this way one will only be able to connect to 
// events embedding the same object type (via connectors). Otherwise the type
// will be set upon the first call of SetValue.
   fType = gROOT->GetClass(classname);
   if (!fType) Error("ctor", "No class named <%s> available.", classname);
}   

//______________________________________________________________________________
AliTrigEventWithObject::AliTrigEventWithObject(const AliTrigEventWithObject &other)
                       :AliTrigEvent(other),
                        fValue(other.fValue),
                        fType(other.fType)
{
// Copy constructor.   
}

//______________________________________________________________________________
AliTrigEventWithObject::operator=(const AliTrigEventWithObject &other)
{
// Assignment operator.
   if (&other == this) return *this;
   AliTrigEvent::operator=(other);
   fValue = other.fValue;
   fType = other.fType;
}

//______________________________________________________________________________
Bool_t AliTrigEventWithObject::ImportData(AliTrigEvent *source)
{
// Import data from likewise signal.
  AliTrigEventWithObject *src = (AliTrigEventWithObject *)source;
  Bool_t done = SetValue(src->GetValue());
  if (!done) Error("ImportData", "Cannot import object <%s> of class <%s> since event type was set to: <%s>",
                   src->GetValue()->GetName(), src->GetValue()->ClassName(), fType->GetName());
  return done;
}

//______________________________________________________________________________
Bool_t AliTrigEventWithObject::SetType(const char *classname)
{
// Set the type of this event. Can be done only once.
  TClass *type = gROOT->GetClass(classname);
  if (!type) {
    Error("SetType", "Unknown class %s", classname);
    return kFALSE;
  }
  if (!fType) fType = type;
  if (fType != type) {
    Error("SetType", "Type %s not matching the one defined for this event <%s>",
          classname, fType->GetName());
    return kFALSE;
  }
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
  TClass *type = value->IsA();
  // Set the type if used for the first time.
  if (!fType) fType = type;
  // Check consistency of the value with event type.  
  if (type != fType) return kFALSE;
  fValue = value;
  return kTRUE;
}  
