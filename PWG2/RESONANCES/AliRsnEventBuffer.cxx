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

//
// Class AliRsnEventBuffer
//
// Implements a temporary buffer of many AliRsnEvent objects
// which is useful for event mixing.
//
// author: Martin Vala (Martin.Vala@cern.ch)
// revised by: A. Pulvirenti (alberto.pulvirenti@ct.infn.it)
//

#include "AliLog.h"

#include "AliRsnCut.h"
#include "AliRsnCutSet.h"
#include "AliRsnEventBuffer.h"

ClassImp(AliRsnEventBuffer)

//_____________________________________________________________________________
AliRsnEventBuffer::AliRsnEventBuffer(Int_t buffSize, Bool_t deleteBufferWhenReset) :
    fDeleteBufferWhenReset(deleteBufferWhenReset),
    fEventsBufferSize(buffSize),
    fEventsBufferIndex(-1)
{
//
// Constructor.
// Initializes to NULL all AliRsnEvent pointers.
//

  Int_t i;
  for (i = 0; i < fEventsBufferSize; i++) fEventsBuffer[i] = 0;
}

//_____________________________________________________________________________
AliRsnEventBuffer::~AliRsnEventBuffer()
{
//
// Destructor.
// Actually does nothing.
//

  //ClearBuffer();
}

//_____________________________________________________________________________
void AliRsnEventBuffer::ClearBuffer()
{
//
// Clears buffer, resetting all pointers.
// If an event is on the HEAP, it is deleted.
//
  Int_t i;
  for (i = 0; i < fEventsBufferSize; i++)
  {
    AliDebug(1, Form("Clearing slot #%p in buffer", fEventsBuffer[i]));
    if (fEventsBuffer[i]->IsOnHeap()) delete fEventsBuffer[i];
    fEventsBuffer[i] = 0;
  }
}

//_____________________________________________________________________________
void AliRsnEventBuffer::ResetIndex()
{
//
// Resets the index for accessing events in buffer
//
  fEventsBufferIndex = -1;
  if (fDeleteBufferWhenReset == kTRUE) ClearBuffer();
}

//_____________________________________________________________________________
void AliRsnEventBuffer::AddEvent(AliRsnEvent * event)
{
//
// Insert a new event in the buffer.
// Since the buffer has a fixed size, when the last index is reached,
// the new event replaces the oldest one in the buffer
//

  if (fEventsBufferIndex >= fEventsBufferSize - 1) ResetIndex();
  fEventsBufferIndex++;
  if (fEventsBuffer[fEventsBufferIndex]) {
    //AliInfo("Replacing event");
    *fEventsBuffer[fEventsBufferIndex] = *event;
  }
  else {
    //AliInfo("New event");
    fEventsBuffer[fEventsBufferIndex] = new AliRsnEvent(*event);
  }
}

//_____________________________________________________________________________
Int_t AliRsnEventBuffer::IndexOf(AliRsnEvent * event)
{
//
// Return position of the event
//

  Int_t i;
  for (i = 0; i < fEventsBufferSize; i++) {
    if (event == fEventsBuffer[i]) return i;
  }
  
  return -1;
}

//_____________________________________________________________________________
AliRsnEvent* AliRsnEventBuffer::GetCurrentEvent()
{
//
// Returns the current event in the buffer
//
  return GetEvent(fEventsBufferIndex);
}

//_____________________________________________________________________________
AliRsnEvent * AliRsnEventBuffer::GetNextEvent()
{
//
// Returns next event in the buffer w.r. to current
//
  if (fEventsBufferIndex == fEventsBufferSize - 1) return GetEvent(0);
  else return GetEvent(fEventsBufferIndex + 1);
}

//_____________________________________________________________________________
AliRsnEvent * AliRsnEventBuffer::GetNextGoodEvent
(Int_t &start, AliRsnCutSet *cuts)
{
//
// Scans the buffer starting from 'start' in order to
// find another event which satisfies the cuts defined in 'cuts'.
// If it finds such an event, returns it and upgrades the 'start' value
// to its position, otherwise returns NULL.
//

  Int_t i = start;
  AliRsnEvent *ref = GetCurrentEvent();
  AliRsnEvent *ev = 0x0;
  for(;;i--) {
    ev = GetEvent(i);
    if (ev == ref) continue;
    if (!ev) break;
    if (!cuts) {
      start = i;
      return ev;
    }
    else if (cuts->IsSelected(AliRsnCut::kMixEvent, ref, ev)) {
      start = i;
      return ev;
    }
  }
  
  return 0x0;
}

//_____________________________________________________________________________
AliRsnEvent * AliRsnEventBuffer::GetEvent(Int_t index)
{
//
// Returns an event in the buffer, provided its index
//
  if (index < 0 || index >= fEventsBufferSize)
  {
    AliDebug(AliLog::kError,Form("Index out of range (index = %d , size = %d)", index, fEventsBufferSize));
    return 0x0;
  }
  return (AliRsnEvent*)fEventsBuffer[index];
}

//_____________________________________________________________________________
Int_t AliRsnEventBuffer::NEmptySlots()
{
//
// Tells if the buffer is completely filled or has empty slots
//
  Int_t i, counter = 0;
  for (i = 0; i < fEventsBufferSize; i++)
  {
    if (!fEventsBuffer[i]) counter++;
  }

  return counter;
}
