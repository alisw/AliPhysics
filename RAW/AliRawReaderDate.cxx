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

///////////////////////////////////////////////////////////////////////////////
//
// This is a class for reading a raw data from a date event and providing
// information about digits
//
///////////////////////////////////////////////////////////////////////////////

#include "AliRawReaderDate.h"
#ifdef ALI_DATE
#include "event.h"
#endif

ClassImp(AliRawReaderDate)


AliRawReaderDate::AliRawReaderDate(
#ifdef ALI_DATE
				   void* event
#else
				   void* /* event */
#endif
)
{
// create an object to read digits from the given date event

#ifdef ALI_DATE
  fEvent = (eventHeaderStruct*) event;
  fSubEvent = NULL;
  fPosition = fEnd = NULL;
#else
  Fatal("AliRawReaderDate", "this class was compiled without DATE");
#endif
}

AliRawReaderDate::AliRawReaderDate(const AliRawReaderDate& rawReader) :
  AliRawReader(rawReader)
{
// copy constructor

  fEvent = rawReader.fEvent;
  fSubEvent = rawReader.fSubEvent;
  fPosition = rawReader.fPosition;
  fEnd = rawReader.fEnd;
}

AliRawReaderDate& AliRawReaderDate::operator = (const AliRawReaderDate& 
						rawReader)
{
// assignment operator

  this->~AliRawReaderDate();
  new(this) AliRawReaderDate(rawReader);
  return *this;
}


UInt_t AliRawReaderDate::GetType() const
{
// get the type from the event header

#ifdef ALI_DATE
  if (!fEvent) return 0;
  return fEvent->eventType;
#else
  return 0;
#endif
}

UInt_t AliRawReaderDate::GetRunNumber() const
{
// get the run number from the event header

#ifdef ALI_DATE
  if (!fEvent) return 0;
  return fEvent->eventRunNb;
#else
  return 0;
#endif
}

const UInt_t* AliRawReaderDate::GetEventId() const
{
// get the event id from the event header

#ifdef ALI_DATE
  if (!fEvent) return NULL;
  return fEvent->eventId;
#else
  return NULL;
#endif
}

const UInt_t* AliRawReaderDate::GetTriggerPattern() const
{
// get the trigger pattern from the event header

#ifdef ALI_DATE
  if (!fEvent) return NULL;
  return fEvent->eventTriggerPattern;
#else
  return NULL;
#endif
}

const UInt_t* AliRawReaderDate::GetDetectorPattern() const
{
// get the detector pattern from the event header

#ifdef ALI_DATE
  if (!fEvent) return NULL;
  return fEvent->eventDetectorPattern;
#else
  return NULL;
#endif
}

const UInt_t* AliRawReaderDate::GetAttributes() const
{
// get the type attributes from the event header

#ifdef ALI_DATE
  if (!fEvent) return NULL;
  return fEvent->eventTypeAttribute;
#else
  return NULL;
#endif
}

UInt_t AliRawReaderDate::GetGDCId() const
{
// get the GDC Id from the event header

#ifdef ALI_DATE
  if (!fEvent) return 0;
  return fEvent->eventGdcId;
#else
  return 0;
#endif
}


Bool_t AliRawReaderDate::ReadMiniHeader()
{
// read a mini header at the current position
// returns kFALSE if the mini header could not be read

#ifdef ALI_DATE
  if (!fEvent) return kFALSE;
  do {
    if (fCount > 0) fPosition += fCount;      // skip payload if event was not selected
    if (!fSubEvent || (fPosition >= fEnd)) {  // new sub event
      if (fPosition >= ((UChar_t*)fEvent)+fEvent->eventSize) return kFALSE;  // end of data
      if (fSubEvent) {
	fSubEvent = (eventHeaderStruct*) (((UChar_t*)fSubEvent) + 
					  fSubEvent->eventSize);
      } else {
	if (fEvent->eventSize <= fEvent->eventHeadSize) return kFALSE; // no sub events
	fSubEvent = (eventHeaderStruct*) (((UChar_t*)fEvent) + 
					  fEvent->eventHeadSize);
      }
      fCount = 0;
      fPosition = ((UChar_t*)fSubEvent) + fSubEvent->eventHeadSize + 
	sizeof(equipmentHeaderStruct);
      fEnd = ((UChar_t*)fSubEvent) + fSubEvent->eventSize;
    }
    if (fPosition >= fEnd) continue;          // no data left in the payload
    if (fPosition + sizeof(AliMiniHeader) > fEnd) {
      Error("ReadMiniHeader", "could not read data!");
      return kFALSE;
    }
    fMiniHeader = (AliMiniHeader*) fPosition;
    fPosition += sizeof(AliMiniHeader);
    if (!CheckMiniHeader()) {
      Warning("ReadMiniHeader", "skipping %d bytes\n"
	      " run: %d  event: %d %d  LDC: %d  GDC: %d\n", 
	      fEnd - fPosition, fSubEvent->eventRunNb, 
	      fSubEvent->eventId[0], fSubEvent->eventId[1],
	      fSubEvent->eventLdcId, fSubEvent->eventGdcId);
      fCount = 0;
      fPosition = fEnd;
      continue;
    }
    fCount = fMiniHeader->fSize;
    if (fPosition + fCount > fEnd) {  // check data size in mini header and sub event
      Error("ReadMiniHeader", "size in mini header exceeds event size!");
      Warning("ReadMiniHeader", "skipping %d bytes\n"
	      " run: %d  event: %d %d  LDC: %d  GDC: %d\n", 
	      fEnd - fPosition, fSubEvent->eventRunNb, 
	      fSubEvent->eventId[0], fSubEvent->eventId[1],
	      fSubEvent->eventLdcId, fSubEvent->eventGdcId);
      fCount = 0;
      fPosition = fEnd;
      continue;
    }
  } while (!IsSelected());
  return kTRUE;
#else
  return kFALSE;
#endif
}

Bool_t AliRawReaderDate::ReadNextData(UChar_t*& data)
{
// reads the next payload at the current position
// returns kFALSE if the data could not be read

  while (fCount == 0) {
    if (!ReadMiniHeader()) return kFALSE;
  }
  data = fPosition;
  fPosition += fCount;  
  fCount = 0;
  return kTRUE;
}

Bool_t AliRawReaderDate::ReadNext(UChar_t* data, Int_t size)
{
// reads the next block of data at the current position
// returns kFALSE if the data could not be read

  if (fPosition + size > fEnd) {
    Error("ReadNext", "could not read data!");
    return kFALSE;
  }
  memcpy(data, fPosition, size);
  fPosition += size;
  fCount -= size;
  return kTRUE;
}


Bool_t AliRawReaderDate::Reset()
{
// reset the current position to the beginning of the event

#ifdef ALI_DATE
  fSubEvent = NULL;
#endif
  fCount = 0;
  fPosition = fEnd = NULL;
  return kTRUE;
}
