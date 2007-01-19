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

///////////////////////////////////////////////////////////////////////////////
///
/// This is a class for reading raw data from a date file or event.
///
/// The AliRawReaderDate is constructed either with a pointer to a
/// date event or with a file name and an event number.
///
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
				   ) :
  fFile(NULL),
  fEvent(NULL),
  fSubEvent(NULL),
  fEquipment(NULL),
  fPosition(NULL),
  fEnd(NULL)
{
// create an object to read digits from the given date event

#ifdef ALI_DATE
  fEvent = (eventHeaderStruct*) event;
#else
  Fatal("AliRawReaderDate", "this class was compiled without DATE");
#endif
}

AliRawReaderDate::AliRawReaderDate(
#ifdef ALI_DATE
				   const char* fileName, Int_t eventNumber
#else
				   const char* /*fileName*/, 
				   Int_t /*eventNumber*/
#endif
				   ) :
  fFile(NULL),
  fEvent(NULL),
  fSubEvent(NULL),
  fEquipment(NULL),
  fPosition(NULL),
  fEnd(NULL)
{
// create an object to read digits from the given date event

#ifdef ALI_DATE
  fFile = fopen(fileName, "rb");
  if (!fFile) {
    Error("AliRawReaderDate", "could not open file %s", fileName);
    return;
  }
  if (eventNumber < 0) return;

  eventHeaderStruct header;
  UInt_t headerSize = sizeof(eventHeaderStruct);
  while (fread(&header, 1, headerSize, fFile) == headerSize) {
    if (eventNumber == 0) {
      UChar_t* buffer = new UChar_t[header.eventSize];
      fseek(fFile, -headerSize, SEEK_CUR);
      if (fread(buffer, 1, header.eventSize, fFile) != header.eventSize) break;
      fEvent = (eventHeaderStruct*) buffer;
      break;
    }
    fseek(fFile, header.eventSize-headerSize, SEEK_CUR);
    eventNumber--;
  }

#else
  Fatal("AliRawReaderDate", "this class was compiled without DATE");
#endif
}

AliRawReaderDate::AliRawReaderDate(const AliRawReaderDate& rawReader) :
  AliRawReader(rawReader),
  fFile(rawReader.fFile),
  fEvent(rawReader.fEvent),
  fSubEvent(rawReader.fSubEvent),
  fEquipment(rawReader.fEquipment),
  fPosition(rawReader.fPosition),
  fEnd(rawReader.fEnd)

{
// copy constructor

  Fatal("AliRawReaderDate", "copy constructor not implemented");
}

AliRawReaderDate& AliRawReaderDate::operator = (const AliRawReaderDate& 
						/*rawReader*/)
{
// assignment operator

  Fatal("operator =", "assignment operator not implemented");
  return *this;
}

AliRawReaderDate::~AliRawReaderDate()
{
// destructor

#ifdef ALI_DATE
  if (fFile) {
    delete[] fEvent;
    fclose(fFile);
  }
#endif
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

const UInt_t* AliRawReaderDate::GetSubEventAttributes() const
{
// get the type attributes from the sub event header

#ifdef ALI_DATE
  if (!fSubEvent) return NULL;
  return fSubEvent->eventTypeAttribute;
#else
  return NULL;
#endif
}

UInt_t AliRawReaderDate::GetLDCId() const
{
// get the LDC Id from the event header

#ifdef ALI_DATE
  if (!fSubEvent) return 0;
  return fSubEvent->eventLdcId;
#else
  return 0;
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


Int_t AliRawReaderDate::GetEquipmentSize() const
{
// get the size of the equipment (including the header)

#ifdef ALI_DATE
  if (!fEquipment) return 0;
  if (fSubEvent->eventVersion <= 0x00030001) {
    return fEquipment->equipmentSize + sizeof(equipmentHeaderStruct);
  } else {
    return fEquipment->equipmentSize;
  }
#else
  return 0;
#endif
}

Int_t AliRawReaderDate::GetEquipmentType() const
{
// get the type from the equipment header

#ifdef ALI_DATE
  if (!fEquipment) return -1;
  return fEquipment->equipmentType;
#else
  return 0;
#endif
}

Int_t AliRawReaderDate::GetEquipmentId() const
{
// get the ID from the equipment header

#ifdef ALI_DATE
  if (!fEquipment) return -1;
  return fEquipment->equipmentId;
#else
  return 0;
#endif
}

const UInt_t* AliRawReaderDate::GetEquipmentAttributes() const
{
// get the attributes from the equipment header

#ifdef ALI_DATE
  if (!fEquipment) return NULL;
  return fEquipment->equipmentTypeAttribute;
#else
  return 0;
#endif
}

Int_t AliRawReaderDate::GetEquipmentElementSize() const
{
// get the basic element size from the equipment header

#ifdef ALI_DATE
  if (!fEquipment) return 0;
  return fEquipment->equipmentBasicElementSize;
#else
  return 0;
#endif
}

Int_t AliRawReaderDate::GetEquipmentHeaderSize() const
{
  // Get the equipment header size
  // 28 bytes by default
#ifdef ALI_DATE
  return sizeof(equipmentHeaderStruct);
#else
  return 0;
#endif
}

Bool_t AliRawReaderDate::ReadHeader()
{
// read a data header at the current position
// returns kFALSE if the data header could not be read

  fErrorCode = 0;

#ifdef ALI_DATE
  fHeader = NULL;
  if (!fEvent) return kFALSE;
  // check whether there are sub events
  if (fEvent->eventSize <= fEvent->eventHeadSize) return kFALSE;

  do {
    // skip payload (if event was not selected)
    if (fCount > 0) fPosition += fCount;

    // get the first or the next equipment if at the end of an equipment
    if (!fEquipment || (fPosition >= fEnd)) {
      fEquipment = NULL;

      // get the first or the next sub event if at the end of a sub event
      if (!fSubEvent || 
	  (fPosition >= ((UChar_t*)fSubEvent) + fSubEvent->eventSize)) {

	// check for end of event data
	if (fPosition >= ((UChar_t*)fEvent)+fEvent->eventSize) return kFALSE;
        if (!TEST_SYSTEM_ATTRIBUTE(fEvent->eventTypeAttribute, 
                                   ATTR_SUPER_EVENT)) {
	  fSubEvent = fEvent;   // no super event
	} else if (fSubEvent) {
	  fSubEvent = (eventHeaderStruct*) (((UChar_t*)fSubEvent) + 
					    fSubEvent->eventSize);
	} else {
	  fSubEvent = (eventHeaderStruct*) (((UChar_t*)fEvent) + 
					    fEvent->eventHeadSize);
	}

	// check the magic word of the sub event
	if (fSubEvent->eventMagic != EVENT_MAGIC_NUMBER) {
	  Error("ReadHeader", "wrong magic number in sub event!\n"
		" run: %d  event: %d %d  LDC: %d  GDC: %d\n", 
		fSubEvent->eventRunNb, 
		fSubEvent->eventId[0], fSubEvent->eventId[1],
		fSubEvent->eventLdcId, fSubEvent->eventGdcId);
	  fErrorCode = kErrMagic;
	  return kFALSE;
	}

	// continue if no data in the subevent
	if (fSubEvent->eventSize == fSubEvent->eventHeadSize) {
	  fPosition = fEnd = ((UChar_t*)fSubEvent) + fSubEvent->eventSize;
	  fCount = 0;
	  continue;
	}

	fEquipment = (equipmentHeaderStruct*)
	  (((UChar_t*)fSubEvent) + fSubEvent->eventHeadSize);

      } else {
	fEquipment = (equipmentHeaderStruct*) fEnd;
      }

      fCount = 0;
      fPosition = ((UChar_t*)fEquipment) + sizeof(equipmentHeaderStruct);
      if (fSubEvent->eventVersion <= 0x00030001) {
        fEnd = fPosition + fEquipment->equipmentSize;
      } else {
        fEnd = ((UChar_t*)fEquipment) + fEquipment->equipmentSize;
      }
    }

    // continue with the next sub event if no data left in the payload
    if (fPosition >= fEnd) continue;

    if (fRequireHeader) {
      // check that there are enough bytes left for the data header
      if (fPosition + sizeof(AliRawDataHeader) > fEnd) {
	Error("ReadHeader", "could not read data header data!");
	Warning("ReadHeader", "skipping %d bytes\n"
		" run: %d  event: %d %d  LDC: %d  GDC: %d\n", 
		fEnd - fPosition, fSubEvent->eventRunNb, 
		fSubEvent->eventId[0], fSubEvent->eventId[1],
		fSubEvent->eventLdcId, fSubEvent->eventGdcId);
	fCount = 0;
	fPosition = fEnd;
	fErrorCode = kErrNoDataHeader;
	continue;
      }

      // "read" the data header
      fHeader = (AliRawDataHeader*) fPosition;
      if ((fPosition + fHeader->fSize) != fEnd) {
	Warning("ReadHeader",
		"raw data size found in the header is wrong (%d != %d)! Using the equipment size instead !",
		fHeader->fSize, fEnd - fPosition);
	fHeader->fSize = fEnd - fPosition;
      }
      fPosition += sizeof(AliRawDataHeader);
    }

    if (fHeader && (fHeader->fSize != 0xFFFFFFFF)) {
      fCount = fHeader->fSize - sizeof(AliRawDataHeader);

      // check consistency of data size in the header and in the sub event
      if (fPosition + fCount > fEnd) {
	Error("ReadHeader", "size in data header exceeds event size!");
	Warning("ReadHeader", "skipping %d bytes\n"
		" run: %d  event: %d %d  LDC: %d  GDC: %d\n", 
		fEnd - fPosition, fSubEvent->eventRunNb, 
		fSubEvent->eventId[0], fSubEvent->eventId[1],
		fSubEvent->eventLdcId, fSubEvent->eventGdcId);
	fCount = 0;
	fPosition = fEnd;
	fErrorCode = kErrSize;
	continue;
      }

    } else {
      fCount = fEnd - fPosition;
    }

  } while (!fEquipment || !IsSelected());

  return kTRUE;
#else
  return kFALSE;
#endif
}

Bool_t AliRawReaderDate::ReadNextData(UChar_t*& data)
{
// reads the next payload at the current position
// returns kFALSE if the data could not be read

  fErrorCode = 0;
  while (fCount == 0) {
    if (!ReadHeader()) return kFALSE;
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

  fErrorCode = 0;
  if (fPosition + size > fEnd) {
    Error("ReadNext", "could not read data!");
    fErrorCode = kErrOutOfBounds;
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
  fEquipment = NULL;
#endif
  fCount = 0;
  fPosition = fEnd = NULL;
  return kTRUE;
}


Bool_t AliRawReaderDate::NextEvent()
{
// go to the next event in the date file

#ifdef ALI_DATE
  if (!fFile) return kFALSE;

  Reset();
  eventHeaderStruct header;
  UInt_t headerSize = sizeof(eventHeaderStruct);
  if (fEvent) delete[] fEvent;
  fEvent = &header;

  while (fread(&header, 1, headerSize, fFile) == headerSize) {
    if (!IsEventSelected()) {
      fseek(fFile, header.eventSize-headerSize, SEEK_CUR);
      continue;
    }
    UChar_t* buffer = new UChar_t[header.eventSize];
    fseek(fFile, -headerSize, SEEK_CUR);
    if (fread(buffer, 1, header.eventSize, fFile) != header.eventSize) {
      Error("NextEvent", "could not read event from file");
      delete[] buffer;
      break;
    }
    fEvent = (eventHeaderStruct*) buffer;
    fEventNumber++;
    return kTRUE;
  };

  fEvent = NULL;
#endif

  return kFALSE;
}

Bool_t AliRawReaderDate::RewindEvents()
{
// go back to the beginning of the date file

  if (!fFile) return kFALSE;

  fseek(fFile, 0, SEEK_SET);
  fEventNumber = -1;
  return Reset();
}


Int_t AliRawReaderDate::CheckData() const
{
// check the consistency of the data

#ifdef ALI_DATE
  if (!fEvent) return 0;
  // check whether there are sub events
  if (fEvent->eventSize <= fEvent->eventHeadSize) return 0;

  eventHeaderStruct* subEvent = NULL;
  UChar_t* position = 0;
  UChar_t* end = 0;
  Int_t result = 0;

  while (kTRUE) {
    // get the first or the next sub event if at the end of a sub event
    if (!subEvent || (position >= end)) {

      // check for end of event data
      if (position >= ((UChar_t*)fEvent)+fEvent->eventSize) return result;
      if (!TEST_SYSTEM_ATTRIBUTE(fEvent->eventTypeAttribute, 
                                 ATTR_SUPER_EVENT)) {
        subEvent = fEvent;   // no super event
      } else if (subEvent) {
	subEvent = (eventHeaderStruct*) (((UChar_t*)subEvent) + 
					 subEvent->eventSize);
      } else {
	subEvent = (eventHeaderStruct*) (((UChar_t*)fEvent) + 
					 fEvent->eventHeadSize);
      }

      // check the magic word of the sub event
      if (subEvent->eventMagic != EVENT_MAGIC_NUMBER) {
	result |= kErrMagic;
	return result;
      }

      position = ((UChar_t*)subEvent) + subEvent->eventHeadSize + 
	sizeof(equipmentHeaderStruct);
      end = ((UChar_t*)subEvent) + subEvent->eventSize;
    }

    // continue with the next sub event if no data left in the payload
    if (position >= end) continue;

    if (fRequireHeader) {
    // check that there are enough bytes left for the data header
      if (position + sizeof(AliRawDataHeader) > end) {
	result |= kErrNoDataHeader;
	position = end;
	continue;
      }

      // check consistency of data size in the data header and in the sub event
      AliRawDataHeader* header = (AliRawDataHeader*) position;
      if ((position + header->fSize) != end) {
	Warning("ReadHeader",
		"raw data size found in the header is wrong (%d != %d)! Using the equipment size instead !",
		header->fSize, end - position);
	header->fSize = end - position;
	result |= kErrSize;
      }
    }
    position = end;
  };

#endif
  return 0;
}
