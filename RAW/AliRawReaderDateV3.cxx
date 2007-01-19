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
/// This is a class for reading raw data from a date file or event (version 3).
///
/// The AliRawReaderDateV3 is constructed either with a pointer to a
/// date event or with a file name and an event number.
///
///////////////////////////////////////////////////////////////////////////////

#include "AliRawReaderDateV3.h"

struct eventHeaderStruct { 
  Int_t  size;          /* size of event in Bytes  */
  UInt_t magic;         /* magic number used for consistency check */
  UInt_t type;          /* event type */
  UInt_t headLen;       /* size of header in bytes */
  UInt_t runNb;         /* run number */
  UInt_t burstNb;       /* burst number */
  UInt_t nbInRun;       /* event number in run */
  UInt_t nbInBurst;     /* event number in burst */
  UInt_t triggerNb;     /* trigger number for this detector */
  UInt_t fileSeqNb;     /* file sequence number for multifiles run */
  UInt_t detectorId[3]; /* detector identification */
  UInt_t time;          /* time in seconds since 0.00 GMT 1.1.1970 */
  UInt_t usec;          /* microseconds */
  UInt_t errorCode;
  UInt_t deadTime;
  UInt_t deadTimeusec;
  UInt_t typeAttribute[2]; /* event type id mask */
};

#define EVENT_MAGIC_NUMBER         ((UInt_t)0xDA1E5AFE)


ClassImp(AliRawReaderDateV3)


AliRawReaderDateV3::AliRawReaderDateV3(void* event) :
  fFile(NULL),
  fEvent(NULL),
  fSubEvent(NULL),
  fPosition(NULL),
  fEnd(NULL)
{
// create an object to read digits from the given date event

  fEvent = (eventHeaderStruct*) event;
}

AliRawReaderDateV3::AliRawReaderDateV3(const char* fileName, 
				       Int_t eventNumber) :
  fFile(NULL),
  fEvent(NULL),
  fSubEvent(NULL),
  fPosition(NULL),
  fEnd(NULL)
{
// create an object to read digits from the given date event

  fFile = fopen(fileName, "rb");
  if (!fFile) {
    Error("AliRawReaderDateV3", "could not open file %s", fileName);
    return;
  }
  if (eventNumber < 0) return;

  eventHeaderStruct header;
  UInt_t headerSize = sizeof(eventHeaderStruct);
  while (fread(&header, 1, headerSize, fFile) == headerSize) {
    if (eventNumber == 0) {
      UChar_t* buffer = new UChar_t[header.size];
      fseek(fFile, -headerSize, SEEK_CUR);
      if (Int_t(fread(buffer, 1, header.size, fFile)) != header.size) break;
      fEvent = (eventHeaderStruct*) buffer;
      break;
    }
    fseek(fFile, header.size-headerSize, SEEK_CUR);
    eventNumber--;
  }
}

AliRawReaderDateV3::AliRawReaderDateV3(const AliRawReaderDateV3& rawReader) :
  AliRawReader(rawReader),
  fFile(rawReader.fFile),
  fEvent(rawReader.fEvent),
  fSubEvent(rawReader.fSubEvent),
  fPosition(rawReader.fPosition),
  fEnd(rawReader.fEnd)

{
// copy constructor

  Fatal("AliRawReaderDateV3", "copy constructor not implemented");
}

AliRawReaderDateV3& AliRawReaderDateV3::operator = (const AliRawReaderDateV3& 
						/*rawReader*/)
{
// assignment operator

  Fatal("operator =", "assignment operator not implemented");
  return *this;
}

AliRawReaderDateV3::~AliRawReaderDateV3()
{
// destructor

  if (fFile) {
    delete[] fEvent;
    fclose(fFile);
  }
}


UInt_t AliRawReaderDateV3::GetType() const
{
// get the type from the event header

  if (!fEvent) return 0;
  return fEvent->type;
}

UInt_t AliRawReaderDateV3::GetRunNumber() const
{
// get the run number from the event header

  if (!fEvent) return 0;
  return fEvent->runNb;
}

const UInt_t* AliRawReaderDateV3::GetEventId() const
{
// get the event id from the event header

  if (!fEvent) return NULL;
  return &(fEvent->nbInRun);
}

const UInt_t* AliRawReaderDateV3::GetTriggerPattern() const
{
// get the trigger pattern from the event header

  return NULL;
}

const UInt_t* AliRawReaderDateV3::GetDetectorPattern() const
{
// get the detector pattern from the event header

  if (!fEvent) return NULL;
  return fEvent->detectorId;
}

const UInt_t* AliRawReaderDateV3::GetAttributes() const
{
// get the type attributes from the event header

  if (!fEvent) return NULL;
  return fEvent->typeAttribute;
}

const UInt_t* AliRawReaderDateV3::GetSubEventAttributes() const
{
// get the type attributes from the sub event header

  if (!fSubEvent) return NULL;
  return fSubEvent->typeAttribute;
}

UInt_t AliRawReaderDateV3::GetLDCId() const
{
// get the LDC Id from the event header

  return UInt_t(-1);
}

UInt_t AliRawReaderDateV3::GetGDCId() const
{
// get the GDC Id from the event header

  return UInt_t(-1);
}


Int_t AliRawReaderDateV3::GetEquipmentSize() const
{
// get the size of the equipment

  if (!fSubEvent) return 0;
  return fSubEvent->size;
}

Int_t AliRawReaderDateV3::GetEquipmentType() const
{
// get the type from the equipment header

  return 0;
}

Int_t AliRawReaderDateV3::GetEquipmentId() const
{
// get the ID from the equipment header

  return 0;
}

const UInt_t* AliRawReaderDateV3::GetEquipmentAttributes() const
{
// get the attributes from the equipment header

  return 0;
}

Int_t AliRawReaderDateV3::GetEquipmentElementSize() const
{
// get the basic element size from the equipment header

  return 0;
}

Int_t AliRawReaderDateV3::GetEquipmentHeaderSize() const
{
// get the size of the equipment header

  return 0;
}

Bool_t AliRawReaderDateV3::ReadHeader()
{
// read a data header at the current position
// returns kFALSE if the data header could not be read

  fErrorCode = 0;

  fHeader = NULL;
  if (!fEvent) return kFALSE;
  // check whether there are sub events
  if (fEvent->size <= Int_t(fEvent->headLen)) return kFALSE;

  do {
    // skip payload (if event was not selected)
    if (fCount > 0) fPosition += fCount;

    // check for end of event data
    if (fPosition >= ((UChar_t*)fEvent)+fEvent->size) return kFALSE;
    if ((fEvent->detectorId[2] & 0x8000) != 0x8000) {
      fSubEvent = fEvent;   // no super event
    } else if (fSubEvent) {
      fSubEvent = (eventHeaderStruct*) (((UChar_t*)fSubEvent) + 
					fSubEvent->size);
    } else {
      fSubEvent = (eventHeaderStruct*) (((UChar_t*)fEvent) + 
					fEvent->headLen);
    }

    // check the magic word of the sub event
    if (fSubEvent->magic != EVENT_MAGIC_NUMBER) {
      Error("ReadHeader", "wrong magic number in sub event!\n"
	    " run: %d  event: %d\n", 
	    fSubEvent->runNb, fSubEvent->nbInRun);
      fErrorCode = kErrMagic;
      return kFALSE;
    }

    // continue if no data in the subevent
    if (fSubEvent->size == Int_t(fSubEvent->headLen)) {
      fPosition = fEnd = ((UChar_t*)fSubEvent) + fSubEvent->size;
      fCount = 0;
      continue;
    }

    fCount = 0;
    fPosition = ((UChar_t*)fSubEvent) + sizeof(eventHeaderStruct);
    fEnd = ((UChar_t*)fSubEvent) + fSubEvent->size;

    // continue with the next sub event if no data left in the payload
    if (fPosition >= fEnd) continue;

    if (fRequireHeader) {
      // check that there are enough bytes left for the data header
      if (fPosition + sizeof(AliRawDataHeader) > fEnd) {
	Error("ReadHeader", "could not read data header data!");
	Warning("ReadHeader", "skipping %d bytes\n"
		" run: %d  event: %d\n", 
		fEnd - fPosition, fSubEvent->runNb, fSubEvent->nbInRun);
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
		" run: %d  event: %d\n", 
		fEnd - fPosition, fSubEvent->runNb, fSubEvent->nbInRun);
	fCount = 0;
	fPosition = fEnd;
	fErrorCode = kErrSize;
	continue;
      }

    } else {
      fCount = fEnd - fPosition;
    }

  } while (!IsSelected());

  return kTRUE;
}

Bool_t AliRawReaderDateV3::ReadNextData(UChar_t*& data)
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

Bool_t AliRawReaderDateV3::ReadNext(UChar_t* data, Int_t size)
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


Bool_t AliRawReaderDateV3::Reset()
{
// reset the current position to the beginning of the event

  fSubEvent = NULL;
  fCount = 0;
  fPosition = fEnd = NULL;
  return kTRUE;
}


Bool_t AliRawReaderDateV3::NextEvent()
{
// go to the next event in the date file

  if (!fFile) return kFALSE;

  Reset();
  eventHeaderStruct header;
  UInt_t headerSize = sizeof(eventHeaderStruct);
  if (fEvent) delete[] fEvent;
  fEvent = &header;

  while (fread(&header, 1, headerSize, fFile) == headerSize) {
    if (!IsEventSelected()) {
      fseek(fFile, header.size-headerSize, SEEK_CUR);
      continue;
    }
    UChar_t* buffer = new UChar_t[header.size];
    fseek(fFile, -headerSize, SEEK_CUR);
    if (Int_t(fread(buffer, 1, header.size, fFile)) != header.size) {
      Error("NextEvent", "could not read event from file");
      delete[] buffer;
      break;
    }
    fEvent = (eventHeaderStruct*) buffer;
    fEventNumber++;
    return kTRUE;
  };

  fEvent = NULL;
  return kFALSE;
}

Bool_t AliRawReaderDateV3::RewindEvents()
{
// go back to the beginning of the date file

  if (!fFile) return kFALSE;

  fseek(fFile, 0, SEEK_SET);
  fEventNumber = -1;
  return Reset();
}


Int_t AliRawReaderDateV3::CheckData() const
{
// check the consistency of the data

  if (!fEvent) return 0;
  // check whether there are sub events
  if (fEvent->size <= Int_t(fEvent->headLen)) return 0;

  eventHeaderStruct* subEvent = NULL;
  UChar_t* position = 0;
  UChar_t* end = 0;
  Int_t result = 0;

  while (kTRUE) {
    // check for end of event data
    if (position >= ((UChar_t*)fEvent)+fEvent->size) return result;
    if ((fEvent->detectorId[2] & 0x8000) != 0x8000) {
      subEvent = fEvent;   // no super event
    } else if (subEvent) {
      subEvent = (eventHeaderStruct*) (((UChar_t*)subEvent) + 
				       subEvent->size);
    } else {
      subEvent = (eventHeaderStruct*) (((UChar_t*)fEvent) + 
				       fEvent->headLen);
    }

    // check the magic word of the sub event
    if (subEvent->magic != EVENT_MAGIC_NUMBER) {
      result |= kErrMagic;
      return result;
    }

    position = ((UChar_t*)subEvent) + subEvent->headLen;
    end = ((UChar_t*)subEvent) + subEvent->size;

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

  return 0;
}
