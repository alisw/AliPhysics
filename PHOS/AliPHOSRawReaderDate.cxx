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
/// The AliPHOSRawReaderDate is constructed either with a pointer to a
/// date event or with a file name and an event number.
///
///////////////////////////////////////////////////////////////////////////////
#include <stdio.h>

#include "AliPHOSRawReaderDate.h"
#include "Riostream.h"

ClassImp(AliPHOSRawReaderDate)

//---------------------------------------------------------
AliPHOSRawReaderDate::AliPHOSRawReaderDate(void* event) :
  fRequireHeader(kTRUE),
  fFile(NULL),
  fEvent(NULL),
  fSubEvent(NULL),
  fEquipment(NULL),
  fPosition(NULL),
  fEnd(NULL)
{
// create an object to read digits from the given date event
  fEvent = (eventHeaderStruct*) event;
}
//---------------------------------------------------------
AliPHOSRawReaderDate::AliPHOSRawReaderDate(const char* fileName, Int_t eventNumber):
  fRequireHeader(kTRUE),
  fFile(NULL),
  fEvent(NULL),
  fSubEvent(NULL),
  fEquipment(NULL),
  fPosition(NULL),
  fEnd(NULL)
{
// create an object to read digits from the given date event
//   char command[256];
//   if(strstr(fileName,".gz") )
//     sprintf(command,"zcat %s",fileName);
//   else
//     sprintf(command,"cat %s",fileName);    
//   printf("Comand %s \n",command) ;
//   fFile  = popen("zcat Run_3186.dat.gz", "rb");

  fFile  = fopen(fileName, "rb");

  if (!fFile) {
    Error("AliPHOSRawReaderDate", "could not open file %s", fileName);
    return;
  }
  if (eventNumber < 0) return;

  eventHeaderStruct header;
  UInt_t headerSize = sizeof(eventHeaderStruct);
  while (fread(&header, 1, headerSize, fFile) == headerSize) {
    if (eventNumber == 0) {
      UChar_t* buffer = new UChar_t[header.eventSize];
      fseek(fFile, -headerSize, SEEK_CUR);
      if (fread(buffer, 1, header.eventSize, fFile) != (UInt_t)header.eventSize) break;
      fEvent = (eventHeaderStruct*) buffer;
      break;
    }
    fseek(fFile, header.eventSize-headerSize, SEEK_CUR);
    eventNumber--;
  }
}
//---------------------------------------------------------
AliPHOSRawReaderDate::AliPHOSRawReaderDate(const AliPHOSRawReaderDate& rawReader) :
  AliRawReader(rawReader),
  fRequireHeader(rawReader.fRequireHeader),
  fFile(rawReader.fFile),
  fEvent(rawReader.fEvent),
  fSubEvent(rawReader.fSubEvent),
  fEquipment(rawReader.fEquipment),
  fPosition(rawReader.fPosition),
  fEnd(rawReader.fEnd)

{
// copy constructor

  Fatal("AliPHOSRawReaderDate", "copy constructor not implemented");
}
//---------------------------------------------------------
AliPHOSRawReaderDate& AliPHOSRawReaderDate::operator = (const AliPHOSRawReaderDate& 
						/*rawReader*/)
{
// assignment operator

  Fatal("operator =", "assignment operator not implemented");
  return *this;
}
//---------------------------------------------------------
AliPHOSRawReaderDate::~AliPHOSRawReaderDate()
{
// destructor

  if (fFile) {
    delete[] fEvent;
    //   pclose(fFile);
    fclose(fFile);
  }
}
//---------------------------------------------------------
Bool_t AliPHOSRawReaderDate::ReadHeader()
{
// read a data header at the current position
// returns kFALSE if the data header could not be read

  fErrorCode = 0;

  fHeader = NULL;
  if (!fEvent) return kFALSE;

  // Analize, if event header good to use. Swap event, if need.
  if (fEvent->eventMagic != EVENT_MAGIC_NUMBER){
    if (fEvent->eventMagic == EVENT_MAGIC_NUMBER_SWAPPED) {
      SwappEvent(fEvent);
    } 
    else {
      Error("ReadHader","Wrong event magic number, MAGIC == %08x (expected %08x).\n",
	    fEvent->eventMagic, EVENT_MAGIC_NUMBER);
      fErrorCode = kErrMagic;
      return kFALSE;
    }
  }

  // check whether there are sub events
  if ((UInt_t)fEvent->eventSize <= fEvent->eventHeadSize) return kFALSE;
  
  do {
    // skip payload (if event was not selected)
    if (fCount > 0) fPosition += fCount;
    
    // get the first or the next equipment if at the end of an equipment
    if (!fEquipment || (fPosition >= fEnd)) {
      fEquipment = NULL;
     
      // printf("No equipment defined \n") ;
 
      // get the first or the next sub event if at the end of a sub event
      if (!fSubEvent || 
	  (fPosition >= ((UChar_t*)fSubEvent) + fSubEvent->eventSize)) {

	// printf("GetNext subevent \n") ;

	// printf("fPosition fPosition=%p fEvent=%p Size=%d Ev+size%p \n",(UChar_t*)fPosition,(UChar_t*)fEvent,fEvent->eventSize,(((UChar_t*)fEvent)+fEvent->eventSize)) ;
	
	// check for end of event data
	if (fPosition >= ((UChar_t*)fEvent)+fEvent->eventSize) return kFALSE;

	// printf("position passed \n") ;

	if(!(fEvent->detectorId[0] & SUPER_EVENT_MASK)){
	  fSubEvent = fEvent;   // no super event
	  // printf("superevent \n") ;
	} else if (fSubEvent) {
	  fSubEvent = (eventHeaderStruct*) (((UChar_t*)fSubEvent) + 
					    fSubEvent->eventSize);
	  // printf("second subevent \n ") ;
	} else {
	  fSubEvent = (eventHeaderStruct*) (((UChar_t*)fEvent) + 
					    fEvent->eventHeadSize);
	  // printf("first subevent \n") ;
	}
	
	// check the magic word of the sub event
	if (fSubEvent->eventMagic != EVENT_MAGIC_NUMBER) {
	  if (fSubEvent->eventMagic == EVENT_MAGIC_NUMBER_SWAPPED) {
	    SwappEvent(fSubEvent);
	  } 
	  else {
	    Error("ReadHeader", "wrong magic number in sub event!\n"
		  " run: %d \n",fSubEvent->eventRunNb) ; 
	    fErrorCode = kErrMagic;
	    return kFALSE;
	  }
	}
	
	// continue if no data in the subevent
	if ((UInt_t)fSubEvent->eventSize == fSubEvent->eventHeadSize) {
	  fPosition = fEnd = ((UChar_t*)fSubEvent) + fSubEvent->eventSize;
	  fCount = 0;
	  // printf("no data in subevent \n") ;
	  continue;
	}
	
	// printf("Set equipment after header of subevent \n") ;
	fEquipment = (equipmentHeaderStruct*)
	  (((UChar_t*)fSubEvent) + fSubEvent->eventHeadSize);
	
      } else {
	fEquipment = (equipmentHeaderStruct*) fEnd;
	// printf("equipment set at fEnd \n" ) ;
      }
      
      fCount = 0;
      fPosition = ((UChar_t*)fEquipment) + sizeof(equipmentHeaderStruct);
      fEnd = fPosition + fEquipment->rawDataLen;
    }
    
    // continue with the next sub event if no data left in the payload
    if (fPosition >= fEnd) continue;
    
    fCount = fEnd - fPosition;
    //    if(IsSelected())
      // printf("Selected \n") ;
    // else
      // printf("Not selected \n") ;
    
  } while (!fEquipment || !IsSelected());
  
  return kTRUE;
}
//---------------------------------------------------------
Bool_t AliPHOSRawReaderDate::ReadNextData(UChar_t*& data)
{
// reads the next payload at the current position
// returns kFALSE if the data could not be read

  fErrorCode = 0;
  while (fCount == 0) {
    if (!ReadHeader()) return kFALSE;
  }
  data = fPosition;
  fPosition += fCount;  
  // printf("fCount %d \n",fCount) ;
  fCount = 0;
  return kTRUE;
}
//---------------------------------------------------------
Bool_t AliPHOSRawReaderDate::ReadNext(UChar_t* data, Int_t size)
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
//---------------------------------------------------------
Bool_t AliPHOSRawReaderDate::Reset()
{
// reset the current position to the beginning of the event

  fSubEvent = NULL;
  fCount = 0;
  fPosition = fEnd = NULL;
  return kTRUE;
}
//---------------------------------------------------------
Bool_t AliPHOSRawReaderDate::NextEvent()
{
// go to the next event in the date file

  if (!fFile) return kFALSE;

  // printf("Old fEvent = %p \n",(UChar_t*)fEvent) ;
  eventHeaderStruct header;
  UInt_t headerSize = sizeof(eventHeaderStruct);
  if (fEvent) delete[] fEvent;
  fEvent = &header;
  // printf("New fEvent = %p \n",(UChar_t*)fEvent) ; 
  
  while (fread(&header, 1, headerSize, fFile) == headerSize) {
    // printf("read \n") ;
    if (!IsEventSelected()) {
      // printf("selected \n") ;
      fseek(fFile, header.eventSize-headerSize, SEEK_CUR);
      continue;
    }
    UChar_t* buffer = new UChar_t[header.eventSize];
    fseek(fFile, -headerSize, SEEK_CUR);
    if (fread(buffer, 1, header.eventSize, fFile) != (UInt_t)header.eventSize) {
      Error("NextEvent", "could not read event from file");
      delete[] buffer;
      break;
    }
    fEvent = (eventHeaderStruct*) buffer;
    // printf("Read fEvent = %p \n",(UChar_t*)fEvent) ;
    fSubEvent = NULL ;
    fEquipment = NULL ;
    fPosition = NULL ;
    fEnd = NULL ;
    return kTRUE;
  };

  fEvent = NULL;
  return kFALSE;
}
//---------------------------------------------------------
Bool_t AliPHOSRawReaderDate::RewindEvents()
{
// go back to the beginning of the date file

  if (!fFile) return kFALSE;

  fseek(fFile, 0, SEEK_SET);
  return Reset();
}
//---------------------------------------------------------
Int_t AliPHOSRawReaderDate::CheckData() const
{
// check the consistency of the data

  if (!fEvent) return 0;
  // check whether there are sub events
  if ((UInt_t)fEvent->eventSize <= fEvent->eventHeadSize) return 0;

  eventHeaderStruct* subEvent = NULL;
  UChar_t* position = 0;
  UChar_t* end = 0;
  Int_t result = 0;

  while (kTRUE) {
    // get the first or the next sub event if at the end of a sub event
    if (!subEvent || (position >= end)) {

      // check for end of event data
      if (position >= ((UChar_t*)fEvent)+fEvent->eventSize) return result;
      //???DP
//       if (!TEST_SYSTEM_ATTRIBUTE(fEvent->eventTypeAttribute, 
//                                  ATTR_SUPER_EVENT)) {
        subEvent = fEvent;   // no super event
//       } else if (subEvent) {
// 	subEvent = (eventHeaderStruct*) (((UChar_t*)subEvent) + 
// 					 subEvent->eventSize);
//       } else {
// 	subEvent = (eventHeaderStruct*) (((UChar_t*)fEvent) + 
// 					 fEvent->eventHeadSize);
//       }

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

    // check that there are enough bytes left for the data header
    if (position + sizeof(AliRawDataHeader) > end) {
      result |= kErrNoDataHeader;
      position = end;
      continue;
    }

    // check consistency of data size in the data header and in the sub event
    AliRawDataHeader* header = (AliRawDataHeader*) position;
    if (header->fSize != 0xFFFFFFFF) {
      if (position + header->fSize > end) {
	result |= kErrSize;
	position = end;
      } else {
	position += header->fSize;
      }
    } else {
      position = end;
    }
  };

  return 0;
}
//---------------------------------------------------------
void AliPHOSRawReaderDate::SwappEvent(eventHeaderStruct * event)
{
  if(event->eventMagic == EVENT_MAGIC_NUMBER_SWAPPED)
    {
      ChangeOrder(event->eventSize);
      ChangeOrder(event->eventType);
      ChangeOrder(event->detectorId[0]);
    }
}
