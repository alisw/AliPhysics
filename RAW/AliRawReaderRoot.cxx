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
// This is a class for reading a raw data from a root file and providing
// information about digits
//
///////////////////////////////////////////////////////////////////////////////

#include <TFile.h>
#include <TTree.h>
#include "AliRawReaderRoot.h"
#include "AliRawEvent.h"
#include "AliRawEventHeader.h"
#include "AliRawEquipmentHeader.h"
#include "AliRawData.h"


ClassImp(AliRawReaderRoot)


AliRawReaderRoot::AliRawReaderRoot(const char* fileName, Int_t eventNumber)
{
// create an object to read digits from the given input file for the
// event with the given number

  fEvent = NULL;
  TDirectory* dir = gDirectory;
  fFile = TFile::Open(fileName);
  dir->cd();
  if (!fFile || !fFile->IsOpen()) {
    Error("AliRawReaderRoot", "could not open file %s", fileName);
    return;
  }
  TTree* tree = (TTree*) fFile->Get("RAW");
  if (!tree) {
    Error("AliRawReaderRoot", "no raw data tree found");
    return;
  }
  TBranch* branch = tree->GetBranch("rawevent");
  if (!branch) {
    Error("AliRawReaderRoot", "no raw data branch found");
    return;
  }

  fEvent = new AliRawEvent;
  branch->SetAddress(&fEvent);
  if (branch->GetEntry(eventNumber) <= 0) {
    Error("AliRawReaderRoot", "no event with number %d found", eventNumber);
    return;
  }
  
  fSubEventIndex = 0;
  fSubEvent = NULL;
  fRawData = NULL;
  fHeader = NULL;

  fCount = 0;
  fPosition = fEnd = NULL;
}

AliRawReaderRoot::AliRawReaderRoot(AliRawEvent* event)
{
// create an object to read digits from the given raw event

  fFile = NULL;
  fEvent = event;
  
  fSubEventIndex = 0;
  fSubEvent = NULL;
  fRawData = NULL;
  fHeader = NULL;

  fCount = 0;
  fPosition = fEnd = NULL;
}

AliRawReaderRoot::AliRawReaderRoot(const AliRawReaderRoot& rawReader) :
  AliRawReader(rawReader)
{
// copy constructor

  fFile = NULL;
  fEvent = rawReader.fEvent;
  
  fSubEventIndex = rawReader.fSubEventIndex;
  fSubEvent = rawReader.fSubEvent;
  fRawData = rawReader.fRawData;
  fHeader = rawReader.fHeader;

  fCount = rawReader.fCount;
  fPosition = rawReader.fPosition;
  fEnd = rawReader.fEnd;
}

AliRawReaderRoot& AliRawReaderRoot::operator = (const AliRawReaderRoot& 
						rawReader)
{
// assignment operator

  this->~AliRawReaderRoot();
  new(this) AliRawReaderRoot(rawReader);
  return *this;
}

AliRawReaderRoot::~AliRawReaderRoot()
{
// delete objects and close root file

  if (fFile) {
    if (fEvent) delete fEvent;
    fFile->Close();
    delete fFile;
  }
}


UInt_t AliRawReaderRoot::GetType() const
{
// get the type from the event header

  if (!fEvent) return 0;
  return fEvent->GetHeader()->GetType();
}

UInt_t AliRawReaderRoot::GetRunNumber() const
{
// get the run number from the event header

  if (!fEvent) return 0;
  return fEvent->GetHeader()->GetRunNumber();
}

const UInt_t* AliRawReaderRoot::GetEventId() const
{
// get the event id from the event header

  if (!fEvent) return NULL;
  return fEvent->GetHeader()->GetId();
}

const UInt_t* AliRawReaderRoot::GetTriggerPattern() const
{
// get the trigger pattern from the event header

  if (!fEvent) return NULL;
  return fEvent->GetHeader()->GetTriggerPattern();
}

const UInt_t* AliRawReaderRoot::GetDetectorPattern() const
{
// get the detector pattern from the event header

  if (!fEvent) return NULL;
  return fEvent->GetHeader()->GetDetectorPattern();
}

const UInt_t* AliRawReaderRoot::GetAttributes() const
{
// get the type attributes from the event header

  if (!fEvent) return NULL;
  return fEvent->GetHeader()->GetTypeAttribute();
}

UInt_t AliRawReaderRoot::GetLDCId() const
{
// get the LDC Id from the event header

  if (!fEvent || !fEvent->GetSubEvent(fSubEventIndex)) return 0;
  return fEvent->GetSubEvent(fSubEventIndex)->GetHeader()->GetLDCId();
}

UInt_t AliRawReaderRoot::GetGDCId() const
{
// get the GDC Id from the event header

  if (!fEvent) return 0;
  return fEvent->GetHeader()->GetGDCId();
}


Int_t AliRawReaderRoot::GetEquipmentSize() const
{
// get the size of the equipment

  if (!fEvent || !fEvent->GetEquipmentHeader()) return 0;
  return fEvent->GetEquipmentHeader()->GetEquipmentSize();
}

Int_t AliRawReaderRoot::GetEquipmentType() const
{
// get the type from the equipment header

  if (!fEvent || !fEvent->GetEquipmentHeader()) return -1;
  return fEvent->GetEquipmentHeader()->GetEquipmentType();
}

Int_t AliRawReaderRoot::GetEquipmentId() const
{
// get the ID from the equipment header

  if (!fEvent || !fEvent->GetEquipmentHeader()) return -1;
  return fEvent->GetEquipmentHeader()->GetId();
}

const UInt_t* AliRawReaderRoot::GetEquipmentAttributes() const
{
// get the attributes from the equipment header

  if (!fEvent || !fEvent->GetEquipmentHeader()) return NULL;
  return fEvent->GetEquipmentHeader()->GetTypeAttribute();
}

Int_t AliRawReaderRoot::GetEquipmentElementSize() const
{
// get the basic element size from the equipment header

  if (!fEvent || !fEvent->GetEquipmentHeader()) return 0;
  return fEvent->GetEquipmentHeader()->GetBasicSizeType();
}


Bool_t AliRawReaderRoot::ReadHeader()
{
// read a data header at the current position
// returns kFALSE if the data header could not be read

  fErrorCode = 0;
  if (!fEvent) return kFALSE;

  do {
    // skip payload (if event was not selected)
    if (fCount > 0) fPosition += fCount;

    // get the first or the next sub event if at the end of a sub event
    if (!fSubEvent || (fPosition >= fEnd)) {

      // check for end of event data
      if (fSubEventIndex >= fEvent->GetNSubEvents()) return kFALSE;
      fSubEvent = fEvent->GetSubEvent(fSubEventIndex++);

      // check the magic word of the sub event
      if (!fSubEvent->GetHeader()->IsValid()) {
	Error("ReadHeader", "wrong magic number in sub event!");
	fSubEvent->GetHeader()->Dump();
	fErrorCode = kErrMagic;
	return kFALSE;
      }

      fRawData = fSubEvent->GetRawData();
      fCount = 0;
      fPosition = (UChar_t*) fRawData->GetBuffer();
      fEnd = ((UChar_t*) fRawData->GetBuffer()) + fRawData->GetSize();
    }

    // continue with the next sub event if no data left in the payload
    if (fPosition >= fEnd) continue;

    // check that there are enough bytes left for the data header
    if (fPosition + sizeof(AliRawDataHeader) > fEnd) {
      Error("ReadHeader", "could not read data header!");
      Warning("ReadHeader", "skipping %d bytes", fEnd - fPosition);
      fSubEvent->GetHeader()->Dump();
      fCount = 0;
      fPosition = fEnd;
      fErrorCode = kErrNoDataHeader;
      continue;
    }

    // "read" the data header
    fHeader = (AliRawDataHeader*) fPosition;
    fPosition += sizeof(AliRawDataHeader);
    if (fHeader->fSize != 0xFFFFFFFF) {
      fCount = fHeader->fSize - sizeof(AliRawDataHeader);
    } else {
      fCount = fEnd - fPosition;
    }

    // check consistency of data size in the header and in the sub event
    if (fPosition + fCount > fEnd) {  
      Error("ReadHeader", "size in data header exceeds event size!");
      Warning("ReadHeader", "skipping %d bytes", fEnd - fPosition);
      fSubEvent->GetHeader()->Dump();
      fCount = 0;
      fPosition = fEnd;
      fErrorCode = kErrSize;
      continue;
    }

  } while (!IsSelected());

  return kTRUE;
}

Bool_t AliRawReaderRoot::ReadNextData(UChar_t*& data)
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

Bool_t AliRawReaderRoot::ReadNext(UChar_t* data, Int_t size)
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


Bool_t AliRawReaderRoot::Reset()
{
// reset the current position to the beginning of the event

  fSubEventIndex = 0;
  fSubEvent = NULL;
  fRawData = NULL;
  fHeader = NULL;

  fCount = 0;
  fPosition = fEnd = NULL;
  return kTRUE;
}


Int_t AliRawReaderRoot::CheckData() const
{
// check the consistency of the data

  if (!fEvent) return 0;

  AliRawEvent* subEvent = NULL;
  Int_t subEventIndex = 0;
  UChar_t* position = 0;
  UChar_t* end = 0;
  Int_t result = 0;

  while (kTRUE) {
    // get the first or the next sub event if at the end of a sub event
    if (!subEvent || (position >= end)) {

      // check for end of event data
      if (subEventIndex >= fEvent->GetNSubEvents()) return result;
      subEvent = fEvent->GetSubEvent(subEventIndex++);

      // check the magic word of the sub event
      if (!fSubEvent->GetHeader()->IsValid()) {
	result |= kErrMagic;
	return result;
      }

      AliRawData* rawData = subEvent->GetRawData();
      position = (UChar_t*) rawData->GetBuffer();
      end = ((UChar_t*) rawData->GetBuffer()) + rawData->GetSize();
    }

    // continue with the next sub event if no data left in the payload
    if (position >= end) continue;

    // check that there are enough bytes left for the data header
    if (position + sizeof(AliRawDataHeader) > end) {
      result |= kErrNoDataHeader;
      position = end;
      continue;
    }

    // check consistency of data size in the header and in the sub event
    AliRawDataHeader* header = (AliRawDataHeader*) position;
    if (fHeader->fSize != 0xFFFFFFFF) {
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

  return result;
}
