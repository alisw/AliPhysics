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

#include "AliRawReaderRoot.h"
#include "AliRawEvent.h"


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
  fMiniHeader = NULL;

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
  fMiniHeader = NULL;

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
  fMiniHeader = rawReader.fMiniHeader;

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

UInt_t AliRawReaderRoot::GetGDCId() const
{
// get the GDC Id from the event header

  if (!fEvent) return 0;
  return fEvent->GetHeader()->GetGDCId();
}


Bool_t AliRawReaderRoot::ReadMiniHeader()
{
// read a mini header at the current position
// returns kFALSE if the mini header could not be read

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
	Error("ReadMiniHeader", "wrong magic number in sub event!");
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

    // check that there are enough bytes left for the mini header
    if (fPosition + sizeof(AliMiniHeader) > fEnd) {
      Error("ReadMiniHeader", "could not read mini header data!");
      Warning("ReadMiniHeader", "skipping %d bytes", fEnd - fPosition);
      fSubEvent->GetHeader()->Dump();
      fCount = 0;
      fPosition = fEnd;
      fErrorCode = kErrNoMiniHeader;
      continue;
    }

    // "read" and check the mini header
    fMiniHeader = (AliMiniHeader*) fPosition;
    fPosition += sizeof(AliMiniHeader);
    if (!CheckMiniHeader()) {
      Error("ReadMiniHeader", "wrong magic word in mini header!");
      Warning("ReadMiniHeader", "skipping %d bytes", fEnd - fPosition);
      fSubEvent->GetHeader()->Dump();
      fCount = 0;
      fPosition = fEnd;
      fErrorCode = kErrMiniMagic;
      continue;
    }
    fCount = fMiniHeader->fSize;

    // check consistency of data size in the mini header and in the sub event
    if (fPosition + fCount > fEnd) {  
      Error("ReadMiniHeader", "size in mini header exceeds event size!");
      Warning("ReadMiniHeader", "skipping %d bytes", fEnd - fPosition);
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
    if (!ReadMiniHeader()) return kFALSE;
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
  fMiniHeader = NULL;

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

    // check that there are enough bytes left for the mini header
    if (position + sizeof(AliMiniHeader) > end) {
      result |= kErrNoMiniHeader;
      position = end;
      continue;
    }

    // "read" and check the mini header
    AliMiniHeader* miniHeader = (AliMiniHeader*) position;
    position += sizeof(AliMiniHeader);
    if (!CheckMiniHeader(miniHeader)) {
      result |= kErrMiniMagic;
      position = end;
      continue;
    }

    // check consistency of data size in the mini header and in the sub event
    if (position + miniHeader->fSize > end) result |= kErrSize;
    position += miniHeader->fSize;
  };

  return result;
}
