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
/// This is a class for reading raw data from a root file.
///
/// The root file is expected to contain a tree of name "RAW" with
/// a branch of name "rawevent" which contains objects of type
/// AliRawEvent.
/// 
/// The file name and the event number are arguments of the constructor
/// of AliRawReaderRoot.
///
///////////////////////////////////////////////////////////////////////////////

#include <TFile.h>
#include <TTree.h>
#include "AliRawReaderRoot.h"
#include "AliRawEvent.h"
#include "AliRawEventHeaderBase.h"
#include "AliRawEquipment.h"
#include "AliRawEquipmentHeader.h"
#include "AliRawData.h"


ClassImp(AliRawReaderRoot)


AliRawReaderRoot::AliRawReaderRoot(const char* fileName, Int_t eventNumber) :
  fFile(NULL),
  fBranch(NULL),
  fEventIndex(eventNumber),
  fEvent(NULL),
  fSubEventIndex(0),
  fSubEvent(NULL),
  fEquipmentIndex(0),
  fEquipment(NULL),
  fRawData(NULL),
  fPosition(NULL),
  fEnd(NULL)
{
// create an object to read digits from the given input file for the
// event with the given number

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
  fBranch = tree->GetBranch("rawevent");
  if (!fBranch) {
    Error("AliRawReaderRoot", "no raw data branch found");
    return;
  }

  fEvent = new AliRawEvent;
  fBranch->SetAddress(&fEvent);
  if (fEventIndex >= 0) {
    if (fBranch->GetEntry(fEventIndex) <= 0) {
      Error("AliRawReaderRoot", "no event with number %d found", fEventIndex);
      return;
    }
  }
}

AliRawReaderRoot::AliRawReaderRoot(AliRawEvent* event) :
  fFile(NULL),
  fBranch(NULL),
  fEventIndex(-1),
  fEvent(event),
  fSubEventIndex(0),
  fSubEvent(NULL),
  fEquipmentIndex(0),
  fEquipment(NULL),
  fRawData(NULL),
  fPosition(NULL),
  fEnd(NULL)
{
// create an object to read digits from the given raw event

}

AliRawReaderRoot::AliRawReaderRoot(const AliRawReaderRoot& rawReader) :
  AliRawReader(rawReader),
  fFile(NULL),
  fBranch(NULL),
  fEventIndex(rawReader.fEventIndex),
  fEvent(NULL),
  fSubEventIndex(rawReader.fSubEventIndex),
  fSubEvent(NULL),
  fEquipmentIndex(rawReader.fEquipmentIndex),
  fEquipment(NULL),
  fRawData(NULL),
  fPosition(NULL),
  fEnd(NULL)
{
// copy constructor

  if (rawReader.fFile) {
    TDirectory* dir = gDirectory;
    fFile = TFile::Open(rawReader.fFile->GetName());
    dir->cd();
    if (!fFile || !fFile->IsOpen()) {
      Error("AliRawReaderRoot", "could not open file %s", 
	    rawReader.fFile->GetName());
      return;
    }
    TTree* tree = (TTree*) fFile->Get("RAW");
    if (!tree) {
      Error("AliRawReaderRoot", "no raw data tree found");
      return;
    }
    fBranch = tree->GetBranch("rawevent");
    if (!fBranch) {
      Error("AliRawReaderRoot", "no raw data branch found");
      return;
    }

    fEvent = new AliRawEvent;
    fBranch->SetAddress(&fEvent);
    if (fEventIndex >= 0) {
      if (fBranch->GetEntry(fEventIndex) <= 0) {
	Error("AliRawReaderRoot", "no event with number %d found", 
	      fEventIndex);
	return;
      }
    }
  } else {
    fEvent = rawReader.fEvent;
  }

  if (fSubEventIndex > 0) {
    fSubEvent = fEvent->GetSubEvent(fSubEventIndex-1);
    fEquipment = fSubEvent->GetEquipment(fEquipmentIndex);
    fRawData = fEquipment->GetRawData();
      fCount = 0;
    fHeader = (AliRawDataHeader*) ((UChar_t*) fRawData->GetBuffer() + 
      ((UChar_t*) rawReader.fHeader - 
       (UChar_t*) rawReader.fRawData->GetBuffer()));
    fPosition = (UChar_t*) fRawData->GetBuffer() + 
      (rawReader.fPosition - (UChar_t*) rawReader.fRawData->GetBuffer());
    fEnd = ((UChar_t*) fRawData->GetBuffer()) + fRawData->GetSize();
  }
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

const AliRawEventHeaderBase* AliRawReaderRoot::GetEventHeader() const
{
  // Get the even header
  // Return NULL in case of failure
  if (!fEvent) return NULL;
  return fEvent->GetHeader();
}

UInt_t AliRawReaderRoot::GetType() const
{
// get the type from the event header

  if (!fEvent) return 0;
  return fEvent->GetHeader()->Get("Type");
}

UInt_t AliRawReaderRoot::GetRunNumber() const
{
// get the run number from the event header

  if (!fEvent) return 0;
  return fEvent->GetHeader()->Get("RunNb");
}

const UInt_t* AliRawReaderRoot::GetEventId() const
{
// get the event id from the event header

  if (!fEvent) return NULL;
  return fEvent->GetHeader()->GetP("Id");
}

const UInt_t* AliRawReaderRoot::GetTriggerPattern() const
{
// get the trigger pattern from the event header

  if (!fEvent) return NULL;
  return fEvent->GetHeader()->GetP("TriggerPattern");
}

const UInt_t* AliRawReaderRoot::GetDetectorPattern() const
{
// get the detector pattern from the event header

  if (!fEvent) return NULL;
  return fEvent->GetHeader()->GetP("DetectorPattern");
}

const UInt_t* AliRawReaderRoot::GetAttributes() const
{
// get the type attributes from the event header

  if (!fEvent) return NULL;
  return fEvent->GetHeader()->GetP("TypeAttribute");
}

const UInt_t* AliRawReaderRoot::GetSubEventAttributes() const
{
// get the type attributes from the sub event header

  if (!fSubEvent) return NULL;
  return fSubEvent->GetHeader()->GetP("TypeAttribute");
}

UInt_t AliRawReaderRoot::GetLDCId() const
{
// get the LDC Id from the event header

  if (!fEvent || !fSubEvent) return 0;
  return fSubEvent->GetHeader()->Get("LdcId");
}

UInt_t AliRawReaderRoot::GetGDCId() const
{
// get the GDC Id from the event header

  if (!fEvent) return 0;
  return fEvent->GetHeader()->Get("GdcId");
}


Int_t AliRawReaderRoot::GetEquipmentSize() const
{
// get the size of the equipment

  if (!fEvent || !fEquipment || !fEquipment->GetEquipmentHeader()) return 0;
  return fEquipment->GetEquipmentHeader()->GetEquipmentSize();
}

Int_t AliRawReaderRoot::GetEquipmentType() const
{
// get the type from the equipment header

  if (!fEvent || !fEquipment || !fEquipment->GetEquipmentHeader()) return -1;
  return fEquipment->GetEquipmentHeader()->GetEquipmentType();
}

Int_t AliRawReaderRoot::GetEquipmentId() const
{
// get the ID from the equipment header

  if (!fEvent || !fEquipment || !fEquipment->GetEquipmentHeader()) return -1;
  return fEquipment->GetEquipmentHeader()->GetId();
}

const UInt_t* AliRawReaderRoot::GetEquipmentAttributes() const
{
// get the attributes from the equipment header

  if (!fEvent || !fEquipment || !fEquipment->GetEquipmentHeader()) return NULL;
  return fEquipment->GetEquipmentHeader()->GetTypeAttribute();
}

Int_t AliRawReaderRoot::GetEquipmentElementSize() const
{
// get the basic element size from the equipment header

  if (!fEvent || !fEquipment || !fEquipment->GetEquipmentHeader()) return 0;
  return fEquipment->GetEquipmentHeader()->GetBasicSizeType();
}

Int_t AliRawReaderRoot::GetEquipmentHeaderSize() const
{
// get the size of the equipment header (28 bytes by default)

  if (!fEvent || !fEquipment || !fEquipment->GetEquipmentHeader()) return 0;
  return fEquipment->GetEquipmentHeader()->HeaderSize();
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

    // get the first or the next equipment if at the end of an equipment
    if (!fEquipment || (fPosition >= fEnd)) {

      // get the first or the next sub event if at the end of a sub event
      if (!fSubEvent || (fEquipmentIndex >= fSubEvent->GetNEquipments())) {

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

	fEquipmentIndex = 0;
	fEquipment = NULL;
	fRawData = NULL;
      }

      // get the next equipment and raw data
      fCount = 0;
      fEquipment = fSubEvent->GetEquipment(fEquipmentIndex++);
      if (!fEquipment) continue;
      fRawData = fEquipment->GetRawData();
      if (!fRawData) {
	fPosition = fEnd;
	continue;
      }
      fPosition = (UChar_t*) fRawData->GetBuffer();
      fEnd = ((UChar_t*) fRawData->GetBuffer()) + fRawData->GetSize();
    }

    // continue with the next equipment if no data left in the payload
    if (fPosition >= fEnd) continue;

    if (fRequireHeader) {
      // check that there are enough bytes left for the data header
      if (fPosition + sizeof(AliRawDataHeader) > fEnd) {
	Error("ReadHeader", "could not read data header!");
	Warning("ReadHeader", "skipping %d bytes", fEnd - fPosition);
	fEquipment->GetEquipmentHeader()->Dump();
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
	Warning("ReadHeader", "skipping %d bytes", fEnd - fPosition);
	fEquipment->GetEquipmentHeader()->Dump();
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
  fEquipmentIndex = 0;
  fEquipment = NULL;
  fRawData = NULL;
  fHeader = NULL;

  fCount = 0;
  fPosition = fEnd = NULL;
  return kTRUE;
}


Bool_t AliRawReaderRoot::NextEvent()
{
// go to the next event in the root file

  if (!fFile) return kFALSE;

  do {
    delete fEvent;
    fEvent = new AliRawEvent;
    fBranch->SetAddress(&fEvent);
    if (fBranch->GetEntry(fEventIndex+1) <= 0)
      return kFALSE;
    fEventIndex++;
  } while (!IsEventSelected());
  return Reset();
}

Bool_t AliRawReaderRoot::RewindEvents()
{
// go back to the beginning of the root file

  if (!fFile) return kFALSE;

  fEventIndex = -1;
  delete fEvent;
  fEvent = new AliRawEvent;
  fBranch->SetAddress(&fEvent);
  return Reset();
}


Int_t AliRawReaderRoot::CheckData() const
{
// check the consistency of the data

  if (!fEvent) return 0;

  AliRawEvent* subEvent = NULL;
  Int_t subEventIndex = 0;
  AliRawEquipment* equipment = NULL;
  Int_t equipmentIndex = 0;
  UChar_t* position = 0;
  UChar_t* end = 0;
  Int_t result = 0;

  while (kTRUE) {
    // get the first or the next sub event if at the end of an equipment
    if (!subEvent || (equipmentIndex >= subEvent->GetNEquipments())) {

      // check for end of event data
      if (subEventIndex >= fEvent->GetNSubEvents()) return result;
      subEvent = fEvent->GetSubEvent(subEventIndex++);

      // check the magic word of the sub event
      if (!fSubEvent->GetHeader()->IsValid()) {
	result |= kErrMagic;
	return result;
      }

      equipmentIndex = 0;
    }

    // get the next equipment and raw data
    equipment = subEvent->GetEquipment(equipmentIndex++);
    if (!equipment) continue;
    AliRawData* rawData = equipment->GetRawData();
    if (!rawData) continue;
    position = (UChar_t*) rawData->GetBuffer();
    end = ((UChar_t*) rawData->GetBuffer()) + rawData->GetSize();

    // continue with the next sub event if no data left in the payload
    if (position >= end) continue;

    if (fRequireHeader) {
    // check that there are enough bytes left for the data header
      if (position + sizeof(AliRawDataHeader) > end) {
	result |= kErrNoDataHeader;
	continue;
      }

      // check consistency of data size in the header and in the equipment
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

  return result;
}
