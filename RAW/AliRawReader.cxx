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
/// This is the base class for reading raw data.
///
/// The derived classes, which operate on concrete raw data formats,
/// should implement
/// - ReadHeader to read the next (data/equipment) header
/// - ReadNextData to read the next raw data block (=1 DDL)
/// - ReadNext to read a given number of bytes
/// - several getters like GetType
///
/// Sequential access to the raw data is provided by the methods
/// ReadHeader, ReadNextData, ReadNextInt, ReadNextShort, ReadNextChar
///
/// If only data from a specific detector (and a given range of DDL numbers)
/// should be read, this can be achieved by the Select method.
/// Several getters provide information about the current event and the
/// current type of raw data.
///
///////////////////////////////////////////////////////////////////////////////

#include <Riostream.h>
#include "AliRawReader.h"
#include "AliDAQ.h"
#include "AliLog.h"

ClassImp(AliRawReader)


AliRawReader::AliRawReader() :
  fEquipmentIdsIn(NULL),
  fEquipmentIdsOut(NULL),
  fRequireHeader(kTRUE),
  fHeader(NULL),
  fCount(0),
  fSelectEquipmentType(-1),
  fSelectMinEquipmentId(-1),
  fSelectMaxEquipmentId(-1),
  fSkipInvalid(kFALSE),
  fSelectEventType(-1),
  fErrorCode(0),
  fEventNumber(-1),
  fErrorLogs("AliRawDataErrorLog",100)
{
// default constructor: initialize data members
}

Bool_t AliRawReader::LoadEquipmentIdsMap(const char *fileName)
{
  // Open the mapping file
  // and load the mapping data
  ifstream input(fileName);
  if (input.is_open()) {
    Warning("AliRawReader","Equipment ID mapping file is found !");
    const Int_t kMaxDDL = 256;
    fEquipmentIdsIn = new TArrayI(kMaxDDL);
    fEquipmentIdsOut = new TArrayI(kMaxDDL);
    Int_t equipIn, equipOut;
    Int_t nIds = 0;
    while (input >> equipIn >> equipOut) {
      if (nIds >= kMaxDDL) {
	Error("AliRawReader","Too many equipment Id mappings found ! Truncating the list !");
	break;
      }
      fEquipmentIdsIn->AddAt(equipIn,nIds); 
      fEquipmentIdsOut->AddAt(equipOut,nIds);
      nIds++;
    }
    fEquipmentIdsIn->Set(nIds);
    fEquipmentIdsOut->Set(nIds);
    input.close();
    return kTRUE;
  }
  else {
    Error("AliRawReader","equipment id map file is not found ! Skipping the mapping !");
    return kFALSE;
  }
}

AliRawReader::AliRawReader(const AliRawReader& rawReader) :
  TObject(rawReader),
  fEquipmentIdsIn(rawReader.fEquipmentIdsIn),
  fEquipmentIdsOut(rawReader.fEquipmentIdsOut),
  fRequireHeader(rawReader.fRequireHeader),
  fHeader(rawReader.fHeader),
  fCount(rawReader.fCount),
  fSelectEquipmentType(rawReader.fSelectEquipmentType),
  fSelectMinEquipmentId(rawReader.fSelectMinEquipmentId),
  fSelectMaxEquipmentId(rawReader.fSelectMaxEquipmentId),
  fSkipInvalid(rawReader.fSkipInvalid),
  fSelectEventType(rawReader.fSelectEventType),
  fErrorCode(0),
  fEventNumber(-1),
  fErrorLogs("AliRawDataErrorLog",100)
{
// copy constructor
}

AliRawReader& AliRawReader::operator = (const AliRawReader& rawReader)
{
// assignment operator
  fEquipmentIdsIn = rawReader.fEquipmentIdsIn;
  fEquipmentIdsOut = rawReader.fEquipmentIdsOut;

  fHeader = rawReader.fHeader;
  fCount = rawReader.fCount;

  fSelectEquipmentType = rawReader.fSelectEquipmentType;
  fSelectMinEquipmentId = rawReader.fSelectMinEquipmentId;
  fSelectMaxEquipmentId = rawReader.fSelectMaxEquipmentId;
  fSkipInvalid = rawReader.fSkipInvalid;
  fSelectEventType = rawReader.fSelectEventType;

  fErrorCode = rawReader.fErrorCode;

  fEventNumber = rawReader.fEventNumber;
  fErrorLogs = *((TClonesArray*)rawReader.fErrorLogs.Clone());

  return *this;
}

AliRawReader::~AliRawReader()
{
  // destructor
  // delete the mapping arrays if
  // initialized
  if (fEquipmentIdsIn) delete fEquipmentIdsIn;
  if (fEquipmentIdsOut) delete fEquipmentIdsOut;
}

Int_t AliRawReader::GetMappedEquipmentId() const
{
  if (!fEquipmentIdsIn || !fEquipmentIdsOut) {
    Error("AliRawReader","equipment Ids mapping is not initialized !");
    return GetEquipmentId();
  }
  Int_t equipmentId = GetEquipmentId();
  for(Int_t iId = 0; iId < fEquipmentIdsIn->GetSize(); iId++) {
    if (equipmentId == fEquipmentIdsIn->At(iId)) {
      equipmentId = fEquipmentIdsOut->At(iId);
      break;
    }
  }
  return equipmentId;
}

Int_t AliRawReader::GetDetectorID() const
{
  // Get the detector ID
  // The list of detector IDs
  // can be found in AliDAQ.h
  Int_t equipmentId;
  if (fEquipmentIdsIn && fEquipmentIdsIn)
    equipmentId = GetMappedEquipmentId();
  else
    equipmentId = GetEquipmentId();

  if (equipmentId >= 0) {
    Int_t ddlIndex;
    return AliDAQ::DetectorIDFromDdlID(equipmentId,ddlIndex);
  }
  else
    return -1;
}

Int_t AliRawReader::GetDDLID() const
{
  // Get the DDL ID (within one sub-detector)
  // The list of detector IDs
  // can be found in AliDAQ.h
  Int_t equipmentId;
  if (fEquipmentIdsIn && fEquipmentIdsIn)
    equipmentId = GetMappedEquipmentId();
  else
    equipmentId = GetEquipmentId();

  if (equipmentId >= 0) {
    Int_t ddlIndex;
    AliDAQ::DetectorIDFromDdlID(equipmentId,ddlIndex);
    return ddlIndex;
  }
  else
    return -1;
}

void AliRawReader::Select(const char *detectorName, Int_t minDDLID, Int_t maxDDLID)
{
// read only data of the detector with the given name and in the given
// range of DDLs (minDDLID <= DDLID <= maxDDLID).
// no selection is applied if a value < 0 is used.
  Int_t detectorID = AliDAQ::DetectorID(detectorName);
  if(detectorID >= 0)
    Select(detectorID,minDDLID,maxDDLID);
}

void AliRawReader::Select(Int_t detectorID, Int_t minDDLID, Int_t maxDDLID)
{
// read only data of the detector with the given ID and in the given
// range of DDLs (minDDLID <= DDLID <= maxDDLID).
// no selection is applied if a value < 0 is used.

  fSelectEquipmentType = -1;

  if (minDDLID < 0)
    fSelectMinEquipmentId = AliDAQ::DdlIDOffset(detectorID);
  else
    fSelectMinEquipmentId = AliDAQ::DdlID(detectorID,minDDLID);

  if (maxDDLID < 0)
    fSelectMaxEquipmentId = AliDAQ::DdlID(detectorID,AliDAQ::NumberOfDdls(detectorID)-1);
  else
    fSelectMaxEquipmentId = AliDAQ::DdlID(detectorID,maxDDLID);
}

void AliRawReader::SelectEquipment(Int_t equipmentType, 
				   Int_t minEquipmentId, Int_t maxEquipmentId)
{
// read only data of the equipment with the given type and in the given
// range of IDs (minEquipmentId <= EquipmentId <= maxEquipmentId).
// no selection is applied if a value < 0 is used.

  fSelectEquipmentType = equipmentType;
  fSelectMinEquipmentId = minEquipmentId;
  fSelectMaxEquipmentId = maxEquipmentId;
}

void AliRawReader::SelectEvents(Int_t type)
{
// read only events with the given type.
// no selection is applied if a value < 0 is used.

  fSelectEventType = type;
}

Bool_t AliRawReader::IsSelected() const
{
// apply the selection (if any)

  if (fSkipInvalid && !IsValid()) return kFALSE;

  if (fSelectEquipmentType >= 0)
    if (GetEquipmentType() != fSelectEquipmentType) return kFALSE;

  Int_t equipmentId;
  if (fEquipmentIdsIn && fEquipmentIdsIn)
    equipmentId = GetMappedEquipmentId();
  else
    equipmentId = GetEquipmentId();

  if ((fSelectMinEquipmentId >= 0) && 
      (equipmentId < fSelectMinEquipmentId))
    return kFALSE;
  if ((fSelectMaxEquipmentId >= 0) && 
      (equipmentId > fSelectMaxEquipmentId))
    return kFALSE;

  return kTRUE;
}

Bool_t AliRawReader::IsEventSelected() const
{
// apply the event selection (if any)

  if (fSelectEventType >= 0) {
    if (GetType() != (UInt_t) fSelectEventType) return kFALSE;
  }

  return kTRUE;
}


Bool_t AliRawReader::ReadNextInt(UInt_t& data)
{
// reads the next 4 bytes at the current position
// returns kFALSE if the data could not be read

  while (fCount == 0) {
    if (!ReadHeader()) return kFALSE;
  }
  if (fCount < (Int_t) sizeof(data)) {
    Error("ReadNextInt", 
	  "too few data left (%d bytes) to read an UInt_t!", fCount);
    return kFALSE;
  }
  if (!ReadNext((UChar_t*) &data, sizeof(data))) {
    Error("ReadNextInt", "could not read data!");
    return kFALSE;
  }
  return kTRUE;
}

Bool_t AliRawReader::ReadNextShort(UShort_t& data)
{
// reads the next 2 bytes at the current position
// returns kFALSE if the data could not be read

  while (fCount == 0) {
    if (!ReadHeader()) return kFALSE;
  }
  if (fCount < (Int_t) sizeof(data)) {
    Error("ReadNextShort", 
	  "too few data left (%d bytes) to read an UShort_t!", fCount);
    return kFALSE;
  }
  if (!ReadNext((UChar_t*) &data, sizeof(data))) {
    Error("ReadNextShort", "could not read data!");
    return kFALSE;
  }
  return kTRUE;
}

Bool_t AliRawReader::ReadNextChar(UChar_t& data)
{
// reads the next 1 byte at the current stream position
// returns kFALSE if the data could not be read

  while (fCount == 0) {
    if (!ReadHeader()) return kFALSE;
  }
  if (!ReadNext((UChar_t*) &data, sizeof(data))) {
    Error("ReadNextChar", "could not read data!");
    return kFALSE;
  }
  return kTRUE;
}


Int_t AliRawReader::CheckData() const
{
// check the consistency of the data
// derived classes should overwrite the default method which returns 0 (no err)

  return 0;
}


void AliRawReader::DumpData(Int_t limit)
{
// print the raw data
// if limit is not negative, only the first and last "limit" lines of raw data
// are printed

  Reset();
  if (!ReadHeader()) {
    Error("DumpData", "no header");
    return;
  }
  printf("header:\n"
	 " type = %d  run = %d  ", GetType(), GetRunNumber());
  if (GetEventId()) {
    printf("event = %8.8x %8.8x\n", GetEventId()[1], GetEventId()[0]);
  } else {
    printf("event = -------- --------\n");
  }
  if (GetTriggerPattern()) {
    printf(" trigger = %8.8x %8.8x  ",
	   GetTriggerPattern()[1], GetTriggerPattern()[0]);
  } else {
    printf(" trigger = -------- --------  ");
  }
  if (GetDetectorPattern()) {
    printf("detector = %8.8x\n", GetDetectorPattern()[0]);
  } else {
    printf("detector = --------\n");
  }
  if (GetAttributes()) {
    printf(" attributes = %8.8x %8.8x %8.8x  ",
	   GetAttributes()[2], GetAttributes()[1], GetAttributes()[0]);
  } else {
    printf(" attributes = -------- -------- --------  ");
  }
  printf("GDC = %d\n", GetGDCId());
  printf("\n");

  do {
    printf("-------------------------------------------------------------------------------\n");
    printf("LDC = %d\n", GetLDCId());

    printf("equipment:\n"
	   " size = %d  type = %d  id = %d\n",
	   GetEquipmentSize(), GetEquipmentType(), GetEquipmentId());
    if (GetEquipmentAttributes()) {
      printf(" attributes = %8.8x %8.8x %8.8x  ", GetEquipmentAttributes()[2],
	     GetEquipmentAttributes()[1], GetEquipmentAttributes()[0]);
    } else {
      printf(" attributes = -------- -------- --------  ");
    }
    printf("element size = %d\n", GetEquipmentElementSize());

    printf("data header:\n"
	   " size = %d  version = %d  valid = %d  compression = %d\n",
	   GetDataSize(), GetVersion(), IsValid(), IsCompressed());

    printf("\n");
    if (limit == 0) continue;

    Int_t size = GetDataSize();
    char line[70];
    for (Int_t i = 0; i < 70; i++) line[i] = ' ';
    line[69] = '\0';
    Int_t pos = 0;
    Int_t max = 16;
    UChar_t byte;

    for (Int_t n = 0; n < size; n++) {
      if (!ReadNextChar(byte)) {
	Error("DumpData", "couldn't read byte number %d\n", n);
	break;
      }
      if (pos >= max) {
	printf("%8.8x  %s\n", n-pos, line);
	for (Int_t i = 0; i < 70; i++) line[i] = ' ';
	line[69] = '\0';
	pos = 0;
	if ((limit > 0) && (n/max == limit)) {
	  Int_t nContinue = ((size-1)/max+1-limit) * max;
	  if (nContinue > n) {
	    printf(" [skipping %d bytes]\n", nContinue-n);
	    n = nContinue-1;
	    continue;
	  }
	}
      }
      Int_t offset = pos/4;
      if ((byte > 0x20) && (byte < 0x7f)) {
	line[pos+offset] = byte;
      } else {
	line[pos+offset] = '.';
      }
      char hex[3];
      sprintf(hex, "%2.2x", byte);
      line[max+max/4+3+2*pos+offset] = hex[0];
      line[max+max/4+4+2*pos+offset] = hex[1];
      pos++;
    }

    if (pos > 0) printf("%8.8x  %s\n", size-pos, line);
    printf("\n");
	   
  } while (ReadHeader());
}

void AliRawReader::AddErrorLog(AliRawDataErrorLog::ERawDataErrorType type,
			       const char *message)
{
  // Add a raw data error message to the list
  // of raw-data decoding errors
  if (fEventNumber < 0) {
    AliError("No events have read so far! Impossible to add a raw data error log!");
    return;
  }
  Int_t ddlId = GetDDLID();
  if (ddlId < 0) {
    AliError("No ddl raw data have been read so far! Impossible to add a raw data error log!");
    return;
  }

  new (fErrorLogs[fErrorLogs.GetEntriesFast()])
    AliRawDataErrorLog(fEventNumber,
		       ddlId,
		       type,
		       message);
}
