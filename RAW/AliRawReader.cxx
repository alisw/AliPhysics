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
// This is the base class for reading raw data and providing
// information about digits
//
// The derived classes, which operate on concrete raw data formats,
// should implement
// - ReadHeader to read the next (mini or equipment) header
// - ReadNextData to read the next raw data block (=1 DDL)
// - ReadNext to read a given number of bytes
// - several getters like GetType
//
// Sequential access to the raw data is provided by the methods
// ReadHeader, ReadNextData, ReadNextInt, ReadNextShort, ReadNextChar
//
// If only data from a specific detector (and a given range of DDL numbers)
// should be read, this can be achieved by the Select method.
// Several getter provide information about the current event and the
// current type of raw data.
//
///////////////////////////////////////////////////////////////////////////////

#include "AliRawReader.h"


ClassImp(AliRawReader)


AliRawReader::AliRawReader() :
  fMiniHeader(NULL),
  fCount(0),
  fSelectDetectorID(-1),
  fSelectMinDDLID(-1),
  fSelectMaxDDLID(-1),
  fSelectEquipmentType(-1),
  fSelectMinEquipmentId(-1),
  fSelectMaxEquipmentId(-1),
  fErrorCode(0)
{
// default constructor: initialize data members

}

AliRawReader::AliRawReader(const AliRawReader& rawReader) :
  TObject(rawReader),
  fMiniHeader(rawReader.fMiniHeader),
  fCount(rawReader.fCount),
  fSelectDetectorID(rawReader.fSelectDetectorID),
  fSelectMinDDLID(rawReader.fSelectMinDDLID),
  fSelectMaxDDLID(rawReader.fSelectMaxDDLID),
  fSelectEquipmentType(rawReader.fSelectEquipmentType),
  fSelectMinEquipmentId(rawReader.fSelectMinEquipmentId),
  fSelectMaxEquipmentId(rawReader.fSelectMaxEquipmentId),
  fErrorCode(0)
{
// copy constructor

}

AliRawReader& AliRawReader::operator = (const AliRawReader& rawReader)
{
// assignment operator

  fMiniHeader = rawReader.fMiniHeader;
  fCount = rawReader.fCount;

  fSelectDetectorID = rawReader.fSelectDetectorID;
  fSelectMinDDLID = rawReader.fSelectMinDDLID;
  fSelectMaxDDLID = rawReader.fSelectMaxDDLID;

  fErrorCode = rawReader.fErrorCode;

  return *this;
}


void AliRawReader::Select(Int_t detectorID, Int_t minDDLID, Int_t maxDDLID)
{
// read only data of the detector with the given ID and in the given
// range of DDLs (minDDLID <= DDLID <= maxDDLID).
// no selection is applied if a value < 0 is used.

  fSelectDetectorID = detectorID;
  fSelectMinDDLID = minDDLID;
  fSelectMaxDDLID = maxDDLID;
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

Bool_t AliRawReader::IsSelected() const
{
// apply the selection (if any)

  if (fSelectDetectorID >= 0) {
    if (!fMiniHeader) return kFALSE;
    if (fMiniHeader->fDetectorID != fSelectDetectorID) return kFALSE;
    if ((fSelectMinDDLID >= 0) && (fMiniHeader->fDDLID < fSelectMinDDLID))
      return kFALSE;
    if ((fSelectMaxDDLID >= 0) && (fMiniHeader->fDDLID > fSelectMaxDDLID))
      return kFALSE;
  }

  if (fSelectEquipmentType >= 0) {
    if (GetEquipmentType() != fSelectEquipmentType) return kFALSE;
    if ((fSelectMinEquipmentId >= 0) && 
	(GetEquipmentId() < fSelectMinEquipmentId))
      return kFALSE;
    if ((fSelectMaxEquipmentId >= 0) && 
	(GetEquipmentId() > fSelectMaxEquipmentId))
      return kFALSE;
  }

  return kTRUE;
}


Bool_t AliRawReader::CheckMiniHeader(AliMiniHeader* miniHeader) const
{
// check the magic number of the mini header

  if (!miniHeader) miniHeader = fMiniHeader;
  if (!miniHeader) return kFALSE;
  UInt_t magicWord = miniHeader->fMagicWord[2]*65536 +
    miniHeader->fMagicWord[1]*256 + miniHeader->fMagicWord[0];
  if (magicWord != AliMiniHeader::kMagicWord) {
    return kFALSE;
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

    printf("mini header:\n"
	   " size = %d  detector = %d  DDL = %d\n"
	   " version = %d  compression = %d\n",
	   GetDataSize(), GetDetectorID(), GetDDLID(),
	   GetVersion(), IsCompressed());

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
