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
// range of DDLs (minDDLID <= DDLID < maxDDLID).
// no selection is applied if a value < 0 is used.

  fSelectDetectorID = detectorID;
  fSelectMinDDLID = minDDLID;
  fSelectMaxDDLID = maxDDLID;
}

Bool_t AliRawReader::IsSelected() const
{
// apply the selection (if any)

  if (fSelectDetectorID >= 0) {
    if (!fMiniHeader) return kFALSE;
    if (fMiniHeader->fDetectorID != fSelectDetectorID) return kFALSE;
    if ((fSelectMinDDLID >= 0) && (fMiniHeader->fDDLID < fSelectMinDDLID))
      return kFALSE;
    if ((fSelectMaxDDLID >= 0) && (fMiniHeader->fDDLID >= fSelectMaxDDLID))
      return kFALSE;
  }
  return kTRUE;
}


Bool_t AliRawReader::CheckMiniHeader(AliMiniHeader* miniHeader) const
{
// check the magic number of the mini header

  if (!miniHeader) miniHeader = fMiniHeader;
  if (!miniHeader) return kFALSE;
  if ((miniHeader->fMagicWord[2] != 0x12) ||
      (miniHeader->fMagicWord[1] != 0x34) ||
      (miniHeader->fMagicWord[0] != 0x56)) {
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

