// @(#) $Id$

// Author: Constantin Loizides <mailto:loizides@ikf.uni-frankfurt.de>
//*-- Copyright &copy ALICE HLT Group

#include "AliHLTRootTypes.h"
#include "AliHLTStandardIncludes.h"
#include "AliHLTLogging.h"

#include "AliHLTDDLRawReader.h"

#if __GNUC__ >= 3
using namespace std;
#endif

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

/** \class AliHLTDDLRawReader
<pre>
//_____________________________________________________________
// AliHLTDDLRawReader (taken from the offline AliROOT code,
// original authors: D.Favretto and A.K.Mohanty)
//
// This is the base class for reading ddl raw data 
// and providing information about digits
</pre>
*/

ClassImp(AliHLTDDLRawReader)

AliHLTDDLRawReader::AliHLTDDLRawReader()
{
  fMiniHeader = NULL;
  fCount = 0;
  fSelectDetectorID = -1;
  fSelectMinDDLID = -1;
  fSelectMaxDDLID = -1;
}

AliHLTDDLRawReader::~AliHLTDDLRawReader()
{
}

void AliHLTDDLRawReader::Select(Int_t detectorID, Int_t minDDLID, Int_t maxDDLID)
{
  // read only data of the detector with the given ID and in the given
  // range of DDLs (minDDLID <= DDLID < maxDDLID).
  // no selection is applied if a value < 0 is used.

  fSelectDetectorID = detectorID;
  fSelectMinDDLID = minDDLID;
  fSelectMaxDDLID = maxDDLID;
}

Bool_t AliHLTDDLRawReader::IsSelected() const
{
  // apply the selection (if any)

  if (fSelectDetectorID >= 0) {
    if (fMiniHeader->fDetectorID != fSelectDetectorID) return kFALSE;
    if ((fSelectMinDDLID >= 0) && (fMiniHeader->fDDLID < fSelectMinDDLID))
      return kFALSE;
    if ((fSelectMaxDDLID >= 0) && (fMiniHeader->fDDLID >= fSelectMaxDDLID))
      return kFALSE;
  }
  return kTRUE;
}

Bool_t AliHLTDDLRawReader::CheckMiniHeader() const
{
  // check the magic number of the mini header

  if ((fMiniHeader->fMagicWord[2] != 0x12) ||
      (fMiniHeader->fMagicWord[1] != 0x34) ||
      (fMiniHeader->fMagicWord[0] != 0x56)) {
    LOG(AliHLTLog::kError,"AliHLTDDLRawReader::CheckMiniHeader","MH")
      <<"DDL mini header has wrong magic word!"<<ENDLOG;
    return kFALSE;
  }
  return kTRUE;
}

Bool_t AliHLTDDLRawReader::ReadNextInt(UInt_t& data)
{
  // reads the next 4 bytes at the current position
  // returns kFALSE if the data could not be read

  while (fCount == 0) {
    if (!ReadMiniHeader()) return kFALSE;
  }
  if (fCount < (Int_t) sizeof(data)) {
    LOG(AliHLTLog::kError,"AliHLTDDLRawReader::ReadNextInt","Data")
      <<AliHLTLog::kDec<<"Too few data left ("<<fCount<<") to read UInt_t!"<<ENDLOG;
    return kFALSE;
  }
  if (!ReadNext((UChar_t*) &data, sizeof(data))) {
    LOG(AliHLTLog::kError,"AliHLTDDLRawReader::ReadNextInt","Data")
      <<"Could not read data."<<ENDLOG;
    return kFALSE;
  }
  return kTRUE;
}

Bool_t AliHLTDDLRawReader::ReadNextShort(UShort_t& data)
{
  // reads the next 2 bytes at the current position
  // returns kFALSE if the data could not be read

  while (fCount == 0) {
    if (!ReadMiniHeader()) return kFALSE;
  }
  if (fCount < (Int_t) sizeof(data)) {
    LOG(AliHLTLog::kError,"AliHLTDDLRawReader::ReadNextShort","Data")
      <<AliHLTLog::kDec<<"Too few data left ("<<fCount<<") to read UShort_t!"<<ENDLOG;
    return kFALSE;
  }
  if (!ReadNext((UChar_t*) &data, sizeof(data))) {
    LOG(AliHLTLog::kError,"AliHLTDDLRawReader::ReadNextShort","Data")
      <<"Could not read data."<<ENDLOG;
    return kFALSE;
  }
  return kTRUE;
}

Bool_t AliHLTDDLRawReader::ReadNextChar(UChar_t& data)
{
  // reads the next 1 byte at the current stream position
  // returns kFALSE if the data could not be read

  while (fCount == 0) {
    if (!ReadMiniHeader()) return kFALSE;
  }
  if (!ReadNext((UChar_t*) &data, sizeof(data))) {
    LOG(AliHLTLog::kError,"AliHLTDDLRawReader::ReadNextChar","Data")
      <<"Could not read data."<<ENDLOG;
    return kFALSE;
  }
  return kTRUE;
}

