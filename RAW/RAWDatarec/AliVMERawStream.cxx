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
/// This class provides access to VME data in test beam raw data.
///
///////////////////////////////////////////////////////////////////////////////

#include "AliVMERawStream.h"
#include "AliRawReader.h"

ClassImp(AliVMERawStream)



AliVMERawStream::AliVMERawStream(AliRawReader* rawReader) :
  fRawReader(rawReader),
  fData(0),
  fNChannels(-1),
  fBlock(0),
  fNumber(0),
  fChannel(0),
  fValue(0),
  fTime(0),
  fTimeMuSec(0)
{
// create an object to read VME raw digits

  fRawReader = rawReader;

  ReadTDC();
  ReadTime();

  fRawReader->Reset();
  fRawReader->SelectEquipment(551, 38, 38);
}

Bool_t AliVMERawStream::Next()
{
// read the next raw digit
// returns kFALSE if there is no digit left

  // V551 string
  if (fNChannels == -1) {
    if (!fRawReader->ReadNextInt(fData)) return kFALSE;
    if (!CheckString("V551")) return kFALSE;
    fNChannels = 0;
  }

  while (fNChannels == 0) {
    // V550 or v551 string
    if (!fRawReader->ReadNextInt(fData)) {
      Error("Next", "incomplete equipment");
      return kFALSE;
    }
    // check for v551 string (end of data)
    const char* v551 = "v551";
    if (fData == *((UInt_t*) v551)) return kFALSE;
    if (!CheckString("V550")) return kFALSE;

    // block
    if (!fRawReader->ReadNextShort(fBlock)) {
      Error("Next", "incomplete equipment");
      return kFALSE;
    }

    // serial number
    if (!fRawReader->ReadNextShort(fNumber)) {
      Error("Next", "incomplete equipment");
      return kFALSE;
    }

    // number of channels
    if (!fRawReader->ReadNextInt((UInt_t&) fNChannels)) {
      Error("Next", "incomplete equipment");
      return kFALSE;
    }
  }

  if (!fRawReader->ReadNextInt(fData)) {
    Error("Next", "incomplete equipment");
    return kFALSE;
  }
  fChannel = (fData >> 12) & 0x03ff;
  fValue = fData & 0x0fff;
  fNChannels--;

  return kTRUE;
}



Bool_t AliVMERawStream::CheckString(const char* str) const
{
// check fData to be equal to the given string

  if (fData != *((UInt_t*) str)) {
    char strData[5];
    memcpy(strData, &fData, 4);
    strData[4] = 0;
    Error("CheckString", "invalid %s string (%s)", str, strData);
    return kFALSE;
  }
  return kTRUE;
}

Bool_t AliVMERawStream::ReadTDC()
{
// read the TDC information

  fRawReader->Reset();
  fRawReader->SelectEquipment(775, 72, 72);

  // V775 string
  if (!fRawReader->ReadNextInt(fData)) return kFALSE;
  if (!CheckString("V775")) return kFALSE;

  // header
  if (!fRawReader->ReadNextInt(fData)) {
    Error("ReadTDC", "incomplete TDC equipment");
    return kFALSE;
  }
  if ((fData & 0x02000000) == 0) {
    Error("ReadTDC", "invalid header: 0x%x", fData);
    return kFALSE;
  }

  // check the array size
  Int_t nTDC = fRawReader->GetDataSize() / 4 - 4; // - V775,header,counter,v775
  if ( nTDC != 3 ) {
    Error("ReadTDC", "wrong number of TDC channels: %d", nTDC);
    return kFALSE;
  }

  // TDC data
  for (Int_t i = 0; i < fgkNTDC; i++) {
    if (!fRawReader->ReadNextInt(fData)) {
      Error("ReadTDC", "incomplete TDC equipment");
      return kFALSE;
    }
    if (fData & 0x07000000) {
      Warning("ReadTDC", "bad TDC data: %x", fData);
    }
    if ((fData & 0x00004000) == 0) {
      Warning("ReadTDC", "TDC data not valid: %x", fData);
    }
    if (fData & 0x00002000) {
      Warning("ReadTDC", "TDC data underflow: %x", fData);
    }
    if (fData & 0x00001000) {
      Warning("ReadTDC", "TDC data overflow: %x", fData);
    }
    fTDCChannel[i]  = (fData >> 16) & 0x1f;
    fTDCValue[i] = fData & 0x0fff;
  }

  // counter
  if (!fRawReader->ReadNextInt(fData)) {
    Error("ReadTDC", "incomplete TDC equipment");
    return kFALSE;
  }
  if ((fData & 0x04000000) == 0) {
    Error("ReadTDC", "invalid counter: 0x%x", fData);
    return kFALSE;
  }

  // v775 string
  if (!fRawReader->ReadNextInt(fData)) {
    Error("ReadTDC", "incomplete TDC equipment");
    return kFALSE;
  }
  if (!CheckString("v775")) return kFALSE;

  return kTRUE;
}

Bool_t AliVMERawStream::ReadTime()
{
// read the time information

  fRawReader->Reset();
  fRawReader->SelectEquipment(1970, 0x12345678, 0x12345678);

  // TIME string
  if (!fRawReader->ReadNextInt(fData)) return kFALSE;
  if (!CheckString("TIME")) return kFALSE;

  // time value
  if (!fRawReader->ReadNextInt(fTime)) {
    Error("ReadTime", "incomplete time equipment");
    return kFALSE;
  }

  // micro seconds value
  if (!fRawReader->ReadNextInt(fTimeMuSec)) {
    Error("ReadTime", "incomplete time equipment");
    return kFALSE;
  }

  // time string
  if (!fRawReader->ReadNextInt(fData)) {
    Error("ReadTime", "incomplete time equipment");
    return kFALSE;
  }
  if (!CheckString("time")) return kFALSE;

  return kTRUE;
}
