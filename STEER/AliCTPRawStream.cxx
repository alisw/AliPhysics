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
///
/// This class provides access to CTP DDL raw data.
///
/// The raw data format is taken form the trigger TDR.
/// The meaning of the trigger class and cluster masks
/// are given in the trigger description file (in /data)
/// and in the AliCentralTrigger class.
///
///////////////////////////////////////////////////////////////////////////////

#include "AliCTPRawStream.h"
#include "AliRawReader.h"
#include "AliLog.h"

ClassImp(AliCTPRawStream)

//_____________________________________________________________________________
AliCTPRawStream::AliCTPRawStream(AliRawReader* rawReader) :
  fClassMask(0),
  fClusterMask(0),
  fRawReader(rawReader)
{
  // create an object to read CTP raw data
  //
  // select the raw data corresponding to
  // the CTP detector id
  fRawReader->Reset();
  AliDebug(1,Form("Selecting raw data for detector %d",kCTPIndex));
  fRawReader->Select(kCTPIndex);
}

//_____________________________________________________________________________
AliCTPRawStream::AliCTPRawStream(const AliCTPRawStream& stream) :
  TObject(stream),
  fClassMask(0),
  fClusterMask(0),
  fRawReader(NULL)
{
  // Copy constructor
  AliFatal("Copy constructor not implemented");
}

//_____________________________________________________________________________
AliCTPRawStream& AliCTPRawStream::operator = (const AliCTPRawStream& 
					      /* stream */)
{
  AliFatal("Assignment operator not implemented");
  return *this;
}

//_____________________________________________________________________________
AliCTPRawStream::~AliCTPRawStream()
{
  // destructor
}

//_____________________________________________________________________________
void AliCTPRawStream::Reset()
{
  // reset raw stream params

  fClassMask = fClusterMask = 0;

  if (fRawReader) fRawReader->Reset();
}

//_____________________________________________________________________________
Bool_t AliCTPRawStream::Next()
{
  // read the whole CTP raw data stream
  // return kFALSE in case of error

  UChar_t *data = NULL;

  if (!fRawReader->ReadNextData(data)) return kFALSE;

  if (fRawReader->GetDataSize() != 32) {
    AliError(Form("Wrong CTP raw data size: %d",fRawReader->GetDataSize()));
    return kFALSE;
  }

  fClusterMask = data[12] >> 2;

  fClassMask =  ((ULong64_t)data[12] & 0x3) << 48;

  fClassMask |= (ULong64_t)data[16] << 36;
  fClassMask |= ((ULong64_t)data[17] & 0xF) << 44;

  fClassMask |= (ULong64_t)data[20] << 24;
  fClassMask |= ((ULong64_t)data[21] & 0xF) << 32;

  fClassMask |= (ULong64_t)data[24] << 12;
  fClassMask |= ((ULong64_t)data[25] & 0xF) << 20;

  fClassMask |= (ULong64_t)data[28];
  fClassMask |= ((ULong64_t)data[29] & 0xF) << 8;

  return kTRUE;
}

