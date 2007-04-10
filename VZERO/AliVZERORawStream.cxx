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
/// This is a class for reading the VZERO DDL raw data
/// The format of the raw data corresponds to the one
/// implemented in AliVZEROBuffer class.
///
///////////////////////////////////////////////////////////////////////////////

#include "AliVZERORawStream.h"
#include "AliRawReader.h"
#include "AliLog.h"
#include "AliDAQ.h"

ClassImp(AliVZERORawStream)

//_____________________________________________________________________________
AliVZERORawStream::AliVZERORawStream(AliRawReader* rawReader) :
  fCell(-1),
  fADC(-1),
  fTime(-1),
  fPosition(-1),
  fRawReader(rawReader),
  fData(NULL)
{
  // create an object to read VZERO raw data
  //
  // select the raw data corresponding to
  // the VZERO detector id
  fRawReader->Reset();
  AliDebug(1,Form("Selecting raw data for detector %d",AliDAQ::DetectorID("VZERO")));
  fRawReader->Select("VZERO");
}

//_____________________________________________________________________________
AliVZERORawStream::~AliVZERORawStream()
{
  // destructor
}

//_____________________________________________________________________________
void AliVZERORawStream::Reset()
{
  // reset raw stream params

  fCell = fADC = fTime = fPosition = -1;

  if (fRawReader) fRawReader->Reset();
}

//_____________________________________________________________________________
Bool_t AliVZERORawStream::Next()
{
  // read next digit from the VZERO raw data stream
  // return kFALSE in case of error or no digits left

  if (fPosition < 0) {
    if (!fRawReader->ReadNextData(fData)) return kFALSE;
    if (fRawReader->GetDataSize()%3 != 0) {
      fRawReader->AddFatalErrorLog(kRawDataSizeErr,Form("size %d != n*12",fRawReader->GetDataSize()));
      AliWarning(Form("Wrong VZERO raw data size: %d, expected n*12 bytes!",fRawReader->GetDataSize()));
      return kFALSE;
    }
    fPosition = 0;
  }

  if (fPosition >= fRawReader->GetDataSize()) return kFALSE;

  fCell = GetNextWord();
  fADC  = GetNextWord();
  fTime = GetNextWord();

  return kTRUE;
}

//_____________________________________________________________________________
Int_t AliVZERORawStream::GetNextWord()
{
  // This method returns the next 32 bit word
  // inside the raw data payload.
  // The method is supposed to be endian (platform)
  // independent.
  if (!fData || fPosition < 0) AliFatal("Raw data payload buffer is not yet initialized !");

  Int_t word = 0;
  word |= fData[fPosition++];
  word |= fData[fPosition++] << 8;
  word |= fData[fPosition++] << 16;
  word |= fData[fPosition++] << 24;

  return word;
}
