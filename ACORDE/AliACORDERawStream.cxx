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
//                                                                           //
//  Reads ACORDE DDL raw data from raw data stream                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliACORDERawStream.h"
#include "AliRawReader.h"
#include "AliLog.h"
#include "AliDAQ.h"
#include "AliRawReaderRoot.h"

ClassImp(AliACORDERawStream)

//_____________________________________________________________________________
AliACORDERawStream::AliACORDERawStream(AliRawReader* rawReader) :
  fRawReader(rawReader),
  fPosition(-1),
  fData(NULL),
  fDataSize(0)
{
  //
  // Create an object to read ACORDE raw data
  //
  // Created:      04 Feb 2008  Mario Sitta
  //

  fWord[0] = fWord[1] = fWord[2] = fWord[3] = 0;

  // Select the raw data corresponding to the ACORDE detector id
//  fRawReader->Reset();
  AliDebug(1,Form("Selecting raw data for detector %d",AliDAQ::DetectorID("ACORDE")));
  fRawReader->Select("ACORDE");

}

//_____________________________________________________________________________
AliACORDERawStream::AliACORDERawStream(const AliACORDERawStream &r) :
  TObject(),
  fRawReader(r.fRawReader),
  fPosition(-1),
  fData(NULL),
  fDataSize(0)
{
  // Simple copy constructor
  ((AliACORDERawStream &) r).Copy(*this);
}

//_____________________________________________________________________________
AliACORDERawStream::~AliACORDERawStream()
{
  // Default destructor
}

//_____________________________________________________________________________
AliACORDERawStream &AliACORDERawStream::operator=(const AliACORDERawStream &r)
{
  // Simple operator=
  if (this != &r)  ((AliACORDERawStream &) r).Copy(*this);
  return *this;
}

//_____________________________________________________________________________
void AliACORDERawStream::Reset()
{
  //
  // Reset the raw stream parameters
  //
  // Input:
  //
  // Output:
  //
  // Created:      04 Feb 2008  Mario Sitta
  //

  fPosition = -1;
  fData = NULL;

  if (fRawReader) fRawReader->Reset();
}

//_____________________________________________________________________________
Bool_t AliACORDERawStream::Next()
{
  //
  // Read next digit from the ACORDE raw data stream;
  // return kFALSE in case of error or no digits left
  //
  // Input:
  //
  // Output:
  //
  // Created:      04 Feb 2008  Mario Sitta
  //

  if (fPosition >= 0) return kFALSE;

  if (!fRawReader->ReadNextData(fData)) return kFALSE;
  if (fRawReader->GetDataSize() == 0) return kFALSE;

  fDataSize = fRawReader->GetDataSize();
  if (fDataSize != 16) {
    fRawReader->AddFatalErrorLog(kRawDataSizeErr,Form("size %d != 16",fDataSize));
    AliWarning(Form("Wrong ACORDE raw data size: %d, expected 16 bytes!",fDataSize));
    return kFALSE;
  }

  fPosition = 0;

  for (Int_t i=0; i<4; i++)
    fWord[i] = GetNextWord();

  return kTRUE;
}

//_____________________________________________________________________________
UInt_t AliACORDERawStream::GetWord(Int_t index) const
{
  //
  // Returns the ``index'' word from ACORDE raw data.
  //
  // Input:
  //         index : the index of the requested word
  // Output:
  //         word  : the 32 bit ``index'' word
  //
  // Created:      12 Feb 2008  Mario Sitta
  //

  if (index < 0 || index > 3) {
    AliWarning(Form("Wrong word index %d, returning 0",index));
    return 0;
  } else {
    return fWord[index];
  }
  
}

//_____________________________________________________________________________
UInt_t AliACORDERawStream::GetNextWord()
{
  //
  // Returns the next 32 bit word inside the raw data payload.
  // The method is supposed to be endian (platform) independent.
  //
  // Input:
  //
  // Output:
  //         word : a 32 bit word containing the data
  //
  // Created:      04 Feb 2008  Mario Sitta
  //

  if (!fData || fPosition < 0)
    AliFatal("Raw data payload buffer is not yet initialized !");

  UInt_t word = 0;
  word |= fData[fPosition++];
  word |= fData[fPosition++] << 8;
  word |= fData[fPosition++] << 16;
  word |= fData[fPosition++] << 24;

  return word;
}

//_____________________________________________________________________________

Int_t AliACORDERawStream::GetNEvents(char* fileName) 
{
	// Returns the Total Number of Events recorded by ACORDE 
	// Note: it may be a better way to do it !!
	// Input: fileName to Analyze
	// Output: Number of Total Events (fNEvents) in fileName
	// Created: 25 March 2008
	// Author: Mario Rodriguez Cahuantzi <mrodrigu@mail.cern.ch>
	
	AliRawReader* rCount = new AliRawReaderRoot(fileName);
	Int_t DyM=0;
	Int_t fNEvents=0;
	while(DyM==0)
  	{
  	if (!rCount->NextEvent()) DyM=1;
	else fNEvents++;
  	}
	delete rCount;
	return fNEvents;
}

//____________________________________________________________________________
