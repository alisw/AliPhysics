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

// Authors:
//	 Mario Rodriguez Cahuantzi <mrodrigu@mail.cern.ch>
//	 Arturo Fernandez Tellez <afernan@mail.cern.ch>
//____________________________________________________________________
//                                                                          
// ACORDE 
// Class for reading ACORDE RAW data in TOF data format
//
#include "AliACORDERawReader.h"
#include "AliBitPacking.h"
#include "TBits.h"

#include <Riostream.h>
#include "TMath.h"
#include "TH1F.h"
#include "TArrayI.h"
#include "AliLog.h"
 
ClassImp(AliACORDERawReader)
  
AliACORDERawReader::AliACORDERawReader (AliRawReader *rawReader, Bool_t isOnline):
       TNamed("ACORDERawReader","read raw ACORDE data"),
       fRawReader(rawReader),
       fData(NULL),
       fPosition(0),
       fIsOnline(isOnline),
	fDataSize(0)
{

	fWord[0] = fWord[1] = fWord[2] = fWord[3] = 0;

}
//_____________________________________________________________________________
 AliACORDERawReader::~AliACORDERawReader ()
{
}
//_____________________________________________________________________________
Bool_t  AliACORDERawReader::Next()
{

// Read next digit from the ACORDE raw data stream;
// return kFALSE in case of error or no digits left

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
Int_t AliACORDERawReader::GetPosition()
{
  // Sets the position in the
  // input stream
  if (((fRawReader->GetDataSize() * 8) % 32) != 0)
    AliFatal(Form("Incorrect raw data size ! %d words are found !",fRawReader->GetDataSize()));
  return (fRawReader->GetDataSize() * 8) / 32;
}
//_____________________________________________________________________________
UInt_t AliACORDERawReader::GetNextWord()
{

  // Returns the next 32 bit word inside the raw data payload.
  // The method is supposed to be endian (platform) independent.


 if (!fData || fPosition < 0)
    AliFatal("Raw data payload buffer is not yet initialized !");

  UInt_t word = 0;
  word |= fData[fPosition++];
  word |= fData[fPosition++] << 8;
  word |= fData[fPosition++] << 16;
  word |= fData[fPosition++] << 24;

  return word;

 
}

