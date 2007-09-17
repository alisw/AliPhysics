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
/// This class provides fast access to Altro raw data.
///
/// It loops over all Altro payload in the raw data given by the AliRawReader.
///
///
///////////////////////////////////////////////////////////////////////////////

#include "AliAltroRawStreamFast.h"
#include "AliRawReader.h"
#include "AliLog.h"

ClassImp(AliAltroRawStreamFast)


//_____________________________________________________________________________
AliAltroRawStreamFast::AliAltroRawStreamFast(AliRawReader* rawReader) :
  fDDLNumber(-1),
  fPrevDDLNumber(-1),
  fDecoder(),
  fData(),
  fBunch(),
  fRawReader(rawReader)
{
  // Default constructor
  fRawReader->Reset();
}

AliAltroRawStreamFast::~AliAltroRawStreamFast()
{
  // Destructor
}

//_____________________________________________________________________________
void AliAltroRawStreamFast::SelectRawData(Int_t detId)
{
  // Select the raw data for specific
  // detector id
  AliDebug(1,Form("Selecting raw data for detector %d",detId));
  fRawReader->Select(detId);
}

//_____________________________________________________________________________
void AliAltroRawStreamFast::SelectRawData(const char *detName)
{
  // Select the raw data for specific
  // detector name
  AliDebug(1,Form("Selecting raw data for detector %s",detName));
  fRawReader->Select(detName);
}
//_____________________________________________________________________________
void AliAltroRawStreamFast::Reset()
{
    // reset altro raw stream params

    fDDLNumber = fPrevDDLNumber = -1;
    if (fRawReader) fRawReader->Reset();

}
//_____________________________________________________________________________
Bool_t AliAltroRawStreamFast::NextDDL()
{
  // Read next DDL raw-data payload and
  // runs the altro decoder over the payload

  UChar_t *dataPtr = NULL;
  do {
    if (!fRawReader->ReadNextData(dataPtr)) return kFALSE;
  } while (fRawReader->GetDataSize() == 0);
  fPrevDDLNumber = fDDLNumber;
  fDDLNumber = fRawReader->GetDDLID();

  // Temporary solution while Per Thomas is
  // changing the decoder code
  dataPtr -= sizeof(AliRawDataHeader);
  Int_t length = fRawReader->GetDataSize() + sizeof(AliRawDataHeader);

  fDecoder.SetMemory(dataPtr, length);
  fDecoder.Decode();

  return kTRUE;
}
//_____________________________________________________________________________
Bool_t AliAltroRawStreamFast::NextChannel()
{
  // Get the data for the next altro channel
  do {
    if (!fDecoder.NextChannel(&fData)) return kFALSE;
  }  while (fData.GetDataSize() == 0);
  return kTRUE;
}
//_____________________________________________________________________________
Bool_t AliAltroRawStreamFast::NextBunch()
{
  // Get the data for the next altro bunch
  return fData.NextBunch(&fBunch);
}
