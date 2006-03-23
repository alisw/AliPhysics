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
/// This class provides access to TPC digits in raw data.
///
/// It loops over all TPC digits in the raw data given by the AliRawReader.
/// The Next method goes to the next digit. If there are no digits left
/// it returns kFALSE.
/// Several getters provide information about the current digit.
///
///////////////////////////////////////////////////////////////////////////////

#include <TSystem.h>

#include "AliTPCRawStream.h"
#include "AliRawReader.h"
#include "AliLog.h"
#include "AliTPCAltroMapping.h"

ClassImp(AliTPCRawStream)

//_____________________________________________________________________________
AliTPCRawStream::AliTPCRawStream(AliRawReader* rawReader) :
  AliAltroRawStream(rawReader)
{
  // create an object to read TPC raw digits

  fRawReader->Select(0);

  TString path = gSystem->Getenv("ALICE_ROOT");
  path += "/TPC/mapping/Patch";
  TString path2;
  for(Int_t i = 0; i < 6; i++) {
    path2 = path;
    path2 += i;
    path2 += ".data";
    fMapping[i] = new AliTPCAltroMapping(path2.Data());
  }

  fNoAltroMapping = kFALSE;
}

//_____________________________________________________________________________
AliTPCRawStream::AliTPCRawStream(const AliTPCRawStream& stream) :
  AliAltroRawStream(stream)
{
  Fatal("AliTPCRawStream", "copy constructor not implemented");
}

//_____________________________________________________________________________
AliTPCRawStream& AliTPCRawStream::operator = (const AliTPCRawStream& 
					      /* stream */)
{
  Fatal("operator =", "assignment operator not implemented");
  return *this;
}

//_____________________________________________________________________________
AliTPCRawStream::~AliTPCRawStream()
{
// destructor

  for(Int_t i = 0; i < 6; i++) delete fMapping[i];
}

//_____________________________________________________________________________
void AliTPCRawStream::Reset()
{
  // reset tpc raw stream params
  AliAltroRawStream::Reset();
}

//_____________________________________________________________________________
void AliTPCRawStream::ApplyAltroMapping()
{
  // Reads the DDL index, loads
  // the corresponding altro mapping
  // object and fills the sector,row and pad indeces
  Int_t ddlNumber = fRawReader->GetDDLID();
  Int_t patchIndex;
  if (ddlNumber < 72) {
    fSector = ddlNumber / 2;
    patchIndex = ddlNumber % 2;
  }
  else {
    fSector = (ddlNumber - 72) / 4 + 36;
    patchIndex = (ddlNumber - 72) % 4 + 2;
  }

  fRow = fMapping[patchIndex]->GetPadRow(fHWAddress);
  fPad = fMapping[patchIndex]->GetPad(fHWAddress);

}
