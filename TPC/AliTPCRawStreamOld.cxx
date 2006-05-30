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
/// This class provides access to TPC digits in raw data.
///
/// It loops over all TPC digits in the raw data given by the AliRawReader.
/// The Next method goes to the next digit. If there are no digits left
/// it returns kFALSE.
/// Several getters provide information about the current digit.
///
///////////////////////////////////////////////////////////////////////////////

#include <TSystem.h>

#include "AliTPCRawStreamOld.h"
#include "AliRawReader.h"
#include "AliLog.h"
#include "AliTPCAltroMapping.h"

ClassImp(AliTPCRawStreamOld)

//_____________________________________________________________________________
AliTPCRawStreamOld::AliTPCRawStreamOld(AliRawReader* rawReader) :
  AliAltroRawStreamOld(rawReader),
  fSector(-1),
  fPrevSector(-1),
  fRow(-1),
  fPrevRow(-1),
  fPad(-1),
  fPrevPad(-1)
{
  // create an object to read TPC raw digits

  SelectRawData(0);

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
AliTPCRawStreamOld::AliTPCRawStreamOld(const AliTPCRawStreamOld& stream) :
  AliAltroRawStreamOld(stream),
  fSector(-1),
  fPrevSector(-1),
  fRow(-1),
  fPrevRow(-1),
  fPad(-1),
  fPrevPad(-1)
{
  Fatal("AliTPCRawStreamOld", "copy constructor not implemented");
}

//_____________________________________________________________________________
AliTPCRawStreamOld& AliTPCRawStreamOld::operator = (const AliTPCRawStreamOld& 
					      /* stream */)
{
  Fatal("operator =", "assignment operator not implemented");
  return *this;
}

//_____________________________________________________________________________
AliTPCRawStreamOld::~AliTPCRawStreamOld()
{
// destructor

  for(Int_t i = 0; i < 6; i++) delete fMapping[i];
}

//_____________________________________________________________________________
void AliTPCRawStreamOld::Reset()
{
  // reset tpc raw stream params
  AliAltroRawStreamOld::Reset();
  fSector = fPrevSector = fRow = fPrevRow = fPad = fPrevPad = -1;
}

//_____________________________________________________________________________
Bool_t AliTPCRawStreamOld::Next()
{
  // Read next TPC signal
  // Apply the TPC altro mapping to get
  // the sector,pad-row and pad indeces
  fPrevSector = fSector;
  fPrevRow = fRow;
  fPrevPad = fPad;
  if (AliAltroRawStreamOld::Next()) {
    if (IsNewHWAddress())
      ApplyAltroMapping();
    return kTRUE;
  }
  else
    return kFALSE;
}

//_____________________________________________________________________________
void AliTPCRawStreamOld::ApplyAltroMapping()
{
  // Take the DDL index, load
  // the corresponding altro mapping
  // object and fill the sector,row and pad indeces
  Int_t ddlNumber = GetDDLNumber();
  Int_t patchIndex;
  if (ddlNumber < 72) {
    fSector = ddlNumber / 2;
    patchIndex = ddlNumber % 2;
  }
  else {
    fSector = (ddlNumber - 72) / 4 + 36;
    patchIndex = (ddlNumber - 72) % 4 + 2;
  }

  Short_t hwAddress = GetHWAddress();
  fRow = fMapping[patchIndex]->GetPadRow(hwAddress);
  fPad = fMapping[patchIndex]->GetPad(hwAddress);

}
