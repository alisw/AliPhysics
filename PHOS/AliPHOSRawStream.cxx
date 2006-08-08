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
/// This class provides access to PHOS digits in raw data.
///
/// It loops over all PHOS digits in the raw data given by the AliRawReader.
/// The Next method goes to the next digit. If there are no digits left
/// it returns kFALSE.
/// Several getters provide information about the current digit.
/// usage: 
/// root > AliRawReaderFile rawReader ; 
/// root > AliPHOSRawStream input(&rawReader) ; 
/// root > while (input.Next()) ..... 
///////////////////////////////////////////////////////////////////////////////

#include <TString.h>
#include <TSystem.h>

#include "AliPHOSRawStream.h"
#include "AliRawReader.h"
#include "AliPHOSAltroMapping.h"

ClassImp(AliPHOSRawStream)


//_____________________________________________________________________________
AliPHOSRawStream::AliPHOSRawStream(AliRawReader* rawReader) :
  AliAltroRawStream(rawReader),
  fModule(-1),
  fPrevModule(-1),
  fRow(-1),
  fPrevRow(-1),
  fColumn(-1),
  fPrevColumn(-1),
  fGain(0)
{
// create an object to read PHOS raw digits

  SelectRawData("PHOS");

  TString path = gSystem->Getenv("ALICE_ROOT");
  path += "/PHOS/mapping/RCU";
  TString path2;
  for(Int_t i = 0; i < 4; i++) {
    path2 = path;
    path2 += i;
    path2 += ".data";
    fMapping[i] = new AliPHOSAltroMapping(path2.Data());
  }

  SetNoAltroMapping(kFALSE);
}

//_____________________________________________________________________________
AliPHOSRawStream::AliPHOSRawStream(const AliPHOSRawStream& stream) :
  AliAltroRawStream(stream),
  fModule(-1),
  fPrevModule(-1),
  fRow(-1),
  fPrevRow(-1),
  fColumn(-1),
  fPrevColumn(-1),
  fGain(0)
{  
  Fatal("AliPHOSRawStream", "copy constructor not implemented");
}

//_____________________________________________________________________________
AliPHOSRawStream& AliPHOSRawStream::operator = (const AliPHOSRawStream& 
					      /* stream */)
{
  Fatal("operator =", "assignment operator not implemented");
  return *this;
}

//_____________________________________________________________________________
AliPHOSRawStream::~AliPHOSRawStream()
{
// destructor

  for(Int_t i = 0; i < 4; i++) delete fMapping[i];
}

//_____________________________________________________________________________
void AliPHOSRawStream::Reset()
{
  // reset phos raw stream params
  AliAltroRawStream::Reset();
  fModule = fPrevModule = fRow = fPrevRow = fColumn = fPrevColumn = -1;
  fGain = 0;
}

//_____________________________________________________________________________
Bool_t AliPHOSRawStream::Next()
{
  // Read next PHOS signal
  // Apply the PHOS altro mapping to get
  // the module,row and column indeces
  fPrevModule = fModule;
  fPrevRow = fRow;
  fPrevColumn = fColumn;
  if (AliAltroRawStream::Next()) {
    if (IsNewHWAddress())
      ApplyAltroMapping();
    return kTRUE;
  }
  else
    return kFALSE;
}

//_____________________________________________________________________________
void AliPHOSRawStream::ApplyAltroMapping()
{
  // Take the DDL index, load
  // the corresponding altro mapping
  // object and fill the sector,row and pad indeces
  Int_t ddlNumber = GetDDLNumber();
  fModule = ddlNumber / 4;

  Int_t rcuIndex = ddlNumber % 4;

  Short_t hwAddress = GetHWAddress();
  fRow = fMapping[rcuIndex]->GetPadRow(hwAddress);
  fColumn = fMapping[rcuIndex]->GetPad(hwAddress);
  fGain = fMapping[rcuIndex]->GetSector(hwAddress);

}
