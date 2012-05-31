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

/* $Id: AliTPCRawStreamV3.cxx 22777 2007-12-05 17:37:33Z marian $ */

///////////////////////////////////////////////////////////////////////////////
///
/// This class provides access to TPC digits in raw data.
///
/// It loops over all TPC digits in the raw data given by the AliRawReader.
/// The NextChannel method loads the data for the next pad. If there is no pad left
/// it returns kFALSE.
///
///////////////////////////////////////////////////////////////////////////////

#include <TSystem.h>

#include "AliTPCRawStreamV3.h"
#include "AliRawReader.h"
#include "AliLog.h"
#include "AliTPCAltroMapping.h"

ClassImp(AliTPCRawStreamV3)

//_____________________________________________________________________________
AliTPCRawStreamV3::AliTPCRawStreamV3(AliRawReader* rawReader, AliAltroMapping **mapping) :
  AliAltroRawStreamV3(rawReader),
  fSector(-1),
  fPrevSector(-1),
  fRow(-1),
  fPrevRow(-1),
  fPad(-1),
  fPrevPad(-1),
  fPatchIndex(-1),
  fIsMapOwner(kFALSE)
{
  // create an object to read TPC raw digits

  SelectRawData("TPC");

  if (mapping == NULL) {
    TString path = gSystem->Getenv("ALICE_ROOT");
    path += "/TPC/mapping/Patch";
    TString path2;
    for(Int_t i = 0; i < 6; i++) {
      path2 = path;
      path2 += i;
      path2 += ".data";
      fMapping[i] = new AliTPCAltroMapping(path2.Data());
    }
    fIsMapOwner = kTRUE;
  }
  else {
    for(Int_t i = 0; i < 6; i++)
      fMapping[i] = mapping[i];
  }


  //fNoAltroMapping = kFALSE;
}
//_____________________________________________________________________________
AliTPCRawStreamV3::~AliTPCRawStreamV3()
{
// destructor

  if (fIsMapOwner)
    for(Int_t i = 0; i < 6; i++) delete fMapping[i];
}

//_____________________________________________________________________________
void AliTPCRawStreamV3::Reset()
{
  // reset tpc raw stream params
  AliAltroRawStreamV3::Reset();
  fSector = fPrevSector = fRow = fPrevRow = fPad = fPrevPad = fPatchIndex = -1;
}

//_____________________________________________________________________________
Bool_t AliTPCRawStreamV3::NextChannel()
{
  // Read next TPC Channel
  // Apply the TPC altro mapping to get
  // the pad-row and pad indeces

/*
  fPrevSector = fSector;
  fPrevRow = fRow;
  fPrevPad = fPad;
  if (AliAltroRawStreamV3::NextChannel()) {
      //    if (IsNewHWAddress())
    if ( GetHWAddress() > -1 )
      ApplyAltroMapping();
    return kTRUE;
  }
  else
    return kFALSE;
  */
  
  fPrevRow = fRow;
  fPrevPad = fPad;
  fRow = -1;
  fPad = -1;
  if (!AliAltroRawStreamV3::NextChannel()) return kFALSE;

  Short_t hwAddress = GetHWAddress();
  if (hwAddress>-1){
    fRow = fMapping[fPatchIndex]->GetPadRow(hwAddress);
    fPad = fMapping[fPatchIndex]->GetPad(hwAddress);
  }
  return kTRUE;
}
//_____________________________________________________________________________
Bool_t AliTPCRawStreamV3::NextDDL()
{
  // Take the DDL index,
  // calculate the patch number and
  // set the sector number
//   return AliAltroRawStreamV3::NextDDL();
  fPrevSector = fSector;
  fSector     = -1;
  if (!AliAltroRawStreamV3::NextDDL()) return kFALSE;
  
  Int_t ddlNumber = GetDDLNumber();
  if (ddlNumber < 72) {
    fSector = ddlNumber / 2;
    fPatchIndex = ddlNumber % 2;
  }
  else {
    fSector = (ddlNumber - 72) / 4 + 36;
    fPatchIndex = (ddlNumber - 72) % 4 + 2;
  }
  return kTRUE;
}
//_____________________________________________________________________________
void AliTPCRawStreamV3::ApplyAltroMapping()
{
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
  

//  if ((fRow < 0) || (fPad < 0))
//    AddMappingErrorLog(Form("hw=%d",hwAddress));
}
