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

// This class handles the mapping of the Altro channels in the TPC
// The mapping is read from an external mapping files
// Author: C.Cheshkov

#include "AliTPCAltroMapping.h"
#include "AliLog.h"
#include <Riostream.h>
//#include <stdlib.h>


ClassImp(AliTPCAltroMapping)

//_____________________________________________________________________________
AliTPCAltroMapping::AliTPCAltroMapping(const char *mappingFile):
  AliAltroMapping(mappingFile),
  fMinPadRow(0),
  fMaxPadRow(0),
  fMaxPad(0),
  fMapping(NULL),
  fInvMapping(NULL)
{
  // Constructor
  ReadMapping();
  CloseMappingFile();
}

//_____________________________________________________________________________
AliTPCAltroMapping::~AliTPCAltroMapping()
{
  // destructor
  DeleteMappingArrays();
}

//_____________________________________________________________________________
AliTPCAltroMapping::AliTPCAltroMapping(const AliTPCAltroMapping& mapping):
  AliAltroMapping(mapping),
  fMinPadRow(mapping.fMinPadRow),
  fMaxPadRow(mapping.fMaxPadRow),
  fMaxPad(mapping.fMaxPad),
  fMapping(mapping.fMapping),
  fInvMapping(mapping.fInvMapping)
{
// Copy Constructor

  Fatal("AliTPCAltroMapping", "copy constructor not implemented");
}

//_____________________________________________________________________________
AliTPCAltroMapping& AliTPCAltroMapping::operator = (const AliTPCAltroMapping& /*mapping*/)
{
//Assigment operator

  Fatal("operator =", "assignment operator not implemented");
  return *this;
}

//_____________________________________________________________________________
Bool_t AliTPCAltroMapping::ReadMapping()
{
  // Initalizes the ALTRO mapping from a file
  // Look at the TPC module for the format of
  // the mapping file
  if (!fIn) {
    AliFatal("Mapping file has not been opened !");
    return kFALSE;
  }

  fMinPadRow = 0x7fffffff;
  fMaxPadRow = 0;
  fMaxPad = 0;
  fMapping = new Short_t*[fMaxHWAddress+1];
  for (Int_t i = 0; i <= fMaxHWAddress; i++) {
    fMapping[i] = new Short_t[2];
    fMapping[i][0] = fMapping[i][1] = -1;
  }
 
  for(Int_t i = 0; i < fNumberOfChannels ; i++) { //5504 is size of irorc mapping at ther moment only for irorc
    Int_t hwAddress;
    if (!(*fIn >> hwAddress)) {
      AliFatal("Syntax of the mapping file is wrong !");
      return kFALSE;
    }
    if (hwAddress > fMaxHWAddress) {
      AliFatal(Form("Hardware (ALTRO) adress (%d) outside the range (0 -> %d) !",hwAddress,fMaxHWAddress));
      return kFALSE;
    }
    Int_t padrow,pad;
    if (!(*fIn >> padrow >> pad)) {
      AliFatal("Syntax of the mapping file is wrong !");
      return kFALSE;
    }
 
    fMapping[hwAddress][0] = padrow;
    fMapping[hwAddress][1] = pad;

    if (padrow > fMaxPadRow) fMaxPadRow = padrow;
    if (padrow < fMinPadRow) fMinPadRow = padrow;
    if (pad > fMaxPad) fMaxPad = pad;
  }

  fInvMapping = new Short_t*[fMaxPadRow - fMinPadRow + 1];
  for (Int_t i = 0; i <= (fMaxPadRow - fMinPadRow); i++) {
    fInvMapping[i] = new Short_t[fMaxPad + 1];
    for (Int_t j = 0; j <= fMaxPad; j++) fInvMapping[i][j] = -1;
  }

  for(Int_t i = 0; i <= fMaxHWAddress; i++) {
    Int_t padrow = fMapping[i][0];
    Int_t pad = fMapping[i][1];
    if(padrow != -1 && pad != -1)
      fInvMapping[padrow-fMinPadRow][pad] = i;
  }

  return kTRUE;
}

//_____________________________________________________________________________
Int_t AliTPCAltroMapping::GetHWAddress(Int_t padrow, Int_t pad, Int_t /* sector */) const
{
  // Get the content of the mapping array
  // return -1 in case there is no hardware
  // adress defined for these pad-row and pad
  if (!fInvMapping) {
    AliWarning("Mapping array was not initalized correctly !");
    return -1;
  }
  if (padrow < fMinPadRow || padrow > fMaxPadRow) {
    AliWarning(Form("Index of pad-row (%d) outside the range (%d -> %d) !",padrow,fMinPadRow,fMaxPadRow));
    return -1;
  }
  if (pad > fMaxPad) {
    AliWarning(Form("Index of pad (%d) outside the range (0 -> %d) !",pad,fMaxPad));
    return -1;
  }
  Int_t hwAddress = fInvMapping[padrow-fMinPadRow][pad];
  if (hwAddress == -1)
    AliWarning(Form("Hardware (ALTRO) adress is not defined for these pad-row (%d) and pad (%d) !",padrow,pad));

  return hwAddress;
}

//_____________________________________________________________________________
Int_t AliTPCAltroMapping::GetPadRow(Int_t hwAddress) const
{
  if (!fMapping) {
    AliWarning("Mapping array was not initalized correctly !");
    return -1;
  }
  if (hwAddress > fMaxHWAddress) {
    AliWarning(Form("Hardware (ALTRO) adress (%d) outside the range (0 -> %d) !",hwAddress,fMaxHWAddress));
    return -1;
  }
  Int_t padrow = fMapping[hwAddress][0];
  if (padrow == -1)
    AliWarning(Form("Hardware (ALTRO) adress (%d) is not defined !",hwAddress));

  return padrow;
}

//_____________________________________________________________________________
Int_t AliTPCAltroMapping::GetPad(Int_t hwAddress) const
{
  if (!fMapping) {
    AliWarning("Mapping array was not initalized correctly !");
    return -1;
  }
  if (hwAddress > fMaxHWAddress) {
    AliWarning(Form("Hardware (ALTRO) adress (%d) outside the range (0 -> %d) !",hwAddress,fMaxHWAddress));
    return -1;
  }
  Int_t pad = fMapping[hwAddress][1];
  if (pad == -1)
    AliWarning(Form("Hardware (ALTRO) adress (%d) is not defined !",hwAddress));

  return pad;
}

//_____________________________________________________________________________
Int_t AliTPCAltroMapping::GetSector(Int_t /* hwAddress */) const
{
  AliWarning("Sector index is not contained in the TPC altro mapping !");
  return -1;
}

//_____________________________________________________________________________
void AliTPCAltroMapping::DeleteMappingArrays()
{
  // Deletes the arrays which have been
  // allocated during the reading of the
  // mapping file
  if (fMapping) {
    for (Int_t i = 0; i <= fMaxHWAddress; i++) delete [] fMapping[i];
    delete [] fMapping;
  }

  if (fInvMapping) {
    for (Int_t i = 0; i <= (fMaxPadRow - fMinPadRow); i++)
      delete [] fInvMapping[i];
    delete [] fInvMapping;
  }
}
