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
AliTPCAltroMapping::AliTPCAltroMapping():
  AliAltroMapping(),
  fMinPadRow(0),
  fMaxPadRow(0),
  fMaxPad(0),
  fInvMapping(NULL)
{
  // Default constructor
}

//_____________________________________________________________________________
AliTPCAltroMapping::AliTPCAltroMapping(const char *mappingFile):
  AliAltroMapping(mappingFile),
  fMinPadRow(0),
  fMaxPadRow(0),
  fMaxPad(0),
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
  if (fInvMapping) delete [] fInvMapping;
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
  fMappingSize = 2*(fMaxHWAddress+1);
  fMapping = new Short_t[fMappingSize];
  for (Int_t i = 0; i <= fMaxHWAddress; i++) {
    fMapping[2*i] = fMapping[2*i+1] = -1;
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
 
    fMapping[2*hwAddress] = padrow;
    fMapping[2*hwAddress+1] = pad;

    if (padrow > fMaxPadRow) fMaxPadRow = padrow;
    if (padrow < fMinPadRow) fMinPadRow = padrow;
    if (pad > fMaxPad) fMaxPad = pad;
  }

  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliTPCAltroMapping::CreateInvMapping()
{
  // Create the inverse mapping
  // needed for the simulation of
  // raw data
  if (fInvMapping) return kTRUE;

  if (!fMapping) {
    AliWarning("Mapping array was not initalized correctly ! Impossible to create the inverse mapping !");
    return kFALSE;
  }

  Int_t nRows = fMaxPadRow - fMinPadRow + 1;
  Int_t nPads = fMaxPad + 1;
  Int_t invMappingSize = nRows*nPads;

  fInvMapping = new Short_t[invMappingSize];
  for (Int_t i = 0; i <= (fMaxPadRow - fMinPadRow); i++) {
    for (Int_t j = 0; j <= fMaxPad; j++) fInvMapping[nPads*i+j] = -1;
  }

  for(Int_t i = 0; i <= fMaxHWAddress; i++) {
    Int_t padrow = fMapping[2*i];
    Int_t pad = fMapping[2*i+1];
    if(padrow != -1 && pad != -1)
      fInvMapping[nPads*(padrow-fMinPadRow)+pad] = i;
  }

  return kTRUE;
}

//_____________________________________________________________________________
Int_t AliTPCAltroMapping::GetHWAddress(Int_t padrow, Int_t pad, Int_t /* sector */)
{
  // Get the content of the mapping array
  // return -1 in case there is no hardware
  // adress defined for these pad-row and pad
  if (!fInvMapping) {
    if (!CreateInvMapping()) return -1;
  }
  if (padrow < fMinPadRow || padrow > fMaxPadRow) {
    AliWarning(Form("Index of pad-row (%d) outside the range (%d -> %d) !",padrow,fMinPadRow,fMaxPadRow));
    return -1;
  }
  if (pad > fMaxPad) {
    AliWarning(Form("Index of pad (%d) outside the range (0 -> %d) !",pad,fMaxPad));
    return -1;
  }
  Int_t hwAddress = fInvMapping[(fMaxPad+1)*(padrow-fMinPadRow)+pad];
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
  Int_t padrow = fMapping[2*hwAddress];
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
  Int_t pad = fMapping[2*hwAddress+1];
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
