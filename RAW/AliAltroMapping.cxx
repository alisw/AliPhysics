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

// This class handles the mapping of the Altro channels
// The mapping is read from an external mapping files
// The class is used by TPC,PHOS and FMD
// Author: C.Cheshkov

#include "AliAltroMapping.h"
#include "AliLog.h"
#include <Riostream.h>
//#include <stdlib.h>


ClassImp(AliAltroMapping)

//_____________________________________________________________________________
AliAltroMapping::AliAltroMapping(const char *mappingFile):
  fNumberOfChannels(0),
  fMaxHWAdress(0),
  fMinPadRow(0),
  fMaxPadRow(0),
  fMaxPad(0),
  fMapping(NULL),
  fInvMapping(NULL)
{
  // Constructor
  // Reads the mapping from an external file
  if (mappingFile)
    ReadMapping(mappingFile);
  else
    AliFatal("Mapping file not specified !");
}

//_____________________________________________________________________________
AliAltroMapping::~AliAltroMapping()
{
  // destructor
  if (fMapping) {
    for (Int_t i = 0; i <= fMaxHWAdress; i++) delete [] fMapping[i];
    delete [] fMapping;
  }

  if (fInvMapping) {
    for (Int_t i = 0; i <= (fMaxPadRow - fMinPadRow); i++)
      delete [] fInvMapping[i];
    delete [] fInvMapping;
  }

}

//_____________________________________________________________________________
AliAltroMapping::AliAltroMapping(const AliAltroMapping& mapping):
  TObject(mapping),
  fNumberOfChannels(mapping.fNumberOfChannels),
  fMaxHWAdress(mapping.fMaxHWAdress),
  fMinPadRow(mapping.fMinPadRow),
  fMaxPadRow(mapping.fMaxPadRow),
  fMaxPad(mapping.fMaxPad),
  fMapping(mapping.fMapping),
  fInvMapping(mapping.fInvMapping)
{
// Copy Constructor

  Fatal("AliAltroMapping", "copy constructor not implemented");
}

//_____________________________________________________________________________
AliAltroMapping& AliAltroMapping::operator = (const AliAltroMapping& /*mapping*/)
{
//Assigment operator

  Fatal("operator =", "assignment operator not implemented");
  return *this;
}

//_____________________________________________________________________________
Bool_t AliAltroMapping::ReadMapping(const char *mappingFile)
{
  // Initalizes the ALTRO mapping from a file
  // Look at the TPC module for the format of
  // the mapping file
  ifstream in(mappingFile);
  if (!in) {
    AliFatal(Form("Missing mapping file (%s) !",mappingFile));
    return kFALSE;
  }
  if (!(in >> fNumberOfChannels)) {
    AliFatal(Form("Syntax of the mapping file is wrong (%s) !",mappingFile));
    return kFALSE;
  }
  if (!(in >> fMaxHWAdress)) {
    AliFatal(Form("Syntax of the mapping file is wrong (%s) !",mappingFile));
    return kFALSE;
  }
//   if (!(in >> fMinPadRow >> fMaxPadRow)) {
//     AliFatal(Form("Syntax of the mapping file is wrong (%s) !",mappingFile));
//     return kFALSE;
//   }
//   if (!(in >> fMaxPad)) {
//     AliFatal(Form("Syntax of the mapping file is wrong (%s) !",mappingFile));
//     return kFALSE;
//   }

 
  fMinPadRow = 0x7fffffff;
  fMaxPadRow = 0;
  fMaxPad = 0;
  fMapping = new Short_t*[fMaxHWAdress+1];
  for (Int_t i = 0; i <= fMaxHWAdress; i++) {
    fMapping[i] = new Short_t[2];
    fMapping[i][0] = fMapping[i][1] = -1;
  }
 
  for(Int_t i = 0; i < fNumberOfChannels ; i++) { //5504 is size of irorc mapping at ther moment only for irorc
    Int_t hwAdress;
    if (!(in >> hwAdress)) {
      AliFatal(Form("Syntax of the mapping file is wrong (%s) !",mappingFile));
      return kFALSE;
    }
    if (hwAdress > fMaxHWAdress) {
      AliFatal(Form("Hardware (ALTRO) adress (%d) outside the range (0 -> %d) !",hwAdress,fMaxHWAdress));
      return kFALSE;
    }
    Int_t padrow,pad;
    if (!(in >> padrow >> pad)) {
      AliFatal(Form("Syntax of the mapping file is wrong (%s) !",mappingFile));
      return kFALSE;
    }
 
    fMapping[hwAdress][0] = padrow;
    fMapping[hwAdress][1] = pad;

    if (padrow > fMaxPadRow) fMaxPadRow = padrow;
    if (padrow < fMinPadRow) fMinPadRow = padrow;
    if (pad > fMaxPad) fMaxPad = pad;
  }

  fInvMapping = new Short_t*[fMaxPadRow - fMinPadRow + 1];
  for (Int_t i = 0; i <= (fMaxPadRow - fMinPadRow); i++) {
    fInvMapping[i] = new Short_t[fMaxPad + 1];
    for (Int_t j = 0; j <= fMaxPad; j++) fInvMapping[i][j] = -1;
  }

  for(Int_t i = 0; i <= fMaxHWAdress; i++) {
    Int_t padrow = fMapping[i][0];
    Int_t pad = fMapping[i][1];
    if(padrow != -1 && pad != -1)
      fInvMapping[padrow-fMinPadRow][pad] = i;
  }

  in.close();
  return kTRUE;
}

//_____________________________________________________________________________
Int_t AliAltroMapping::GetHWAdress(Int_t padrow, Int_t pad) const
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
  Int_t hwAdress = fInvMapping[padrow-fMinPadRow][pad];
  if (hwAdress == -1)
    AliWarning(Form("Hardware (ALTRO) adress is not defined for these pad-row (%d) and pad (%d) !",padrow,pad));

  return hwAdress;
}

//_____________________________________________________________________________
Int_t AliAltroMapping::GetPadRow(Int_t hwAdress) const
{
  if (!fMapping) {
    AliWarning("Mapping array was not initalized correctly !");
    return -1;
  }
  if (hwAdress > fMaxHWAdress) {
    AliWarning(Form("Hardware (ALTRO) adress (%d) outside the range (0 -> %d) !",hwAdress,fMaxHWAdress));
    return -1;
  }
  Int_t padrow = fMapping[hwAdress][0];
  if (padrow == -1)
    AliWarning(Form("Hardware (ALTRO) adress (%d) is not defined !",hwAdress));

  return padrow;
}

//_____________________________________________________________________________
Int_t AliAltroMapping::GetPad(Int_t hwAdress) const
{
  if (!fMapping) {
    AliWarning("Mapping array was not initalized correctly !");
    return -1;
  }
  if (hwAdress > fMaxHWAdress) {
    AliWarning(Form("Hardware (ALTRO) adress (%d) outside the range (0 -> %d) !",hwAdress,fMaxHWAdress));
    return -1;
  }
  Int_t pad = fMapping[hwAdress][1];
  if (pad == -1)
    AliWarning(Form("Hardware (ALTRO) adress (%d) is not defined !",hwAdress));

  return pad;
}
