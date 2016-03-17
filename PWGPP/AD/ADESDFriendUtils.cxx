// -*- C++ -*-
// $Id$

/**************************************************************************
 * Author: C. Mayer                                                       *
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

#include "AliLog.h"
#include "AliESDADfriend.h"
#include "AliADCalibData.h"
#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliCDBStorage.h"
#include "ADESDFriendUtils.h"
#include "AliESDADfriend.h"

ClassImp(ADESDFriendUtils);

ADESDFriendUtils::ADESDFriendUtils()
  : TObject()
  , fCalibData(NULL) {
  for (Int_t ch=0; ch<AliADRawStream::kNChannels; ++ch)
    for (Int_t bc=0; bc<AliADRawStream::kNEvOfInt; ++bc)
      fADCPedSub[ch][bc] = 0.0f;
}

ADESDFriendUtils::~ADESDFriendUtils() {
}

AliCDBManager* ADESDFriendUtils::Init(Int_t runNumber) {
  // (1) set up the OCDB
  AliCDBManager *man = AliCDBManager::Instance();
  if (NULL == man) {
    AliFatal("CDB manager not found");
    return NULL;
  }

  if (!man->IsDefaultStorageSet())
    man->SetDefaultStorage("raw://");

  man->SetRun(runNumber);

  // (2) Get the AD calibration data OCDB object (pedestals)
  AliCDBEntry *entry = man->Get("AD/Calib/Data");
  if (NULL == entry) {
    AliFatal("AD/Calib/Data not found");
    return NULL;
  }
  fCalibData = dynamic_cast<AliADCalibData*>(entry->GetObject());
  if (NULL == fCalibData) {
    AliFatal("No calibration data from calibration database");
    return NULL;
  }
  return man;
}

void ADESDFriendUtils::Update(const AliESDADfriend* esdADfriend) {
  if (NULL == esdADfriend) {
    AliError("NULL == esdADfriend");
    return;
  }  
  for (Int_t ch=0; ch<AliADRawStream::kNChannels; ++ch) {    
    for (Int_t bc=0; bc<AliADRawStream::kNEvOfInt; ++bc) {
      fADCPedSub[ch][bc]  = Float_t(esdADfriend->GetPedestal(ch, bc));
      fADCPedSub[ch][bc] -= fCalibData->GetPedestal(ch + 16*esdADfriend->GetIntegratorFlag(ch, bc));
    }
  }
}

Bool_t ADESDFriendUtils::IsPileUp(Int_t ch, Float_t threshold) const {  
  Bool_t isPileUp = kFALSE;
  for (Int_t bc=13; bc<AliADRawStream::kNEvOfInt-1 && !isPileUp; ++bc)  
    isPileUp |= (fADCPedSub[ch][bc+1] > fADCPedSub[ch][bc] + threshold);
  
  return isPileUp;
}
