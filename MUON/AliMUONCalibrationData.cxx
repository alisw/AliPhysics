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

// $Id$

#include "AliMUONCalibrationData.h"

#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliMUONCalibParam.h"
#include "AliLog.h"
#include "AliMUONV3DStore.h"
#include "Riostream.h"

ClassImp(AliMUONCalibrationData)

//_____________________________________________________________________________
AliMUONCalibrationData::AliMUONCalibrationData(Int_t runNumber, 
                                               Bool_t deferredInitialization) 
: TObject(), 
fIsValid(kTRUE),
fRunNumber(runNumber), 
fGains(0x0), 
fPedestals(0x0)
{
  if ( deferredInitialization == kFALSE )
  {
    Gains();
    Pedestals();
  }
}


//_____________________________________________________________________________
AliMUONCalibrationData::~AliMUONCalibrationData()
{
  delete fPedestals;
  delete fGains;
}

//_____________________________________________________________________________
AliCDBEntry*
AliMUONCalibrationData::GetEntry(const char* path) const
{
  AliInfo(Form("Fetching %s from Condition DataBase for run %d",path,fRunNumber));
  
  AliCDBManager* man = AliCDBManager::Instance();
  if (!man->IsDefaultStorageSet())
  {
    AliError("No default CDB storage set !");
    fIsValid = kFALSE;
    return 0;
  }
  
  AliCDBStorage* storage = man->GetDefaultStorage();
  
  AliCDBEntry* entry = storage->Get(path,fRunNumber);
  return entry;
}

//_____________________________________________________________________________
AliMUONCalibParam*
AliMUONCalibrationData::Gain(Int_t detElemId, 
                             Int_t manuId, Int_t manuChannel) const
{
  AliMUONCalibParam* gain = 
  static_cast<AliMUONCalibParam*>(Gains()->Get(detElemId,manuId,manuChannel));
  if (!gain)
  {
    AliError(Form("Could not get gain for detElemId=%d manuId=%d "
                  "manuChannel=%d",detElemId,manuId,manuChannel));
  }
  return gain;
}

//_____________________________________________________________________________
AliMUONV3DStore*
AliMUONCalibrationData::Gains() const
{
  if (!fGains)
  {
    AliCDBEntry* entry = GetEntry("MUON/Calib/Gains");
    if (entry)
    {
      fGains = dynamic_cast<AliMUONV3DStore*>(entry->GetObject());
      if (!fGains)
      {
        AliError("Gains not of the expected type !!!");
      }
    }
    else
    {
      AliError("Could not get gains !");
    }
  }
  return fGains;
}

//_____________________________________________________________________________
Bool_t
AliMUONCalibrationData::IsValid() const
{
  return fIsValid;
}

//_____________________________________________________________________________
AliMUONV3DStore*
AliMUONCalibrationData::Pedestals() const
{
  if (!fPedestals)
  {
    AliCDBEntry* entry = GetEntry("MUON/Calib/Pedestals");
    if (entry)
    {
      fPedestals = dynamic_cast<AliMUONV3DStore*>(entry->GetObject());
      if (!fPedestals)
      {
        AliError("fPedestals not of the expected type !!!");
      }
    }
    else
    {
      AliError("Could not get pedestals !");
    }
  }
  return fPedestals;
}

//_____________________________________________________________________________
void
AliMUONCalibrationData::Print(Option_t*) const
{
  cout << "RunNumber " << RunNumber()
    << " fGains=" << fGains
  << " fPedestals=" << fPedestals
  << endl;
}

//_____________________________________________________________________________
Int_t
AliMUONCalibrationData::RunNumber() const
{
  return fRunNumber;
}

//_____________________________________________________________________________
AliMUONCalibParam*
AliMUONCalibrationData::Pedestal(Int_t detElemId, 
                                 Int_t manuId, Int_t manuChannel) const
{
  AliMUONCalibParam* ped = 
    static_cast<AliMUONCalibParam*>(Pedestals()->Get(detElemId,manuId,manuChannel));
  if (!ped)
  {
    AliError(Form("Could not get pedestal for detElemId=%d manuId=%d "
                  "manuChannel=%d",detElemId,manuId,manuChannel));
  }
  return ped;
}


