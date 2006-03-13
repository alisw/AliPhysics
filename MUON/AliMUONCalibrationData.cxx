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
#include "AliLog.h"
#include "AliMUONV2DStore.h"
#include "AliMUONVCalibParam.h"
#include "Riostream.h"

ClassImp(AliMUONCalibrationData)

//_____________________________________________________________________________
AliMUONCalibrationData::AliMUONCalibrationData(Int_t runNumber, 
                                               Bool_t deferredInitialization) 
: TObject(), 
fIsValid(kTRUE),
fRunNumber(runNumber), 
fGains(0x0), 
fPedestals(0x0),
fDeadChannels(0x0)
{
  //
  // Default ctor.
  // If deferredInitialization is false, we read *all* calibrations
  // at once.
  // So when using this class to access only one kind of calibrations (e.g.
  // only pedestals), you should put deferredInitialization to kTRUE, which
  // will instruct this object to fetch the data only when neeeded.
  //
  if ( deferredInitialization == kFALSE )
  {
    Gains();
    Pedestals();
    DeadChannels();
  }
}

//______________________________________________________________________________
AliMUONCalibrationData::AliMUONCalibrationData(const AliMUONCalibrationData& right) 
  : TObject(right) 
{  
/// Protected copy constructor (not implemented)

  AliFatal("Copy constructor not provided.");
}

//_____________________________________________________________________________
AliMUONCalibrationData::~AliMUONCalibrationData()
{
  //
  // dtor. Note that we're the owner of our pointers.
  //
  delete fPedestals;
  delete fGains;
  delete fDeadChannels;
}

//______________________________________________________________________________
AliMUONCalibrationData& 
AliMUONCalibrationData::operator=(const AliMUONCalibrationData& right)
{
/// Protected assignement operator (not implemented)

  // check assignement to self
  if (this == &right) return *this;

  AliFatal("Assignement operator not provided.");
    
  return *this;  
}    

//_____________________________________________________________________________
AliMUONVCalibParam*
AliMUONCalibrationData::DeadChannel(Int_t detElemId, Int_t manuId) const
{
  //
  // Return the calibration for a given (detElemId, manuId) pair
  // Note that for DeadChannel, it's "legal" to return 0x0 (e.g. if a manu
  // is perfect, we might simply forget it in the store).
  //
  return
  static_cast<AliMUONVCalibParam*>(DeadChannels()->Get(detElemId,manuId));
}

//_____________________________________________________________________________
AliMUONV2DStore*
AliMUONCalibrationData::DeadChannels() const
{
  //
  // Create (if needed) and return the internal store for DeadChannels.
  //
  if (!fDeadChannels)
  {
    AliCDBEntry* entry = GetEntry("MUON/Calib/DeadChannels");
    if (entry)
    {
      fDeadChannels = dynamic_cast<AliMUONV2DStore*>(entry->GetObject());
      if (!fDeadChannels)
      {
        AliError("fDeadChannels not of the expected type !!!");
      }
    }
    else
    {
      AliError("Could not get dead channels !");
    }
  }
  return fDeadChannels;
}

//_____________________________________________________________________________
AliCDBEntry*
AliMUONCalibrationData::GetEntry(const char* path) const
{
  //
  // Access the CDB for a given path (e.g. MUON/Calib/Pedestals),
  // and return the corresponding CDBEntry.
  //
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
AliMUONVCalibParam*
AliMUONCalibrationData::Gain(Int_t detElemId, Int_t manuId) const
{
  //
  // Return the gains for a given (detElemId, manuId) pair
  // Note that, unlike the DeadChannel case, if the result is 0x0, that's an
  // error (meaning that we should get gains for all channels).
  //
  AliMUONVCalibParam* gain = 
    static_cast<AliMUONVCalibParam*>(Gains()->Get(detElemId,manuId));
  if (!gain)
  {
    AliError(Form("Could not get gain for detElemId=%d manuId=%d ",
                    detElemId,manuId));
  }
  return gain;
}

//_____________________________________________________________________________
AliMUONV2DStore*
AliMUONCalibrationData::Gains() const
{
  //
  // Create (if needed) and return the internal store for gains.
  //
  if (!fGains)
  {
    AliCDBEntry* entry = GetEntry("MUON/Calib/Gains");
    if (entry)
    {
      fGains = dynamic_cast<AliMUONV2DStore*>(entry->GetObject());
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
AliMUONV2DStore*
AliMUONCalibrationData::Pedestals() const
{
  //
  // Create (if needed) and return the internal storage for pedestals.
  //
  if (!fPedestals)
  {
    AliCDBEntry* entry = GetEntry("MUON/Calib/Pedestals");
    if (entry)
    {
      fPedestals = dynamic_cast<AliMUONV2DStore*>(entry->GetObject());
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
  //
  // A very basic dump of our guts.
  //  
  cout << "RunNumber " << RunNumber()
    << " fGains=" << fGains
    << " fPedestals=" << fPedestals
    << " fDeadChannels=" << fDeadChannels
  << endl;
}


//_____________________________________________________________________________
AliMUONVCalibParam*
AliMUONCalibrationData::Pedestal(Int_t detElemId, Int_t manuId) const
{
  //
  // Return the pedestals for a given (detElemId, manuId) pair.
  // A return value of 0x0 is considered an error, meaning we should get
  // pedestals for all channels.
  //
  AliMUONVCalibParam* ped = 
    static_cast<AliMUONVCalibParam*>(Pedestals()->Get(detElemId,manuId));
  if (!ped)
  {
    AliError(Form("Could not get pedestal for detElemId=%d manuId=%d ",
                  detElemId,manuId));
  }
  return ped;
}


