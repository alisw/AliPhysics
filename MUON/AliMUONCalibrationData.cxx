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
#include "AliLog.h"
#include "AliMUONTriggerEfficiencyCells.h"
#include "AliMUONTriggerLut.h"
#include "AliMUONV1DStore.h"
#include "AliMUONV2DStore.h"
#include "AliMUONVCalibParam.h"
#include "Riostream.h"

///
/// This class is to be considered as a convenience class.
/// Its aim is to ease retrieval of calibration data from the 
/// condition database.
///
/// It acts as a "facade" to a bunch of underlying 
/// containers/calibration classes.
///

ClassImp(AliMUONCalibrationData)

//_____________________________________________________________________________
AliMUONCalibrationData::AliMUONCalibrationData(Int_t runNumber, 
                                               Bool_t deferredInitialization) 
: TObject(), 
fIsValid(kTRUE),
fRunNumber(runNumber), 
fGains(0x0), 
fPedestals(0x0),
fDeadChannels(0x0),
fLocalTriggerBoardMasks(0x0),
fRegionalTriggerBoardMasks(0x0),
fGlobalTriggerBoardMasks(0x0),
fTriggerLut(0x0),
fTriggerEfficiency(0x0)
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
    OnDemandGains();
    OnDemandPedestals();
    OnDemandDeadChannels();
    OnDemandLocalTriggerBoardMasks();
    OnDemandRegionalTriggerBoardMasks();
    OnDemandGlobalTriggerBoardMasks();
    OnDemandTriggerLut();
    OnDemandTriggerEfficiency();
  }
}

//_____________________________________________________________________________
AliMUONCalibrationData::AliMUONCalibrationData(const AliMUONCalibrationData&)
: TObject()
{
  AliFatal("Implement me if needed");
}

//_____________________________________________________________________________
AliMUONCalibrationData&
AliMUONCalibrationData::operator=(const AliMUONCalibrationData&)
{
  AliFatal("Implement me if needed");
  return *this;
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
  delete fLocalTriggerBoardMasks;
  delete fRegionalTriggerBoardMasks;
  delete fGlobalTriggerBoardMasks;
  delete fTriggerLut;
  delete fTriggerEfficiency;
}

//_____________________________________________________________________________
AliMUONVCalibParam*
AliMUONCalibrationData::DeadChannels(Int_t detElemId, Int_t manuId) const
{
  //
  // Return the calibration for a given (detElemId, manuId) pair
  // Note that for DeadChannel, it's "legal" to return 0x0 (e.g. if a manu
  // is perfect, we might simply forget it in the store).
  //
  return
  static_cast<AliMUONVCalibParam*>(OnDemandDeadChannels()->Get(detElemId,manuId));
}

//_____________________________________________________________________________
AliMUONV2DStore*
AliMUONCalibrationData::OnDemandDeadChannels() const
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
  return AliCDBManager::Instance()->Get(path,fRunNumber);
}

//_____________________________________________________________________________
AliMUONVCalibParam*
AliMUONCalibrationData::Gains(Int_t detElemId, Int_t manuId) const
{
  //
  // Return the gains for a given (detElemId, manuId) pair
  // Note that, unlike the DeadChannel case, if the result is 0x0, that's an
  // error (meaning that we should get gains for all channels).
  //
  AliMUONVCalibParam* gain = 
    static_cast<AliMUONVCalibParam*>(OnDemandGains()->Get(detElemId,manuId));
  if (!gain)
  {
    AliError(Form("Could not get gain for detElemId=%d manuId=%d ",
                    detElemId,manuId));
  }
  return gain;
}

//_____________________________________________________________________________
AliMUONV2DStore*
AliMUONCalibrationData::OnDemandGains() const
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
AliMUONVCalibParam* 
AliMUONCalibrationData::GlobalTriggerBoardMasks() const
{
  //
  // Return the masks for the global trigger board.
  //
  return OnDemandGlobalTriggerBoardMasks();
}

//_____________________________________________________________________________
AliMUONVCalibParam*
AliMUONCalibrationData::OnDemandGlobalTriggerBoardMasks() const
{
  //
  // Create (if needed) and return the internal store for GlobalTriggerBoardMasks.
  //
  if (!fGlobalTriggerBoardMasks)
  {
    AliCDBEntry* entry = GetEntry("MUON/Calib/GlobalTriggerBoardMasks");
    if (entry)
    {
      fGlobalTriggerBoardMasks = dynamic_cast<AliMUONVCalibParam*>(entry->GetObject());
      if (!fGlobalTriggerBoardMasks)
      {
        AliError("fGlobalTriggerBoardMasks not of the expected type !!!");
      }
    }
    else
    {
      AliError("Could not get global trigger board masks !");
    }
  }
  return fGlobalTriggerBoardMasks;
}

//_____________________________________________________________________________
AliMUONVCalibParam* 
AliMUONCalibrationData::LocalTriggerBoardMasks(Int_t localBoardNumber) const
{
  //
  // Return the masks for a given trigger local board.
  //
  AliMUONVCalibParam* ltbm = 
  static_cast<AliMUONVCalibParam*>(OnDemandLocalTriggerBoardMasks()->Get(localBoardNumber));
  if (!ltbm)
  {
    AliError(Form("Could not get mask for localBoardNumber=%d",localBoardNumber));
  }
  return ltbm;  
}

//_____________________________________________________________________________
AliMUONV1DStore*
AliMUONCalibrationData::OnDemandLocalTriggerBoardMasks() const
{
  //
  // Create (if needed) and return the internal store for LocalTriggerBoardMasks.
  //
  if (!fLocalTriggerBoardMasks)
  {
    AliCDBEntry* entry = GetEntry("MUON/Calib/LocalTriggerBoardMasks");
    if (entry)
    {
      fLocalTriggerBoardMasks = dynamic_cast<AliMUONV1DStore*>(entry->GetObject());
      if (!fLocalTriggerBoardMasks)
      {
        AliError("fLocalTriggerBoardMasks not of the expected type !!!");
      }
    }
    else
    {
      AliError("Could not get local trigger board masks !");
    }
  }
  return fLocalTriggerBoardMasks;
}

//_____________________________________________________________________________
AliMUONV2DStore*
AliMUONCalibrationData::OnDemandPedestals() const
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
  << " fLocalTriggerBoardMasks=" << fLocalTriggerBoardMasks
  << " fRegionalTriggerBoardMasks=" << fRegionalTriggerBoardMasks
  << " fGlobalTriggerBoardMasks=" << fGlobalTriggerBoardMasks
  << " fTriggerLut=" << fTriggerLut
  << endl;
}


//_____________________________________________________________________________
AliMUONVCalibParam*
AliMUONCalibrationData::Pedestals(Int_t detElemId, Int_t manuId) const
{
  //
  // Return the pedestals for a given (detElemId, manuId) pair.
  // A return value of 0x0 is considered an error, meaning we should get
  // pedestals for all channels.
  //
  AliMUONVCalibParam* ped = 
    static_cast<AliMUONVCalibParam*>(OnDemandPedestals()->Get(detElemId,manuId));
  if (!ped)
  {
    AliError(Form("Could not get pedestal for detElemId=%d manuId=%d ",
                  detElemId,manuId));
  }
  return ped;
}

//_____________________________________________________________________________
AliMUONVCalibParam* 
AliMUONCalibrationData::RegionalTriggerBoardMasks(Int_t index) const
{
  //
  // Return the masks for a given trigger regional board.
  //
  AliMUONVCalibParam* rtbm = 
  static_cast<AliMUONVCalibParam*>(OnDemandRegionalTriggerBoardMasks()->Get(index));
  if (!rtbm)
  {
    AliError(Form("Could not get mask for regionalBoard index=%d",index));
  }
  return rtbm;  
}

//_____________________________________________________________________________
AliMUONV1DStore*
AliMUONCalibrationData::OnDemandRegionalTriggerBoardMasks() const
{
  //
  // Create (if needed) and return the internal store for RegionalTriggerBoardMasks.
  //
  if (!fRegionalTriggerBoardMasks)
  {
    AliCDBEntry* entry = GetEntry("MUON/Calib/RegionalTriggerBoardMasks");
    if (entry)
    {
      fRegionalTriggerBoardMasks = dynamic_cast<AliMUONV1DStore*>(entry->GetObject());
      if (!fRegionalTriggerBoardMasks)
      {
        AliError("fRegionalTriggerBoardMasks not of the expected type !!!");
      }
    }
    else
    {
      AliError("Could not get regional trigger board masks !");
    }
  }
  return fRegionalTriggerBoardMasks;
}

//_____________________________________________________________________________
AliMUONTriggerEfficiencyCells*
AliMUONCalibrationData::TriggerEfficiency() const
{
  //
  // Return the trigger efficiency.
  //
  return OnDemandTriggerEfficiency();
}

//_____________________________________________________________________________
AliMUONTriggerEfficiencyCells* 
AliMUONCalibrationData::OnDemandTriggerEfficiency() const
{
  //
  //
  //
  if (!fTriggerEfficiency)
  {
    AliCDBEntry* entry = GetEntry("MUON/Calib/TriggerEfficiency");
    if (entry)
    {
      fTriggerEfficiency = dynamic_cast<AliMUONTriggerEfficiencyCells*>(entry->GetObject());
      if (!fTriggerEfficiency)
      {
        AliError("fTriggerEfficiency not of the expected type !!!");
      }
    }
    else
    {
      AliError("Could not get trigger efficiency !");
    }
  }
  return fTriggerEfficiency;
}

//_____________________________________________________________________________
AliMUONTriggerLut*
AliMUONCalibrationData::TriggerLut() const
{
  //
  // Return the trigger look up table.
  //
  return OnDemandTriggerLut();
}

//_____________________________________________________________________________
AliMUONTriggerLut* 
AliMUONCalibrationData::OnDemandTriggerLut() const
{
  //
  //
  //
  if (!fTriggerLut)
  {
    AliCDBEntry* entry = GetEntry("MUON/Calib/TriggerLut");
    if (entry)
    {
      fTriggerLut = dynamic_cast<AliMUONTriggerLut*>(entry->GetObject());
      if (!fTriggerLut)
      {
        AliError("fTriggerLut not of the expected type !!!");
      }
    }
    else
    {
      AliError("Could not get trigger lut !");
    }
  }
  return fTriggerLut;
}



