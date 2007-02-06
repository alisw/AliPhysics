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
#include "TMap.h"

/// \class AliMUONCalibrationData
///
/// For the moment, this class stores pedestals, gains, hv (for tracker)
/// and lut, masks and efficiencies (for trigger) that are fetched from the CDB.
///
/// This class is to be considered as a convenience class.
/// Its aim is to ease retrieval of calibration data from the 
/// condition database.
///
/// It acts as a "facade" to a bunch of underlying 
/// containers/calibration classes.
///
/// \author Laurent Aphecetche

/// \cond CLASSIMP
ClassImp(AliMUONCalibrationData)
/// \endcond

//_____________________________________________________________________________
AliMUONCalibrationData::AliMUONCalibrationData(Int_t runNumber, 
                                               Bool_t deferredInitialization) 
: TObject(), 
fIsValid(kTRUE),
fRunNumber(runNumber), 
fGains(0x0), 
fPedestals(0x0),
fHV(0x0),
fLocalTriggerBoardMasks(0x0),
fRegionalTriggerBoardMasks(0x0),
fGlobalTriggerBoardMasks(0x0),
fTriggerLut(0x0),
fTriggerEfficiency(0x0)
{
/// Default ctor.

  // If deferredInitialization is false, we read *all* calibrations
  // at once.
  // So when using this class to access only one kind of calibrations (e.g.
  // only pedestals), you should put deferredInitialization to kTRUE, which
  // will instruct this object to fetch the data only when neeeded.

  if ( deferredInitialization == kFALSE )
  {
    OnDemandGains();
    OnDemandPedestals();
    OnDemandHV();
    OnDemandLocalTriggerBoardMasks();
    OnDemandRegionalTriggerBoardMasks();
    OnDemandGlobalTriggerBoardMasks();
    OnDemandTriggerLut();
    OnDemandTriggerEfficiency();
  }
}

//_____________________________________________________________________________
AliMUONCalibrationData::~AliMUONCalibrationData()
{
/// Destructor. Note that we're the owner of our pointers.

  delete fPedestals;
  delete fGains;
  delete fHV;
  delete fLocalTriggerBoardMasks;
  delete fRegionalTriggerBoardMasks;
  delete fGlobalTriggerBoardMasks;
  delete fTriggerLut;
  delete fTriggerEfficiency;
}

//_____________________________________________________________________________
TMap*
AliMUONCalibrationData::HV() const
{
/// Return the calibration for a given (detElemId, manuId) pair

  return OnDemandHV();
}

//_____________________________________________________________________________
TMap*
AliMUONCalibrationData::OnDemandHV() const
{
/// Create (if needed) and return the internal store for DeadChannels.

  if (!fHV)
  {
    AliCDBEntry* entry = GetEntry("MUON/Calib/HV");
    if (entry)
    {
      fHV = dynamic_cast<TMap*>(entry->GetObject());
      if (!fHV)
      {
        AliError("fHV not of the expected type !!!");
      }
    }
    else
    {
      AliError("Could not get HV values !");
    }
  }
  return fHV;
}

//_____________________________________________________________________________
AliCDBEntry*
AliMUONCalibrationData::GetEntry(const char* path) const
{
/// Access the CDB for a given path (e.g. MUON/Calib/Pedestals),
/// and return the corresponding CDBEntry.

  return AliCDBManager::Instance()->Get(path,fRunNumber);
}

//_____________________________________________________________________________
AliMUONVCalibParam*
AliMUONCalibrationData::Gains(Int_t detElemId, Int_t manuId) const
{
/// Return the gains for a given (detElemId, manuId) pair
/// Note that, unlike the DeadChannel case, if the result is 0x0, that's an
/// error (meaning that we should get gains for all channels).

  AliMUONV2DStore* gains = Gains();
  if (!gains)
  {
    AliError("Could not get gains");
    return 0x0;
  }
  
  AliMUONVCalibParam* g = 
    static_cast<AliMUONVCalibParam*>(gains->Get(detElemId,manuId));
  if (!g)
  {
    AliError(Form("Could not get gain for detElemId=%d manuId=%d ",
                    detElemId,manuId));
  }
  return g;
}

//_____________________________________________________________________________
AliMUONV2DStore*
AliMUONCalibrationData::Gains() const
{
  /// Create (if needed) and return the internal store for gains.
  return OnDemandGains();
}

//_____________________________________________________________________________
AliMUONV2DStore*
AliMUONCalibrationData::OnDemandGains() const
{
/// Create (if needed) and return the internal store for gains.

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
/// Return the masks for the global trigger board.

  return OnDemandGlobalTriggerBoardMasks();
}

//_____________________________________________________________________________
AliMUONVCalibParam*
AliMUONCalibrationData::OnDemandGlobalTriggerBoardMasks() const
{
/// Create (if needed) and return the internal store for GlobalTriggerBoardMasks.

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
/// Return the masks for a given trigger local board.

  AliMUONV1DStore* store = OnDemandLocalTriggerBoardMasks();
  if (!store)
  {
    AliError("Could not get LocalTriggerBoardMasks");
    return 0x0;
  }
  
  AliMUONVCalibParam* ltbm = 
    static_cast<AliMUONVCalibParam*>(store->Get(localBoardNumber));
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
/// Create (if needed) and return the internal store for LocalTriggerBoardMasks.

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
/// Create (if needed) and return the internal storage for pedestals.

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
/// A very basic dump of our guts.

  cout << "RunNumber " << RunNumber()
  << " fGains=" << fGains
  << " fPedestals=" << fPedestals
  << " fHV=" << fHV
  << " fLocalTriggerBoardMasks=" << fLocalTriggerBoardMasks
  << " fRegionalTriggerBoardMasks=" << fRegionalTriggerBoardMasks
  << " fGlobalTriggerBoardMasks=" << fGlobalTriggerBoardMasks
  << " fTriggerLut=" << fTriggerLut
  << endl;
}

//_____________________________________________________________________________
AliMUONV2DStore*
AliMUONCalibrationData::Pedestals() const
{
  /// Return pedestals
  return OnDemandPedestals();
}

//_____________________________________________________________________________
AliMUONVCalibParam*
AliMUONCalibrationData::Pedestals(Int_t detElemId, Int_t manuId) const
{
/// Return the pedestals for a given (detElemId, manuId) pair.
/// A return value of 0x0 is considered an error, meaning we should get
/// pedestals for all channels.

  AliMUONV2DStore* pedestals = OnDemandPedestals();
  if (!pedestals) 
  {
    AliError(Form("Did not find pedestals for DE %d manuId %d",
                  detElemId,manuId));
    return 0x0;
  }
  
  AliMUONVCalibParam* ped = 
    static_cast<AliMUONVCalibParam*>(pedestals->Get(detElemId,manuId));
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
/// Return the masks for a given trigger regional board.

  AliMUONV1DStore* store = OnDemandRegionalTriggerBoardMasks();
  
  if (!store)
  {
    AliError("Could not get RegionalTriggerBoardMasks");
    return 0x0;
  }
  
  AliMUONVCalibParam* rtbm = 
    static_cast<AliMUONVCalibParam*>(store->Get(index));
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
/// Create (if needed) and return the internal store for RegionalTriggerBoardMasks.

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
/// Return the trigger efficiency.

  return OnDemandTriggerEfficiency();
}

//_____________________________________________________________________________
AliMUONTriggerEfficiencyCells* 
AliMUONCalibrationData::OnDemandTriggerEfficiency() const
{
/// \todo: add comment

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
/// Return the trigger look up table.

  return OnDemandTriggerLut();
}

//_____________________________________________________________________________
AliMUONTriggerLut* 
AliMUONCalibrationData::OnDemandTriggerLut() const
{
/// \todo: add comment

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



