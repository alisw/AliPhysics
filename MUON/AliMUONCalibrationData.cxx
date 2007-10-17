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
#include "AliCodeTimer.h"
#include "AliLog.h"
#include "AliMUONTriggerEfficiencyCells.h"
#include "AliMUONTriggerLut.h"
#include "AliMUONVStore.h"
#include "AliMUONVStore.h"
#include "AliMUONVCalibParam.h"
#include <Riostream.h>
#include <TClass.h>
#include <TMap.h>

//-----------------------------------------------------------------------------
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
//-----------------------------------------------------------------------------

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
fTriggerEfficiency(0x0),
fCapacitances(0x0),
fNeighbours(0x0)
{
/// Default ctor.

  // If deferredInitialization is false, we read *all* calibrations
  // at once.
  // So when using this class to access only one kind of calibrations (e.g.
  // only pedestals), you should put deferredInitialization to kTRUE, which
  // will instruct this object to fetch the data only when neeeded.

  if ( deferredInitialization == kFALSE )
  {
    Gains();
    Pedestals();
    HV();
    LocalTriggerBoardMasks(0);
    RegionalTriggerBoardMasks(0);
    GlobalTriggerBoardMasks();
    TriggerLut();
    TriggerEfficiency();
    Capacitances();
    Neighbours();
  }
}

//_____________________________________________________________________________
AliMUONCalibrationData::~AliMUONCalibrationData()
{
  /// Destructor. Note that we're the owner of our pointers.
  Reset();
}

//_____________________________________________________________________________
AliMUONVStore*
AliMUONCalibrationData::Capacitances() const
{
  /// Create (if needed) and return the internal store for capacitances.
  
  if (!fCapacitances)
  {
    fCapacitances = CreateCapacitances(fRunNumber);
  }
  return fCapacitances;
}

//_____________________________________________________________________________
AliMUONVStore*
AliMUONCalibrationData::CreateCapacitances(Int_t runNumber)
{
  /// Create capa store from OCDB for a given run
  
  return dynamic_cast<AliMUONVStore*>(CreateObject(runNumber,"MUON/Calib/Capacitances"));
}

//_____________________________________________________________________________
AliMUONVStore*
AliMUONCalibrationData::CreateGains(Int_t runNumber)
{
  /// Create a new gain store from the OCDB for a given run
  return dynamic_cast<AliMUONVStore*>(CreateObject(runNumber,"MUON/Calib/Gains"));
}

//_____________________________________________________________________________
AliMUONVCalibParam*
AliMUONCalibrationData::CreateGlobalTriggerBoardMasks(Int_t runNumber)
{
  /// Create the internal store for GlobalTriggerBoardMasks from OCDB
  
  return dynamic_cast<AliMUONVCalibParam*>(CreateObject(runNumber,"MUON/Calib/GlobalTriggerBoardMasks"));
}

//_____________________________________________________________________________
TMap*
AliMUONCalibrationData::CreateHV(Int_t runNumber)
{
  /// Create a new HV map from the OCDB for a given run
  return dynamic_cast<TMap*>(CreateObject(runNumber,"MUON/Calib/HV"));
}

//_____________________________________________________________________________
AliMUONVStore*
AliMUONCalibrationData::CreateLocalTriggerBoardMasks(Int_t runNumber)
{
  /// Get the internal store for LocalTriggerBoardMasks from OCDB
  
  return dynamic_cast<AliMUONVStore*>(CreateObject(runNumber,"MUON/Calib/LocalTriggerBoardMasks"));
}

//_____________________________________________________________________________
AliMUONVStore*
AliMUONCalibrationData::CreateNeighbours(Int_t runNumber)
{
  /// Create a neighbour store from the OCDB for a given run
  return dynamic_cast<AliMUONVStore*>(CreateObject(runNumber,"MUON/Calib/Neighbours"));
}

//_____________________________________________________________________________
TObject*
AliMUONCalibrationData::CreateObject(Int_t runNumber, const char* path)
{
  /// Access the CDB for a given path (e.g. MUON/Calib/Pedestals),
  /// and return the corresponding TObject.
  
  AliCodeTimerAutoClass(Form("%d : %s",runNumber,path));
  
  AliCDBManager* man = AliCDBManager::Instance();
  
  Bool_t undefStorage(kFALSE);
  
  if ( !man->IsDefaultStorageSet() )
  {
    TString storage("local://$ALICE_ROOT");
    AliInfoClass(Form("CDB Storage not set. Will use %s for MUON stuff",storage.Data()));
    man->SetDefaultStorage(storage.Data());
    undefStorage = kTRUE;
  }
  
  Bool_t cacheStatus = man->GetCacheFlag();
  
  man->SetCacheFlag(kFALSE);
  
  AliCDBEntry* entry =  AliCDBManager::Instance()->Get(path,runNumber);
  
  man->SetCacheFlag(cacheStatus);
  
  if (entry)
  {
    TObject* object = entry->GetObject();
    entry->SetOwner(kFALSE);
    delete entry;
    return object;
  }
  
  if ( undefStorage ) 
  {
    man->UnsetDefaultStorage();
  }
  
  return 0x0;
}

//_____________________________________________________________________________
AliMUONVStore*
AliMUONCalibrationData::CreatePedestals(Int_t runNumber)
{
  /// Create a new pedestal store from the OCDB for a given run
  return dynamic_cast<AliMUONVStore*>(CreateObject(runNumber,"MUON/Calib/Pedestals"));
}

//_____________________________________________________________________________
AliMUONVStore*
AliMUONCalibrationData::CreateRegionalTriggerBoardMasks(Int_t runNumber)
{
  /// Create the internal store for RegionalTriggerBoardMasks from OCDB
  
  return dynamic_cast<AliMUONVStore*>(CreateObject(runNumber,"MUON/Calib/RegionalTriggerBoardMasks"));
}

//_____________________________________________________________________________
AliMUONTriggerEfficiencyCells* 
AliMUONCalibrationData::CreateTriggerEfficiency(Int_t runNumber)
{
  /// Create trigger efficiency object from OCBD
  
  return dynamic_cast<AliMUONTriggerEfficiencyCells*>(CreateObject(runNumber,"MUON/Calib/TriggerEfficiency"));
}

//_____________________________________________________________________________
AliMUONTriggerLut* 
AliMUONCalibrationData::CreateTriggerLut(Int_t runNumber)
{
  /// Create trigger LUT from OCDB
  
  return dynamic_cast<AliMUONTriggerLut*>(CreateObject(runNumber,"MUON/Calib/TriggerLut"));
}

//_____________________________________________________________________________
AliMUONVStore*
AliMUONCalibrationData::Gains() const
{
  /// Create (if needed) and return the internal store for gains.
  if (!fGains)
  {
    fGains = CreateGains(fRunNumber);
  }
  return fGains;
}

//_____________________________________________________________________________
AliMUONVCalibParam*
AliMUONCalibrationData::Gains(Int_t detElemId, Int_t manuId) const
{
/// Return the gains for a given (detElemId, manuId) pair
/// Note that, unlike the DeadChannel case, if the result is 0x0, that's an
/// error (meaning that we should get gains for all channels).

  AliMUONVStore* gains = Gains();
  if (!gains)
  {
    return 0x0;
  }
  
  return static_cast<AliMUONVCalibParam*>(gains->FindObject(detElemId,manuId));
}

//_____________________________________________________________________________
AliMUONVCalibParam* 
AliMUONCalibrationData::GlobalTriggerBoardMasks() const
{
  /// Return the masks for the global trigger board.
  
  if (!fGlobalTriggerBoardMasks)
  {
    fGlobalTriggerBoardMasks = CreateGlobalTriggerBoardMasks(fRunNumber);
  }
  return fGlobalTriggerBoardMasks;
}

//_____________________________________________________________________________
TMap*
AliMUONCalibrationData::HV() const
{
  /// Return the calibration for a given (detElemId, manuId) pair
  
  if (!fHV)
  {
    fHV = CreateHV(fRunNumber);
  }
  return fHV;
}

//_____________________________________________________________________________
AliMUONVStore*
AliMUONCalibrationData::Neighbours() const
{
  /// Create (if needed) and return the internal store for neighbours.
  if (!fNeighbours)
  {
    fNeighbours = CreateNeighbours(fRunNumber);
  }
  return fNeighbours;
}

//_____________________________________________________________________________
AliMUONVCalibParam* 
AliMUONCalibrationData::LocalTriggerBoardMasks(Int_t localBoardNumber) const
{
/// Return the masks for a given trigger local board.

  if (!fLocalTriggerBoardMasks)
  {
    fLocalTriggerBoardMasks = CreateLocalTriggerBoardMasks(fRunNumber);
  }

  if ( fLocalTriggerBoardMasks ) 
  {
    AliMUONVCalibParam* ltbm = 
      static_cast<AliMUONVCalibParam*>(fLocalTriggerBoardMasks->FindObject(localBoardNumber));
    if (!ltbm)
    {
      AliError(Form("Could not get mask for localBoardNumber=%d",localBoardNumber));
    }
    return ltbm;  
  }
  return 0x0;
}

//_____________________________________________________________________________
AliMUONVStore*
AliMUONCalibrationData::Pedestals() const
{
  /// Return pedestals
  if (!fPedestals)
  {
    fPedestals = CreatePedestals(fRunNumber);
  }
  return fPedestals;
}

//_____________________________________________________________________________
AliMUONVCalibParam*
AliMUONCalibrationData::Pedestals(Int_t detElemId, Int_t manuId) const
{
  /// Return the pedestals for a given (detElemId, manuId) pair.
  /// A return value of 0x0 is considered an error, meaning we should get
  /// pedestals for all channels.
  
  AliMUONVStore* pedestals = Pedestals();
  if (!pedestals) 
  {
    return 0x0;
  }
  
  return static_cast<AliMUONVCalibParam*>(pedestals->FindObject(detElemId,manuId));
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
AliMUONVCalibParam* 
AliMUONCalibrationData::RegionalTriggerBoardMasks(Int_t regionalBoardNumber) const
{
/// Return the masks for a given trigger regional board.

  if ( !fRegionalTriggerBoardMasks ) 
  {
    fRegionalTriggerBoardMasks = CreateRegionalTriggerBoardMasks(fRunNumber);
  }
  
  if ( fRegionalTriggerBoardMasks ) 
  {
    AliMUONVCalibParam* rtbm = 
      static_cast<AliMUONVCalibParam*>(fRegionalTriggerBoardMasks->FindObject(regionalBoardNumber));
    
    if (!rtbm)
    {
      AliError(Form("Could not get mask for regionalBoard index=%d",regionalBoardNumber));
    }
    return rtbm;  
  }
  return 0x0;
}

//_____________________________________________________________________________
AliMUONTriggerEfficiencyCells*
AliMUONCalibrationData::TriggerEfficiency() const
{
/// Return the trigger efficiency.

  if (!fTriggerEfficiency)
  {
    fTriggerEfficiency = CreateTriggerEfficiency(fRunNumber);
  }
  return fTriggerEfficiency;
}


//_____________________________________________________________________________
AliMUONTriggerLut*
AliMUONCalibrationData::TriggerLut() const
{
/// Return the trigger look up table.

  if (!fTriggerLut)
  {
    fTriggerLut = CreateTriggerLut(fRunNumber);
  }
  return fTriggerLut;
}

//_____________________________________________________________________________
void
AliMUONCalibrationData::Reset()
{
/// Reset all data

  delete fPedestals;
  fPedestals = 0x0;
  delete fGains;
  fGains = 0x0;
  delete fHV;
  fHV = 0x0;
  delete fLocalTriggerBoardMasks;
  fLocalTriggerBoardMasks = 0x0;
  delete fRegionalTriggerBoardMasks;
  fRegionalTriggerBoardMasks = 0x0;
  delete fGlobalTriggerBoardMasks;
  fGlobalTriggerBoardMasks = 0x0;
  delete fTriggerLut;
  fTriggerLut = 0x0;
  delete fTriggerEfficiency;
  fTriggerEfficiency = 0x0;
  delete fCapacitances;
  fCapacitances = 0x0;
  delete fNeighbours;
  fNeighbours = 0x0;
}



