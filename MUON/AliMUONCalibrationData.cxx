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
#include "AliDCSValue.h"
#include "AliLog.h"
#include "AliMpDCSNamer.h"
#include "AliMUONConstants.h"
#include "AliMUONGlobalCrateConfig.h"
#include "AliMUONRegionalTriggerConfig.h"
#include "AliMUONRejectList.h"
#include "AliMUONTriggerEfficiencyCells.h"
#include "AliMUONTriggerLut.h"
#include "AliMUONVCalibParam.h"
#include "AliMUONVStore.h"
#include "AliMUONVStore.h"

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

AliMUONVStore* AliMUONCalibrationData::fgBypassPedestals(0x0);
AliMUONVStore* AliMUONCalibrationData::fgBypassGains(0x0);

//_____________________________________________________________________________
AliMUONCalibrationData::AliMUONCalibrationData(Int_t runNumber, 
                                               Bool_t deferredInitialization) 
: TObject(), 
fIsValid(kTRUE),
fRunNumber(runNumber), 
fGains(0x0), 
fPedestals(0x0),
fHV(0x0),
fTriggerDCS(0x0),
fLocalTriggerBoardMasks(0x0),
fRegionalTriggerConfig(0x0),
fGlobalTriggerCrateConfig(0x0),
fTriggerLut(0x0),
fTriggerEfficiency(0x0),
fCapacitances(0x0),
fNeighbours(0x0),
fOccupancyMap(0x0),
fRejectList(0x0),
fConfig(0x0)
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
    OccupancyMap();
    RejectList();
    HV();
    TriggerDCS();
    LocalTriggerBoardMasks(0);
    RegionalTriggerConfig();
    GlobalTriggerCrateConfig();
    TriggerLut();
    TriggerEfficiency();
    Capacitances();
    Neighbours();
    Config();
  }
}

//_____________________________________________________________________________
AliMUONCalibrationData::~AliMUONCalibrationData()
{
  /// Destructor. Note that we're the owner of our pointers if the OCDB cache
  /// is not set. Otherwise the cache is supposed to take care of them...
  if (!(AliCDBManager::Instance()->GetCacheFlag())) Reset();
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
AliMUONCalibrationData::CreateCapacitances(Int_t runNumber, Int_t* startOfValidity)
{
  /// Create capa store from OCDB for a given run
  
  return dynamic_cast<AliMUONVStore*>(CreateObject(runNumber,"MUON/Calib/Capacitances",startOfValidity));
}

//_____________________________________________________________________________
AliMUONVStore*
AliMUONCalibrationData::CreateGains(Int_t runNumber, Int_t* startOfValidity)
{
  /// Create a new gain store from the OCDB for a given run
  return dynamic_cast<AliMUONVStore*>(CreateObject(runNumber,"MUON/Calib/Gains",startOfValidity));
}

//_____________________________________________________________________________
AliMUONGlobalCrateConfig*
AliMUONCalibrationData::CreateGlobalTriggerCrateConfig(Int_t runNumber, Int_t* startOfValidity)
{
  /// Create the internal store for GlobalTriggerCrateConfig from OCDB
  
  return dynamic_cast<AliMUONGlobalCrateConfig*>(CreateObject(runNumber,"MUON/Calib/GlobalTriggerCrateConfig",startOfValidity));
}


//______________________________________________________________________________
void AliMUONCalibrationData::PatchHVValues(TObjArray& values,
                                           Int_t& nbelowready,
                                           Int_t& noff,
                                           Int_t& ntrips,
                                           Int_t& neor,
                                           TString* msg)
{
  /// We do here a little bit of massaging of the HV values, if needed.
  ///
  /// The main point is to "gather" the values near the end of the run (the last XX seconds)
  /// to avoid the ramp-down before end-of-run syndrom...
  ///
  /// The rest is more of a debug/expert tool to have closer look at trends
  ///
  
  Double_t HVBEAMTUNING(1300); 
  
  Bool_t eorProblem(kFALSE);
  
  UInt_t mergeDelay(300); // in seconds
  
  // First start by removing values within the last mergeDelay seconds, keeping only
  // the last one
  
  AliDCSValue* last = static_cast<AliDCSValue*>(values.At(values.GetLast()));
  
  Int_t* toberemoved = new Int_t[values.GetLast()+1];
  
  memset(toberemoved,0,(values.GetLast()+1)*sizeof(Int_t));
  
  Int_t ntoberemoved(0);
  
  for ( Int_t i = values.GetLast()-1; i > 0; --i ) 
  {
    AliDCSValue* val = static_cast<AliDCSValue*>(values.At(i));
    
    if ( last->GetTimeStamp() - val->GetTimeStamp() < mergeDelay )
    {
      toberemoved[i]=1;
      ++ntoberemoved;
    }
  }
  
  if (ntoberemoved)
  {
    // ok, we have some values within the same mergeDelay seconds
    // we'll "merge" them by taking the last one, except if
    // the last one is below HVBEAMTUNING, in which case we
    // remove that one too (meaning the ramp-down was requesting
    // before the end-of-run)
    
    if ( last->GetFloat() < HVBEAMTUNING ) 
    {
      eorProblem=kTRUE;
      if (msg) *msg = "ERROR RAMP-DOWN BEFORE EOR";
      toberemoved[values.GetLast()]=1;
    }
    
    for ( Int_t i = 0; i <= values.GetLast(); ++i ) 
    {
      if ( toberemoved[i] ) values.RemoveAt(i);
    }
    
    values.Compress();
    
  }
  
  delete[] toberemoved;
  
  if (eorProblem)
  {
    ++neor;
    return;
  }
  
  // now for the rest of the diagnosis
  
  Int_t ntmpoff(0);
  Int_t ntmpready(0);
  
  for ( Int_t i = 0; i <= values.GetLast(); ++i ) 
  {
    AliDCSValue* val = static_cast<AliDCSValue*>(values.At(i));
    
    if ( val->GetFloat() < AliMpDCSNamer::TrackerHVOFF() ) 
    {
      ++ntmpoff;
    }
    else if ( val->GetFloat() < HVBEAMTUNING )
    {
      ++ntmpready;
    }
  }
  
  if ( ntmpoff )
  {
    if ( ntmpoff == values.GetLast()+1 ) 
    {
      if (msg) *msg = "ERROR HV OFF";
      ++noff;
    }
    else
    {
      if (msg) *msg = "ERROR TRIP";
      ++ntrips;    
    }
  }
  
  if ( ntmpready == values.GetLast()+1 ) 
  {
    if (msg) *msg = "ERROR BELOW READY";    
  }      
  
  if (ntmpready) ++nbelowready;
}

//_____________________________________________________________________________
TMap*
AliMUONCalibrationData::CreateHV(Int_t runNumber, Int_t* startOfValidity, Bool_t patched)
{
  /// Create a new HV map from the OCDB for a given run
  TMap* hvMap = dynamic_cast<TMap*>(CreateObject(runNumber,"MUON/Calib/HV",startOfValidity));
  if (patched)
  {
    TIter next(hvMap);
    TObjString* hvChannelName;
    
    while ( ( hvChannelName = static_cast<TObjString*>(next()) ) )
    {
      TString name(hvChannelName->String());
      
      if ( name.Contains("sw") ) continue; // skip switches
      
      TPair* hvPair = static_cast<TPair*>(hvMap->FindObject(name.Data()));
      TObjArray* values = static_cast<TObjArray*>(hvPair->Value());
      if (!values)
      {
        AliErrorClass(Form("Could not get values for alias %s",name.Data()));
      }
      else
      {
        int nbelowready(0);
        int noff(0);
        int ntrips(0);
        int neor(0);
        
        PatchHVValues(*values,nbelowready,noff,ntrips,neor,0x0);
        if (neor)
        {
          nbelowready=noff=ntrips=neor=0;
          PatchHVValues(*values,nbelowready,noff,ntrips,neor,0x0);
          if (neor)
          {
            AliErrorClass("neor is not null after PatchHVValue ! This is serious !");
          }
        }
      }
    }
    
  }
  return hvMap;
}

//_____________________________________________________________________________
TMap*
AliMUONCalibrationData::CreateTriggerDCS(Int_t runNumber, Int_t* startOfValidity)
{
  /// Create a new Trigger HV and curent map from the OCDB for a given run
  return dynamic_cast<TMap*>(CreateObject(runNumber,"MUON/Calib/TriggerDCS",startOfValidity));
}

//_____________________________________________________________________________
AliMUONVStore*
AliMUONCalibrationData::CreateLocalTriggerBoardMasks(Int_t runNumber, Int_t* startOfValidity)
{
  /// Get the internal store for LocalTriggerBoardMasks from OCDB
  
  return dynamic_cast<AliMUONVStore*>(CreateObject(runNumber,"MUON/Calib/LocalTriggerBoardMasks",startOfValidity));
}

//_____________________________________________________________________________
AliMUONVStore*
AliMUONCalibrationData::CreateNeighbours(Int_t runNumber, Int_t* startOfValidity)
{
  /// Create a neighbour store from the OCDB for a given run
  return dynamic_cast<AliMUONVStore*>(CreateObject(runNumber,"MUON/Calib/Neighbours",startOfValidity));
}

//_____________________________________________________________________________
TObject*
AliMUONCalibrationData::CreateObject(Int_t runNumber, const char* path, Int_t* startOfValidity)
{
  /// Access the CDB for a given path (e.g. MUON/Calib/Pedestals),
  /// and return the corresponding TObject.
  
  AliCodeTimerAutoClass(Form("%d : %s",runNumber,path),0);
  
  AliCDBManager* man = AliCDBManager::Instance();
  
  AliCDBEntry* entry =  man->Get(path,runNumber);
  
  if (entry)
  {
		if ( startOfValidity ) *startOfValidity = entry->GetId().GetFirstRun();
		
    TObject* object = entry->GetObject();
    if (!(man->GetCacheFlag()))
    {
      entry->SetOwner(kFALSE);
      delete entry;      
    }
    return object;
  }
	else
	{
		if ( startOfValidity )  *startOfValidity = AliCDBRunRange::Infinity();
  }
	
  {
    
    AliCodeTimerAutoClass(Form("Failed to get %s for run %d",path,runNumber),1);

  }
  
  return 0x0;
}

//_____________________________________________________________________________
AliMUONVStore*
AliMUONCalibrationData::CreateOccupancyMap(Int_t runNumber, Int_t* startOfValidity)
{
  /// Create a new occupancy map store from the OCDB for a given run
  return dynamic_cast<AliMUONVStore*>(CreateObject(runNumber,"MUON/Calib/OccupancyMap",startOfValidity));
}

//_____________________________________________________________________________
AliMUONRejectList*
AliMUONCalibrationData::CreateRejectList(Int_t runNumber, Int_t* startOfValidity)
{
  /// Create a new rejectlist store from the OCDB for a given run
  return dynamic_cast<AliMUONRejectList*>(CreateObject(runNumber,"MUON/Calib/RejectList",startOfValidity));
}

//_____________________________________________________________________________
AliMUONVStore*
AliMUONCalibrationData::CreatePedestals(Int_t runNumber, Int_t* startOfValidity)
{
  /// Create a new pedestal store from the OCDB for a given run
  return dynamic_cast<AliMUONVStore*>(CreateObject(runNumber,"MUON/Calib/Pedestals",startOfValidity));
}

//_____________________________________________________________________________
AliMUONVStore*
AliMUONCalibrationData::CreateConfig(Int_t runNumber, Int_t* startOfValidity)
{
  /// Create a new config store from the OCDB for a given run
  return dynamic_cast<AliMUONVStore*>(CreateObject(runNumber,"MUON/Calib/Config",startOfValidity));
}


//_____________________________________________________________________________
AliMUONRegionalTriggerConfig*
AliMUONCalibrationData::CreateRegionalTriggerConfig(Int_t runNumber, Int_t* startOfValidity)
{
  /// Create the internal store for RegionalTriggerConfig from OCDB
  
  return dynamic_cast<AliMUONRegionalTriggerConfig*>(CreateObject(runNumber,"MUON/Calib/RegionalTriggerConfig",startOfValidity));
}

//_____________________________________________________________________________
AliMUONTriggerEfficiencyCells* 
AliMUONCalibrationData::CreateTriggerEfficiency(Int_t runNumber, Int_t* startOfValidity)
{
  /// Create trigger efficiency object from OCBD
  
  return dynamic_cast<AliMUONTriggerEfficiencyCells*>(CreateObject(runNumber,"MUON/Calib/TriggerEfficiency",startOfValidity));
}

//_____________________________________________________________________________
AliMUONTriggerLut* 
AliMUONCalibrationData::CreateTriggerLut(Int_t runNumber, Int_t* startOfValidity)
{
  /// Create trigger LUT from OCDB
  
  return dynamic_cast<AliMUONTriggerLut*>(CreateObject(runNumber,"MUON/Calib/TriggerLut",startOfValidity));
}

//_____________________________________________________________________________
AliMUONVStore*
AliMUONCalibrationData::Gains() const
{
  /// Create (if needed) and return the internal store for gains.
  if (fgBypassGains) return fgBypassGains;
  
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
AliMUONGlobalCrateConfig* 
AliMUONCalibrationData::GlobalTriggerCrateConfig() const
{
  /// Return the config for the global trigger board.
  
  if (!fGlobalTriggerCrateConfig)
  {
    fGlobalTriggerCrateConfig = CreateGlobalTriggerCrateConfig(fRunNumber);
  }
  return fGlobalTriggerCrateConfig;
}


//_____________________________________________________________________________
TMap*
AliMUONCalibrationData::HV(Bool_t patched) const
{
  /// Return the calibration for a given (detElemId, manuId) pair
  
  if (!fHV)
  {
    fHV = CreateHV(fRunNumber,0,patched);
  }
  return fHV;
}

//_____________________________________________________________________________
TMap*
AliMUONCalibrationData::TriggerDCS() const
{
  /// Return the calibration for a given (detElemId, manuId) pair
  
  if (!fTriggerDCS)
  {
    fTriggerDCS = CreateTriggerDCS(fRunNumber);
  }
  return fTriggerDCS;
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
AliMUONCalibrationData::OccupancyMap() const
{
  /// Get occupancy map
  if (!fOccupancyMap)
  {
    fOccupancyMap = CreateOccupancyMap(fRunNumber);
  }
  return fOccupancyMap;
}

//_____________________________________________________________________________
AliMUONRejectList*
AliMUONCalibrationData::RejectList() const
{
  /// Get reject list
  if (!fRejectList)
  {
    fRejectList = CreateRejectList(fRunNumber);
  }
  return fRejectList;
}

//_____________________________________________________________________________
void
AliMUONCalibrationData::BypassStores(AliMUONVStore* ped, AliMUONVStore* gain)
{
  /// Force the use of those pedestals and gains
  fgBypassPedestals = ped;
  fgBypassGains = gain;
  
}

//_____________________________________________________________________________
AliMUONVStore*
AliMUONCalibrationData::Pedestals() const
{
  /// Return pedestals
  
  if (fgBypassPedestals) return fgBypassPedestals;
  
  if (!fPedestals)
  {
    fPedestals = CreatePedestals(fRunNumber);
  }
  return fPedestals;
}

//_____________________________________________________________________________
AliMUONVStore*
AliMUONCalibrationData::Config() const
{
  /// Return config
  
  if (!fConfig)
  {
    fConfig = CreateConfig(fRunNumber);
  }
  return fConfig;
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
  << " fConfig=" << fConfig
  << " fHV=" << fHV
  << " fTriggerDCS=" << fTriggerDCS
  << " fLocalTriggerBoardMasks=" << fLocalTriggerBoardMasks
  << " fRegionalTriggerConfig=" << fRegionalTriggerConfig
  << " fGlobalTriggerCrateConfig=" << fGlobalTriggerCrateConfig
  << " fTriggerLut=" << fTriggerLut
  << endl;
}


//_____________________________________________________________________________
AliMUONRegionalTriggerConfig* 
AliMUONCalibrationData::RegionalTriggerConfig() const
{
  /// Return the config for the regional trigger board.
  
  if (!fRegionalTriggerConfig)
  {
    fRegionalTriggerConfig = CreateRegionalTriggerConfig(fRunNumber);
    }
  return fRegionalTriggerConfig;
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

  AliCodeTimerAuto("",0);
  
  delete fConfig;
  fConfig = 0x0;
  delete fPedestals;
  fPedestals = 0x0;
  delete fGains;
  fGains = 0x0;
  delete fHV;
  fHV = 0x0;
  delete fTriggerDCS;
  fTriggerDCS = 0x0;
  delete fLocalTriggerBoardMasks;
  fLocalTriggerBoardMasks = 0x0;
  delete fRegionalTriggerConfig;
  fRegionalTriggerConfig = 0x0;
  delete fGlobalTriggerCrateConfig;
  fGlobalTriggerCrateConfig = 0x0;
  
  delete fTriggerLut;
  fTriggerLut = 0x0;
  delete fTriggerEfficiency;
  fTriggerEfficiency = 0x0;
  delete fCapacitances;
  fCapacitances = 0x0;
  delete fNeighbours;
  fNeighbours = 0x0;
}

//_____________________________________________________________________________
void
AliMUONCalibrationData::Check(Int_t runNumber)
{
  /// Self-check to see if we can read all data for a given run 
  /// from the current OCDB...
  
  if ( ! CreateCapacitances(runNumber) )
  {
    AliErrorClass("Could not read capacitances");
  }
  else
  {
    AliInfoClass("Capacitances read OK");
  }

  if ( ! CreateGains(runNumber) ) 
  {
    AliErrorClass("Could not read gains");
  }
  else
  {
    AliInfoClass("Gains read OK");
  }

  if ( ! CreateGlobalTriggerCrateConfig(runNumber) ) 
  {
    AliErrorClass("Could not read Trigger Crate Config");
  }
  else
  {
    AliInfoClass("TriggerBoardMasks read OK");
  }

  if ( !  CreateHV(runNumber) )
  {
    AliErrorClass("Could not read HV");
  }
  else
  {
    AliInfoClass("HV read OK");
  }

  if ( !  CreateTriggerDCS(runNumber) )
  {
    AliErrorClass("Could not read Trigger HV and Currents");
  }
  else
  {
    AliInfoClass("Trigger HV and Currents read OK");
  }

  if ( ! CreateNeighbours(runNumber) )
  {
    AliErrorClass("Could not read Neighbours");
  }
  else
  {
    AliInfoClass("Neighbours read OK");
  }

  if ( !  CreateLocalTriggerBoardMasks(runNumber) )
  {
    AliErrorClass("Could not read LocalTriggerBoardMasks");
  }
  else
  {
    AliInfoClass("LocalTriggerBoardMasks read OK");
  }
  
  if ( ! CreatePedestals(runNumber) )
  {
    AliErrorClass("Could not read pedestals");
  }
  else
  {
    AliInfoClass("Pedestals read OK");
  }

  if ( ! CreateConfig(runNumber) )
  {
    AliErrorClass("Could not read config");
  }
  else
  {
    AliInfoClass("Config read OK");
  }
  
  if ( ! CreateRegionalTriggerConfig(runNumber) )
  {
    AliErrorClass("Could not read RegionalTriggerConfig");
  }
  else
  {
    AliInfoClass("RegionalTriggerBoardMasks read OK");
  }
  
  if ( ! CreateTriggerLut(runNumber) )
  {
    AliErrorClass("Could not read TriggerLut");
  }
  else
  {
    AliInfoClass("TriggerLut read OK");
  }

  if ( ! CreateTriggerEfficiency(runNumber) )
  {
    AliErrorClass("Could not read TriggerEfficiency");
  }
  else    
  {
    AliInfoClass("TriggerEfficiency read OK");
  }
}


