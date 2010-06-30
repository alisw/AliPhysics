/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup calib
/// \class AliMUONCalibrationData
/// \brief Single entry point to access MUON calibration data.
/// 
//  Author Laurent Aphecetche

#ifndef ALIMUONCALIBRATIONDATA_H
#define ALIMUONCALIBRATIONDATA_H

#ifndef ROOT_TObject
#include "TObject.h"
#endif

class AliCDBEntry;
class AliMUONGlobalCrateConfig;
class AliMUONRegionalTriggerConfig;
class AliMUONRejectList;
class AliMUONTriggerEfficiencyCells;
class AliMUONTriggerLut;
class AliMUONVCalibParam;
class AliMUONVStore;
class AliMUONVStore;
class TMap;

class AliMUONCalibrationData : public TObject
{
public:
  /** Constructor.
    * @param runNumber is used as a key to the CDB
    * @param deferredInitialization if kFALSE, all the calibrations are fetched
    * regardless of whether you'll use them or not.
    */
  AliMUONCalibrationData(Int_t runNumber=-1, Bool_t deferredInitialization=kTRUE);
  virtual ~AliMUONCalibrationData();

  AliMUONVStore* Capacitances() const;

  /// Create a capa store (which must be deleted) from OCDB for the given run
  static AliMUONVStore* CreateCapacitances(Int_t runNumber, Int_t* startOfValidity=0);

  /// Create a gain store (which must be deleted) from OCDB for the given run
  static AliMUONVStore* CreateGains(Int_t runNumber, Int_t* startOfValidity=0);

  /// Create a global trigger mask (which must be deleted) from OCDB for the given run
  static AliMUONGlobalCrateConfig* CreateGlobalTriggerCrateConfig(Int_t runNumber, Int_t* startOfValidity=0);
  
  /// Create a hv map (which must be deleted) from OCDB for the given run
  static TMap* CreateHV(Int_t runNumber, Int_t* startOfValidity=0);

  /// Create a Trigger HV and current  map (which must be deleted) from OCDB for the given run
  static TMap* CreateTriggerDCS(Int_t runNumber, Int_t* startOfValidity=0);

  /// Create a neighbours store (which must be deleted) from OCDB for the given run
  static AliMUONVStore* CreateNeighbours(Int_t runNumber, Int_t* startOfValidity=0);
  
  /// Create a local trigger mask store (which must be deleted) for a given run
  static AliMUONVStore* CreateLocalTriggerBoardMasks(Int_t runNumber, Int_t* startOfValidity=0);

  /// Create an occupancy map store (which must be deleted) from OCDB for the given run
  static AliMUONVStore* CreateOccupancyMap(Int_t runNumber, Int_t* startOfValidity=0);

  /// Create a rejectlist store (which must be deleted) from OCDB for the given run
  static AliMUONRejectList* CreateRejectList(Int_t runNumber, Int_t* startOfValidity=0);

  /// Create a pedestal store (which must be deleted) from OCDB for the given run
  static AliMUONVStore* CreatePedestals(Int_t runNumber, Int_t* startOfValidity=0);

  /// Create a configuration store (which must be deleted) from OCDB for the given run
  static AliMUONVStore* CreateConfig(Int_t runNumber, Int_t* startOfValidity=0);

  /// Create a regional trigger mask store (which must be deleted) for a given run
  static AliMUONRegionalTriggerConfig* CreateRegionalTriggerConfig(Int_t runNumber, Int_t* startOfValidity=0);

  /// Create a trigger Look Up Table (which must be deleted) for a given run
  static AliMUONTriggerLut* CreateTriggerLut(Int_t runNumber, Int_t* startOfValidity=0);
  /// Create a trigger efficiency map (which must be deleted) for a given run
  static AliMUONTriggerEfficiencyCells* CreateTriggerEfficiency(Int_t runNumber, Int_t* startOfValidity=0);
  
  /// Get all the gains
  AliMUONVStore* Gains() const;

  /// Get the configuration for the global trigger board.
  AliMUONGlobalCrateConfig* GlobalTriggerCrateConfig() const;
    
  /// Get the Gain calibration object for channels within (detElemId,manuId).
  AliMUONVCalibParam* Gains(Int_t detElemId, Int_t manuId) const;
    
  /// Get the HV values
  TMap* HV() const;

  /// Get the Trigger HV and current values
  TMap* TriggerDCS() const;
    
  /// Whether this object is valid or not (might be invalid if fetching from CDB failed).
  Bool_t IsValid() const { return fIsValid; }
    
  /// Get the mask for a given local trigger board.
  AliMUONVCalibParam* LocalTriggerBoardMasks(Int_t localBoardNumber) const;
    
  /// Get the neighbours store
  AliMUONVStore* Neighbours() const;
  
  /// Get the pedestal store
  AliMUONVStore* Pedestals() const;

  /// Get the config store
  AliMUONVStore* Config() const;
  
  /// Get the occupancy map store
  AliMUONVStore* OccupancyMap() const;

  /// Get the reject list store
  AliMUONRejectList* RejectList() const;

  /// Get the Pedestal calibration object for channels within (detElemId,manuId).
  AliMUONVCalibParam* Pedestals(Int_t detElemId, Int_t manuId) const;
  
  /// Dump to screen.
  virtual void Print(Option_t* opt="") const;

  /// Get the config for regional trigger.
  AliMUONRegionalTriggerConfig* RegionalTriggerConfig() const;


  /// The runnumber used by this object.
  Int_t RunNumber() const { return fRunNumber; }
  
  /// Get the trigger Look Up Table.
  AliMUONTriggerLut* TriggerLut() const;
  
  /// Get the trigger efficiency map
  AliMUONTriggerEfficiencyCells* TriggerEfficiency() const;
  
  void Reset();

  static TObject* CreateObject(Int_t runNumber, const char* path, Int_t* startOfValidity=0x0);
  
  static void Check(Int_t runNumber);
  
  static void BypassStores(AliMUONVStore* ped, AliMUONVStore* gain);
  
protected:
  /// Not implemented
  AliMUONCalibrationData(const AliMUONCalibrationData& other);
  /// Not implemented
  AliMUONCalibrationData& operator=(const AliMUONCalibrationData& other);
  
private:
  mutable Bool_t fIsValid; ///<  Whether we were able to correctly initialize
  Int_t fRunNumber; ///<  The run number for which we hold calibrations
  mutable AliMUONVStore* fGains; //!< Gains
  mutable AliMUONVStore* fPedestals; //!< Pedestals
  mutable TMap* fHV; //!< HV
  mutable TMap* fTriggerDCS; //!< Trigger HV and Currents
  mutable AliMUONVStore* fLocalTriggerBoardMasks; //!< Local trigger board maska  
  mutable AliMUONRegionalTriggerConfig* fRegionalTriggerConfig; //!< Regional trigger config
  mutable AliMUONGlobalCrateConfig* fGlobalTriggerCrateConfig; //!< Global trigger crate config
  
  mutable AliMUONTriggerLut* fTriggerLut; //!< TRigger LUTs
  mutable AliMUONTriggerEfficiencyCells* fTriggerEfficiency; //!< Trigger efficiency cells
  mutable AliMUONVStore* fCapacitances; //!< Manu capacitances
  mutable AliMUONVStore* fNeighbours; //!< list of neighbours for all channels
  
  mutable AliMUONVStore* fOccupancyMap; //!< occupancy map
  
  mutable AliMUONRejectList* fRejectList; //!< reject list

  static AliMUONVStore* fgBypassPedestals;
  static AliMUONVStore* fgBypassGains;
  
  mutable AliMUONVStore* fConfig; //!< configuration of the tracker
  
  ClassDef(AliMUONCalibrationData,13) // Storage for all MUON calibration data.
};

#endif
