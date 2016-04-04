/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id: AliMUONCalibrationData.h 59457 2012-11-06 12:36:48Z laphecet $

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

  /// Create a global trigger mask (which must be deleted) from OCDB for the given run
  static AliMUONGlobalCrateConfig* CreateGlobalTriggerCrateConfig(Int_t runNumber, Int_t* startOfValidity=0);

  /// Create a hv map (which must be deleted) from OCDB for the given run
  static TMap* CreateHV(Int_t runNumber, Int_t* startOfValidity=0, Bool_t patched=kTRUE, TList* messages=0x0, Bool_t dryRun=kFALSE);

  /// Create a Trigger HV and current  map (which must be deleted) from OCDB for the given run
  static TMap* CreateTriggerDCS(Int_t runNumber, Int_t* startOfValidity=0);

  /// Create a MCH LV map (which must be deleted) from OCDB for the given run
  static TMap* CreateLV(Int_t runNumber, Int_t* startOfValidity=0);

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

  /// Get the configuration for the global trigger board.
  AliMUONGlobalCrateConfig* GlobalTriggerCrateConfig() const;

  /// Get the HV values. Use patched=kFALSE to get unprocessed (i.e. "raw") values as they are in the OCDB
  TMap* HV(Bool_t patched=kTRUE) const;

  /// Get the Trigger HV and current values
  TMap* TriggerDCS() const;

  /// Get the MCH LV
  TMap* LV() const;

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

  static void BypassStores(AliMUONVStore* ped);

  static void PatchSt1DCSAliases(TMap& hvMap);

  static Bool_t PatchHVValues(TObjArray& values, TString* msg=0x0, Bool_t dryRun=kFALSE);

  static UInt_t PatchHVDCSAliasesSt1WasAppliedMask() { return fgkPatchHVDCSAliasesSt1WasAppliedMask; }

  static UInt_t PatchHVAllWasAppliedMask() { return fgkPatchHVAllWasAppliedMask; }

protected:
  /// Not implemented
  AliMUONCalibrationData(const AliMUONCalibrationData& other);
  /// Not implemented
  AliMUONCalibrationData& operator=(const AliMUONCalibrationData& other);

  static Bool_t CheckHVGroup(TObjArray& values, Int_t first, Int_t last, Double_t& value,
                             Int_t& slope, TString* msg);

  static void AddToMap(const TMap& sourceMap,
                       TMap& destMap,
                       const TString& key,
                       const char* source,
                       const char* dest);

private:
  mutable Bool_t fIsValid; ///<  Whether we were able to correctly initialize
  Int_t fRunNumber; ///<  The run number for which we hold calibrations
  mutable AliMUONVStore* fPedestals; //!<! Pedestals
  mutable TMap* fHV; //!<! HV
  mutable TMap* fTriggerDCS; //!<! Trigger HV and Currents
  mutable AliMUONVStore* fLocalTriggerBoardMasks; //!<! Local trigger board maska
  mutable AliMUONRegionalTriggerConfig* fRegionalTriggerConfig; //!<! Regional trigger config
  mutable AliMUONGlobalCrateConfig* fGlobalTriggerCrateConfig; //!<! Global trigger crate config

  mutable AliMUONTriggerLut* fTriggerLut; //!<! TRigger LUTs
  mutable AliMUONTriggerEfficiencyCells* fTriggerEfficiency; //!<! Trigger efficiency cells
  mutable AliMUONVStore* fNeighbours; //!<! list of neighbours for all channels

  mutable AliMUONVStore* fOccupancyMap; //!<! occupancy map

  mutable AliMUONRejectList* fRejectList; //!<! reject list

  static AliMUONVStore* fgBypassPedestals;

  mutable AliMUONVStore* fConfig; //!<! configuration of the tracker

  mutable TMap* fLV; //!<! MCH LV
  
  static UInt_t fgkPatchHVDCSAliasesSt1WasAppliedMask; //!<! mask to indicate that the DCS alias naming is not messed up in St1
  static UInt_t fgkPatchHVAllWasAppliedMask; //!<! mask to indicate that HV values were massaged already

  ClassDef(AliMUONCalibrationData,17) // Storage for all MUON calibration data.
};

#endif
