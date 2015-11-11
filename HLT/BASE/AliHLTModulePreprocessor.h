//-*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTMODULEPREPROCESSOR_H
#define ALIHLTMODULEPREPROCESSOR_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

// @file   AliHLTModulePreprocessor.h
// @author Matthias Richter
// @date   2008-01-22
// @brief  Base class for HLT module preprocessors
// 

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "TObject.h"

class AliHLTPreprocessor;
class TMap;
class AliCDBMetaData;
class AliCDBEntry;
class AliHLTShuttleInterface;

/**
 * @class AliHLTModulePreprocessor
 * Implementation of the HLT version for the Shuttle Preprocessor.
 * This class implements the same interface as the AliPreprocessor and
 * allows multiple preprocessors for the HLT. All methods are redirected
 * to the AliHLTPreprocessor, which acts as the HLT preprocessor
 * to the outside, and a container for all module preprocessors to the
 * inside.
 *
 * @author Matthias Richter
 */
class AliHLTModulePreprocessor : public TObject
{
public:
  /** Constructor*/
  AliHLTModulePreprocessor();
  /** Destructor */
  virtual ~AliHLTModulePreprocessor();

  /**
   * Set the container class which is the gateway to the shuttle.
   */
  void SetShuttleInterface(AliHLTShuttleInterface* pInterface);

  /**
   * Initialize the Preprocessor.
   *
   * @param run run number
   * @param startTime start time of data
   * @param endTime end time of data
   */
  virtual void Initialize(Int_t run, UInt_t startTime, UInt_t endTime) = 0;

  /**
   * Function to process data. Inside the preparation and storing to OCDB
   * should be handled.
   *
   * @param dcsAliasMap the map containing aliases and corresponding DCS
   * 			values and timestamps
   *
   * @return 0 on success; error code otherwise
   */
  virtual UInt_t Process(TMap* dcsAliasMap) = 0;

  /** Get the run no */
  Int_t GetRun() const;

  /** Get the start time */
  UInt_t GetStartTime() const;

  /** Get the end time */
  UInt_t GetEndTime() const;

 /**
   * Get the id of the module
   * Each module preprocessor is identified by a unique id.
   * The function is pure virtual and must be implemented by the child class.
   * @return module id (string)
   */
  virtual const char* GetModuleID() {return ClassName();}

 /**
   * Get all detectors the module process
   * Detectors are determined in the bit mask as provided by "detectorMask"
   * implement: GetModuleNumber() {return AliHLTModulePreprocessor::DetetorID("detectorname");}
   * The function is pure virtual and must be implemented by the child class.
   * @return bit mask for active detectors (Int_t)
   */
  virtual Int_t GetModuleNumber() = 0;

  /** Get detector bit mask (bit mask comparable to the detectorMask but for the respective detector only)
   *  out of detector name
   *@param detectorName const char* of the detector
   *@return Int_t for the detector bit mask
   */
  Int_t DetectorBitMask(const char *detectorName);

 /** Get the information whether detector was active (return value: 1 = active, 0 = inactive) 
      @param detectorbitmask bitmask of the detector to be checked (s. DetectorBitMask)
      @return 1 if active 0 if inactive
 */
  Bool_t GetDetectorStatus(Int_t detectorbitmask);

  /** Get detector name out of detector id (bit position in detectorMask of the resp. detector)
   *@param detectorID integer ID of the detector
   *@return const char* name of the detector corresponding to the ID
   */
  const char *DetectorName(Int_t detectorID);

  /// extract object from the alias map
  /// return the last object from the value set
  TObject* GetFromMap(TMap* dcsAliasMap, const char* stringId) const;

  /** number of detectors */
  //static const Int_t kNDetectors;    // Number of detectors

protected:
  // the AliPreprocessor interface, all functions redirected via the
  // AliHLTShuttleInterface to the AliHLTPreprocessor
  Bool_t Store(const char* pathLevel2, const char* pathLevel3, TObject* object,
	       AliCDBMetaData* metaData, Int_t validityStart = 0, Bool_t validityInfinite = kFALSE);
  Bool_t StoreReferenceData(const char* pathLevel2, const char* pathLevel3, TObject* object,
			    AliCDBMetaData* metaData);
  Bool_t StoreReferenceFile(const char* localFile, const char* gridFileName);

  Bool_t StoreRunMetadataFile(const char* localFile, const char* gridFileName);

  const char* GetFile(Int_t system, const char* id, const char* source);

  TList* GetFileSources(Int_t system, const char* id = 0);

  TList* GetFileIDs(Int_t system, const char* source);

  const char* GetRunParameter(const char* param);

  AliCDBEntry* GetFromOCDB(const char* pathLevel2, const char* pathLevel3);

  const char* GetRunType();

  void Log(const char* message);

private:
  /** copy constructor prohibited */
  AliHLTModulePreprocessor(const AliHLTModulePreprocessor& preproc);
  /** assignment operator prohibited */
  AliHLTModulePreprocessor& operator=(const AliHLTModulePreprocessor& rhs);

  /** the interface class which is the gateway to the shuttle */
  AliHLTShuttleInterface* fpInterface;                             //! transient

 /** determine which which detectors were active */
  Int_t fActiveDetectors; // bit array of active detectors

  /** array of detector names */
  //static const char *fgkDetectorName[]; // Detector names

  ClassDef(AliHLTModulePreprocessor, 1);
};
#endif
