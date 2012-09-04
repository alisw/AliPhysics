//-*- Mode: C++ -*-
// $Id: AliHLTPreprocessor.h 23318 2008-01-14 12:43:28Z hristov $

#ifndef ALIHLTPREPROCESSOR_H
#define ALIHLTPREPROCESSOR_H
//* This file is property of and copyright by the                          * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/// @file   AliHLTPreprocessor.h
/// @author Matthias Richter
/// @date   2008-01-22
/// @brief  Container for HLT module preprocessors, acts to the outside as
///         HLT preprocessor used by the Offline Shuttle 
/// 

#include "TList.h"
#include "AliPreprocessor.h"
#include "AliHLTShuttleInterface.h"

/**
 * @class AliHLTPreprocessor
 * Implementation of the HLT version for the Shuttle Preprocessor.
 * Since HLT requires a more modular concept of the pre-processors, this
 * class acts as HLT pre-processor to the outside and container class for
 * the specific HLT module pre-processors to the inside.
 *
 * The base class for HLT module preprocessors is provided by the
 * AliHLTModulePreprocessor class, which implements the same interface as
 * the AliPreprocessor.
 *
 * The main purpose of the container class is to loop over all module
 * preprocessors and to make the AliPreprocessor interface methods
 * publicly available.
 */
class AliHLTPreprocessor : public AliPreprocessor , public AliHLTShuttleInterface
{
 public:
  /**
   * Constructor for AliHLTPreprocessor
   *
   * @param shuttle pointer to the hosting shuttle
   */
  AliHLTPreprocessor(AliShuttleInterface* shuttle);
  /** Destructor */
  virtual ~AliHLTPreprocessor();

  /**
   * Initialize the Preprocessor.
   *
   * @param run run number
   * @param startTime start time of data
   * @param endTime end time of data
   */
  virtual void Initialize(Int_t run, UInt_t startTime, UInt_t endTime);

  /**
   * Function to process data. Inside the preparation and storing to OCDB
   * should be handled.
   *
   * @param dcsAliasMap the map containing aliases and corresponding DCS
   * 			values and timestamps
   *
   * @return 0 on success; a value greater than 0 refers to an error
   */
  virtual UInt_t Process(TMap* dcsAliasMap);

  /**
   * Indicates if DCS data can be processed.
   *
   * @return true if DCS data can be processed, else false. 
   */
  virtual Bool_t ProcessDCS() {return kFALSE;}

  /** Define for name of the HLT Preproc */
  static const char* fgkHLTPreproc; 			// see above

  /** Get the run no which has been previously initialized */
  Int_t PreprocessorGetRun() {return fRun;}

  /** Get the start time no which has been previously initialized */
  UInt_t PreprocessorGetStartTime() {return fStartTime;}

  /** Get the end time no which has been previously initialized */
  UInt_t PreprocessorGetEndTime() {return fEndTime;}

  // AliPreprocessor methods made publicly available
  // the subsequent functions map the AliPreprocessor interface functions in order
  // to be used by the module proprocessors.
  Bool_t PreprocessorStore(const char* pathLevel2, const char* pathLevel3, TObject* object,
	       AliCDBMetaData* metaData, Int_t validityStart = 0, Bool_t validityInfinite = kFALSE) {
    return AliPreprocessor::Store(pathLevel2, pathLevel3, object, metaData, validityStart, validityInfinite);
  }

  Bool_t PreprocessorStoreReferenceData(const char* pathLevel2, const char* pathLevel3, TObject* object,
			    AliCDBMetaData* metaData) {
    return AliPreprocessor::StoreReferenceData(pathLevel2, pathLevel3, object, metaData);
  }

  Bool_t PreprocessorStoreReferenceFile(const char* localFile, const char* gridFileName) {
    return AliPreprocessor::StoreReferenceFile(localFile, gridFileName);
  }

  Bool_t PreprocessorStoreRunMetadataFile(const char* localFile, const char* gridFileName) {
    return AliPreprocessor::StoreRunMetadataFile(localFile, gridFileName);
  }
    
  const char* PreprocessorGetFile(Int_t system, const char* id, const char* source) {
    return AliPreprocessor::GetFile(system, id, source);
  }

  TList* PreprocessorGetFileSources(Int_t system, const char* id = 0) {
    return AliPreprocessor::GetFileSources(system, id);
  }

  TList* PreprocessorGetFileIDs(Int_t system, const char* source) {
    return AliPreprocessor::GetFileIDs(system, source);
  }

  const char* PreprocessorGetRunParameter(const char* param) {
    return AliPreprocessor::GetRunParameter(param);
  }

  AliCDBEntry* PreprocessorGetFromOCDB(const char* pathLevel2, const char* pathLevel3) {
    return AliPreprocessor::GetFromOCDB(pathLevel2, pathLevel3);
  }

  const char* PreprocessorGetRunType() {
    return AliPreprocessor::GetRunType();
  }

  void PreprocessorLog(const char* message) {
    AliPreprocessor::Log(message);
  }

 protected:

 private:
  /** copy constructor prohibited */
  AliHLTPreprocessor(const AliHLTPreprocessor& preproc);
  /** assignment operator prohibited */
  AliHLTPreprocessor& operator=(const AliHLTPreprocessor& rhs);

  /** list of HLT module processors */
  TList fProcessors;                                               //!transient

  /** determine which which detectors were active */
  Int_t fActiveDetectors; // bit array of active detectors
  
  /** array of default libraries */
  static const char* fgkHLTDefaultShuttleLibs[];                   //!transient

  ClassDef(AliHLTPreprocessor, 1);
};
#endif


