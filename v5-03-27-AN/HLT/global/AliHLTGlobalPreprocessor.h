//-*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTGLOBALPREPROCESSOR_H
#define ALIHLTGLOBALPREPROCESSOR_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               */

//  @file   AliHLTGlobalPreprocessor.h
//  @author Matthias Richter
//  @date   2010-08-20
//  @brief  HLT Preprocessor plugin for global HLT
// 

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "AliHLTModulePreprocessor.h"

/**
 * @class AliHLTGlobalPreprocessor
 * HLT preprocessor for global HLT objects
 *
 * <h2>Produced OCDB objects:</h2>
 * - HLT/Calib/Streamerinfo <br>
 *   The streamer info object is produced by the ROOTSchemaEvolutionComponent
 *   See ProcessStreamerInfo() for details.
 *
 * @author Matthias Richter
 */
class AliHLTGlobalPreprocessor : public AliHLTModulePreprocessor
{
 public:
	
  /** Standard Constructor */
  AliHLTGlobalPreprocessor();

  /** Destructor */
  ~AliHLTGlobalPreprocessor();

  /**
   * Initialize the Preprocessor.
   *
   * @param run run number
   * @param startTime start time of data
   * @param endTime end time of data
   */
  void Initialize(Int_t run, UInt_t startTime, UInt_t endTime);

  /**
   * Function to process data. Inside the preparation and storing to OCDB
   * should be handled.
   *
   * @param dcsAliasMap the map containing aliases and corresponding DCS
   * 			values and timestamps
   *
   * @return 0 on success; error code otherwise
   */
  UInt_t Process(TMap* dcsAliasMap);

  /** Define bit mask of the active detectors needed by this preprocessor module */
  Int_t GetModuleNumber();

  /// DCS alias 'StreamerInfo' -> Calib/StreamerInfo
  static const char* fgkStreamerInfoAlias;
  static const char* fgkStreamerInfoName;
  static const char* fgkStreamerInfoType;

 protected:

 private:
  /** copy constructor prohibited */
  AliHLTGlobalPreprocessor(const AliHLTGlobalPreprocessor& preproc);
  /** assignment operator prohibited */
  AliHLTGlobalPreprocessor& operator=(const AliHLTGlobalPreprocessor& rhs);

  int ProcessStreamerInfo(TObject* object);

  ClassDef(AliHLTGlobalPreprocessor, 0);
};
#endif
