//-*- Mode: C++ -*-
// @(#) $Id: AliHLTSamplePreprocessor.h 23318 2008-01-14 12:43:28Z hristov $

#ifndef ALIHLTSAMPLEPREPROCESSOR_H
#define ALIHLTSAMPLEPREPROCESSOR_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               */

/**
 * @file   AliHLTSamplePreprocessor.h
 * @author Kenneth Aamodt, Sebastian Bablok
 * @date   2007-12-06
 * @brief  HLT Preprocessor plugin for the AliHLTSample library
 */

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "AliHLTModulePreprocessor.h"

/**
 * @class AliHLTSamplePreprocessor
 * A sample HLT preprocessor.
 *
 * @date 2008-01-22
 */
class AliHLTSamplePreprocessor : public AliHLTModulePreprocessor
{
 public:
	
  /** Constructor */
  AliHLTSamplePreprocessor();

  /** Destructor */
  ~AliHLTSamplePreprocessor();

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

  /** Define for Temperature Histogram filename */
  static const char* fgkTempHistoFileName; 		// see above

  /** Define module id */
  const char* GetModuleID() {return "AliHLTSamplePreprocessor";};

  Int_t GetModuleNumber() {return AliHLTModulePreprocessor::DetectorBitMask("TPC");};

 protected:

 private:
  /** copy constructor prohibited */
  AliHLTSamplePreprocessor(const AliHLTSamplePreprocessor& preproc);
  /** assignment operator prohibited */
  AliHLTSamplePreprocessor& operator=(const AliHLTSamplePreprocessor& rhs);

  /**
   * Function fetch and prepare a dummy temperature histogram from the 
   * HLT FXS.
   *
   * @return 0 in case of success, else an error code
   */
  UInt_t GetTempHisto();
		
  ClassDef(AliHLTSamplePreprocessor, 0);
};
#endif
