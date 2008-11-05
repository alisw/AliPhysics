//-*- Mode: C++ -*-
// @(#) $Id: AliHLTCompPreprocessor.h 23318 2008-01-14 12:43:28Z hristov $

#ifndef ALIHLTCOMPPREPROCESSOR_H
#define ALIHLTCOMPPREPROCESSOR_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               */

/**
 * @file   AliHLTCompPreprocessor.h
 * @author Jenny Wagner, Matthias Richter
 * @brief  HLT Preprocessor plugin for the AliHLTComp library
 */

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "AliHLTModulePreprocessor.h"
#include "AliHLTPreprocessor.h"

/**
 * @class AliHLTCompPreprocessor
 * HLT preprocessor for the libAliHLTComp module.
 *
 * @author Jenny Wagner, Matthias Richter
 *
 * @date 2008-01-22
 */
class AliHLTCompPreprocessor : public AliHLTModulePreprocessor
{
 public:
	
  /** Standard Constructor */
  AliHLTCompPreprocessor();

  /** Destructor */
  ~AliHLTCompPreprocessor();

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

  /** Define name of huffman tables stored at FXS */
  static const char* fgkHuffmanFileId;			// see above

  /** Define module id */
  const char* GetModuleID() {return "AliHLTCompPreprocessor";};

  /** Define bit mask of the active detectors needed by this preprocessor module */
  Int_t GetModuleNumber() {
    Int_t modulenumber = 0;
    modulenumber = AliHLTModulePreprocessor::DetectorBitMask("TPC") | AliHLTModulePreprocessor::DetectorBitMask("PHOS");
    return modulenumber;
  };

 protected:

 private:
  /** copy constructor prohibited */
  AliHLTCompPreprocessor(const AliHLTCompPreprocessor& preproc);
  /** assignment operator prohibited */
  AliHLTCompPreprocessor& operator=(const AliHLTCompPreprocessor& rhs);

  /**
   * Function fetch and prepare the Huffman tables from the HLT FXS
   *
   * @return 0 in case of success, else an error code
   */
  UInt_t GetHuffmanTables();

  /** mark if TPC was active (1 = active) */
  Bool_t fTPCactive;  // mark if TPC was active
  /** mark if PHOS was active (1 = active) */
  Bool_t fPHOSactive; // makr if PHOS was active
		
  ClassDef(AliHLTCompPreprocessor, 1);
};
#endif
