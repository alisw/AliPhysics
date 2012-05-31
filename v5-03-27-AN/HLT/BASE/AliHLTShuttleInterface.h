//-*- Mode: C++ -*-
// @(#) $Id$

#ifndef ALIHLTSHUTTLEINTERFACE_H
#define ALIHLTSHUTTLEINTERFACE_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/**
 * @file   AliHLTShuttleInterface.h
 * @author Matthias Richter
 * @date   2008-01-22
 * @brief  Pure virtual interface to the HLT shuttle methods
 */

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

/**
 * @class AliHLTShuttleInterface
 * This class implements the redirection of the shuttle access methods for
 * AliHLTModulePreprocessor classes. The AliHLTShuttleInterface has been
 * declared pure virtual to avoid dependencies between the libHLTshuttle and
 * the component libraries. It implements the same interface to the shuttle
 * as the AliPreprocessor.
 *
 * The AliHLTPreprocessor initializes each AliHLTModulePreprocessor with this
 * interface. From the interface methods of AliHLTModulePreprocessor classes,
 * the call is redirected via AliHLTShuttleInterface to the AliHLTPreprocessor,
 * which makes the methods publicly available AliPreprocessor.
 *
 * @author Matthias Richter
 */
class AliHLTShuttleInterface
{
public:
  /** Constructor*/
  AliHLTShuttleInterface();
  /** Destructor */
  virtual ~AliHLTShuttleInterface();

  /** Get the run no */
  virtual Int_t PreprocessorGetRun() = 0;

  /** Get the start time */
  virtual UInt_t PreprocessorGetStartTime() = 0;

  /** Get the end time */
  virtual UInt_t PreprocessorGetEndTime() = 0;

  // the AliPreprocessor interface, all functions redirected via the
  // AliHLTPreprocessor
  virtual Bool_t PreprocessorStore(const char* pathLevel2, const char* pathLevel3, TObject* object,
	       AliCDBMetaData* metaData, Int_t validityStart = 0, Bool_t validityInfinite = kFALSE) = 0;
  virtual Bool_t PreprocessorStoreReferenceData(const char* pathLevel2, const char* pathLevel3, TObject* object,
			    AliCDBMetaData* metaData) = 0;
  virtual Bool_t PreprocessorStoreReferenceFile(const char* localFile, const char* gridFileName) = 0;

  virtual Bool_t PreprocessorStoreRunMetadataFile(const char* localFile, const char* gridFileName) = 0;

  virtual const char* PreprocessorGetFile(Int_t system, const char* id, const char* source) = 0;

  virtual TList* PreprocessorGetFileSources(Int_t system, const char* id = 0) = 0;

  virtual TList* PreprocessorGetFileIDs(Int_t system, const char* source) = 0;

  virtual const char* PreprocessorGetRunParameter(const char* param) = 0;

  virtual AliCDBEntry* PreprocessorGetFromOCDB(const char* pathLevel2, const char* pathLevel3) = 0;

  virtual const char* PreprocessorGetRunType() = 0;

  virtual void PreprocessorLog(const char* message) = 0;

protected:

private:
  /** copy constructor prohibited */
  AliHLTShuttleInterface(const AliHLTShuttleInterface&);
  /** assignment operator prohibited */
  AliHLTShuttleInterface& operator=(const AliHLTShuttleInterface&);

  ClassDef(AliHLTShuttleInterface, 0);
};
#endif
