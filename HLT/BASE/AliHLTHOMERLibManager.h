//-*- Mode: C++ -*-
// @(#) $Id$

#ifndef ALIHLTHOMERLIBMANAGER_H
#define ALIHLTHOMERLIBMANAGER_H
/* This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTHOMERLibManager.h
    @author Matthias Richter
    @date   
    @brief  dynamic HLT HOMER reader/writer generation and destruction

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
                                                                          */

#include "AliHLTLogging.h"

class AliHLTHOMERReader;
class AliHLTHOMERWriter;

/**
 * @class AliHLTHOMERLibManager
 * Handler of HLTOUT data for buffer input.
 */
class AliHLTHOMERLibManager : public AliHLTLogging {
 public:
  /** standard constructor */
  AliHLTHOMERLibManager();
  /** destructor */
  virtual ~AliHLTHOMERLibManager();

  /**
   * Open a HOMER reader.
   * Load HOMER library dynamically and create object working on the provided
   * buffer.
   */
  AliHLTHOMERReader* OpenReader(const AliHLTUInt8_t* pBuffer, int size);

  /**
   * Delete a HOMER reader.
   * Clean-up of the object is done inside the HOMER library.
   */
  int DeleteReader(AliHLTHOMERReader* pReader);

  /**
   * Open a HOMER writer.
   * Load HOMER library dynamically and create object working on the provided
   * buffer.
   */
  AliHLTHOMERWriter* OpenWriter();

  /**
   * Delete a HOMER writer.
   * Clean-up of the object is done inside the HOMER library.
   */
  int DeleteWriter(AliHLTHOMERWriter* pWriter);

 protected:

 private:
  /** copy constructor prohibited */
  AliHLTHOMERLibManager(const AliHLTHOMERLibManager&);
  /** assignment operator prohibited */
  AliHLTHOMERLibManager& operator=(const AliHLTHOMERLibManager&);

  /**
   * Load the HOMER library.
   */
  int LoadHOMERLibrary();

  /** status of the loading of the HOMER library */
  int fLibraryStatus; //!transient

  /** entry in the HOMER library */
  void* fFctCreateReaderFromBuffer; //!transient

  /** entry in the HOMER library */
  void* fFctDeleteReader; //!transient

  /** entry in the HOMER library */
  void* fFctCreateWriter; //!transient

  /** entry in the HOMER library */
  void* fFctDeleteWriter; //!transient

  ClassDef(AliHLTHOMERLibManager, 0)
};
#endif
