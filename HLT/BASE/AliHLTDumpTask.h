//-*- Mode: C++ -*-
// $Id$
#ifndef ALIHLTDUMPTASK_H
#define ALIHLTDUMPTASK_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/** @file   AliHLTDumpTask.h
    @author Matthias Richter
    @date   
    @brief  Base class for data sinks with direct buffer access.
*/

#include "AliHLTOUT.h"
#include "AliHLTTask.h"

/**
 * @class AliHLTDumpTask
 * A base class for data sinks without a processing component. Tasks of this
 * type can be added at the end of the processing chain and allow direct
 * access to the data blocks produced by the last components of the chain.
 * 
 * The constructor takes the chains as a blank separated list of chain ids.
 *
 * @ingroup alihlt_system
 */
class AliHLTDumpTask : public AliHLTTask {
 public:
  /** constructor */
  AliHLTDumpTask(const char* chains=NULL);
  /** standard destructor */
  virtual ~AliHLTDumpTask();

  /**
   * Set the parent reconstruction chains this task is going to subscribe to.
   * @param chains    blank separated list of chain ids
   */
  int SetChains(const char* chains);

  /**
   * Get string of the source chains.
   */
  const char* GetSourceChains() const;

  /**
   * Get the list of data blocks.
   * Subscribe to the publishers if not yet done and return pointer
   * to block list.
   */
  const AliHLTComponentBlockDataList& GetDataBlocks();

  /**
   * Explicitly release data blocks.
   */
  int ReleaseDataBlocks();

 protected:

 private:
  /** copy constructor prohibited */
  AliHLTDumpTask(const AliHLTDumpTask&);
  /** assignment operator prohibited */
  AliHLTDumpTask& operator=(const AliHLTDumpTask&);

  /**
   * Custom initialization for child tasks.
   * Create and init the dummy task.
   */
  int CustomInit(AliHLTComponentHandler* pCH);

  /**
   * Custom clean up for child tasks.
   */
  int CustomCleanup();

  /** a dummy task to pretend existence of a consumer */
  AliHLTTask* fpDummyTask; //!transient

  /** the configuration for the dummy task */
  AliHLTConfiguration* fpDummyConfiguration; //!transient

  /** list of block descriptors of the output */
  AliHLTComponentBlockDataList fBlocks;

  ClassDef(AliHLTDumpTask, 0)
};
#endif
