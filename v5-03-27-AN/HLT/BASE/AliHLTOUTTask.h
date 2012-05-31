//-*- Mode: C++ -*-
// $Id$
#ifndef ALIHLTOUTTASK_H
#define ALIHLTOUTTASK_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/** @file   AliHLTOUTTask.h
    @author Matthias Richter
    @date   
    @brief  A special HLTOUT sibling working as a data sink in chains
*/

#include "AliHLTOUT.h"
#include "AliHLTDumpTask.h"

/**
 * @class AliHLTOUTTask
 * A special HLTOUT sibling implementing AliHLTDataSink functionality in
 * order to be run at the end of a reconstruction chain and generation of
 * an HLTOUT sub-collection.
 * 
 * The constructor takes the chains as a blank separated list of chain ids.
 *
 * @ingroup alihlt_system
 */
class AliHLTOUTTask : public AliHLTOUT, public AliHLTDumpTask {
 public:
  /** constructor */
  AliHLTOUTTask(const char* chains);
  /** standard destructor */
  virtual ~AliHLTOUTTask();

 protected:

 private:
  /** standard constructor prohibited */
  AliHLTOUTTask();
  /** copy constructor prohibited */
  AliHLTOUTTask(const AliHLTOUTTask&);
  /** assignment operator prohibited */
  AliHLTOUTTask& operator=(const AliHLTOUTTask&);

  /**
   * Generate the index of the HLTOUT data.
   * Must be implemented by the child classes.
   */
  int GenerateIndex();

  /**
   * Cleanup and reset the data input.
   */
  int ResetInput();

  /**
   * Get the data buffer
   * @param [in] index    index of the block
   * @param [out] pBuffer buffer of the selected data block
   * @param [out] size    size of the selected data block
   */
  int GetDataBuffer(AliHLTUInt32_t index, const AliHLTUInt8_t* &pBuffer, 
		    AliHLTUInt32_t& size);

  /**
   * Check byte order of data block
   */
  AliHLTOUTByteOrder CheckBlockByteOrder(AliHLTUInt32_t index);

  /**
   * Check alignment of data block
   */
  int CheckBlockAlignment(AliHLTUInt32_t index, AliHLTOUT::AliHLTOUTDataType type);

  ClassDef(AliHLTOUTTask, 1)
};
#endif
