// @(#) $Id$

#ifndef ALIHLTPROCESSOR_H
#define ALIHLTPROCESSOR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTProcessor.h
    @author Matthias Richter, Timm Steinbeck
    @date   
    @brief  Base class declaration for HLT analysis components. */

#include "AliHLTComponent.h"

/**
 * @class AliHLTProcessor
 * Base class of HLT data analysis components.
 * The class provides a common interface for the implementation of HLT data
 * analysis components. The child class must implement the functions:
 * - DoInit (optional)
 * - DoDeinit (optional)
 * - DoEvent
 * - GetComponentID
 * - GetInputDataTypes
 * - GetOutputDataType
 * - GetOutputDataSize
 * - Spawn
 *
 * @ingroup AliHLTbase
 */
class AliHLTProcessor : public AliHLTComponent {
 public:
  /** standard constructor */
  AliHLTProcessor();
  /** standard destructor */
  virtual ~AliHLTProcessor();

  /* depricated */
  int Init( AliHLTComponentEnvironment* environ, void* environ_param, int argc, const char** argv );
  /* depricated */
  int Deinit();

  /**
   * Event processing function.
   * The method is called by the framework to process one event. After 
   * preparation of data structures. The call is redirected to DoEvent.
   * @return neg. error code if failed 
   */
  int ProcessEvent( const AliHLTComponent_EventData& evtData,
		    const AliHLTComponent_BlockData* blocks, 
		    AliHLTComponent_TriggerData& trigData,
		    AliHLTUInt8_t* outputPtr, 
		    AliHLTUInt32_t& size,
		    AliHLTUInt32_t& outputBlockCnt, 
		    AliHLTComponent_BlockData*& outputBlocks,
		    AliHLTComponent_EventDoneData*& edd );

  // Information member functions for registration.

  /**
   * Return @ref AliHLTComponent::kProcessor type as component type.
   * @return component type id
   */
  TComponentType GetComponentType() { return AliHLTComponent::kProcessor;}

 private:
  /**
   * Data processing method for the component.
   * @param evtData       event data structure
   * @param blocks        input data block descriptors
   * @param trigData	  trigger data structure
   * @param outputPtr	  pointer to target buffer
   * @param size	  <i>input</i>: size of target buffer
   *            	  <i>output</i>:size of produced data
   * @param outputBlocks  list to receive output block descriptors
   */
  virtual int DoEvent( const AliHLTComponent_EventData& evtData,
		       const AliHLTComponent_BlockData* blocks, 
		       AliHLTComponent_TriggerData& trigData,
		       AliHLTUInt8_t* outputPtr, 
		       AliHLTUInt32_t& size,
		       vector<AliHLTComponent_BlockData>& outputBlocks ) = 0;

  ClassDef(AliHLTProcessor, 0)
};
#endif
