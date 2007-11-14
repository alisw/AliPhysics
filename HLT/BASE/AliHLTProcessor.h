//-*- Mode: C++ -*-
// @(#) $Id$

#ifndef ALIHLTPROCESSOR_H
#define ALIHLTPROCESSOR_H
/* This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTProcessor.h
    @author Matthias Richter, Timm Steinbeck
    @date   
    @brief  Base class declaration for HLT analysis components. */

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt   

#include "AliHLTComponent.h"

/**
 * @class AliHLTProcessor
 * Base class of HLT data analysis components.
 * The class provides a common interface for the implementation of HLT data
 * analysis components. The child class must implement the functions:
 * - @ref DoInit (optional)
 * - @ref DoDeinit (optional)
 * - @ref DoEvent
 * - @ref GetComponentID
 * - @ref GetInputDataTypes
 * - @ref GetOutputDataType
 * - @ref GetOutputDataSize
 * - @ref Spawn
 *
 * @ingroup alihlt_component
 */
class AliHLTProcessor : public AliHLTComponent {
 public:
  /** standard constructor */
  AliHLTProcessor();
  /** standard destructor */
  virtual ~AliHLTProcessor();

  /* depricated */
  int Init( AliHLTComponentEnvironment* environ, void* environParam, int argc, const char** argv );
  /* depricated */
  int Deinit();

  /**
   * Event processing function.
   * The method is called by the framework to process one event. After 
   * preparation of data structures. The call is redirected to DoEvent.
   * @return neg. error code if failed 
   */
  int DoProcessing( const AliHLTComponentEventData& evtData,
		    const AliHLTComponentBlockData* blocks, 
		    AliHLTComponentTriggerData& trigData,
		    AliHLTUInt8_t* outputPtr, 
		    AliHLTUInt32_t& size,
		    vector<AliHLTComponentBlockData>& outputBlocks,
		    AliHLTComponentEventDoneData*& edd );

  // Information member functions for registration.

  /**
   * Return @ref AliHLTComponent::kProcessor type as component type.
   * @return component type id
   */
  TComponentType GetComponentType() { return AliHLTComponent::kProcessor;}

 protected:
  /**
   * The low-level data processing method for the component.
   * This is the custom processing method and can be overloaded by 
   * the component.
   * @param evtData       event data structure
   * @param blocks        input data block descriptors
   * @param trigData	  trigger data structure
   * @param outputPtr	  pointer to target buffer
   * @param size	  <i>input</i>: size of target buffer
   *            	  <i>output</i>:size of produced data
   * @param outputBlocks  list to receive output block descriptors
   * @return neg. error code if failed                                <br>
   *         -ENOSPC      output buffer too small
   */
  virtual int DoEvent( const AliHLTComponentEventData& evtData,
		       const AliHLTComponentBlockData* blocks, 
		       AliHLTComponentTriggerData& trigData,
		       AliHLTUInt8_t* outputPtr, 
		       AliHLTUInt32_t& size,
		       vector<AliHLTComponentBlockData>& outputBlocks );

  /**
   * The high-level data processing method.
   * This is the default processing method; the method is called
   * if no low level @ref DoEvent method is overloaded by the component.
   * @param evtData       event data structure
   * @param trigData	  trigger data structure
   * @return neg. error code if failed
   */
  virtual int DoEvent( const AliHLTComponentEventData& evtData, AliHLTComponentTriggerData& trigData);

  ClassDef(AliHLTProcessor, 1)
};
#endif
