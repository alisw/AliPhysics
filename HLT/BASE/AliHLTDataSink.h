//-*- Mode: C++ -*-
// @(#) $Id$

#ifndef ALIHLTDATASINK_H
#define ALIHLTDATASINK_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTDataSink.h
    @author Matthias Richter
    @date   
    @brief  Base class declaration for HLT data sink components.
    @note   The class is used in Offline (AliRoot) context
*/

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt   

#include "AliHLTComponent.h"

/**
 * @class AliHLTDataSink
 * Base class of HLT data sink components.
 * The class provides a common interface for the implementation of HLT data
 * sink components. The child class must implement the functions:
 * - @ref DoInit (optional)
 * - @ref DoDeinit (optional)
 * - @ref DumpEvent
 * - @ref GetComponentID
 * - @ref GetInputDataTypes
 * - @ref Spawn
 *
 * @ingroup alihlt_component
 */
class AliHLTDataSink : public AliHLTComponent {
 public:
  /** standard constructor */
  AliHLTDataSink();
  /** standard destructor */
  virtual ~AliHLTDataSink();

  /**
   * Event processing function.
   * The method is called by the framework to process one event. After 
   * preparation of data structures. The call is redirected to DumpEvent.
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
   * Return @ref AliHLTComponent::kSink type as component type.
   * @return component type id
   */
  TComponentType GetComponentType() { return AliHLTComponent::kSink;}

  /**
   * Default implementation for all data sinks.
   * There are no output data types.
   */
  AliHLTComponentDataType GetOutputDataType();

  /**
   * Default implementation for all data sinks.
   * There is no output data.
   */
  void GetOutputDataSize( unsigned long& constBase, double& inputMultiplier );

 protected:

  /**
   * The low-level data processing method for the component.
   * This is the custom processing method and can be overloaded by 
   * the component.
   * @param evtData       event data structure
   * @param blocks        input data block descriptors
   * @param trigData	  trigger data structure
   */
  virtual int DumpEvent( const AliHLTComponentEventData& evtData,
		 const AliHLTComponentBlockData* blocks, 
		 AliHLTComponentTriggerData& trigData );

  /**
   * The high-level data processing method.
   * This is the default processing method; the method is called
   * if no low level @ref DumpEvent method is overloaded by the component.
   * @param evtData       event data structure
   * @param trigData	  trigger data structure
   */
  virtual int DumpEvent( const AliHLTComponentEventData& evtData, AliHLTComponentTriggerData& trigData);


  ClassDef(AliHLTDataSink, 1)
};
#endif
