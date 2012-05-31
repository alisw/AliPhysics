//-*- Mode: C++ -*-
// @(#) $Id$

#ifndef ALIHLTDATASOURCE_H
#define ALIHLTDATASOURCE_H
/* This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTDataSource.h
    @author Matthias Richter
    @date   
    @brief  Base class declaration for HLT data source components.
    @note   The class is used in Offline (AliRoot) context
*/

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt   

#include "AliHLTComponent.h"

/**
 * @class AliHLTDataSource
 * Base class of HLT data source components.
 * The class provides a common interface for the implementation of HLT data
 * source components. The child class must implement the functions:
 * - @ref DoInit (optional)
 * - @ref DoDeinit (optional)
 * - @ref GetEvent
 * - @ref GetComponentID
 * - @ref GetOutputDataType
 * - @ref GetOutputDataSize
 * - @ref Spawn
 *
 * @ingroup alihlt_component
 */
class AliHLTDataSource : public AliHLTComponent {
 public:
  /** standard constructor */
  AliHLTDataSource();
  /** standard destructor */
  virtual ~AliHLTDataSource();

  /**
   * Event processing function.
   * The method is called by the framework to process one event. After 
   * preparation of data structures. The call is redirected to GetEvent.
   * @return neg. error code if failed                                <br>
   *         -ENOSPC      output buffer too small
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
   * Return @ref AliHLTComponent::kSource type as component type.
   * @return component type id
   */
  TComponentType GetComponentType() { return AliHLTComponent::kSource;}

  /**
   * Default implementation for all data sources.
   * There are no input data types.
   */
  void GetInputDataTypes( vector<AliHLTComponentDataType>& list);

protected:

  /**
   * The low-level data processing method for the component.
   * This is the custom processing method and can be overloaded by 
   * the component.
   * @param [in] evtData       event data structure
   * @param [in] trigData	  trigger data structure
   * @param [in] outputPtr	  pointer to target buffer
   * @param [in,out] size	  <i>input</i>: size of target buffer
   *            	  <i>output</i>:size of produced data
   * @param [in] outputBlocks  list to receive output block descriptors
   * @return neg. error code if failed
   */
  virtual int GetEvent( const AliHLTComponentEventData& evtData,
		AliHLTComponentTriggerData& trigData,
		AliHLTUInt8_t* outputPtr, 
		AliHLTUInt32_t& size,
		vector<AliHLTComponentBlockData>& outputBlocks );

  /**
   * The high-level data processing method.
   * This is the default processing method; the method is called
   * if no low level @ref GetEvent method is overloaded by the component.
   * @param evtData       event data structure
   * @param trigData	  trigger data structure
   * @return neg. error code if failed
   */
  virtual int GetEvent( const AliHLTComponentEventData& evtData, AliHLTComponentTriggerData& trigData);

private:

  ClassDef(AliHLTDataSource, 3)
};
#endif
