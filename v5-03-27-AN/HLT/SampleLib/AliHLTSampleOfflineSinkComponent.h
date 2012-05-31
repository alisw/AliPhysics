//-*- Mode: C++ -*-
// @(#) $Id$

#ifndef ALIHLTSAMPLEOFFLINESINKCOMPONENT_H
#define ALIHLTSAMPLEOFFLINESINKCOMPONENT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx sink for full Copyright notice                               */

/** @file   AliHLTSampleOfflineSinkComponent.h
    @author Matthias Richter
    @date   
    @brief  This is a sample offline interface component.
*/

#include "AliHLTOfflineDataSink.h"

class AliLoader;

/**
 * @class AliHLTSampleOfflineSinkComponent
 * This is a sample offline interface component.
 *
 * @ingroup alihlt_system
 */
class AliHLTSampleOfflineSinkComponent : public AliHLTOfflineDataSink {
 public:
  /**
   * standard constructor. The agent is automatically registered in the
   * global agent manager
   */
  AliHLTSampleOfflineSinkComponent();
  /** destructor */
  virtual ~AliHLTSampleOfflineSinkComponent();

  /**
   * Get the id of the component.
   * Each component is identified by a unique id.
   * The function is pure virtual and must be implemented by the child class.
   * @return component id (string)
   */
  const char* GetComponentID();

  /**
   * Get the input data types of the component.
   * The function is pure virtual and must be implemented by the child class.
   * @return list of data types in the vector reference
   */
  virtual void GetInputDataTypes( vector<AliHLTComponentDataType>& );

  /**
   * Spawn function.
   * Each component must implement a spawn function to create a new instance of 
   * the class. Basically the function must return <i>new <b>my_class_name</b></i>.
   * @return new class instance
   */
  virtual AliHLTComponent* Spawn();

  /**
   * Fill ESD for one event.
   * Fill the ESD with the previously reconstructed data. Interface method called by
   * the AliHLTOfflineInterface framework.
   * @param eventNo       event No. \em Note: this is an internal enumeration of the
   *                      processed events.
   * @param runLoader     the AliRoot runloader
   * @param esd           an AliESDEvent instance
   * @return neg. error code if failed 
   */
  int FillESD(int eventNo, AliRunLoader* runLoader, AliESDEvent* esd);

 protected:
  /**
   * Init method.
   */
  int DoInit( int argc, const char** argv );

  /**
   * Deinit method.
   */
  int DoDeinit();

  /**
   * Data sink method.
   * @param evtData       event data structure
   * @param trigData	  trigger data structure
   * @return
   */
  int DumpEvent(const AliHLTComponentEventData& evtData,
		AliHLTComponentTriggerData& trigData);

  using AliHLTOfflineDataSink::DumpEvent;

 private:
  ClassDef(AliHLTSampleOfflineSinkComponent, 0);
};

#endif
