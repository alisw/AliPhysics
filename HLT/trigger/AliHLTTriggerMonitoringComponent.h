//-*- Mode: C++ -*-
// $Id: AliHLTTriggerMonitoringComponent.h 24328 2008-03-06 13:26:00Z richterm $
#ifndef ALIHLTTRIGGERMONITORINGCOMPONENT_H
#define ALIHLTTRIGGERMONITORINGCOMPONENT_H

/* This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTTriggerMonitoringComponent.h
    @author Jochen Thaeder
    @date   
    @brief  Produces a monitoring Trigger for AliEve
*/

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt   

#include "AliHLTProcessor.h"
#include "AliHLTEventSummary.h"

/**
 * @class  AliHLTTriggerMonitoringComponent
 * @brief  Produces a monitoring Trigger for AliEve
 *
 * This class produces a trigger accoring to events with more than x tracks, 
 * with associated clusters more the y and sends it out for every event.
 * x and y can be set as option
 *
 * @ingroup alihlt_trigger
 */

class AliHLTTriggerMonitoringComponent : public AliHLTProcessor {
  
public:
  
  /*
   * ---------------------------------------------------------------------------------
   *                            Constructor / Destructor
   * ---------------------------------------------------------------------------------
   */

  /** constructor */
  AliHLTTriggerMonitoringComponent();

  /** destructor */
  virtual ~AliHLTTriggerMonitoringComponent();

  /*
   * ---------------------------------------------------------------------------------
   * Public functions to implement AliHLTComponent's interface.
   * These functions are required for the registration process
   * ---------------------------------------------------------------------------------
   */
  
  /** interface function, see AliHLTComponent for description */
  const char* GetComponentID();

  /** interface function, see AliHLTComponent for description */
  void GetInputDataTypes( AliHLTComponentDataTypeList& list );

  /** interface function, see AliHLTComponent for description */
  AliHLTComponentDataType GetOutputDataType();

  /** interface function, see AliHLTComponent for description */
  virtual void GetOutputDataSize( unsigned long& constBase, double& inputMultiplier );

  /** interface function, see AliHLTComponent for description */
  AliHLTComponent* Spawn();

 protected:

  using AliHLTProcessor::DoEvent;

  /*
   * ---------------------------------------------------------------------------------
   * Protected functions to implement AliHLTComponent's interface.
   * These functions provide initialization as well as the actual processing
   * capabilities of the component. 
   * ---------------------------------------------------------------------------------
   */
  
  /** Initialization */
  Int_t DoInit( int argc, const char** argv );

  /** DeInitialization */
  Int_t DoDeinit();
  
  /** EventLoop */
  Int_t DoEvent( const AliHLTComponentEventData& evtData, AliHLTComponentTriggerData& trigData);

  ///////////////////////////////////////////////////////////////////////////////////

private:
 
  /** copy constructor prohibited */
  AliHLTTriggerMonitoringComponent (const AliHLTTriggerMonitoringComponent&);

  /** assignment operator prohibited */
  AliHLTTriggerMonitoringComponent& operator= (const AliHLTTriggerMonitoringComponent&);

  // ------------------------------------------------------------------------------------------

  /** Process trigger */
  void Trigger();

  /*
   * ---------------------------------------------------------------------------------
   *                             Members - private
   * ---------------------------------------------------------------------------------
   */

  /** Cut on total number of tracks */
  Int_t fTotalTrackCut;

  /** Cut on number of long tracks */
  Int_t fLongTrackCut;
  
  /** Event summary class*/
  AliHLTEventSummary* fEventSummary;       //! transient
  
  ClassDef(AliHLTTriggerMonitoringComponent, 0);

};
#endif

