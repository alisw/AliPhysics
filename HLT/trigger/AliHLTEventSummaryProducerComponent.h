 //-*- Mode: C++ -*-
// $Id$
#ifndef ALIHLTEVENTSUMMARYPRODUCERCOMPONENT_H
#define ALIHLTEVENTSUMMARYPRODUCERCOMPONENT_H

/* This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTEventSummaryProducerComponent.h
    @author Jochen Thaeder
    @date   
    @brief  Produces a event summary as @see AliHLTEventSummary
*/

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt   

#include "AliHLTProcessor.h"
#include "AliHLTEventSummary.h"

/**
 * @class  AliHLTEventSummaryProducerComponent
 * @brief  Produces a event summary as @see AliHLTEventSummary
 *
 * This class produces a event summary, updating informations for all subdetectors 
 * and sends it out for every event.
 *
 * @ingroup alihlt_run_statistics alihlt_trigger
 */

class AliHLTEventSummaryProducerComponent : public AliHLTProcessor {
  
public:
  
  /** constructor */
  AliHLTEventSummaryProducerComponent();
  /** destructor */
  virtual ~AliHLTEventSummaryProducerComponent();

  // Public functions to implement AliHLTComponent's interface.
  // These functions are required for the registration process
  
  const char* GetComponentID();
  void GetInputDataTypes( AliHLTComponentDataTypeList& list );

  AliHLTComponentDataType GetOutputDataType();

  virtual void GetOutputDataSize( unsigned long& constBase, double& inputMultiplier );
  AliHLTComponent* Spawn();

 protected:

  using AliHLTProcessor::DoEvent;

  // Protected functions to implement AliHLTComponent's interface.
  // These functions provide initialization as well as the actual processing
  // capabilities of the component. 
  
  /** Initialize the component. */
  Int_t DoInit( int argc, const char** argv );

  /** DeInitialize the component. */
  Int_t DoDeinit();
  
  /** Process the data in the component */
  Int_t DoEvent( const AliHLTComponentEventData& evtData, AliHLTComponentTriggerData& trigData);

  // ------------------------------------------------------------------------------------------

  /** Process trigger data block 
   *  @param  trigData to @see AliHLTComponentTriggerData
   */
  void ProcessTriggerData( AliHLTComponentTriggerData& trigData );

private:
 
  /** copy constructor prohibited */
  AliHLTEventSummaryProducerComponent (const AliHLTEventSummaryProducerComponent&);

  /** assignment operator prohibited */
  AliHLTEventSummaryProducerComponent& operator= (const AliHLTEventSummaryProducerComponent&);

  /** Event summary class*/
  AliHLTEventSummary* fEventSummary;                                                                  //! transient

  ClassDef(AliHLTEventSummaryProducerComponent, 0);

};
#endif

