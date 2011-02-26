//-*- Mode: C++ -*-
// $Id$
#ifndef ALIHLTRUNSUMMARYPRODUCERCOMPONENT_H
#define ALIHLTRUNSUMMARYPRODUCERCOMPONENT_H

/* This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTRunSummaryProducerComponent.h
    @author Jochen Thaeder
    @date   
    @brief  Produces a run summary
*/

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt   

#include "AliHLTProcessor.h"
#include "AliHLTRunSummary.h"
#include "AliHLTEventSummary.h"

/**
 * @class  AliHLTRunSummaryProducerComponent
 * @brief  Produces a run summary
 *
 * This class produces a run summary, updating informations for whole run 
 * and sends it out for every event, in order to have it displayed with AliEve
 *
 * @ingroup alihlt_run_statistics alihlt_trigger
 */

class AliHLTRunSummaryProducerComponent : public AliHLTProcessor {
  
public:
  
  /** constructor */
  AliHLTRunSummaryProducerComponent();
  /** destructor */
  virtual ~AliHLTRunSummaryProducerComponent();

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
  /** Process event summary data block 
   *  @param  eventSummary to @see AliHLTEventSummary
   */
  void ProcessEventSummary( AliHLTEventSummary* eventSummary );

  /** Process trigger data block 
   *  @param  trigData to @see AliHLTComponentTriggerData
   */
  void ProcessTriggerData( AliHLTComponentTriggerData& trigData );

private:
 
  /** copy constructor prohibited */
  AliHLTRunSummaryProducerComponent (const AliHLTRunSummaryProducerComponent&);

  /** assignment operator prohibited */
  AliHLTRunSummaryProducerComponent& operator= (const AliHLTRunSummaryProducerComponent&);

  /** Run summary class*/
  AliHLTRunSummary* fRunSummary;                                                                     //! transient
 
  ClassDef(AliHLTRunSummaryProducerComponent, 0);

};
#endif

