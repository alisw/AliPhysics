//-*- Mode: C++ -*-
// $Id$
#ifndef ALIHLTTPCRUNSTATISTICSPRODUCERCOMPONENT_H
#define ALIHLTTPCRUNSTATISTICSPRODUCERCOMPONENT_H

/* This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTTPCRunStatisticsProducerComponent.h
 *  @author Jochen Thaeder
 *  @date   
 *  @brief  Component for producing the @see AliHLTTPCRunStatistics
 */

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt   

#include "AliHLTProcessor.h"

#include "AliHLTTPCRunStatistics.h"
#include "AliHLTTPCEventStatistics.h"

/**
 * @class  AliHLTTPCRunStatisticsProducerComponent
 * @brief  Component for producing the @see AliHLTTPCRunStatistics
 *
 * @ingroup alihlt_run_statistics alihlt_trigger
 */

class AliHLTTPCRunStatisticsProducerComponent : public AliHLTProcessor {
  
public:
  
  /** constructor */
  AliHLTTPCRunStatisticsProducerComponent();
  /** destructor */
  virtual ~AliHLTTPCRunStatisticsProducerComponent();

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
  
  /** Initialize the trigger component. */
  Int_t DoInit( int argc, const char** argv );

  /** DeInitialize the trigger component. */
  Int_t DoDeinit();
  
  /** Process the data in the trigger component */
  Int_t DoEvent( const AliHLTComponentEventData& evtData, AliHLTComponentTriggerData& trigData);

  // ------------------------------------------------------------------------------------------

  /** Get Process event statistics coming from @see AliHLTTPCEventStatistics
   *  @param evStat  event statistics as @see AliHLTTPCEventStatistics
   */
  void ProcessEventStatistics( AliHLTTPCEventStatistics* evStat );

  /** Get ptr to @see AliHLTTPCRunStatistics, is a TObject */
  AliHLTTPCRunStatistics* GetRunStatistics() { return fRunStat; }  

private:
 
  /** copy constructor prohibited */
  AliHLTTPCRunStatisticsProducerComponent (const AliHLTTPCRunStatisticsProducerComponent&);

  /** assignment operator prohibited */
  AliHLTTPCRunStatisticsProducerComponent& operator= (const AliHLTTPCRunStatisticsProducerComponent&);

  /** Event Statistics class*/
  AliHLTTPCRunStatistics* fRunStat;                                                                  //! transient

  /** If run header is set */
  Bool_t fIsHeader;                                                                                  // see above
  
  ClassDef(AliHLTTPCRunStatisticsProducerComponent, 0);

};
#endif

