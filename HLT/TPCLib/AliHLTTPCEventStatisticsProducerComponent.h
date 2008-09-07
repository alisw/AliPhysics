//-*- Mode: C++ -*-
// $Id$
#ifndef ALIHLTTPCEVENTSTATISTICSPRODUCERCOMPONENT_H
#define ALIHLTTPCEVENTSTATISTICSPRODUCERCOMPONENT_H

/* This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTTPCEventStatisticsProducerComponent.h
    @author Jochen Thaeder
    @date   
    @brief  Component for the @see AliHLTTPCEventStatisticsProducer class
*/

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt   

#include "AliHLTProcessor.h"

#include "AliHLTTPCEventStatistics.h"
#include "AliHLTTPCTrackArray.h"
#include "AliESDEvent.h"
#include "TTree.h"

/**
 * @class  AliHLTTPCEventStatisticsProducerComponent
 * @brief  Component for the @see AliHLTTPCEventStatisticsProducer class
 *
 *
 * @ingroup alihlt_run_statistics alihlt_trigger
 */

class AliHLTTPCEventStatisticsProducerComponent : public AliHLTProcessor {
  
public:
  
  /** constructor */
  AliHLTTPCEventStatisticsProducerComponent();
  /** destructor */
  virtual ~AliHLTTPCEventStatisticsProducerComponent();

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

  /** Initialize for next event */
  void InitializeEvent();                  

  /** Add new cluster block 
   * @param ptr   Pointer to data block
   * @param size  Size of data block
   * @param slice Slice
   * @param patch Patch
   */
  void AddClusters( void* ptr, Int_t slice, Int_t patch );

  /** Add new tracks block
   * @param ptr   Pointer to data block
   * @param size  Size of data block
   * @param slice Slice, default is -1 for GlobalMerger tracks
   */
  void AddTracks( void* ptr, Int_t slice=-1 );     

  /** Add ESD
   * @param esdEvent Pointer to AliESDEvent
   */
  void AddESD( TTree* esdTree );     

  /** Process even -> get process statistics */
  void ProcessEvent();                          

  /** Fill tracks per slice ( @see fNTracksPerSlice ), if global tracks */
  void FillTracksPerSlice();    

  /** Get ptr to @see AliHLTTPCEventStatistics, is a TObject */
  AliHLTTPCEventStatistics* GetEventStatistics() { return fEvStat; }    

private:
 
  /** copy constructor prohibited */
  AliHLTTPCEventStatisticsProducerComponent (const AliHLTTPCEventStatisticsProducerComponent&);

  /** assignment operator prohibited */
  AliHLTTPCEventStatisticsProducerComponent& operator= (const AliHLTTPCEventStatisticsProducerComponent&);

  /** Threshold for long tracks */
  Int_t fClusterThreshold;                                                                           // see above

  /** Event Statistics class*/
  AliHLTTPCEventStatistics* fEvStat;                                                                 //! transient

  /** Tracks per slice */
  Int_t fNTracksPerSlice[36];                                                                        // see above

  /** Track container of class @see AliHLTTPCTrackArray */
  AliHLTTPCTrackArray *fTracks;                                                                      //! transient

  /** If tracks in global coordinates @see fTracksPerSlice has still to be filled */
  Bool_t fGlobalTracks;                                                                              // see above

  /** Number of slices with tracks */
  Int_t fNSlice;                                                                                     // see above

  ClassDef(AliHLTTPCEventStatisticsProducerComponent, 0);

};
#endif

