//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/** @file   AliHLTGlobalTrackMatcherComponent.h
    @author Sergey Gorbunov
    @brief  Component for monitor V0 physics
*/



#ifndef ALIHLTGLOBALTRACKMATCHERCOMPONENT_H
#define ALIHLTGLOBALTRACKMATCHERCOMPONENT_H
#include "AliHLTComponentBenchmark.h"

class AliHLTProcessor;
class AliHLTGlobalTrackMatcher;
class AliHLTCaloClusterReader;
class TObjArray;
/**
 * @class AliHLTTPCV0HistoComponent
 * Component for monitor V0 physics 
 */
class AliHLTGlobalTrackMatcherComponent : public AliHLTProcessor
{
public:
  /** default constructor */
  AliHLTGlobalTrackMatcherComponent();
  /** destructor */
  virtual ~AliHLTGlobalTrackMatcherComponent();

  // Public functions to implement AliHLTComponent's interface.
  // These functions are required for the registration process

  /** interface function, see AliHLTComponent for description */
  const char* GetComponentID();
  /** interface function, see AliHLTComponent for description */
  void GetInputDataTypes(AliHLTComponentDataTypeList& list);
  /** interface function, see AliHLTComponent for description */
  AliHLTComponentDataType GetOutputDataType();
  /** interface function, see AliHLTComponent for description */
  virtual void GetOutputDataSize( unsigned long& constBase, double& inputMultiplier );
  /** interface function, see AliHLTComponent for description */
  AliHLTComponent* Spawn();

    ///Inherited from AliHLTComponent: Get list of OCDB objects
  void GetOCDBObjectDescription( TMap* const targetMap); //Methods
  
protected:

  // Protected functions to implement AliHLTComponent's interface.
  // These functions provide initialization as well as the actual processing
  // capabilities of the component. 

  /** interface function, see AliHLTComponent for description */
  int DoInit( int argc, const char** argv );
  /** interface function, see AliHLTComponent for description */
  int DoDeinit();
  /** interface function, see AliHLTComponent for description */
  
  //FOR METHOD extrapolation:
    /// inherited from AliHLTComponent: handle re-configuration event
  int Reconfigure(const char* cdbEntry, const char* chainId);

  /// inherited from AliHLTComponent, scan one argument and
  /// its parameters
  int ScanConfigurationArgument(int argc, const char** argv);
  //FOR METHOD
  
    /// the default configuration entry for this component
  const char* fOCDBEntry; // Method for TrackMatcher
  
  int DoEvent( const AliHLTComponentEventData& /*evtData*/, AliHLTComponentTriggerData& trigData );

  //int Reconfigure(const char* cdbEntry, const char* chainId);

  using AliHLTProcessor::DoEvent;
  
  Int_t fMethod; //TString for method choice for extrapolation
  
private:
  /** copy constructor prohibited */
  AliHLTGlobalTrackMatcherComponent(const AliHLTGlobalTrackMatcherComponent&);
  /** assignment operator prohibited */
  AliHLTGlobalTrackMatcherComponent& operator=(const AliHLTGlobalTrackMatcherComponent&);
  /**
   * Configure the component.
   * Parse a string for the configuration arguments and set the component
   * properties.
   */ 
  //  int Configure(const char* arguments);
  

  AliHLTGlobalTrackMatcher * fTrackMatcher;   //Instance of the track matcher class

  Int_t fNEvents;    //Number of events processed

  Double_t fBz;   //Magnetic field of event

  AliHLTCaloClusterReader * fClusterReader;   //Instance of helper class to read calorimeter structs

  TObjArray * fTrackArray;

  ClassDef(AliHLTGlobalTrackMatcherComponent, 0);

};
#endif
