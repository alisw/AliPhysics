#ifndef ALIHLTMUONFULLTRACKERCOMPONENT_H
#define ALIHLTMUONFULLTRACKERCOMPONENT_H

//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               */


/// \ingroup libAliHLTMUON.so
/// \class AliHLTMUONFullTrackerComponent
/// \brief Component to handle full tracker
///
//  Indranil Das, email : indra.das@saha.ac.in | indra.ehep@gmail.com , Saha Institute of Nuclear Physics

#include "AliHLTMUONProcessor.h"
#include "AliHLTMUONFullTracker.h"
#include <cassert>


class AliHLTMUONFullTrackerComponent : public AliHLTMUONProcessor {
 public:
  /// Constructor
  AliHLTMUONFullTrackerComponent();
  /// Destructor
  virtual ~AliHLTMUONFullTrackerComponent();

  /// Get the component ID
  const char* GetComponentID();
  /// Get Output Data Types
  void GetInputDataTypes( vector<AliHLTComponentDataType>& list);
  /// Get Output Data Types
  AliHLTComponentDataType GetOutputDataType();
  /// Output Data Types
  virtual int GetOutputDataTypes(AliHLTComponentDataTypeList& list);
  /// Output Data Sizes
  virtual void GetOutputDataSize( unsigned long& constBase, double& inputMultiplier );
  /// Spawn
  AliHLTComponent* Spawn();

 protected:
  // AliHLTComponent interface functions
  int DoInit( int argc, const char** argv );
  /// DeInitialisation
  int DoDeinit();
  /// Process events
  int DoEvent( const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks, 
	       AliHLTComponentTriggerData& trigData, AliHLTUInt8_t* outputPtr, 
	       AliHLTUInt32_t& size, vector<AliHLTComponentBlockData>& outputBlocks );
  /// Reconfigure
  int Reconfigure(const char* cdbEntry, const char* chainId);
  /// Read preprocessor values
  int ReadPreprocessorValues(const char* modules);

  using AliHLTProcessor::DoEvent;

 private:

  // Do not allow copying of this class.
  /// Not implemented.
  AliHLTMUONFullTrackerComponent(const AliHLTMUONFullTrackerComponent& /*obj*/);
  /// Not implemented.
  AliHLTMUONFullTrackerComponent& operator = (const AliHLTMUONFullTrackerComponent& /*obj*/);

  /**
   * Configure the component.
   * Parse a string for the configuration arguments and set the component
   * properties.
   *
   * This function illustrates the scanning of an argument string. The string
   * was presumably fetched from the CDB.
   */
  int Configure(const char* arguments);
  
  AliHLTMUONFullTracker *fTracker; /// Pointer to the full tracker
  
  ClassDef(AliHLTMUONFullTrackerComponent, 0)
    };
#endif
