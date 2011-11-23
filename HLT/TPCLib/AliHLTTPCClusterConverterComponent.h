// XEmacs -*-C++-*-
// $Id$

#ifndef ALIHLTTPCCLUSTERCONVERTERCOMPONENT_H
#define ALIHLTTPCCLUSTERCONVERTERCOMPONENT_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/// @file   AliHLTTPCClusterConverterComponent.h
/// @author Kalliopi Kanaki
/// @date   
/// @brief  The TPC cluster format conversion component.
///

#include "AliHLTProcessor.h"

struct AliHLTTPCSpacePointData;
struct AliHLTTPCTrackSegmentData;
class TObjArray;

/**
 * @class AliHLTTPCClusterConverterComponent
 * The TPC cluster format conversion component
 *
 * The purpose of the component is publishing the 
 * USED clusters of the tracks, translating them to
 * the offline format. For this reason, it subscribes
 * to the cluster finder, the slice tracker and the
 * global merger.
 *
 * @ingroup alihlt_tpc_components
 */
 
class AliHLTTPCClusterConverterComponent : public AliHLTProcessor{

public:
  /** default constructor */
  AliHLTTPCClusterConverterComponent();
  /** destructor */
  virtual ~AliHLTTPCClusterConverterComponent();

  // Public functions to implement AliHLTComponent's interface.
  // These functions are required for the registration process

  /** interface function, see AliHLTComponent for description */
  const char* GetComponentID();
  /** interface function, see AliHLTComponent for description */
  void GetInputDataTypes(AliHLTComponentDataTypeList& list);
  /** interface function, see AliHLTComponent for description */
  AliHLTComponentDataType GetOutputDataType();
  /** interface function, see AliHLTComponent for description */
  int GetOutputDataTypes(AliHLTComponentDataTypeList& tgtList);
  /** interface function, see AliHLTComponent for description */
  virtual void GetOutputDataSize( unsigned long& constBase, double& inputMultiplier );
  /** interface function, see AliHLTComponent for description */
  AliHLTComponent* Spawn();

protected:

  /** interface function, see AliHLTComponent for description */
  int DoInit( int argc, const char** argv );
  /** interface function, see AliHLTComponent for description */
  int DoDeinit();
  /** interface function, see AliHLTComponent for description */
  int DoEvent(const AliHLTComponentEventData& /*evtData*/, AliHLTComponentTriggerData& trigData);

  int Reconfigure(const char* cdbEntry, const char* chainId);

  using AliHLTProcessor::DoEvent;
  
private:
  /** copy constructor prohibited */
  AliHLTTPCClusterConverterComponent(const AliHLTTPCClusterConverterComponent&);
  /** assignment operator prohibited */
  AliHLTTPCClusterConverterComponent& operator=(const AliHLTTPCClusterConverterComponent&);
  /**
   * Configure the component.
   * Parse a string for the configuration arguments and set the component
   * properties.
   */  
  int Configure(const char* arguments);
  
  vector<AliHLTTPCSpacePointData>   fClusters; //! transient
  vector<AliHLTTPCTrackSegmentData> fTracks;   //!transient
  TObjArray *fOffArray;                        //!transient
  
  ClassDef(AliHLTTPCClusterConverterComponent, 0);

};
#endif
