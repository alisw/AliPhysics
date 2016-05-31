// -*- Mode: C++ -*-

#ifndef ALIHLTTPCCLUSTERTRANSFORMATIONCOMPONENT_H
#define ALIHLTTPCCLUSTERTRANSFORMATIONCOMPONENT_H

//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/** @file   AliHLTTPCClusterTransformationComponent.h
    @author Sergey Gorbunov
    @date   
    @brief
*/

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "AliHLTProcessor.h"
#include "AliHLTComponentBenchmark.h"

class AliHLTTPCClusterTransformation;


/**
 * @class AliHLTTPCClusterTransformationComponent
 * @ingroup alihlt_tpc_components
 */

class AliHLTTPCClusterTransformationComponent : public AliHLTProcessor {
    
public:

  /** standard constructor */    
  AliHLTTPCClusterTransformationComponent();           
  /** destructor */
  virtual ~AliHLTTPCClusterTransformationComponent();

  // Public functions to implement AliHLTComponent's interface.
  // These functions are required for the registration process
      
  /** interface function, see AliHLTComponent for description */
  const char* GetComponentID();							     
  /** interface function, see AliHLTComponent for description */
  void GetInputDataTypes( vector<AliHLTComponentDataType>& list);			     
  /** interface function, see AliHLTComponent for description */
  AliHLTComponentDataType GetOutputDataType();					     
  /** interface function, see AliHLTComponent for description */
  int GetOutputDataTypes(AliHLTComponentDataTypeList& tgtList);			   
  /** interface function, see AliHLTComponent for description */
  virtual void GetOutputDataSize( unsigned long& constBase, double& inputMultiplier ); 
  /** interface function, see AliHLTComponent for description */
  AliHLTComponent* Spawn(); 							   
  /** interface function, see @ref AliHLTComponent for description */
  void GetOCDBObjectDescription( TMap* const targetMap);
  
protected:
	
  // Protected functions to implement AliHLTComponent's interface.
  // These functions provide initialization as well as the actual processing capabilities of the component. 

  int DoInit( int argc, const char** argv );
  int DoEvent( const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks, 
		     AliHLTComponentTriggerData& trigData, AliHLTUInt8_t* outputPtr, 
		     AliHLTUInt32_t& size, vector<AliHLTComponentBlockData>& outputBlocks );
  int DoDeinit();
  
  int Reconfigure(const char* cdbEntry, const char* chainId);  

  int ScanConfigurationArgument(int argc, const char** argv);

  using AliHLTProcessor::DoEvent;
  
private:
   
          
  /** copy constructor prohibited */
  AliHLTTPCClusterTransformationComponent(const AliHLTTPCClusterTransformationComponent&);

  /** assignment operator prohibited */
  AliHLTTPCClusterTransformationComponent& operator=(const AliHLTTPCClusterTransformationComponent&);

  static const char* fgkOCDBEntryClusterTransformation; //!transient
  static AliHLTTPCClusterTransformation fgTransform; //!transient
  static Bool_t fgTimeInitialisedFromEvent; //!transient
  bool fOfflineMode; //Run in offline mode (for hlt simulation in offline)
  int fInitializeByObjectInDoEvent;	//Do not initialize the transformation but wait for transformation objects arriving in DoEvent loop
  bool fInitialized;	//Are we initialized?
  bool fTPCPresent;	//Is TPC present in GRP, if not skip init

  AliHLTComponentBenchmark fBenchmark; // benchmarks

  ClassDef(AliHLTTPCClusterTransformationComponent, 0)
};

#endif
