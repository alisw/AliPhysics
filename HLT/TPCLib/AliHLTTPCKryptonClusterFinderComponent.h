// $Id$

#ifndef ALIHLTTPCKRYPTONCLUSTERFINDERCOMPONENT_H
#define ALIHLTTPCKRYPTONCLUSTERFINDERCOMPONENT_H

//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/** @file   AliHLTTPCKryptonClusterFinderComponent.h
    @author Kenneth Aamodt, Kalliopi Kanaki
    @date   
    @brief  The TPC krypton cluster finder component.
*/

#include "AliHLTProcessor.h"

class AliHLTTPCKryptonClusterFinder;
class AliHLTTPCDigitReader;

/**
 * @class AliHLTTPCKryptonClusterFinderComponent
 * Component for the krypton ClusterFinder
 * @ingroup alihlt_tpc_components
 */
class AliHLTTPCKryptonClusterFinderComponent : public AliHLTProcessor
    {
    public:

        /**
         * constructor 
         */
	AliHLTTPCKryptonClusterFinderComponent();
	/** destructor */
	virtual ~AliHLTTPCKryptonClusterFinderComponent();

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

    protected:
	
	// Protected functions to implement AliHLTComponent's interface.
	// These functions provide initialization as well as the actual processing
	// capabilities of the component. 

	int DoInit( int argc, const char** argv );
	int DoDeinit();
	int DoEvent( const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks, 
		     AliHLTComponentTriggerData& trigData, AliHLTUInt8_t* outputPtr, 
		     AliHLTUInt32_t& size, vector<AliHLTComponentBlockData>& outputBlocks );

	int Reconfigure(const char* cdbEntry, const char* chainId);
	
	int Configure(const char* arguments);
      
	using AliHLTProcessor::DoEvent;

    private:
	/** copy constructor prohibited */
	AliHLTTPCKryptonClusterFinderComponent(const AliHLTTPCKryptonClusterFinderComponent&);
	/** assignment operator prohibited */
	AliHLTTPCKryptonClusterFinderComponent& operator=(const AliHLTTPCKryptonClusterFinderComponent&);
	/** the cluster finder object */
	AliHLTTPCKryptonClusterFinder* fKryptonClusterFinder;                               //!transient
	/** the reader object for data decoding */
	AliHLTTPCDigitReader* fReader;                                                      //!transient

	AliHLTUInt32_t fSpecification;

	ClassDef(AliHLTTPCKryptonClusterFinderComponent, 1)
};
#endif
