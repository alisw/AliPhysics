// XEmacs -*-C++-*-
// @(#) $Id$

#ifndef ALIHLTTPCVERTEXFINDERCOMPONENT_H
#define ALIHLTTPCVERTEXFINDERCOMPONENT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "AliHLTProcessor.h"

class AliHLTTPCVertexFinder;

/**
 * @class AliHLTTPCVertexFinderComponent
 * A vertex finder component for the TPC.
 * This component has never been tested in the new framework and needs certainly
 * some investigation.
 */
class AliHLTTPCVertexFinderComponent : public AliHLTProcessor
    {
    public:
	AliHLTTPCVertexFinderComponent();
	virtual ~AliHLTTPCVertexFinderComponent();

	// Public functions to implement AliHLTComponent's interface.
	// These functions are required for the registration process

	const char* GetComponentID();
	void GetInputDataTypes( vector<AliHLTComponentDataType>& list);
	AliHLTComponentDataType GetOutputDataType();
	virtual void GetOutputDataSize( unsigned long& constBase, double& inputMultiplier );
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
	
    private:
      /** copy constructor prohibited */
      AliHLTTPCVertexFinderComponent(const AliHLTTPCVertexFinderComponent&);
      /** assignment operator prohibited */
      AliHLTTPCVertexFinderComponent& operator=(const AliHLTTPCVertexFinderComponent&);

      AliHLTTPCVertexFinder* fVertexFinder; //! transient

      ClassDef(AliHLTTPCVertexFinderComponent, 0);
    };
#endif
