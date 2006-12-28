// XEmacs -*-C++-*-
// @(#) $Id$

#ifndef ALIHLTTPCVERTEXFINDERCOMPONENT_H
#define ALIHLTTPCVERTEXFINDERCOMPONENT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* AliHLTTPCVertexFinderComponent
 */

#include "AliHLTProcessor.h"
#include "AliHLTTPCDefinitions.h"

class AliHLTTPCVertexFinder;

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

	AliHLTTPCVertexFinder* fVertexFinder;

	ClassDef(AliHLTTPCVertexFinderComponent, 0)

    };
#endif
