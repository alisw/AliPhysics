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

class AliL3VertexFinder;

class AliHLTTPCVertexFinderComponent : public AliHLTProcessor
    {
    public:
	AliHLTTPCVertexFinderComponent();
	virtual ~AliHLTTPCVertexFinderComponent();

	// Public functions to implement AliHLTComponent's interface.
	// These functions are required for the registration process

	const char* GetComponentID();
	void GetInputDataTypes( vector<AliHLTComponent_DataType>& list);
	AliHLTComponent_DataType GetOutputDataType();
	virtual void GetOutputDataSize( unsigned long& constBase, double& inputMultiplier );
	AliHLTComponent* Spawn();
	
    protected:
	
	// Protected functions to implement AliHLTComponent's interface.
	// These functions provide initialization as well as the actual processing
	// capabilities of the component. 

	int DoInit( int argc, const char** argv );
	int DoDeinit();
	int DoEvent( const AliHLTComponent_EventData& evtData, const AliHLTComponent_BlockData* blocks, 
		     AliHLTComponent_TriggerData& trigData, AliHLTUInt8_t* outputPtr, 
		     AliHLTUInt32_t& size, vector<AliHLTComponent_BlockData>& outputBlocks );
	
    private:

	AliL3VertexFinder* fVertexFinder;

	ClassDef(AliHLTTPCVertexFinderComponent, 0)

    };
#endif
