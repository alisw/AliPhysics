// XEMacs -*-C++-*-
// @(#) $Id$

#ifndef ALIHLTDUMMYCOMPONENT_H
#define ALIHLTDUMMYCOMPONENT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* AliHLTDummyComponent
** An implementiation of a dummy component that just copies its input data
** as a test, demonstration, and example of the new HLT component system.
 */

#include "AliHLTProcessor.h"

class AliHLTDummyComponent : public AliHLTProcessor
    {
    public:
	AliHLTDummyComponent();
	virtual ~AliHLTDummyComponent();

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

	// The size of the output data produced, as a percentage of the input data's size.
	// Can be greater than 100 (%)
	unsigned fOutputPercentage;
	
	ClassDef(AliHLTDummyComponent, 0)

    };
#endif
