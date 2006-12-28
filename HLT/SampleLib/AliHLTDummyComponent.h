// $Id$

#ifndef ALIHLTDUMMYCOMPONENT_H
#define ALIHLTDUMMYCOMPONENT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTDummyComponent.h
    @author Timm Steinbeck, Matthias Richter
    @date   
    @brief  Declaration of a dummy component. */


#include "AliHLTProcessor.h"

/**
 * @class AliHLTDummyComponent
 * @brief A dummy HLT processing component. 
 *
 * An implementiation of a dummy component that just copies its input data
 * as a test, demonstration, and example of the HLT component scheme.
 * @ingroup alihlt_tutorial
 */
class AliHLTDummyComponent : public AliHLTProcessor
    {
    public:
	AliHLTDummyComponent();
	virtual ~AliHLTDummyComponent();

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

	// The size of the output data produced, as a percentage of the input data's size.
	// Can be greater than 100 (%)
	unsigned fOutputPercentage;
	
	ClassDef(AliHLTDummyComponent, 0)

    };
#endif
