// XEmacs -*-C++-*-
// @(#) $Id$

#ifndef ALIHLTTPCGLOBALMERGERCOMPONENT_H
#define ALIHLTTPCGLOBALMERGERCOMPONENT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* AliHLTTPCGlobalMergerComponent
 */

#include "AliHLTProcessor.h"
#include "AliHLTTPCDefinitions.h"

class AliL3GlobalMerger;
class AliL3Vertex;

class AliHLTTPCGlobalMergerComponent : public AliHLTProcessor
    {
    public:
	AliHLTTPCGlobalMergerComponent();
	virtual ~AliHLTTPCGlobalMergerComponent();

	// Public functions to implement AliHLTComponent's interface.
	// These functions are required for the registration process

	const char* GetComponentID();
	void GetInputDataTypes( vector<AliHLTComponent_DataType>& list);
	AliHLTComponent_DataType GetOutputDataType();
	virtual void GetOutputDataSize( unsigned long& constBase, double& inputMultiplier );
	AliHLTComponent* Spawn();
	
    protected:
	
	void SetMergerParameters(Double_t maxy=2.0,Double_t maxz=3.0,Double_t maxkappa=0.003,
				 Double_t maxpsi=0.1,Double_t maxtgl=0.05);

	// Protected functions to implement AliHLTComponent's interface.
	// These functions provide initialization as well as the actual processing
	// capabilities of the component. 

	int DoInit( int argc, const char** argv );
	int DoDeinit();
	int DoEvent( const AliHLTComponent_EventData& evtData, const AliHLTComponent_BlockData* blocks, 
		     AliHLTComponent_TriggerData& trigData, AliHLTUInt8_t* outputPtr, 
		     AliHLTUInt32_t& size, vector<AliHLTComponent_BlockData>& outputBlocks );
	
    private:

	AliL3GlobalMerger* fGlobalMerger;
	AliL3Vertex* fVertex;

	struct SliceData
	    {
		int fSlice;
		const AliHLTComponent_BlockData* fVertexBlock;
		unsigned fVertexBlockIndex;
		const AliHLTComponent_BlockData* fTrackletBlock;
		unsigned fTrackletBlockIndex;
	    };

	ClassDef(AliHLTTPCGlobalMergerComponent, 0)

    };
#endif
