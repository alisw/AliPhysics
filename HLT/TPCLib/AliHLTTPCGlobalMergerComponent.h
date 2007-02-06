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

class AliHLTTPCGlobalMerger;
class AliHLTTPCVertex;

class AliHLTTPCGlobalMergerComponent : public AliHLTProcessor
    {
    public:
	AliHLTTPCGlobalMergerComponent();
	virtual ~AliHLTTPCGlobalMergerComponent();

	// Public functions to implement AliHLTComponent's interface.
	// These functions are required for the registration process

	const char* GetComponentID();
	void GetInputDataTypes( vector<AliHLTComponentDataType>& list);
	AliHLTComponentDataType GetOutputDataType();
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
	int DoEvent( const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks, 
		     AliHLTComponentTriggerData& trigData, AliHLTUInt8_t* outputPtr, 
		     AliHLTUInt32_t& size, vector<AliHLTComponentBlockData>& outputBlocks );
	
    private:

      AliHLTTPCGlobalMerger* fGlobalMerger; //!
      AliHLTTPCVertex* fVertex; //!

	struct SliceData
	    {
		int fSlice;
		const AliHLTComponentBlockData* fVertexBlock;
		unsigned fVertexBlockIndex;
		const AliHLTComponentBlockData* fTrackletBlock;
		unsigned fTrackletBlockIndex;
	    };

	ClassDef(AliHLTTPCGlobalMergerComponent, 0)

    };
#endif
