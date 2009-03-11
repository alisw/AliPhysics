#ifndef ALIHLTTRDTRACKERV1COMPONENT_H
#define ALIHLTTRDTRACKERV1COMPONENT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTTRDTrackerV1Component.h
    @author Timm Steinbeck, Matthias Richter
    @date   
    @brief  Declaration of a TRDTracker component. */

#include "AliHLTProcessor.h"

class TFile;
class TTree;

class TGeoManager;
class AliCDBManager;
class AliMagF;
class AliTRDtrackerV1;
class AliTRDrecoParam;
class AliTRDReconstructor;

/**
 * @class AliHLTTRDTrackerV1Component
 * @brief A TRDTrackerV1 HLT processing component. 
 *
 * Uses the second generation TRD tracker AliTRDtrackerV1
*/

class AliHLTTRDTrackerV1Component : public AliHLTProcessor
    {
    public:
	AliHLTTRDTrackerV1Component();
	virtual ~AliHLTTRDTrackerV1Component();

	// Public functions to implement AliHLTComponent's interface.
	// These functions are required for the registration process

	const char* GetComponentID();
	void GetInputDataTypes( vector<AliHLTComponent_DataType>& list);
	AliHLTComponent_DataType GetOutputDataType();
	virtual void GetOutputDataSize( unsigned long& constBase, double& inputMultiplier );
	AliHLTComponent* Spawn();
	
    protected:
	AliHLTUInt32_t TransportTracks(TClonesArray *inTracksArray, AliHLTUInt8_t* output,
				       vector<AliHLTComponent_BlockData>& outputBlocks, AliHLTUInt32_t inOffset, AliHLTUInt32_t inSpec);
	
	// Protected functions to implement AliHLTComponent's interface.
	// These functions provide initialization as well as the actual processing
	// capabilities of the component. 

	int DoInit( int argc, const char** argv );
	int DoDeinit();
	int DoEvent( const AliHLTComponentEventData& evtData, 
		     const AliHLTComponentBlockData* blocks, 
		     AliHLTComponent_TriggerData& /*trigData*/, 
		     AliHLTUInt8_t* outputPtr, 
		     AliHLTUInt32_t& size, 
		     vector<AliHLTComponent_BlockData>& outputBlocks );
	using AliHLTProcessor::DoEvent;
	
    private:
	/** copy constructor prohibited */
	AliHLTTRDTrackerV1Component(const AliHLTTRDTrackerV1Component&);
	/** assignment operator prohibited */
	AliHLTTRDTrackerV1Component& operator=(const AliHLTTRDTrackerV1Component&);

	// The size of the output data produced, as a percentage of the input data's size.
	// Can be greater than 100 (%)
	unsigned fOutputPercentage; // Output volume in percentage of the input

	string fStrorageDBpath; // Default path for OCDB
	AliCDBManager *fCDB; //! Pointer to OCDB

	string fGeometryFileName; // Path to geometry file 
	
	AliTRDtrackerV1 *fTracker;//! Offline-pure/HLT tracker V1
	AliTRDrecoParam *fRecoParam; //! Offline reco params
	AliTRDReconstructor * fReconstructor;
	
	ClassDef(AliHLTTRDTrackerV1Component, 0)

    };
#endif
