// $Id$

#ifndef ALIHLTTRDTRACKERCOMPONENT_H
#define ALIHLTTRDTRACKERCOMPONENT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTTRDTrackerComponent.h
    @author Timm Steinbeck, Matthias Richter
    @date   
    @brief  Declaration of a TRDTracker component. */


#include "AliHLTProcessor.h"
class AliCDBManager;
class TFile;
class TGeoManager;
//class AliTRDtrackerHLT;
class AliTRDtracker;
class AliMagF;

/**
 * @class AliHLTTRDTrackerComponent
 * @brief A TRDTracker HLT processing component. 
 *
 * An implementiation of a TRDTracker component that just copies its input data
 * as a test, demonstration, and example of the HLT component scheme.
 * @ingroup alihlt_tutorial
 */
class AliHLTTRDTrackerComponent : public AliHLTProcessor
    {
    public:
	AliHLTTRDTrackerComponent();
	virtual ~AliHLTTRDTrackerComponent();

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
/* 	int DoEvent( const AliHLTComponent_EventData& evtData, const AliHLTComponent_BlockData* blocks,  */
/* 		     AliHLTComponent_TriggerData& trigData, AliHLTUInt8_t* outputPtr,  */
/* 		     AliHLTUInt32_t& size, vector<AliHLTComponent_BlockData>& outputBlocks ); */
	int DoEvent( const AliHLTComponentEventData & evtData,
		     AliHLTComponentTriggerData & trigData );

	using AliHLTProcessor::DoEvent;
	
    private:
	/** copy constructor prohibited */
	AliHLTTRDTrackerComponent(const AliHLTTRDTrackerComponent&);
	/** assignment operator prohibited */
	AliHLTTRDTrackerComponent& operator=(const AliHLTTRDTrackerComponent&);

	// The size of the output data produced, as a percentage of the input data's size.
	// Can be greater than 100 (%)
	unsigned fOutputPercentage; // Output volume in percentage of the input

	string fStrorageDBpath; // Default path for OCDB
	AliCDBManager *fCDB; //! Pointer to OCDB

	string fGeometryFileName; // Path to geometry file 
	TFile *fGeometryFile; //! // Pointer to the geom root file
	TGeoManager *fGeoManager; //! Pointer to geometry manager 

	//AliTRDtrackerHLT *fTracker;//! Offline-like/HLT tracker
	AliTRDtracker *fTracker;//! Offline-pure/HLT tracker

	ClassDef(AliHLTTRDTrackerComponent, 0)

    };
#endif
