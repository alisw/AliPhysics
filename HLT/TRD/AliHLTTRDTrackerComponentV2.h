// $Id$

#ifndef ALIHLTTRDTRACKERCOMPONENTV2_H
#define ALIHLTTRDTRACKERCOMPONENTV2_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTTRDTrackerComponentV2.h
    @author Timm Steinbeck, Matthias Richter
    @date   
    @brief  Declaration of a TRDTracker component. */


#include "AliHLTProcessor.h"
class AliCDBManager;
class TFile;
class TGeoManager;
class AliTRDtrackerHLT;
class AliTRDtracker;
class AliMagFMaps;

/**
 * @class AliHLTTRDTrackerComponentV2
 * @brief A TRDTracker HLT processing component. 
 *
 * An implementiation of a TRDTracker component that just copies its input data
 * as a test, demonstration, and example of the HLT component scheme.
 * @ingroup alihlt_tutorial
 */
class AliHLTTRDTrackerComponentV2 : public AliHLTProcessor
    {
    public:
	AliHLTTRDTrackerComponentV2();
	virtual ~AliHLTTRDTrackerComponentV2();

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
	
    private:

	// The size of the output data produced, as a percentage of the input data's size.
	// Can be greater than 100 (%)
	unsigned fOutputPercentage; // Output volume in percentage of the input

	string fStrorageDBpath; // Default path for OCDB
	AliCDBManager *fCDB; //! Pointer to OCDB

	AliMagFMaps* fField; //! magn. field settings

	string fGeometryFileName; // Path to geometry file 
	TFile *fGeometryFile; //! // Pointer to the geom root file
	TGeoManager *fGeoManager; //! Pointer to geometry manager 

	AliTRDtrackerHLT *fTracker;//! Offline-like/HLT tracker
	//AliTRDtracker *fTracker;//! Offline-like/HLT tracker

	ClassDef(AliHLTTRDTrackerComponentV2, 0)

    };
#endif
