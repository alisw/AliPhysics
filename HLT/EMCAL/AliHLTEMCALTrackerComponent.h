#ifndef ALIHLTEMCALTRACKERCOMPONENT_H
#define ALIHLTEMCALTRACKERCOMPONENT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTEMCALTrackerComponent.h
    @author Mateusz Ploskon
    @date   
    @brief  Declaration of a EMCALTracker component. */


#include "AliHLTProcessor.h"
class TFolder;
class TFile;
class TGeoManager;

class AliCDBManager;
class AliMagFMaps;
class AliEMCALTracker;

/**
 * @class AliHLTEMCALTrackerComponent
 * @brief A EMCALTracker HLT processing component. 
 *
 * An implementiation of a EMCALTracker component that just copies its input data
 * as a test, demonstration, and example of the HLT component scheme.
 * @ingroup alihlt_tutorial
 */
class AliHLTEMCALTrackerComponent : public AliHLTProcessor
    {
    public:
	AliHLTEMCALTrackerComponent();
	virtual ~AliHLTEMCALTrackerComponent();

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
	AliHLTEMCALTrackerComponent(const AliHLTEMCALTrackerComponent&);
	/** assignment operator prohibited */
	AliHLTEMCALTrackerComponent& operator=(const AliHLTEMCALTrackerComponent&);

	// The size of the output data produced, as a percentage of the input data's size.
	// Can be greater than 100 (%)
	unsigned fOutputPercentage; // Output volume in percentage of the input

	string fStrorageDBpath; // Default path for OCDB
	AliCDBManager *fCDB; //! Pointer to OCDB

	AliMagFMaps* fField; //! magn. field settings

	string fGeometryFileName; // Path to geometry file 
	TFile *fGeometryFile; //! // Pointer to the geom root file
	TGeoManager *fGeoManager; //! Pointer to geometry manager 

	AliEMCALTracker *fTracker;//! Offline emcal tracker

	TFolder         *fInputFolder;//! input objects - convenient handling

	ClassDef(AliHLTEMCALTrackerComponent, 0)
    };
#endif
