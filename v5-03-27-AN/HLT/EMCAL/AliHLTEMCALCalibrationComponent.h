#ifndef ALIHLTEMCALCALIBRATIONCOMPONENT_H
#define ALIHLTEMCALCALIBRATIONCOMPONENT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTEMCALCalibrationComponent.h
    @author Timm Steinbeck, Matthias Richter
    @date   
    @brief  Declaration of a EMCALCalibration component. */


#include "AliHLTCalibrationProcessor.h"
class AliCDBManager;

/**
 * @class AliHLTEMCALCalibrationComponent
 * @brief A EMCALCalibration HLT processing component. 
 *
 * An implementiation of a EMCALCalibration component that just copies its input data
 * as a test, demonstration, and example of the HLT component scheme.
 * - @ref InitCalibration (optional)
 * - @ref ScanArgument (optional)
 * - @ref DeinitCalibration (optional)
 * - @ref ProcessCalibration
 * - @ref ShipDataToFXS
 * - @ref GetComponentID
 * - @ref GetInputDataTypes
 * - @ref GetOutputDataType
 * - @ref GetOutputDataSize
 * - @ref Spawn
 * @ingroup alihlt_tutorial
 */
class AliHLTEMCALCalibrationComponent : public AliHLTCalibrationProcessor
    {
    public:
	AliHLTEMCALCalibrationComponent();
	virtual ~AliHLTEMCALCalibrationComponent();

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

	virtual Int_t InitCalibration();
	virtual Int_t ScanArgument(int argc, const char** argv);
	virtual Int_t DeinitCalibration();
/* 	virtual Int_t ProcessCalibration(const AliHLTComponent_EventData& evtData, */
/* 					 const AliHLTComponent_BlockData* blocks, */
/* 					 AliHLTComponent_TriggerData& trigData, AliHLTUInt8_t* outputPtr, */
/* 					 AliHLTUInt32_t& size, */
/* 					 vector<AliHLTComponent_BlockData>& outputBlocks); */
	virtual Int_t ProcessCalibration( const AliHLTComponentEventData& evtData, AliHLTComponentTriggerData& trigData);
	virtual Int_t ShipDataToFXS( const AliHLTComponentEventData& evtData, AliHLTComponentTriggerData& trigData);

	using AliHLTCalibrationProcessor::ProcessCalibration;
	using AliHLTCalibrationProcessor::ShipDataToFXS;
	
    private:
	/** copy constructor prohibited */
	AliHLTEMCALCalibrationComponent(const AliHLTEMCALCalibrationComponent&);
	/** assignment operator prohibited */
	AliHLTEMCALCalibrationComponent& operator=(const AliHLTEMCALCalibrationComponent&);

	// The size of the output data produced, as a percentage of the input data's size.
	// Can be greater than 100 (%)
	
	unsigned fOutputPercentage; // Output volume in percentage of the input
	string fStrorageDBpath; // Default path for OCDB
	AliCDBManager *fCDB; //! Pointer to OCDB
	
	ClassDef(AliHLTEMCALCalibrationComponent, 0)

    };
#endif
