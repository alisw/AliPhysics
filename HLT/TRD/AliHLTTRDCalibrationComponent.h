#ifndef ALIHLTTRDCALIBRATIONCOMPONENT_H
#define ALIHLTTRDCALIBRATIONCOMPONENT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTTRDCalibrationComponent.h
    @author Timm Steinbeck, Matthias Richter
    @date   
    @brief  Declaration of a TRDCalibration component. */


#include "AliHLTCalibrationProcessor.h"
class AliCDBManager;
class AliTRDCalibraFillHisto;

/**
 * @class AliHLTTRDCalibrationComponent
 * @brief A TRDCalibration HLT processing component. 
 *
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
class AliHLTTRDCalibrationComponent : public AliHLTCalibrationProcessor
    {
    public:
	AliHLTTRDCalibrationComponent();
	virtual ~AliHLTTRDCalibrationComponent();

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
	
	AliTRDCalibraFillHisto *fTRDCalibraFillHisto;
	
	virtual Int_t InitCalibration();
	virtual Int_t ScanArgument(int argc, const char** argv);
	virtual Int_t DeinitCalibration();
	virtual Int_t ProcessCalibration(const AliHLTComponent_EventData& evtData,
					 const AliHLTComponent_BlockData* blocks,
					 AliHLTComponent_TriggerData& trigData, AliHLTUInt8_t* outputPtr,
					 AliHLTUInt32_t& size,
					 vector<AliHLTComponent_BlockData>& outputBlocks);
/* 	virtual Int_t ProcessCalibration( const AliHLTComponentEventData& evtData, AliHLTComponentTriggerData& trigData); */
	virtual Int_t ShipDataToFXS( const AliHLTComponentEventData& evtData, AliHLTComponentTriggerData& trigData);
	
	using AliHLTCalibrationProcessor::ProcessCalibration;
	using AliHLTCalibrationProcessor::ShipDataToFXS;

	Bool_t fUseHLTTracks;
	Int_t ReadTracks(TClonesArray *outArray, void* inputPtr, AliHLTUInt32_t size);
	
    private:
	/** copy constructor prohibited */
	AliHLTTRDCalibrationComponent(const AliHLTTRDCalibrationComponent&);
	/** assignment operator prohibited */
	AliHLTTRDCalibrationComponent& operator=(const AliHLTTRDCalibrationComponent&);
	TObjArray* FormOutput();
	// The size of the output data produced, as a percentage of the input data's size.
	// Can be greater than 100 (%)
	
	unsigned fOutputPercentage; // Output volume in percentage of the input
	string fStrorageDBpath; // Default path for OCDB
	AliCDBManager *fCDB; //! Pointer to OCDB
	
	ClassDef(AliHLTTRDCalibrationComponent, 0)

    };
#endif
