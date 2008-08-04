#ifndef ALIHLTMUONCLUSTERFINDERCOMPONENT_H
#define ALIHLTMUONCLUSTERFINDERCOMPONENT_H
/* This file is property of and copyright by the ALICE HLT Project        *
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/* $Id: $ */

///
///  @file   AliHLTMUONClusterFinderComponent.h
///  @author Artur Szostak <artursz@iafrica.com>
///  @date   28 May 2007
///  @brief  Cluster finder component for the dimuon HLT derived from offline code.
///

#include "AliHLTMUONProcessor.h"

#if __GNUC__ && __GNUC__ < 3
#define std
#endif

class AliRawReaderMemory;
class AliMUONDigitMaker;
class AliMUONGeometryTransformer;
class AliMUONCalibrationData;
class AliMUONDigitCalibrator;
class AliMUONVClusterFinder;
class AliMUONSimpleClusterServer;
class AliMUONRecoParam;

/**
 * @class AliHLTMUONClusterFinderComponent
 * @brief Cluster finding component for the dHLT tracker DDL raw data.
 *
 * This cluster finder component runs offline algorithms online within HLT.
 * It processes the raw DDL data from dimuon spectrometer tracker stations
 * and returns the cluster information in the offline cluster store format or
 * dHLT internal hit coordinate format.
 *
 * <h2>General properties:</h2>
 *
 * Component ID: \b MUONClusterFinder <br>
 * Library: \b libAliHLTMUON.so <br>
 * Input Data Types: ('DDL_RAW ', 'MUON') <br>
 * Output Data Types: ('CLUSTORE', 'MUON'), ('RECHITS ', 'MUON') <br>
 *
 * <h2>Mandatory arguments:</h2>
 * None <br>
 *
 * <h2>Optional arguments:</h2>
 * \li -delaysetup <br>
 *      Specifying this option causes the component to initialise from CDB only after
 *      receiving the first (normally Start-of-Run) event to process in DoEvent. <br>
 * \li -cdbpath <i>path</i> <br>
 *      Specifies the CDB path to use, given by <i>path</i>. This option will override
 *      the CDB path automatically set by the HLT framework. <br>
 * \li -run <i>number</i> <br>
 *      Specifies the run number to use, given by <i>number</i>. This option will
 *      override the current run number automatically set by the HLT framework. <br>
 * \li -tryrecover <br>
 *      This is a special option to the raw data decoder to turn on logic which will
 *      try and recover from corrupt raw DDL data. This is off by default. <br>
 *
 * <h2>Standard configuration:</h2>
 * This component should normally be configured with no extra options in the XML
 * configuration. <br>
 *
 * <h2>Default CDB entries:</h2>
 * The component loads electronics mapping and calibration information from the MUON
 * subdirectory in the CDB, MUON/Calib and MUON/Align.
 *
 * <h2>Performance:</h2>
 * A few tens of Hertz.
 *
 * <h2>Memory consumption:</h2>
 * A few MBytes.
 *
 * <h2>Output size:</h2>
 * Output size is about 25% of incoming raw input data for nominal p+p events.
 *
 * @ingroup alihlt_dimuon_component
 */
class AliHLTMUONClusterFinderComponent : public AliHLTMUONProcessor
{
public:
	AliHLTMUONClusterFinderComponent();
	virtual ~AliHLTMUONClusterFinderComponent();

	// Public functions to implement AliHLTComponent's interface.
	// These functions are required for the registration process

	virtual const char* GetComponentID();
	virtual void GetInputDataTypes(AliHLTComponentDataTypeList& list);
	virtual AliHLTComponentDataType GetOutputDataType();
	virtual int GetOutputDataTypes(AliHLTComponentDataTypeList& list);
	virtual void GetOutputDataSize(unsigned long& constBase, double& inputMultiplier);
	virtual AliHLTComponent* Spawn();
	
protected:
	
	// Protected functions to implement AliHLTComponent's interface.
	// These functions provide initialization as well as the actual processing
	// capabilities of the component.

	virtual int DoInit(int argc, const char** argv);
	virtual int Reconfigure(const char* cdbEntry, const char* componentId);
	virtual int ReadPreprocessorValues(const char* modules);
	virtual int DoDeinit();
	virtual int DoEvent(
			const AliHLTComponentEventData& evtData,
			const AliHLTComponentBlockData* blocks,
			AliHLTComponentTriggerData& trigData,
			AliHLTUInt8_t* outputPtr,
			AliHLTUInt32_t& size,
			AliHLTComponentBlockDataList& outputBlocks
		);
	
	using AliHLTProcessor::DoEvent;
	
private:

	// Do not allow copying of this class.
	/// Not implemented.
	AliHLTMUONClusterFinderComponent(const AliHLTMUONClusterFinderComponent& /*obj*/);
	/// Not implemented.
	AliHLTMUONClusterFinderComponent& operator = (const AliHLTMUONClusterFinderComponent& /*obj*/);
	
	void FreeObjects();

	int ReadConfigFromCDB(
			bool loadParams = true,
			bool loadMapping = true,
			bool loadGeom = true,
			bool loadCalib = true
		);
	
	AliRawReaderMemory* fRawReader;  ///< Raw reader for decoding input data for offline code.
	AliMUONDigitMaker* fDigitMaker;  ///< Digit maker to convert raw DDL data to digits.
	AliMUONGeometryTransformer* fTransformer;  ///< Transformer for geometry information.
	AliMUONCalibrationData* fCalibrationData;  ///< Calibration data object for calibrating digits.
	AliMUONDigitCalibrator* fDigitCalibrator;  ///< The object which performs calibration of digits.
	AliMUONVClusterFinder* fClusterFinder;  ///< Cluster finder implementing offline algorithm.
	AliMUONSimpleClusterServer* fClusterServer;  ///< Object for driving the cluster finder.
	AliMUONRecoParam* fRecoParam;  ///< Reconstruction parameters for offline code.
	
	bool fMakeClusterStore;  ///< Indicates if the cluster store object will be generated.
	bool fMakeRecHits;  ///< Indicates if the output will be in dHLT internal reconstructed hit format.
	
	ClassDef(AliHLTMUONClusterFinderComponent, 0) // Cluster finder component for dHLT running offline algorithms.
};

#endif // ALIHLTMUONCLUSTERFINDERCOMPONENT_H

